#![forbid(unsafe_code)]

//! Squeeze module implementing optimal LZ77 parsing.

use crate::lz77::{
    find_longest_match, lz77_greedy, Lz77Store, Lz77StoreView,
    UninitializedLz77Store, ZopfliBlockState,
};
use crate::symbols::{
    get_dist_extra_bits, get_dist_symbol, get_dist_symbol_extra_bits, get_length_extra_bits,
    get_length_symbol, get_length_symbol_extra_bits,
};
use crate::tree::{calculate_bit_lengths, calculate_entropy};
use crate::util::{
    LARGE_FLOAT, MAX_CL_BIT_LENGTH, NUM_D, NUM_LL, RLE_CODE_COPY, RLE_CODE_ZERO_11_138,
    RLE_CODE_ZERO_3_10, RLE_COPY_MAX, RLE_COPY_MIN, RLE_COPY_THRESHOLD, RLE_ZERO_17_MAX,
    RLE_ZERO_17_MIN, RLE_ZERO_18_MAX, RLE_ZERO_18_MIN,
};
use std::cmp;

/// Statistics of literal/length and distance symbols.
#[derive(Clone, Debug)]
pub struct SymbolStats {
    /// Frequencies of literal and length symbols.
    pub litlens: [usize; NUM_LL],
    /// Frequencies of distance symbols.
    pub dists: [usize; NUM_D],
    /// Calculated cost (in bits) for each literal/length symbol.
    pub ll_symbols: [f64; NUM_LL],
    /// Calculated cost (in bits) for each distance symbol.
    pub d_symbols: [f64; NUM_D],
}

impl SymbolStats {
    /// Creates a new `SymbolStats` with all frequencies and costs initialized to zero.
    pub fn new() -> Self {
        Self {
            litlens: [0; NUM_LL],
            dists: [0; NUM_D],
            ll_symbols: [0.0; NUM_LL],
            d_symbols: [0.0; NUM_D],
        }
    }

    /// Resets the literal/length and distance frequencies to zero.
    pub fn clear_freqs(&mut self) {
        self.litlens.fill(0);
        self.dists.fill(0);
    }

    /// Copies frequencies and costs from a source `SymbolStats`.
    pub fn copy_stats(&mut self, source: &SymbolStats) {
        self.litlens.copy_from_slice(&source.litlens);
        self.dists.copy_from_slice(&source.dists);
        self.ll_symbols.copy_from_slice(&source.ll_symbols);
        self.d_symbols.copy_from_slice(&source.d_symbols);
    }

    /// Adds the weighed frequencies of another `SymbolStats` to this one.
    pub fn add_weighed(&mut self, other: &Self, w2: f64) {
        for i in 0..NUM_LL {
            self.litlens[i] = (self.litlens[i] as f64 + other.litlens[i] as f64 * w2) as usize;
        }
        for i in 0..NUM_D {
            self.dists[i] = (self.dists[i] as f64 + other.dists[i] as f64 * w2) as usize;
        }
        self.litlens[256] = 1; /* End symbol. */
    }

    /// Calculates entropy/costs for the current symbol frequencies.
    pub fn calculate_statistics(&mut self) {
        let _ = calculate_entropy(&self.litlens, &mut self.ll_symbols);
        let _ = calculate_entropy(&self.dists, &mut self.d_symbols);
    }

    /// Collects symbol frequencies from the given `Lz77Store` and calculates statistics.
    pub fn get_statistics(&mut self, store: &Lz77StoreView<'_>) {
        for i in 0..store.litlens.len() {
            let dist = store.dists[i] as usize;
            let litlen = store.litlens[i] as usize;
            if dist == 0 {
                self.litlens[litlen] += 1;
            } else {
                self.litlens[get_length_symbol(litlen) as usize] += 1;
                self.dists[get_dist_symbol(dist) as usize] += 1;
            }
        }
        self.litlens[256] = 1; /* End symbol. */
        self.calculate_statistics();
    }

    /// Randomizes frequencies using the provided random state to avoid local minima.
    pub fn randomize_freqs(&mut self, ran_state: &mut RanState) {
        let n_ll = self.litlens.len();
        for i in 0..n_ll {
            if (ran_state.ran() >> 4) % 3 == 0 {
                let rand_idx = ran_state.ran() as usize % n_ll;
                self.litlens[i] = self.litlens[rand_idx];
            }
        }
        let n_d = self.dists.len();
        for i in 0..n_d {
            if (ran_state.ran() >> 4) % 3 == 0 {
                let rand_idx = ran_state.ran() as usize % n_d;
                self.dists[i] = self.dists[rand_idx];
            }
        }
        self.litlens[256] = 1; /* End symbol. */
    }
}

impl Default for SymbolStats {
    fn default() -> Self {
        Self::new()
    }
}

/// Random state container using a simple multiply-with-carry generator.
#[derive(Debug, Clone, Copy)]
pub struct RanState {
    /// Internal state parameter W.
    pub m_w: u32,
    /// Internal state parameter Z.
    pub m_z: u32,
}

impl RanState {
    /// Creates a new `RanState` with default seed values.
    pub fn new() -> Self {
        Self { m_w: 1, m_z: 2 }
    }

    /// Generates a pseudo-random 32-bit unsigned integer and updates the internal state.
    pub fn ran(&mut self) -> u32 {
        self.m_z = 36969 * (self.m_z & 65535) + (self.m_z >> 16);
        self.m_w = 18000 * (self.m_w & 65535) + (self.m_w >> 16);
        (self.m_z << 16).wrapping_add(self.m_w)
    }
}

impl Default for RanState {
    fn default() -> Self {
        Self::new()
    }
}

/// Cost model representing how the bit cost of a literal, length, or distance is calculated.
#[derive(Debug, Clone, Copy)]
pub enum CostModel<'a> {
    /// Fixed Huffman tree cost model.
    Fixed,
    /// Dynamic cost model based on custom symbol statistics.
    Stat(&'a SymbolStats),
}

impl<'a> CostModel<'a> {
    /// Calculates the bit cost of a literal, length, or distance.

    pub fn get_cost(&self, litlen: usize, dist: usize) -> f64 {
        match self {
            CostModel::Fixed => {
                const FIXED_LITLEN_LIMIT: usize = 143;
                const FIXED_LITLEN_SHORT_COST: f64 = 8.0;
                const FIXED_LITLEN_LONG_COST: f64 = 9.0;
                const FIXED_LSYM_LIMIT: usize = 279;
                const FIXED_LSYM_SHORT_COST: f64 = 7.0;
                const FIXED_LSYM_LONG_COST: f64 = 8.0;
                const FIXED_DIST_SYMBOL_COST: f64 = 5.0;

                if dist == 0 {
                    if litlen <= FIXED_LITLEN_LIMIT {
                        FIXED_LITLEN_SHORT_COST
                    } else {
                        FIXED_LITLEN_LONG_COST
                    }
                } else {
                    let dbits = get_dist_extra_bits(dist) as f64;
                    let lbits = get_length_extra_bits(litlen) as f64;
                    let lsym = get_length_symbol(litlen) as usize;
                    let mut cost = 0.0;
                    if lsym <= FIXED_LSYM_LIMIT {
                        cost += FIXED_LSYM_SHORT_COST;
                    } else {
                        cost += FIXED_LSYM_LONG_COST;
                    }
                    cost += FIXED_DIST_SYMBOL_COST; /* Every dist symbol has length 5. */
                    cost + dbits + lbits
                }
            }
            CostModel::Stat(stats) => {
                if dist == 0 {
                    stats.ll_symbols[litlen]
                } else {
                    let lsym = get_length_symbol(litlen) as usize;
                    let lbits = get_length_extra_bits(litlen) as f64;
                    let dsym = get_dist_symbol(dist) as usize;
                    let dbits = get_dist_extra_bits(dist) as f64;
                    lbits + dbits + stats.ll_symbols[lsym] + stats.d_symbols[dsym]
                }
            }
        }
    }

    /// Returns the minimum cost among all literal/length/distance choices.
    pub fn get_min_cost(&self) -> f64 {
        let mut mincost: f64;
        let mut bestlength = 0;
        let mut bestdist = 0;

        const D_SYMBOLS: [usize; 30] = [
            1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025,
            1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577,
        ];

        mincost = LARGE_FLOAT;
        for i in 3..259 {
            let c = self.get_cost(i, 1);
            if c < mincost {
                bestlength = i;
                mincost = c;
            }
        }

        mincost = LARGE_FLOAT;
        for i in 0..30 {
            let c = self.get_cost(3, D_SYMBOLS[i]);
            if c < mincost {
                bestdist = D_SYMBOLS[i];
                mincost = c;
            }
        }

        self.get_cost(bestlength, bestdist)
    }
}

/// Adjusts the distance code lengths to satisfy the requirements of some older or buggy deflate
/// decoders
/// that require at least two distance codes to be defined (even if only one distance code is used).
pub fn patch_distance_codes_for_buggy_decoders(d_lengths: &mut [u32]) {
    let mut num_dist_codes = 0;
    for i in 0..30 {
        if d_lengths[i] > 0 {
            num_dist_codes += 1;
        }
        if num_dist_codes >= 2 {
            return;
        }
    }

    if num_dist_codes == 0 {
        d_lengths[0] = 1;
        d_lengths[1] = 1;
    } else if num_dist_codes == 1 {
        let idx = if d_lengths[0] > 0 { 1 } else { 0 };
        d_lengths[idx] = 1;
    }
}

fn encode_tree(
    ll_lengths: &[u32],
    d_lengths: &[u32],
    use_16: bool,
    use_17: bool,
    use_18: bool,
) -> usize {
    let mut hlit = 29;
    while hlit > 0 && ll_lengths[257 + hlit - 1] == 0 {
        hlit -= 1;
    }
    let mut hdist = 29;
    while hdist > 0 && d_lengths[1 + hdist - 1] == 0 {
        hdist -= 1;
    }
    let hlit2 = hlit + 257;
    let lld_total = hlit2 + hdist + 1;

    let mut clcounts = [0usize; 19];
    let mut i = 0;
    while i < lld_total {
        let symbol = if i < hlit2 { ll_lengths[i] } else { d_lengths[i - hlit2] };
        let mut count = 1;
        if use_16 || (symbol == 0 && (use_17 || use_18)) {
            let mut j = i + 1;
            while j < lld_total {
                let next_symbol = if j < hlit2 { ll_lengths[j] } else { d_lengths[j - hlit2] };
                if symbol == next_symbol {
                    count += 1;
                } else {
                    break;
                }
                j += 1;
            }
        }
        i += count;

        let mut count = count;
        // Repetitions of zeroes
        if symbol == 0 && count >= RLE_ZERO_17_MIN {
            if use_18 {
                while count >= RLE_ZERO_18_MIN {
                    let count2 = cmp::min(RLE_ZERO_18_MAX, count);
                    clcounts[RLE_CODE_ZERO_11_138] += 1;
                    count -= count2;
                }
            }
            if use_17 {
                while count >= RLE_ZERO_17_MIN {
                    let count2 = cmp::min(RLE_ZERO_17_MAX, count);
                    clcounts[RLE_CODE_ZERO_3_10] += 1;
                    count -= count2;
                }
            }
        }

        // Repetitions of any symbol
        if use_16 && count >= RLE_COPY_THRESHOLD {
            count -= 1; // Since the first one is hardcoded.
            clcounts[symbol as usize] += 1;
            while count >= RLE_COPY_MIN {
                let count2 = cmp::min(RLE_COPY_MAX, count);
                clcounts[RLE_CODE_COPY] += 1;
                count -= count2;
            }
        }

        clcounts[symbol as usize] += count;
    }

    let mut clcl = [0u32; 19];
    let _ = calculate_bit_lengths(&clcounts, MAX_CL_BIT_LENGTH, &mut clcl);

    let mut hclen = 15;
    const ORDER: [usize; 19] = [16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15];
    while hclen > 0 && clcounts[ORDER[hclen + 4 - 1]] == 0 {
        hclen -= 1;
    }

    let mut result_size = 14; // hlit, hdist, hclen bits
    result_size += (hclen + 4) * 3; // clcl bits
    for j in 0..19 {
        result_size += (clcl[j] as usize) * clcounts[j];
    }
    result_size += clcounts[16] * 2;
    result_size += clcounts[17] * 3;
    result_size += clcounts[18] * 7;

    result_size
}

/// Calculates the size of the dynamic Huffman tree in bits.
pub fn calculate_tree_size(ll_lengths: &[u32], d_lengths: &[u32]) -> usize {
    let mut result = 0;
    for i in 0..8 {
        let size = encode_tree(ll_lengths, d_lengths, (i & 1) != 0, (i & 2) != 0, (i & 4) != 0);
        if result == 0 || size < result {
            result = size;
        }
    }
    result
}

/// Helper to calculate the bit size of literal/length/distance symbols in a range
/// using a simple loop. Used for small inputs or as a fallback.
pub fn calculate_block_symbol_size_small(
    ll_lengths: &[u32],
    d_lengths: &[u32],
    store: &Lz77StoreView<'_>,
    lstart: usize,
    lend: usize,
) -> usize {
    let mut result = 0;
    for i in lstart..lend {
        let dist = store.dists[i];
        let litlen = store.litlens[i];
        if dist == 0 {
            let idx = litlen as usize;
            if idx < NUM_LL {
                result += ll_lengths[idx] as usize;
            }
        } else {
            let ll_symbol = get_length_symbol(litlen as usize) as usize;
            let d_symbol = get_dist_symbol(dist as usize) as usize;
            result += ll_lengths[ll_symbol] as usize;
            result += d_lengths[d_symbol] as usize;
            result += get_length_symbol_extra_bits(ll_symbol) as usize;
            result += get_dist_symbol_extra_bits(d_symbol) as usize;
        }
    }
    result += ll_lengths[256] as usize; // end symbol
    result
}

/// Calculates the bit size of literal/length/distance symbols in a range,
/// using precalculated histogram counts when the size of the range is large.
pub fn calculate_block_symbol_size_given_counts(
    ll_counts: &[usize],
    d_counts: &[usize],
    ll_lengths: &[u32],
    d_lengths: &[u32],
    store: &Lz77StoreView<'_>,
    lstart: usize,
    lend: usize,
) -> usize {
    if lstart + NUM_LL * 3 > lend {
        calculate_block_symbol_size_small(ll_lengths, d_lengths, store, lstart, lend)
    } else {
        let mut result = 0;
        for i in 0..256 {
            result += (ll_lengths[i] as usize) * ll_counts[i];
        }
        for i in 257..286 {
            result += (ll_lengths[i] as usize) * ll_counts[i];
            result += get_length_symbol_extra_bits(i) as usize * ll_counts[i];
        }
        for i in 0..30 {
            result += (d_lengths[i] as usize) * d_counts[i];
            result += get_dist_symbol_extra_bits(i) as usize * d_counts[i];
        }
        result += ll_lengths[256] as usize; // end symbol
        result
    }
}

/// Calculates the bit size of literal/length/distance symbols in a range,
/// dynamically computing a histogram if the range is large.
pub fn calculate_block_symbol_size(
    ll_lengths: &[u32],
    d_lengths: &[u32],
    store: &Lz77StoreView<'_>,
    lstart: usize,
    lend: usize,
) -> usize {
    if lstart + NUM_LL * 3 > lend {
        calculate_block_symbol_size_small(ll_lengths, d_lengths, store, lstart, lend)
    } else {
        let mut ll_counts = [0usize; NUM_LL];
        let mut d_counts = [0usize; NUM_D];
        store.get_histogram(lstart, lend, &mut ll_counts, &mut d_counts);
        calculate_block_symbol_size_given_counts(
            &ll_counts, &d_counts, ll_lengths, d_lengths, store, lstart, lend,
        )
    }
}

/// Smooths out symbol counts/frequencies to make them more suitable for RLE encoding.
pub fn optimize_huffman_for_rle(counts: &mut [usize]) {
    const RLE_ZERO_STRIDE_THRESHOLD: usize = 5;
    const RLE_NONZERO_STRIDE_THRESHOLD: usize = 7;
    let mut length = counts.len();
    while length > 0 {
        if counts[length - 1] != 0 {
            break;
        }
        length -= 1;
    }
    if length == 0 {
        return;
    }

    let mut good_for_rle = vec![false; length];
    let mut symbol = counts[0];
    let mut stride = 0;
    for i in 0..length + 1 {
        if i == length || counts[i] != symbol {
            if (symbol == 0 && stride >= RLE_ZERO_STRIDE_THRESHOLD)
                || (symbol != 0 && stride >= RLE_NONZERO_STRIDE_THRESHOLD)
            {
                for k in 0..stride {
                    good_for_rle[i - k - 1] = true;
                }
            }
            stride = 1;
            if i != length {
                symbol = counts[i];
            }
        } else {
            stride += 1;
        }
    }

    stride = 0;
    let mut limit = counts[0];
    let mut sum = 0;
    for i in 0..length + 1 {
        let abs_diff = if i == length {
            0
        } else if counts[i] > limit {
            counts[i] - limit
        } else {
            limit - counts[i]
        };

        if i == length || good_for_rle[i] || abs_diff >= 4 {
            if stride >= 4 || (stride >= 3 && sum == 0) {
                let mut count = (sum + stride / 2) / stride;
                if count < 1 {
                    count = 1;
                }
                if sum == 0 {
                    count = 0;
                }
                for k in 0..stride {
                    counts[i - k - 1] = count;
                }
            }
            stride = 0;
            sum = 0;
            if i < length.saturating_sub(3) {
                limit = (counts[i] + counts[i + 1] + counts[i + 2] + counts[i + 3] + 2) / 4;
            } else if i < length {
                limit = counts[i];
            } else {
                limit = 0;
            }
        }
        stride += 1;
        if i != length {
            sum += counts[i];
        }
    }
}

/// Attempts to optimize the Huffman trees for RLE encoding and returns the new size.
///
/// If the optimized tree yields a smaller total size (data + tree representation),
/// the `ll_lengths` and `d_lengths` arrays are updated with the new lengths.
pub fn try_optimize_huffman_for_rle(
    store: &Lz77StoreView<'_>,
    lstart: usize,
    lend: usize,
    ll_counts: &[usize],
    d_counts: &[usize],
    ll_lengths: &mut [u32],
    d_lengths: &mut [u32],
) -> f64 {
    let mut ll_counts2 = [0usize; NUM_LL];
    let mut d_counts2 = [0usize; NUM_D];
    let mut ll_lengths2 = [0u32; NUM_LL];
    let mut d_lengths2 = [0u32; NUM_D];

    let treesize = calculate_tree_size(ll_lengths, d_lengths) as f64;
    let datasize = calculate_block_symbol_size_given_counts(
        ll_counts, d_counts, ll_lengths, d_lengths, store, lstart, lend,
    ) as f64;

    ll_counts2.copy_from_slice(ll_counts);
    d_counts2.copy_from_slice(d_counts);
    optimize_huffman_for_rle(&mut ll_counts2);
    optimize_huffman_for_rle(&mut d_counts2);

    let _ = calculate_bit_lengths(&ll_counts2, 15, &mut ll_lengths2);
    let _ = calculate_bit_lengths(&d_counts2, 15, &mut d_lengths2);
    patch_distance_codes_for_buggy_decoders(&mut d_lengths2);

    let treesize2 = calculate_tree_size(&ll_lengths2, &d_lengths2) as f64;
    let datasize2 = calculate_block_symbol_size_given_counts(
        ll_counts,
        d_counts,
        &ll_lengths2,
        &d_lengths2,
        store,
        lstart,
        lend,
    ) as f64;

    if treesize2 + datasize2 < treesize + datasize {
        ll_lengths.copy_from_slice(&ll_lengths2);
        d_lengths.copy_from_slice(&d_lengths2);
        treesize2 + datasize2
    } else {
        treesize + datasize
    }
}

/// Calculates the Huffman bit lengths using dynamic trees and returns the total size in bits.
pub fn get_dynamic_lengths(
    store: &Lz77StoreView<'_>,
    lstart: usize,
    lend: usize,
    ll_lengths: &mut [u32],
    d_lengths: &mut [u32],
) -> f64 {
    let mut ll_counts = [0usize; NUM_LL];
    let mut d_counts = [0usize; NUM_D];

    store.get_histogram(lstart, lend, &mut ll_counts, &mut d_counts);
    ll_counts[256] = 1; /* End symbol. */

    let _ = calculate_bit_lengths(&ll_counts, 15, ll_lengths);
    let _ = calculate_bit_lengths(&d_counts, 15, d_lengths);
    patch_distance_codes_for_buggy_decoders(d_lengths);

    try_optimize_huffman_for_rle(store, lstart, lend, &ll_counts, &d_counts, ll_lengths, d_lengths)
}

/// Calculates the size of a block in bits for a given block type (uncompressed, fixed, or dynamic).
pub fn calculate_block_size(
    store: &Lz77StoreView<'_>,
    lstart: usize,
    lend: usize,
    btype: i32,
) -> f64 {
    let mut ll_lengths = [0u32; NUM_LL];
    let mut d_lengths = [0u32; NUM_D];

    let mut result = 3.0; // bfinal and btype bits

    if btype == 0 {
        let length = lz77_get_byte_range(store, lstart, lend);
        let rem = length % 65535;
        let blocks = length / 65535 + if rem > 0 { 1 } else { 0 };
        return (blocks * 5 * 8 + length * 8) as f64;
    } else if btype == 1 {
        get_fixed_tree(&mut ll_lengths, &mut d_lengths);
        result += calculate_block_symbol_size(&ll_lengths, &d_lengths, store, lstart, lend) as f64;
    } else {
        result += get_dynamic_lengths(store, lstart, lend, &mut ll_lengths, &mut d_lengths);
    }

    result
}

/// Calculates the byte range spanned by the LZ77 commands in the store.
pub fn lz77_get_byte_range(store: &Lz77StoreView<'_>, lstart: usize, lend: usize) -> usize {
    if lstart == lend {
        return 0;
    }
    let l = lend - 1;
    let elem_len = if store.dists[l] == 0 { 1 } else { store.litlens[l] as usize };
    (store.pos[l] as usize + elem_len) - store.pos[lstart] as usize
}

/// Populates the bit length arrays for a fixed Huffman tree.
pub fn get_fixed_tree(ll_lengths: &mut [u32], d_lengths: &mut [u32]) {
    ll_lengths[0..144].fill(8);
    ll_lengths[144..256].fill(9);
    ll_lengths[256..280].fill(7);
    ll_lengths[280..288].fill(8);
    d_lengths[0..32].fill(5);
}

/// Finds the best path of literal/length/distance choices using the given cost model,
/// updating `length_array` and returning the cost of the optimal path.
pub fn get_best_lengths<'a>(
    s: &mut ZopfliBlockState<'a>,
    in_data: &[u8],
    instart: usize,
    inend: usize,
    cost_model: &CostModel,
    length_array: &mut [u16],
) -> f64 {
    let blocksize = inend - instart;
    if instart == inend {
        return 0.0;
    }

    let mut costs = vec![LARGE_FLOAT as f32; blocksize + 1];
    costs[0] = 0.0;
    length_array[0] = 0;

    let window_size = crate::util::WINDOW_SIZE;
    let mut hash = crate::hash::ZopfliHash::new(window_size);
    let windowstart = if instart > window_size { instart - window_size } else { 0 };

    hash.warmup(in_data, windowstart, inend);
    for i in windowstart..instart {
        hash.update(in_data, i, inend);
    }

    let mincost = cost_model.get_min_cost();

    let mut i = instart;
    while i < inend {
        let mut j = i - instart;
        hash.update(in_data, i, inend);

        // Shortcut long repetitions
        let same_idx = i & crate::util::WINDOW_MASK;
        if hash.same[same_idx] > (crate::util::MAX_MATCH as u16) * 2
            && i > instart + crate::util::MAX_MATCH + 1
            && i + crate::util::MAX_MATCH * 2 + 1 < inend
            && hash.same[(i - crate::util::MAX_MATCH) & crate::util::WINDOW_MASK]
                > crate::util::MAX_MATCH as u16
        {
            let symbolcost = cost_model.get_cost(crate::util::MAX_MATCH, 1);
            for _k in 0..crate::util::MAX_MATCH {
                let next_j = j + crate::util::MAX_MATCH;
                costs[next_j] = (costs[j] as f64 + symbolcost) as f32;
                length_array[next_j] = crate::util::MAX_MATCH as u16;
                i += 1;
                j += 1;
                hash.update(in_data, i, inend);
            }
        }

        let mut sublen = [0u16; 259];
        let res =
            find_longest_match(s, &hash, in_data, i, crate::util::MAX_MATCH, Some(&mut sublen));
        let _dist = res.0;
        let leng = res.1;

        /* Literal. */
        if i + 1 <= inend {
            let new_cost = cost_model.get_cost(in_data[i] as usize, 0) + costs[j] as f64;
            if new_cost < costs[j + 1] as f64 {
                costs[j + 1] = new_cost as f32;
                length_array[j + 1] = 1;
            }
        }

        /* Lengths. */
        let kend = cmp::min(leng as usize, inend - i);
        let mincostaddcostj = mincost + costs[j] as f64;
        for k in 3..kend + 1 {
            if costs[j + k] as f64 <= mincostaddcostj {
                continue;
            }

            let new_cost = cost_model.get_cost(k, sublen[k] as usize) + costs[j] as f64;
            if new_cost < costs[j + k] as f64 {
                debug_assert!(k <= crate::util::MAX_MATCH);
                costs[j + k] = new_cost as f32;
                length_array[j + k] = k as u16;
            }
        }

        i += 1;
    }

    debug_assert!(costs[blocksize] >= 0.0);
    costs[blocksize] as f64
}

/// Traces backward from the end of the block using the `length_array` to reconstruct the path.
pub fn trace_backwards(size: usize, length_array: &[u16]) -> Vec<u16> {
    let mut path = Vec::new();
    if size == 0 {
        return path;
    }
    let mut index = size;
    loop {
        let length = length_array[index];
        path.push(length);
        debug_assert!(length as usize <= index);
        debug_assert!(length as usize <= crate::util::MAX_MATCH);
        debug_assert!(length != 0);
        index -= length as usize;
        if index == 0 {
            break;
        }
    }
    path.reverse();
    path
}

/// Follows the reconstructed path of lengths to generate the actual LZ77 commands
/// and append them to the store.
pub fn follow_path<'a>(
    s: &mut ZopfliBlockState<'a>,
    in_data: &[u8],
    instart: usize,
    inend: usize,
    path: &[u16],
    store: &mut Lz77Store<'a>,
) {
    if instart == inend {
        return;
    }

    let window_size = crate::util::WINDOW_SIZE;
    let mut hash = crate::hash::ZopfliHash::new(window_size);
    let windowstart = if instart > window_size { instart - window_size } else { 0 };

    hash.warmup(in_data, windowstart, inend);
    for i in windowstart..instart {
        hash.update(in_data, i, inend);
    }

    let mut pos = instart;
    for &length in path {
        debug_assert!(pos < inend);
        hash.update(in_data, pos, inend);

        let mut final_length = length;

        if length >= crate::util::MIN_MATCH as u16 {
            let res = find_longest_match(s, &hash, in_data, pos, length as usize, None);
            let dist = res.0;
            let sink_length = res.1;
            debug_assert!(!(sink_length != length && length > 2 && sink_length > 2));

            #[cfg(debug_assertions)]
            {
                let d = dist as usize;
                let l = length as usize;
                debug_assert!(pos >= d);
                debug_assert!(pos + l <= in_data.len());
                debug_assert_eq!(&in_data[pos - d..pos - d + l], &in_data[pos..pos + l]);
            }

            store.store_lit_len_dist(length, dist, pos);
        } else {
            final_length = 1;
            store.store_lit_len_dist(in_data[pos] as u16, 0, pos);
        }

        debug_assert!(pos + final_length as usize <= inend);
        for j in 1..final_length as usize {
            hash.update(in_data, pos + j, inend);
        }

        pos += final_length as usize;
    }
}

/// Performs optimal LZ77 parse using a fixed cost model.
pub fn lz77_optimal_fixed<'a>(
    s: &mut ZopfliBlockState<'a>,
    in_data: &'a [u8],
    instart: usize,
    inend: usize,
    store: &mut Lz77Store<'a>,
) {
    let blocksize = inend - instart;
    let mut length_array = vec![0u16; blocksize + 1];

    s.blockstart = instart;
    s.blockend = inend;

    let cost_model = CostModel::Fixed;
    let _cost = get_best_lengths(s, in_data, instart, inend, &cost_model, &mut length_array);
    let path = trace_backwards(blocksize, &length_array);
    follow_path(s, in_data, instart, inend, &path, store);
}

/// Performs optimal LZ77 parse using a dynamic cost model and iterative statistical refinement.
pub fn lz77_optimal<'a>(
    s: &mut ZopfliBlockState<'a>,
    in_data: &'a [u8],
    instart: usize,
    inend: usize,
    numiterations: usize,
    store: &mut Lz77Store<'a>,
) {
    s.blockstart = instart;
    s.blockend = inend;

    let blocksize = inend - instart;
    let mut length_array = vec![0u16; blocksize + 1];

    let mut ran_state = RanState::default();
    let mut stats = SymbolStats::default();
    let mut beststats = SymbolStats::default();

    let mut currentstore = UninitializedLz77Store::new().initialize(in_data);

    /* Initial run. */
    lz77_greedy(s, in_data, instart, inend, &mut currentstore);
    stats.get_statistics(&currentstore.as_view());

    let mut bestcost = LARGE_FLOAT;
    let mut lastcost = 0.0;
    let mut lastrandomstep = -1;

    for i in 0..numiterations {
        // Clear current store
        currentstore.clear();

        let cost_model = CostModel::Stat(&stats);
        let _best_len_cost =
            get_best_lengths(s, in_data, instart, inend, &cost_model, &mut length_array);
        let path = trace_backwards(blocksize, &length_array);
        follow_path(s, in_data, instart, inend, &path, &mut currentstore);

        let cost = calculate_block_size(
            &currentstore.as_view(),
            /* lstart= */ 0,
            currentstore.litlens.len(),
            /* btype= */ 2,
        );

        if s.options.verbose_more || (s.options.verbose && cost < bestcost) {
            eprintln!("Iteration {}: {} bit", i, cost as i32);
        }

        if cost < bestcost {
            /* Copy to the output store. */
            *store = currentstore.clone();
            beststats.copy_stats(&stats);
            bestcost = cost;
        }

        let laststats = stats.clone();
        stats.clear_freqs();
        stats.get_statistics(&currentstore.as_view());

        if lastrandomstep != -1 {
            /* This makes it converge slower but better. Do it only once the
            randomness kicks in so that if the user does few iterations, it gives a
            better result sooner. */
            stats.add_weighed(&laststats, 0.5);
            stats.calculate_statistics();
        }

        if i > 5 && cost == lastcost {
            stats.copy_stats(&beststats);
            stats.randomize_freqs(&mut ran_state);
            stats.calculate_statistics();
            lastrandomstep = i as i32;
        }

        lastcost = cost;
    }
}

/// Estimates the optimal block size for a chunk of LZ77 data.
/// Automatically evaluates uncompressed, fixed, and dynamic block costs.
pub fn calculate_block_size_auto_type(
    store: &Lz77StoreView<'_>,
    lstart: usize,
    lend: usize,
) -> f64 {
    let uncompressedcost = calculate_block_size(store, lstart, lend, /* btype= */ 0);
    // Don't do the expensive fixed cost calculation for larger blocks that are
    // unlikely to use it.
    let fixedcost = if store.litlens.len() > 1000 {
        uncompressedcost
    } else {
        calculate_block_size(store, lstart, lend, /* btype= */ 1)
    };
    let dyncost = calculate_block_size(store, lstart, lend, /* btype= */ 2);

    if uncompressedcost < fixedcost && uncompressedcost < dyncost {
        uncompressedcost
    } else if fixedcost < dyncost {
        fixedcost
    } else {
        dyncost
    }
}
