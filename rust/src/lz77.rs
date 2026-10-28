//! LZ77 structures and algorithms.

use crate::cache::ZopfliLongestMatchCache;
use crate::symbols::{get_dist_symbol, get_length_symbol};
use crate::util::{SafeOptions, MIN_MATCH, NUM_D, NUM_LL};
use std::cmp;
use std::iter;

/// LZ77 store using parallel arrays for efficiency.
#[derive(Debug, Clone)]
pub struct Lz77Store<'a> {
    /// Literal/length symbols.
    pub litlens: Vec<u16>,
    /// Back-reference distances.
    pub dists: Vec<u16>,
    /// Positions in original data where each command starts.
    pub pos: Vec<u32>,
    /// Reference to the original input data.
    pub data: &'a [u8],
    /// Cumulative counts of literal/length symbols.
    pub ll_counts: Vec<u32>,
    /// Cumulative counts of distance symbols.
    pub d_counts: Vec<u32>,
}

impl<'a> Lz77Store<'a> {
    /// Creates a new, empty LZ77 store for the given data.
    fn new(data: &'a [u8]) -> Self {
        Self {
            litlens: Vec::new(),
            dists: Vec::new(),
            pos: Vec::new(),
            data,
            ll_counts: Vec::new(),
            d_counts: Vec::new(),
        }
    }

    /// Borrows this store as an immutable `Lz77StoreView`.
    pub fn as_view(&self) -> Lz77StoreView<'_> {
        Lz77StoreView {
            litlens: &self.litlens,
            dists: &self.dists,
            pos: &self.pos,
            data: self.data,
            ll_counts: &self.ll_counts,
            d_counts: &self.d_counts,
        }
    }

    /// Clears all internal vectors in the store.
    pub fn clear(&mut self) {
        self.litlens.clear();
        self.dists.clear();
        self.pos.clear();
        self.ll_counts.clear();
        self.d_counts.clear();
    }

    /// Appends a literal, length, or distance command to the store.
    pub fn store_lit_len_dist(&mut self, length: u16, dist: u16, pos: usize) {
        let origsize = self.litlens.len();
        let llstart = NUM_LL * (origsize / NUM_LL);
        let dstart = NUM_D * (origsize / NUM_D);

        if origsize % NUM_LL == 0 {
            if origsize == 0 {
                self.ll_counts.extend(iter::repeat(0).take(NUM_LL));
            } else {
                let prev_block_start = origsize - NUM_LL;
                let range = prev_block_start..prev_block_start + NUM_LL;
                self.ll_counts.extend_from_within(range);
            }
        }

        if origsize % NUM_D == 0 {
            if origsize == 0 {
                self.d_counts.extend(iter::repeat(0).take(NUM_D));
            } else {
                let prev_block_start = origsize - NUM_D;
                let range = prev_block_start..prev_block_start + NUM_D;
                self.d_counts.extend_from_within(range);
            }
        }

        self.litlens.push(length);
        self.dists.push(dist);
        self.pos.push(pos as u32);

        if dist == 0 {
            self.ll_counts[llstart + length as usize] += 1;
        } else {
            self.ll_counts[llstart + get_length_symbol(length as usize) as usize] += 1;
            self.d_counts[dstart + get_dist_symbol(dist as usize) as usize] += 1;
        }
    }

    /// Extracts the histogram of symbols in the range `[lstart, lend)`.
    pub fn get_histogram(
        &self,
        lstart: usize,
        lend: usize,
        ll_counts: &mut [usize; NUM_LL],
        d_counts: &mut [usize; NUM_D],
    ) {
        self.as_view().get_histogram(lstart, lend, ll_counts, d_counts);
    }
}

/// Borrowed view of LZ77 store data for calculations without allocations.
#[derive(Debug, Clone, Copy)]
pub struct Lz77StoreView<'a> {
    /// Literal/length symbols.
    pub litlens: &'a [u16],
    /// Back-reference distances.
    pub dists: &'a [u16],
    /// Positions in original data where each command starts.
    pub pos: &'a [u32],
    /// Reference to the original input data.
    pub data: &'a [u8],
    /// Cumulative counts of literal/length symbols.
    pub ll_counts: &'a [u32],
    /// Cumulative counts of distance symbols.
    pub d_counts: &'a [u32],
}

impl<'a> Lz77StoreView<'a> {
    /// Creates an empty LZ77 store view.
    pub fn empty() -> Self {
        Self {
            litlens: &[],
            dists: &[],
            pos: &[],
            data: &[],
            ll_counts: &[],
            d_counts: &[],
        }
    }

    /// Number of LZ77 commands in this view.
    pub fn len(&self) -> usize {
        self.litlens.len()
    }

    /// Returns true if this view has no commands.
    pub fn is_empty(&self) -> bool {
        self.litlens.is_empty()
    }

    /// Extracts the histogram of symbols in the range `[lstart, lend)`.
    pub fn get_histogram(
        &self,
        lstart: usize,
        lend: usize,
        ll_counts: &mut [usize; NUM_LL],
        d_counts: &mut [usize; NUM_D],
    ) {
        if lstart == lend {
            ll_counts.fill(0);
            d_counts.fill(0);
            return;
        }

        if lstart.wrapping_add(NUM_LL * 3) > lend || self.ll_counts.is_empty() || self.d_counts.is_empty() {
            ll_counts.fill(0);
            d_counts.fill(0);
            for (&dist, &litlen) in self.dists[lstart..lend].iter().zip(&self.litlens[lstart..lend])
            {
                if dist == 0 {
                    let idx = litlen as usize;
                    if idx < NUM_LL {
                        ll_counts[idx] += 1;
                    }
                } else {
                    ll_counts[get_length_symbol(litlen as usize) as usize] += 1;
                    d_counts[get_dist_symbol(dist as usize) as usize] += 1;
                }
            }
        } else {
            self.get_histogram_at(lend - 1, ll_counts, d_counts);
            if lstart > 0 {
                let mut ll_counts2 = [0; NUM_LL];
                let mut d_counts2 = [0; NUM_D];
                self.get_histogram_at(lstart - 1, &mut ll_counts2, &mut d_counts2);

                for (count, &count2) in ll_counts.iter_mut().zip(&ll_counts2) {
                    *count = count.saturating_sub(count2);
                }
                for (count, &count2) in d_counts.iter_mut().zip(&d_counts2) {
                    *count = count.saturating_sub(count2);
                }
            }
        }
    }

    fn get_histogram_at(
        &self,
        lpos: usize,
        ll_counts: &mut [usize; NUM_LL],
        d_counts: &mut [usize; NUM_D],
    ) {
        let llpos = NUM_LL * (lpos / NUM_LL);
        let dpos = NUM_D * (lpos / NUM_D);

        for (count, &val) in ll_counts.iter_mut().zip(&self.ll_counts[llpos..llpos + NUM_LL]) {
            *count = val as usize;
        }
        let end_idx = cmp::min(llpos + NUM_LL, self.litlens.len());
        for i in lpos + 1..end_idx {
            let symbol = if self.dists[i] != 0 {
                get_length_symbol(self.litlens[i] as usize) as usize
            } else {
                self.litlens[i] as usize
            };
            if symbol < NUM_LL {
                ll_counts[symbol] = ll_counts[symbol].saturating_sub(1);
            }
        }

        for (count, &val) in d_counts.iter_mut().zip(&self.d_counts[dpos..dpos + NUM_D]) {
            *count = val as usize;
        }
        let end_idx_d = cmp::min(dpos + NUM_D, self.litlens.len());
        for i in lpos + 1..end_idx_d {
            if self.dists[i] != 0 {
                let symbol = get_dist_symbol(self.dists[i] as usize) as usize;
                d_counts[symbol] = d_counts[symbol].saturating_sub(1);
            }
        }
    }
}

/// An uninitialized LZ77 store that does not yet borrow the input data.
///
/// This struct follows the typestate pattern, requiring initialization with input data
/// before it can be used for compression.
#[derive(Debug, Clone, Copy, Default)]
pub struct UninitializedLz77Store;

impl UninitializedLz77Store {
    /// Creates a new uninitialized LZ77 store.
    pub fn new() -> Self {
        UninitializedLz77Store
    }

    /// Initializes the store with the given input data, transitioning it to the `Lz77Store` state.
    pub fn initialize<'a>(self, data: &'a [u8]) -> Lz77Store<'a> {
        Lz77Store::new(data)
    }
}

/// State information for compressing a block.
#[derive(Debug, Clone)]
pub struct ZopfliBlockState<'a> {
    /// Safe options configurations.
    pub options: &'a SafeOptions,
    /// Optional longest match cache.
    pub lmc: Option<ZopfliLongestMatchCache>,
    /// Block start offset in input array.
    pub blockstart: usize,
    /// Block end offset in input array.
    pub blockend: usize,
}

impl<'a> ZopfliBlockState<'a> {
    fn try_get_from_longest_match_cache(
        &self,
        pos: usize,
        limit: &mut usize,
        sublen: Option<&mut [u16]>,
        distance: &mut u16,
        length: &mut u16,
    ) -> bool {
        let lmc = match &self.lmc {
            Some(cache) => cache,
            None => return false,
        };

        let lmcpos = pos - self.blockstart;
        if lmcpos >= lmc.length.len() {
            return false;
        }

        let cache_available = lmc.length[lmcpos] == 0 || lmc.dist[lmcpos] != 0;
        let limit_ok_for_cache = cache_available
            && (*limit == crate::util::MAX_MATCH
                || lmc.length[lmcpos] as usize <= *limit
                || (sublen.is_some()
                    && lmc.max_cached_sublen(lmcpos, lmc.length[lmcpos] as usize) as usize
                        >= *limit));

        if limit_ok_for_cache && cache_available {
            let cached_len = lmc.length[lmcpos] as usize;
            let mut sublen_valid = true;
            if sublen.is_some() {
                if cached_len > lmc.max_cached_sublen(lmcpos, cached_len) as usize {
                    sublen_valid = false;
                }
            }
            if sublen_valid {
                let mut len = cached_len;
                if len > *limit {
                    len = *limit;
                }
                *length = len as u16;
                if let Some(sub) = sublen {
                    lmc.cache_to_sublen(lmcpos, len, sub);
                    *distance = sub[len];
                } else {
                    *distance = lmc.dist[lmcpos];
                }
                return true;
            }
            *limit = cached_len;
        }

        false
    }

    fn store_in_longest_match_cache(
        &mut self,
        pos: usize,
        limit: usize,
        sublen: Option<&[u16]>,
        distance: u16,
        length: u16,
    ) {
        let lmc = match &mut self.lmc {
            Some(cache) => cache,
            None => return,
        };

        let lmcpos = pos - self.blockstart;
        if lmcpos >= lmc.length.len() {
            return;
        }

        let cache_available = lmc.length[lmcpos] == 0 || lmc.dist[lmcpos] != 0;

        if limit == crate::util::MAX_MATCH && sublen.is_some() && !cache_available {
            let d_val = if (length as usize) < crate::util::MIN_MATCH { 0 } else { distance };
            let l_val = if (length as usize) < crate::util::MIN_MATCH { 0 } else { length };
            lmc.dist[lmcpos] = d_val;
            lmc.length[lmcpos] = l_val;
            if let Some(sub) = sublen {
                lmc.sublen_to_cache(sub, lmcpos, length as usize);
            }
        }
    }
}

fn get_match(array: &[u8], mut scan_pos: usize, mut match_pos: usize, end_pos: usize) -> usize {
    assert!(end_pos <= array.len());
    assert!(scan_pos >= match_pos);

    if let Some(safe_end) = end_pos.checked_sub(8) {
        while scan_pos < safe_end {
            // SAFETY: Since `scan_pos < safe_end` and `safe_end == end_pos - 8`, we have
            // `scan_pos + 8 <= end_pos <= array.len()`.
            // Since `scan_pos >= match_pos`, we also have `match_pos + 8 <= array.len()`.
            // Thus, reading 8 bytes from both `scan_pos` and `match_pos` is within bounds.
            let (scan_chunk, match_chunk) = unsafe {
                (
                    u64::from_le_bytes(*(array.as_ptr().add(scan_pos) as *const [u8; 8])),
                    u64::from_le_bytes(*(array.as_ptr().add(match_pos) as *const [u8; 8])),
                )
            };
            if scan_chunk == match_chunk {
                scan_pos += 8;
                match_pos += 8;
            } else {
                break;
            }
        }
    }
    while scan_pos < end_pos {
        // SAFETY: `scan_pos < end_pos <= array.len()`.
        // Also, `scan_pos >= match_pos` guarantees `match_pos < end_pos <= array.len()`.
        // Therefore, both indices are strictly less than `array.len()`.
        let eq = unsafe { *array.get_unchecked(scan_pos) == *array.get_unchecked(match_pos) };
        if eq {
            scan_pos += 1;
            match_pos += 1;
        } else {
            break;
        }
    }
    scan_pos
}

/// Finds the longest match for the given position in the array.
/// Returns (distance, length).
pub fn find_longest_match(
    state: &mut ZopfliBlockState<'_>,
    h: &crate::hash::ZopfliHash,
    array: &[u8],
    pos: usize,
    mut limit: usize,
    mut sublen: Option<&mut [u16]>,
) -> (u16, u16) {
    let mut distance = 0;
    let mut length = 0;

    if state.try_get_from_longest_match_cache(
        pos,
        &mut limit,
        sublen.as_deref_mut(),
        &mut distance,
        &mut length,
    ) {
        return (distance, length);
    }

    if limit > crate::util::MAX_MATCH {
        limit = crate::util::MAX_MATCH;
    }
    assert!(state.blockend <= array.len());

    if pos >= state.blockend {
        return (0, 0);
    }

    if state.blockend - pos < crate::util::MIN_MATCH {
        return (0, 0);
    }

    if pos + limit > state.blockend {
        limit = state.blockend - pos;
    }

    // Static assertions to elide bounds checks

    assert_eq!(h.head.len(), 65536);
    assert_eq!(h.prev.len(), 32768);
    assert_eq!(h.prev2.len(), 32768);
    assert_eq!(h.same.len(), 32768);
    assert_eq!(h.hashval2.len(), 32768);

    let hpos = pos & crate::util::WINDOW_MASK;
    let mut bestdist = 0;
    // Invariant: `bestlength` is bounded by `1 <= bestlength <= limit` (`limit <= MAX_MATCH`).
    // This bounds `bestlength` and prevents any integer overflow when calculating index offsets.
    let mut bestlength = 1;

    let mut use_secondary = false;

    let val_idx = (h.val as usize) & 65535;
    let pp = h.head[val_idx];
    if pp == -1 {
        return (0, 0);
    }
    let pp_idx = (pp as usize) & 32767;
    let mut p = h.prev[pp_idx];

    let p_idx = (p as usize) & 32767;
    let mut dist = if p_idx < pp_idx { pp_idx - p_idx } else { (32768 - p_idx) + pp_idx };

    let mut chain_counter = crate::util::MAX_CHAIN_HITS;

    while dist < crate::util::WINDOW_SIZE {
        let mut currentlength = 0;

        if dist > 0 {
            let match_pos = pos - dist;

            if pos + bestlength >= state.blockend
                || unsafe {
                    // SAFETY: Left side of `||` short-circuits if `pos + bestlength >= blockend`.
                    // Thus, this block is evaluated only when `pos + bestlength < state.blockend`.
                    // Since `state.blockend <= array.len()`, `pos + bestlength < array.len()`.
                    // Since `match_pos = pos - dist` and `dist > 0`, we have:
                    // `match_pos + bestlength < pos + bestlength < array.len()`.
                    // Therefore, both indices are within bounds.
                    *array.get_unchecked(pos + bestlength)
                        == *array.get_unchecked(match_pos + bestlength)
                }
            {
                let mut same_scan = pos;
                let mut same_match = match_pos;
                let same0 = h.same[pos & crate::util::WINDOW_MASK] as usize;
                // SAFETY: `pos < state.blockend` is verified by the early return on entry,
                // and `state.blockend <= array.len()` is asserted on entry. Since neither `pos`
                // nor `state` is mutated, `pos < array.len()` holds.
                // Since `match_pos = pos - dist` and `dist > 0`, `match_pos < pos < array.len()`.
                // Both indices are guaranteed to be within bounds.
                if same0 >= MIN_MATCH
                    && unsafe { *array.get_unchecked(pos) == *array.get_unchecked(match_pos) }
                {
                    let same1 = h.same[(pos - dist) & crate::util::WINDOW_MASK] as usize;
                    let mut same = if same0 < same1 { same0 } else { same1 };
                    if same > limit {
                        same = limit;
                    }
                    same_scan += same;
                    same_match += same;
                }

                let final_scan = get_match(array, same_scan, same_match, pos + limit);
                currentlength = final_scan - pos;
            }

            if currentlength > bestlength {
                if let Some(sub) = sublen.as_deref_mut() {
                    let end_idx = cmp::min(currentlength + 1, sub.len());
                    if bestlength + 1 < end_idx {
                        sub[bestlength + 1..end_idx].fill(dist as u16);
                    }
                }
                bestdist = dist;
                bestlength = currentlength;
                if currentlength >= limit {
                    break;
                }
            }
        }

        let curr_p_idx = (p as usize) & 32767;
        if !use_secondary && bestlength >= h.same[hpos] as usize && h.val2 == h.hashval2[curr_p_idx]
        {
            use_secondary = true;
        }

        let pp_next_idx = curr_p_idx;
        let p_next = if use_secondary { h.prev2[curr_p_idx] } else { h.prev[curr_p_idx] };

        if p_next == p {
            break;
        }
        p = p_next;

        let next_p_idx = (p as usize) & 32767;
        dist += if next_p_idx < pp_next_idx {
            pp_next_idx - next_p_idx
        } else {
            (32768 - next_p_idx) + pp_next_idx
        };

        chain_counter = chain_counter.saturating_sub(1);
        if chain_counter == 0 {
            break;
        }
    }

    state.store_in_longest_match_cache(
        pos,
        limit,
        sublen.as_ref().map(|s| &**s),
        bestdist as u16,
        bestlength as u16,
    );

    (bestdist as u16, bestlength as u16)
}

fn verify_len_dist(data: &[u8], pos: usize, dist: usize, length: usize) {
    debug_assert!(pos + length <= data.len());
    debug_assert_eq!(&data[pos - dist..pos - dist + length], &data[pos..pos + length]);
}

fn get_length_score(length: usize, distance: usize) -> i32 {
    const DISTANCE_PENALTY_THRESHOLD: usize = 1024;
    let score =
        if distance > DISTANCE_PENALTY_THRESHOLD { length.saturating_sub(1) } else { length };
    score as i32
}

/// Greedy and lazy matching LZ77 parsing.
pub fn lz77_greedy(
    state: &mut ZopfliBlockState<'_>,
    in_data: &[u8],
    instart: usize,
    inend: usize,
    store: &mut Lz77Store<'_>,
) {
    state.blockstart = instart;
    state.blockend = inend;

    let mut prev_length = 0;
    let mut prev_match = 0;
    let mut match_available = false;

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

    let mut i = instart;
    while i < inend {
        hash.update(in_data, i, inend);

        let mut temp_sublen = [0u16; 259];

        let res = find_longest_match(
            state,
            &hash,
            in_data,
            i,
            crate::util::MAX_MATCH,
            Some(&mut temp_sublen),
        );
        let mut dist = res.0;
        let mut leng = res.1;

        let lengthscore = get_length_score(leng as usize, dist as usize);
        let prevlengthscore = get_length_score(prev_length as usize, prev_match as usize);

        if match_available {
            match_available = false;
            if lengthscore > prevlengthscore + 1 {
                store.store_lit_len_dist(in_data[i - 1] as u16, 0, i - 1);
                if lengthscore >= crate::util::MIN_MATCH as i32
                    && (leng as usize) < crate::util::MAX_MATCH
                {
                    match_available = true;
                    prev_length = leng;
                    prev_match = dist;
                    i += 1;
                    continue;
                }
            } else {
                leng = prev_length;
                dist = prev_match;
                verify_len_dist(in_data, i - 1, dist as usize, leng as usize);
                store.store_lit_len_dist(leng, dist, i - 1);
                let mut j = 2;
                while j < leng {
                    debug_assert!(i < inend);
                    i += 1;
                    hash.update(in_data, i, inend);
                    j += 1;
                }
                i += 1;
                continue;
            }
        } else if lengthscore >= crate::util::MIN_MATCH as i32
            && (leng as usize) < crate::util::MAX_MATCH
        {
            match_available = true;
            prev_length = leng;
            prev_match = dist;
            i += 1;
            continue;
        }

        if lengthscore >= crate::util::MIN_MATCH as i32 {
            verify_len_dist(in_data, i, dist as usize, leng as usize);
            store.store_lit_len_dist(leng, dist, i);
        } else {
            leng = 1;
            store.store_lit_len_dist(in_data[i] as u16, 0, i);
        }

        let mut j = 1;
        while j < leng {
            debug_assert!(i < inend);
            i += 1;
            hash.update(in_data, i, inend);
            j += 1;
        }

        i += 1;
    }
}
