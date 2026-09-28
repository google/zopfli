#![forbid(unsafe_code)]

//! DEFLATE compression core.
//!
//! This module implements the DEFLATE compression algorithm (RFC 1951)
//! in 100% pure safe Rust with exact equivalence to the original C-Zopfli logic.

use crate::blocksplitter;
use crate::cache::ZopfliLongestMatchCache;
use crate::error::Error;
use crate::lz77::{Lz77Store, UninitializedLz77Store, ZopfliBlockState};
use crate::squeeze;
use crate::symbols::{
    get_dist_extra_bits, get_dist_extra_bits_value, get_dist_symbol, get_length_extra_bits,
    get_length_extra_bits_value, get_length_symbol,
};
use crate::tree;
use crate::util::{
    BlockType, SafeOptions, MACRO_BLOCK_SIZE, MIN_LENGTH_SYMBOL, NUM_D, NUM_LL, RLE_CODE_COPY,
    RLE_CODE_ZERO_11_138, RLE_CODE_ZERO_3_10, RLE_COPY_EXTRA_BITS, RLE_COPY_MAX, RLE_COPY_MIN,
    RLE_COPY_THRESHOLD, RLE_ZERO_17_EXTRA_BITS, RLE_ZERO_17_MAX, RLE_ZERO_17_MIN,
    RLE_ZERO_18_EXTRA_BITS, RLE_ZERO_18_MAX, RLE_ZERO_18_MIN,
};
use std::cmp;

/// Appends a single bit to the output bitstream.
fn add_bit(bit: u8, bp: &mut u8, out: &mut Vec<u8>) {
    if *bp == 0 || out.is_empty() {
        out.push(0);
    }
    let last_idx = out.len() - 1;
    out[last_idx] |= bit << *bp;
    *bp = (*bp + 1) & 7;
}

/// Appends multiple bits to the output bitstream (LSB first).
fn add_bits(symbol: u32, length: u32, bp: &mut u8, out: &mut Vec<u8>) {
    for i in 0..length {
        let bit = ((symbol >> i) & 1) as u8;
        add_bit(bit, bp, out);
    }
}

/// Appends multiple bits to the output bitstream (MSB first, as used for Huffman codes).
fn add_huffman_bits(symbol: u32, length: u32, bp: &mut u8, out: &mut Vec<u8>) {
    for i in 0..length {
        let bit = ((symbol >> (length - i - 1)) & 1) as u8;
        add_bit(bit, bp, out);
    }
}

/// Encodes the Huffman tree and returns how many bits its encoding takes.
/// If `out` is `None` (size-only mode), we only compute and return the bit length.
fn encode_tree(
    ll_lengths: &[u32],
    d_lengths: &[u32],
    use_16: bool,
    use_17: bool,
    use_18: bool,
    bp: &mut u8,
    out: Option<&mut Vec<u8>>,
) -> usize {
    let mut hlit = 29usize; // 286 - 257
    let mut hdist = 29usize; // 32 - 1, but gzip does not like hdist > 29.

    // Trim zeros
    while hlit > 0 && ll_lengths[MIN_LENGTH_SYMBOL + hlit - 1] == 0 {
        hlit -= 1;
    }
    while hdist > 0 && d_lengths[1 + hdist - 1] == 0 {
        hdist -= 1;
    }
    let hlit2 = hlit + 257;
    let lld_total = hlit2 + hdist + 1;

    let size_only = out.is_none();
    let mut rle = Vec::new();
    let mut rle_bits = Vec::new();
    let mut clcounts = [0usize; 19];

    let mut i = 0;
    while i < lld_total {
        let symbol = if i < hlit2 { ll_lengths[i] as u8 } else { d_lengths[i - hlit2] as u8 };
        let mut count = 1;
        if use_16 || (symbol == 0 && (use_17 || use_18)) {
            let mut j = i + 1;
            while j < lld_total {
                let next_symbol =
                    if j < hlit2 { ll_lengths[j] as u8 } else { d_lengths[j - hlit2] as u8 };
                if symbol == next_symbol {
                    count += 1;
                    j += 1;
                } else {
                    break;
                }
            }
        }
        i += count;

        // Repetitions of zeroes
        if symbol == 0 && count >= RLE_ZERO_17_MIN {
            if use_18 {
                while count >= RLE_ZERO_18_MIN {
                    let count2 = cmp::min(RLE_ZERO_18_MAX, count);
                    if !size_only {
                        rle.push(RLE_CODE_ZERO_11_138 as u32);
                        rle_bits.push((count2 - RLE_ZERO_18_MIN) as u32);
                    }
                    clcounts[RLE_CODE_ZERO_11_138] += 1;
                    count -= count2;
                }
            }
            if use_17 {
                while count >= RLE_ZERO_17_MIN {
                    let count2 = cmp::min(RLE_ZERO_17_MAX, count);
                    if !size_only {
                        rle.push(RLE_CODE_ZERO_3_10 as u32);
                        rle_bits.push((count2 - RLE_ZERO_17_MIN) as u32);
                    }
                    clcounts[RLE_CODE_ZERO_3_10] += 1;
                    count -= count2;
                }
            }
        }

        // Repetitions of any symbol
        if use_16 && count >= RLE_COPY_THRESHOLD {
            count -= 1; // Since the first one is hardcoded.
            clcounts[symbol as usize] += 1;
            if !size_only {
                rle.push(symbol as u32);
                rle_bits.push(0);
            }
            while count >= RLE_COPY_MIN {
                let count2 = cmp::min(RLE_COPY_MAX, count);
                if !size_only {
                    rle.push(RLE_CODE_COPY as u32);
                    rle_bits.push((count2 - RLE_COPY_MIN) as u32);
                }
                clcounts[RLE_CODE_COPY] += 1;
                count -= count2;
            }
        }

        // No or insufficient repetition
        clcounts[symbol as usize] += count;
        while count > 0 {
            if !size_only {
                rle.push(symbol as u32);
                rle_bits.push(0);
            }
            count -= 1;
        }
    }

    let mut clcl = [0u32; 19];
    let mut clsymbols = [0u32; 19];
    let _ = tree::calculate_bit_lengths(&clcounts, 7, &mut clcl);
    if !size_only {
        let _ = tree::lengths_to_symbols(&clcl, 7, &mut clsymbols);
    }

    const ORDER: [usize; 19] = [16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15];

    let mut hclen = 15usize;
    while hclen > 0 && clcounts[ORDER[hclen + 4 - 1]] == 0 {
        hclen -= 1;
    }

    let mut result_size = 14; // hlit, hdist, hclen bits
    result_size += (hclen + 4) * 3; // clcl bits
    for k in 0..19 {
        result_size += clcl[k] as usize * clcounts[k];
    }
    result_size += clcounts[16] * 2;
    result_size += clcounts[17] * 3;
    result_size += clcounts[18] * 7;

    if let Some(out_vec) = out {
        add_bits(hlit as u32, 5, bp, out_vec);
        add_bits(hdist as u32, 5, bp, out_vec);
        add_bits(hclen as u32, 4, bp, out_vec);

        for k in 0..(hclen + 4) {
            add_bits(clcl[ORDER[k]], 3, bp, out_vec);
        }

        for k in 0..rle.len() {
            let rle_val = rle[k] as usize;
            let symbol = clsymbols[rle_val];
            add_huffman_bits(symbol, clcl[rle_val], bp, out_vec);
            if rle_val == RLE_CODE_COPY {
                add_bits(rle_bits[k], RLE_COPY_EXTRA_BITS, bp, out_vec);
            } else if rle_val == RLE_CODE_ZERO_3_10 {
                add_bits(rle_bits[k], RLE_ZERO_17_EXTRA_BITS, bp, out_vec);
            } else if rle_val == RLE_CODE_ZERO_11_138 {
                add_bits(rle_bits[k], RLE_ZERO_18_EXTRA_BITS, bp, out_vec);
            }
        }
    }

    result_size
}

/// Matches AddDynamicTree in C. Finds the best permutation of rle options and writes the tree.
fn add_dynamic_tree(ll_lengths: &[u32], d_lengths: &[u32], bp: &mut u8, out: &mut Vec<u8>) {
    let mut best = 0;
    let mut bestsize = 0usize;

    for i in 0..8 {
        let size = encode_tree(
            ll_lengths,
            d_lengths,
            (i & 1) != 0,
            (i & 2) != 0,
            (i & 4) != 0,
            &mut 0,
            None,
        );
        if bestsize == 0 || size < bestsize {
            bestsize = size;
            best = i;
        }
    }

    let _ = encode_tree(
        ll_lengths,
        d_lengths,
        (best & 1) != 0,
        (best & 2) != 0,
        (best & 4) != 0,
        bp,
        Some(out),
    );
}

/// Writes Huffman-encoded LZ77 codes.
fn add_lz77_data(
    lz77: &Lz77Store,
    lstart: usize,
    lend: usize,
    ll_symbols: &[u32],
    ll_lengths: &[u32],
    d_symbols: &[u32],
    d_lengths: &[u32],
    bp: &mut u8,
    out: &mut Vec<u8>,
) {
    for i in lstart..lend {
        let dist = lz77.dists[i] as usize;
        let litlen = lz77.litlens[i] as usize;
        if dist == 0 {
            add_huffman_bits(ll_symbols[litlen], ll_lengths[litlen], bp, out);
        } else {
            let lls = get_length_symbol(litlen) as usize;
            let ds = get_dist_symbol(dist) as usize;
            add_huffman_bits(ll_symbols[lls], ll_lengths[lls], bp, out);
            add_bits(get_length_extra_bits_value(litlen), get_length_extra_bits(litlen), bp, out);
            add_huffman_bits(d_symbols[ds], d_lengths[ds], bp, out);
            add_bits(get_dist_extra_bits_value(dist), get_dist_extra_bits(dist), bp, out);
        }
    }
}

/// Appends a non-compressed block (BTYPE=00) to the output bitstream.
fn add_non_compressed_block(
    final_block: bool,
    in_data: &[u8],
    instart: usize,
    inend: usize,
    bp: &mut u8,
    out: &mut Vec<u8>,
) {
    let mut pos = instart;
    loop {
        let mut blocksize = 65535;
        if pos + blocksize > inend {
            blocksize = inend - pos;
        }
        let currentfinal = pos + blocksize >= inend;
        let nlen = !blocksize as u16;

        add_bit(if final_block && currentfinal { 1 } else { 0 }, bp, out);
        add_bit(0, bp, out);
        add_bit(0, bp, out);

        // Any bits of input up to the next byte boundary are ignored.
        *bp = 0;

        out.push((blocksize & 0xFF) as u8);
        out.push(((blocksize >> 8) & 0xFF) as u8);
        out.push((nlen & 0xFF) as u8);
        out.push(((nlen >> 8) & 0xFF) as u8);

        for i in 0..blocksize {
            out.push(in_data[pos + i]);
        }

        if currentfinal {
            break;
        }
        pos += blocksize;
    }
}

/// Compresses a single block.
pub fn add_lz77_block(
    _options: &SafeOptions,
    btype: BlockType,
    final_block: bool,
    lz77: &Lz77Store,
    lstart: usize,
    lend: usize,
    bp: &mut u8,
    out: &mut Vec<u8>,
) -> Result<(), Error> {
    if btype == BlockType::Uncompressed {
        let length = squeeze::lz77_get_byte_range(&lz77.as_view(), lstart, lend);
        let pos = if lstart == lend { 0 } else { lz77.pos[lstart] as usize };
        add_non_compressed_block(final_block, lz77.data, pos, pos + length, bp, out);
        return Ok(());
    }

    add_bit(if final_block { 1 } else { 0 }, bp, out);
    let btype_bits = match btype {
        BlockType::FixedTree => 1,
        BlockType::DynamicTree => 2,
        _ => unreachable!(),
    };
    add_bit((btype_bits & 1) as u8, bp, out);
    add_bit(((btype_bits & 2) >> 1) as u8, bp, out);

    let mut ll_lengths = [0u32; NUM_LL];
    let mut d_lengths = [0u32; NUM_D];

    if btype == BlockType::FixedTree {
        squeeze::get_fixed_tree(&mut ll_lengths, &mut d_lengths);
    } else {
        squeeze::get_dynamic_lengths(
            &lz77.as_view(),
            lstart,
            lend,
            &mut ll_lengths,
            &mut d_lengths,
        );
        add_dynamic_tree(&ll_lengths, &d_lengths, bp, out);
    }

    let mut ll_symbols = [0u32; NUM_LL];
    let mut d_symbols = [0u32; NUM_D];
    tree::lengths_to_symbols(&ll_lengths, 15, &mut ll_symbols)?;
    tree::lengths_to_symbols(&d_lengths, 15, &mut d_symbols)?;

    add_lz77_data(lz77, lstart, lend, &ll_symbols, &ll_lengths, &d_symbols, &d_lengths, bp, out);
    add_huffman_bits(ll_symbols[256], ll_lengths[256], bp, out);

    Ok(())
}

/// Chooses the cheapest block type, and calls add_lz77_block.
fn add_lz77_block_auto_type(
    options: &SafeOptions,
    final_block: bool,
    lz77: &Lz77Store,
    lstart: usize,
    lend: usize,
    bp: &mut u8,
    out: &mut Vec<u8>,
) -> Result<(), Error> {
    let uncompressedcost =
        squeeze::calculate_block_size(&lz77.as_view(), lstart, lend, /* btype= */ 0);
    let mut fixedcost =
        squeeze::calculate_block_size(&lz77.as_view(), lstart, lend, /* btype= */ 1);
    let dyncost =
        squeeze::calculate_block_size(&lz77.as_view(), lstart, lend, /* btype= */ 2);

    let expensivefixed = (lz77.litlens.len() < 1000) || fixedcost <= dyncost * 1.1;

    if lstart == lend {
        add_bits(if final_block { 1 } else { 0 }, 1, bp, out);
        add_bits(1, 2, bp, out); // btype 01
        add_bits(0, 7, bp, out); // end symbol 256 code
        return Ok(());
    }

    let mut fixedstore = UninitializedLz77Store::new().initialize(lz77.data);
    if expensivefixed {
        let instart = lz77.pos[lstart] as usize;
        let inend = instart + squeeze::lz77_get_byte_range(&lz77.as_view(), lstart, lend);
        let mut s = ZopfliBlockState {
            options,
            lmc: Some(ZopfliLongestMatchCache::new(inend - instart)),
            blockstart: instart,
            blockend: inend,
        };
        squeeze::lz77_optimal_fixed(&mut s, lz77.data, instart, inend, &mut fixedstore);
        fixedcost =
            squeeze::calculate_block_size(
                &fixedstore.as_view(),
                0,
                fixedstore.litlens.len(),
                1,
            );
    }

    if uncompressedcost < fixedcost && uncompressedcost < dyncost {
        add_lz77_block(options, BlockType::Uncompressed, final_block, lz77, lstart, lend, bp, out)?;
    } else if fixedcost < dyncost {
        if expensivefixed {
            add_lz77_block(
                options,
                BlockType::FixedTree,
                final_block,
                &fixedstore,
                0,
                fixedstore.litlens.len(),
                bp,
                out,
            )?;
        } else {
            add_lz77_block(
                options,
                BlockType::FixedTree,
                final_block,
                lz77,
                lstart,
                lend,
                bp,
                out,
            )?;
        }
    } else {
        add_lz77_block(options, BlockType::DynamicTree, final_block, lz77, lstart, lend, bp, out)?;
    }

    Ok(())
}

/// Matches ZopfliDeflatePart in C. Coordinates block splitting and optimal LZ77 runs.
pub fn deflate_part(
    options: &SafeOptions,
    btype: BlockType,
    final_block: bool,
    in_data: &[u8],
    instart: usize,
    inend: usize,
    bp: &mut u8,
    out: &mut Vec<u8>,
) -> Result<(), Error> {
    if btype == BlockType::Uncompressed {
        add_non_compressed_block(final_block, in_data, instart, inend, bp, out);
        return Ok(());
    } else if btype == BlockType::FixedTree {
        let mut store = UninitializedLz77Store::new().initialize(in_data);
        let mut s = ZopfliBlockState {
            options,
            lmc: Some(ZopfliLongestMatchCache::new(inend - instart)),
            blockstart: instart,
            blockend: inend,
        };
        squeeze::lz77_optimal_fixed(&mut s, in_data, instart, inend, &mut store);
        add_lz77_block(
            options,
            btype,
            final_block,
            &store,
            /* lstart= */ 0,
            store.litlens.len(),
            bp,
            out,
        )?;
        return Ok(());
    }

    // Dynamic block splitting and optimal LZ77 storing
    let mut splitpoints_uncompressed = Vec::new();
    if options.blocksplitting {
        splitpoints_uncompressed = blocksplitter::split(
            options,
            in_data,
            instart,
            inend,
            options.blocksplittingmax as usize,
        );
    }
    let npoints = splitpoints_uncompressed.len();
    let mut splitpoints = vec![0usize; npoints];
    let mut totalcost = 0.0;
    let mut master_lz77 = UninitializedLz77Store::new().initialize(in_data);

    for i in 0..=npoints {
        let start = if i == 0 { instart } else { splitpoints_uncompressed[i - 1] };
        let end = if i == npoints { inend } else { splitpoints_uncompressed[i] };
        let mut s = ZopfliBlockState {
            options,
            lmc: Some(ZopfliLongestMatchCache::new(end - start)),
            blockstart: start,
            blockend: end,
        };
        let mut store = UninitializedLz77Store::new().initialize(in_data);
        squeeze::lz77_optimal(
            &mut s,
            in_data,
            start,
            end,
            options.numiterations as usize,
            &mut store,
        );
        totalcost += squeeze::calculate_block_size_auto_type(
            &store.as_view(),
            0,
            store.litlens.len(),
        );

        // Append to master_lz77
        for k in 0..store.litlens.len() {
            master_lz77.store_lit_len_dist(store.litlens[k], store.dists[k], store.pos[k] as usize);
        }
        if i < npoints {
            splitpoints[i] = master_lz77.litlens.len();
        }
    }

    // Second block splitting attempt on LZ77 store
    if options.blocksplitting && npoints > 1 {
        let splitpoints2 = blocksplitter::split_lz77(
            options,
            &master_lz77,
            options.blocksplittingmax as usize,
        );

        let mut totalcost2 = 0.0;
        let npoints2 = splitpoints2.len();
        for i in 0..=npoints2 {
            let start = if i == 0 { 0 } else { splitpoints2[i - 1] };
            let end = if i == npoints2 { master_lz77.litlens.len() } else { splitpoints2[i] };
            totalcost2 += squeeze::calculate_block_size_auto_type(
                &master_lz77.as_view(),
                start,
                end,
            );
        }

        if totalcost2 < totalcost {
            splitpoints = splitpoints2;
        }
    }

    let final_npoints = splitpoints.len();
    for i in 0..=final_npoints {
        let start = if i == 0 { 0 } else { splitpoints[i - 1] };
        let end = if i == final_npoints { master_lz77.litlens.len() } else { splitpoints[i] };
        let is_final = i == final_npoints && final_block;
        add_lz77_block_auto_type(options, is_final, &master_lz77, start, end, bp, out)?;
    }

    Ok(())
}

/// Splits the input into macroblocks of size 256KB and compresses each block using deflate_part.
pub fn deflate(
    options: &SafeOptions,
    btype: BlockType,
    final_block: bool,
    in_data: &[u8],
    bp: &mut u8,
    out: &mut Vec<u8>,
) -> Result<(), Error> {
    if in_data.is_empty() {
        deflate_part(
            options,
            btype,
            final_block,
            in_data,
            /* instart= */ 0,
            /* inend= */ 0,
            bp,
            out,
        )?;
        return Ok(());
    }

    let mut i = 0usize;
    while i < in_data.len() {
        let macrofinal = i + MACRO_BLOCK_SIZE >= in_data.len();
        let chunk_final = final_block && macrofinal;
        let size = if macrofinal { in_data.len() - i } else { MACRO_BLOCK_SIZE };
        deflate_part(options, btype, chunk_final, in_data, i, i + size, bp, out)?;
        i += size;
    }

    Ok(())
}
