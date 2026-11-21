//! Huffman tree building and symbol generation.

#![forbid(unsafe_code)]

use crate::error::Error;
use std::cmp::Ordering;
use std::collections::BinaryHeap;

#[derive(Debug, Clone, Copy, Eq, PartialEq)]
struct QueueElement {
    weight: usize,
    index: usize,
    all_index: usize,
}

impl Ord for QueueElement {
    fn cmp(&self, other: &Self) -> Ordering {
        other.weight.cmp(&self.weight).then_with(|| other.index.cmp(&self.index))
    }
}

impl PartialOrd for QueueElement {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[allow(dead_code)]
#[derive(Debug, Clone, Copy, Default)]
struct HuffmanNode {
    parent: Option<usize>,
    index: usize,
    weight: usize,
}

/// Converts a series of Huffman tree bitlengths, to the bit values of the symbols.
pub fn lengths_to_symbols(lengths: &[u32], maxbits: u32, symbols: &mut [u32]) -> Result<(), Error> {
    const MAX_HUFFMAN_BITS: u32 = 32;
    if maxbits > MAX_HUFFMAN_BITS {
        return Err(Error::TooFewMaxBits);
    }
    let n = lengths.len();
    if symbols.len() != n {
        return Err(Error::SliceLengthMismatch);
    }

    for val in symbols.iter_mut() {
        *val = 0;
    }

    let mut bl_count = vec![0; (maxbits + 1) as usize];
    let mut next_code = vec![0; (maxbits + 1) as usize];

    for &len in lengths {
        if len > maxbits {
            return Err(Error::TooFewMaxBits);
        }
        bl_count[len as usize] += 1;
    }

    let mut code = 0;
    bl_count[0] = 0;
    for bits in 1..=maxbits as usize {
        code = (code + bl_count[bits - 1]) << 1;
        next_code[bits] = code;
    }

    for i in 0..n {
        let len = lengths[i] as usize;
        if len != 0 {
            symbols[i] = next_code[len];
            next_code[len] += 1;
        }
    }

    Ok(())
}

/// Calculates the entropy of each symbol, based on the counts of each symbol.
pub fn calculate_entropy(count: &[usize], bitlengths: &mut [f64]) -> Result<(), Error> {
    let n = count.len();
    if bitlengths.len() != n {
        return Err(Error::SliceLengthMismatch);
    }

    let k_inv_log2 = 1.4426950408889;
    const EPSILON: f64 = -1e-5;
    let sum: usize = count.iter().sum();

    let log2sum =
        if sum == 0 { (n as f64).ln() * k_inv_log2 } else { (sum as f64).ln() * k_inv_log2 };

    for i in 0..n {
        if count[i] == 0 {
            bitlengths[i] = log2sum;
        } else {
            bitlengths[i] = log2sum - (count[i] as f64).ln() * k_inv_log2;
        }
        if bitlengths[i] < 0.0 && bitlengths[i] > EPSILON {
            bitlengths[i] = 0.0;
        }
        debug_assert!(bitlengths[i] >= 0.0);
    }

    Ok(())
}

/// Helper function to calculate standard unlimited Huffman tree bitlengths.
pub fn unlimited_code_lengths(histogram: &[usize], bitlengths: &mut [u32]) -> Result<(), Error> {
    let size = histogram.len();
    if bitlengths.len() != size {
        return Err(Error::SliceLengthMismatch);
    }

    if size == 0 {
        return Ok(());
    }
    if size == 1 {
        bitlengths[0] = if histogram[0] > 0 { 1 } else { 0 };
        return Ok(());
    }

    let mut all = Vec::with_capacity(size * 2);
    let mut pq = BinaryHeap::new();
    let mut count = 0;

    for i in 0..size {
        let weight = histogram[i];
        all.push(HuffmanNode { parent: None, index: i, weight });
        if weight > 0 {
            pq.push(QueueElement { weight, index: i, all_index: count });
        }
        count += 1;
    }

    while pq.len() > 1 {
        let left = pq.pop().unwrap();
        let right = pq.pop().unwrap();

        let node_idx = count;
        all.push(HuffmanNode {
            parent: None,
            index: count,
            weight: left.weight.saturating_add(right.weight),
        });
        count += 1;

        all[left.all_index].parent = Some(node_idx);
        all[right.all_index].parent = Some(node_idx);

        pq.push(QueueElement {
            weight: left.weight.saturating_add(right.weight),
            index: node_idx,
            all_index: node_idx,
        });
    }

    if count > size {
        all[count - 1].index = 0;
        for i in (size..=(count - 2)).rev() {
            if let Some(parent_idx) = all[i].parent {
                all[i].index = all[parent_idx].index + 1;
            }
        }
        for i in 0..size {
            if all[i].parent.is_none() {
                bitlengths[i] = 0;
            } else {
                if let Some(parent_idx) = all[i].parent {
                    bitlengths[i] = (all[parent_idx].index + 1) as u32;
                }
            }
        }
    } else {
        for i in 0..size {
            bitlengths[i] = 0;
        }
        if let Some(el) = pq.pop() {
            bitlengths[el.all_index] = 1;
        }
    }

    Ok(())
}

/// Calculates the bitlengths for the Huffman tree, based on the counts of each symbol.
pub fn calculate_bit_lengths(
    count: &[usize],
    maxbits: i32,
    bitlengths: &mut [u32],
) -> Result<(), Error> {
    let n = count.len();
    if bitlengths.len() != n {
        return Err(Error::SliceLengthMismatch);
    }

    crate::katajainen::length_limited_code_lengths(count, maxbits, bitlengths)?;

    Ok(())
}
