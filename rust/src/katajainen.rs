//! Katajainen length-limited Huffman code length calculation.

#![forbid(unsafe_code)]

use crate::error::Error;
use std::mem;

#[derive(Clone, Copy, Debug, Default)]
struct Leaf {
    weight: usize,
    count: usize, // symbol index
}

#[derive(Clone, Copy, Debug, Default)]
struct Node {
    weight: usize,
    tail: Option<usize>,
    count: i32,
    in_use: bool,
}

#[derive(Debug, Clone, Default)]
struct NodePool {
    nodes: Vec<Node>,
    next: usize,
}

fn get_free_node(
    mut lists: Option<&mut [[Option<usize>; 2]]>,
    maxbits: usize,
    pool: &mut NodePool,
) -> usize {
    loop {
        if pool.next >= pool.nodes.len() {
            // Garbage collection
            for node in &mut pool.nodes {
                node.in_use = false;
            }
            if let Some(ref mut lists_ref) = lists {
                // Ensure all nodes referred in lists are marked in_use
                for i in 0..(maxbits * 2) {
                    let mut curr = lists_ref[i / 2][i % 2];
                    while let Some(idx) = curr {
                        pool.nodes[idx].in_use = true;
                        curr = pool.nodes[idx].tail;
                    }
                }
            }
            pool.next = 0;
        }
        if !pool.nodes[pool.next].in_use {
            break;
        }
        pool.next += 1;
    }
    let res = pool.next;
    pool.nodes[res].in_use = true;
    pool.next += 1;
    res
}

fn init_node(pool: &mut NodePool, idx: usize, weight: usize, count: i32, tail: Option<usize>) {
    pool.nodes[idx].weight = weight;
    pool.nodes[idx].count = count;
    pool.nodes[idx].tail = tail;
    pool.nodes[idx].in_use = true;
}

fn boundary_pm(
    lists: &mut [[Option<usize>; 2]],
    maxbits: usize,
    leaves: &[Leaf],
    numsymbols: usize,
    pool: &mut NodePool,
    index: usize,
    is_final: bool,
) {
    let old_chain_idx = lists[index][1].expect("lists element should be initialized");
    let lastcount = pool.nodes[old_chain_idx].count as usize;

    if index == 0 && lastcount >= numsymbols {
        return;
    }

    let new_chain_idx = get_free_node(Some(lists), maxbits, pool);
    let old_chain_idx = lists[index][1].expect("lists element should be initialized");
    let old_chain_tail = pool.nodes[old_chain_idx].tail;

    lists[index][0] = Some(old_chain_idx);
    lists[index][1] = Some(new_chain_idx);

    if index == 0 {
        init_node(pool, new_chain_idx, leaves[lastcount].weight, (lastcount + 1) as i32, None);
    } else {
        let left_child_idx = lists[index - 1][0].expect("left child should be initialized");
        let right_child_idx = lists[index - 1][1].expect("right child should be initialized");
        let sum = pool.nodes[left_child_idx].weight + pool.nodes[right_child_idx].weight;

        if lastcount < numsymbols && sum > leaves[lastcount].weight {
            init_node(
                pool,
                new_chain_idx,
                leaves[lastcount].weight,
                (lastcount + 1) as i32,
                old_chain_tail,
            );
        } else {
            init_node(pool, new_chain_idx, sum, lastcount as i32, Some(right_child_idx));
            if !is_final {
                boundary_pm(lists, maxbits, leaves, numsymbols, pool, index - 1, false);
                boundary_pm(lists, maxbits, leaves, numsymbols, pool, index - 1, false);
            }
        }
    }
}

fn extract_bit_lengths(
    chain_idx: Option<usize>,
    pool: &NodePool,
    leaves: &[Leaf],
    bitlengths: &mut [u32],
) {
    let mut curr = chain_idx;
    while let Some(idx) = curr {
        let count = pool.nodes[idx].count as usize;
        for i in 0..count {
            bitlengths[leaves[i].count] += 1;
        }
        curr = pool.nodes[idx].tail;
    }
}

/// Outputs minimum-redundancy length-limited code bitlengths for symbols with the
/// given counts. The bitlengths are limited by maxbits.
///
/// Returns proper Result error or Ok(()).
pub fn length_limited_code_lengths(
    frequencies: &[usize],
    maxbits: i32,
    bitlengths: &mut [u32],
) -> Result<(), Error> {
    let n = frequencies.len();
    if bitlengths.len() != n {
        return Err(Error::SliceLengthMismatch);
    }

    for val in bitlengths.iter_mut() {
        *val = 0;
    }

    let mut leaves = Vec::with_capacity(n);
    for (i, &freq) in frequencies.iter().enumerate() {
        if freq > 0 {
            leaves.push(Leaf { weight: freq, count: i });
        }
    }

    let numsymbols = leaves.len();

    if maxbits < 0
        || maxbits > 32
        || (1usize.checked_shl(maxbits as u32).unwrap_or(usize::MAX)) < numsymbols
    {
        return Err(Error::TooFewMaxBits);
    }
    if numsymbols == 0 {
        return Ok(());
    }
    if numsymbols == 1 {
        bitlengths[leaves[0].count] = 1;
        return Ok(());
    }

    for leaf in &mut leaves {
        let limit = 1usize << (mem::size_of::<usize>() * 8 - 9);
        if leaf.weight >= limit {
            return Err(Error::KatajainenCountOverflow);
        }
        leaf.weight = (leaf.weight << 9) | leaf.count;
    }

    leaves.sort_by_key(|l| l.weight);

    for leaf in &mut leaves {
        leaf.weight >>= 9;
    }

    let maxbits = maxbits as usize;
    let pool_size = 2 * maxbits * (maxbits + 1);
    let mut pool = NodePool {
        nodes: vec![Node { weight: 0, tail: None, count: 0, in_use: false }; pool_size],
        next: 0,
    };

    let mut lists = vec![[None, None]; maxbits];

    let node0 = get_free_node(None, maxbits, &mut pool);
    let node1 = get_free_node(None, maxbits, &mut pool);
    init_node(&mut pool, node0, leaves[0].weight, 1, None);
    init_node(&mut pool, node1, leaves[1].weight, 2, None);

    for i in 0..maxbits {
        lists[i][0] = Some(node0);
        lists[i][1] = Some(node1);
    }

    let num_boundary_pm_runs = 2 * numsymbols - 4;
    for i in 0..num_boundary_pm_runs {
        let is_final = i == num_boundary_pm_runs - 1;
        boundary_pm(&mut lists, maxbits, &leaves, numsymbols, &mut pool, maxbits - 1, is_final);
    }

    extract_bit_lengths(lists[maxbits - 1][1], &pool, &leaves, bitlengths);

    Ok(())
}
