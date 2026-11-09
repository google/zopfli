//! Zopfli utility functions and option structures.

#![forbid(unsafe_code)]

/// Supported output formats.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    /// Gzip container format (RFC 1952)
    Gzip,
    /// Zlib container format (RFC 1950)
    Zlib,
    /// Raw deflate stream format (RFC 1951)
    Deflate,
}

/// Options for Zopfli compression.
#[derive(Debug, Clone, Copy)]
pub struct SafeOptions {
    /// Whether to print verbose output.
    pub verbose: bool,
    /// Whether to print even more verbose output.
    pub verbose_more: bool,
    /// Number of iterations.
    pub numiterations: i32,
    /// Whether to enable block splitting.
    pub blocksplitting: bool,
    /// Left for backwards compatibility.
    pub blocksplittinglast: bool,
    /// Maximum number of blocks to split into.
    pub blocksplittingmax: i32,
}

impl Default for SafeOptions {
    fn default() -> Self {
        Self {
            verbose: false,
            verbose_more: false,
            numiterations: 15,
            blocksplitting: true,
            blocksplittinglast: false,
            blocksplittingmax: 15,
        }
    }
}

/// Minimum and maximum length that can be encoded in deflate.
pub const MAX_MATCH: usize = 258;
/// Minimum length that can be encoded in deflate.
pub const MIN_MATCH: usize = 3;

/// Minimum length symbol in DEFLATE.
pub const MIN_LENGTH_SYMBOL: usize = 257;

/// Number of distinct literal/length symbols in DEFLATE.
pub const NUM_LL: usize = 288;
/// Number of distinct distance symbols in DEFLATE.
pub const NUM_D: usize = 32;

// RLE constants used by deflate and squeeze.
pub const RLE_CODE_COPY: usize = 16;
pub const RLE_CODE_ZERO_3_10: usize = 17;
pub const RLE_CODE_ZERO_11_138: usize = 18;

pub const RLE_COPY_MIN: usize = 3;
pub const RLE_COPY_MAX: usize = 6;
pub const RLE_COPY_THRESHOLD: usize = 4; // RLE_COPY_MIN + 1
pub const RLE_ZERO_17_MIN: usize = 3;
pub const RLE_ZERO_17_MAX: usize = 10;
pub const RLE_ZERO_18_MIN: usize = 11;
pub const RLE_ZERO_18_MAX: usize = 138;

pub const RLE_COPY_EXTRA_BITS: u32 = 2;
pub const RLE_ZERO_17_EXTRA_BITS: u32 = 3;
pub const RLE_ZERO_18_EXTRA_BITS: u32 = 7;

/// Maximum code length for the code length Huffman tree in DEFLATE.
pub const MAX_CL_BIT_LENGTH: i32 = 7;

/// The window size for deflate. Must be a power of two.
pub const WINDOW_SIZE: usize = 32768;
/// The window mask used to wrap indices into the window.
pub const WINDOW_MASK: usize = WINDOW_SIZE - 1;

/// A block structure size to divide the input into.
pub const MACRO_BLOCK_SIZE: usize = 1000000;

/// Used to initialize costs.
pub const LARGE_FLOAT: f64 = 1e30;

/// For longest match cache length.
pub const CACHE_LENGTH: usize = 8;

/// Limit the max hash chain hits.
pub const MAX_CHAIN_HITS: usize = 8192;

/// Represent the deflate block type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum BlockType {
    /// Non-compressed block.
    Uncompressed = 0,
    /// Block with fixed Huffman tree.
    FixedTree = 1,
    /// Block with dynamic Huffman tree.
    DynamicTree = 2,
}
