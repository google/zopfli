//! C-compatible FFI bindings for Zopfli.

#![allow(non_camel_case_types)]

use std::cmp;
use std::ffi::c_int;
use std::mem::ManuallyDrop;
use std::panic;
use std::process;
use std::ptr;
use std::slice;

use crate::deflate::{deflate, deflate_part};
use crate::gzip_container::gzip_compress;
use crate::lz77::Lz77StoreView;
use crate::squeeze;
use crate::util::{BlockType, SafeOptions, NUM_D, NUM_LL};
use crate::zlib_container::zlib_compress;

/// C-compatible representation of Zopfli options.
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct CZopfliOptions {
    /// Verbose option.
    pub verbose: c_int,
    /// Verbose more option.
    pub verbose_more: c_int,
    /// Number of iterations.
    pub numiterations: c_int,
    /// Enable block splitting.
    pub blocksplitting: c_int,
    /// Left for backwards compatibility.
    pub blocksplittinglast: c_int,
    /// Maximum blocks.
    pub blocksplittingmax: c_int,
}

impl CZopfliOptions {
    /// Safe converter from C struct to safe Rust Options structure.
    pub fn to_safe(&self) -> SafeOptions {
        SafeOptions {
            verbose: self.verbose != 0,
            verbose_more: self.verbose_more != 0,
            numiterations: cmp::max(0, self.numiterations),
            blocksplitting: self.blocksplitting != 0,
            blocksplittinglast: self.blocksplittinglast != 0,
            blocksplittingmax: cmp::max(0, self.blocksplittingmax),
        }
    }
}

impl Default for CZopfliOptions {
    fn default() -> Self {
        Self {
            verbose: 0,
            verbose_more: 0,
            numiterations: 15,
            blocksplitting: 1,
            blocksplittinglast: 0,
            blocksplittingmax: 15,
        }
    }
}

/// C-compatible representation of output formats.
#[repr(i32)]
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum CZopfliFormat {
    /// Gzip output format.
    ZOPFLI_FORMAT_GZIP = 0,
    /// Zlib output format.
    ZOPFLI_FORMAT_ZLIB = 1,
    /// Deflate output format.
    ZOPFLI_FORMAT_DEFLATE = 2,
}

/// C-compatible representation of `ZopfliLZ77Store`.
#[repr(C)]
#[derive(Debug, Copy, Clone)]
pub struct CZopfliLZ77Store {
    /// Literal symbols or length values.
    pub litlens: *mut u16,
    /// Distances (0 indicates literal).
    pub dists: *mut u16,
    /// Size of `litlens` and `dists` arrays.
    pub size: libc::size_t,
    /// Reference to the original data.
    pub data: *const u8,
    /// Position in data where each LZ77 command begins.
    pub pos: *mut libc::c_uint,
    /// Cumulative counts of literal/length symbols.
    pub ll_counts: *mut libc::c_uint,
    /// Cumulative counts of distance symbols.
    pub d_counts: *mut libc::c_uint,
}

impl CZopfliLZ77Store {
    /// Borrows the C-compatible `CZopfliLZ77Store` as a safe `Lz77StoreView` without copying.
    ///
    /// # Safety
    ///
    /// - `self.litlens`, `self.dists`, and `self.pos` must be valid for reads of `self.size`
    ///   elements.
    /// - If non-null, `self.ll_counts` and `self.d_counts` must be valid for reads of their
    ///   respective cumulative block sizes (`NUM_LL` / `NUM_D` rounded up).
    pub unsafe fn as_view<'a>(&self) -> Lz77StoreView<'a> {
        if self.size == 0 {
            return Lz77StoreView::empty();
        }

        // SAFETY: Pointer validity and aliasing requirements are guaranteed by caller contract.
        let (litlens, dists, pos) = unsafe {
            (
                slice_from_c_raw_parts(self.litlens, self.size),
                slice_from_c_raw_parts(self.dists, self.size),
                slice_from_c_raw_parts(self.pos, self.size),
            )
        };

        let ll_num = NUM_LL;
        let ll_size = ll_num * self.size.div_ceil(ll_num);
        let d_num = NUM_D;
        let d_size = d_num * self.size.div_ceil(d_num);

        let ll_counts = if !self.ll_counts.is_null() && ll_size > 0 {
            // SAFETY: Guaranteed by caller contract for non-null `self.ll_counts`.
            unsafe { slice_from_c_raw_parts(self.ll_counts, ll_size) }
        } else {
            &[]
        };

        let d_counts = if !self.d_counts.is_null() && d_size > 0 {
            // SAFETY: Guaranteed by caller contract for non-null `self.d_counts`.
            unsafe { slice_from_c_raw_parts(self.d_counts, d_size) }
        } else {
            &[]
        };

        Lz77StoreView {
            litlens,
            dists,
            pos,
            data: &[],
            ll_counts,
            d_counts,
        }
    }
}

/// A growable buffer wrapper managing C-allocated (`libc::malloc`/`libc::realloc`) memory.
///
/// Note: Does not implement `Drop` to preserve caller-owned memory on early exit or panic;
/// ownership must be explicitly returned to C callers via `into_raw`.
struct CBuffer {
    ptr: *mut u8,
    len: usize,
    capacity: usize,
}

impl CBuffer {
    /// Takes ownership of an existing C-allocated buffer and initialized length.
    ///
    /// # Safety
    ///
    /// If `ptr` is non-null, it must have been allocated via `libc::malloc`/`libc::realloc`
    /// and be valid for reads and writes of `len` bytes.
    unsafe fn from_raw(ptr: *mut u8, len: usize) -> Self {
        // Conservative initial capacity: assumes only `len` bytes are allocated.
        // Any append operation will immediately reallocate to a known power-of-two capacity.
        Self { ptr, len, capacity: len }
    }

    /// Appends slice bytes, reallocating via `libc::realloc` when capacity is exceeded.
    fn extend_from_slice(&mut self, data: &[u8]) {
        if data.is_empty() {
            return;
        }
        let new_len = self.len.checked_add(data.len()).unwrap_or_else(|| process::abort());
        if new_len > self.capacity {
            let new_cap = new_len.checked_next_power_of_two().unwrap_or_else(|| process::abort());
            // SAFETY: `self.ptr` is either null or was allocated via malloc/realloc,
            // and `new_cap > 0`.
            let new_ptr = unsafe {
                libc::realloc(self.ptr as *mut libc::c_void, new_cap) as *mut u8
            };
            if new_ptr.is_null() {
                process::abort();
            }
            self.ptr = new_ptr;
            self.capacity = new_cap;
        }
        // SAFETY: `self.ptr` has capacity for `new_len`, source and dest do not overlap.
        unsafe {
            ptr::copy_nonoverlapping(data.as_ptr(), self.ptr.add(self.len), data.len());
        }
        self.len = new_len;
    }

    /// Returns a shared slice of the initialized bytes.
    fn as_slice(&self) -> &[u8] {
        if self.len == 0 || self.ptr.is_null() {
            &[]
        } else {
            // SAFETY: `self.ptr` is valid for reads of `self.len` bytes.
            unsafe { slice::from_raw_parts(self.ptr, self.len) }
        }
    }

    /// Returns a mutable slice of the initialized bytes.
    fn as_mut_slice(&mut self) -> &mut [u8] {
        if self.len == 0 || self.ptr.is_null() {
            &mut []
        } else {
            // SAFETY: `self.ptr` is valid for writes of `self.len` bytes, and `&mut self`
            // is exclusive.
            unsafe { slice::from_raw_parts_mut(self.ptr, self.len) }
        }
    }

    /// Length in bytes.
    fn len(&self) -> usize {
        self.len
    }

    /// Relinquishes ownership back to the C caller pointers.
    fn into_raw(self, out: &mut *mut u8, outsize: &mut libc::size_t) {
        let this = ManuallyDrop::new(self);
        *out = this.ptr;
        *outsize = this.len;
    }
}

/// Initializes Zopfli options with standard default values.
///
/// For C callers:
/// - If `options` is non-null, it must point to a valid, aligned writeable memory region for
///   `CZopfliOptions`.
#[unsafe(no_mangle)]
pub extern "C" fn Rust_ZopfliInitOptions(options: Option<&mut CZopfliOptions>) {
    if let Some(options) = options {
        *options = CZopfliOptions::default();
    }
}

/// Initializes a `CZopfliLZ77Store`.
///
/// For C callers:
/// - If `store` is non-null, it must point to a valid, aligned writeable memory region for
///   `CZopfliLZ77Store`.
#[unsafe(no_mangle)]
pub extern "C" fn Rust_ZopfliInitLZ77Store(
    data: *const u8,
    store: Option<&mut CZopfliLZ77Store>,
) {
    if let Some(store) = store {
        *store = CZopfliLZ77Store {
            litlens: ptr::null_mut(),
            dists: ptr::null_mut(),
            size: 0,
            data,
            pos: ptr::null_mut(),
            ll_counts: ptr::null_mut(),
            d_counts: ptr::null_mut(),
        };
    }
}

/// Frees a pointer allocated via `libc::malloc`/`libc::realloc` and sets it to null.
///
/// # Safety
///
/// If `*ptr` is non-null, it must have been allocated via `libc::malloc`/`realloc`
/// and not yet freed.
unsafe fn free_c_ptr<T>(ptr: &mut *mut T) {
    if !(*ptr).is_null() {
        // SAFETY: `*ptr` was allocated via `libc::malloc`/`realloc` per caller contract.
        unsafe {
            libc::free(*ptr as *mut libc::c_void);
        }
        *ptr = ptr::null_mut();
    }
}

/// Frees internal dynamically allocated memory in a `CZopfliLZ77Store`.
///
/// # Safety
///
/// - `store` must point to a valid `CZopfliLZ77Store` whose internal buffer pointers were
///   allocated via `libc::malloc`/`libc::realloc` (or are null).
/// - The pointers must not have been previously freed.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn Rust_ZopfliCleanLZ77Store(store: Option<&mut CZopfliLZ77Store>) {
    if let Some(store) = store {
        // SAFETY: Pointers were allocated via malloc/realloc per caller contract.
        unsafe {
            free_c_ptr(&mut store.litlens);
            free_c_ptr(&mut store.dists);
            free_c_ptr(&mut store.pos);
            free_c_ptr(&mut store.ll_counts);
            free_c_ptr(&mut store.d_counts);
        }
        store.size = 0;
    }
}

/// Helper to handle safe conversion, raw slice access, and output allocation for compression.
///
/// # Safety
///
/// - `in_data` must be valid for reads of `insize` bytes.
/// - If `out` is `Some`, `*out` must be null or point to a buffer of at least `*outsize` bytes
///   allocated via `libc::malloc`/`libc::realloc`.
/// - If `outsize` is `Some`, `*outsize` must match the byte length of `*out` (or be 0 if `*out` is
///   null).
#[inline(always)]
unsafe fn compress_helper<F>(
    options: Option<&CZopfliOptions>,
    in_data: *const u8,
    insize: libc::size_t,
    out: Option<&mut *mut u8>,
    outsize: Option<&mut libc::size_t>,
    compress_fn: F,
) where
    F: FnOnce(&SafeOptions, &[u8], &mut Vec<u8>) -> Result<(), crate::error::Error>
        + std::panic::UnwindSafe,
{
    let (Some(options), Some(out), Some(outsize)) = (options, out, outsize) else {
        return;
    };
    if insize > 0 && in_data.is_null() {
        return;
    }
    if (*out).is_null() && *outsize > 0 {
        return;
    }

    // SAFETY: Input slice validity guaranteed by caller safety contract.
    let data_slice = unsafe { slice_from_c_raw_parts(in_data, insize) };
    // SAFETY: Output buffer validity guaranteed by caller safety contract.
    let mut c_buffer = unsafe { CBuffer::from_raw(*out, *outsize) };

    let caught = panic::catch_unwind(panic::AssertUnwindSafe(|| -> Result<Vec<u8>, ()> {
        let safe_opts = options.to_safe();
        let mut out_buffer = Vec::new();
        compress_fn(&safe_opts, data_slice, &mut out_buffer).map_err(|_| ())?;
        Ok(out_buffer)
    }));

    if let Ok(Ok(out_buffer)) = caught {
        c_buffer.extend_from_slice(&out_buffer);
        c_buffer.into_raw(out, outsize);
    }
}

/// Compresses per deflate specification and appends the compressed result to the output.
///
/// # Safety
///
/// - `in_data` must be valid for reads of `insize` bytes.
/// - If `out` is `Some`, `*out` must be null or point to a buffer of at least `*outsize` bytes
///   allocated via `libc::malloc`/`libc::realloc`.
/// - If `outsize` is `Some`, `*outsize` must match the byte length of `*out` (or be 0 if `*out` is
///   null).
#[unsafe(no_mangle)]
pub unsafe extern "C" fn Rust_ZopfliCompress(
    options: Option<&CZopfliOptions>,
    output_type: libc::c_int,
    in_data: *const u8,
    insize: libc::size_t,
    out: Option<&mut *mut u8>,
    outsize: Option<&mut libc::size_t>,
) {
    let compress_fn = |opts: &SafeOptions, data: &[u8], buffer: &mut Vec<u8>| {
        match output_type {
            0 => gzip_compress(opts, data, buffer),
            1 => zlib_compress(opts, data, buffer),
            2 => {
                let mut bp = 0u8;
                deflate(
                    opts,
                    BlockType::DynamicTree,
                    /* final_block= */ true,
                    data,
                    &mut bp,
                    buffer,
                )
            }
            _ => std::process::abort(),
        }
    };

    // SAFETY: Preconditions forwarded directly to `compress_helper`.
    unsafe {
        compress_helper(options, in_data, insize, out, outsize, compress_fn);
    }
}

/// Compresses the input data into the gzip format and appends it to the output.
///
/// # Safety
///
/// - `in_data` must be valid for reads of `insize` bytes.
/// - If `out` is `Some`, `*out` must be null or point to a buffer of at least `*outsize` bytes
///   allocated via `libc::malloc`/`libc::realloc`.
/// - If `outsize` is `Some`, `*outsize` must match the byte length of `*out` (or be 0 if `*out` is
///   null).
#[unsafe(no_mangle)]
pub unsafe extern "C" fn Rust_ZopfliGzipCompress(
    options: Option<&CZopfliOptions>,
    in_data: *const u8,
    insize: libc::size_t,
    out: Option<&mut *mut u8>,
    outsize: Option<&mut libc::size_t>,
) {
    // SAFETY: Preconditions forwarded directly to `compress_helper`.
    unsafe {
        compress_helper(options, in_data, insize, out, outsize, gzip_compress);
    }
}

/// Compresses the input data into the zlib format and appends it to the output.
///
/// # Safety
///
/// - `in_data` must be valid for reads of `insize` bytes.
/// - If `out` is `Some`, `*out` must be null or point to a buffer of at least `*outsize` bytes
///   allocated via `libc::malloc`/`libc::realloc`.
/// - If `outsize` is `Some`, `*outsize` must match the byte length of `*out` (or be 0 if `*out` is
///   null).
#[unsafe(no_mangle)]
pub unsafe extern "C" fn Rust_ZopfliZlibCompress(
    options: Option<&CZopfliOptions>,
    in_data: *const u8,
    insize: libc::size_t,
    out: Option<&mut *mut u8>,
    outsize: Option<&mut libc::size_t>,
) {
    // SAFETY: Preconditions forwarded directly to `compress_helper`.
    unsafe {
        compress_helper(options, in_data, insize, out, outsize, zlib_compress);
    }
}

/// Helper to handle bit-level deflate appending and output buffer reallocation across FFI.
///
/// # Safety
///
/// - `in_data` must be valid for reads of `in_len` bytes.
/// - If `out` is `Some`, `*out` must be null or point to a buffer of at least `*outsize` bytes
///   allocated via `libc::malloc`/`libc::realloc`.
/// - If `outsize` is `Some`, `*outsize` must match the byte length of `*out` (or be 0 if `*out` is
///   null).
/// - If `bp` is `Some`, `*bp` must contain a bit offset in `0..=7`.
#[allow(clippy::too_many_arguments)]
#[inline(always)]
unsafe fn deflate_helper<F>(
    options: Option<&CZopfliOptions>,
    btype: libc::c_int,
    in_data: *const u8,
    in_len: libc::size_t,
    bp: Option<&mut u8>,
    out: Option<&mut *mut u8>,
    outsize: Option<&mut libc::size_t>,
    deflate_fn: F,
) where
    F: FnOnce(
            &SafeOptions,
            BlockType,
            &[u8],
            &mut u8,
            &mut Vec<u8>,
        ) -> Result<(), crate::error::Error>
        + std::panic::UnwindSafe,
{
    let (Some(options), Some(bp), Some(out), Some(outsize)) = (options, bp, out, outsize) else {
        return;
    };
    if in_len > 0 && in_data.is_null() {
        return;
    }
    if (*out).is_null() && *outsize > 0 {
        return;
    }

    let block_type = match btype {
        0 => BlockType::Uncompressed,
        1 => BlockType::FixedTree,
        2 => BlockType::DynamicTree,
        _ => return,
    };

    // SAFETY: Input slice validity guaranteed by caller safety contract.
    let in_slice = unsafe { slice_from_c_raw_parts(in_data, in_len) };
    // SAFETY: Output buffer validity guaranteed by caller safety contract.
    let mut c_buffer = unsafe { CBuffer::from_raw(*out, *outsize) };
    let old_bp = *bp & 7;

    let caught = panic::catch_unwind(panic::AssertUnwindSafe(|| -> Result<(Vec<u8>, u8), ()> {
        let safe_opts = options.to_safe();
        let mut local_bp = old_bp;
        let mut out_buffer = Vec::new();

        // If the output already has data and bp != 0, initialize the working buffer
        // with the existing last byte so that bits are appended into it.
        if c_buffer.len() > 0 && old_bp != 0 {
            let last_byte = c_buffer.as_slice()[c_buffer.len() - 1];
            out_buffer.push(last_byte);
        }

        let res = deflate_fn(&safe_opts, block_type, in_slice, &mut local_bp, &mut out_buffer);
        if res.is_err() {
            return Err(());
        }

        Ok((out_buffer, local_bp))
    }));

    if let Ok(Ok((out_buffer, new_bp))) = caught {
        if c_buffer.len() > 0 && old_bp != 0 {
            if !out_buffer.is_empty() {
                let last_idx = c_buffer.len() - 1;
                c_buffer.as_mut_slice()[last_idx] = out_buffer[0];
                if out_buffer.len() > 1 {
                    c_buffer.extend_from_slice(&out_buffer[1..]);
                }
            }
        } else if !out_buffer.is_empty() {
            c_buffer.extend_from_slice(&out_buffer);
        }
        *bp = new_bp;
        c_buffer.into_raw(out, outsize);
    }
}

/// Deflates data per the deflate specification and appends the result to the output.
///
/// # Safety
///
/// - `in_data` must be valid for reads of `insize` bytes.
/// - If `out` is `Some`, `*out` must be null or point to a buffer of at least `*outsize` bytes
///   allocated via `libc::malloc`/`libc::realloc`.
/// - If `outsize` is `Some`, `*outsize` must match the byte size of `*out` (or be 0 if `*out` is
///   null).
/// - If `bp` is `Some`, `*bp` must contain a bit offset in `0..=7`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn Rust_ZopfliDeflate(
    options: Option<&CZopfliOptions>,
    btype: libc::c_int,
    final_block: libc::c_int,
    in_data: *const u8,
    insize: libc::size_t,
    bp: Option<&mut u8>,
    out: Option<&mut *mut u8>,
    outsize: Option<&mut libc::size_t>,
) {
    let deflate_fn = |opts: &SafeOptions,
                      bt: BlockType,
                      data: &[u8],
                      local_bp: &mut u8,
                      buffer: &mut Vec<u8>| {
        deflate(opts, bt, final_block != 0, data, local_bp, buffer)
    };

    // SAFETY: Preconditions forwarded to `deflate_helper`.
    unsafe {
        deflate_helper(options, btype, in_data, insize, bp, out, outsize, deflate_fn);
    }
}

/// Deflates a part of the input data, allowing specification of start and end byte with `instart`
/// and `inend`.
///
/// # Safety
///
/// - `in_data` must be valid for reads of `inend` bytes.
/// - `instart` must be `<= inend`.
/// - If `out` is `Some`, `*out` must be null or point to a buffer of at least `*outsize` bytes
///   allocated via `libc::malloc`/`libc::realloc`.
/// - If `outsize` is `Some`, `*outsize` must match the byte size of `*out` (or be 0 if `*out` is
///   null).
/// - If `bp` is `Some`, `*bp` must contain a bit offset in `0..=7`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn Rust_ZopfliDeflatePart(
    options: Option<&CZopfliOptions>,
    btype: libc::c_int,
    final_block: libc::c_int,
    in_data: *const u8,
    instart: libc::size_t,
    inend: libc::size_t,
    bp: Option<&mut u8>,
    out: Option<&mut *mut u8>,
    outsize: Option<&mut libc::size_t>,
) {
    if instart > inend {
        return;
    }

    let deflate_fn = |opts: &SafeOptions,
                      bt: BlockType,
                      data: &[u8],
                      local_bp: &mut u8,
                      buffer: &mut Vec<u8>| {
        deflate_part(opts, bt, final_block != 0, data, instart, inend, local_bp, buffer)
    };

    // SAFETY: Preconditions forwarded to `deflate_helper`.
    unsafe {
        deflate_helper(options, btype, in_data, inend, bp, out, outsize, deflate_fn);
    }
}

/// Calculates block size in bits.
///
/// # Safety
///
/// - If `lz77` is `Some`, its internal buffer pointers must satisfy the safety invariants
///   required by `CZopfliLZ77Store::as_view` (valid for reads of `lz77.size` elements).
/// - `lstart` and `lend` must satisfy `lstart <= lend <= lz77.size`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn Rust_ZopfliCalculateBlockSize(
    lz77: Option<&CZopfliLZ77Store>,
    lstart: libc::size_t,
    lend: libc::size_t,
    btype: libc::c_int,
) -> libc::c_double {
    let Some(lz77) = lz77 else {
        return 0.0;
    };
    if lstart > lend || lend > lz77.size {
        return 0.0;
    }
    let block_type = match btype {
        0..=2 => btype,
        _ => return 0.0,
    };
    let caught = panic::catch_unwind(panic::AssertUnwindSafe(|| -> f64 {
        // SAFETY: `lz77` validity is guaranteed by the caller safety contract.
        let view = unsafe { lz77.as_view() };
        squeeze::calculate_block_size(&view, lstart, lend, block_type)
    }));
    caught.unwrap_or(0.0)
}

/// Calculates block size in bits, automatically using the best btype.
///
/// # Safety
///
/// - If `lz77` is `Some`, its internal buffer pointers must satisfy the safety invariants
///   required by `CZopfliLZ77Store::as_view` (valid for reads of `lz77.size` elements).
/// - `lstart` and `lend` must satisfy `lstart <= lend <= lz77.size`.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn Rust_ZopfliCalculateBlockSizeAutoType(
    lz77: Option<&CZopfliLZ77Store>,
    lstart: libc::size_t,
    lend: libc::size_t,
) -> libc::c_double {
    let Some(lz77) = lz77 else {
        return 0.0;
    };
    if lstart > lend || lend > lz77.size {
        return 0.0;
    }
    let caught = panic::catch_unwind(panic::AssertUnwindSafe(|| -> f64 {
        // SAFETY: `lz77` validity is guaranteed by the caller safety contract.
        let view = unsafe { lz77.as_view() };
        squeeze::calculate_block_size_auto_type(&view, lstart, lend)
    }));
    caught.unwrap_or(0.0)
}

/// Calculates the CRC32 checksum of the provided data buffer.
///
/// # Safety
///
/// `data` must be valid for reads of `size` bytes.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn Rust_ZopfliCRC32(data: *const u8, size: libc::size_t) -> libc::c_ulong {
    // SAFETY: `data` validity guaranteed by caller safety contract.
    let slice = unsafe { slice_from_c_raw_parts(data, size) };
    if slice.is_empty() {
        return 0;
    }
    crate::checksum::crc32(slice) as libc::c_ulong
}

/// Calculates the Adler32 checksum of the provided data buffer.
///
/// # Safety
///
/// `data` must be valid for reads of `size` bytes.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn Rust_ZopfliAdler32(data: *const u8, size: libc::size_t) -> libc::c_uint {
    // SAFETY: `data` validity guaranteed by caller safety contract.
    let slice = unsafe { slice_from_c_raw_parts(data, size) };
    if slice.is_empty() {
        return 1;
    }
    crate::checksum::adler32(slice) as libc::c_uint
}

/// Frees a buffer previously allocated by Zopfli FFI functions.
///
/// # Safety
///
/// `buffer` must be null or previously allocated by `libc::malloc`/`libc::realloc` and not freed.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn zopfli_free_buffer(buffer: *mut u8) {
    if !buffer.is_null() {
        // SAFETY: `buffer` was allocated via malloc/realloc per caller contract.
        unsafe {
            libc::free(buffer as *mut libc::c_void);
        }
    }
}

/// Helper to safely create a slice from a raw C pointer and size.
/// Handles null pointers and zero sizes by returning an empty slice.
///
/// # Safety
///
/// - `ptr` must be properly aligned and valid for reads of `size` elements of type `T`.
/// - `size * size_of::<T>()` must not exceed `isize::MAX as usize`.
pub unsafe fn slice_from_c_raw_parts<'a, T>(ptr: *const T, size: libc::size_t) -> &'a [T] {
    if size == 0 {
        &[]
    } else {
        assert!(!ptr.is_null());
        // SAFETY: `ptr` is non-null (asserted) and valid for reads of `size` elements per contract.
        // Total size does not wrap isize::MAX as it fits in addressable memory allocated by caller.
        unsafe { slice::from_raw_parts(ptr, size) }
    }
}
