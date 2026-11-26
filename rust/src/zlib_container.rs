//! Zlib container compression.

#![forbid(unsafe_code)]

use crate::deflate::OutBuffer;
use crate::error::Error;
use crate::util::{BlockType, SafeOptions};

/// Compresses the data according to the zlib specification, RFC 1950.
pub fn zlib_compress(
    options: &SafeOptions,
    in_data: &[u8],
    out: &mut impl OutBuffer,
) -> Result<(), Error> {
    let checksum = crate::checksum::adler32(in_data);
    let mut bp = 0u8;

    let cmfflg = 256 * 120 + 3 * 64;
    let fcheck = 31 - (cmfflg % 31);
    let final_header = cmfflg + fcheck;

    out.push((final_header / 256) as u8);
    out.push((final_header % 256) as u8);

    crate::deflate::deflate(
        options,
        BlockType::DynamicTree,
        /* final_block= */ true,
        in_data,
        &mut bp,
        out,
    )?;

    // Append Adler32 checksum in BIG-ENDIAN format
    out.push((checksum >> 24) as u8);
    out.push((checksum >> 16) as u8);
    out.push((checksum >> 8) as u8);
    out.push(checksum as u8);

    if options.verbose {
        // Like C, `outsize` includes output that was in `out` before this call, the subtraction
        // wraps when the output is larger than the input, and an empty input gives "inf".
        let insize = in_data.len();
        let outsize = out.len();
        eprintln!(
            "Original Size: {}, Zlib: {}, Compression: {:.6}% Removed",
            insize as i32,
            outsize as i32,
            100.0 * insize.wrapping_sub(outsize) as f64 / insize as f64
        );
    }

    Ok(())
}
