//! DEFLATE symbol helper functions.

use crate::util::MIN_LENGTH_SYMBOL;

const DIST_DIRECT_LIMIT: usize = 5;

/// Gets the amount of extra bits for the given distance, according to the DEFLATE spec.
pub const fn get_dist_extra_bits(dist: usize) -> u32 {
    if dist < DIST_DIRECT_LIMIT {
        0
    } else {
        (31 ^ ((dist - 1) as u32).leading_zeros()) - 1
    }
}

/// Gets value of the extra bits for the given distance, according to the DEFLATE spec.
pub const fn get_dist_extra_bits_value(dist: usize) -> u32 {
    if dist < DIST_DIRECT_LIMIT {
        0
    } else {
        let l = 31 ^ ((dist - 1) as u32).leading_zeros() as i32;
        let val = (dist as i32 - (1 + (1 << l))) & ((1 << (l - 1)) - 1);
        val as u32
    }
}

/// Gets the symbol for the given distance, according to the DEFLATE spec.
pub const fn get_dist_symbol(dist: usize) -> u32 {
    if dist < DIST_DIRECT_LIMIT {
        (dist as u32).saturating_sub(1)
    } else {
        let l = 31 ^ ((dist - 1) as u32).leading_zeros() as i32;
        let r = (((dist - 1) as i32) >> (l - 1)) & 1;
        (l * 2 + r) as u32
    }
}

/// Gets the amount of extra bits for the given length, according to the DEFLATE spec.
pub const fn get_length_extra_bits(l: usize) -> u32 {
    const TABLE: [u32; 259] = {
        let mut t = [0; 259];
        // Populate the lookup table exactly matching the original Zopfli table.
        t[11] = 1;
        t[12] = 1;
        t[13] = 1;
        t[14] = 1;
        t[15] = 1;
        t[16] = 1;
        t[17] = 1;
        t[18] = 1;
        let mut i = 19;
        while i < 35 {
            t[i] = 2;
            i += 1;
        }
        while i < 67 {
            t[i] = 3;
            i += 1;
        }
        while i < 131 {
            t[i] = 4;
            i += 1;
        }
        while i < 258 {
            t[i] = 5;
            i += 1;
        }
        t[258] = 0;
        t
    };
    if l >= TABLE.len() {
        0
    } else {
        TABLE[l]
    }
}

/// Gets value of the extra bits for the given length, according to the DEFLATE spec.
pub const fn get_length_extra_bits_value(l: usize) -> u32 {
    const TABLE: [u32; 259] = {
        let mut t = [0; 259];
        // Populate the table manually with identical structure and values as legacy C.
        t[12] = 1;
        t[14] = 1;
        t[16] = 1;
        t[18] = 1;
        t[20] = 1;
        t[21] = 2;
        t[22] = 3;
        t[24] = 1;
        t[25] = 2;
        t[26] = 3;
        t[28] = 1;
        t[29] = 2;
        t[30] = 3;
        t[32] = 1;
        t[33] = 2;
        t[34] = 3;

        let mut i = 36;
        while i < 43 {
            t[i] = (i - 35) as u32;
            i += 1;
        }
        i = 44;
        while i < 51 {
            t[i] = (i - 43) as u32;
            i += 1;
        }
        i = 52;
        while i < 59 {
            t[i] = (i - 51) as u32;
            i += 1;
        }
        i = 60;
        while i < 67 {
            t[i] = (i - 59) as u32;
            i += 1;
        }

        i = 68;
        while i < 83 {
            t[i] = (i - 67) as u32;
            i += 1;
        }
        i = 84;
        while i < 99 {
            t[i] = (i - 83) as u32;
            i += 1;
        }
        i = 100;
        while i < 115 {
            t[i] = (i - 99) as u32;
            i += 1;
        }
        i = 116;
        while i < 131 {
            t[i] = (i - 115) as u32;
            i += 1;
        }

        i = 132;
        while i < 163 {
            t[i] = (i - 131) as u32;
            i += 1;
        }
        i = 164;
        while i < 195 {
            t[i] = (i - 163) as u32;
            i += 1;
        }
        i = 196;
        while i < 227 {
            t[i] = (i - 195) as u32;
            i += 1;
        }
        i = 228;
        while i < 258 {
            t[i] = (i - 227) as u32;
            i += 1;
        }
        t[258] = 0;
        t
    };
    if l >= TABLE.len() {
        0
    } else {
        TABLE[l]
    }
}

/// Gets the symbol for the given length, according to the DEFLATE spec.
/// Returns the symbol in the range [257-285] (inclusive).
pub const fn get_length_symbol(l: usize) -> u32 {
    const TABLE: [u32; 259] = {
        let mut t = [0; 259];
        t[3] = 257;
        t[4] = 258;
        t[5] = 259;
        t[6] = 260;
        t[7] = 261;
        t[8] = 262;
        t[9] = 263;
        t[10] = 264;
        t[11] = 265;
        t[12] = 265;
        t[13] = 266;
        t[14] = 266;
        t[15] = 267;
        t[16] = 267;
        t[17] = 268;
        t[18] = 268;

        let mut i = 19;
        while i < 23 {
            t[i] = 269;
            i += 1;
        }
        while i < 27 {
            t[i] = 270;
            i += 1;
        }
        while i < 31 {
            t[i] = 271;
            i += 1;
        }
        while i < 35 {
            t[i] = 272;
            i += 1;
        }

        while i < 43 {
            t[i] = 273;
            i += 1;
        }
        while i < 51 {
            t[i] = 274;
            i += 1;
        }
        while i < 59 {
            t[i] = 275;
            i += 1;
        }
        while i < 67 {
            t[i] = 276;
            i += 1;
        }

        while i < 83 {
            t[i] = 277;
            i += 1;
        }
        while i < 99 {
            t[i] = 278;
            i += 1;
        }
        while i < 115 {
            t[i] = 279;
            i += 1;
        }
        while i < 131 {
            t[i] = 280;
            i += 1;
        }

        while i < 163 {
            t[i] = 281;
            i += 1;
        }
        while i < 195 {
            t[i] = 282;
            i += 1;
        }
        while i < 227 {
            t[i] = 283;
            i += 1;
        }
        while i < 258 {
            t[i] = 284;
            i += 1;
        }
        t[258] = 285;
        t
    };
    if l >= TABLE.len() {
        0
    } else {
        TABLE[l]
    }
}

/// Gets the amount of extra bits for the given length symbol.
pub fn get_length_symbol_extra_bits(s: usize) -> u32 {
    const TABLE: [u32; 29] =
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0];
    // SAFETY: reproduces C's `table[s - 257]` with no bounds check; an s outside the
    // table reads past it. Memory-unsafe; see the PR note.
    unsafe { *TABLE.get_unchecked(s.wrapping_sub(MIN_LENGTH_SYMBOL)) }
}

/// Gets the amount of extra bits for the given distance symbol.
pub fn get_dist_symbol_extra_bits(s: usize) -> u32 {
    const TABLE: [u32; 30] = [
        0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12,
        13, 13,
    ];
    // SAFETY: reproduces C's `table[s]` with no bounds check; an s past the table reads
    // past it. Memory-unsafe; see the PR note.
    unsafe { *TABLE.get_unchecked(s) }
}
