/*
Copyright 2011 Google Inc. All Rights Reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Author: lode.vandevenne@gmail.com (Lode Vandevenne)
Author: jyrki.alakuijala@gmail.com (Jyrki Alakuijala)
*/

#include "deflate.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include "blocksplitter.h"
#include "lz77.h"
#include "squeeze.h"
#include "tree.h"

static void AddBit(int bit,
                   unsigned char* bp, unsigned char** out, size_t* outsize) {
  if (((*bp) & 7) == 0) ZOPFLI_APPEND_DATA(0, out, outsize);
  (*out)[*outsize - 1] |= bit << ((*bp) & 7);
  (*bp)++;
}

static void AddBits(unsigned symbol, unsigned length,
                    unsigned char* bp, unsigned char** out, size_t* outsize) {
  /* TODO(lode): make more efficient (add more bits at once). */
  unsigned i;
  for (i = 0; i < length; i++) {
    unsigned bit = (symbol >> i) & 1;
    if (((*bp) & 7) == 0) ZOPFLI_APPEND_DATA(0, out, outsize);
    (*out)[*outsize - 1] |= bit << ((*bp) & 7);
    (*bp)++;
  }
}

/*
Adds bits, like AddBits, but the order is inverted. The deflate specification
uses both orders in one standard.
*/
static void AddHuffmanBits(unsigned symbol, unsigned length,
                           unsigned char* bp, unsigned char** out,
                           size_t* outsize) {
  /* TODO(lode): make more efficient (add more bits at once). */
  unsigned i;
  for (i = 0; i < length; i++) {
    unsigned bit = (symbol >> (length - i - 1)) & 1;
    if (((*bp) & 7) == 0) ZOPFLI_APPEND_DATA(0, out, outsize);
    (*out)[*outsize - 1] |= bit << ((*bp) & 7);
    (*bp)++;
  }
}

/*
Ensures there are at least 2 distance codes to support buggy decoders.
Zlib 1.2.1 and below have a bug where it fails if there isn't at least 1
distance code (with length > 0), even though it's valid according to the
deflate spec to have 0 distance codes. On top of that, some mobile phones
require at least two distance codes. To support these decoders too (but
potentially at the cost of a few bytes), add dummy code lengths of 1.
References to this bug can be found in the changelog of
Zlib 1.2.2 and here: http://www.jonof.id.au/forum/index.php?topic=515.0.

d_lengths: the 32 lengths of the distance codes.
*/
static void PatchDistanceCodesForBuggyDecoders(unsigned* d_lengths) {
  int num_dist_codes = 0; /* Amount of non-zero distance codes */
  int i;
  for (i = 0; i < 30 /* Ignore the two unused codes from the spec */; i++) {
    if (d_lengths[i]) num_dist_codes++;
    if (num_dist_codes >= 2) return; /* Two or more codes is fine. */
  }

  if (num_dist_codes == 0) {
    d_lengths[0] = d_lengths[1] = 1;
  } else if (num_dist_codes == 1) {
    d_lengths[d_lengths[0] ? 1 : 0] = 1;
  }
}

static void AddDynamicTree(const unsigned* ll_lengths,
                           const unsigned* d_lengths,
                           unsigned char* bp,
                           unsigned char** out, size_t* outsize) {
  unsigned* lld_lengths = 0;  /* All litlen and dist lengthts with ending zeros
      trimmed together in one array. */
  unsigned lld_total;  /* Size of lld_lengths. */
  unsigned* rle = 0;  /* Runlength encoded version of lengths of litlen and dist
      trees. */
  unsigned* rle_bits = 0;  /* Extra bits for rle values 16, 17 and 18. */
  size_t rle_size = 0;  /* Size of rle array. */
  size_t rle_bits_size = 0;  /* Should have same value as rle_size. */
  unsigned hlit = 29; /* 286 - 257 */
  unsigned hdist = 29;  /* 32 - 1, but gzip does not like hdist > 29.*/
  unsigned hclen;
  size_t i, j;
  size_t clcounts[19];
  unsigned clcl[19];  /* Code length code lengths. */
  unsigned clsymbols[19];
  /* The order in which code length code lengths are encoded as per deflate. */
  unsigned order[19] = {
    16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15
  };

  /* Trim zeros. */
  while (hlit > 0 && ll_lengths[257 + hlit - 1] == 0) hlit--;
  while (hdist > 0 && d_lengths[1 + hdist - 1] == 0) hdist--;

  lld_total = hlit + 257 + hdist + 1;
  lld_lengths = (unsigned*)malloc(sizeof(*lld_lengths) * lld_total);
  if (!lld_lengths) exit(-1); /* Allocation failed. */

  for (i = 0; i < lld_total; i++) {
    lld_lengths[i] = i < 257 + hlit
        ? ll_lengths[i] : d_lengths[i - 257 - hlit];
    assert(lld_lengths[i] < 16);
  }

  for (i = 0; i < lld_total; i++) {
    size_t count = 0;
    for (j = i; j < lld_total && lld_lengths[i] == lld_lengths[j]; j++) {
      count++;
    }
    if (count >= 4 || (count >= 3 && lld_lengths[i] == 0)) {
      if (lld_lengths[i] == 0) {
        if (count > 10) {
          if (count > 138) count = 138;
          ZOPFLI_APPEND_DATA(18, &rle, &rle_size);
          ZOPFLI_APPEND_DATA(count - 11, &rle_bits, &rle_bits_size);
        } else {
          ZOPFLI_APPEND_DATA(17, &rle, &rle_size);
          ZOPFLI_APPEND_DATA(count - 3, &rle_bits, &rle_bits_size);
        }
      } else {
        unsigned repeat = count - 1;  /* Since the first one is hardcoded. */
        ZOPFLI_APPEND_DATA(lld_lengths[i], &rle, &rle_size);
        ZOPFLI_APPEND_DATA(0, &rle_bits, &rle_bits_size);
        while (repeat >= 6) {
          ZOPFLI_APPEND_DATA(16, &rle, &rle_size);
          ZOPFLI_APPEND_DATA(6 - 3, &rle_bits, &rle_bits_size);
          repeat -= 6;
        }
        if (repeat >= 3) {
          ZOPFLI_APPEND_DATA(16, &rle, &rle_size);
          ZOPFLI_APPEND_DATA(3 - 3, &rle_bits, &rle_bits_size);
          repeat -= 3;
        }
        while (repeat != 0) {
          ZOPFLI_APPEND_DATA(lld_lengths[i], &rle, &rle_size);
          ZOPFLI_APPEND_DATA(0, &rle_bits, &rle_bits_size);
          repeat--;
        }
      }

      i += count - 1;
    } else {
      ZOPFLI_APPEND_DATA(lld_lengths[i], &rle, &rle_size);
      ZOPFLI_APPEND_DATA(0, &rle_bits, &rle_bits_size);
    }
    assert(rle[rle_size - 1] <= 18);
  }

  for (i = 0; i < 19; i++) {
    clcounts[i] = 0;
  }
  for (i = 0; i < rle_size; i++) {
    clcounts[rle[i]]++;
  }

  ZopfliCalculateBitLengths(clcounts, 19, 7, clcl);
  ZopfliLengthsToSymbols(clcl, 19, 7, clsymbols);

  hclen = 15;
  /* Trim zeros. */
  while (hclen > 0 && clcounts[order[hclen + 4 - 1]] == 0) hclen--;

  AddBits(hlit, 5, bp, out, outsize);
  AddBits(hdist, 5, bp, out, outsize);
  AddBits(hclen, 4, bp, out, outsize);

  for (i = 0; i < hclen + 4; i++) {
    AddBits(clcl[order[i]], 3, bp, out, outsize);
  }

  for (i = 0; i < rle_size; i++) {
    unsigned symbol = clsymbols[rle[i]];
    AddHuffmanBits(symbol, clcl[rle[i]], bp, out, outsize);
    /* Extra bits. */
    if (rle[i] == 16) AddBits(rle_bits[i], 2, bp, out, outsize);
    else if (rle[i] == 17) AddBits(rle_bits[i], 3, bp, out, outsize);
    else if (rle[i] == 18) AddBits(rle_bits[i], 7, bp, out, outsize);
  }

  free(lld_lengths);
  free(rle);
  free(rle_bits);
}

/*
Gives the exact size of the tree, in bits, as it will be encoded in DEFLATE.
*/
static size_t CalculateTreeSize(const unsigned* ll_lengths,
                                const unsigned* d_lengths,
                                size_t* ll_counts, size_t* d_counts) {
  unsigned char* dummy = 0;
  size_t dummysize = 0;
  unsigned char bp = 0;

  (void)ll_counts;
  (void)d_counts;

  AddDynamicTree(ll_lengths, d_lengths, &bp, &dummy, &dummysize);
  free(dummy);

  return dummysize * 8 + (bp & 7);
}

/*
Adds all lit/len and dist codes from the lists as huffman symbols. Does not add
end code 256. expected_data_size is the uncompressed block size, used for
assert, but you can set it to 0 to not do the assertion.
*/
static void AddLZ77Data(const unsigned short* litlens,
                        const unsigned short* dists,
                        size_t lstart, size_t lend,
                        size_t expected_data_size,
                        const unsigned* ll_symbols, const unsigned* ll_lengths,
                        const unsigned* d_symbols, const unsigned* d_lengths,
                        unsigned char* bp,
                        unsigned char** out, size_t* outsize) {
  size_t testlength = 0;
  size_t i;

  for (i = lstart; i < lend; i++) {
    unsigned dist = dists[i];
    unsigned litlen = litlens[i];
    if (dist == 0) {
      assert(litlen < 256);
      assert(ll_lengths[litlen] > 0);
      AddHuffmanBits(ll_symbols[litlen], ll_lengths[litlen], bp, out, outsize);
      testlength++;
    } else {
      unsigned lls = ZopfliGetLengthSymbol(litlen);
      unsigned ds = ZopfliGetDistSymbol(dist);
      assert(litlen >= 3 && litlen <= 288);
      assert(ll_lengths[lls] > 0);
      assert(d_lengths[ds] > 0);
      AddHuffmanBits(ll_symbols[lls], ll_lengths[lls], bp, out, outsize);
      AddBits(ZopfliGetLengthExtraBitsValue(litlen),
              ZopfliGetLengthExtraBits(litlen),
              bp, out, outsize);
      AddHuffmanBits(d_symbols[ds], d_lengths[ds], bp, out, outsize);
      AddBits(ZopfliGetDistExtraBitsValue(dist),
              ZopfliGetDistExtraBits(dist),
              bp, out, outsize);
      testlength += litlen;
    }
  }
  assert(expected_data_size == 0 || testlength == expected_data_size);
}

static void GetFixedTree(unsigned* ll_lengths, unsigned* d_lengths) {
  size_t i;
  for (i = 0; i < 144; i++) ll_lengths[i] = 8;
  for (i = 144; i < 256; i++) ll_lengths[i] = 9;
  for (i = 256; i < 280; i++) ll_lengths[i] = 7;
  for (i = 280; i < 288; i++) ll_lengths[i] = 8;
  for (i = 0; i < 32; i++) d_lengths[i] = 5;
}

/*
Calculates size of the part after the header and tree of an LZ77 block, in bits.
*/
static size_t CalculateBlockSymbolSize(const unsigned* ll_lengths,
                                       const unsigned* d_lengths,
                                       const unsigned short* litlens,
                                       const unsigned short* dists,
                                       size_t lstart, size_t lend) {
  size_t result = 0;
  size_t i;
  for (i = lstart; i < lend; i++) {
    if (dists[i] == 0) {
      result += ll_lengths[litlens[i]];
    } else {
      result += ll_lengths[ZopfliGetLengthSymbol(litlens[i])];
      result += d_lengths[ZopfliGetDistSymbol(dists[i])];
      result += ZopfliGetLengthExtraBits(litlens[i]);
      result += ZopfliGetDistExtraBits(dists[i]);
    }
  }
  result += ll_lengths[256]; /*end symbol*/
  return result;
}

double ZopfliCalculateBlockSize(const unsigned short* litlens,
                                const unsigned short* dists,
                                size_t lstart, size_t lend, int btype) {
  size_t ll_counts[288];
  size_t d_counts[32];

  unsigned ll_lengths[288];
  unsigned d_lengths[32];

  double result = 3; /*bfinal and btype bits*/

  assert(btype == 1 || btype == 2); /* This is not for uncompressed blocks. */

  if(btype == 1) {
    GetFixedTree(ll_lengths, d_lengths);
  } else {
    ZopfliLZ77Counts(litlens, dists, lstart, lend, ll_counts, d_counts);
    ZopfliCalculateBitLengths(ll_counts, 288, 15, ll_lengths);
    ZopfliCalculateBitLengths(d_counts, 32, 15, d_lengths);
    PatchDistanceCodesForBuggyDecoders(d_lengths);
    result += CalculateTreeSize(ll_lengths, d_lengths, ll_counts, d_counts);
  }

  result += CalculateBlockSymbolSize(
      ll_lengths, d_lengths, litlens, dists, lstart, lend);

  return result;
}

/*
Adds a deflate block with the given LZ77 data to the output.
options: global program options
btype: the block type, must be 1 or 2
final: whether to set the "final" bit on this block, must be the last block
litlens: literal/length array of the LZ77 data, in the same format as in
    ZopfliLZ77Store.
dists: distance array of the LZ77 data, in the same format as in
    ZopfliLZ77Store.
lstart: where to start in the LZ77 data
lend: where to end in the LZ77 data (not inclusive)
expected_data_size: the uncompressed block size, used for assert, but you can
  set it to 0 to not do the assertion.
bp: output bit pointer
out: dynamic output array to append to
outsize: dynamic output array size
*/
static void AddLZ77Block(const ZopfliOptions* options, int btype, int final,
                         const unsigned short* litlens,
                         const unsigned short* dists,
                         size_t lstart, size_t lend,
                         size_t expected_data_size,
                         unsigned char* bp, unsigned char** out, size_t* outsize) {
  size_t ll_counts[288];
  size_t d_counts[32];
  unsigned ll_lengths[288];
  unsigned d_lengths[32];
  unsigned ll_symbols[288];
  unsigned d_symbols[32];
  size_t detect_block_size = *outsize;
  size_t compressed_size;
  size_t uncompressed_size = 0;
  size_t i;

  AddBit(final, bp, out, outsize);
  AddBit(btype & 1, bp, out, outsize);
  AddBit((btype & 2) >> 1, bp, out, outsize);

  if (btype == 1) {
    /* Fixed block. */
    GetFixedTree(ll_lengths, d_lengths);
  } else {
    /* Dynamic block. */
    unsigned detect_tree_size;
    assert(btype == 2);
    ZopfliLZ77Counts(litlens, dists, lstart, lend, ll_counts, d_counts);
    ZopfliCalculateBitLengths(ll_counts, 288, 15, ll_lengths);
    ZopfliCalculateBitLengths(d_counts, 32, 15, d_lengths);
    PatchDistanceCodesForBuggyDecoders(d_lengths);
    detect_tree_size = *outsize;
    AddDynamicTree(ll_lengths, d_lengths, bp, out, outsize);
    if (options->verbose) {
      fprintf(stderr, "treesize: %d\n", (int)(*outsize - detect_tree_size));
    }

    /* Assert that for every present symbol, the code length is non-zero. */
    /* TODO(lode): remove this in release version. */
    for (i = 0; i < 288; i++) assert(ll_counts[i] == 0 || ll_lengths[i] > 0);
    for (i = 0; i < 32; i++) assert(d_counts[i] == 0 || d_lengths[i] > 0);
  }

  ZopfliLengthsToSymbols(ll_lengths, 288, 15, ll_symbols);
  ZopfliLengthsToSymbols(d_lengths, 32, 15, d_symbols);

  detect_block_size = *outsize;
  AddLZ77Data(litlens, dists, lstart, lend, expected_data_size,
              ll_symbols, ll_lengths, d_symbols, d_lengths,
              bp, out, outsize);
  /* End symbol. */
  AddHuffmanBits(ll_symbols[256], ll_lengths[256], bp, out, outsize);

  for (i = lstart; i < lend; i++) {
    uncompressed_size += dists[i] == 0 ? 1 : litlens[i];
  }
  compressed_size = *outsize - detect_block_size;
  if (options->verbose) {
    fprintf(stderr, "compressed block size: %d (%dk) (unc: %d)\n",
           (int)compressed_size, (int)(compressed_size / 1024),
           (int)(uncompressed_size));
  }
}

static void DeflateDynamicBlock(const ZopfliOptions* options, int final,
                                const unsigned char* in,
                                size_t instart, size_t inend,
                                unsigned char* bp,
                                unsigned char** out, size_t* outsize) {
  ZopfliBlockState s;
  size_t blocksize = inend - instart;
  ZopfliLZ77Store store;
  int btype = 2;

  ZopfliInitLZ77Store(&store);

  s.options = options;
  s.blockstart = instart;
  s.blockend = inend;
#ifdef ZOPFLI_LONGEST_MATCH_CACHE
  s.lmc = (ZopfliLongestMatchCache*)malloc(sizeof(ZopfliLongestMatchCache));
  ZopfliInitCache(blocksize, s.lmc);
#endif

  ZopfliLZ77Optimal(&s, in, instart, inend, &store);

  /* For small block, encoding with fixed tree can be smaller. For large block,
  don't bother doing this expensive test, dynamic tree will be better.*/
  if (store.size < 1000) {
    double dyncost, fixedcost;
    ZopfliLZ77Store fixedstore;
    ZopfliInitLZ77Store(&fixedstore);
    ZopfliLZ77OptimalFixed(&s, in, instart, inend, &fixedstore);
    dyncost = ZopfliCalculateBlockSize(store.litlens, store.dists,
        0, store.size, 2);
    fixedcost = ZopfliCalculateBlockSize(fixedstore.litlens, fixedstore.dists,
        0, fixedstore.size, 1);
    if (fixedcost < dyncost) {
      btype = 1;
      ZopfliCleanLZ77Store(&store);
      store = fixedstore;
    } else {
      ZopfliCleanLZ77Store(&fixedstore);
    }
  }

  AddLZ77Block(s.options, btype, final,
               store.litlens, store.dists, 0, store.size,
               blocksize, bp, out, outsize);

#ifdef ZOPFLI_LONGEST_MATCH_CACHE
  ZopfliCleanCache(s.lmc);
  free(s.lmc);
#endif
  ZopfliCleanLZ77Store(&store);
}

static void DeflateFixedBlock(const ZopfliOptions* options, int final,
                              const unsigned char* in,
                              size_t instart, size_t inend,
                              unsigned char* bp,
                              unsigned char** out, size_t* outsize) {
  ZopfliBlockState s;
  size_t blocksize = inend - instart;
  ZopfliLZ77Store store;

  ZopfliInitLZ77Store(&store);

  s.options = options;
  s.blockstart = instart;
  s.blockend = inend;
#ifdef ZOPFLI_LONGEST_MATCH_CACHE
  s.lmc = (ZopfliLongestMatchCache*)malloc(sizeof(ZopfliLongestMatchCache));
  ZopfliInitCache(blocksize, s.lmc);
#endif

  ZopfliLZ77OptimalFixed(&s, in, instart, inend, &store);

  AddLZ77Block(s.options, 1, final, store.litlens, store.dists, 0, store.size,
               blocksize, bp, out, outsize);

#ifdef ZOPFLI_LONGEST_MATCH_CACHE
  ZopfliCleanCache(s.lmc);
  free(s.lmc);
#endif
  ZopfliCleanLZ77Store(&store);
}

static void DeflateNonCompressedBlock(const ZopfliOptions* options, int final,
                                      const unsigned char* in, size_t instart,
                                      size_t inend,
                                      unsigned char* bp,
                                      unsigned char** out, size_t* outsize) {
  size_t i;
  size_t blocksize = inend - instart;
  unsigned short nlen = ~blocksize;

  (void)options;
  assert(blocksize < 65536);  /* Non compressed blocks are max this size. */

  AddBit(final, bp, out, outsize);
  /* BTYPE 00 */
  AddBit(0, bp, out, outsize);
  AddBit(0, bp, out, outsize);

  /* Any bits of input up to the next byte boundary are ignored. */
  *bp = 0;

  ZOPFLI_APPEND_DATA(blocksize % 256, out, outsize);
  ZOPFLI_APPEND_DATA((blocksize / 256) % 256, out, outsize);
  ZOPFLI_APPEND_DATA(nlen % 256, out, outsize);
  ZOPFLI_APPEND_DATA((nlen / 256) % 256, out, outsize);

  for (i = instart; i < inend; i++) {
    ZOPFLI_APPEND_DATA(in[i], out, outsize);
  }
}

static void DeflateBlock(const ZopfliOptions* options,
                         int btype, int final,
                         const unsigned char* in, size_t instart, size_t inend,
                         unsigned char* bp,
                         unsigned char** out, size_t* outsize) {
  if (btype == 0) {
    DeflateNonCompressedBlock(
        options, final, in, instart, inend, bp, out, outsize);
  } else if (btype == 1) {
     DeflateFixedBlock(options, final, in, instart, inend, bp, out, outsize);
  } else {
    assert (btype == 2);
    DeflateDynamicBlock(options, final, in, instart, inend, bp, out, outsize);
  }
}

/*
Does squeeze strategy where first block splitting is done, then each block is
squeezed.
Parameters: see description of the ZopfliDeflate function.
*/
static void DeflateSplittingFirst(const ZopfliOptions* options,
                                  int btype, int final,
                                  const unsigned char* in,
                                  size_t instart, size_t inend,
                                  unsigned char* bp,
                                  unsigned char** out, size_t* outsize) {
  size_t i;
  size_t* splitpoints = 0;
  size_t npoints = 0;
  if (btype == 0) {
    ZopfliBlockSplitSimple(in, instart, inend, 65535, &splitpoints, &npoints);
  } else if (btype == 1) {
    /* If all blocks are fixed tree, splitting into separate blocks only
    increases the total size. Leave npoints at 0, this represents 1 block. */
  } else {
    ZopfliBlockSplit(options, in, instart, inend,
                     options->blocksplittingmax, &splitpoints, &npoints);
  }

  for (i = 0; i <= npoints; i++) {
    size_t start = i == 0 ? instart : splitpoints[i - 1];
    size_t end = i == npoints ? inend : splitpoints[i];
    DeflateBlock(options, btype, i == npoints && final, in, start, end,
                 bp, out, outsize);
  }

  free(splitpoints);
}

/*
Does squeeze strategy where first the best possible lz77 is done, and then based
on that data, block splitting is done.
Parameters: see description of the ZopfliDeflate function.
*/
static void DeflateSplittingLast(const ZopfliOptions* options,
                                 int btype, int final,
                                 const unsigned char* in,
                                 size_t instart, size_t inend,
                                 unsigned char* bp,
                                 unsigned char** out, size_t* outsize) {
  size_t i;
  ZopfliBlockState s;
  ZopfliLZ77Store store;
  size_t* splitpoints = 0;
  size_t npoints = 0;

  if (btype == 0) {
    /* This function only supports LZ77 compression. DeflateSplittingFirst
       supports the special case of noncompressed data. Punt it to that one. */
    DeflateSplittingFirst(options, btype, final,
                          in, instart, inend,
                          bp, out, outsize);
  }
  assert(btype == 1 || btype == 2);

  ZopfliInitLZ77Store(&store);

  s.options = options;
  s.blockstart = instart;
  s.blockend = inend;
#ifdef ZOPFLI_LONGEST_MATCH_CACHE
  s.lmc = (ZopfliLongestMatchCache*)malloc(sizeof(ZopfliLongestMatchCache));
  ZopfliInitCache(inend - instart, s.lmc);
#endif

  if (btype == 2) {
    ZopfliLZ77Optimal(&s, in, instart, inend, &store);
  } else {
    assert (btype == 1);
    ZopfliLZ77OptimalFixed(&s, in, instart, inend, &store);
  }

  if (btype == 1) {
    /* If all blocks are fixed tree, splitting into separate blocks only
    increases the total size. Leave npoints at 0, this represents 1 block. */
  } else {
    ZopfliBlockSplitLZ77(options, store.litlens, store.dists, store.size,
                         options->blocksplittingmax, &splitpoints, &npoints);
  }

  for (i = 0; i <= npoints; i++) {
    size_t start = i == 0 ? 0 : splitpoints[i - 1];
    size_t end = i == npoints ? store.size : splitpoints[i];
    AddLZ77Block(options, btype, i == npoints && final,
                 store.litlens, store.dists, start, end, 0,
                 bp, out, outsize);
  }

#ifdef ZOPFLI_LONGEST_MATCH_CACHE
  ZopfliCleanCache(s.lmc);
  free(s.lmc);
#endif

  ZopfliCleanLZ77Store(&store);
}

/*
Deflate a part, to allow ZopfliDeflate() to use multiple master blocks if
needed.
It is possible to call this function multiple times in a row, shifting
instart and inend to next bytes of the data. If instart is larger than 0, then
previous bytes are used as the initial dictionary for LZ77.
This function will usually output multiple deflate blocks. If final is 1, then
the final bit will be set on the last block.
*/
void ZopfliDeflatePart(const ZopfliOptions* options, int btype, int final,
                       const unsigned char* in, size_t instart, size_t inend,
                       unsigned char* bp, unsigned char** out,
                       size_t* outsize) {
  if (options->blocksplitting) {
    if (options->blocksplittinglast) {
      DeflateSplittingLast(options, btype, final, in, instart, inend,
                           bp, out, outsize);
    } else {
      DeflateSplittingFirst(options, btype, final, in, instart, inend,
                            bp, out, outsize);
    }
  } else {
    DeflateBlock(options, btype, final, in, instart, inend, bp, out, outsize);
  }
}

void ZopfliDeflate(const ZopfliOptions* options, int btype, int final,
                   const unsigned char* in, size_t insize,
                   unsigned char* bp, unsigned char** out, size_t* outsize) {
#if ZOPFLI_MASTER_BLOCK_SIZE == 0
  ZopfliDeflatePart(options, btype, final, in, 0, insize, bp, out, outsize);
#else
  size_t i = 0;
  while (i < insize) {
    int masterfinal = (i + ZOPFLI_MASTER_BLOCK_SIZE >= insize);
    int final2 = final && masterfinal;
    size_t size = masterfinal ? insize - i : ZOPFLI_MASTER_BLOCK_SIZE;
    ZopfliDeflatePart(options, btype, final2,
                      in, i, i + size, bp, out, outsize);
    i += size;
  }
#endif
  if (options->verbose) {
    fprintf(stderr,
            "Original Size: %d, Deflate: %d, Compression: %f%% Removed\n",
            (int)insize, (int)*outsize,
            100.0 * (double)(insize - *outsize) / (double)insize);
  }
}
