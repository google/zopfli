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

/*
Functions for basic LZ77 compression and utilities for the "squeeze" LZ77
compression.
*/

#ifndef ZOPFLI_LZ77_H_
#define ZOPFLI_LZ77_H_

#include <stdlib.h>

#include "cache.h"
#include "hash.h"
#include "zopfli.h"

/*
Stores lit/length and dist pairs for LZ77.
litlens: Contains the literal symbols or length values.
dists: Indicates the distance, or 0 to indicate that there is no distance and
litlens contains a literal instead of a length.
litlens and dists both have the same size.
*/
typedef struct ZopfliLZ77Store {
  unsigned short* litlens;  /* Lit or len. */
  unsigned short* dists;  /* If 0: indicates literal in corresponding litlens,
      if > 0: length in corresponding litlens, this is the distance. */
  size_t size;
} ZopfliLZ77Store;

void ZopfliInitLZ77Store(ZopfliLZ77Store* store);
void ZopfliCleanLZ77Store(ZopfliLZ77Store* store);
void ZopfliCopyLZ77Store(const ZopfliLZ77Store* source, ZopfliLZ77Store* dest);
void ZopfliStoreLitLenDist(unsigned short length, unsigned short dist,
                           ZopfliLZ77Store* store);

/*
Some state information for compressing a block.
This is currently a bit under-used (with mainly only the longest match cache),
but is kept for easy future expansion.
*/
typedef struct ZopfliBlockState {
  const ZopfliOptions* options;

#ifdef ZOPFLI_LONGEST_MATCH_CACHE
  /* Cache for length/distance pairs found so far. */
  ZopfliLongestMatchCache* lmc;
#endif

  /* The start (inclusive) and end (not inclusive) of the current block. */
  size_t blockstart;
  size_t blockend;
} ZopfliBlockState;

/*
Finds the longest match (length and corresponding distance) for LZ77
compression.
Even when not using "sublen", it can be more efficient to provide an array,
because only then the caching is used.
array: the data
pos: position in the data to find the match for
size: size of the data
limit: limit length to maximum this value (default should be 258). This allows
    finding a shorter dist for that length (= less extra bits). Must be
    in the range [ZOPFLI_MIN_MATCH, ZOPFLI_MAX_MATCH].
sublen: output array of 259 elements, or null. Has, for each length, the
    smallest distance required to reach this length. Only 256 of its 259 values
    are used, the first 3 are ignored (the shortest length is 3. It is purely
    for convenience that the array is made 3 longer).
*/
void ZopfliFindLongestMatch(
    ZopfliBlockState *s, const ZopfliHash* h, const unsigned char* array,
    size_t pos, size_t size, size_t limit,
    unsigned short* sublen, unsigned short* distance, unsigned short* length);

/*
Verifies if length and dist are indeed valid, only used for assertion.
*/
void ZopfliVerifyLenDist(const unsigned char* data, size_t datasize, size_t pos,
                         unsigned short dist, unsigned short length);

/*
Counts the number of literal, length and distance symbols in the given lz77
arrays.
litlens: lz77 lit/lengths
dists: ll77 distances
start: where to begin counting in litlens and dists
end: where to stop counting in litlens and dists (not inclusive)
ll_count: count of each lit/len symbol, must have size 288 (see deflate
    standard)
d_count: count of each dist symbol, must have size 32 (see deflate standard)
*/
void ZopfliLZ77Counts(const unsigned short* litlens,
                      const unsigned short* dists,
                      size_t start, size_t end,
                      size_t* ll_count, size_t* d_count);

/*
Does LZ77 using an algorithm similar to gzip, with lazy matching, rather than
with the slow but better "squeeze" implementation.
The result is placed in the ZopfliLZ77Store.
If instart is larger than 0, it uses values before instart as starting
dictionary.
*/
void ZopfliLZ77Greedy(ZopfliBlockState* s, const unsigned char* in,
                      size_t instart, size_t inend,
                      ZopfliLZ77Store* store);

#endif  /* ZOPFLI_LZ77_H_ */
