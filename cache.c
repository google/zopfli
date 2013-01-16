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

#include "cache.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef USE_LONGEST_MATCH_CACHE

void InitLongestMatchCache(size_t blocksize, LongestMatchCache* lmc) {
  size_t i;
  lmc->length = (unsigned short*)malloc(sizeof(unsigned short) * blocksize);
  lmc->dist = (unsigned short*)malloc(sizeof(unsigned short) * blocksize);
  /* Rather large amount of memory. */
  lmc->sublen = (unsigned char*)malloc(NUM_CACHED_LENGTHS * 3 * blocksize);

  /* length > 0 and dist 0 is invalid combination, which indicates on purpose
  that this cache value is not filled in yet. */
  for (i = 0; i < blocksize; i++) lmc->length[i] = 1;
  for (i = 0; i < blocksize; i++) lmc->dist[i] = 0;
  for (i = 0; i < NUM_CACHED_LENGTHS * blocksize * 3; i++) lmc->sublen[i] = 0;
}

void CleanLongestMatchCache(LongestMatchCache* lmc) {
  free(lmc->length);
  free(lmc->dist);
  free(lmc->sublen);
}

void SublenToCache(const unsigned short* sublen, size_t pos, size_t length,
                   LongestMatchCache* lmc) {
  size_t i;
  size_t j = 0;
  unsigned bestlength = 0;
  unsigned char* cache;

#if NUM_CACHED_LENGTHS == 0
  return;
#endif

  cache = &lmc->sublen[NUM_CACHED_LENGTHS * pos * 3];
  if (length < 3) return;
  for (i = 3; i <= length; i++) {
    if (i == length || sublen[i] != sublen[i + 1]) {
      cache[j * 3] = i - 3;
      cache[j * 3 + 1] = sublen[i] % 256;
      cache[j * 3 + 2] = (sublen[i] >> 8) % 256;
      bestlength = i;
      j++;
      if (j >= NUM_CACHED_LENGTHS) break;
    }
  }
  if (j < NUM_CACHED_LENGTHS) {
    assert(bestlength == length);
    cache[(NUM_CACHED_LENGTHS - 1) * 3] = bestlength - 3;
  } else {
    assert(bestlength <= length);
  }
  assert(bestlength == MaxCachedSublen(lmc, pos, length));
}

void CacheToSublen(const LongestMatchCache* lmc, size_t pos, size_t length,
                   unsigned short* sublen) {
  size_t i, j;
  unsigned maxlength = MaxCachedSublen(lmc, pos, length);
  unsigned prevlength = 0;
  unsigned char* cache;
#if NUM_CACHED_LENGTHS == 0
  return;
#endif
  if (length < 3) return;
  cache = &lmc->sublen[NUM_CACHED_LENGTHS * pos * 3];
  for (j = 0; j < NUM_CACHED_LENGTHS; j++) {
    unsigned length = cache[j * 3] + 3;
    unsigned dist = cache[j * 3 + 1] + 256 * cache[j * 3 + 2];
    for (i = prevlength; i <= length; i++) {
      sublen[i] = dist;
    }
    if (length == maxlength) break;
    prevlength = length + 1;
  }
}

/*
Returns the length up to which could be stored in the cache.
*/
unsigned MaxCachedSublen(const LongestMatchCache* lmc,
                         size_t pos, size_t length) {
  unsigned char* cache;
#if NUM_CACHED_LENGTHS == 0
  return 0;
#endif
  cache = &lmc->sublen[NUM_CACHED_LENGTHS * pos * 3];
  (void)length;
  if (cache[1] == 0 && cache[2] == 0) return 0;  /* No sublen cached. */
  return cache[(NUM_CACHED_LENGTHS - 1) * 3] + 3;
}

#endif  /* USE_LONGEST_MATCH_CACHE */
