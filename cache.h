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
The cache that speeds up FindLongestMatch of lz77.c.
*/

#ifndef ZOPFLI_CACHE_H_
#define ZOPFLI_CACHE_H_

#include "util.h"

#ifdef USE_LONGEST_MATCH_CACHE

/*
Cache used by FindLongestMatch to remember previously found length/dist values.
This is needed because the squeeze runs will ask these values multiple times for
the same position.
Uses large amounts of memory, since it has to remember the distance belonging
to every possible shorter-than-the-best length (the so called "sublen" array).
*/
typedef struct LongestMatchCache {
  unsigned short* length;
  unsigned short* dist;
  unsigned char* sublen; /* For each length, the distance */
} LongestMatchCache;

/* Initializes the LongestMatchCache. */
void InitLongestMatchCache(size_t blocksize, LongestMatchCache* lmc);

/* Frees up the memory of the LongestMatchCache. */
void CleanLongestMatchCache(LongestMatchCache* lmc);

/* Stores sublen array in the cache. */
void SublenToCache(const unsigned short* sublen, size_t pos, size_t length,
                   LongestMatchCache* lmc);

/* Extracts sublen array from the cache. */
void CacheToSublen(const LongestMatchCache* lmc, size_t pos, size_t length,
                   unsigned short* sublen);
/* Returns the length up to which could be stored in the cache. */
unsigned MaxCachedSublen(const LongestMatchCache* lmc,
                         size_t pos, size_t length);

#endif  /* USE_LONGEST_MATCH_CACHE */

#endif  /* ZOPFLI_CACHE_H_ */
