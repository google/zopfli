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

#include "hash.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#define HASH_SHIFT 5
#define HASH_MASK 32767

void ZopfliAllocHash(size_t window_size, ZopfliHash* h) {
  h->head = (int*)malloc(sizeof(*h->head) * 65536);
  h->prev = (unsigned short*)malloc(sizeof(*h->prev) * window_size);
  h->hashval = (int*)malloc(sizeof(*h->hashval) * window_size);

#ifdef ZOPFLI_HASH_SAME
  h->same = (unsigned short*)malloc(sizeof(*h->same) * window_size);
#endif

#ifdef ZOPFLI_HASH_SAME_HASH
  h->head2 = (int*)malloc(sizeof(*h->head2) * 65536);
  h->prev2 = (unsigned short*)malloc(sizeof(*h->prev2) * window_size);
  h->hashval2 = (int*)malloc(sizeof(*h->hashval2) * window_size);
#endif
}

void ZopfliResetHash(size_t window_size, ZopfliHash* h) {
  size_t i;

  h->val = 0;
  /* -1 indicates no head so far. */
  memset(h->head, -1, sizeof(*h->head) * 65536);
  memset(h->hashval, -1, sizeof(*h->hashval) * window_size);
  for (i = 0; i < window_size; ++i) {
    h->prev[i] = i; /* If prev[j] == j, then prev[j] is uninitialized. */
  }

#ifdef ZOPFLI_HASH_SAME
  memset(h->same, 0, sizeof(*h->same) * window_size);
#endif

#ifdef ZOPFLI_HASH_SAME_HASH
  h->val2 = 0;
  memset(h->head2, -1, sizeof(*h->head2) * 65536);
  memset(h->hashval2, -1, sizeof(*h->hashval2) * window_size);
  for (i = 0; i < window_size; ++i) {
    h->prev2[i] = i;
  }
#endif
}

void ZopfliCleanHash(ZopfliHash* h) {
  free(h->head);
  free(h->prev);
  free(h->hashval);

#ifdef ZOPFLI_HASH_SAME_HASH
  free(h->head2);
  free(h->prev2);
  free(h->hashval2);
#endif

#ifdef ZOPFLI_HASH_SAME
  free(h->same);
#endif
}

/*
Update the sliding hash value with the given byte. All calls to this function
must be made on consecutive input characters. Since the hash value exists out
of multiple input bytes, a few warmups with this function are needed initially.
*/
static void UpdateHashValue(ZopfliHash* h, unsigned char c) {
  h->val = (((h->val) << HASH_SHIFT) ^ (c)) & HASH_MASK;
}

void ZopfliUpdateHash(const unsigned char* array, size_t pos, size_t end,
                ZopfliHash* h) {
  unsigned short hpos = pos & ZOPFLI_WINDOW_MASK;
#ifdef ZOPFLI_HASH_SAME
  size_t amount = 0;
#endif

  UpdateHashValue(h, pos + ZOPFLI_MIN_MATCH <= end ?
      array[pos + ZOPFLI_MIN_MATCH - 1] : 0);
  h->hashval[hpos] = h->val;
  if (h->head[h->val] != -1 && h->hashval[h->head[h->val]] == h->val) {
    h->prev[hpos] = h->head[h->val];
  }
  else h->prev[hpos] = hpos;
  h->head[h->val] = hpos;

#ifdef ZOPFLI_HASH_SAME
  /* Update "same". */
  if (h->same[(pos - 1) & ZOPFLI_WINDOW_MASK] > 1) {
    amount = h->same[(pos - 1) & ZOPFLI_WINDOW_MASK] - 1;
  }
  while (pos + amount + 1 < end &&
      array[pos] == array[pos + amount + 1] && amount < (unsigned short)(-1)) {
    amount++;
  }
  h->same[hpos] = amount;
#endif

#ifdef ZOPFLI_HASH_SAME_HASH
  h->val2 = ((h->same[hpos] - ZOPFLI_MIN_MATCH) & 255) ^ h->val;
  h->hashval2[hpos] = h->val2;
  if (h->head2[h->val2] != -1 && h->hashval2[h->head2[h->val2]] == h->val2) {
    h->prev2[hpos] = h->head2[h->val2];
  }
  else h->prev2[hpos] = hpos;
  h->head2[h->val2] = hpos;
#endif
}

void ZopfliWarmupHash(const unsigned char* array, size_t pos, size_t end,
                ZopfliHash* h) {
  UpdateHashValue(h, array[pos + 0]);
  if (pos + 1 < end) UpdateHashValue(h, array[pos + 1]);
}
