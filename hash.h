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
The hash for FindLongestMatch of lz77.c.
*/

#ifndef ZOPFLI_HASH_H_
#define ZOPFLI_HASH_H_

#include "util.h"

typedef struct Hash {
  int* head;  /* Hash value to index of its most recent occurance. */
  unsigned short* prev;  /* Index to index of prev. occurance of same hash. */
  int* hashval;  /* Index to hash value at this index. */
  int val;  /* Current hash value. */

#ifdef USE_HASH_SAME_HASH
  /* Fields with similar purpose as the above hash, but for the second hash with
  a value that is calculated differently.  */
  int* head2;  /* Hash value to index of its most recent occurance. */
  unsigned short* prev2;  /* Index to index of prev. occurance of same hash. */
  int* hashval2;  /* Index to hash value at this index. */
  int val2;  /* Current hash value. */
#endif

#ifdef USE_HASH_SAME
  unsigned short* same;  /* Amount of repetitions of same byte after this .*/
#endif
} Hash;

/* Allocates and initializes all fields of Hash. */
void InitHash(size_t window_size, Hash* h);

/* Frees all fields of Hash. */
void CleanHash(Hash* h);

/*
Updates the hash values based on the current position in the array. All calls
to this must be made for consecutive bytes.
*/
void UpdateHash(const unsigned char* array, size_t pos, size_t end, Hash* h);

/*
Prepopulates hash:
Fills in the initial values in the hash, before UpdateHash can be used
correctly.
*/
void WarmupHash(const unsigned char* array, size_t pos, size_t end, Hash* h);

#endif  /* ZOPFLI_HASH_H_ */
