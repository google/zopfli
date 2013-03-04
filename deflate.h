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
Modified by madler@alumni.caltech.edu (Mark Adler)
Exposed DeflatePart() as an external function.
*/

#ifndef ZOPFLI_DEFLATE_H_
#define ZOPFLI_DEFLATE_H_

/*
Functions to compress compatible with the deflate specification.
*/

#include "util.h"

/*
Compresses according to the deflate specification and append the compressed
result to the output.
This function will usually output multiple deflate blocks. If final is 1, then
the final bit will be set on the last block.

options: global program options
btype: the deflate block type. Use 2 for best compression.
  -0: non compressed blocks (00)
  -1: blocks with fixed tree (01)
  -2: blocks with dynamic tree (10)
final: whether this is the last section of the input, sets the final bit to the
  last deflate block.
in: the input bytes
insize (Deflate() only): number of input bytes
instart (DeflatePart() only): offset of the start of the data to compress at in
   if instart is not zero, then the data preceding instart will be used as the
   LZ77 dictionary
inend (DeflatePart() only): offset + 1 of the end of the data to compress at in
bp: bit pointer for the output array. This must initially be 0, and for
  consecutive calls must be reused (it can have values from 0-7). This is
  because deflate appends blocks as bit-based data, rather than on byte
  boundaries.
out: pointer to the dynamic output array to which the result is appended. Must
  be freed after use.
outsize: pointer to the dynamic output array size.
*/
void Deflate(const Options* options, int btype, int final,
             const unsigned char* in, size_t insize,
             unsigned char* bp, unsigned char** out, size_t* outsize);
void DeflatePart(const Options* options, int btype, int final,
                 const unsigned char* in, size_t instart, size_t inend,
                 unsigned char* bp, unsigned char** out,
                 size_t* outsize);

/*
Outputs the tree to a dynamic block (btype 10) according to the deflate
specification.
*/
void AddDynamicTree(const unsigned* ll_lengths, const unsigned* d_lengths,
                    unsigned char* bp, unsigned char** out, size_t* outsize);

/*
Calculates block size in bits.
litlens: lz77 lit/lengths
dists: ll77 distances
lstart: start of block
lend: end of block (not inclusive)
*/
double CalculateBlockSize(
    const unsigned short* litlens, const unsigned short* dists,
    size_t lstart, size_t lend, int btype);
#endif  /* ZOPFLI_DEFLATE_H_ */
