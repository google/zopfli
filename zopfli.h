/*
Copyright 2013 Google Inc. All Rights Reserved.
Author: lode@google.com (Lode Vandevenne)
*/

#ifndef UTIL_COMPRESSION_ZOPFLI_INTERNAL_ZOPFLI_H_
#define UTIL_COMPRESSION_ZOPFLI_INTERNAL_ZOPFLI_H_

#include <stdlib.h> /* for size_t */

/*
Options used throughout the program.
*/
typedef struct ZopfliOptions {
  /* Whether to print output */
  int verbose;

  /*
  Maximum amount of times to rerun forward and backward pass to optimize LZ77
  compression cost. Good values: 10, 15 for small files, 5 for files over
  several MB in size or it will be too slow.
  */
  int numiterations;

  /*
  If true, splits the data in multiple deflate blocks with optimal choice
  for the block boundaries. Block splitting gives better compression. Default:
  true (1).
  */
  int blocksplitting;

  /*
  If true, chooses the optimal block split points only after doing the iterative
  LZ77 compression. If false, chooses the block split points first, then does
  iterative LZ77 on each individual block. Depending on the file, either first
  or last gives the best compression. Default: false (0).
  */
  int blocksplittinglast;

  /*
  Maximum amount of blocks to split into (0 for unlimited, but this can give
  extreme results that hurt compression on some files). Default value: 15.
  */
  int blocksplittingmax;
} ZopfliOptions;

/* Initializes options with default values. */
void ZopfliInitOptions(ZopfliOptions* options);

/* Output format */
typedef enum {
  ZOPFLI_FORMAT_GZIP,
  ZOPFLI_FORMAT_ZLIB,
  ZOPFLI_FORMAT_DEFLATE
} ZopfliFormat;

/*
Compresses according to the given output format and appends the result to the
output.

options: global program options
output_type: the output format to use
out: pointer to the dynamic output array to which the result is appended. Must
  be freed after use
outsize: pointer to the dynamic output array size
*/
void ZopfliCompress(const ZopfliOptions* options, ZopfliFormat output_type,
                    const unsigned char* in, size_t insize,
                    unsigned char** out, size_t* outsize);

#endif  /* UTIL_COMPRESSION_ZOPFLI_INTERNAL_ZOPFLI_H_ */
