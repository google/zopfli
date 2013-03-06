/*
Copyright 2013 Google Inc. All Rights Reserved.
Author: lode@google.com (Lode Vandevenne)
*/

#include "zopfli.h"

#include "deflate.h"
#include "gzip_container.h"
#include "zlib_container.h"

#include <assert.h>

void ZopfliInitOptions(ZopfliOptions* options) {
  options->verbose = 0;
  options->numiterations = 15;
  options->blocksplitting = 1;
  options->blocksplittinglast = 0;
  options->blocksplittingmax = 15;
}

void ZopfliCompress(const ZopfliOptions* options, ZopfliFormat output_type,
                    const unsigned char* in, size_t insize,
                    unsigned char** out, size_t* outsize)
{
  if (output_type == ZOPFLI_FORMAT_GZIP) {
    ZopfliGzipCompress(options, in, insize, out, outsize);
  } else if (output_type == ZOPFLI_FORMAT_ZLIB) {
    ZopfliZlibCompress(options, in, insize, out, outsize);
  } else if (output_type == ZOPFLI_FORMAT_DEFLATE) {
    unsigned char bp = 0;
    ZopfliDeflate(options, 2 /* Dynamic block */, 1,
                  in, insize, &bp, out, outsize);
  } else {
    assert(0);
  }
}
