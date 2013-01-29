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
Example program using the zopfli compression algorithm library. It is a tool
which behaves like a subset of gzip which can only compress, not decompress.
*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "deflate.h"
#include "util.h"

/*
Loads a file into a memory array.
*/
static void LoadFile(const char* filename,
                     unsigned char** out, size_t* outsize) {
  FILE* file;

  *out = 0;
  *outsize = 0;
  file = fopen(filename, "rb");
  if (!file) return;

  fseek(file , 0 , SEEK_END);
  *outsize = ftell(file);
  rewind(file);

  *out = (unsigned char*)malloc(*outsize);

  if (*outsize && (*out)) {
    size_t testsize = fread(*out, 1, *outsize, file);
    if (testsize != *outsize) {
      /* It could be a directory */
      free(*out);
      *out = 0;
      *outsize = 0;
    }
  }

  assert(!(*outsize) || out);  /* If size is not zero, out must be allocated. */
  fclose(file);
}

/*
Saves a file from a memory array, overwriting the file if it existed.
*/
static void SaveFile(const char* filename,
                     const unsigned char* in, size_t insize) {
  FILE* file = fopen(filename, "wb" );
  assert(file);
  fwrite((char*)in, 1, insize, file);
  fclose(file);
}

/* Table of CRCs of all 8-bit messages. */
static unsigned long crc_table[256];

/* Flag: has the table been computed? Initially false. */
static int crc_table_computed = 0;

/* Makes the table for a fast CRC. */
void MakeCRCTable() {
  unsigned long c;
  int n, k;
  for (n = 0; n < 256; n++) {
    c = (unsigned long) n;
    for (k = 0; k < 8; k++) {
      if (c & 1) {
        c = 0xedb88320L ^ (c >> 1);
      } else {
        c = c >> 1;
      }
    }
    crc_table[n] = c;
  }
  crc_table_computed = 1;
}

/*
Updates a running crc with the bytes buf[0..len-1] and returns
the updated crc. The crc should be initialized to zero.
*/
unsigned long UpdateCRC(unsigned long crc,
                        const unsigned char *buf, size_t len) {
  unsigned long c = crc ^ 0xffffffffL;
  unsigned n;

  if (!crc_table_computed)
    MakeCRCTable();
  for (n = 0; n < len; n++) {
    c = crc_table[(c ^ buf[n]) & 0xff] ^ (c >> 8);
  }
  return c ^ 0xffffffffL;
}

/* Returns the CRC of the bytes buf[0..len-1]. */
unsigned long CRC(const unsigned char* buf, int len) {
  return UpdateCRC(0L, buf, len);
}

/*
Compresses the data according to the gzip specification.
*/
void Gzip(const Options* options,
          const unsigned char* in, size_t insize,
          unsigned char** out, size_t* outsize) {
  unsigned long crcvalue = CRC(in, insize);
  unsigned char bp = 0;

  APPEND_DATA(31, out, outsize);  /* ID1 */
  APPEND_DATA(139, out, outsize);  /* ID2 */
  APPEND_DATA(8, out, outsize);  /* CM */
  APPEND_DATA(0, out, outsize);  /* FLG */
  /* MTIME */
  APPEND_DATA(0, out, outsize);
  APPEND_DATA(0, out, outsize);
  APPEND_DATA(0, out, outsize);
  APPEND_DATA(0, out, outsize);

  APPEND_DATA(2, out, outsize);  /* XFL, 2 indicates best compression. */
  APPEND_DATA(3, out, outsize);  /* OS follows Unix conventions. */

  Deflate(options, 2 /* Dynamic block */, 1, in, insize, &bp, out, outsize);

  /* CRC */
  APPEND_DATA(crcvalue % 256, out, outsize);
  APPEND_DATA((crcvalue >> 8) % 256, out, outsize);
  APPEND_DATA((crcvalue >> 16) % 256, out, outsize);
  APPEND_DATA((crcvalue >> 24) % 256, out, outsize);

  /* ISIZE */
  APPEND_DATA(insize % 256, out, outsize);
  APPEND_DATA((insize >> 8) % 256, out, outsize);
  APPEND_DATA((insize >> 16) % 256, out, outsize);
  APPEND_DATA((insize >> 24) % 256, out, outsize);

  if (options->verbose) {
    fprintf(stderr,
            "Original Size: %d, Compressed: %d, Compression: %f%% Removed\n",
            (int)insize, (int)*outsize,
            100.0f * (float)(insize - *outsize) / (float)insize);
  }
}

/*
outfilename: filename to write output to, or 0 to write to stdout instead
*/
void GzipFile(const Options* options,
              const char* infilename,
              const char* outfilename) {
  unsigned char* in;
  size_t insize;
  unsigned char* out = 0;
  size_t outsize = 0;
  LoadFile(infilename, &in, &insize);
  if (insize == 0) {
    fprintf(stderr, "Invalid filename: %s\n", infilename);
    return;
  }
  Gzip(options, in, insize, &out, &outsize);
  if (outfilename) {
    SaveFile(outfilename, out, outsize);
  } else {
    size_t i;
    for (i = 0; i < outsize; i++) {
      /* Works only if terminal does not convert newlines. */
      printf("%c", out[i]);
    }
  }

  free(out);
  free(in);
}

/*
Add two strings together. Size does not matter. Result must be freed.
*/
static char* AddStrings(const char* str1, const char* str2) {
  size_t len = strlen(str1) + strlen(str2);
  char* result = (char*)malloc(len + 1);
  if (!result) exit(-1); /* Allocation failed. */
  strcpy(result, str1);
  strcat(result, str2);
  return result;
}

static char StringsEqual(const char* str1, const char* str2) {
  return strcmp(str1, str2) == 0;
}

int main(int argc, char* argv[]) {
  Options options;
  const char* filename = 0;
  int output_to_stdout = 0;
  int i;

  InitOptions(&options);

  for (i = 1; i < argc; i++) {
    if (StringsEqual(argv[i], "-v")) options.verbose = 1;
    else if (StringsEqual(argv[i], "-c")) output_to_stdout = 1;
    else if (StringsEqual(argv[i], "-i5")) options.numiterations = 5;
    else if (StringsEqual(argv[i], "-i10")) options.numiterations = 10;
    else if (StringsEqual(argv[i], "-i15")) options.numiterations = 15;
    else if (StringsEqual(argv[i], "-i25")) options.numiterations = 25;
    else if (StringsEqual(argv[i], "-i50")) options.numiterations = 50;
    else if (StringsEqual(argv[i], "-i100")) options.numiterations = 100;
    else if (StringsEqual(argv[i], "-i250")) options.numiterations = 250;
    else if (StringsEqual(argv[i], "-i500")) options.numiterations = 500;
    else if (StringsEqual(argv[i], "-i1000")) options.numiterations = 1000;
    else if (StringsEqual(argv[i], "-h")) {
      fprintf(stderr, "Usage: zopfli [OPTION]... FILE\n"
          "  -h    gives this help\n"
          "  -c    write the result on standard output, instead of disk"
          " filename + '.gz'\n"
          "  -v    verbose mode\n"
          "  -i5 less compression, but faster\n"
          "  -i10 less compression, but faster\n"
          "  -i15 default compression, 15 iterations\n"
          "  -i25 more compression, but slower\n"
          "  -i50 more compression, but slower\n"
          "  -i100 more compression, but slower\n"
          "  -i250 more compression, but slower\n"
          "  -i500 more compression, but slower\n"
          "  -i1000 more compression, but slower\n");
      return 0;
    }
  }

  for (i = 1; i < argc; i++) {
    if (argv[i][0] != '-') {
      char* outfilename;
      filename = argv[i];
      outfilename = output_to_stdout ?
          0 : AddStrings(filename, ".gz");
      if (options.verbose && outfilename) {
        fprintf(stderr, "Saving to: %s\n", outfilename);
      }
      GzipFile(&options, filename, outfilename);
      free(outfilename);
    }
  }

  if (!filename) {
    fprintf(stderr,
            "Please provide filename\nFor help, type: %s -h\n", argv[0]);
  }

  return 0;
}
