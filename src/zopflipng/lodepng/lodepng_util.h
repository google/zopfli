/*
LodePNG Utils

Copyright (c) 2005-2012 Lode Vandevenne

This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

    1. The origin of this software must not be misrepresented; you must not
    claim that you wrote the original software. If you use this software
    in a product, an acknowledgment in the product documentation would be
    appreciated but is not required.

    2. Altered source versions must be plainly marked as such, and must not be
    misrepresented as being the original software.

    3. This notice may not be removed or altered from any source
    distribution.
*/

/*
Extra C++ utilities for LodePNG, for convenience.
*/

#include <string>
#include <vector>
#include "lodepng.h"

#pragma once

namespace lodepng
{

/*
Returns info from the header of the PNG by value, purely for convenience.
Does NOT check for errors. Returns bogus info if the PNG has an error.
Does not require cleanup of allocated memory because no palette or text chunk
info is in the LodePNGInfo object after checking only the header of the PNG.
*/
LodePNGInfo getPNGHeaderInfo(const std::vector<unsigned char>& png);

/*
Get the names and sizes of all chunks in the PNG file.
Returns 0 if ok, non-0 if error happened.
*/
unsigned getChunkInfo(std::vector<std::string>& names, std::vector<size_t>& sizes,
                      const std::vector<unsigned char>& png);

/*
Returns the names and full chunks (including the name and everything else that
makes up the chunk) for all chunks except IHDR, PLTE, IDAT and IEND.
It separates the chunks into 3 separate lists, representing the chunks between
certain critical chunks: 0: IHDR-PLTE, 1: PLTE-IDAT, 2: IDAT-IEND
Returns 0 if ok, non-0 if error happened.
*/
unsigned getChunks(std::vector<std::string> names[3],
                   std::vector<std::vector<unsigned char> > chunks[3],
                   const std::vector<unsigned char>& png);

/*
Inserts chunks into the given png file. The chunks must be fully encoded,
including length, type, content and CRC.
The array index determines where it goes:
0: between IHDR and PLTE, 1: between PLTE and IDAT, 2: between IDAT and IEND.
They're appended at the end of those locations within the PNG.
Returns 0 if ok, non-0 if error happened.
*/
unsigned insertChunks(std::vector<unsigned char>& png,
                      const std::vector<std::vector<unsigned char> > chunks[3]);

/*
Get the filtertypes of each scanline in this PNG file.
Returns 0 if ok, 1 if PNG decoding error happened.

For a non-interlaced PNG, it returns one filtertype per scanline, in order.

For interlaced PNGs, it returns a result as if it's not interlaced. It returns
one filtertype per scanline, in order. The values match pass 6 and 7 of the
Adam7 interlacing, alternating between the two, so that the values correspond
the most to their scanlines.
*/
unsigned getFilterTypes(std::vector<unsigned char>& filterTypes, const std::vector<unsigned char>& png);

/*
Get the filtertypes of each scanline in every interlace pass this PNG file.
Returns 0 if ok, 1 if PNG decoding error happened.

For a non-interlaced PNG, it returns one filtertype per scanline, in order, in
a single std::vector in filterTypes.

For an interlaced PNG, it returns 7 std::vectors in filterTypes, one for each
Adam7 pass. The amount of values per pass can be calculated as follows, where
w and h are the size of the image and all divisions are integer divisions:
pass 1: (h + 7) / 8
pass 2: w <= 4 ? 0 : (h + 7) / 8
pass 3: h <= 4 ? 0 : (h + 7) / 8
pass 4: w <= 2 ? 0 : (h + 3) / 4
pass 5: h <= 2 ? 0 : (h + 3) / 4
pass 6: w <= 1 ? 0 : (h + 1) / 2
pass 7: h <= 1 ? 0 : (h + 1) / 2
*/
unsigned getFilterTypesInterlaced(std::vector<std::vector<unsigned char> >& filterTypes,
                                  const std::vector<unsigned char>& png);

/*
Returns the value of the i-th pixel in an image with 1, 2, 4 or 8-bit color.
E.g. if bits is 4 and i is 5, it returns the 5th nibble (4-bit group), which
is the second half of the 3th byte, in big endian (PNG's endian order).
*/
int getPaletteValue(const unsigned char* data, size_t i, int bits);

/*
The information for extractZlibInfo.
*/
struct ZlibBlockInfo
{
  int btype; //block type (0-2)
  size_t compressedbits; //size of compressed block in bits
  size_t uncompressedbytes; //size of uncompressed block in bytes

  // only filled in for block type 2
  size_t treebits; //encoded tree size in bits
  int hlit; //the HLIT value that was filled in for this tree
  int hdist; //the HDIST value that was filled in for this tree
  int hclen; //the HCLEN value that was filled in for this tree
  std::vector<int> clcl; //19 code length code lengths (compressed tree's tree)
  std::vector<int> treecodes; //N tree codes, with values 0-18. Values 17 or 18 are followed by the repetition value.
  std::vector<int> litlenlengths; //288 code lengths for lit/len symbols
  std::vector<int> distlengths; //32 code lengths for dist symbols

  // only filled in for block types 1 or 2
  std::vector<int> lz77_lcode; //LZ77 codes. 0-255: literals. 256: end symbol. 257-285: length code of length/dist pairs
  // the next vectors have the same size as lz77_lcode, but an element only has meaningful value if lz77_lcode contains a length code.
  std::vector<int> lz77_dcode;
  std::vector<int> lz77_lbits;
  std::vector<int> lz77_dbits;
  std::vector<int> lz77_lvalue;
  std::vector<int> lz77_dvalue;
  size_t numlit; //number of lit codes in this block
  size_t numlen; //number of len codes in this block
};

//Extracts all info needed from a PNG file to reconstruct the zlib compression exactly.
void extractZlibInfo(std::vector<ZlibBlockInfo>& zlibinfo, const std::vector<unsigned char>& in);

} // namespace lodepng
