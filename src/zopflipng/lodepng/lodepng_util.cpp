/*
LodePNG Utils

Copyright (c) 2005-2014 Lode Vandevenne

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

#include "lodepng_util.h"
#include <iostream>

namespace lodepng
{

LodePNGInfo getPNGHeaderInfo(const std::vector<unsigned char>& png)
{
  unsigned w, h;
  lodepng::State state;
  lodepng_inspect(&w, &h, &state, &png[0], png.size());
  return state.info_png;
}

unsigned getChunkInfo(std::vector<std::string>& names, std::vector<size_t>& sizes,
                      const std::vector<unsigned char>& png)
{
  // Listing chunks is based on the original file, not the decoded png info.
  const unsigned char *chunk, *begin, *end, *next;
  end = &png.back() + 1;
  begin = chunk = &png.front() + 8;

  while(chunk + 8 < end && chunk >= begin)
  {
    char type[5];
    lodepng_chunk_type(type, chunk);
    if(std::string(type).size() != 4) return 1;

    unsigned length = lodepng_chunk_length(chunk);
    if(chunk + length + 12 > end) return 1;
    names.push_back(type);
    sizes.push_back(length);

    next = lodepng_chunk_next_const(chunk);
    if (next <= chunk) return 1; // integer overflow
    chunk = next;
  }
  return 0;
}

unsigned getChunks(std::vector<std::string> names[3],
                   std::vector<std::vector<unsigned char> > chunks[3],
                   const std::vector<unsigned char>& png)
{
  const unsigned char *chunk, *next, *begin, *end;
  end = &png.back() + 1;
  begin = chunk = &png.front() + 8;

  int location = 0;

  while(chunk + 8 < end && chunk >= begin)
  {
    char type[5];
    lodepng_chunk_type(type, chunk);
    std::string name(type);
    if(name.size() != 4) return 1;

    next = lodepng_chunk_next_const(chunk);
    if (next <= chunk) return 1; // integer overflow

    if(name == "IHDR")
    {
      location = 0;
    }
    else if(name == "PLTE")
    {
      location = 1;
    }
    else if(name == "IDAT")
    {
      location = 2;
    }
    else if(name == "IEND")
    {
      break; // anything after IEND is not part of the PNG or the 3 groups here.
    }
    else
    {
      if(next > end) return 1; // invalid chunk, content too far
      names[location].push_back(name);
      chunks[location].push_back(std::vector<unsigned char>(chunk, next));
    }

    chunk = next;
  }
  return 0;
}


unsigned insertChunks(std::vector<unsigned char>& png,
                      const std::vector<std::vector<unsigned char> > chunks[3])
{
  const unsigned char *chunk, *next, *begin, *end;
  end = &png.back() + 1;
  begin = chunk = &png.front() + 8;

  long l0 = 0; //location 0: IHDR-l0-PLTE (or IHDR-l0-l1-IDAT)
  long l1 = 0; //location 1: PLTE-l1-IDAT (or IHDR-l0-l1-IDAT)
  long l2 = 0; //location 2: IDAT-l2-IEND

  while(chunk + 8 < end && chunk >= begin)
  {
    char type[5];
    lodepng_chunk_type(type, chunk);
    std::string name(type);
    if(name.size() != 4) return 1;

    next = lodepng_chunk_next_const(chunk);
    if (next <= chunk) return 1; // integer overflow

    if(name == "PLTE")
    {
      if(l0 == 0) l0 = chunk - begin + 8;
    }
    else if(name == "IDAT")
    {
      if(l0 == 0) l0 = chunk - begin + 8;
      if(l1 == 0) l1 = chunk - begin + 8;
    }
    else if(name == "IEND")
    {
      if(l2 == 0) l2 = chunk - begin + 8;
    }

    chunk = next;
  }

  std::vector<unsigned char> result;
  result.insert(result.end(), png.begin(), png.begin() + l0);
  for(size_t i = 0; i < chunks[0].size(); i++) result.insert(result.end(), chunks[0][i].begin(), chunks[0][i].end());
  result.insert(result.end(), png.begin() + l0, png.begin() + l1);
  for(size_t i = 0; i < chunks[1].size(); i++) result.insert(result.end(), chunks[1][i].begin(), chunks[1][i].end());
  result.insert(result.end(), png.begin() + l1, png.begin() + l2);
  for(size_t i = 0; i < chunks[2].size(); i++) result.insert(result.end(), chunks[2][i].begin(), chunks[2][i].end());
  result.insert(result.end(), png.begin() + l2, png.end());

  png = result;
  return 0;
}

unsigned getFilterTypesInterlaced(std::vector<std::vector<unsigned char> >& filterTypes,
                                  const std::vector<unsigned char>& png)
{
  //Get color type and interlace type
  lodepng::State state;
  unsigned w, h;
  unsigned error;
  error = lodepng_inspect(&w, &h, &state, &png[0], png.size());

  if(error) return 1;

  //Read literal data from all IDAT chunks
  const unsigned char *chunk, *begin, *end, *next;
  end = &png.back() + 1;
  begin = chunk = &png.front() + 8;

  std::vector<unsigned char> zdata;

  while(chunk + 8 < end && chunk >= begin)
  {
    char type[5];
    lodepng_chunk_type(type, chunk);
    if(std::string(type).size() != 4) return 1; //Probably not a PNG file

    if(std::string(type) == "IDAT")
    {
      const unsigned char* cdata = lodepng_chunk_data_const(chunk);
      unsigned clength = lodepng_chunk_length(chunk);
      if(chunk + clength + 12 > end || clength > png.size() || chunk + clength + 12 < begin) {
        // corrupt chunk length
        return 1;
      }

      for(unsigned i = 0; i < clength; i++)
      {
        zdata.push_back(cdata[i]);
      }
    }

    next = lodepng_chunk_next_const(chunk);
    if (next <= chunk) return 1; // integer overflow
    chunk = next;
  }

  //Decompress all IDAT data
  std::vector<unsigned char> data;
  error = lodepng::decompress(data, &zdata[0], zdata.size());

  if(error) return 1;

  if(state.info_png.interlace_method == 0)
  {
    filterTypes.resize(1);

    //A line is 1 filter byte + all pixels
    size_t linebytes = 1 + lodepng_get_raw_size(w, 1, &state.info_png.color);

    for(size_t i = 0; i < data.size(); i += linebytes)
    {
      filterTypes[0].push_back(data[i]);
    }
  }
  else
  {
    //Interlaced
    filterTypes.resize(7);
    static const unsigned ADAM7_IX[7] = { 0, 4, 0, 2, 0, 1, 0 }; /*x start values*/
    static const unsigned ADAM7_IY[7] = { 0, 0, 4, 0, 2, 0, 1 }; /*y start values*/
    static const unsigned ADAM7_DX[7] = { 8, 8, 4, 4, 2, 2, 1 }; /*x delta values*/
    static const unsigned ADAM7_DY[7] = { 8, 8, 8, 4, 4, 2, 2 }; /*y delta values*/
    size_t pos = 0;
    for(size_t j = 0; j < 7; j++)
    {
      unsigned w2 = (w - ADAM7_IX[j] + ADAM7_DX[j] - 1) / ADAM7_DX[j];
      unsigned h2 = (h - ADAM7_IY[j] + ADAM7_DY[j] - 1) / ADAM7_DY[j];
      if(ADAM7_IX[j] >= w) w2 = 0;
      if(ADAM7_IY[j] >= h) h2 = 0;
      size_t linebytes = 1 + lodepng_get_raw_size(w2, 1, &state.info_png.color);
      for(size_t i = 0; i < h2; i++)
      {
        filterTypes[j].push_back(data[pos]);
        pos += linebytes;
      }
    }
  }
  return 0; /* OK */
}


unsigned getFilterTypes(std::vector<unsigned char>& filterTypes, const std::vector<unsigned char>& png)
{
  std::vector<std::vector<unsigned char> > passes;
  unsigned error = getFilterTypesInterlaced(passes, png);
  if(error) return error;

  if(passes.size() == 1)
  {
    filterTypes.swap(passes[0]);
  }
  else
  {
    lodepng::State state;
    unsigned w, h;
    lodepng_inspect(&w, &h, &state, &png[0], png.size());
    /*
    Interlaced. Simplify it: put pass 6 and 7 alternating in the one vector so
    that one filter per scanline of the uninterlaced image is given, with that
    filter corresponding the closest to what it would be for non-interlaced
    image.
    */
    for(size_t i = 0; i < h; i++)
    {
      filterTypes.push_back(i % 2 == 0 ? passes[5][i / 2] : passes[6][i / 2]);
    }
  }
  return 0; /* OK */
}

int getPaletteValue(const unsigned char* data, size_t i, int bits)
{
  if(bits == 8) return data[i];
  else if(bits == 4) return (data[i / 2] >> ((i % 2) * 4)) & 15;
  else if(bits == 2) return (data[i / 4] >> ((i % 4) * 2)) & 3;
  else if(bits == 1) return (data[i / 8] >> (i % 8)) & 1;
  else return 0;
}

//This uses a stripped down version of picoPNG to extract detailed zlib information while decompressing.
static const unsigned long LENBASE[29] =
    {3,4,5,6,7,8,9,10,11,13,15,17,19,23,27,31,35,43,51,59,67,83,99,115,131,163,195,227,258};
static const unsigned long LENEXTRA[29] =
    {0,0,0,0,0,0,0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4,  4,  5,  5,  5,  5,  0};
static const unsigned long DISTBASE[30] =
    {1,2,3,4,5,7,9,13,17,25,33,49,65,97,129,193,257,385,513,769,1025,1537,2049,3073,4097,6145,8193,12289,16385,24577};
static const unsigned long DISTEXTRA[30] =
    {0,0,0,0,1,1,2, 2, 3, 3, 4, 4, 5, 5,  6,  6,  7,  7,  8,  8,   9,   9,  10,  10,  11,  11,  12,   12,   13,   13};
static const unsigned long CLCL[19] =
    {16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15}; //code length code lengths

struct ExtractZlib // Zlib decompression and information extraction
{
  std::vector<ZlibBlockInfo>* zlibinfo;
  ExtractZlib(std::vector<ZlibBlockInfo>* info) : zlibinfo(info) {};
  int error;

  unsigned long readBitFromStream(size_t& bitp, const unsigned char* bits)
  {
    unsigned long result = (bits[bitp >> 3] >> (bitp & 0x7)) & 1;
    bitp++;
    return result;
  }

  unsigned long readBitsFromStream(size_t& bitp, const unsigned char* bits, size_t nbits)
  {
    unsigned long result = 0;
    for(size_t i = 0; i < nbits; i++) result += (readBitFromStream(bitp, bits)) << i;
    return result;
  }

  struct HuffmanTree
  {
    int makeFromLengths(const std::vector<unsigned long>& bitlen, unsigned long maxbitlen)
    { //make tree given the lengths
      unsigned long numcodes = (unsigned long)(bitlen.size()), treepos = 0, nodefilled = 0;
      std::vector<unsigned long> tree1d(numcodes), blcount(maxbitlen + 1, 0), nextcode(maxbitlen + 1, 0);
      //count number of instances of each code length
      for(unsigned long bits = 0; bits < numcodes; bits++) blcount[bitlen[bits]]++;
      for(unsigned long bits = 1; bits <= maxbitlen; bits++)
      {
        nextcode[bits] = (nextcode[bits - 1] + blcount[bits - 1]) << 1;
      }
      //generate all the codes
      for(unsigned long n = 0; n < numcodes; n++) if(bitlen[n] != 0) tree1d[n] = nextcode[bitlen[n]]++;
      tree2d.clear(); tree2d.resize(numcodes * 2, 32767); //32767 here means the tree2d isn't filled there yet
      for(unsigned long n = 0; n < numcodes; n++) //the codes
      for(unsigned long i = 0; i < bitlen[n]; i++) //the bits for this code
      {
        unsigned long bit = (tree1d[n] >> (bitlen[n] - i - 1)) & 1;
        if(treepos > numcodes - 2) return 55;
        if(tree2d[2 * treepos + bit] == 32767) //not yet filled in
        {
          if(i + 1 == bitlen[n])
          {
            //last bit
            tree2d[2 * treepos + bit] = n;
            treepos = 0;
          }
          else
          {
            //addresses are encoded as values > numcodes
            tree2d[2 * treepos + bit] = ++nodefilled + numcodes;
            treepos = nodefilled;
          }
        }
        else treepos = tree2d[2 * treepos + bit] - numcodes; //subtract numcodes from address to get address value
      }
      return 0;
    }
    int decode(bool& decoded, unsigned long& result, size_t& treepos, unsigned long bit) const
    { //Decodes a symbol from the tree
      unsigned long numcodes = (unsigned long)tree2d.size() / 2;
      if(treepos >= numcodes) return 11; //error: you appeared outside the codetree
      result = tree2d[2 * treepos + bit];
      decoded = (result < numcodes);
      treepos = decoded ? 0 : result - numcodes;
      return 0;
    }
    //2D representation of a huffman tree: one dimension is "0" or "1", the other contains all nodes and leaves.
    std::vector<unsigned long> tree2d;
  };

  void inflate(std::vector<unsigned char>& out, const std::vector<unsigned char>& in, size_t inpos = 0)
  {
    size_t bp = 0, pos = 0; //bit pointer and byte pointer
    error = 0;
    unsigned long BFINAL = 0;
    while(!BFINAL && !error)
    {
      size_t uncomprblockstart = pos;
      size_t bpstart = bp;
      if(bp >> 3 >= in.size()) { error = 52; return; } //error, bit pointer will jump past memory
      BFINAL = readBitFromStream(bp, &in[inpos]);
      unsigned long BTYPE = readBitFromStream(bp, &in[inpos]); BTYPE += 2 * readBitFromStream(bp, &in[inpos]);
      zlibinfo->resize(zlibinfo->size() + 1);
      zlibinfo->back().btype = BTYPE;
      if(BTYPE == 3) { error = 20; return; } //error: invalid BTYPE
      else if(BTYPE == 0) inflateNoCompression(out, &in[inpos], bp, pos, in.size());
      else inflateHuffmanBlock(out, &in[inpos], bp, pos, in.size(), BTYPE);
      size_t uncomprblocksize = pos - uncomprblockstart;
      zlibinfo->back().compressedbits = bp - bpstart;
      zlibinfo->back().uncompressedbytes = uncomprblocksize;
    }
  }

  void generateFixedTrees(HuffmanTree& tree, HuffmanTree& treeD) //get the tree of a deflated block with fixed tree
  {
    std::vector<unsigned long> bitlen(288, 8), bitlenD(32, 5);;
    for(size_t i = 144; i <= 255; i++) bitlen[i] = 9;
    for(size_t i = 256; i <= 279; i++) bitlen[i] = 7;
    tree.makeFromLengths(bitlen, 15);
    treeD.makeFromLengths(bitlenD, 15);
  }

  //the code tree for Huffman codes, dist codes, and code length codes
  HuffmanTree codetree, codetreeD, codelengthcodetree;
  unsigned long huffmanDecodeSymbol(const unsigned char* in, size_t& bp, const HuffmanTree& tree, size_t inlength)
  {
    //decode a single symbol from given list of bits with given code tree. return value is the symbol
    bool decoded; unsigned long ct;
    for(size_t treepos = 0;;)
    {
      if((bp & 0x07) == 0 && (bp >> 3) > inlength) { error = 10; return 0; } //error: end reached without endcode
      error = tree.decode(decoded, ct, treepos, readBitFromStream(bp, in));
      if(error) return 0; //stop, an error happened
      if(decoded) return ct;
    }
  }

  void getTreeInflateDynamic(HuffmanTree& tree, HuffmanTree& treeD,
                             const unsigned char* in, size_t& bp, size_t inlength)
  {
    size_t bpstart = bp;
    //get the tree of a deflated block with dynamic tree, the tree itself is also Huffman compressed with a known tree
    std::vector<unsigned long> bitlen(288, 0), bitlenD(32, 0);
    if(bp >> 3 >= inlength - 2) { error = 49; return; } //the bit pointer is or will go past the memory
    size_t HLIT =  readBitsFromStream(bp, in, 5) + 257; //number of literal/length codes + 257
    size_t HDIST = readBitsFromStream(bp, in, 5) + 1; //number of dist codes + 1
    size_t HCLEN = readBitsFromStream(bp, in, 4) + 4; //number of code length codes + 4
    zlibinfo->back().hlit = HLIT - 257;
    zlibinfo->back().hdist = HDIST - 1;
    zlibinfo->back().hclen = HCLEN - 4;
    std::vector<unsigned long> codelengthcode(19); //lengths of tree to decode the lengths of the dynamic tree
    for(size_t i = 0; i < 19; i++) codelengthcode[CLCL[i]] = (i < HCLEN) ? readBitsFromStream(bp, in, 3) : 0;
    //code length code lengths
    for(size_t i = 0; i < codelengthcode.size(); i++) zlibinfo->back().clcl.push_back(codelengthcode[i]);
    error = codelengthcodetree.makeFromLengths(codelengthcode, 7); if(error) return;
    size_t i = 0, replength;
    while(i < HLIT + HDIST)
    {
      unsigned long code = huffmanDecodeSymbol(in, bp, codelengthcodetree, inlength); if(error) return;
      zlibinfo->back().treecodes.push_back(code); //tree symbol code
      if(code <= 15)  { if(i < HLIT) bitlen[i++] = code; else bitlenD[i++ - HLIT] = code; } //a length code
      else if(code == 16) //repeat previous
      {
        if(bp >> 3 >= inlength) { error = 50; return; } //error, bit pointer jumps past memory
        replength = 3 + readBitsFromStream(bp, in, 2);
        unsigned long value; //set value to the previous code
        if((i - 1) < HLIT) value = bitlen[i - 1];
        else value = bitlenD[i - HLIT - 1];
        for(size_t n = 0; n < replength; n++) //repeat this value in the next lengths
        {
          if(i >= HLIT + HDIST) { error = 13; return; } //error: i is larger than the amount of codes
          if(i < HLIT) bitlen[i++] = value; else bitlenD[i++ - HLIT] = value;
        }
      }
      else if(code == 17) //repeat "0" 3-10 times
      {
        if(bp >> 3 >= inlength) { error = 50; return; } //error, bit pointer jumps past memory
        replength = 3 + readBitsFromStream(bp, in, 3);
        zlibinfo->back().treecodes.push_back(replength); //tree symbol code repetitions
        for(size_t n = 0; n < replength; n++) //repeat this value in the next lengths
        {
          if(i >= HLIT + HDIST) { error = 14; return; } //error: i is larger than the amount of codes
          if(i < HLIT) bitlen[i++] = 0; else bitlenD[i++ - HLIT] = 0;
        }
      }
      else if(code == 18) //repeat "0" 11-138 times
      {
        if(bp >> 3 >= inlength) { error = 50; return; } //error, bit pointer jumps past memory
        replength = 11 + readBitsFromStream(bp, in, 7);
        zlibinfo->back().treecodes.push_back(replength); //tree symbol code repetitions
        for(size_t n = 0; n < replength; n++) //repeat this value in the next lengths
        {
          if(i >= HLIT + HDIST) { error = 15; return; } //error: i is larger than the amount of codes
          if(i < HLIT) bitlen[i++] = 0; else bitlenD[i++ - HLIT] = 0;
        }
      }
      else { error = 16; return; } //error: somehow an unexisting code appeared. This can never happen.
    }
    if(bitlen[256] == 0) { error = 64; return; } //the length of the end code 256 must be larger than 0
    error = tree.makeFromLengths(bitlen, 15);
    if(error) return; //now we've finally got HLIT and HDIST, so generate the code trees, and the function is done
    error = treeD.makeFromLengths(bitlenD, 15);
    if(error) return;
    zlibinfo->back().treebits = bp - bpstart;
    //lit/len/end symbol lengths
    for(size_t j = 0; j < bitlen.size(); j++) zlibinfo->back().litlenlengths.push_back(bitlen[j]);
    //dist lengths
    for(size_t j = 0; j < bitlenD.size(); j++) zlibinfo->back().distlengths.push_back(bitlenD[j]);
  }

  void inflateHuffmanBlock(std::vector<unsigned char>& out,
                           const unsigned char* in, size_t& bp, size_t& pos, size_t inlength, unsigned long btype)
  {
    size_t numcodes = 0, numlit = 0, numlen = 0; //for logging
    if(btype == 1) { generateFixedTrees(codetree, codetreeD); }
    else if(btype == 2) { getTreeInflateDynamic(codetree, codetreeD, in, bp, inlength); if(error) return; }
    for(;;)
    {
      unsigned long code = huffmanDecodeSymbol(in, bp, codetree, inlength); if(error) return;
      numcodes++;
      zlibinfo->back().lz77_lcode.push_back(code); //output code
      zlibinfo->back().lz77_dcode.push_back(0);
      zlibinfo->back().lz77_lbits.push_back(0);
      zlibinfo->back().lz77_dbits.push_back(0);
      zlibinfo->back().lz77_lvalue.push_back(0);
      zlibinfo->back().lz77_dvalue.push_back(0);

      if(code == 256) break; //end code
      else if(code <= 255) //literal symbol
      {
        out.push_back((unsigned char)(code));
        pos++;
        numlit++;
      }
      else if(code >= 257 && code <= 285) //length code
      {
        size_t length = LENBASE[code - 257], numextrabits = LENEXTRA[code - 257];
        if((bp >> 3) >= inlength) { error = 51; return; } //error, bit pointer will jump past memory
        length += readBitsFromStream(bp, in, numextrabits);
        unsigned long codeD = huffmanDecodeSymbol(in, bp, codetreeD, inlength); if(error) return;
        if(codeD > 29) { error = 18; return; } //error: invalid dist code (30-31 are never used)
        unsigned long dist = DISTBASE[codeD], numextrabitsD = DISTEXTRA[codeD];
        if((bp >> 3) >= inlength) { error = 51; return; } //error, bit pointer will jump past memory
        dist += readBitsFromStream(bp, in, numextrabitsD);
        size_t start = pos, back = start - dist; //backwards
        for(size_t i = 0; i < length; i++)
        {
          out.push_back(out[back++]);
          pos++;
          if(back >= start) back = start - dist;
        }
        numlen++;
        zlibinfo->back().lz77_dcode.back() = codeD; //output distance code
        zlibinfo->back().lz77_lbits.back() = numextrabits; //output length extra bits
        zlibinfo->back().lz77_dbits.back() = numextrabitsD; //output dist extra bits
        zlibinfo->back().lz77_lvalue.back() = length; //output length
        zlibinfo->back().lz77_dvalue.back() = dist; //output dist
      }
    }
    zlibinfo->back().numlit = numlit; //output number of literal symbols
    zlibinfo->back().numlen = numlen; //output number of length symbols
  }

  void inflateNoCompression(std::vector<unsigned char>& out,
                            const unsigned char* in, size_t& bp, size_t& pos, size_t inlength)
  {
    while((bp & 0x7) != 0) bp++; //go to first boundary of byte
    size_t p = bp / 8;
    if(p >= inlength - 4) { error = 52; return; } //error, bit pointer will jump past memory
    unsigned long LEN = in[p] + 256u * in[p + 1], NLEN = in[p + 2] + 256u * in[p + 3]; p += 4;
    if(LEN + NLEN != 65535) { error = 21; return; } //error: NLEN is not one's complement of LEN
    if(p + LEN > inlength) { error = 23; return; } //error: reading outside of in buffer
    for(unsigned long n = 0; n < LEN; n++)
    {
      out.push_back(in[p++]); //read LEN bytes of literal data
      pos++;
    }
    bp = p * 8;
  }

  int decompress(std::vector<unsigned char>& out, const std::vector<unsigned char>& in) //returns error value
  {
    if(in.size() < 2) { return 53; } //error, size of zlib data too small
    //error: 256 * in[0] + in[1] must be a multiple of 31, the FCHECK value is supposed to be made that way
    if((in[0] * 256 + in[1]) % 31 != 0) { return 24; }
    unsigned long CM = in[0] & 15, CINFO = (in[0] >> 4) & 15, FDICT = (in[1] >> 5) & 1;
    //error: only compression method 8: inflate with sliding window of 32k is supported by the PNG spec
    if(CM != 8 || CINFO > 7) { return 25; }
    //error: the PNG spec says about the zlib stream: "The additional flags shall not specify a preset dictionary."
    if(FDICT != 0) { return 26; }
    inflate(out, in, 2);
    return error; //note: adler32 checksum was skipped and ignored
  }
};

struct ExtractPNG //PNG decoding and information extraction
{
  std::vector<ZlibBlockInfo>* zlibinfo;
  ExtractPNG(std::vector<ZlibBlockInfo>* info) : zlibinfo(info) {};
  int error;
  void decode(const unsigned char* in, size_t size)
  {
    error = 0;
    if(size == 0 || in == 0) { error = 48; return; } //the given data is empty
    readPngHeader(&in[0], size); if(error) return;
    size_t pos = 33; //first byte of the first chunk after the header
    std::vector<unsigned char> idat; //the data from idat chunks
    bool IEND = false;
    //loop through the chunks, ignoring unknown chunks and stopping at IEND chunk.
    //IDAT data is put at the start of the in buffer
    while(!IEND)
    {
      //error: size of the in buffer too small to contain next chunk
      if(pos + 8 >= size) { error = 30; return; }
      size_t chunkLength = read32bitInt(&in[pos]); pos += 4;
      if(chunkLength > 2147483647) { error = 63; return; }
      //error: size of the in buffer too small to contain next chunk
      if(pos + chunkLength >= size) { error = 35; return; }
      //IDAT chunk, containing compressed image data
      if(in[pos + 0] == 'I' && in[pos + 1] == 'D' && in[pos + 2] == 'A' && in[pos + 3] == 'T')
      {
        idat.insert(idat.end(), &in[pos + 4], &in[pos + 4 + chunkLength]);
        pos += (4 + chunkLength);
      }
      else if(in[pos + 0] == 'I' && in[pos + 1] == 'E' && in[pos + 2] == 'N' && in[pos + 3] == 'D')
      {
          pos += 4;
          IEND = true;
      }
      else //it's not an implemented chunk type, so ignore it: skip over the data
      {
        pos += (chunkLength + 4); //skip 4 letters and uninterpreted data of unimplemented chunk
      }
      pos += 4; //step over CRC (which is ignored)
    }
    std::vector<unsigned char> out; //now the out buffer will be filled
    ExtractZlib zlib(zlibinfo); //decompress with the Zlib decompressor
    error = zlib.decompress(out, idat);
    if(error) return; //stop if the zlib decompressor returned an error
  }

  //read the information from the header and store it in the Info
  void readPngHeader(const unsigned char* in, size_t inlength)
  {
    if(inlength < 29) { error = 27; return; } //error: the data length is smaller than the length of the header
    if(in[0] != 137 || in[1] != 80 || in[2] != 78 || in[3] != 71
    || in[4] != 13 || in[5] != 10 || in[6] != 26 || in[7] != 10) { error = 28; return; } //no PNG signature
    //error: it doesn't start with a IHDR chunk!
    if(in[12] != 'I' || in[13] != 'H' || in[14] != 'D' || in[15] != 'R') { error = 29; return; }
  }

  unsigned long readBitFromReversedStream(size_t& bitp, const unsigned char* bits)
  {
    unsigned long result = (bits[bitp >> 3] >> (7 - (bitp & 0x7))) & 1;
    bitp++;
    return result;
  }

  unsigned long readBitsFromReversedStream(size_t& bitp, const unsigned char* bits, unsigned long nbits)
  {
    unsigned long result = 0;
    for(size_t i = nbits - 1; i < nbits; i--) result += ((readBitFromReversedStream(bitp, bits)) << i);
    return result;
  }

  void setBitOfReversedStream(size_t& bitp, unsigned char* bits, unsigned long bit)
  {
    bits[bitp >> 3] |=  (bit << (7 - (bitp & 0x7))); bitp++;
  }

  unsigned long read32bitInt(const unsigned char* buffer)
  {
    return (unsigned int)((buffer[0] << 24u) | (buffer[1] << 16u) | (buffer[2] << 8u) | buffer[3]);
  }
};

void extractZlibInfo(std::vector<ZlibBlockInfo>& zlibinfo, const std::vector<unsigned char>& in)
{
  ExtractPNG decoder(&zlibinfo);
  decoder.decode(&in[0], in.size());

  if(decoder.error) std::cout << "extract error: " << decoder.error << std::endl;
}

} // namespace lodepng
