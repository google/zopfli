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

#include "squeeze.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "blocksplitter.h"
#include "deflate.h"
#include "symbols.h"
#include "tree.h"
#include "util.h"

typedef struct SymbolStats {
  /* The literal and length symbols. */
  size_t litlens[ZOPFLI_NUM_LL];
  /* The 32 unique dist symbols, not the 32768 possible dists. */
  size_t dists[ZOPFLI_NUM_D];

  /* Length of each lit/len symbol in bits. */
  double ll_symbols[ZOPFLI_NUM_LL];
  /* Length of each dist symbol in bits. */
  double d_symbols[ZOPFLI_NUM_D];
} SymbolStats;

/* Sets everything to 0. */
static void InitStats(SymbolStats* stats) {
  memset(stats->litlens, 0, ZOPFLI_NUM_LL * sizeof(stats->litlens[0]));
  memset(stats->dists, 0, ZOPFLI_NUM_D * sizeof(stats->dists[0]));

  memset(stats->ll_symbols, 0, ZOPFLI_NUM_LL * sizeof(stats->ll_symbols[0]));
  memset(stats->d_symbols, 0, ZOPFLI_NUM_D * sizeof(stats->d_symbols[0]));
}

static void CopyStats(SymbolStats* source, SymbolStats* dest) {
  memcpy(dest->litlens, source->litlens,
         ZOPFLI_NUM_LL * sizeof(dest->litlens[0]));
  memcpy(dest->dists, source->dists, ZOPFLI_NUM_D * sizeof(dest->dists[0]));

  memcpy(dest->ll_symbols, source->ll_symbols,
         ZOPFLI_NUM_LL * sizeof(dest->ll_symbols[0]));
  memcpy(dest->d_symbols, source->d_symbols,
         ZOPFLI_NUM_D * sizeof(dest->d_symbols[0]));
}

/* Adds the bit lengths. */
static void AddWeighedStatFreqs(const SymbolStats* stats1, double w1,
                                const SymbolStats* stats2, double w2,
                                SymbolStats* result) {
  size_t i;
  for (i = 0; i < ZOPFLI_NUM_LL; i++) {
    result->litlens[i] =
        (size_t) (stats1->litlens[i] * w1 + stats2->litlens[i] * w2);
  }
  for (i = 0; i < ZOPFLI_NUM_D; i++) {
    result->dists[i] =
        (size_t) (stats1->dists[i] * w1 + stats2->dists[i] * w2);
  }
  result->litlens[256] = 1;  /* End symbol. */
}

typedef struct RanState {
  unsigned int m_w, m_z;
} RanState;

static void InitRanState(RanState* state) {
  state->m_w = 1;
  state->m_z = 2;
}

/* Get random number: "Multiply-With-Carry" generator of G. Marsaglia */
static unsigned int Ran(RanState* state) {
  state->m_z = 36969 * (state->m_z & 65535) + (state->m_z >> 16);
  state->m_w = 18000 * (state->m_w & 65535) + (state->m_w >> 16);
  return (state->m_z << 16) + state->m_w;  /* 32-bit result. */
}

static void RandomizeFreqs(RanState* state, size_t* freqs, int n) {
  int i;
  for (i = 0; i < n; i++) {
    if ((Ran(state) >> 4) % 3 == 0) freqs[i] = freqs[Ran(state) % n];
  }
}

static void RandomizeStatFreqs(RanState* state, SymbolStats* stats) {
  RandomizeFreqs(state, stats->litlens, ZOPFLI_NUM_LL);
  RandomizeFreqs(state, stats->dists, ZOPFLI_NUM_D);
  stats->litlens[256] = 1;  /* End symbol. */
}

static void ClearStatFreqs(SymbolStats* stats) {
  size_t i;
  for (i = 0; i < ZOPFLI_NUM_LL; i++) stats->litlens[i] = 0;
  for (i = 0; i < ZOPFLI_NUM_D; i++) stats->dists[i] = 0;
}

/*
Function that calculates a cost based on a model for the given LZ77 symbol.
litlen: means literal symbol if dist is 0, length otherwise.
*/
typedef double CostModelFun(unsigned litlen, unsigned dist, void* context);

/*
Cost model which should exactly match fixed tree.
type: CostModelFun
*/
static double GetCostFixed(unsigned litlen, unsigned dist, void* unused) {
  (void)unused;
  if (dist == 0) {
    if (litlen <= 143) return 8;
    else return 9;
  } else {
    int dbits = ZopfliGetDistExtraBits(dist);
    int lbits = ZopfliGetLengthExtraBits(litlen);
    int lsym = ZopfliGetLengthSymbol(litlen);
    int cost = 0;
    if (lsym <= 279) cost += 7;
    else cost += 8;
    cost += 5;  /* Every dist symbol has length 5. */
    return cost + dbits + lbits;
  }
}

/*
Finds the minimum possible cost this cost model can return for valid length and
distance symbols.
*/
static double GetCostModelMinCost(CostModelFun* costmodel, void* costcontext) {
  double mincost;
  int bestlength = 0; /* length that has lowest cost in the cost model */
  int bestdist = 0; /* distance that has lowest cost in the cost model */
  int i;
  /*
  Table of distances that have a different distance symbol in the deflate
  specification. Each value is the first distance that has a new symbol. Only
  different symbols affect the cost model so only these need to be checked.
  See RFC 1951 section 3.2.5. Compressed blocks (length and distance codes).
  */
  static const int dsymbols[30] = {
    1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513,
    769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577
  };

  mincost = ZOPFLI_LARGE_FLOAT;
  for (i = 3; i < 259; i++) {
    double c = costmodel(i, 1, costcontext);
    if (c < mincost) {
      bestlength = i;
      mincost = c;
    }
  }

  mincost = ZOPFLI_LARGE_FLOAT;
  for (i = 0; i < 30; i++) {
    double c = costmodel(3, dsymbols[i], costcontext);
    if (c < mincost) {
      bestdist = dsymbols[i];
      mincost = c;
    }
  }

  return costmodel(bestlength, bestdist, costcontext);
}

static size_t zopfli_min(size_t a, size_t b) {
  return a < b ? a : b;
}

/*
Performs the forward pass for "squeeze". Gets the most optimal length to reach
every byte from a previous byte, using cost calculations.
s: the ZopfliBlockState
in: the input data array
instart: where to start
inend: where to stop (not inclusive)
costmodel: function to calculate the cost of some lit/len/dist pair.
costcontext: abstract context for the costmodel function
length_array: output array of size (inend - instart) which will receive the best
    length to reach this byte from a previous byte.
returns the cost that was, according to the costmodel, needed to get to the end.
*/
static double GetBestLengths(ZopfliBlockState *s,
                             const unsigned char* in,
                             size_t instart, size_t inend,
                             CostModelFun* costmodel, void* costcontext,
                             unsigned short* length_array,
                             ZopfliHash* h, float* costs) {
  /* Best cost to get here so far. */
  size_t blocksize = inend - instart;
  size_t i = 0, k, kend;
  unsigned short leng;
  unsigned short dist;
  unsigned short sublen[259];
  size_t windowstart = instart > ZOPFLI_WINDOW_SIZE
      ? instart - ZOPFLI_WINDOW_SIZE : 0;
  double result;
  double mincost = GetCostModelMinCost(costmodel, costcontext);
  double mincostaddcostj;

  if (instart == inend) return 0;

  ZopfliResetHash(ZOPFLI_WINDOW_SIZE, h);
  ZopfliWarmupHash(in, windowstart, inend, h);
  for (i = windowstart; i < instart; i++) {
    ZopfliUpdateHash(in, i, inend, h);
  }

  for (i = 1; i < blocksize + 1; i++) costs[i] = ZOPFLI_LARGE_FLOAT;
  costs[0] = 0;  /* Because it's the start. */
  length_array[0] = 0;

  for (i = instart; i < inend; i++) {
    size_t j = i - instart;  /* Index in the costs array and length_array. */
    ZopfliUpdateHash(in, i, inend, h);

#ifdef ZOPFLI_SHORTCUT_LONG_REPETITIONS
    /* If we're in a long repetition of the same character and have more than
    ZOPFLI_MAX_MATCH characters before and after our position. */
    if (h->same[i & ZOPFLI_WINDOW_MASK] > ZOPFLI_MAX_MATCH * 2
        && i > instart + ZOPFLI_MAX_MATCH + 1
        && i + ZOPFLI_MAX_MATCH * 2 + 1 < inend
        && h->same[(i - ZOPFLI_MAX_MATCH) & ZOPFLI_WINDOW_MASK]
            > ZOPFLI_MAX_MATCH) {
      double symbolcost = costmodel(ZOPFLI_MAX_MATCH, 1, costcontext);
      /* Set the length to reach each one to ZOPFLI_MAX_MATCH, and the cost to
      the cost corresponding to that length. Doing this, we skip
      ZOPFLI_MAX_MATCH values to avoid calling ZopfliFindLongestMatch. */
      for (k = 0; k < ZOPFLI_MAX_MATCH; k++) {
        costs[j + ZOPFLI_MAX_MATCH] = costs[j] + symbolcost;
        length_array[j + ZOPFLI_MAX_MATCH] = ZOPFLI_MAX_MATCH;
        i++;
        j++;
        ZopfliUpdateHash(in, i, inend, h);
      }
    }
#endif

    ZopfliFindLongestMatch(s, h, in, i, inend, ZOPFLI_MAX_MATCH, sublen,
                           &dist, &leng);

    /* Literal. */
    if (i + 1 <= inend) {
      double newCost = costmodel(in[i], 0, costcontext) + costs[j];
      assert(newCost >= 0);
      if (newCost < costs[j + 1]) {
        costs[j + 1] = newCost;
        length_array[j + 1] = 1;
      }
    }
    /* Lengths. */
    kend = zopfli_min(leng, inend-i);
    mincostaddcostj = mincost + costs[j];
    for (k = 3; k <= kend; k++) {
      double newCost;

      /* Calling the cost model is expensive, avoid this if we are already at
      the minimum possible cost that it can return. */
     if (costs[j + k] <= mincostaddcostj) continue;

      newCost = costmodel(k, sublen[k], costcontext) + costs[j];
      assert(newCost >= 0);
      if (newCost < costs[j + k]) {
        assert(k <= ZOPFLI_MAX_MATCH);
        costs[j + k] = newCost;
        length_array[j + k] = k;
      }
    }
  }

  assert(costs[blocksize] >= 0);
  result = costs[blocksize];

  return result;
}

/*
Calculates the optimal path of lz77 lengths to use, from the calculated
length_array. The length_array must contain the optimal length to reach that
byte. The path will be filled with the lengths to use, so its data size will be
the amount of lz77 symbols.
*/
static void TraceBackwards(size_t size, const unsigned short* length_array,
                           unsigned short** path, size_t* pathsize) {
  size_t index = size;
  if (size == 0) return;
  for (;;) {
    ZOPFLI_APPEND_DATA(length_array[index], path, pathsize);
    assert(length_array[index] <= index);
    assert(length_array[index] <= ZOPFLI_MAX_MATCH);
    assert(length_array[index] != 0);
    index -= length_array[index];
    if (index == 0) break;
  }

  /* Mirror result. */
  for (index = 0; index < *pathsize / 2; index++) {
    unsigned short temp = (*path)[index];
    (*path)[index] = (*path)[*pathsize - index - 1];
    (*path)[*pathsize - index - 1] = temp;
  }
}

static void FollowPath(ZopfliBlockState* s,
                       const unsigned char* in, size_t instart, size_t inend,
                       unsigned short* path, size_t pathsize,
                       ZopfliLZ77Store* store, ZopfliHash *h) {
  size_t i, j, pos = 0;
  size_t windowstart = instart > ZOPFLI_WINDOW_SIZE
      ? instart - ZOPFLI_WINDOW_SIZE : 0;

  size_t total_length_test = 0;

  if (instart == inend) return;

  ZopfliResetHash(ZOPFLI_WINDOW_SIZE, h);
  ZopfliWarmupHash(in, windowstart, inend, h);
  for (i = windowstart; i < instart; i++) {
    ZopfliUpdateHash(in, i, inend, h);
  }

  pos = instart;
  for (i = 0; i < pathsize; i++) {
    unsigned short length = path[i];
    unsigned short dummy_length;
    unsigned short dist;
    assert(pos < inend);

    ZopfliUpdateHash(in, pos, inend, h);

    /* Add to output. */
    if (length >= ZOPFLI_MIN_MATCH) {
      /* Get the distance by recalculating longest match. The found length
      should match the length from the path. */
      ZopfliFindLongestMatch(s, h, in, pos, inend, length, 0,
                             &dist, &dummy_length);
      assert(!(dummy_length != length && length > 2 && dummy_length > 2));
      ZopfliVerifyLenDist(in, inend, pos, dist, length);
      ZopfliStoreLitLenDist(length, dist, pos, store);
      total_length_test += length;
    } else {
      length = 1;
      ZopfliStoreLitLenDist(in[pos], 0, pos, store);
      total_length_test++;
    }


    assert(pos + length <= inend);
    for (j = 1; j < length; j++) {
      ZopfliUpdateHash(in, pos + j, inend, h);
    }

    pos += length;
  }
}

/* Calculates the entropy of the statistics */
static void CalculateStatistics(SymbolStats* stats) {
  ZopfliCalculateEntropy(stats->litlens, ZOPFLI_NUM_LL, stats->ll_symbols);
  ZopfliCalculateEntropy(stats->dists, ZOPFLI_NUM_D, stats->d_symbols);
}

/* Appends the symbol statistics from the store. */
static void GetStatistics(const ZopfliLZ77Store* store, SymbolStats* stats) {
  size_t i;
  for (i = 0; i < store->size; i++) {
    if (store->dists[i] == 0) {
      stats->litlens[store->litlens[i]]++;
    } else {
      stats->litlens[ZopfliGetLengthSymbol(store->litlens[i])]++;
      stats->dists[ZopfliGetDistSymbol(store->dists[i])]++;
    }
  }
  stats->litlens[256] = 1;  /* End symbol. */

  CalculateStatistics(stats);
}

/*
Does a single run for ZopfliLZ77Optimal. For good compression, repeated runs
with updated statistics should be performed.
s: the block state
in: the input data array
instart: where to start
inend: where to stop (not inclusive)
path: pointer to dynamically allocated memory to store the path
pathsize: pointer to the size of the dynamic path array
length_array: array of size (inend - instart) used to store lengths
costmodel: function to use as the cost model for this squeeze run
costcontext: abstract context for the costmodel function
store: place to output the LZ77 data
returns the cost that was, according to the costmodel, needed to get to the end.
    This is not the actual cost.
*/
static double LZ77OptimalRun(ZopfliBlockState* s,
    const unsigned char* in, size_t instart, size_t inend,
    unsigned short** path, size_t* pathsize,
    unsigned short* length_array, CostModelFun* costmodel,
    void* costcontext, ZopfliLZ77Store* store,
    ZopfliHash* h, float* costs) {
  double cost = GetBestLengths(s, in, instart, inend, costmodel,
                costcontext, length_array, h, costs);
  free(*path);
  *path = 0;
  *pathsize = 0;
  TraceBackwards(inend - instart, length_array, path, pathsize);
  FollowPath(s, in, instart, inend, *path, *pathsize, store, h);
  assert(cost < ZOPFLI_LARGE_FLOAT);
  return cost;
}

/* ---- Distance-class aware optimal parsing ---- */

/*
Per block position, the longest match for each distance symbol, sorted by
length (descending). Packed to keep memory near that of the longest match
cache it replaces: one count byte per position, three bytes per entry.
*/
typedef struct ClassCache {
  unsigned char* cnt;     /* entries at this position; bit 7: long run */
  unsigned char* len;     /* match length - ZOPFLI_MIN_MATCH */
  unsigned short* dist;
  size_t n, cap;
} ClassCache;

#define CLASS_CACHE_RUN 0x80

static void ClassCachePush(ClassCache* c, unsigned short l, unsigned short d) {
  if (c->n == c->cap) {
    c->cap = c->cap + c->cap / 2 + 1024;
    c->len = (unsigned char*)realloc(c->len, c->cap);
    c->dist = (unsigned short*)realloc(c->dist, c->cap * sizeof(*c->dist));
    if (!c->len || !c->dist) exit(-1); /* Allocation failed. */
  }
  c->len[c->n] = (unsigned char)(l - ZOPFLI_MIN_MATCH);
  c->dist[c->n] = d;
  c->n++;
}

static void BuildClassCache(ZopfliBlockState* s,
                            const unsigned char* in, size_t instart,
                            size_t inend, ZopfliHash* h, ClassCache* c) {
  unsigned short sublen[259];
  size_t blocksize = inend - instart, i;
  size_t windowstart = instart > ZOPFLI_WINDOW_SIZE
      ? instart - ZOPFLI_WINDOW_SIZE : 0;
  unsigned short clen[30], cdist[30];
  memset(c, 0, sizeof(*c));
  c->cnt = (unsigned char*)malloc(blocksize + 1);
  if (!c->cnt) exit(-1); /* Allocation failed. */
  ZopfliResetHash(ZOPFLI_WINDOW_SIZE, h);
  ZopfliWarmupHash(in, windowstart, inend, h);
  for (i = windowstart; i < instart; i++) ZopfliUpdateHash(in, i, inend, h);
  for (i = instart; i < inend; i++) {
    int k, m = 0, run = 0;
    unsigned char order[30];
    ZopfliUpdateHash(in, i, inend, h);
    /* Same condition as ZOPFLI_SHORTCUT_LONG_REPETITIONS in GetBestLengths. */
    if (h->same[i & ZOPFLI_WINDOW_MASK] > ZOPFLI_MAX_MATCH * 2
        && i > instart + ZOPFLI_MAX_MATCH + 1
        && i + ZOPFLI_MAX_MATCH * 2 + 1 < inend
        && h->same[(i - ZOPFLI_MAX_MATCH) & ZOPFLI_WINDOW_MASK]
            > ZOPFLI_MAX_MATCH) {
      run = 1;
      for (k = 0; k < 30; k++) { clen[k] = 0; cdist[k] = 0; }
    } else {
      ZopfliFindClassMatches(h, in, i, inend, clen, cdist);
    }
    if (inend - i >= ZOPFLI_MIN_MATCH) {
      /* Deep search (cached): nearest distance for every length. */
      unsigned short d0, l0;
      ZopfliFindLongestMatch(s, h, in, i, inend, ZOPFLI_MAX_MATCH, sublen,
                             &d0, &l0);
      for (k = 3; k <= l0; k++) {
        int ds = ZopfliGetDistSymbol(sublen[k]);
        if (k > clen[ds]) { clen[ds] = k; cdist[ds] = sublen[k]; }
      }
    }
    for (k = 0; k < 30; k++) if (clen[k]) order[m++] = (unsigned char)k;
    /* insertion sort by length descending */
    for (k = 1; k < m; k++) {
      unsigned char t = order[k]; int q = k - 1;
      while (q >= 0 && clen[order[q]] < clen[t]) { order[q + 1] = order[q]; q--; }
      order[q + 1] = t;
    }
    for (k = 0; k < m; k++) {
      ClassCachePush(c, clen[order[k]], cdist[order[k]]);
    }
    c->cnt[i - instart] = (unsigned char)(m | (run ? CLASS_CACHE_RUN : 0));
  }
}

/* Longest match at a position from the class cache, nearest distance first,
as ZopfliFindLongestMatch would return it. */
static void ClassCacheLongest(const ClassCache* c, size_t e, size_t n,
                              unsigned short* leng, unsigned short* dist) {
  size_t k;
  *leng = 0;
  *dist = 0;
  if (n == 0) return;
  *leng = c->len[e] + ZOPFLI_MIN_MATCH;
  *dist = c->dist[e];
  for (k = 1; k < n && c->len[e + k] == c->len[e]; k++) {
    if (c->dist[e + k] < *dist) *dist = c->dist[e + k];
  }
}

static int GreedyLengthScore(int length, int distance) {
  return distance > 1024 ? length - 1 : length;
}

/*
The same lazy greedy parse as ZopfliLZ77Greedy, but reading matches from the
class cache instead of searching the hash chains again.
*/
static void GreedyFromClassCache(const unsigned char* in, size_t instart,
                                 size_t inend, const ClassCache* c,
                                 ZopfliLZ77Store* store) {
  size_t blocksize = inend - instart, i, j;
  size_t* eoff = (size_t*)malloc(sizeof(size_t) * (blocksize + 1));
  unsigned short leng, dist;
  unsigned prev_length = 0, prev_match = 0;
  int lengthscore, prevlengthscore, match_available = 0;
  if (!eoff) exit(-1); /* Allocation failed. */
  eoff[0] = 0;
  for (i = 0; i < blocksize; i++) {
    eoff[i + 1] = eoff[i] + (c->cnt[i] & ~CLASS_CACHE_RUN);
  }
  for (i = 0; i < blocksize; i++) {
    size_t pos = instart + i;
    ClassCacheLongest(c, eoff[i], eoff[i + 1] - eoff[i], &leng, &dist);
    lengthscore = GreedyLengthScore(leng, dist);
    prevlengthscore = GreedyLengthScore(prev_length, prev_match);
    if (match_available) {
      match_available = 0;
      if (lengthscore > prevlengthscore + 1) {
        ZopfliStoreLitLenDist(in[pos - 1], 0, pos - 1, store);
        if (lengthscore >= ZOPFLI_MIN_MATCH && leng < ZOPFLI_MAX_MATCH) {
          match_available = 1;
          prev_length = leng;
          prev_match = dist;
          continue;
        }
      } else {
        leng = prev_length;
        dist = prev_match;
        ZopfliVerifyLenDist(in, inend, pos - 1, dist, leng);
        ZopfliStoreLitLenDist(leng, dist, pos - 1, store);
        for (j = 2; j < leng; j++) i++;
        continue;
      }
    } else if (lengthscore >= ZOPFLI_MIN_MATCH && leng < ZOPFLI_MAX_MATCH) {
      match_available = 1;
      prev_length = leng;
      prev_match = dist;
      continue;
    }
    if (lengthscore >= ZOPFLI_MIN_MATCH) {
      ZopfliVerifyLenDist(in, inend, pos, dist, leng);
      ZopfliStoreLitLenDist(leng, dist, pos, store);
    } else {
      leng = 1;
      ZopfliStoreLitLenDist(in[pos], 0, pos, store);
    }
    for (j = 1; j < leng; j++) i++;
  }
  free(eoff);
}

static void CleanClassCache(ClassCache* c) {
  free(c->cnt); free(c->len); free(c->dist);
}

/* Forward DP using per-class matches. Fills length_array and dist_array. */
static void GetBestLengthsClass(const unsigned char* in, size_t instart,
                                size_t inend, const SymbolStats* st,
                                const ClassCache* c,
                                unsigned short* length_array,
                                unsigned short* dist_array, float* costs) {
  size_t blocksize = inend - instart, j;
  double lcost[ZOPFLI_MAX_MATCH + 1];
  double dcost[30];
  double minl = ZOPFLI_LARGE_FLOAT, mind = ZOPFLI_LARGE_FLOAT, mincost;
  size_t e0 = 0;
  int k;
  for (k = 3; k <= ZOPFLI_MAX_MATCH; k++) {
    lcost[k] = st->ll_symbols[ZopfliGetLengthSymbol(k)]
        + ZopfliGetLengthExtraBits(k);
    if (lcost[k] < minl) minl = lcost[k];
  }
  for (k = 0; k < 30; k++) {
    dcost[k] = st->d_symbols[k] + ZopfliGetDistSymbolExtraBits(k);
    if (dcost[k] < mind) mind = dcost[k];
  }
  mincost = minl + mind;
  for (j = 1; j <= blocksize; j++) costs[j] = ZOPFLI_LARGE_FLOAT;
  costs[0] = 0;
  length_array[0] = 0;
  for (j = 0; j < blocksize; j++) {
    double cj;
    double newCost;
    size_t e, eend;
    if (c->cnt[j] & CLASS_CACHE_RUN) {
      /* Mirror ZOPFLI_SHORTCUT_LONG_REPETITIONS: step by max-length matches
      of distance 1 instead of trying every length. Entries of the skipped
      positions are passed over too. */
      double symbolcost = lcost[ZOPFLI_MAX_MATCH] + dcost[0];
      size_t q;
      for (q = 0; q < ZOPFLI_MAX_MATCH && j + ZOPFLI_MAX_MATCH <= blocksize; q++) {
        costs[j + ZOPFLI_MAX_MATCH] = costs[j] + symbolcost;
        length_array[j + ZOPFLI_MAX_MATCH] = ZOPFLI_MAX_MATCH;
        dist_array[j + ZOPFLI_MAX_MATCH] = 1;
        e0 += c->cnt[j] & ~CLASS_CACHE_RUN;
        j++;
      }
      if (j >= blocksize) break;
    }
    e = e0;
    eend = e0 + (c->cnt[j] & ~CLASS_CACHE_RUN);
    e0 = eend;
    cj = costs[j];
    newCost = cj + st->ll_symbols[in[instart + j]];
    if (newCost < costs[j + 1]) {
      costs[j + 1] = newCost;
      length_array[j + 1] = 1;
      dist_array[j + 1] = 0;
    }
    if (e < eend) {
      int kmax = c->len[e] + ZOPFLI_MIN_MATCH;
      double best = ZOPFLI_LARGE_FLOAT;
      unsigned short bestd = 0;
      double floorc = cj + mincost;
      if ((size_t)kmax > blocksize - j) kmax = (int)(blocksize - j);
      for (k = kmax; k >= ZOPFLI_MIN_MATCH; k--) {
        while (e < eend && c->len[e] + ZOPFLI_MIN_MATCH >= k) {
          double dc = dcost[ZopfliGetDistSymbol(c->dist[e])];
          if (dc < best) { best = dc; bestd = c->dist[e]; }
          e++;
        }
        if (costs[j + k] <= floorc) continue;
        newCost = cj + best + lcost[k];
        if (newCost < costs[j + k]) {
          costs[j + k] = newCost;
          length_array[j + k] = k;
          dist_array[j + k] = bestd;
        }
      }
    }
  }
}

static void TraceAndStoreClass(const unsigned char* in, size_t instart,
                               size_t inend, const unsigned short* length_array,
                               const unsigned short* dist_array,
                               ZopfliLZ77Store* store) {
  size_t size = inend - instart, idx = size, n = 0, i;
  unsigned* steps;
  if (size == 0) return;
  steps = (unsigned*)malloc(sizeof(unsigned) * (size + 1));
  if (!steps) exit(-1);
  while (idx > 0) {
    steps[n++] = (unsigned)idx;
    idx -= length_array[idx];
  }
  for (i = n; i-- > 0;) {
    size_t end = steps[i];
    unsigned short l = length_array[end];
    size_t pos = instart + end - l;
    if (l >= ZOPFLI_MIN_MATCH) {
      ZopfliVerifyLenDist(in, inend, pos, dist_array[end], l);
      ZopfliStoreLitLenDist(l, dist_array[end], pos, store);
    } else {
      ZopfliStoreLitLenDist(in[pos], 0, pos, store);
    }
  }
  free(steps);
}

void ZopfliLZ77Optimal(ZopfliBlockState *s,
                       const unsigned char* in, size_t instart, size_t inend,
                       int numiterations,
                       ZopfliLZ77Store* store) {
  /* Dist to get to here with smallest cost. */
  size_t blocksize = inend - instart;
  unsigned short* length_array =
      (unsigned short*)malloc(sizeof(unsigned short) * (blocksize + 1));
  unsigned short* dist_array =
      (unsigned short*)malloc(sizeof(unsigned short) * (blocksize + 1));
  ClassCache ccache;
  ZopfliLZ77Store currentstore;
  ZopfliHash hash;
  ZopfliHash* h = &hash;
  SymbolStats stats, beststats, laststats;
  int i;
  float* costs = (float*)malloc(sizeof(float) * (blocksize + 1));
  double cost;
  double bestcost = ZOPFLI_LARGE_FLOAT;
  double lastcost = 0;
  /* Try randomizing the costs a bit once the size stabilizes. */
  RanState ran_state;
  int lastrandomstep = -1;

  if (!costs) exit(-1); /* Allocation failed. */
  if (!length_array) exit(-1); /* Allocation failed. */
  if (!dist_array) exit(-1); /* Allocation failed. */

  InitRanState(&ran_state);
  InitStats(&stats);
  ZopfliInitLZ77Store(in, &currentstore);
  ZopfliAllocHash(ZOPFLI_WINDOW_SIZE, h);

  /* Do regular deflate, then loop multiple shortest path runs, each time using
  the statistics of the previous run. */

  /* Collect, once, the best match per distance symbol at every position. */
  BuildClassCache(s, in, instart, inend, h, &ccache);

  /* Initial run: the same lazy greedy parse as ZopfliLZ77Greedy, read from
  the cache instead of searching again. */
  GreedyFromClassCache(in, instart, inend, &ccache, &currentstore);
  GetStatistics(&currentstore, &stats);

  /* Repeat statistics with each time the cost model from the previous stat
  run. */
  for (i = 0; i < numiterations; i++) {
    ZopfliCleanLZ77Store(&currentstore);
    ZopfliInitLZ77Store(in, &currentstore);
    GetBestLengthsClass(in, instart, inend, &stats, &ccache,
                        length_array, dist_array, costs);
    TraceAndStoreClass(in, instart, inend, length_array, dist_array,
                       &currentstore);
    cost = ZopfliCalculateBlockSize(&currentstore, 0, currentstore.size, 2);
    if (s->options->verbose_more || (s->options->verbose && cost < bestcost)) {
      fprintf(stderr, "Iteration %d: %d bit\n", i, (int) cost);
    }
    if (cost < bestcost) {
      /* Copy to the output store. */
      ZopfliCopyLZ77Store(&currentstore, store);
      CopyStats(&stats, &beststats);
      bestcost = cost;
    }
    CopyStats(&stats, &laststats);
    ClearStatFreqs(&stats);
    GetStatistics(&currentstore, &stats);
    if (lastrandomstep != -1) {
      /* This makes it converge slower but better. Do it only once the
      randomness kicks in so that if the user does few iterations, it gives a
      better result sooner. */
      AddWeighedStatFreqs(&stats, 1.0, &laststats, 0.5, &stats);
      CalculateStatistics(&stats);
    }
    if (i > 5 && cost == lastcost) {
      CopyStats(&beststats, &stats);
      RandomizeStatFreqs(&ran_state, &stats);
      CalculateStatistics(&stats);
      lastrandomstep = i;
    }
    lastcost = cost;
  }

  CleanClassCache(&ccache);
  free(length_array);
  free(dist_array);
  free(costs);
  ZopfliCleanLZ77Store(&currentstore);
  ZopfliCleanHash(h);
}

void ZopfliLZ77OptimalFixed(ZopfliBlockState *s,
                            const unsigned char* in,
                            size_t instart, size_t inend,
                            ZopfliLZ77Store* store)
{
  /* Dist to get to here with smallest cost. */
  size_t blocksize = inend - instart;
  unsigned short* length_array =
      (unsigned short*)malloc(sizeof(unsigned short) * (blocksize + 1));
  unsigned short* path = 0;
  size_t pathsize = 0;
  ZopfliHash hash;
  ZopfliHash* h = &hash;
  float* costs = (float*)malloc(sizeof(float) * (blocksize + 1));

  if (!costs) exit(-1); /* Allocation failed. */
  if (!length_array) exit(-1); /* Allocation failed. */

  ZopfliAllocHash(ZOPFLI_WINDOW_SIZE, h);

  s->blockstart = instart;
  s->blockend = inend;

  /* Shortest path for fixed tree This one should give the shortest possible
  result for fixed tree, no repeated runs are needed since the tree is known. */
  LZ77OptimalRun(s, in, instart, inend, &path, &pathsize,
                 length_array, GetCostFixed, 0, store, h, costs);

  free(length_array);
  free(path);
  free(costs);
  ZopfliCleanHash(h);
}
