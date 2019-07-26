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
Bounded package merge algorithm, based on the paper
"A Fast and Space-Economical Algorithm for Length-Limited Coding
Jyrki Katajainen, Alistair Moffat, Andrew Turpin".
*/

#ifdef __cplusplus
#include <algorithm>
extern "C" {
#endif

#include "katajainen.h"
#include <assert.h>
#include <stdlib.h>
#include <limits.h>

#ifdef __cplusplus
}
#endif

typedef struct Node Node;

/*
Nodes forming chains. Also used to represent leaves.
*/
struct Node {
  size_t weight;  /* Total weight (symbol count) of this chain. */
  Node* tail;  /* Previous node(s) of this chain, or 0 if none. */
  int count;  /* Leaf symbol index, or number of leaves before this chain. */
};

/*
Initializes a chain node with the given values and marks it as in use.
*/
static void InitNode(size_t weight, int count, Node* tail, Node* node) {
  node->weight = weight;
  node->count = count;
  node->tail = tail;
}

static void BoundaryPMFinal(Node* (*lists)[2],
    Node* leaves, int numsymbols, Node* pool, int index) {
  int lastcount = lists[index][1]->count;  /* Count of last chain of list. */

  size_t sum = lists[index - 1][0]->weight + lists[index - 1][1]->weight;

  if (lastcount < numsymbols && sum > leaves[lastcount].weight) {
    Node* oldchain = lists[index][1]->tail;

    lists[index][1] = pool;
    pool->count = lastcount + 1;
    pool->tail = oldchain;
  } else {
    lists[index][1]->tail = lists[index - 1][1];
  }
}

/*
Initializes each list with as lookahead chains the two leaves with lowest
weights.
*/
static void InitLists(
    Node* pool, const Node* leaves, int maxbits, Node* (*lists)[2]) {
  int i;
  Node* node0 = pool;
  Node* node1 = pool + 1;
  InitNode(leaves[0].weight, 1, 0, node0);
  InitNode(leaves[1].weight, 2, 0, node1);
  for (i = 0; i < maxbits; i++) {
    lists[i][0] = node0;
    lists[i][1] = node1;
  }
}

/*
Converts result of boundary package-merge to the bitlengths. The result in the
last chain of the last list contains the amount of active leaves in each list.
chain: Chain to extract the bit length from (last chain from last list).
*/
static void ExtractBitLengths(Node* chain, Node* leaves, unsigned* bitlengths) {
  int counts[16] = {0};
  unsigned end = 16;
  unsigned ptr = 15;
  unsigned value = 1;
  Node* node;
  int val;

  for (node = chain; node; node = node->tail) {
    counts[--end] = node->count;
  }

  val = counts[15];
  while (ptr >= end) {
    for (; val > counts[ptr - 1]; val--) {
      bitlengths[leaves[val - 1].count] = value;
    }
    ptr--;
    value++;
  }
}

#ifndef __cplusplus
/*
Comparator for sorting the leaves. Has the function signature for qsort.
*/
static int LeafComparator(const void* a, const void* b) {
  return ((const Node*)a)->weight - ((const Node*)b)->weight;
}
#else
struct {
  bool operator()(const Node a, const Node b) {
    return (a.weight < b.weight);
  }
} cmp;
#endif

#ifdef __cplusplus
extern "C"
#endif
int ZopfliLengthLimitedCodeLengths(
    const size_t* frequencies, int n, int maxbits, unsigned* bitlengths) {
  Node* pool;
  int i;
  int numsymbols = 0;  /* Amount of symbols with frequency > 0. */
  int numBoundaryPMRuns;
  Node* nodes;
  unsigned char stack[16];

  /* Array of lists of chains. Each list requires only two lookahead chains at
  a time, so each list is a array of two Node*'s. */
  Node* (*lists)[2];

  /* One leaf per symbol. Only numsymbols leaves will be used. */
  Node* leaves = (Node*)malloc(n * sizeof(*leaves));

  /* Initialize all bitlengths at 0. */
  for (i = 0; i < n; i++) {
    bitlengths[i] = 0;
  }

  /* Count used symbols and place them in the leaves. */
  for (i = 0; i < n; i++) {
    if (frequencies[i]) {
      leaves[numsymbols].weight = frequencies[i];
      leaves[numsymbols].count = i;  /* Index of symbol this leaf represents. */
      numsymbols++;
    }
  }

  /* Check special cases and error conditions. */
  if ((1 << maxbits) < numsymbols) {
    free(leaves);
    return 1;  /* Error, too few maxbits to represent symbols. */
  }
  if (numsymbols == 0) {
    free(leaves);
    return 0;  /* No symbols at all. OK. */
  }
  if (numsymbols == 1) {
    bitlengths[leaves[0].count] = 1;
    free(leaves);
    return 0;  /* Only one symbol, give it bitlength 1, not 0. OK. */
  }
  if (numsymbols == 2) {
    bitlengths[leaves[0].count]++;
    bitlengths[leaves[1].count]++;
    free(leaves);
    return 0;
  }

  /* Sort the leaves from lightest to heaviest. Add count into the same
  variable for stable sorting. */
  for (i = 0; i < numsymbols; i++) {
    if (leaves[i].weight >=
        ((size_t)1 << (sizeof(leaves[0].weight) * CHAR_BIT - 9))) {
      free(leaves);
      return 1;  /* Error, we need 9 bits for the count. */
    }
    leaves[i].weight = (leaves[i].weight << 9) | leaves[i].count;
  }
#ifdef __cplusplus
  std::sort(leaves, leaves + numsymbols, cmp);
#else
  qsort(leaves, numsymbols, sizeof(Node), LeafComparator);
#endif
  for (i = 0; i < numsymbols; i++) {
    leaves[i].weight >>= 9;
  }

  if (numsymbols - 1 < maxbits) {
    maxbits = numsymbols - 1;
  }

  /* Initialize node memory pool. */
  nodes = (Node*)malloc(maxbits * 2 * numsymbols * sizeof(Node));
  pool = nodes;

  lists = (Node* (*)[2])malloc(maxbits * sizeof(*lists));
  InitLists(pool, leaves, maxbits, lists);
  pool += 2;

  /* In the last list, 2 * numsymbols - 2 active chains need to be created. Two
  are already created in the initialization. Each BoundaryPM run creates one. */
  numBoundaryPMRuns = 2 * numsymbols - 4;
  for (i = 0; i < numBoundaryPMRuns - 1; i++) {
    /*
    Performs a Boundary Package-Merge step. Puts a new chain in the given list. The
    new chain is, depending on the weights, a leaf or a combination of two chains
    from the previous list.
    */
    unsigned stackpos;
    stack[0] = maxbits - 1;

    for (stackpos = 0; ;) {
      unsigned char index = stack[stackpos];

      int lastcount = lists[index][1]->count;  /* Count of last chain of list. */

      Node* newchain = pool++;
      Node* oldchain = lists[index][1];
      size_t sum;

      /* These are set up before the recursive calls below, so that there is a list
      pointing to the new node, to let the garbage collection know it's in use. */
      lists[index][0] = oldchain;
      lists[index][1] = newchain;

      sum = lists[index - 1][0]->weight + lists[index - 1][1]->weight;

      if (lastcount < numsymbols && sum > leaves[lastcount].weight) {
        /* New leaf inserted in list, so count is incremented. */
        InitNode(leaves[lastcount].weight, lastcount + 1, oldchain->tail, newchain);
      } else {
        InitNode(sum, lastcount, lists[index - 1][1], newchain);
        /* Two lookahead chains of previous list used up, create new ones. */
        if (index == 1) {
          if (lists[0][1]->count < numsymbols) {
            lastcount = lists[0][1]->count;
            lists[0][0] = lists[0][1];
            lists[0][1] = pool++;
            InitNode(leaves[lastcount].weight, lastcount + 1, 0, lists[0][1]);
            lastcount++;
            if(lastcount < numsymbols){
              lists[0][0] = lists[0][1];
              lists[0][1] = pool++;
              InitNode(leaves[lastcount].weight, lastcount + 1, 0, lists[0][1]);
            }
          }
        }
        else {
          stack[stackpos++] = index - 1;
          stack[stackpos++] = index - 1;
        }
      }
      if (!stackpos--) {
        break;
      }
    }
  }
  BoundaryPMFinal(lists, leaves, numsymbols, pool, maxbits - 1);

  ExtractBitLengths(lists[maxbits - 1][1], leaves, bitlengths);

  free(lists);
  free(leaves);
  free(nodes);
  return 0;  /* OK. */
}
