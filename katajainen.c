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

#include "katajainen.h"
#include <assert.h>
#include <stdlib.h>

typedef struct Node Node;

/*
Nodes forming chains. Also used to represent leaves.
*/
struct Node {
  size_t weight;  /* Total weight (symbol count) of this chain. */
  Node* tail;  /* Previous node(s) of this chain, or 0 if none. */
  int count;  /* Leaf symbol index, or number of leaves before this chain. */
  char inuse;  /* Tracking for garbage collection. */
};

/*
Memory pool for nodes.
*/
typedef struct NodePool {
  Node* nodes;  /* The pool. */
  Node* next;  /* Pointer to a possibly free node in the pool. */
  int size;  /* Size of the memory pool. */
} NodePool;

/*
Initializes a chain node with the given values and marks it as in use.
*/
static void InitNode(size_t weight, int count, Node* tail, Node* node) {
  node->weight = weight;
  node->count = count;
  node->tail = tail;
  node->inuse = 1;
}

/*
Finds a free location in the memory pool. Performs garbage collection if needed.
lists: If given, used to mark in-use nodes during garbage collection.
maxbits: Size of lists.
pool: Memory pool to get free node from.
*/
static Node* GetFreeNode(Node* (*lists)[2], int maxbits, NodePool* pool) {
  for (;;) {
    if (pool->next >= &pool->nodes[pool->size]) {
      /* Garbage collection. */
      int i;
      for (i = 0; i < pool->size; i++) {
        pool->nodes[i].inuse = 0;
      }
      if (lists) {
        for (i = 0; i < maxbits * 2; i++) {
          Node* node;
          for (node = lists[i / 2][i % 2]; node; node = node->tail) {
            node->inuse = 1;
          }
        }
      }
      pool->next = &pool->nodes[0];
    }
    if (!pool->next->inuse) break;  /* Found one. */
    pool->next++;
  }
  return pool->next++;
}


/*
Performs a Boundary Package-Merge step. Puts a new chain in the given list. The
new chain is, depending on the weights, a leaf or a combination of two chains
from the previous list.
lists: The lists of chains.
maxbits: Number of lists.
leaves: The leaves, one per symbol.
numsymbols: Number of leaves.
pool: the node memory pool.
index: The index of the list in which a new chain or leaf is required.
final: Whether this is the last time this function is called. If it is then it
  is no more needed to recursively call self.
*/
static void BoundaryPM(Node* (*lists)[2], int maxbits,
    Node* leaves, int numsymbols, NodePool* pool, int index, char final) {
  Node* newchain;
  Node* oldchain;
  int lastcount = lists[index][1]->count;  /* Count of last chain of list. */

  if (index == 0 && lastcount >= numsymbols) return;

  newchain = GetFreeNode(lists, maxbits, pool);
  oldchain = lists[index][1];

  /* These are set up before the recursive calls below, so that there is a list
  pointing to the new node, to let the garbage collection know it's in use. */
  lists[index][0] = oldchain;
  lists[index][1] = newchain;

  if (index == 0) {
    /* New leaf node in list 0. */
    InitNode(leaves[lastcount].weight, lastcount + 1, 0, newchain);
  } else {
    size_t sum = lists[index - 1][0]->weight + lists[index - 1][1]->weight;
    if (lastcount < numsymbols && sum > leaves[lastcount].weight) {
      /* New leaf inserted in list, so count is incremented. */
      InitNode(leaves[lastcount].weight, lastcount + 1, oldchain->tail,
          newchain);
    } else {
      InitNode(sum, lastcount, lists[index - 1][1], newchain);
      if (!final) {
        /* Two lookahead chains of previous list used up, create new ones. */
        BoundaryPM(lists, maxbits, leaves, numsymbols, pool, index - 1, 0);
        BoundaryPM(lists, maxbits, leaves, numsymbols, pool, index - 1, 0);
      }
    }
  }
}

/*
Initializes each list with as lookahead chains the two leaves with lowest
weights.
*/
static void InitLists(
    NodePool* pool, const Node* leaves, int maxbits, Node* (*lists)[2]) {
  int i;
  Node* node0 = GetFreeNode(0, maxbits, pool);
  Node* node1 = GetFreeNode(0, maxbits, pool);
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
  Node* node;
  for (node = chain; node; node = node->tail) {
    int i;
    for (i = 0; i < node->count; i++) {
      bitlengths[leaves[i].count]++;
    }
  }
}

/*
Comparator for sorting the leaves. Has the function signature for qsort.
*/
static int LeafComparator(const void* a, const void* b) {
  return ((const Node*)a)->weight - ((const Node*)b)->weight;
}

int ZopfliLengthLimitedCodeLengths(
    const size_t* frequencies, int n, int maxbits, unsigned* bitlengths) {
  NodePool pool;
  int i;
  int numsymbols = 0;  /* Amount of symbols with frequency > 0. */
  int numBoundaryPMRuns;

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

  /* Sort the leaves from lightest to heaviest. */
  qsort(leaves, numsymbols, sizeof(Node), LeafComparator);

  /* Initialize node memory pool. */
  pool.size = 2 * maxbits * (maxbits + 1);
  pool.nodes = (Node*)malloc(pool.size * sizeof(*pool.nodes));
  pool.next = pool.nodes;
  for (i = 0; i < pool.size; i++) {
    pool.nodes[i].inuse = 0;
  }

  lists = (Node* (*)[2])malloc(maxbits * sizeof(*lists));
  InitLists(&pool, leaves, maxbits, lists);

  /* In the last list, 2 * numsymbols - 2 active chains need to be created. Two
  are already created in the initialization. Each BoundaryPM run creates one. */
  numBoundaryPMRuns = 2 * numsymbols - 4;
  for (i = 0; i < numBoundaryPMRuns; i++) {
    char final = i == numBoundaryPMRuns - 1;
    BoundaryPM(lists, maxbits, leaves, numsymbols, &pool, maxbits - 1, final);
  }

  ExtractBitLengths(lists[maxbits - 1][1], leaves, bitlengths);

  free(lists);
  free(leaves);
  free(pool.nodes);
  return 0;  /* OK. */
}
