// Copyright 2013 Google Inc. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// Author: lode.vandevenne@gmail.com (Lode Vandevenne)
// Author: jyrki.alakuijala@gmail.com (Jyrki Alakuijala)

// Library to recompress and optimize PNG images. Uses Zopfli as the compression
// backend, chooses optimal PNG color model, and tries out several PNG filter
// strategies.

#ifndef ZOPFLIPNG_LIB_H_
#define ZOPFLIPNG_LIB_H_

#include <string>
#include <vector>

enum ZopfliPNGFilterStrategy {
  kStrategyZero = 0,
  kStrategyOne = 1,
  kStrategyTwo = 2,
  kStrategyThree = 3,
  kStrategyFour = 4,
  kStrategyMinSum,
  kStrategyEntropy,
  kStrategyPredefined,
  kStrategyBruteForce,
  kNumFilterStrategies /* Not a strategy but used for the size of this enum */
};

struct ZopfliPNGOptions {
  ZopfliPNGOptions();

  // Allow altering hidden colors of fully transparent pixels
  bool lossy_transparent;
  // Convert 16-bit per channel images to 8-bit per channel
  bool lossy_8bit;

  // Filter strategies to try
  std::vector<ZopfliPNGFilterStrategy> filter_strategies;

  // Automatically choose filter strategy using less good compression
  bool auto_filter_strategy;

  // PNG chunks to keep
  // chunks to literally copy over from the original PNG to the resulting one
  std::vector<std::string> keepchunks;

  // Use Zopfli deflate compression
  bool use_zopfli;

  // Zopfli number of iterations
  int num_iterations;

  // Zopfli number of iterations on large images
  int num_iterations_large;

  // 0=none, 1=first, 2=last, 3=both
  int block_split_strategy;
};

// Returns 0 on success, error code otherwise.
// If verbose is true, it will print some info while working.
int ZopfliPNGOptimize(const std::vector<unsigned char>& origpng,
    const ZopfliPNGOptions& png_options,
    bool verbose,
    std::vector<unsigned char>* resultpng);

#endif  // ZOPFLIPNG_LIB_H_
