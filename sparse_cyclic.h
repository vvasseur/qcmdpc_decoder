/*
   Copyright (c) 2019 Valentin Vasseur

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to
   deal in the Software without restriction, including without limitation the
   rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
   sell copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
   IN THE SOFTWARE
*/
#ifndef SPARSE_CYCLIC_H
#define SPARSE_CYCLIC_H
#include "types.h"
#include "xoroshiro128plus.h"

index_t insert_sorted(sparse_t array, index_t value, index_t max_i);
void insert_sorted_noinc(sparse_t array, index_t value, index_t max_i);

sparse_t sparse_new(index_t weight);
void sparse_free(sparse_t array);
void sparse_rand(sparse_t array, index_t length, index_t weight, prng_t prng);

void transpose(sparse_t dst, const sparse_t src, index_t weight,
               index_t length);

void multiply_xor_mod2(dense_t z, const sparse_t x, const dense_t y,
                       index_t block_weight, index_t block_length);
void multiply_add(dense_t z, const sparse_t x, const dense_t y,
                  index_t block_weight, index_t block_length);
#ifdef AVX
void multiply_xor_mod2_avx2(dense_t z, const sparse_t x, const dense_t y,
                            index_t block_weight, index_t block_length);
void multiply_avx2(dense_t z, const sparse_t x, const dense_t y,
                   index_t block_weight, index_t block_length);
#endif
#endif
