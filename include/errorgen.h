/*
   Copyright (c) 2020-2021 Valentin Vasseur

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
#pragma once
#include "types.h"
#include "xoshiro256plusplus.h"

void generate_random_error(sparse_t e_block, index_t weight, prng_t prng);
void generate_random_syndrome_error(sparse_t e_block, index_t weight,
                                    prng_t prng);
void generate_around_word(sparse_t e_dst, index_t weight_dst, sparse_t e_src,
                          index_t weight_src, index_t intersections,
                          prng_t prng);
void generate_near_codeword(sparse_t e_block, code_t *H, prng_t prng);
void generate_near_codeword2(sparse_t e_block, code_t *H, prng_t prng);
void generate_codeword(sparse_t e_block, code_t *H, prng_t prng);
