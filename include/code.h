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
#include "xoroshiro128plus.h"

void transpose_columns(code_t *H);
void transpose_rows(code_t *H);

void compute_codeword(cw_t codeword, code_t *H, msg_t message);
void compute_syndrome(syndrome_t *syndrome, code_t *H, e_t *e_dense);
void compute_counters(counters_t counters, bit_t *syndrome, code_t *H);

void error_sparse_to_dense(e_t *e_dense, const sparse_t e_sparse,
                           index_t weight);
void syndrome_add_sparse_error(syndrome_t *syndrome, const sparse_t e_sparse,
                               index_t weight);

void generate_random_message(msg_t message, prng_t prng);
