/*
   Copyright (c) 2023 Valentin Vasseur

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

#include <stdatomic.h>

typedef struct {
    int n_threads;
    int max_iter;
    atomic_int run;
    long int *n_test;
    long int *n_success;
    long int **n_iter;
} decoding_results_t;

void init_decoding_results(decoding_results_t *res, int n_threads,
                           int max_iter);
void clear_decoding_results(decoding_results_t *res);
void sum_decoding_results(long int *test_total, long int *success_total,
                          long int *iter_total, const decoding_results_t *res);
void decoder_loop(decoding_results_t *results, int n_threads, long int r);
void decoder_stop(decoding_results_t *res);
