/*
   Copyright (c) 2019-2021 Valentin Vasseur

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
#include <stdbool.h>
#include <stdint.h>

#include "log2.h"
#include "param.h"

/* Round relevant arrays size to the next multiple of 16 * 256 bits (to use the
 * 16 ymm AVX registers). */
#define ROUND_UP(N, S) ((((N) + (S)-1) / (S)) * (S))
#define SIZE_AVX (ROUND_UP(BLOCK_LENGTH * 8 * sizeof(bit_t), 256 * 16) / 8)

#define SWAP(X, Y)                                                             \
    do {                                                                       \
        __typeof__(X) _T = X;                                                  \
        X = Y;                                                                 \
        Y = _T;                                                                \
    } while (0)
#define SATURATE(X, SAT) ((X) > (SAT)) ? (SAT) : (((X) > (SAT)) ? -(SAT) : (X))

typedef int64_t index_t;
typedef index_t *sparse_t;
typedef float llr_t;

typedef uint8_t bit_t;
typedef bit_t *dense_t;

typedef struct array a_t;
typedef struct flip_list fl_t;
typedef struct decoder *decoder_t;
typedef struct decoder_bp *decoder_bp_t;

/* Double linked list to store previous flips */
struct flip_list {
    index_t first;
    uint8_t tod[INDEX * BLOCK_LENGTH];
    index_t prev[INDEX * BLOCK_LENGTH];
    index_t next[INDEX * BLOCK_LENGTH];
    index_t length;
};

/* Array to store black and gray positions */
struct array {
    uint8_t index[INDEX * BLOCK_LENGTH];
    index_t position[INDEX * BLOCK_LENGTH];
    index_t counter[INDEX * BLOCK_LENGTH];
    index_t length;
};

typedef struct {
    index_t index;
    index_t position;
    bit_t counter;
} pos_counter_t;

typedef struct {
    index_t columns[INDEX][BLOCK_WEIGHT];
    index_t rows[INDEX][BLOCK_WEIGHT];
} code_t;

typedef struct {
    bit_t vec[INDEX][2 * SIZE_AVX] __attribute__((aligned(32)));
    index_t weight;
} e_t;

typedef struct {
    bit_t vec[2 * SIZE_AVX] __attribute__((aligned(32)));
    index_t weight;
} syndrome_t;

typedef bit_t msg_t[2 * SIZE_AVX] __attribute__((aligned(32)));
typedef bit_t cw_t[INDEX][2 * SIZE_AVX] __attribute__((aligned(32)));
typedef bit_t bits_t[INDEX][BLOCK_LENGTH];
typedef bit_t counters_t[INDEX][2 * SIZE_AVX] __attribute__((aligned(32)));

/* State of the decoder */
struct decoder {
    code_t *H;
    syndrome_t *syndrome;
    e_t *e;
    bits_t bits;
    counters_t counters;
    index_t iter;
    bool blocked;
#if (ALGO == BACKFLIP) || (ALGO == BACKFLIP2)
    fl_t fl;
#endif
#if (ALGO == GRAY_BGF) || (ALGO == GRAY_BGB) || (ALGO == GRAY_B) ||            \
    (ALGO == GRAY_BG)
    a_t gray;
    a_t black;
#endif
#if (ALGO == SORT)
    pos_counter_t sorted_counters[GRAY_SIZE];
#endif
};

/* State of the decoder for belief propagation */
struct decoder_bp {
    code_t *H;
    syndrome_t *syndrome;
    e_t *e;
    e_t bits;
    index_t iter;
    msg_t message;
    cw_t codeword;
    llr_t r[INDEX][BLOCK_LENGTH];
    llr_t v_to_c[INDEX][BLOCK_WEIGHT][BLOCK_LENGTH];
    llr_t c_to_v[INDEX][BLOCK_WEIGHT][BLOCK_LENGTH];
    llr_t tree[1 << LOG2(INDEX * BLOCK_LENGTH)];
};
