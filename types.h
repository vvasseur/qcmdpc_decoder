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
#ifndef TYPES_H
#define TYPES_H
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

typedef int_fast32_t index_t;
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

/* State of the decoder */
struct decoder {
    index_t Hcolumns[INDEX][BLOCK_WEIGHT];
    index_t Hrows[INDEX][BLOCK_WEIGHT];
    bit_t bits[INDEX][BLOCK_LENGTH];
    bit_t syndrome[2 * SIZE_AVX] __attribute__((aligned(32)));
    bit_t e[INDEX][2 * SIZE_AVX] __attribute__((aligned(32)));
    bit_t counters[INDEX][2 * SIZE_AVX] __attribute__((aligned(32)));
    index_t syndrome_weight;
    index_t error_weight;
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
    index_t Hcolumns[INDEX][BLOCK_WEIGHT];
    index_t Hrows[INDEX][BLOCK_WEIGHT];
    bit_t message[2 * SIZE_AVX] __attribute__((aligned(32)));
    bit_t codeword[INDEX][2 * SIZE_AVX] __attribute__((aligned(32)));
    bit_t bits[INDEX][2 * SIZE_AVX] __attribute__((aligned(32)));
    bit_t syndrome[2 * SIZE_AVX] __attribute__((aligned(32)));
    llr_t r[INDEX][BLOCK_LENGTH];
    llr_t v_to_c[INDEX][BLOCK_WEIGHT][BLOCK_LENGTH];
    llr_t c_to_v[INDEX][BLOCK_WEIGHT][BLOCK_LENGTH];
    llr_t tree[1 << LOG2(INDEX * BLOCK_LENGTH)];
    index_t iter;
};
#endif
