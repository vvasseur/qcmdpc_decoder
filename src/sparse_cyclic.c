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
#ifdef AVX
#include <immintrin.h>
#endif
#include <stdlib.h>

#include "sparse_cyclic.h"

sparse_t sparse_new(index_t weight) {
    sparse_t h = (sparse_t)malloc(weight * sizeof(index_t));

    return h;
}

void sparse_free(sparse_t array) { free(array); }

/* Insert in place. */
void insert_sorted_noinc(sparse_t array, index_t value, index_t max_i) {
    index_t i;
    for (i = 0; i < max_i && array[i] <= value; i++)
        ;
    for (index_t j = max_i; j > i; j--)
        array[j] = array[j - 1];
    array[i] = value;
}

/* Insert in place. */
index_t insert_sorted(sparse_t array, index_t value, index_t max_i) {
    index_t i;
    for (i = 0; i < max_i && array[i] <= value; i++, value++)
        ;
    for (index_t j = max_i; j > i; j--)
        array[j] = array[j - 1];
    array[i] = value;
    return value;
}

/* Pick a random (sparse) binary block h of weight 'weight' in a previously
 * allocated block. */
void sparse_rand(sparse_t array, index_t weight, index_t length, prng_t prng) {
    /* Get an ordered list of positions for which the bit should be set to 1. */
    for (index_t i = 0; i < weight; i++) {
        index_t rand = prng->random_lim(--length, &prng->s0, &prng->s1);
        insert_sorted(array, rand, i);
    }
}

void transpose(sparse_t dst, const sparse_t src, index_t weight,
               index_t length) {
    index_t l = 0;
    if (src[0] == 0) {
        dst[0] = 0;
        l = 1;
    }
    else {
        dst[0] = -src[weight - 1] + length;
    }
    for (index_t k = 1; k < weight; ++k) {
        dst[k] = -src[weight + l - 1 - k] + length;
    }
}

void transpose_columns(code_t *H) {
    for (index_t i = 0; i < INDEX; ++i) {
        transpose(H->rows[i], H->columns[i], BLOCK_WEIGHT, BLOCK_LENGTH);
    }
}

void transpose_rows(code_t *H) {
    for (index_t i = 0; i < INDEX; ++i) {
        transpose(H->columns[i], H->rows[i], BLOCK_WEIGHT, BLOCK_LENGTH);
    }
}

struct mult_t {
    dense_t y;
    dense_t z;
    index_t len;
};

/* Multiply modulo 2 the sparse vector 'x' of weight 'block_weight' by the dense
 * vector 'y' of length 'block_length' and xor the result in 'z'. */
void multiply_xor_mod2(dense_t restrict z, const sparse_t x,
                       const dense_t restrict y, index_t block_weight,
                       index_t block_length) {
    /* Avoid doing a modulo operation by precomputing the wrapping around. */
    struct mult_t queue[2 * block_weight];
    struct mult_t *restrict queue1 = queue;
    struct mult_t *restrict queue2 = queue + block_weight;
    for (index_t k = 0; k < block_weight; ++k) {
        queue1[k].z = z + x[k];
        queue1[k].y = y;
        queue1[k].len = block_length - x[k];
        queue2[k].z = z;
        queue2[k].y = y + block_length - x[k];
        queue2[k].len = x[k];
    }
    for (index_t k = 0; k < 2 * block_weight; ++k) {
        dense_t restrict y = queue[k].y;
        dense_t restrict z = queue[k].z;
        index_t len = queue[k].len;

        for (index_t i = 0; i < len; ++i) {
            z[i] ^= y[i];
        }
    }
}

/* Multiply the sparse vector 'x' of weight 'block_weight' by the dense
 * vector 'y' of length 'block_length' and add the result in 'z'. */
void multiply_add(dense_t restrict z, const sparse_t x,
                  const dense_t restrict y, index_t block_weight,
                  index_t block_length) {
    /* Avoid doing a modulo operation by precomputing the wrapping around. */
    struct mult_t queue[2 * block_weight];
    struct mult_t *restrict queue1 = queue;
    struct mult_t *restrict queue2 = queue + block_weight;
    for (index_t k = 0; k < block_weight; ++k) {
        queue1[k].z = z + x[k];
        queue1[k].y = y;
        queue1[k].len = block_length - x[k];
        queue2[k].z = z;
        queue2[k].y = y + block_length - x[k];
        queue2[k].len = x[k];
    }
    for (index_t k = 0; k < 2 * block_weight; ++k) {
        dense_t restrict y = queue[k].y;
        dense_t restrict z = queue[k].z;
        index_t len = queue[k].len;

        for (index_t i = 0; i < len; ++i) {
            z[i] += y[i];
        }
    }
}

#ifdef AVX
#define BUFF_LEN 8
/* Multiply modulo 2 the sparse vector 'x' of weight 'block_weight' by the dense
 * vector 'y' of length 'block_length' and xor the result in 'z'. */
void multiply_xor_mod2_avx2(dense_t restrict z, const sparse_t x,
                            const dense_t restrict y, index_t block_weight,
                            index_t block_length) {
    __m256i x_buff[BUFF_LEN];
    for (index_t i = 0; i < block_length / 32; i += BUFF_LEN) {
        for (index_t k = 0; k < BUFF_LEN; ++k)
            x_buff[k] = _mm256_load_si256((__m256i *)(z + 32 * (i + k)));

        for (index_t j = 0; j < block_weight; ++j) {
            index_t off = x[j] + 32 * i;
            for (index_t k = 0; k < BUFF_LEN; ++k)
                asm("vpxor   %[y], %[x], %[x]"
                    : [x] "+x"(x_buff[k])
                    : [y] "m"(*(__m256i *)&y[off + 32 * k])
                    :);
        }

        for (index_t k = 0; k < BUFF_LEN; ++k)
            _mm256_store_si256((__m256i *)(z + 32 * (i + k)), x_buff[k]);
    }
}

/* Multiply the sparse vector 'x' of weight 'block_weight' by the dense
 * vector 'y' of length 'block_length' and store the result in 'z'. */
void multiply_avx2(dense_t restrict z, const sparse_t x,
                   const dense_t restrict y, index_t block_weight,
                   index_t block_length) {
    __m256i x_buff[BUFF_LEN];
    for (index_t i = 0; i < block_length / 32; i += BUFF_LEN) {
        for (index_t k = 0; k < BUFF_LEN; ++k)
            x_buff[k] = _mm256_setzero_si256();

        for (index_t j = 0; j < block_weight; ++j) {
            index_t off = x[j] + 32 * i;
            for (index_t k = 0; k < BUFF_LEN; ++k)
                asm("vpaddb  %[y], %[x], %[x]"
                    : [x] "+x"(x_buff[k])
                    : [y] "m"(*(__m256i *)&y[off + 32 * k])
                    :);
        }

        for (index_t k = 0; k < BUFF_LEN; ++k)
            _mm256_store_si256((__m256i *)(z + 32 * (i + k)), x_buff[k]);
    }
}
#endif
