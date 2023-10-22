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
#include "errorgen.h"
#include "param.h"
#include "sparse_cyclic.h"

void generate_random_error(sparse_t e_block, index_t weight, prng_t prng) {
    sparse_rand(e_block, weight, INDEX * BLOCK_LENGTH, prng);
}

void generate_random_syndrome_error(sparse_t e_block, index_t weight,
                                    prng_t prng) {
    sparse_rand(e_block, weight, BLOCK_LENGTH, prng);
}

/* Generate an error pattern with a fixed number intersections and specific
 * total weight. */
void generate_around_word(sparse_t e_dst, index_t weight_dst, sparse_t e_src,
                          index_t weight_src, index_t intersections,
                          prng_t prng) {
    bit_t e[INDEX * BLOCK_LENGTH] = {0};
    bit_t h[INDEX * BLOCK_LENGTH] = {0};

    /* Pick a near-codeword. */
    index_t shift = prng->random_lim(BLOCK_LENGTH - 1, &prng->s0, &prng->s1);
    for (index_t l = 0; l < weight_src; ++l) {
        index_t k = e_src[l] >= BLOCK_LENGTH;
        index_t i =
            k ? (e_src[l] - BLOCK_LENGTH + shift) % BLOCK_LENGTH + BLOCK_LENGTH
              : (e_src[l] + shift) % BLOCK_LENGTH;
        h[i] = 1;
    }

    /* Pick `intersections` common positions with the near-codeword. */
    index_t error_weight = 0;
    while (error_weight < intersections) {
        index_t lrand = prng->random_lim(weight_src - 1, &prng->s0, &prng->s1);
        index_t k = e_src[lrand] >= BLOCK_LENGTH;
        index_t i = k ? (e_src[lrand] - BLOCK_LENGTH + shift) % BLOCK_LENGTH +
                            BLOCK_LENGTH
                      : (e_src[lrand] + shift) % BLOCK_LENGTH;
        if (!e[i]) {
            e[i] = 1;
            insert_sorted_noinc(e_dst, i, error_weight++);
        }
    }

    /* Complete the error pattern. */
    while (error_weight < weight_dst) {
        index_t jrand =
            prng->random_lim(INDEX * BLOCK_LENGTH - 1, &prng->s0, &prng->s1);
        if (!e[jrand] && !h[jrand]) {
            e[jrand] = 1;
            insert_sorted_noinc(e_dst, jrand, error_weight++);
        }
    }
}

/* Generate an error pattern with ERROR_FLOOR_P intersection with a
 * (BLOCK_WEIGHT, BLOCK_WEIGHT) near-codeword. */
void generate_near_codeword(sparse_t e_block, code_t *H, prng_t prng) {
    index_t e_src[BLOCK_WEIGHT];

    /* Pick a near-codeword. */
    index_t k = prng->random_lim(INDEX - 1, &prng->s0, &prng->s1);

    for (index_t l = 0; l < BLOCK_WEIGHT; ++l) {
        e_src[l] = k * BLOCK_LENGTH + H->columns[k][l];
    }

    generate_around_word(e_block, ERROR_WEIGHT, e_src, BLOCK_WEIGHT,
                         ERROR_FLOOR_P, prng);
}

/* Generate an error pattern with ERROR_FLOOR_P intersection with a
 * (2 * BLOCK_WEIGHT, ~2 * BLOCK_WEIGHT) near-codeword. */
void generate_near_codeword2(sparse_t e_block, code_t *H, prng_t prng) {
    index_t e_src[BLOCK_WEIGHT];

    for (index_t l = 0; l < BLOCK_WEIGHT; ++l) {
        e_src[l] = H->columns[0][l];
    }

    /* Pick a near-codeword. */
    index_t shift = prng->random_lim(BLOCK_LENGTH - 1, &prng->s0, &prng->s1);
    for (index_t l = 0; l < BLOCK_WEIGHT; ++l) {
        index_t i = H->columns[1][l] + shift;
        i = (i <= BLOCK_LENGTH) ? i : (i - BLOCK_LENGTH);
        e_src[BLOCK_WEIGHT + l] = BLOCK_LENGTH + i;
    }

    generate_around_word(e_block, ERROR_WEIGHT, e_src, 2 * BLOCK_WEIGHT,
                         ERROR_FLOOR_P, prng);
}

/* Generate an error pattern with ERROR_FLOOR_P intersection with a
 * codeword of weight 2 * BLOCK_WEIGHT. */
void generate_codeword(sparse_t e_block, code_t *H, prng_t prng) {
    index_t e_src[BLOCK_WEIGHT];

    for (index_t k = 0; k < INDEX; ++k) {
        index_t l;
        for (l = 0; l < BLOCK_WEIGHT; ++l) {
            index_t i = H->columns[INDEX - 1 - k][l];
            if (i >= BLOCK_LENGTH)
                break;
            e_src[k * BLOCK_WEIGHT + l] = i + k * BLOCK_LENGTH;
        }
        for (; l < BLOCK_WEIGHT; ++l) {
            index_t i = H->columns[INDEX - 1 - k][l];
            e_src[k * BLOCK_WEIGHT + l] = i + (k - 1) * BLOCK_LENGTH;
        }
    }

    generate_around_word(e_block, ERROR_WEIGHT, e_src, 2 * BLOCK_WEIGHT,
                         ERROR_FLOOR_P, prng);
}
