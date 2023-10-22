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
#include <stdlib.h>
#include <string.h>

#include "param.h"
#include "sparse_cyclic.h"

void generate_random_error(sparse_t e_block, index_t weight, prng_t prng) {
    sparse_rand(e_block, weight, INDEX * BLOCK_LENGTH, prng);
}

void generate_random_syndrome_error(sparse_t e_block, index_t weight,
                                    prng_t prng) {
    sparse_rand(e_block, weight, BLOCK_LENGTH, prng);
}

/* Generate an error pattern with ERROR_FLOOR_P intersection with a
 * (BLOCK_WEIGHT, BLOCK_WEIGHT) near-codeword. */
void generate_near_codeword(sparse_t e_block, code_t *H, prng_t prng) {
    bit_t e[INDEX * BLOCK_LENGTH];
    bit_t h[INDEX * BLOCK_LENGTH];

    index_t max_inter;
    do {
        max_inter = 0;
        memset(e, 0, INDEX * BLOCK_LENGTH * sizeof(bit_t));
        memset(h, 0, INDEX * BLOCK_LENGTH * sizeof(bit_t));

        /* Pick a near-codeword. */
        index_t shift =
            prng->random_lim(BLOCK_LENGTH - 1, &prng->s0, &prng->s1);
        index_t k = prng->random_lim(INDEX - 1, &prng->s0, &prng->s1);
        {
            index_t l;
            for (l = 0; l < BLOCK_WEIGHT; ++l) {
                index_t i = H->columns[k][l] + shift;
                if (i >= BLOCK_LENGTH)
                    break;
                h[i + k * BLOCK_LENGTH] = 1;
            }
            for (; l < BLOCK_WEIGHT; ++l) {
                index_t i = H->columns[k][l] + shift;
                h[i + (k - 1) * BLOCK_LENGTH] = 1;
            }
        }

        /* Pick ERROR_FLOOR_P common positions with the near-codeword. */
        index_t error_weight = 0;
        while (error_weight < ERROR_FLOOR_P) {
            index_t lrand =
                prng->random_lim(BLOCK_WEIGHT - 1, &prng->s0, &prng->s1);
            index_t i = (H->columns[k][lrand] + shift) % BLOCK_LENGTH +
                        k * BLOCK_LENGTH;
            if (!e[i]) {
                e[i] = 1;
                insert_sorted_noinc(e_block, i, error_weight++);
            }
        }

        /* Complete the error pattern. */
        while (error_weight < ERROR_WEIGHT) {
            index_t jrand = prng->random_lim(INDEX * BLOCK_LENGTH - 1,
                                             &prng->s0, &prng->s1);
            if (!e[jrand] && !h[jrand]) {
                e[jrand] = 1;
                insert_sorted_noinc(e_block, jrand, error_weight++);
            }
        }

        /* Check that no other near-codeword is closer. */
        index_t block_inter[INDEX] = {0};
        for (index_t l1 = 0; l1 < ERROR_WEIGHT; ++l1) {
            index_t index = e_block[l1] / BLOCK_LENGTH;
            for (index_t l2 = 0; l2 < BLOCK_WEIGHT; ++l2) {
                index_t inter = 0;
                for (index_t l3 = 0; l3 < BLOCK_WEIGHT; ++l3) {
                    index_t j =
                        (e_block[l1] + BLOCK_LENGTH - H->columns[index][l2] +
                         H->columns[index][l3]) %
                            BLOCK_LENGTH +
                        index * BLOCK_LENGTH;
                    inter += e[j];
                }
                if (inter > block_inter[index])
                    block_inter[index] = inter;
            }
        }
        max_inter =
            (block_inter[0] > block_inter[1]) ? block_inter[0] : block_inter[1];
    } while (max_inter != ERROR_FLOOR_P);
}

/* Generate an error pattern with ERROR_FLOOR_P intersection with a
 * (2 * BLOCK_WEIGHT, ~2 * BLOCK_WEIGHT) near-codeword. */
void generate_near_codeword2(sparse_t e_block, code_t *H, prng_t prng) {
    bit_t e[INDEX * BLOCK_LENGTH];
    bit_t h[INDEX * BLOCK_LENGTH];

    index_t max_inter;
    do {
        max_inter = 0;
        memset(e, 0, INDEX * BLOCK_LENGTH * sizeof(bit_t));
        memset(h, 0, INDEX * BLOCK_LENGTH * sizeof(bit_t));

        /* Pick a near-codeword. */
        index_t shift[INDEX];
        for (index_t k = 0; k < INDEX; ++k) {
            shift[k] = prng->random_lim(BLOCK_LENGTH - 1, &prng->s0, &prng->s1);
            index_t l;
            for (l = 0; l < BLOCK_WEIGHT; ++l) {
                index_t i = H->columns[k][l] + shift[k];
                if (i >= BLOCK_LENGTH)
                    break;
                h[i + k * BLOCK_LENGTH] = 1;
            }
            for (; l < BLOCK_WEIGHT; ++l) {
                index_t i = H->columns[k][l] + shift[k];
                h[i + (k - 1) * BLOCK_LENGTH] = 1;
            }
        }

        /* Pick ERROR_FLOOR_P common positions with the near-codeword. */
        index_t error_weight = 0;
        while (error_weight < ERROR_FLOOR_P) {
            index_t krand = prng->random_lim(INDEX - 1, &prng->s0, &prng->s1);
            index_t lrand =
                prng->random_lim(BLOCK_WEIGHT - 1, &prng->s0, &prng->s1);
            index_t i =
                (H->columns[krand][lrand] + shift[krand]) % BLOCK_LENGTH +
                krand * BLOCK_LENGTH;
            if (!e[i]) {
                e[i] = 1;
                insert_sorted_noinc(e_block, i, error_weight++);
            }
        }

        /* Complete the error pattern. */
        while (error_weight < ERROR_WEIGHT) {
            index_t jrand = prng->random_lim(INDEX * BLOCK_LENGTH - 1,
                                             &prng->s0, &prng->s1);
            if (!e[jrand] && !h[jrand]) {
                e[jrand] = 1;
                insert_sorted_noinc(e_block, jrand, error_weight++);
            }
        }

        /* Check that no other near-codeword is closer. */
        index_t block_inter[INDEX] = {0};
        for (index_t l1 = 0; l1 < ERROR_WEIGHT; ++l1) {
            index_t index = e_block[l1] / BLOCK_LENGTH;
            for (index_t l2 = 0; l2 < BLOCK_WEIGHT; ++l2) {
                index_t inter = 0;
                for (index_t l3 = 0; l3 < BLOCK_WEIGHT; ++l3) {
                    index_t j =
                        (e_block[l1] + BLOCK_LENGTH - H->columns[index][l2] +
                         H->columns[index][l3]) %
                            BLOCK_LENGTH +
                        index * BLOCK_LENGTH;
                    inter += e[j];
                }
                if (inter > block_inter[index])
                    block_inter[index] = inter;
            }
        }
        for (index_t i = 0; i < INDEX; ++i)
            max_inter += block_inter[i];
    } while (max_inter != ERROR_FLOOR_P);
}

/* Generate an error pattern with ERROR_FLOOR_P intersection with a
 * codeword of weight 2 * BLOCK_WEIGHT. */
void generate_codeword(sparse_t e_block, code_t *H, prng_t prng) {
    bit_t e[INDEX * BLOCK_LENGTH];
    bit_t h[INDEX * BLOCK_LENGTH];

    index_t max_inter;
    do {
        max_inter = 0;
        memset(e, 0, INDEX * BLOCK_LENGTH * sizeof(bit_t));
        memset(h, 0, INDEX * BLOCK_LENGTH * sizeof(bit_t));

        /* Pick a codeword. */
        index_t shift =
            prng->random_lim(BLOCK_LENGTH - 1, &prng->s0, &prng->s1);
        for (index_t k = 0; k < INDEX; ++k) {
            index_t l;
            for (l = 0; l < BLOCK_WEIGHT; ++l) {
                index_t i = H->columns[INDEX - 1 - k][l] + shift;
                if (i >= BLOCK_LENGTH)
                    break;
                h[i + k * BLOCK_LENGTH] = 1;
            }
            for (; l < BLOCK_WEIGHT; ++l) {
                index_t i = H->columns[INDEX - 1 - k][l] + shift;
                h[i + (k - 1) * BLOCK_LENGTH] = 1;
            }
        }

        /* Pick ERROR_FLOOR_P common positions with the codeword. */
        index_t error_weight = 0;
        while (error_weight < ERROR_FLOOR_P) {
            index_t krand = prng->random_lim(INDEX - 1, &prng->s0, &prng->s1);
            index_t lrand =
                prng->random_lim(BLOCK_WEIGHT - 1, &prng->s0, &prng->s1);
            index_t i = H->columns[krand][lrand] + (1 - krand) * BLOCK_LENGTH;
            if (!e[i]) {
                e[i] = 1;
                insert_sorted_noinc(e_block, i, error_weight++);
            }
        }

        /* Complete the error pattern. */
        while (error_weight < ERROR_WEIGHT) {
            index_t jrand = prng->random_lim(INDEX * BLOCK_LENGTH - 1,
                                             &prng->s0, &prng->s1);
            if (!e[jrand] && !h[jrand]) {
                e[jrand] = 1;
                insert_sorted_noinc(e_block, jrand, error_weight++);
            }
        }

        /* Check that no other codeword is closer. */
        for (index_t l1 = 0; l1 < ERROR_WEIGHT; ++l1) {
            index_t index = INDEX - 1 - e_block[l1] / BLOCK_LENGTH;
            for (index_t l2 = 0; l2 < BLOCK_WEIGHT; ++l2) {
                index_t inter = 0;
                for (index_t l3 = 0; l3 < BLOCK_WEIGHT; ++l3) {
                    index_t j =
                        (e_block[l1] + BLOCK_LENGTH - H->columns[index][l2] +
                         H->columns[index][l3]) %
                            BLOCK_LENGTH +
                        (INDEX - 1 - index) * BLOCK_LENGTH;
                    inter += e[j];
                }
                for (index_t l3 = 0; l3 < BLOCK_WEIGHT; ++l3) {
                    index_t j =
                        (e_block[l1] + BLOCK_LENGTH - H->columns[index][l2] +
                         H->columns[INDEX - 1 - index][l3]) %
                            BLOCK_LENGTH +
                        index * BLOCK_LENGTH;
                    inter += e[j];
                }
                if (inter > max_inter)
                    max_inter = inter;
            }
        }
    } while (max_inter != ERROR_FLOOR_P);
}
