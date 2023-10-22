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
#include "param.h"
#if (ALGO == BP)
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "decoder_bp.h"
#include "sparse_cyclic.h"
#include "threshold.h"

static void to_binary(decoder_bp_t dec);
static void compute_codeword(decoder_bp_t dec);
static void compute_syndrome(decoder_bp_t dec);
static size_t get_error_weight(decoder_bp_t dec);
static size_t get_syndrome_weight(decoder_bp_t dec);
static void to_binary(decoder_bp_t dec);
static void extrinsic(llr_t tree[1 << LOG2(INDEX * BLOCK_LENGTH)], size_t size,
                      llr_t (*op)(llr_t, llr_t));

void random_message(decoder_bp_t dec, prng_t prng) {
    for (index_t i = 0; i < BLOCK_LENGTH; i += 64) {
        index_t rand = prng->random_uint64_t(&prng->s0, &prng->s1);
        for (index_t j = 0; j < 64 && i + j < BLOCK_LENGTH; ++j) {
            bit_t b = (rand >> j) & 1L;
            dec->message[i + j] = b;
            dec->message[BLOCK_LENGTH + i + j] = b;
        }
    }
}

static void compute_codeword(decoder_bp_t dec) {
    memset(dec->codeword, 0, INDEX * 2 * SIZE_AVX * sizeof(bit_t));
#ifndef AVX
    for (index_t k = 0; k < INDEX; ++k) {
        multiply_xor_mod2(dec->codeword[k], dec->H.columns[INDEX - 1 - k],
                          dec->message, BLOCK_WEIGHT, BLOCK_LENGTH);
    }
#else
    for (index_t k = 0; k < INDEX; ++k) {
        multiply_xor_mod2_avx2(dec->codeword[k], dec->H.rows[INDEX - 1 - k],
                               dec->message, BLOCK_WEIGHT, SIZE_AVX);
    }
#endif
}

static void compute_syndrome(decoder_bp_t dec) {
    memset(dec->syndrome, 0, 2 * SIZE_AVX * sizeof(bit_t));
#ifndef AVX
    for (index_t k = 0; k < INDEX; ++k) {
        multiply_xor_mod2(dec->syndrome, dec->H.columns[k], dec->bits[k],
                          BLOCK_WEIGHT, BLOCK_LENGTH);
    }
#else
    for (index_t k = 0; k < INDEX; ++k) {
        multiply_xor_mod2_avx2(dec->syndrome, dec->H.rows[k], dec->bits[k],
                               BLOCK_WEIGHT, SIZE_AVX);
    }
#endif
}

static size_t get_error_weight(decoder_bp_t dec) {
    size_t weight = 0;
    for (index_t k = 0; k < INDEX; ++k)
        for (index_t j = 0; j < BLOCK_LENGTH; ++j)
            weight += dec->codeword[k][j] ^ dec->bits[k][j];
    return weight;
}

static size_t get_syndrome_weight(decoder_bp_t dec) {
    size_t weight = 0;
    for (index_t j = 0; j < BLOCK_LENGTH; ++j)
        weight += dec->syndrome[j];
    return weight;
}

static void to_binary(decoder_bp_t dec) {
    memset(dec->bits, 0, INDEX * 2 * SIZE_AVX * sizeof(bit_t));
    for (index_t k = 0; k < INDEX; ++k)
        for (index_t j = 0; j < BLOCK_LENGTH; ++j) {
            llr_t val = dec->r[k][j];
            for (index_t l = 0; l < BLOCK_WEIGHT; ++l)
                val += dec->v_to_c[k][l][j];
            dec->bits[k][j] = (val < 0);
            dec->bits[k][BLOCK_LENGTH + j] = (val < 0);
        }
}

void init_decoder_error(decoder_bp_t dec, const sparse_t e_block,
                        const sparse_t e2_block) {
    compute_codeword(dec);

    const llr_t proba_init =
        log((llr_t)(INDEX * BLOCK_LENGTH - ERROR_WEIGHT) / ERROR_WEIGHT);
    for (index_t k = 0; k < INDEX; ++k)
        for (index_t j = 0; j < BLOCK_LENGTH; ++j) {
            dec->r[k][j] = (1 - 2 * dec->codeword[k][j]) * proba_init;
        }

    {
        index_t k;
        for (k = 0; k < ERROR_WEIGHT; ++k) {
            index_t j = e_block[k];
            if (j >= BLOCK_LENGTH)
                break;
            dec->r[0][j] = -dec->r[0][j];
        }
        for (; k < ERROR_WEIGHT; ++k) {
            index_t j = e_block[k] - BLOCK_LENGTH;
            dec->r[1][j] = -dec->r[1][j];
        }
    }

    for (index_t k = 0; k < INDEX; ++k) {
        for (index_t j = 0; j < BLOCK_LENGTH; ++j) {
            for (index_t l = 0; l < BLOCK_WEIGHT; ++l) {
                dec->v_to_c[k][l][j] = SATURATE(dec->r[k][j], BP_SATURATE);
            }
        }
    }
    (void)e2_block;
}

void reset_decoder(decoder_bp_t dec) { (void)dec; }

/* Operation functions to use with the extrinsic function. */
static llr_t add(llr_t a, llr_t b) { return a + b; }

static llr_t mult(llr_t a, llr_t b) { return a * b; }

/* Use a tree that is traversed from bottom to top and then from top to
 * bottom to compute all the sums (or products) of (size - 1) values among
 * 'size'. */
static void extrinsic(llr_t tree[1 << LOG2(INDEX * BLOCK_LENGTH)], size_t size,
                      llr_t (*op)(llr_t, llr_t)) {
    const size_t lsize = LOG2(size);

    size_t sizes[lsize + 1];
    char carries[lsize + 1];
    sizes[lsize] = size;
    carries[lsize] = 0;
    for (size_t k = lsize; k-- > 1;) {
        size_t n = (sizes[k + 1] + carries[k + 1]) >> 1;
        char c = (sizes[k + 1] + carries[k + 1]) & 1;
        sizes[k] = n;
        carries[k] = c;
        llr_t *restrict tree_current = tree + (1 << k);
        llr_t *restrict tree_below = tree + (1 << (k + 1));
        for (size_t i = 0; i < n; ++i)
            tree_current[i] = op(tree_below[2 * i], tree_below[2 * i + 1]);
        if (c)
            tree_current[n] = tree_below[2 * n];
    }

    SWAP(tree[2], tree[3]);
    for (size_t k = 2; k <= lsize; ++k) {
        size_t n = sizes[k - 1];
        char c = carries[k - 1];
        llr_t *restrict tree_current = tree + (1 << k);
        llr_t *restrict tree_above = tree + (1 << (k - 1));
        for (size_t i = 0; i < n; ++i) {
            SWAP(tree_current[2 * i], tree_current[2 * i + 1]);
            tree_current[2 * i] = op(tree_current[2 * i], tree_above[i]);
            tree_current[2 * i + 1] =
                op(tree_current[2 * i + 1], tree_above[i]);
        }
        if (c)
            tree_current[2 * n] = tree_above[n];
    }
}

int qcmdpc_decode(decoder_bp_t dec, int max_iter) {
    llr_t *messages_v = dec->tree + (1L << LOG2(BLOCK_WEIGHT));
    llr_t(*messages_c)[BLOCK_WEIGHT] = (llr_t(*)[BLOCK_WEIGHT])(
        dec->tree + (1L << LOG2(INDEX * BLOCK_WEIGHT)));

    dec->iter = 0;
    while (dec->iter < max_iter) {
        ++dec->iter;
        for (index_t i = 0; i < BLOCK_LENGTH; ++i) {
            for (index_t k = 0; k < INDEX; ++k)
                for (index_t l = 0; l < BLOCK_WEIGHT; ++l) {
                    index_t j =
                        ((i > dec->H.columns[k][l]) ? 0 : BLOCK_LENGTH) + i -
                        dec->H.columns[k][l];
                    messages_c[k][l] = tanh(dec->v_to_c[k][l][j] / 2);
                }
            extrinsic(dec->tree, INDEX * BLOCK_WEIGHT, mult);
            for (index_t k = 0; k < INDEX; ++k)
                for (index_t l = 0; l < BLOCK_WEIGHT; ++l) {
                    index_t j =
                        ((i > dec->H.columns[k][l]) ? 0 : BLOCK_LENGTH) + i -
                        dec->H.columns[k][l];
                    dec->c_to_v[k][l][j] = SATURATE(
                        2 * atanh(messages_c[k][l]) * BP_SCALE, BP_SATURATE);
                }
        }
        for (index_t k = 0; k < INDEX; ++k)
            for (index_t j = 0; j < BLOCK_LENGTH; ++j) {
                for (index_t l = 0; l < BLOCK_WEIGHT; ++l)
                    messages_v[l] = dec->c_to_v[k][l][j];

                extrinsic(dec->tree, BLOCK_WEIGHT, add);
                for (index_t l = 0; l < BLOCK_WEIGHT; ++l) {
                    dec->v_to_c[k][l][j] =
                        SATURATE(dec->r[k][j] + messages_v[l], BP_SATURATE);
                }
            }
        to_binary(dec);
        compute_syndrome(dec);
        if (get_syndrome_weight(dec) == SYNDROME_STOP)
            break;
    }
    return !get_error_weight(dec);
}
#endif
