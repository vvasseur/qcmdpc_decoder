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
#include "codegen.h"
#include "param.h"
#include "sparse_cyclic.h"

void generate_random_code(code_t *H, prng_t prng) {
    for (index_t i = 0; i < INDEX; ++i) {
        sparse_rand(H->columns[i], BLOCK_WEIGHT, BLOCK_LENGTH, prng);
    }
    transpose_columns(H);
}

void generate_weak_type1(code_t *H, prng_t prng) {
    index_t k = prng->random_lim(INDEX, prng->s);
    index_t delta = 1 + prng->random_lim(BLOCK_LENGTH / 2, prng->s);
    index_t shift = prng->random_lim(BLOCK_LENGTH, prng->s);

    index_t length_left = BLOCK_LENGTH;
    for (index_t i = 0; i < WEAK_P; i++) {
        uint32_t a = (delta * (i + shift)) % BLOCK_LENGTH;
        insert_sorted_noinc(H->columns[k], a, i);
    }
    length_left -= WEAK_P;
    for (index_t i = WEAK_P; i < BLOCK_WEIGHT; i++) {
        uint32_t rand = prng->random_lim(length_left--, prng->s);
        insert_sorted(H->columns[k], rand, i);
    }

    sparse_rand(H->columns[INDEX - 1 - k], BLOCK_WEIGHT, BLOCK_LENGTH, prng);
    transpose_columns(H);
}

/* Generate a polynomial with a multiplicity of WEAK_P using the stars and bars
 * principle. */
void generate_weak_type2(code_t *H, prng_t prng) {
    index_t k = prng->random_lim(INDEX, prng->s);
    index_t delta = 1 + prng->random_lim(BLOCK_LENGTH / 2, prng->s);

    const index_t s = BLOCK_WEIGHT - WEAK_P;
    /* First the ois and zis represent the "bars". */
    index_t ois[s + 1];
    index_t zis[s + 1];
    ois[0] = 0;
    ois[s] = BLOCK_WEIGHT;
    index_t left = BLOCK_WEIGHT - 1;
    for (index_t i = 1; i < s; i++) {
        uint32_t rand = prng->random_lim(left--, prng->s);
        insert_sorted(ois, rand, i);
    }
    zis[0] = 0;
    zis[s] = BLOCK_LENGTH - BLOCK_WEIGHT;
    left = BLOCK_LENGTH - BLOCK_WEIGHT - 1;
    for (index_t i = 1; i < s; i++) {
        uint32_t rand = prng->random_lim(left--, prng->s);
        insert_sorted(zis, rand, i);
    }
    /* Convert ois and zis to run-length. */
    for (index_t i = 0; i < s; ++i) {
        ois[i] = ois[i + 1] - ois[i];
        zis[i] = zis[i + 1] - zis[i];
    }

    index_t shift = prng->random_lim(ois[0] + zis[0], prng->s);
    index_t current_pos = (BLOCK_LENGTH - shift) % BLOCK_LENGTH;
    index_t i = 0;
    for (index_t l1 = 0; l1 < s; ++l1) {
        current_pos = (current_pos + zis[l1]) % BLOCK_LENGTH;
        for (index_t l2 = 0; l2 < ois[l1]; ++l2) {
            uint32_t a = (delta * (current_pos + l2)) % BLOCK_LENGTH;
            insert_sorted_noinc(H->columns[k], a, i++);
        }
        current_pos = (current_pos + ois[l1]) % BLOCK_LENGTH;
    }

    sparse_rand(H->columns[INDEX - 1 - k], BLOCK_WEIGHT, BLOCK_LENGTH, prng);
    transpose_columns(H);
}

void generate_weak_type3(code_t *H, prng_t prng) {
    index_t length = BLOCK_LENGTH;
    index_t shift = prng->random_lim(BLOCK_LENGTH, prng->s);

    /* Choose WEAK_P common values. */
    for (index_t i = 0; i < WEAK_P; i++) {
        index_t rand = prng->random_lim(length--, prng->s);
        rand = insert_sorted(H->columns[0], rand, i);
        uint32_t a = (rand + shift) % BLOCK_LENGTH;
        insert_sorted_noinc(H->columns[1], a, i);
    }

    /* Complete H->columns[0]. */
    for (index_t i = WEAK_P; i < BLOCK_WEIGHT; i++) {
        index_t rand = prng->random_lim(length--, prng->s);
        insert_sorted(H->columns[0], rand, i);
    }
    length += BLOCK_WEIGHT - WEAK_P;

    /* Complete H->columns[1] without increasing the number of intersections
     * with H->columns[0].
     */
    for (index_t i = WEAK_P; i < BLOCK_WEIGHT; i++) {
    gen : {
        index_t rand = prng->random_lim(length, prng->s);

        for (index_t j = 0; j < i && H->columns[1][j] <= rand; j++, rand++)
            ;
        for (index_t k = 0; k < BLOCK_WEIGHT; k++) {
            uint32_t a = (H->columns[0][k] + shift) % BLOCK_LENGTH;
            if (a == rand)
                goto gen;
        }
        insert_sorted_noinc(H->columns[1], rand, i);
        --length;
    }
    }
    transpose_columns(H);
}
