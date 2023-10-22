#include <string.h>

#include "code.h"
#include "sparse_cyclic.h"

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

void compute_codeword(cw_t codeword, code_t *H, msg_t message) {
    memset(codeword, 0, INDEX * 2 * SIZE_AVX * sizeof(bit_t));
#ifndef AVX
    for (index_t k = 0; k < INDEX; ++k) {
        multiply_xor_mod2(codeword[k], H->columns[INDEX - 1 - k], message,
                          BLOCK_WEIGHT, BLOCK_LENGTH);
    }
#else
    for (index_t k = 0; k < INDEX; ++k) {
        multiply_xor_mod2_avx2(codeword[k], H->rows[INDEX - 1 - k], message,
                               BLOCK_WEIGHT, SIZE_AVX);
    }
#endif
}

void compute_syndrome(syndrome_t *syndrome, code_t *H, e_t *e_dense) {
    memset(syndrome->vec, 0, 2 * SIZE_AVX * sizeof(bit_t));
#ifndef AVX
    for (index_t i = 0; i < INDEX; ++i) {
        multiply_xor_mod2(syndrome->vec, H->columns[i], e_dense->vec[i],
                          BLOCK_WEIGHT, BLOCK_LENGTH);
    }
#else
    for (index_t i = 0; i < INDEX; ++i) {
        multiply_xor_mod2_avx2(syndrome->vec, H->rows[i], e_dense->vec[i],
                               BLOCK_WEIGHT, SIZE_AVX);
    }
#endif

    syndrome->weight = 0;

    for (index_t j = 0; j < BLOCK_LENGTH; ++j) {
        syndrome->weight += syndrome->vec[j];
    }
}

/* Computing all the counters at is more efficient if we consider the
 * quasi-cyclic structure. */
void compute_counters(counters_t counters, bit_t *syndrome, code_t *H) {
    memcpy(syndrome + BLOCK_LENGTH, syndrome, BLOCK_LENGTH * sizeof(bit_t));
    for (index_t i = 0; i < INDEX; ++i) {
#ifndef AVX
        memset(counters[i], 0, BLOCK_LENGTH * sizeof(bit_t));
        multiply_add(counters[i], H->rows[i], syndrome, BLOCK_WEIGHT,
                     BLOCK_LENGTH);
#else
        multiply_avx2(counters[i], H->columns[i], syndrome, BLOCK_WEIGHT,
                      SIZE_AVX);
#endif
    }
}

void error_sparse_to_dense(e_t *e_dense, const sparse_t e_sparse,
                           index_t weight) {
    memset(e_dense->vec, 0, INDEX * 2 * SIZE_AVX * sizeof(bit_t));

    index_t k;
    for (k = 0; k < weight; ++k) {
        index_t j = e_sparse[k];
        if (j >= BLOCK_LENGTH)
            break;
        e_dense->vec[0][j] ^= 1;
    }
    for (; k < weight; ++k) {
        index_t j = e_sparse[k] - BLOCK_LENGTH;
        e_dense->vec[1][j] ^= 1;
    }
    for (index_t i = 0; i < INDEX; ++i) {
        memcpy(e_dense->vec[i] + BLOCK_LENGTH, e_dense->vec[i],
               BLOCK_LENGTH * sizeof(bit_t));
    }

    e_dense->weight = weight;
}

void syndrome_add_sparse_error(syndrome_t *syndrome, const sparse_t e_sparse,
                               index_t weight) {
    for (index_t k = 0; k < weight; ++k) {
        syndrome->vec[e_sparse[k]] ^= 1;
    }
}

void generate_random_message(msg_t message, prng_t prng) {
    for (index_t i = 0; i < BLOCK_LENGTH; i += 64) {
        index_t rand = prng->random_uint64_t(&prng->s0, &prng->s1);
        for (index_t j = 0; j < 64 && i + j < BLOCK_LENGTH; ++j) {
            bit_t b = (rand >> j) & 1L;
            message[i + j] = b;
            message[BLOCK_LENGTH + i + j] = b;
        }
    }
}
