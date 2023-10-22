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
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "code.h"
#include "codegen.h"
#include "errorgen.h"
#include "param.h"
#include "qcmdpc_decoder.h"
#include "types.h"
#include "xoshiro256plusplus.h"

#if (ALGO == BP)
#include "decoder_bp.h"
#else
#include "decoder.h"
#endif

void init_decoding_results(decoding_results_t *res, int n_threads,
                           int max_iter) {
    res->n_threads = n_threads;
    res->max_iter = max_iter;
    res->run = 0;

    res->n_test = calloc(n_threads, sizeof(long int));
    res->n_success = calloc(n_threads, sizeof(long int));
    res->n_iter = malloc(n_threads * sizeof(long int *));
    for (index_t i = 0; i < n_threads; ++i) {
        res->n_iter[i] = calloc(max_iter + 1, sizeof(long int));
    }
}

void clear_decoding_results(decoding_results_t *res) {
    free(res->n_test);
    free(res->n_success);
    for (index_t i = 0; i < res->n_threads; ++i) {
        free(res->n_iter[i]);
    }
    free(res->n_iter);
}

void sum_decoding_results(long int *test_total, long int *success_total,
                          long int *iter_total, const decoding_results_t *res) {
    *test_total = 0;
    *success_total = 0;
    memset(iter_total, 0, (res->max_iter + 1) * sizeof(long int));

    for (int i = 0; i < res->n_threads; ++i) {
        *test_total += res->n_test[i];
        *success_total += res->n_success[i];
    }

    for (int i = 0; i < res->n_threads; ++i) {
        for (int it = 0; it <= res->max_iter; ++it) {
            iter_total[it] += res->n_iter[i][it];
        }
    }
}

struct process_args {
    /* Number of test rounds */
    long int r;
    /* PRNG seeds */
    uint64_t s[4];

    int id;

    long int thread_iter;

    decoding_results_t *results;
};

void *process(void *arg) {
    struct process_args *args = arg;

    long int iter = args->thread_iter;
    decoding_results_t *results = args->results;
    int tid = args->id;

    code_t H;
    e_t e __attribute__((aligned(32)));
    syndrome_t syndrome __attribute__((aligned(32)));

    /* Error pattern */
    index_t error_sparse[ERROR_WEIGHT];

    /* Error pattern on the syndrome (for Ouroboros) */
#if OUROBOROS
    index_t syndrome_error_sparse[ERROR_WEIGHT / 2];
#endif

#if (ALGO == BP)
    decoder_bp_t dec = aligned_alloc(32, sizeof(struct decoder_bp));
#else
    decoder_t dec = aligned_alloc(32, sizeof(struct decoder));
#endif

    struct PRNG prng;
    memcpy(prng.s, args->s, 4 * sizeof(uint64_t));
    prng.random_lim = random_lim;
    prng.random_uint64_t = random_uint64_t;

    for (int i = 0; i < tid; ++i) {
        jump(prng.s);
    }
    init_decoder(dec, &H, &e, &syndrome);

    ++results->run;
    while (results->run && (iter == -1 || results->n_test[tid] < iter)) {
#if WEAK == 1
        generate_weak_type1(&H, &prng);
#elif WEAK == 2
        generate_weak_type2(&H, &prng);
#elif WEAK == 3
        generate_weak_type3(&H, &prng);
#else
        generate_random_code(&H, &prng);
#endif

#if ERROR_FLOOR == 1
        generate_near_codeword(error_sparse, &H, &prng);
#elif ERROR_FLOOR == 2
        generate_near_codeword2(error_sparse, &H, &prng);
#elif ERROR_FLOOR == 3
        generate_codeword(error_sparse, &H, &prng);
#else
        generate_random_error(error_sparse, ERROR_WEIGHT, &prng);
#endif

        reset_decoder(dec);
        error_sparse_to_dense(&e, error_sparse, ERROR_WEIGHT);

#if (ALGO == BP)
        init_bp(dec, &prng);
#else
        compute_syndrome(&syndrome, &H, &e);
#endif

#if OUROBOROS
        generate_random_syndrome_error(syndrome_error_sparse, SYNDROME_STOP,
                                       &prng);
        syndrome_add_sparse_error(&syndrome, syndrome_error_sparse,
                                  SYNDROME_STOP);
#endif

#if (ALGO == SBS) || (ALGO == SORT)
        if (qcmdpc_decode(dec, results->max_iter, &prng)) {
#else
        if (qcmdpc_decode(dec, results->max_iter)) {
#endif
            results->n_success[tid]++;
            results->n_iter[tid][dec->iter]++;
        }

        results->n_test[tid]++;
    }
    if (results->run)
        --results->run;

    free(dec);

    return NULL;
}

void decoder_loop(decoding_results_t *results, int n_threads, long int r) {
    /* PRNG seeds */
    uint64_t s[4] = {0};
    seed_random(s);

    struct process_args args[n_threads];
    for (int i = 0; i < n_threads; i++) {
        args[i].s[0] = s[0];
        args[i].s[1] = s[1];
        args[i].id = i;

        args[i].thread_iter = (r < 0) ? -1 : (i + r) / n_threads;

        args[i].results = results;
    }

    pthread_t threads[n_threads];

    for (int i = 0; i < n_threads; i++)
        pthread_create(&threads[i], NULL, process, (void *)&args[i]);

    for (int i = 0; i < n_threads; i++)
        pthread_join(threads[i], NULL);
}

void decoder_stop(decoding_results_t *res) { res->run = 0; }
