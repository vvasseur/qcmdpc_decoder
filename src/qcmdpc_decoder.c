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
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include "param.h"

#include "cli.h"
#if (ALGO == BP)
#include "decoder_bp.h"
#else
#include "decoder.h"
#endif
#include "error_floor.h"
#include "sparse_cyclic.h"
#include "threshold.h"
#include "weak.h"

static void print_parameters(FILE *f);
static void print_stats(FILE *f, long int *n_test, long int *n_success);
static void inthandler(int signo);

static long int *n_test = NULL;
static long int *n_success = NULL;
static long int **n_iter = NULL;
static int n_threads = 1;
static int max_iter = 100;
static int stop = 0;

static void print_parameters(FILE *f) {
    const char *algo[] = {"CLASSIC",  "BACKFLIP", "BACKFLIP2", "SBS", "GRAY_B",
                          "GRAY_BGF", "GRAY_BGB", "GRAY_BG",   "BP",  "SORT"};

    fprintf(f,
            "-DINDEX=%d "
            "-DBLOCK_LENGTH=%d "
            "-DBLOCK_WEIGHT=%d "
            "-DERROR_WEIGHT=%d "
            "-DOUROBOROS=%d "
            "-DWEAK=%d "
            "-DWEAK_P=%d "
            "-DERROR_FLOOR=%d "
            "-DERROR_FLOOR_P=%d "
#if (ALGO == BP)
            "-DBP_SCALE=%lg "
            "-DBP_SATURATE=%lg "
#elif (ALGO == GRAY_B) || (ALGO == GRAY_BGF) || (ALGO == GRAY_BGB) ||          \
    (ALGO == GRAY_BG)
            "-DTHRESHOLD_C0=%lg "
            "-DTHRESHOLD_C1=%lg "
#elif (ALGO == BACKFLIP2)
            "-DTHRESHOLD_A0=%lg "
            "-DTHRESHOLD_A1=%lg "
            "-DTHRESHOLD_A2=%lg "
            "-DTHRESHOLD_A3=%lg "
            "-DTHRESHOLD_A4=%lg "
#elif (ALGO == BACKFLIP)
            "-DTTL_C0=%lg "
            "-DTTL_C1=%lg "
#elif (ALGO == BACKFLIP2) || (ALGO == BACKFLIP)
            "-DTTL_SATURATE=%d "
#elif (ALGO == SORT)
            "-DGRAY_SIZE=%d "
#endif
            "-DALGO=%s\n",
            INDEX, BLOCK_LENGTH, BLOCK_WEIGHT, ERROR_WEIGHT, OUROBOROS, WEAK,
            WEAK_P, ERROR_FLOOR, ERROR_FLOOR_P,
#if (ALGO == BP)
            BP_SCALE, BP_SATURATE,
#elif (ALGO == GRAY_B) || (ALGO == GRAY_BGF) || (ALGO == GRAY_BGB) ||          \
    (ALGO == GRAY_BG)
            THRESHOLD_C0, THRESHOLD_C1,
#elif (ALGO == BACKFLIP2)
            THRESHOLD_A0, THRESHOLD_A1, THRESHOLD_A2, THRESHOLD_A3,
            THRESHOLD_A4,
#elif (ALGO == BACKFLIP)
            TTL_C0, TTL_C1,
#elif (ALGO == BACKFLIP2) || (ALGO == BACKFLIP)
            TTL_SATURATE,
#elif (ALGO == SORT)
            GRAY_SIZE,
#endif
            algo[ALGO]);
    fflush(f);
}

static void print_stats(FILE *f, long int *n_test, long int *n_success) {
    if (!n_test && !n_success)
        return;
    long int n_test_total = 0;
    long int n_success_total = 0;
    for (int i = 0; i < n_threads; ++i) {
        n_test_total += n_test[i];
        n_success_total += n_success[i];
    }
    long int *n_iter_total = calloc(max_iter + 1, sizeof(long int));

    for (int i = 0; i < n_threads; ++i) {
        for (int it = 0; it <= max_iter; ++it) {
            n_iter_total[it] += n_iter[i][it];
        }
    }

    fprintf(f, "%ld", n_test_total);
    for (int it = 0; it <= max_iter; ++it) {
        if (n_iter_total[it])
            fprintf(f, " %d:%ld", it, n_iter_total[it]);
    }
    if (n_success_total != n_test_total)
        fprintf(f, " >%d:%ld", max_iter, n_test_total - n_success_total);
    fprintf(f, "\n");
    fflush(f);
    free(n_iter_total);
}

static void inthandler(int signo) {
    (void)signo;
    stop = 1;
}
static void huphandler(int signo) {
    (void)signo;
    print_stats(stdout, n_test, n_success);
}

struct process_args {
    /* Number of test rounds */
    long int r;
    /* PRNG seeds */
    uint64_t s[2];

    int id;
};

void *print(void *arg) {
    (void)arg;
    while (!stop) {
        sleep(TIME_BETWEEN_PRINTS);
        print_stats(stdout, n_test, n_success);
    }

    return NULL;
}

void *process(void *arg) {
    struct process_args *args = arg;

    /* Error pattern */
    sparse_t e_block = sparse_new(ERROR_WEIGHT);

    /* Error pattern on the syndrome (for Ouroboros) */
#if !OUROBOROS
    sparse_t e2_block = NULL;
#else
    sparse_t e2_block = sparse_new(ERROR_WEIGHT / 2);
#endif

#if (ALGO == BP)
    struct decoder_bp *dec = aligned_alloc(32, sizeof(struct decoder_bp));
#else
    struct decoder *dec = aligned_alloc(32, sizeof(struct decoder));
#endif

    prng_t prng = malloc(sizeof(struct PRNG));
    prng->s0 = args->s[0];
    prng->s1 = args->s[1];
    prng->random_lim = random_lim;
    prng->random_uint64_t = random_uint64_t;

    for (int i = 0; i < args->id; ++i) {
        jump(&prng->s0, &prng->s1);
    }

    long int thread_total_tests = (args->id + args->r) / n_threads;
    while (!stop && (args->r == -1 || n_test[args->id] < thread_total_tests)) {
#if (ALGO == BP)
        random_message(dec, prng);
#endif
#if WEAK == 1
        generate_weak_type1(dec->Hcolumns, prng);
#elif WEAK == 2
        generate_weak_type2(dec->Hcolumns, prng);
#elif WEAK == 3
        generate_weak_type3(dec->Hcolumns, prng);
#else
        for (index_t i = 0; i < INDEX; ++i) {
            sparse_rand(dec->Hcolumns[i], BLOCK_WEIGHT, BLOCK_LENGTH, prng);
        }
#endif

#if ERROR_FLOOR == 1
        near_codeword(e_block, dec->Hcolumns, prng);
#elif ERROR_FLOOR == 2
        near_codeword2(e_block, dec->Hcolumns, prng);
#elif ERROR_FLOOR == 3
        codeword(e_block, dec->Hcolumns, prng);
#else
        sparse_rand(e_block, ERROR_WEIGHT, INDEX * BLOCK_LENGTH, prng);
#endif
#if OUROBOROS
        sparse_rand(e2_block, SYNDROME_STOP, BLOCK_LENGTH, prng);
#endif
        reset_decoder(dec);
        init_decoder_error(dec, e_block, e2_block);

#if (ALGO == SBS) || (ALGO == SORT)
        if (qcmdpc_decode(dec, max_iter, prng)) {
#else
        if (qcmdpc_decode(dec, max_iter)) {
#endif
            n_success[args->id]++;
            n_iter[args->id][dec->iter]++;
        }

        n_test[args->id]++;
    }
    free(prng);
    sparse_free(e_block);
    if (e2_block) {
        sparse_free(e2_block);
    }
    free(dec);

    return NULL;
}

int main(int argc, char *argv[]) {
    struct sigaction action_hup;
    action_hup.sa_handler = huphandler;
    sigemptyset(&action_hup.sa_mask);
    action_hup.sa_flags = 0;

    struct sigaction action_int;
    action_int.sa_handler = inthandler;
    sigemptyset(&action_int.sa_mask);
    action_int.sa_flags = 0;

    sigaction(SIGINT, &action_int, NULL);
    sigaction(SIGTERM, &action_int, NULL);
    sigaction(SIGHUP, &action_hup, NULL);

    /* Number of test rounds */
    long int r = -1;
    /* PRNG seeds */
    uint64_t s[2];
    int quiet = 0;

    seed_random(&s[0], &s[1]);

    parse_arguments(argc, argv, &max_iter, &r, &n_threads, &quiet);
    print_parameters(stdout);

    struct process_args args[n_threads];
    for (int i = 0; i < n_threads; i++) {
        args[i].r = r;
        args[i].s[0] = s[0];
        args[i].s[1] = s[1];
        args[i].id = i;
    }

    /* Keep independent statistics for all threads. */
    n_test = calloc(n_threads, sizeof(long int));
    n_success = calloc(n_threads, sizeof(long int));
    n_iter = malloc(n_threads * sizeof(long int *));
    for (index_t i = 0; i < n_threads; ++i) {
        n_iter[i] = calloc(max_iter + 1, sizeof(long int));
    }

    pthread_t threads[n_threads + 1];

    for (int i = 0; i < n_threads; i++)
        pthread_create(&threads[i], NULL, process, (void *)&args[i]);
    if (!quiet)
        pthread_create(&threads[n_threads], NULL, print, NULL);

    for (int i = 0; i < n_threads; i++)
        pthread_join(threads[i], NULL);
    pthread_cancel(threads[n_threads]);

    print_stats(stdout, n_test, n_success);

    free(n_test);
    free(n_success);
    for (index_t i = 0; i < n_threads; ++i) {
        free(n_iter[i]);
    }
    free(n_iter);
    exit(EXIT_SUCCESS);
}
