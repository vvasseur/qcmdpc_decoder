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
#include <getopt.h>
#include <pthread.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "param.h"
#include "qcmdpc_decoder.h"

decoding_results_t *current_results = NULL;
pthread_t *print_thread = NULL;

#define _GNU_SOURCE

static void print_parameters(FILE *f);
static void print_usage(FILE *f, char *arg0);
static void print_stats(FILE *f);
static void inthandler(int signo);
static void huphandler(int signo);
static void parse_arguments(int argc, char *argv[], int *max_iter, long int *N,
                            int *threads, int *quiet);

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
#endif
#if (ALGO == GRAY_B) || (ALGO == GRAY_BGF) || (ALGO == GRAY_BGB) ||            \
    (ALGO == GRAY_BG)
            "-DTHRESHOLD_C0=%lg "
            "-DTHRESHOLD_C1=%lg "
#endif
#if (ALGO == BACKFLIP2)
            "-DTHRESHOLD_A0=%lg "
            "-DTHRESHOLD_A1=%lg "
            "-DTHRESHOLD_A2=%lg "
            "-DTHRESHOLD_A3=%lg "
            "-DTHRESHOLD_A4=%lg "
#endif
#if (ALGO == BACKFLIP)
            "-DTTL_C0=%lg "
            "-DTTL_C1=%lg "
#endif
#if (ALGO == BACKFLIP2) || (ALGO == BACKFLIP)
            "-DTTL_SATURATE=%d "
#endif
#if (ALGO == SORT)
            "-DGRAY_SIZE=%d "
#endif
            "-DALGO=%s\n",
            INDEX, BLOCK_LENGTH, BLOCK_WEIGHT, ERROR_WEIGHT, OUROBOROS, WEAK,
            WEAK_P, ERROR_FLOOR, ERROR_FLOOR_P,
#if (ALGO == BP)
            BP_SCALE, BP_SATURATE,
#endif
#if (ALGO == GRAY_B) || (ALGO == GRAY_BGF) || (ALGO == GRAY_BGB) ||            \
    (ALGO == GRAY_BG)
            THRESHOLD_C0, THRESHOLD_C1,
#endif
#if (ALGO == BACKFLIP2)
            THRESHOLD_A0, THRESHOLD_A1, THRESHOLD_A2, THRESHOLD_A3,
            THRESHOLD_A4,
#endif
#if (ALGO == BACKFLIP)
            TTL_C0, TTL_C1,
#endif
#if (ALGO == BACKFLIP2) || (ALGO == BACKFLIP)
            TTL_SATURATE,
#endif
#if (ALGO == SORT)
            GRAY_SIZE,
#endif
            algo[ALGO]);
    fflush(f);
}

static void print_usage(FILE *f, char *arg0) {
    fprintf(f,
            "usage: %s [OPTIONS]\n"
            "\n"
            "-i, --max-iter         maximum number of iterations\n"
            "-N, --rounds           number of rounds to perform\n"
            "-T, --threads          number of threads to use\n"
            "-q, --quiet            do not regularly output results (only on "
            "SIGHUP)\n",
            arg0);
    exit(2);
}

static void print_stats(FILE *f) {
    if (!current_results->n_test && !current_results->n_success)
        return;
    long int n_test_total;
    long int n_success_total;
    long int n_iter_total[current_results->max_iter + 1];

    sum_decoding_results(&n_test_total, &n_success_total, n_iter_total,
                         current_results);

    fprintf(f, "%ld", n_test_total);
    for (int it = 0; it <= current_results->max_iter; ++it) {
        if (n_iter_total[it])
            fprintf(f, " %d:%ld", it, n_iter_total[it]);
    }
    if (n_success_total != n_test_total)
        fprintf(f, " >%d:%ld", current_results->max_iter,
                n_test_total - n_success_total);
    fprintf(f, "\n");
    fflush(f);
}

static void inthandler(int signo) {
    (void)signo;
    if (print_thread)
        pthread_cancel(*print_thread);
    decoder_stop(current_results);
}

static void huphandler(int signo) {
    (void)signo;
    print_stats(stdout);
}

static void parse_arguments(int argc, char *argv[], int *max_iter, long int *N,
                            int *threads, int *quiet) {
    const char *options = "i:N:T:q";
    static struct option longopts[] = {{"max-iter", required_argument, 0, 'i'},
                                       {"rounds", required_argument, 0, 'N'},
                                       {"threads", required_argument, 0, 'T'},
                                       {"quiet", no_argument, 0, 'q'},
                                       {NULL, 0, 0, 0}};

    int ch;
    while ((ch = getopt_long(argc, argv, options, longopts, NULL)) != -1) {
        switch (ch) {
        case 'i':
            *max_iter = atol(optarg);
            if (*max_iter < 0)
                print_usage(stderr, argv[0]);
            break;
        case 'N':
            *N = atol(optarg);
            if (*N < 1)
                print_usage(stderr, argv[0]);
            break;
        case 'T':
            *threads = atoi(optarg);
            if (*threads <= 0)
                print_usage(stderr, argv[0]);
            break;
        case 'q':
            *quiet = 1;
            break;
        default:
            print_usage(stderr, argv[0]);
            break;
        }
    }
}

void *print(void *arg) {
    (void)arg;
    do {
        sleep(TIME_BETWEEN_PRINTS);
        print_stats(stdout);
    } while (current_results->run);

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
    int quiet = 0;
    int n_threads = 1;
    int max_iter = 100;
    decoding_results_t results;
    current_results = &results;

    parse_arguments(argc, argv, &max_iter, &r, &n_threads, &quiet);
    print_parameters(stdout);

    /* Keep independent statistics for all threads. */
    init_decoding_results(&results, n_threads, max_iter);

    if (!quiet) {
        print_thread = malloc(sizeof(pthread_t));
        pthread_create(print_thread, NULL, print, (void *)NULL);
    }

    decoder_loop(&results, n_threads, r);

    if (!quiet)
        pthread_cancel(*print_thread);

    print_stats(stdout);

    clear_decoding_results(&results);

    exit(EXIT_SUCCESS);
}
