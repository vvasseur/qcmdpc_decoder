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
#include <stdio.h>
#include <stdlib.h>

#include "cli.h"

#define _GNU_SOURCE

void print_usage(char *arg0) {
    fprintf(stderr,
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

void parse_arguments(int argc, char *argv[], int *max_iter, long int *N,
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
                print_usage(argv[0]);
            break;
        case 'N':
            *N = atol(optarg);
            if (*N < 1)
                print_usage(argv[0]);
            break;
        case 'T':
            *threads = atoi(optarg);
            if (*threads <= 0)
                print_usage(argv[0]);
            break;
        case 'q':
            *quiet = 1;
            break;
        default:
            print_usage(argv[0]);
            break;
        }
    }
}
