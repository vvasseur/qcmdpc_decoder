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
#include <math.h>

#include "param.h"
#include "threshold.h"

static double lnbino(unsigned n, unsigned t);
static double xlny(double x, double y);
static double lnbinomialpmf(unsigned n, unsigned k, double p);
static double logE(unsigned t, unsigned i);
static double X_val(unsigned t);
static double counters_C0(unsigned S, unsigned t, double x);
static double counters_C1(unsigned S, unsigned t, double x);

static double lnbino(unsigned n, unsigned t) {
    return ((t == 0) || (n == t))
               ? 0.
               : lgamma(n + 1) - lgamma(t + 1) - lgamma(n - t + 1);
}

static double xln1my(double x, double y) {
    return (x == 0.) ? 0. : x * log1p(-y);
}

static double xlny(double x, double y) { return (x == 0.) ? 0. : x * log(y); }

/* Log of the probability mass function of a binomial distribution */
static double lnbinomialpmf(unsigned n, unsigned k, double p) {
    return lnbino(n, k) + xlny(k, p) + xln1my(n - k, p);
}

static double logE(unsigned t, unsigned i) {
    return lnbino(INDEX * BLOCK_WEIGHT, i) +
           lnbino(INDEX * (BLOCK_LENGTH - BLOCK_WEIGHT), t - i) -
           lnbino(INDEX * BLOCK_LENGTH, t);
}

/* X = sum((l - 1) * E_l, l odd) */
static double X_val(unsigned t) {
    unsigned i;
    double x;
    double denom = 0.;

    /* logE(n, w, t, i) decreases fast when 'i' varies.
     * For i >= 10 it is very likely to be negligible. */
    for (x = 0, i = 1; (i < 10) && (i < t); i += 2) {
        x += (i - 1) * exp(logE(t, i));
        denom += exp(logE(t, i));
    }

    return (denom == 0.) ? 0. : x / denom;
}

/* Probability for a bit of the syndrome to be zero, knowing the syndrome
 * weight 'S' and 'X'. */
static double counters_C0(unsigned S, unsigned t, double x) {
    return ((INDEX * BLOCK_WEIGHT - 1) * S - x) / (INDEX * BLOCK_LENGTH - t) /
           BLOCK_WEIGHT;
}

/* Probability for a bit of the syndrome to be non-zero, knowing the syndrome
 * weight 'S' and 'X'. */
static double counters_C1(unsigned S, unsigned t, double x) {
    return (S + x) / t / BLOCK_WEIGHT;
}

unsigned compute_threshold(unsigned S, unsigned t) {
    double p, q;

    double x = X_val(t) * S;
    p = counters_C0(S, t, x);
    q = counters_C1(S, t, x);

    double lq = log(q);
    double lp = log(p);
    double lqbar = log1p(-q);
    double lpbar = log1p(-p);
    double lt = log(t);
    double ltbar = log(INDEX * BLOCK_LENGTH - t);

    unsigned threshold;
    threshold = ceil((BLOCK_WEIGHT * (lqbar - lpbar) + lt - ltbar) /
                     (lp - lq + lqbar - lpbar));

    if (threshold > BLOCK_WEIGHT || q >= 1)
        threshold = BLOCK_WEIGHT;
    if (threshold < (BLOCK_WEIGHT + 1) / 2)
        threshold = (BLOCK_WEIGHT + 1) / 2;

    return threshold;
}

unsigned compute_threshold_alpha(unsigned S, unsigned t, double alpha) {
    double p;

    double x = X_val(t) * S;
    p = counters_C0(S, t, x);
    p = (p <= 1) ? p : 1;

    unsigned threshold;
    if (p >= 1.0) {
        threshold = BLOCK_WEIGHT;
        return threshold;
    }
    else {
        threshold = BLOCK_WEIGHT + 1;
        double diff = 0.;
        do {
            threshold--;
            diff = -exp(lnbinomialpmf(BLOCK_WEIGHT, threshold, p)) +
                   alpha / BLOCK_LENGTH;
        } while (diff >= 0. && threshold > (BLOCK_WEIGHT + 1) / 2);
        threshold = threshold < BLOCK_WEIGHT ? (threshold + 1) : BLOCK_WEIGHT;

        return threshold;
    }
}

#if (ALGO == GRAY_BGF) || (ALGO == GRAY_BGB) || (ALGO == GRAY_B) ||            \
    (ALGO == GRAY_BG)
unsigned compute_threshold_affine(unsigned S) {
    unsigned thr = THRESHOLD_C0 + THRESHOLD_C1 * S;
    return (thr > (BLOCK_WEIGHT + 1) / 2) ? thr : ((BLOCK_WEIGHT + 1) / 2);
}
#endif
