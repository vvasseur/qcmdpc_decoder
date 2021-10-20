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
#ifndef PARAM_H
#define PARAM_H

/* In seconds */
#define TIME_BETWEEN_PRINTS 60

#define CLASSIC 0
#define BACKFLIP 1
#define BACKFLIP2 2
#define SBS 3
#define GRAY_B 4
#define GRAY_BGF 5
#define GRAY_BGB 6
#define GRAY_BG 7
#define BP 8
#define SORT 9

#if !defined(PRESET_CPA) && !defined(PRESET_CCA) &&                            \
    !(defined(INDEX) && defined(BLOCK_LENGTH) && defined(BLOCK_WEIGHT) &&      \
      defined(ERROR_WEIGHT))
#define PRESET_CCA 128
#endif
#ifndef OUROBOROS
#define OUROBOROS 0
#endif

#ifdef PRESET_CPA
#if PRESET_CPA == 128 && !OUROBOROS
#ifndef INDEX
#define INDEX 2
#endif
#ifndef BLOCK_LENGTH
#define BLOCK_LENGTH 10163
#endif
#ifndef BLOCK_WEIGHT
#define BLOCK_WEIGHT 71
#endif
#ifndef ERROR_WEIGHT
#define ERROR_WEIGHT 134
#endif
#elif PRESET_CPA == 192 && !OUROBOROS
#ifndef INDEX
#define INDEX 2
#endif
#ifndef BLOCK_LENGTH
#define BLOCK_LENGTH 19853
#endif
#ifndef BLOCK_WEIGHT
#define BLOCK_WEIGHT 103
#endif
#ifndef ERROR_WEIGHT
#define ERROR_WEIGHT 199
#endif
#elif PRESET_CPA == 256 && !OUROBOROS
#ifndef INDEX
#define INDEX 2
#endif
#ifndef BLOCK_LENGTH
#define BLOCK_LENGTH 32749
#endif
#ifndef BLOCK_WEIGHT
#define BLOCK_WEIGHT 137
#endif
#ifndef ERROR_WEIGHT
#define ERROR_WEIGHT 264
#endif
#elif PRESET_CPA == 128 && OUROBOROS
#ifndef INDEX
#define INDEX 2
#endif
#ifndef BLOCK_LENGTH
#define BLOCK_LENGTH 11027
#endif
#ifndef BLOCK_WEIGHT
#define BLOCK_WEIGHT 67
#endif
#ifndef ERROR_WEIGHT
#define ERROR_WEIGHT 154
#endif
#elif PRESET_CPA == 192 && OUROBOROS
#ifndef INDEX
#define INDEX 2
#endif
#ifndef BLOCK_LENGTH
#define BLOCK_LENGTH 21683
#endif
#ifndef BLOCK_WEIGHT
#define BLOCK_WEIGHT 99
#endif
#ifndef ERROR_WEIGHT
#define ERROR_WEIGHT 226
#endif
#elif PRESET_CPA == 256 && OUROBOROS
#ifndef INDEX
#define INDEX 2
#endif
#ifndef BLOCK_LENGTH
#define BLOCK_LENGTH 36131
#endif
#ifndef BLOCK_WEIGHT
#define BLOCK_WEIGHT 133
#endif
#ifndef ERROR_WEIGHT
#define ERROR_WEIGHT 300
#endif
#endif
#endif

#ifdef PRESET_CCA
#if PRESET_CCA == 128
#ifndef INDEX
#define INDEX 2
#endif
#ifndef BLOCK_LENGTH
#define BLOCK_LENGTH 12323
#endif
#ifndef BLOCK_WEIGHT
#define BLOCK_WEIGHT 71
#endif
#ifndef ERROR_WEIGHT
#define ERROR_WEIGHT 134
#endif
#elif PRESET_CCA == 192
#ifndef INDEX
#define INDEX 2
#endif
#ifndef BLOCK_LENGTH
#define BLOCK_LENGTH 24659
#endif
#ifndef BLOCK_WEIGHT
#define BLOCK_WEIGHT 103
#endif
#ifndef ERROR_WEIGHT
#define ERROR_WEIGHT 199
#endif
#elif PRESET_CCA == 256
#ifndef INDEX
#define INDEX 2
#endif
#ifndef BLOCK_LENGTH
#define BLOCK_LENGTH 40973
#endif
#ifndef BLOCK_WEIGHT
#define BLOCK_WEIGHT 137
#endif
#ifndef ERROR_WEIGHT
#define ERROR_WEIGHT 264
#endif
#endif
#endif

#if (!OUROBOROS && BLOCK_WEIGHT == 71 && ERROR_WEIGHT == 134)
#ifndef THRESHOLD_C0
#define THRESHOLD_C0 13.530
#endif
#ifndef THRESHOLD_C1
#define THRESHOLD_C1 0.0069722
#endif
#ifndef TTL_C0
#define TTL_C0 1.1
#endif
#ifndef TTL_C1
#define TTL_C1 0.45
#endif
#elif (!OUROBOROS && BLOCK_WEIGHT == 103 && ERROR_WEIGHT == 199)
#ifndef THRESHOLD_C0
#define THRESHOLD_C0 15.2588
#endif
#ifndef THRESHOLD_C1
#define THRESHOLD_C1 0.005265
#endif
#ifndef TTL_C0
#define TTL_C0 1.41
#endif
#ifndef TTL_C1
#define TTL_C1 0.36
#endif
#elif (!OUROBOROS && BLOCK_WEIGHT == 137 && ERROR_WEIGHT == 264)
#ifndef THRESHOLD_C0
#define THRESHOLD_C0 17.8785
#endif
#ifndef THRESHOLD_C1
#define THRESHOLD_C1 0.00402312
#endif
#ifndef TTL_C0
#define TTL_C0 1.
#endif
#ifndef TTL_C1
#define TTL_C1 0.45
#endif
#elif (OUROBOROS && BLOCK_WEIGHT == 67 && ERROR_WEIGHT == 154)
#ifndef THRESHOLD_C0
#define THRESHOLD_C0 13.209
#endif
#ifndef THRESHOLD_C1
#define THRESHOLD_C1 0.0060515
#endif
#ifndef TTL_C0
#define TTL_C0 1.16
#endif
#ifndef TTL_C1
#define TTL_C1 0.46
#endif
#elif (OUROBOROS && BLOCK_WEIGHT == 99 && ERROR_WEIGHT == 226)
#ifndef THRESHOLD_C0
#define THRESHOLD_C0 15.561
#endif
#ifndef THRESHOLD_C1
#define THRESHOLD_C1 0.0046692
#endif
#ifndef TTL_C0
#define TTL_C0 1.4
#endif
#ifndef TTL_C1
#define TTL_C1 0.4
#endif
#elif (OUROBOROS && BLOCK_WEIGHT == 133 && ERROR_WEIGHT == 300)
#ifndef THRESHOLD_C0
#define THRESHOLD_C0 17.061
#endif
#ifndef THRESHOLD_C1
#define THRESHOLD_C1 0.0038460
#endif
#ifndef TTL_C0
#define TTL_C0 0.9
#endif
#ifndef TTL_C1
#define TTL_C1 0.44
#endif
#endif

#if OUROBOROS && !defined(SYNDROME_STOP)
#define SYNDROME_STOP ((ERROR_WEIGHT) / 2)
#else
#define SYNDROME_STOP 0
#endif

#ifndef WEAK
#define WEAK 0
#endif
#ifndef WEAK_P
#define WEAK_P 0
#endif

#ifndef ERROR_FLOOR
#define ERROR_FLOOR 0
#endif
#ifndef ERROR_FLOOR_P
#define ERROR_FLOOR_P 0
#endif

#if !defined(ALGO)
#define ALGO GRAY_BGF
#endif

#ifndef BP_SCALE
#define BP_SCALE 0.4
#endif
#ifndef BP_SATURATE
#define BP_SATURATE 1000.
#endif

#ifndef GRAY_DELTA
#define GRAY_DELTA 3
#endif

#ifndef GRAY_SIZE
#define GRAY_SIZE (((INDEX) * (BLOCK_LENGTH)) / 10)
#endif

#ifndef THRESHOLD_A0
#define THRESHOLD_A0 4.
#endif
#ifndef THRESHOLD_A1
#define THRESHOLD_A1 1.
#endif
#ifndef THRESHOLD_A2
#define THRESHOLD_A2 0.25
#endif
#ifndef THRESHOLD_A3
#define THRESHOLD_A3 0.0625
#endif
#ifndef THRESHOLD_A4
#define THRESHOLD_A4 0.015625
#endif
#ifndef THRESHOLD_C0
#define THRESHOLD_C0 13.530
#endif
#ifndef THRESHOLD_C1
#define THRESHOLD_C1 0.0069722
#endif

#ifndef TTL_C0
#define TTL_C0 1.1
#endif
#ifndef TTL_C1
#define TTL_C1 0.45
#endif
#ifndef TTL_SATURATE
#define TTL_SATURATE 5
#endif

#if INDEX != 2
#error "INDEX != 2: Not implemented"
#endif
#if BLOCK_WEIGHT > 255
#error "BLOCK_WEIGHT > 255: Not implemented"
#endif
#if BLOCK_LENGTH > 65536
#error "BLOCK_LENGTH > 65536: Not implemented"
#endif
#if ALGO == BP && OUROBOROS
#error "Ouroboros with belief propagation decoding: Not implemented"
#endif
#endif
