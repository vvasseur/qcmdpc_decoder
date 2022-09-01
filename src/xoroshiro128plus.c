/*  Written in 2016-2018 by David Blackman and Sebastiano Vigna (vigna@acm.org)

    To the extent possible under law, the author has dedicated all copyright
    and related and neighboring rights to this software to the public domain
    worldwide. This software is distributed without any warranty.

    See <http://creativecommons.org/publicdomain/zero/1.0/>. */

#include <stdint.h>
#include <stdio.h>

#include "xoroshiro128plus.h"

/* This is the successor to xorshift128+. It is the fastest full-period
   generator passing BigCrush without systematic failures, but due to the
   relatively short period it is acceptable only for applications with a
   mild amount of parallelism; otherwise, use a xorshift1024* generator.

   Beside passing BigCrush, this generator passes the PractRand test suite
   up to (and included) 16TB, with the exception of binary rank tests,
   which fail due to the lowest bit being an LFSR; all other bits pass all
   tests. We suggest to use a sign test to extract a random Boolean value.

   Note that the generator uses a simulated rotate operation, which most C
   compilers will turn into a single instruction. In Java, you can use
   Long.rotateLeft(). In languages that do not make low-level rotation
   instructions accessible xorshift128+ could be faster.

   The state must be seeded so that it is not everywhere zero. If you have
   a 64-bit seed, we suggest to seed a splitmix64 generator and use its
   output to fill s. */

static inline uint64_t rotl(const uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

uint64_t random_uint64_t(uint64_t *S0, uint64_t *S1) {
    const uint64_t s0 = *S0;
    uint64_t s1 = *S1;
    const uint64_t result = s0 + s1;

    s1 ^= s0;
    *S0 = rotl(s0, 24) ^ s1 ^ (s1 << 16); // a, b
    *S1 = rotl(s1, 37);                   // c

    return result;
}

int seed_random(uint64_t *S0, uint64_t *S1) {
    uint64_t new_s[2];
    FILE *urandom_fp;

    urandom_fp = fopen("/dev/urandom", "r");
    if (urandom_fp == NULL)
        return 0;
    if (fread(&new_s, 8, 2, urandom_fp) != 2)
        return 0;
    fclose(urandom_fp);

    *S0 = new_s[0];
    *S1 = new_s[1];

    return 1;
}

uint64_t random_lim(uint64_t limit, uint64_t *S0, uint64_t *S1) {
    uint64_t divisor = 0xffffffffffffffffUL / (limit + 1);
    uint64_t retval;

    do {
        retval = random_uint64_t(S0, S1) / divisor;
    } while (retval > limit);

    return retval;
}

/* This is the jump function for the generator. It is equivalent
   to 2^64 calls to next(); it can be used to generate 2^64
   non-overlapping subsequences for parallel computations. */

void jump(uint64_t *S0, uint64_t *S1) {
    static const uint64_t JUMP[] = {0xdf900294d8f554a5, 0x170865df4b3201fc};

    uint64_t s0 = 0;
    uint64_t s1 = 0;
    for (size_t i = 0; i < sizeof(JUMP) / sizeof(*JUMP); i++)
        for (int b = 0; b < 64; b++) {
            if (JUMP[i] & UINT64_C(1) << b) {
                s0 ^= *S0;
                s1 ^= *S1;
            }
            random_uint64_t(S0, S1);
        }

    *S0 = s0;
    *S1 = s1;
}
