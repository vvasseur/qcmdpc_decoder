/*
   Written in 2019 by David Blackman and Sebastiano Vigna (vigna@acm.org)

   To the extent possible under law, the author has dedicated all copyright and
   related and neighboring rights to this software to the public domain
   worldwide. This software is distributed without any warranty.

   See <http://creativecommons.org/publicdomain/zero/1.0/>.
*/
#include <stdint.h>
#include <stdio.h>

#include "xoshiro256plusplus.h"

/*
   This is xoshiro256++ 1.0, one of our all-purpose, rock-solid generators. It
   has excellent (sub-ns) speed, a state (256 bits) that is large enough for any
   parallel application, and it passes all tests we are aware of.

   For generating just floating-point numbers, xoshiro256+ is even faster.

   The state must be seeded so that it is not everywhere zero. If you have a
   64-bit seed, we suggest to seed a splitmix64 generator and use its output to
   fill s.
*/

static inline uint64_t rotl(const uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

uint64_t random_uint64_t(uint64_t *s) {
    const uint64_t result = rotl(s[0] + s[3], 23) + s[0];

    const uint64_t t = s[1] << 17;

    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];

    s[2] ^= t;

    s[3] = rotl(s[3], 45);

    return result;
}

int seed_random(uint64_t *s) {
    FILE *urandom_fp;

    urandom_fp = fopen("/dev/urandom", "r");
    if (urandom_fp == NULL)
        return 0;
    if (fread(s, 8, 4, urandom_fp) != 4)
        return 0;
    fclose(urandom_fp);

    return 1;
}

/* See
 * <https://lemire.me/blog/2019/06/06/nearly-divisionless-random-integer-generation-on-various-systems/>.
 */
uint64_t random_lim(uint64_t limit, uint64_t *s) {
    uint64_t x = random_uint64_t(s);
    __uint128_t m = (__uint128_t)x * (__uint128_t)limit;
    uint64_t l = (uint64_t)m;
    if (l < limit) {
        uint64_t t = -limit % limit;
        while (l < t) {
            x = random_uint64_t(s);
            m = (__uint128_t)x * (__uint128_t)limit;
            l = (uint64_t)m;
        }
    }
    return m >> 64;
}

/* This is the jump function for the generator. It is equivalent to 2^128 calls
 * to next(); it can be used to generate 2^128 non-overlapping subsequences for
 * parallel computations. */
void jump(uint64_t *s) {
    static const uint64_t JUMP[] = {0x180ec6d33cfd0aba, 0xd5a61266f0c9392c,
                                    0xa9582618e03fc9aa, 0x39abdc4529b1661c};

    uint64_t s0 = 0;
    uint64_t s1 = 0;
    uint64_t s2 = 0;
    uint64_t s3 = 0;
    for (size_t i = 0; i < sizeof(JUMP) / sizeof(*JUMP); i++)
        for (size_t b = 0; b < 64; b++) {
            if (JUMP[i] & UINT64_C(1) << b) {
                s0 ^= s[0];
                s1 ^= s[1];
                s2 ^= s[2];
                s3 ^= s[3];
            }
            random_uint64_t(s);
        }

    s[0] = s0;
    s[1] = s1;
    s[2] = s2;
    s[3] = s3;
}
