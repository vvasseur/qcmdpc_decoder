#pragma once

#include <stdint.h>

uint64_t random_uint64_t(uint64_t *S0, uint64_t *S1);
int seed_random(uint64_t *S0, uint64_t *S1);
uint64_t random_lim(uint64_t limit, uint64_t *S0, uint64_t *S1);
void jump(uint64_t *S0, uint64_t *S1);

struct PRNG {
    uint64_t s0;
    uint64_t s1;
    uint64_t (*random_lim)(uint64_t limit, uint64_t *s0, uint64_t *s1);
    uint64_t (*random_uint64_t)(uint64_t *s0, uint64_t *s1);
};

typedef struct PRNG *prng_t;
