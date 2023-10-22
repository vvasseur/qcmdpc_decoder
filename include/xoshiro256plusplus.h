#pragma once
#include <stdint.h>

uint64_t random_uint64_t(uint64_t *s);
int seed_random(uint64_t *s);
uint64_t random_lim(uint64_t limit, uint64_t *s);
void jump(uint64_t *s);

struct PRNG {
    uint64_t s[4];
    uint64_t (*random_lim)(uint64_t limit, uint64_t *s);
    uint64_t (*random_uint64_t)(uint64_t *s);
};

typedef struct PRNG *prng_t;
