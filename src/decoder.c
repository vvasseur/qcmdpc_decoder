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
#include "param.h"
#if (ALGO != BP)
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "code.h"
#include "decoder.h"
#include "param.h"
#include "threshold.h"

static bit_t get_counter(decoder_t dec, index_t index, index_t position);
static void flip_column(decoder_t dec, index_t index, index_t position);
static void single_flip(decoder_t dec, index_t index, index_t position);

static bit_t get_counter(decoder_t dec, index_t index, index_t position) {
    bit_t counter = 0;
    index_t offset = position;

    index_t l;
    for (l = 0; l < BLOCK_WEIGHT; ++l) {
        index_t i = offset + dec->H->columns[index][l];
        if (i >= BLOCK_LENGTH) {
            offset -= BLOCK_LENGTH;
            break;
        }
        counter += dec->syndrome->vec[i];
    }
    for (; l < BLOCK_WEIGHT; ++l) {
        index_t i = offset + dec->H->columns[index][l];
        counter += dec->syndrome->vec[i];
    }
    return counter;
}

static void flip_column(decoder_t dec, index_t index, index_t position) {
    index_t offset = position;

    index_t l;
    for (l = 0; l < BLOCK_WEIGHT; ++l) {
        index_t i = offset + dec->H->columns[index][l];
        if (i >= BLOCK_LENGTH) {
            offset -= BLOCK_LENGTH;
            break;
        }
        dec->syndrome->vec[i] ^= 1;
    }
    for (; l < BLOCK_WEIGHT; ++l) {
        index_t i = offset + dec->H->columns[index][l];
        dec->syndrome->vec[i] ^= 1;
    }
}

static void single_flip(decoder_t dec, index_t index, index_t position) {
    bit_t counter = get_counter(dec, index, position);
    flip_column(dec, index, position);
    dec->bits[index][position] ^= 1;
    dec->syndrome->weight += BLOCK_WEIGHT - 2 * counter;
    dec->e->weight +=
        2 * (dec->bits[index][position] ^ dec->e->vec[index][position]) - 1;
}

void init_decoder(decoder_t dec, code_t *H, e_t *e, syndrome_t *syndrome) {
    dec->H = H;
    dec->e = e;
    dec->syndrome = syndrome;
}

void reset_decoder(decoder_t dec) {
    memset(dec->bits, 0, INDEX * BLOCK_LENGTH * sizeof(bit_t));
#if (ALGO == BACKFLIP) || (ALGO == BACKFLIP2)
    dec->fl.first = -1;
    dec->fl.length = 0;
#endif
    dec->iter = 0;
}

#if (ALGO == CLASSIC)
int qcmdpc_decode(decoder_t dec, int max_iter) {
    dec->blocked = false;
    while (dec->iter < max_iter && dec->syndrome->weight != SYNDROME_STOP &&
           !dec->blocked) {
        ++dec->iter;

        compute_counters(dec->counters, dec->syndrome->vec, dec->H);

        unsigned threshold =
            compute_threshold(dec->syndrome->weight, dec->e->weight);

        dec->blocked = true;
        for (index_t k = 0; k < INDEX; ++k)
            for (index_t j = 0; j < BLOCK_LENGTH; ++j)
                if (dec->counters[k][j] >= threshold) {
                    single_flip(dec, k, j);
                    dec->blocked = false;
                }
    }

    return !dec->e->weight;
}
#endif

#if (ALGO == BACKFLIP) || (ALGO == BACKFLIP2)
/* Double linked list implementation, used to store the flip list. */
static void fl_add(fl_t *fl, index_t pos) {
    fl->next[pos] = fl->first;
    fl->prev[pos] = -1;
    if (fl->first != -1)
        fl->prev[fl->first] = pos;
    fl->first = pos;
    ++fl->length;
}

static void fl_remove(fl_t *fl, index_t pos) {
    index_t next = fl->next[pos];
    index_t prev = fl->prev[pos];
    if (next != -1) {
        fl->prev[next] = prev;
    }
    if (prev != -1) {
        fl->next[prev] = next;
    }
    else {
        fl->first = next;
    }
    --fl->length;
}

#if (ALGO == BACKFLIP)
static inline int affine_ttl(int diff) {
    int ttl = (int)(TTL_C0 + TTL_C1 * diff);

    ttl = (ttl < 1) ? 1 : ttl;
    return (ttl > TTL_SATURATE) ? TTL_SATURATE : ttl;
}
#endif

int qcmdpc_decode(decoder_t dec, int max_iter) {
    unsigned threshold = 0;
#if (ALGO == BACKFLIP2)
    unsigned threshold2 = 0;
    unsigned threshold3 = 0;
    unsigned threshold4 = 0;
    unsigned threshold5 = 0;
#endif
    /* Only recompute the threshold when necessary */
    dec->blocked = false;
    while (dec->iter < max_iter && dec->syndrome->weight != SYNDROME_STOP) {
        ++dec->iter;
        compute_counters(dec->counters, dec->syndrome->vec, dec->H);

        if (!dec->blocked) {
            int t = (ERROR_WEIGHT > dec->fl.length)
                        ? ERROR_WEIGHT - dec->fl.length
                        : 1;
#if (ALGO == BACKFLIP)
            threshold = compute_threshold(dec->syndrome->weight, t);
#else // (ALGO == BACKFLIP2)
            threshold =
                compute_threshold_alpha(dec->syndrome->weight, t, THRESHOLD_A0);
            threshold2 =
                compute_threshold_alpha(dec->syndrome->weight, t, THRESHOLD_A1);
            threshold3 =
                compute_threshold_alpha(dec->syndrome->weight, t, THRESHOLD_A2);
            threshold4 =
                compute_threshold_alpha(dec->syndrome->weight, t, THRESHOLD_A3);
            threshold5 =
                compute_threshold_alpha(dec->syndrome->weight, t, THRESHOLD_A4);
#endif
        }

        dec->blocked = true;
        for (index_t k = 0; k < INDEX; ++k) {
            for (index_t j = 0; j < BLOCK_LENGTH; ++j) {
                if (dec->counters[k][j] >= threshold) {
                    if (dec->bits[k][j])
                        fl_remove(&dec->fl, k * BLOCK_LENGTH + j);
                    else {
#if (ALGO == BACKFLIP)
                        uint8_t ttl =
                            affine_ttl(dec->counters[k][j] - threshold);
#else // (ALGO == BACKFLIP2)
                        uint8_t ttl =
                            (dec->counters[k][j] < threshold2)
                                ? 1
                                : ((dec->counters[k][j] < threshold3)
                                       ? 2
                                       : ((dec->counters[k][j] < threshold4)
                                              ? 3
                                              : ((dec->counters[k][j] <
                                                  threshold5)
                                                     ? 4
                                                     : 5)));
#endif

                        fl_add(&dec->fl, k * BLOCK_LENGTH + j);
                        dec->fl.tod[k * BLOCK_LENGTH + j] =
                            (dec->iter + ttl) % (TTL_SATURATE + 1);
                    }
                    single_flip(dec, k, j);
                    dec->blocked = false;
                }
            }
        }
        if (dec->syndrome->weight != SYNDROME_STOP && dec->fl.length) {
            uint8_t current_iter = dec->iter % (TTL_SATURATE + 1);
            index_t fl_pos = dec->fl.first;
            while (fl_pos != -1) {
                if (dec->fl.tod[fl_pos] == current_iter) {
                    index_t k = (fl_pos < BLOCK_LENGTH) ? 0 : 1;
                    index_t j =
                        fl_pos - (fl_pos < BLOCK_LENGTH ? 0 : BLOCK_LENGTH);

                    fl_remove(&dec->fl, fl_pos);

                    single_flip(dec, k, j);
                    dec->blocked = false;
                }
                fl_pos = dec->fl.next[fl_pos];
            }
        }
    }

    return !dec->e->weight;
}
#endif

#if (ALGO == GRAY_BGF) || (ALGO == GRAY_BGB) || (ALGO == GRAY_B) ||            \
    (ALGO == GRAY_BG)
int qcmdpc_decode(decoder_t dec, int max_iter) {
    unsigned threshold = 0;
    /* Only recompute the threshold when necessary */
    dec->blocked = false;
    while (dec->iter < max_iter && dec->syndrome->weight != SYNDROME_STOP &&
           !dec->blocked) {
        ++dec->iter;
        compute_counters(dec->counters, dec->syndrome->vec, dec->H);

        if (!dec->blocked)
            threshold = compute_threshold_affine(dec->syndrome->weight);

        dec->blocked = true;

        dec->gray.length = 0;
        dec->black.length = 0;
        for (index_t k = 0; k < INDEX; ++k)
            for (index_t j = 0; j < BLOCK_LENGTH; ++j) {
                if (dec->counters[k][j] >= threshold) {
                    single_flip(dec, k, j);
                    dec->blocked = false;

                    index_t curr = dec->black.length++;
                    dec->black.index[curr] = k;
                    dec->black.position[curr] = j;
                }
                else if (dec->counters[k][j] + (unsigned)GRAY_DELTA >=
                         threshold) {
                    index_t curr = dec->gray.length++;
                    dec->gray.index[curr] = k;
                    dec->gray.position[curr] = j;
                }
            }
/* We count each black or gray step as an iteration. */
#if (ALGO == GRAY_BGF)
        if (dec->iter < 2)
#endif
        {
            ++dec->iter;
            for (index_t i = 0; i < dec->black.length; ++i) {
                index_t k = dec->black.index[i];
                index_t j = dec->black.position[i];
                dec->black.counter[i] = get_counter(dec, k, j);
            }
            for (index_t i = 0; i < dec->black.length; ++i) {
                if (dec->black.counter[i] >= (BLOCK_WEIGHT + 1) / 2 + 1) {
                    index_t k = dec->black.index[i];
                    index_t j = dec->black.position[i];

                    single_flip(dec, k, j);
                    dec->blocked = false;
                }
            }
#if (ALGO == GRAY_BGB)
            if (dec->iter < 3)
#elif (ALGO == GRAY_B)
            if (0)
#endif
            {
                ++dec->iter;
                for (index_t i = 0; i < dec->gray.length; ++i) {
                    index_t k = dec->gray.index[i];
                    index_t j = dec->gray.position[i];
                    dec->gray.counter[i] = get_counter(dec, k, j);
                }
                for (index_t i = 0; i < dec->gray.length; ++i) {
                    if (dec->gray.counter[i] >= (BLOCK_WEIGHT + 1) / 2 + 1) {
                        index_t k = dec->gray.index[i];
                        index_t j = dec->gray.position[i];

                        single_flip(dec, k, j);
                        dec->blocked = false;
                    }
                }
            }
        }
    }

    return !dec->e->weight;
}
#endif

#if (ALGO == SBS)
int qcmdpc_decode(decoder_t dec, int max_iter, prng_t prng) {
#elif (ALGO == SORT)
/* Step-by-step is used as the final step of the sorted gray decoder. */
static int qcmdpc_decode_sbs(decoder_t dec, int max_iter, prng_t prng) {
#endif
#if (ALGO == SBS) || (ALGO == SORT)
    unsigned threshold = 0;
    /* Only recompute the threshold when necessary */
    dec->blocked = false;
    unsigned long missed = 0;
    while (dec->iter < max_iter && dec->syndrome->weight != SYNDROME_STOP) {
        ++dec->iter;
        if (missed > INDEX * BLOCK_LENGTH) {
            compute_counters(dec->counters, dec->syndrome->vec, dec->H);
            bool found = false;
            for (index_t k = 0; k < INDEX && !found; ++k) {
                for (index_t j = 0; j < BLOCK_LENGTH && !found; ++j) {
                    found |= (dec->counters[k][j] >= threshold);
                }
            }
            if (!found)
                break;
            missed = 0;
        }
        if (!dec->blocked)
            threshold =
                compute_threshold(dec->syndrome->weight, dec->e->weight);

        dec->blocked = true;
        /* Randomly pick a 1 in the syndrome */
        int i;
        do {
            i = prng->random_lim(BLOCK_LENGTH, prng->s);
        } while (!dec->syndrome->vec[i]);
        int k;
        k = prng->random_lim(INDEX, prng->s);
        int l;
        l = prng->random_lim(BLOCK_WEIGHT, prng->s);

        int j = i + dec->H->rows[k][l];
        j = (j >= BLOCK_LENGTH) ? (j - BLOCK_LENGTH) : j;

        bit_t counter = get_counter(dec, k, j);
        if (counter >= threshold) {
            single_flip(dec, k, j);
            dec->blocked = false;
        }
        else
            ++missed;
    }

    return !dec->e->weight;
}
#endif

#if (ALGO == SORT)
/* Heap implementation to extract the GRAY_SIZE greater counters. */
#define LCHILD(i) 2 * (i) + 1
#define RCHILD(i) 2 * (i) + 2
#define PARENT(i) ((i)-1) / 2

void swap(pos_counter_t *array, index_t i1, index_t i2) {
    pos_counter_t tmp = array[i1];
    array[i1] = array[i2];
    array[i2] = tmp;
}

void sift_down(pos_counter_t *array, index_t start, index_t end) {
    index_t root = start;

    while (LCHILD(root) <= end) {
        index_t child = LCHILD(root);
        index_t to_swap = root;

        if (array[to_swap].counter > array[child].counter)
            to_swap = child;
        if (child + 1 <= end &&
            array[to_swap].counter > array[child + 1].counter)
            to_swap = child + 1;
        if (to_swap == root)
            return;
        else {
            swap(array, root, to_swap);
            root = to_swap;
        }
    }
}

void heapify(pos_counter_t *array, index_t size) {
    for (index_t start = PARENT(size - 1); start >= 0; --start)
        sift_down(array, start, size - 1);
}

void insert(pos_counter_t *array, pos_counter_t pc, index_t size) {
    if (array[0].counter < pc.counter) {
        array[0] = pc;
        sift_down(array, 0, size - 1);
    }
}

static int cmp_pos_counter_t(const void *a, const void *b) {
    return (((pos_counter_t *)b)->counter - ((pos_counter_t *)a)->counter);
}

static void sort_counters(pos_counter_t sorted_counters[GRAY_SIZE],
                          const bit_t counters[INDEX][2 * SIZE_AVX],
                          index_t size) {
    index_t j = 0;
    for (index_t k = 0; k < INDEX; ++k)
        for (index_t i = 0; i < BLOCK_LENGTH; ++i) {
            pos_counter_t pc;
            pc.index = k;
            pc.position = i;
            pc.counter = counters[k][i];

            if (j < size)
                sorted_counters[j] = pc;
            else {
                if (j == size)
                    heapify(sorted_counters, size);
                insert(sorted_counters, pc, size);
            }
            ++j;
        }
    qsort(sorted_counters, size, sizeof(pos_counter_t), cmp_pos_counter_t);
}

int qcmdpc_decode(decoder_t dec, int max_iter, prng_t prng) {
    compute_counters(dec->counters, dec->syndrome->vec, dec->H);
    sort_counters(dec->sorted_counters, dec->counters, GRAY_SIZE);

    const unsigned threshold = (BLOCK_WEIGHT + 1) / 2;

    index_t i = 0;
    while (dec->iter < INDEX * BLOCK_LENGTH && dec->iter < max_iter) {
        ++dec->iter;
        index_t k = dec->sorted_counters[i].index;
        index_t j = dec->sorted_counters[i].position;

        bit_t counter = get_counter(dec, k, j);
        if (counter >= threshold)
            single_flip(dec, k, j);
        i = (++i == GRAY_SIZE) ? 0 : i;
    }

    qcmdpc_decode_sbs(dec, max_iter, prng);

    return !dec->e->weight;
}
#endif
#endif
