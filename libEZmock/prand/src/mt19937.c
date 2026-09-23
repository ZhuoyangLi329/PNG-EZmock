/*******************************************************************************
* mt19937.c: this file is part of the prand library.
 
* prand: parallel random number generator.

* Github repository:
        https://github.com/cheng-zhao/prand

* Copyright (c) 2020 Cheng Zhao <zhaocheng03@gmail.com>
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.

*******************************************************************************/

#include "mt19937.h"
#include "mt19937_jump.h"
#include <stdlib.h>
#include <string.h>

/*******************************************************************************
  Implementation of the Mersenne Twister 19937 random number generator.
  ref: https://doi.org/10.1145%2F272991.272995

  The initialisation with a seed is performed following the 2002 version:
  http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/emt19937ar.html

  This variant has been written from scratch, but produces the same
  sequences as the version provided by the authors of the algorithm:
  http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
*******************************************************************************/

/*============================================================================*\
                            Definitions of constants
\*============================================================================*/

#ifndef MT19937_N
  #define N               624
#else
  #define N               MT19937_N
#endif

#define M               397
#define K               19937
#define MA              0x9908b0dfUL
#define UPPER_MASK(x)   (0x80000000UL & x)      /* most significant w-r bits */
#define LOWER_MASK(x)   (0x7fffffffUL & x)      /* least significant r bits */

/* Normalisation for sampling a float-point number in the range [0,1). */
#define NORM            0x1p-32                 /* 2^{-32} */
/* Normalisation for sampling a float-point number in the range (0,1). */
#define NORM_POS        0x1.fffffffep-33        /* 1 / (2^{32} + 1) */

#define DEFAULT_SEED    1

#ifndef MT19937_UNROLL
  #define MAGIC(y)      (((y) & 1) ? MA : 0)
#else
/******************************************************************************
  The following trick is taken from
  https://github.com/cslarsen/mersenne-twister
  which is distributed under the modified BSD license.
  It is supposed to work faster with `-ftree-vectorize`, especially for
  machines with SSE/AVX.
******************************************************************************/
  #define MAGIC(y)      (((((int32_t)y) << 31) >> 31) & MA)
#endif

/*============================================================================*\
                            Definition of the state
\*============================================================================*/

typedef struct {
  uint32_t mt[N];
  int idx;
} mt19937_state_t;


/*============================================================================*\
                     Functions for random number generation
\*============================================================================*/

/******************************************************************************
Function `mt19937_seed`:
  Initialise the state with an integer.
Arguments:
  * `state`:    the state to be intialised;
  * `seed`:     a positive integer for the initialisation.
******************************************************************************/
static void mt19937_seed(void *state, uint64_t seed) {
  mt19937_state_t *stat = (mt19937_state_t *) state;
  int i;
  /* Initialise the state with LCG.
   * Bit truncation with 0xffffffffUL is not necessary for 32-bit states. */
  stat->mt[0] = seed;   /* validation of seed is done in `mt19937_init' */
  for (i = 1; i < N; i++)
    stat->mt[i] = 1812433253UL * (stat->mt[i-1] ^ (stat->mt[i-1] >> 30)) + i;

  stat->idx = i;
}

/******************************************************************************
Function `mt19937_get`:
  Generate an integer and update the state.
Arguments:
  * `state`:    the state for the generator.
Return:
  A pseudo-random integer.
******************************************************************************/
static uint64_t mt19937_get(void *state) {
  mt19937_state_t *stat = (mt19937_state_t *) state;
  uint32_t y;

  if (stat->idx >= N) {         /* generate N words at one time */
    int k;
    for (k = 0; k < N - M; k++) {
      y = UPPER_MASK(stat->mt[k]) | LOWER_MASK(stat->mt[k+1]);
      stat->mt[k] = stat->mt[k+M] ^ (y >> 1) ^ MAGIC(y);
    }
    for (; k < N - 1; k++) {
      y = UPPER_MASK(stat->mt[k]) | LOWER_MASK(stat->mt[k+1]);
      stat->mt[k] = stat->mt[k+M-N] ^ (y >> 1) ^ MAGIC(y);
    }
    y = UPPER_MASK(stat->mt[N-1]) | LOWER_MASK(stat->mt[0]);
    stat->mt[N-1] = stat->mt[M-1] ^ (y >> 1) ^ MAGIC(y);
    stat->idx = 0;
  }

  /* tempering */
  y = stat->mt[stat->idx];

  y ^= (y >> 11);
  y ^= (y << 7) & 0x9d2c5680UL;
  y ^= (y << 15) & 0xefc60000UL;
  y ^= (y >> 18);

  stat->idx += 1;
  return y;
}

/******************************************************************************
Function `mt19937_get_double`:
  Generate a double-precision floating-point number in the range [0,1).
Arguments:
  * `state`:    the state for the generator.
Return:
  A pseudo-random floating-point number.
******************************************************************************/
static double mt19937_get_double(void *state) {
  return mt19937_get(state) * NORM;
}

/******************************************************************************
Function `mt19937_get_double_pos`:
  Generate a double-precision floating-point number in the range (0,1).
Arguments:
  * `state`:    the state for the generator.
Return:
  A pseudo-random floating-point number.
******************************************************************************/
static double mt19937_get_double_pos(void *state) {
  return (mt19937_get(state) + 1) * NORM_POS;
}


/*============================================================================*\
                         Functions for multiple streams
\*============================================================================*/

/******************************************************************************
  ref: https://doi.org/10.1007/978-3-540-85912-3_26

  The algorithm for reconstructing the state from a polynomial is taken from
  the boost library (https://www.boost.org, random/mersenne_twister.hpp),
  which is distributed under the Boost Software License:

  Copyright Jens Maurer 2000-2001
  Copyright Steven Watanabe 2010
  Distributed under the Boost Software License, Version 1.0.
  (See copy at https://www.boost.org/LICENSE_1_0.txt)

******************************************************************************/

/******************************************************************************
Function `copy_state`:
  Copy the state.
Arguments:
  * `dst`:      the destination;
  * `src`:      the state to be copied.
******************************************************************************/
static inline void copy_state(mt19937_state_t *dst,
    const mt19937_state_t *src) {
  if (dst == src) return;
  for (int i = 0; i < N; i++) dst->mt[i] = src->mt[i];
  dst->idx = src->idx;
}

/******************************************************************************
Function `next_state`:
  Update the state and return the next element.
Arguments:
  * `s`:        the state for the generator.
Return:
  The next element (indicated by the state pointer) of the state array.
******************************************************************************/
static uint32_t next_state(mt19937_state_t *s) {
  if (s->idx >= N) {            /* generate N words at one time */
    int k;
    uint32_t y;
    for (k = 0; k < N - M; k++) {
      y = UPPER_MASK(s->mt[k]) | LOWER_MASK(s->mt[k+1]);
      s->mt[k] = s->mt[k+M] ^ (y >> 1) ^ MAGIC(y);
    }
    for (; k < N - 1; k++) {
      y = UPPER_MASK(s->mt[k]) | LOWER_MASK(s->mt[k+1]);
      s->mt[k] = s->mt[k+M-N] ^ (y >> 1) ^ MAGIC(y);
    }
    y = UPPER_MASK(s->mt[N-1]) | LOWER_MASK(s->mt[0]);
    s->mt[N-1] = s->mt[M-1] ^ (y >> 1) ^ MAGIC(y);
    s->idx = 0;
  }
  return s->mt[s->idx++];
}

/******************************************************************************
Function `recover_state`:
  Reconstruct the state from a polynomial.
Arguments:
  * `s`:        the state to be reconstructed;
  * `poly`:     the polynomial.
******************************************************************************/
static void recover_state(mt19937_state_t *s, const uint32_t *poly) {
  uint32_t y0 = 0;
  for (int i = K - N + 1; i <= K; i++) s->mt[i % N] = COEF(poly, i);
  for (int i = K + 1; i >= N - 1; i--) {
    uint32_t y1 = s->mt[i % N] ^ s->mt[(i + M) % N];
    if (COEF(poly, i - N + 1))
      y1 = ((y1 ^ MA) << 1) | UINT32_C(1);
    else
      y1 <<= 1;
    s->mt[(i + 1) % N] = UPPER_MASK(y0) | LOWER_MASK(y1);
    y0 = y1;
  }
  s->idx = 0;
}

/******************************************************************************
Function `get_poly`:
  Compute the polynomial for a given skipping step from pre-computed values.
Arguments:
  * `step`:     the number of steps to be skipped.
Return:
  The pointer to the evaluated polynomial.
******************************************************************************/
static uint32_t *get_poly(const uint64_t step) {
  uint32_t *poly = calloc(N * 11, sizeof(uint32_t));
  if (!poly) {
    return NULL;
  }

  uint32_t *pm = poly + N;      /* 2N words for the result of multiplication */
  uint32_t *tmp = pm + (N << 1);        /* temporary array for multiplication */

  /* split n with base 8 */
  int i, init;
  uint64_t n = step;
  i = init = 0;
  while (n) {
    int j = n & 0x07UL;
    if (j) {
      if (!init) {      /* initialise the polynomial */
        memcpy(poly, mt19937_poly[i][j-1], sizeof(uint32_t) * N);
        init = 1;
      }
      else {            /* polynomial multiplication and modular reduction */
        poly_mul(pm, poly, mt19937_poly[i][j-1], N, tmp);
        poly_mod_phi(pm, tmp);
        memcpy(poly, pm, sizeof(uint32_t) * N);
      }
    }
    i += 1;
    n >>= 3;
  }

  if (!init) memcpy(poly, mt19937_poly[0][0], sizeof(uint32_t) * N);
  return poly;
}

/******************************************************************************
Function `state_forward`:
  Update the state with polynomial multiplication.
Arguments:
  * `out`:      the pointer to the resulting state;
  * `in`:       the pointer to the initial state;
  * `poly`:     the pre-computed power of minimal polynomial.
******************************************************************************/
static void state_forward(mt19937_state_t *out, const mt19937_state_t *in,
    uint32_t *poly) {
  uint32_t *pt, *pm, *ph, *tmp;

  pt = poly;            /* N words for storing the pre-computed polynomial */
  pm = poly + N;        /* 2N words for the polynomial from advancing states */
  ph = pm + (N << 1);   /* 3N words for the multiplication */
  tmp = ph + N * 3;     /* 5N words for temporary array */

  /* Compute the polynomial pm by advancing the generator for 2K steps. */
  copy_state(out, in);
  memset(pm, 0, sizeof(uint32_t) * N * 2);
  for (int i = 2 * K - 1; i >= 0; i--) {
    pm[DIV_NBIT(i)] |= (next_state(out) & 1UL) << MOD_NBIT(i);
  }

  /* Compute ph = pt * pm, and extract coefficients from 2K-1 to K. */
  poly_mul_ub(ph, pm, pt, N, tmp);
  memset(pm, 0, sizeof(uint32_t) * N);
  for (int i = 0; i <= K; i++) {
    int j = 2 * K - 1 - i;
    pm[DIV_NBIT(i)] |= COEF(ph, j) << MOD_NBIT(i);
  }

  /* Recover state from the polynomial. */
  recover_state(out, pm);
}

/******************************************************************************
Function `mt19937_jump_seq`:
  Jump ahead in sequence for each stream.
  It is supposed to be called only at initialisation.
Arguments:
  * `state`:            the state array for multiple streams;
  * `init_state`:       the pointer to the initial state;
  * `nstream`:          total number of streams;
  * `step`:             step size for jumping ahead;
  * `err`:              an integer for storing the error message.
******************************************************************************/
static void mt19937_jump_seq(void **state, const void *init_state,
    const unsigned int nstream, const uint64_t step, int *err) {
  mt19937_state_t **stat = (mt19937_state_t **) state;
  mt19937_state_t *istat = (mt19937_state_t *) init_state;
  uint32_t *poly;

  if (PRAND_IS_ERROR(*err)) return;
  /* fill state[0] with init_state */
  copy_state(stat[0], istat);

  /* check of the upper limit should be done before calling this function */
  if (!step) {
    for (unsigned int i = 1; i < nstream; i++) copy_state(stat[i], stat[i - 1]);
    return;
  }

  /* jump-ahead polynomial */
  poly = get_poly(step);
  if (!poly) {
    *err = PRAND_ERR_MEMORY_JUMP;
    return;
  }

  /* Advance states with the polynomial. */
  for (unsigned int i = 1; i < nstream; i++)
    state_forward(stat[i], stat[i - 1], poly);

  free(poly);
}

/******************************************************************************
Function `mt19937_jump`:
  Jump ahead for one stream.
Arguments:
  * `state`:    the current state (to be over-written);
  * `step`:     step size for jumping ahead;
  * `err`:      an integer for storing the error message.
******************************************************************************/
static void mt19937_jump(void *state, const uint64_t step, int *err) {
  mt19937_state_t *stat = (mt19937_state_t *) state;
  if (PRAND_IS_ERROR(*err)) return;

  if (!step) return;
  else if (step > MT19937_MAX_STEP) {
    *err = PRAND_ERR_STEP;
    return;
  }

  /* jump-ahead polynomial */
  uint32_t *poly = get_poly(step);
  if (!poly) {
    *err = PRAND_ERR_MEMORY_JUMP;
    return;
  }

  /* Advance states with the polynomial. */
  state_forward(stat, stat, poly);

  free(poly);
}

/******************************************************************************
Function `mt19937_jump_all`:
  Jump ahead the same number of steps for all streams.
Arguments:
  * `rng`:      the random number generator interface;
  * `step`:     step size for jumping ahead;
  * `err`:      an integer for storing the error message.
******************************************************************************/
static void mt19937_jump_all(prand_t *rng, const uint64_t step, int *err) {
  if (PRAND_IS_ERROR(*err)) return;

  if (!step) return;
  else if (step > MT19937_MAX_STEP) {
    *err = PRAND_ERR_STEP;
    return;
  }

  /* jump-ahead polynomial */
  uint32_t *poly = get_poly(step);
  if (!poly) {
    *err = PRAND_ERR_MEMORY_JUMP;
    return;
  }

  /* Advance states with the polynomial. */
  for (int i = 0; i < rng->nstream; i++)
    state_forward(((mt19937_state_t **) (rng->state_stream))[i],
        ((mt19937_state_t **) (rng->state_stream))[i], poly);

  free(poly);
}

/******************************************************************************
Function `mt19937_reset`:
  Reset the state for one stream, with a given seed and number of skip steps.
Arguments:
  * `state`:    the current state (to be over-written);
  * `seed`:     an integer for initalisation the generator;
  * `step`:     step size for jumping ahead;
  * `err`:      an integer for storing the error message.
******************************************************************************/
static void mt19937_reset(void *state, const uint64_t seed,
    const uint64_t step, int *err) {
  mt19937_state_t *stat = (mt19937_state_t *) state;
  if (PRAND_IS_ERROR(*err)) return;

  if (seed == 0) {
    *err = PRAND_WARN_SEED;
    mt19937_seed(stat, DEFAULT_SEED);
  }
  else mt19937_seed(stat, seed);

  mt19937_jump(stat, step, err);
}

/******************************************************************************
Function `mt19937_reset_all`:
  Reset the state for all streams, with a given seed and number of skip steps.
Arguments:
  * `rng`:      the random number generator interface;
  * `seed`:     an integer for initalisation the generator;
  * `step`:     step size for jumping ahead;
  * `err`:      an integer for storing the error message.
******************************************************************************/
static void mt19937_reset_all(prand_t *rng, const uint64_t seed,
    const uint64_t step, int *err) {
  mt19937_state_t *stat = (mt19937_state_t *) (rng->state);
  if (PRAND_IS_ERROR(*err)) return;

  if (seed == 0) {
    *err = PRAND_WARN_SEED;
    mt19937_seed(stat, DEFAULT_SEED);
  }
  else mt19937_seed(stat, seed);

  if (!step) {
    for (int i = 1; i < rng->nstream; i++)
      memcpy(rng->state_stream[i], stat, sizeof(mt19937_state_t));
    return;
  }
  else if (step > MT19937_MAX_STEP) {
    *err = PRAND_ERR_STEP;
    return;
  }

  if (rng->nstream <= 1) mt19937_jump(stat, step, err);
  else mt19937_jump_seq(rng->state_stream, stat, rng->nstream, step, err);
}


/*============================================================================*\
                          Interface for initialisation
\*============================================================================*/

/******************************************************************************
Function `mt19937_init`:
  Initialisation of the MT19937 generator, with the universal API.
Arguments:
  * `seed`:     an integer for initalisation the generator;
  * `nstream`:  total number of streams;
  * `step`:     step size for jumping ahead;
  * `err`:      an integer for storing the error message.
Return:
  A universal instance of the random number generator.
******************************************************************************/
prand_t *mt19937_init(const uint64_t seed, const unsigned int nstream,
    const uint64_t step, int *err) {
  /* `step` should not be larger than the pre-computed length. */
  if (step > MT19937_MAX_STEP) {
    *err = PRAND_ERR_STEP;
    return NULL;
  }

  prand_t *rng = malloc(sizeof(prand_t));
  if (!rng) {
    *err = PRAND_ERR_MEMORY;
    return NULL;
  }

  unsigned int numstr = (nstream == 0) ? 1 : nstream;
  rng->state_stream = malloc(sizeof(mt19937_state_t *) * numstr);
  if (!rng->state_stream) {
    free(rng);
    *err = PRAND_ERR_MEMORY;
    return NULL;
  }

  mt19937_state_t *states =  malloc(sizeof(mt19937_state_t) * numstr);
  if (!states) {
    free(rng->state_stream);
    free(rng);
    *err = PRAND_ERR_MEMORY;
    return NULL;
  }
  for (unsigned int i = 0; i < numstr; i++)
    rng->state_stream[i] = states + i;

  rng->state = rng->state_stream[0];
  rng->nstream = numstr;
  rng->type = PRAND_RNG_MT19937;
  rng->min = 0;
  rng->max = 0xffffffffUL;      /* 2^32 - 1 */

  rng->get = &mt19937_get;
  rng->get_double = &mt19937_get_double;
  rng->get_double_pos = &mt19937_get_double_pos;
  rng->reset = &mt19937_reset;
  rng->reset_all = &mt19937_reset_all;
  rng->jump = &mt19937_jump;
  rng->jump_all = &mt19937_jump_all;

  if (seed == 0) {
    *err = PRAND_WARN_SEED;
    mt19937_seed(rng->state, DEFAULT_SEED);
  }
  else mt19937_seed(rng->state, seed);

  if (nstream <= 1)
    mt19937_jump(rng->state, step, err);
  else
    mt19937_jump_seq(rng->state_stream, rng->state, nstream, step, err);

  return rng;
}

