/*******************************************************************************
* mrg32k3a.c: this file is part of the prand library.
 
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

#include "mrg32k3a.h"
#include "mrg32k3a_jump.h"
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

/*******************************************************************************
  Implementation of the MRG32k3a random number generator.
  ref: https://doi.org/10.1287/opre.47.1.159
*******************************************************************************/

/*============================================================================*\
                            Definitions of constants
\*============================================================================*/

#define m1              4294967087LL            /* 2^32 - 209 */
#define m2              4294944443LL            /* 2^32 - 22853 */
#define a12             1403580
#define a13             (-810728)
#define a21             527612
#define a23             (-1370589)

/* The following two numbers ensure postive p1 and p2. */
#define add1            3482050076509336LL      /* m1 * a13 */
#define add2            5886603609186927LL      /* m2 * a23 */

/* Normalisation for sampling a float-point number in the range [0,1). */
#define norm            0x1.000000d00000bp-32   /* 1 / (1 + m1) */
/* Normalisation for sampling a float-point number in the range (0,1). */
#define norm_pos        0x1.000000cf0000ap-32   /* 1 / (2 + m1) */


/*============================================================================*\
                         Macros for the initialisation
\*============================================================================*/

/* LCG for initialising the states. */
#define PRAND_LCG(n) ((69069 * n + 1) & 0xffffffffUL)

#define DEFAULT_SEED    1


/*============================================================================*\
                            Definition of the state
\*============================================================================*/

typedef struct {
  int64_t s10, s11, s12;
  int64_t s20, s21, s22;
} mrg32k3a_state_t;


/*============================================================================*\
                     Functions for random number generation
\*============================================================================*/

/******************************************************************************
Function `mrg32k3a_seed`:
  Initialise the state with an integer.
Arguments:
  * `state`:    the state to be intialised;
  * `seed`:     a positive integer for the initialisation.
******************************************************************************/
static void mrg32k3a_seed(void *state, uint64_t seed) {
  mrg32k3a_state_t *stat = (mrg32k3a_state_t *) state;

  /* Initialise states with LCG.
   * The validation of the seed is done in `mrg32k3a_init`. */
  seed = PRAND_LCG(seed);
  stat->s10 = seed % m1;
  seed = PRAND_LCG(seed);
  stat->s11 = seed % m1;
  seed = PRAND_LCG(seed);
  stat->s12 = seed % m1;

  seed = PRAND_LCG(seed);
  stat->s20 = seed % m2;
  seed = PRAND_LCG(seed);
  stat->s21 = seed % m2;
  seed = PRAND_LCG(seed);
  stat->s22 = seed % m2;
}

/******************************************************************************
Function `mrg32k3a_get`:
  Generate an integer and update the state.
Arguments:
  * `state`:    the state for the generator.
Return:
  A pseudo-random integer.
******************************************************************************/
static uint64_t mrg32k3a_get(void *state) {
  mrg32k3a_state_t *stat = (mrg32k3a_state_t *) state;

  /* Component 1 */
  int64_t p1 = (a12 * stat->s11 + a13 * stat->s10 + add1) % m1;
  stat->s10 = stat->s11;
  stat->s11 = stat->s12;
  stat->s12 = p1;

  /* Component 2 */
  int64_t p2 = (a21 * stat->s22 + a23 * stat->s20 + add2) % m2;
  stat->s20 = stat->s21;
  stat->s21 = stat->s22;
  stat->s22 = p2;

  /* Combination */
  if (p1 <= p2) return (p1 - p2 + m1);
  else return (p1 - p2);
}

/******************************************************************************
Function `mrg32k3a_get_double`:
  Generate a double-precision floating-point number in the range [0,1).
Arguments:
  * `state`:    the state for the generator.
Return:
  A pseudo-random floating-point number.
******************************************************************************/
static double mrg32k3a_get_double(void *state) {
  return mrg32k3a_get(state) * norm;
}

/******************************************************************************
Function `mrg32k3a_get_double_pos`:
  Generate a double-precision floating-point number in the range (0,1).
Arguments:
  * `state`:    the state for the generator.
Return:
  A pseudo-random floating-point number.
******************************************************************************/
static double mrg32k3a_get_double_pos(void *state) {
  return (mrg32k3a_get(state) + 1) * norm_pos;
}


/*============================================================================*\
                         Functions for multiple streams
\*============================================================================*/

/******************************************************************************
  ref: https://doi.org/10.1016/B978-0-12-384988-5.00016-4
******************************************************************************/

/******************************************************************************
Function `matrix_dot`:
  Inplace 3x3 matrix manipulation: C = (A * B) % m.
Arguments:
  * `C`:        the output matrix, can be the same as either of the multipliers;
  * `A`, `B`:   the multipliers;
  * `m`:        the divisor for modular reduction.
******************************************************************************/
static inline void matrix_dot(uint64_t C[9], const uint64_t A[9],
    const uint64_t B[9], const int64_t m) {
  uint64_t out[9];

  /* Perform modular reduction at all steps to prevent overflow. */
  out[0] = ((A[0] * B[0]) % m + (A[1] * B[3]) % m + (A[2] * B[6]) % m) % m;
  out[1] = ((A[0] * B[1]) % m + (A[1] * B[4]) % m + (A[2] * B[7]) % m) % m;
  out[2] = ((A[0] * B[2]) % m + (A[1] * B[5]) % m + (A[2] * B[8]) % m) % m;
  out[3] = ((A[3] * B[0]) % m + (A[4] * B[3]) % m + (A[5] * B[6]) % m) % m;
  out[4] = ((A[3] * B[1]) % m + (A[4] * B[4]) % m + (A[5] * B[7]) % m) % m;
  out[5] = ((A[3] * B[2]) % m + (A[4] * B[5]) % m + (A[5] * B[8]) % m) % m;
  out[6] = ((A[6] * B[0]) % m + (A[7] * B[3]) % m + (A[8] * B[6]) % m) % m;
  out[7] = ((A[6] * B[1]) % m + (A[7] * B[4]) % m + (A[8] * B[7]) % m) % m;
  out[8] = ((A[6] * B[2]) % m + (A[7] * B[5]) % m + (A[8] * B[8]) % m) % m;

  C[0] = out[0];
  C[1] = out[1];
  C[2] = out[2];
  C[3] = out[3];
  C[4] = out[4];
  C[5] = out[5];
  C[6] = out[6];
  C[7] = out[7];
  C[8] = out[8];
}

/******************************************************************************
Function `matrix_pow`:
  Compute the power of the two state advancing matrices for MRG32k3a,
  with modular reduction:
    A1 = A1^(step) % m1, A2 = A2^(step) % m2
Arguments:
  * `A1`, `A2`: the two state advancing matrices for the MRG32k3a generator;
  * `step`:     the number of steps to be skipped.
******************************************************************************/
static void matrix_pow(uint64_t A1[9], uint64_t A2[9], const uint64_t step) {
  uint64_t n = step;

  int i = 0;
  bool init = false;
  while (n) {
    int j = n & 0x07UL;
    if (j) {
      if (!init) {      /* Initialise A1 and A2. */
        for (int k = 0; k < 9; k++) {
          A1[k] = A1gj8i[i][j-1][k];
          A2[k] = A2gj8i[i][j-1][k];
          init = 1;
        }
      }
      else {            /* Compute A1 and A2 by matrix production. */
        matrix_dot(A1, A1gj8i[i][j-1], A1, m1);
        matrix_dot(A2, A2gj8i[i][j-1], A2, m2);
      }
    }
    i += 1;
    n >>= 3;
  }

  if (!init) {          /* step == 0 */
    for (int k = 0; k < 9; k++) {
      A1[k] = A1gj8i[0][0][k];
      A2[k] = A2gj8i[0][0][k];
    }
  }
}

/******************************************************************************
Function `copy_state`:
  Copy the state.
Arguments:
  * `dst`:      the destination;
  * `src`:      the state to be copied.
******************************************************************************/
static inline void copy_state(mrg32k3a_state_t *dst,
    const mrg32k3a_state_t *src) {
  if (dst == src) return;
  dst->s10 = src->s10;
  dst->s11 = src->s11;
  dst->s12 = src->s12;
  dst->s20 = src->s20;
  dst->s21 = src->s21;
  dst->s22 = src->s22;
}

/******************************************************************************
Function `state_forward`:
  Update the state with matrix production.
Arguments:
  * `out`:      the pointer to the resulting state;
  * `in`:       the pointer to the initial state;
  * `A1`, `A2`: the conversion matrices.
******************************************************************************/
static inline void state_forward(mrg32k3a_state_t *out,
    const mrg32k3a_state_t *in, const uint64_t A1[9], const uint64_t A2[9]) {
  uint64_t s0, s1, s2;
  s0 = in->s10;
  s1 = in->s11;
  s2 = in->s12;
  out->s10 = ((A1[0]*s0) % m1 + (A1[1]*s1) % m1 + (A1[2]*s2) % m1) % m1;
  out->s11 = ((A1[3]*s0) % m1 + (A1[4]*s1) % m1 + (A1[5]*s2) % m1) % m1;
  out->s12 = ((A1[6]*s0) % m1 + (A1[7]*s1) % m1 + (A1[8]*s2) % m1) % m1;

  s0 = in->s20;
  s1 = in->s21;
  s2 = in->s22;
  out->s20 = ((A2[0]*s0) % m2 + (A2[1]*s1) % m2 + (A2[2]*s2) % m2) % m2;
  out->s21 = ((A2[3]*s0) % m2 + (A2[4]*s1) % m2 + (A2[5]*s2) % m2) % m2;
  out->s22 = ((A2[6]*s0) % m2 + (A2[7]*s1) % m2 + (A2[8]*s2) % m2) % m2;
}

/******************************************************************************
Function `mrg32k3a_jump_seq`:
  Jump ahead in sequence for each stream.
  It is supposed to be called only at initialisation.
Arguments:
  * `state`:            the state array for multiple streams;
  * `init_state`:       the pointer to the initial state;
  * `nstream`:          total number of streams;
  * `step`:             step size for jumping ahead;
  * `err`:              an integer for storing the error message.
******************************************************************************/
static void mrg32k3a_jump_seq(void **state, const void *init_state,
    const unsigned int nstream, const uint64_t step, int *err) {
  mrg32k3a_state_t **stat = (mrg32k3a_state_t **) state;
  mrg32k3a_state_t *istat = (mrg32k3a_state_t *) init_state;
  uint64_t A1[9], A2[9];

  if (PRAND_IS_ERROR(*err)) return;
  /* fill state[0] with init_state */
  copy_state(stat[0], istat);

  /* check of the upper limit should be done before calling this function */
  if (!step) {
    for (unsigned int i = 1; i < nstream; i++) copy_state(stat[i], stat[i - 1]);
    return;
  }

  /* jump-ahead matrices */
  matrix_pow(A1, A2, step);

  /* Advance states with the matrices. */
  for (unsigned int i = 1; i < nstream; i++)
    state_forward(stat[i], stat[i - 1], A1, A2);
}

/******************************************************************************
Function `mrg32k3a_jump`:
  Jump ahead for one stream.
Arguments:
  * `state`:    the current state (to be over-written);
  * `step`:     step size for jumping ahead;
  * `err`:      an integer for storing the error message.
******************************************************************************/
static void mrg32k3a_jump(void *state, const uint64_t step, int *err) {
  mrg32k3a_state_t *stat = (mrg32k3a_state_t *) state;
  uint64_t A1[9], A2[9];
  if (PRAND_IS_ERROR(*err)) return;

  if (!step) return;
  else if (step > MRG32K3A_MAX_STEP) {
    *err = PRAND_ERR_STEP;
    return;
  }

  /* jump-ahead matrices */
  matrix_pow(A1, A2, step);

  /* Advance states with the matrices. */
  state_forward(stat, stat, A1, A2);
}

/******************************************************************************
Function `mrg32k3a_jump_all`:
  Jump ahead the same number of steps for all streams.
Arguments:
  * `rng`:      the random number generator interface;
  * `step`:     step size for jumping ahead;
  * `err`:      an integer for storing the error message.
******************************************************************************/
static void mrg32k3a_jump_all(prand_t *rng, const uint64_t step, int *err) {
  uint64_t A1[9], A2[9];
  if (PRAND_IS_ERROR(*err)) return;

  if (!step) return;
  else if (step > MRG32K3A_MAX_STEP) {
    *err = PRAND_ERR_STEP;
    return;
  }

  /* jump-ahead matrices */
  matrix_pow(A1, A2, step);

  /* Advance states with the matrices. */
  for (int i = 0; i < rng->nstream; i++)
    state_forward(((mrg32k3a_state_t **) (rng->state_stream))[i],
        ((mrg32k3a_state_t **) (rng->state_stream))[i], A1, A2);
}

/******************************************************************************
Function `mrg32k3a_reset`:
  Reset the state for one stream, with a given seed and number of skip steps.
Arguments:
  * `state`:    the current state (to be over-written);
  * `seed`:     an integer for initalisation the generator;
  * `step`:     step size for jumping ahead;
  * `err`:      an integer for storing the error message.
******************************************************************************/
static void mrg32k3a_reset(void *state, const uint64_t seed,
    const uint64_t step, int *err) {
  mrg32k3a_state_t *stat = (mrg32k3a_state_t *) state;
  if (PRAND_IS_ERROR(*err)) return;

  if (seed == 0) {
    *err = PRAND_WARN_SEED;
    mrg32k3a_seed(stat, DEFAULT_SEED);
  }
  else mrg32k3a_seed(stat, seed);

  mrg32k3a_jump(stat, step, err);
}

/******************************************************************************
Function `mrg32k3a_reset_all`:
  Reset the state for all streams, with a given seed and number of skip steps.
Arguments:
  * `rng`:      the random number generator interface;
  * `seed`:     an integer for initalisation the generator;
  * `step`:     step size for jumping ahead;
  * `err`:      an integer for storing the error message.
******************************************************************************/
static void mrg32k3a_reset_all(prand_t *rng, const uint64_t seed,
    const uint64_t step, int *err) {
  mrg32k3a_state_t *stat = (mrg32k3a_state_t *) (rng->state);
  if (PRAND_IS_ERROR(*err)) return;

  if (seed == 0) {
    *err = PRAND_WARN_SEED;
    mrg32k3a_seed(stat, DEFAULT_SEED);
  }
  else mrg32k3a_seed(stat, seed);

  if (!step) {
    for (int i = 1; i < rng->nstream; i++)
      memcpy(rng->state_stream[i], stat, sizeof(mrg32k3a_state_t));
    return;
  }
  else if (step > MRG32K3A_MAX_STEP) {
    *err = PRAND_ERR_STEP;
    return;
  }

  if (rng->nstream <= 1) mrg32k3a_jump(stat, step, err);
  else mrg32k3a_jump_seq(rng->state_stream, stat, rng->nstream, step, err);
}


/*============================================================================*\
                          Interface for initialisation
\*============================================================================*/

/******************************************************************************
Function `mrg32k3a_init`:
  Initialisation of the MRG32k3a generator, with the universal API.
Arguments:
  * `seed`:     an integer for initalisation the generator;
  * `nstream`:  total number of streams;
  * `step`:     step size for jumping ahead;
  * `err`:      an integer for storing the error message.
Return:
  A universal instance of the random number generator.
******************************************************************************/
prand_t *mrg32k3a_init(const uint64_t seed, const unsigned int nstream,
    const uint64_t step, int *err) {
  /* `step` should not be larger than the pre-computed length. */
  if (step > MRG32K3A_MAX_STEP) {
    *err = PRAND_ERR_STEP;
    return NULL;
  }

  prand_t *rng = malloc(sizeof(prand_t));
  if (!rng) {
    *err = PRAND_ERR_MEMORY;
    return NULL;
  }

  unsigned int numstr = (nstream == 0) ? 1 : nstream;
  rng->state_stream = malloc(sizeof(mrg32k3a_state_t *) * numstr);
  if (!rng->state_stream) {
    free(rng);
    *err = PRAND_ERR_MEMORY;
    return NULL;
  }

  mrg32k3a_state_t *states = malloc(sizeof(mrg32k3a_state_t) * numstr);
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
  rng->type = PRAND_RNG_MRG32K3A;
  rng->min = 0;
  rng->max = m1;

  rng->get = &mrg32k3a_get;
  rng->get_double = &mrg32k3a_get_double;
  rng->get_double_pos = &mrg32k3a_get_double_pos;
  rng->reset = &mrg32k3a_reset;
  rng->reset_all = &mrg32k3a_reset_all;
  rng->jump = &mrg32k3a_jump;
  rng->jump_all = &mrg32k3a_jump_all;

  if (seed == 0) {
    *err = PRAND_WARN_SEED;
    mrg32k3a_seed(rng->state, DEFAULT_SEED);
  }
  else mrg32k3a_seed(rng->state, seed);

  if (nstream == 0)
    mrg32k3a_jump(rng->state, step, err);
  else
    mrg32k3a_jump_seq(rng->state_stream, rng->state, nstream, step, err);

  return rng;
}

