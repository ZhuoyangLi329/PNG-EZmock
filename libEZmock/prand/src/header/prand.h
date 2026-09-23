/*******************************************************************************
* prand.h: this file is part of the prand library.
 
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

#ifndef __PRAND_H__
#define __PRAND_H__

#include <stdint.h>

/*============================================================================*\
                   List of available random number generators
\*============================================================================*/
typedef enum {
  PRAND_RNG_MRG32K3A = 0,
  PRAND_RNG_MT19937 = 1
} prand_rng_enum;


/*============================================================================*\
                           Definitions of error codes
\*============================================================================*/
#define PRAND_ERR_MEMORY                (-1)
#define PRAND_ERR_MEMORY_JUMP           (-2)
#define PRAND_ERR_STEP                  (-3)
#define PRAND_ERR_UNDEF_RNG             (-4)
#define PRAND_WARN_SEED                 1

#define PRAND_IS_ERROR(err)             ((err) < 0)
#define PRAND_IS_WARN(err)              ((err) > 0)

/*============================================================================*\
              Universal interface of the random number generators
\*============================================================================*/

typedef struct prand_struct {
  void *state;                  /* the state for single stream */
  void **state_stream;          /* states for multiple streams */
  int nstream;                  /* number of random streams */
  prand_rng_enum type;          /* type of the random number generator */
  int64_t min;                  /* minimum value of the random integer */
  int64_t max;                  /* maximum value of the random integer */
  /* function pointers for sampling numbers */
  uint64_t (*get) (void *);
  double (*get_double) (void *);
  double (*get_double_pos) (void *);
  /* function pointers for reseting states with seed and skipping steps */
  void (*reset) (void *, const uint64_t, const uint64_t, int *);
  void (*reset_all) (struct prand_struct *, const uint64_t, const uint64_t,
      int *);
  /* function pointers for jumping ahead */
  void (*jump) (void *, const uint64_t, int *);
  void (*jump_all) (struct prand_struct *, const uint64_t, int *);
} prand_t;

/******************************************************************************
Function `prand_init`:
  Initialisation of the interface for the selected random number generator.
Arguments:
  * `type`:     the ID of the pre-defined random number generator;
  * `seed`:     an integer for initalisation the generator;
  * `nstream`:  total number of streams;
  * `step`:     step size for jumping ahead;
  * `err`:      an integer for storing the error message.
Return:
  A universal interface of the random number generator.
******************************************************************************/
prand_t *prand_init(const prand_rng_enum type, const uint64_t seed,
    const unsigned int nstream, const uint64_t step, int *err);

/******************************************************************************
Function `prand_errmsg`:
  Produce error message for a given error code.
Arguments:
  * `err`:      the error code
Return:
  A string with error message.
******************************************************************************/
const char *prand_errmsg(const int err);

/******************************************************************************
Function `prand_destroy`:
  Release memory allocated for the random number generator interface.
Arguments:
  * `rng`:      the instance of the random number generator.
******************************************************************************/
void prand_destroy(prand_t *rng);

#endif

