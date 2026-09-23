/*******************************************************************************
* mt19937.h: this file is part of the prand library.
 
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

#ifndef __MT19937_H__
#define __MT19937_H__

#include "prand.h"

/*============================================================================*\
      Definitions for the representation of polynomials with 32-bit words
\*============================================================================*/

#define MT19937_N       624
#define WORD_SIZE       32
#define MUL_NBIT(x)     ((x) << 5)
#define DIV_NBIT(x)     ((x) >> 5)
#define MOD_NBIT(x)     ((x) & 0x1fU)
#define COEF(x,i)       (((x[DIV_NBIT(i)]) >> MOD_NBIT(i)) & 1)


/*============================================================================*\
    Low-level functions for polynomial multiplication and modular reduction
\*============================================================================*/

/******************************************************************************
Function `poly_mul`:
  Compute r = a * b, where a and b contain both n words,
    using a pre-defined algorithm.
Arguments:
  * `r`:        pointer to the result, with at least 2 * `n` words;
  * `a`, `b`:   pointers to the two multipliers, each with `n` words;
  * `n`:        the length of `a` and `b` (in words);
  * `tmp`:      a temporary array with a rough requirement of 4 * `n` words.
TODO:
  Include the toom-3 algorithm for long polynomials.
******************************************************************************/
void poly_mul(uint32_t *r, const uint32_t *a, const uint32_t *b,
    const unsigned int n, uint32_t *tmp);

/******************************************************************************
Function `poly_mul_ub`:
  Compute r = a * b, where a and b contain 2n and n words, respectively.
Arguments:
  * `r`:        pointer to the result, with at least 3 * `n` words;
  * `a`, `b`:   pointers to the two multipliers;
  * `n`:        the length of `b` (in words);
  * `tmp`:      a temporary array with a rough requirement of 5 * `n` words.
******************************************************************************/
void poly_mul_ub(uint32_t *r, const uint32_t *a, const uint32_t *b,
    const unsigned int n, uint32_t *tmp);

/******************************************************************************
Function `poly_mod_phi`:
  Compute r %= phi, where phi is the minimal polynomial for MT19937.
Arguments:
  * `r`:        pointer to the result;
  * `tmp`:      a temporary array, with at least the size of one pre-defined
                block.
******************************************************************************/
void poly_mod_phi(uint32_t *r, uint32_t *tmp);


/*============================================================================*\
                            Initialisation function
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
    const uint64_t step, int *err);

#endif

