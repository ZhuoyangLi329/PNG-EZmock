/*******************************************************************************
* mt19937_poly.c: this file is part of the prand library.
 
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
#include <string.h>

/*============================================================================*\
           Functions for polynomial multiplication with 32-bit words
\*============================================================================*/

/* Maximum number of words for the expanded multiplication functions. */
#define EXPD_MUL_THRES  6

/* Lookup table for the grade-school multiplication algorithm. */
static const uint32_t WORD_MASK[2] = {UINT32_C(0), UINT32_C(0xffffffff)};

/******************************************************************************
Function `poly_mul1`:
  Compute r = a * b, using the grade-school multiplication algorithm.
Arguments:
  * `r`:        pointer to the result, with at least 2 words;
  * `a`, `b`:   the two multipliers.
TODO:
  Process more bits at once with a larger lookup table.
******************************************************************************/
static inline void poly_mul1(uint32_t *r, const uint32_t a, const uint32_t b) {
  uint32_t t;
  r[0] = r[1] = 0;
  t = b & WORD_MASK[a & 1]; r[0] ^= t;
  t = b & WORD_MASK[(a >> 1) & 1]; r[0] ^= t << 1; r[1] ^= t >> 31;
  t = b & WORD_MASK[(a >> 2) & 1]; r[0] ^= t << 2; r[1] ^= t >> 30;
  t = b & WORD_MASK[(a >> 3) & 1]; r[0] ^= t << 3; r[1] ^= t >> 29;
  t = b & WORD_MASK[(a >> 4) & 1]; r[0] ^= t << 4; r[1] ^= t >> 28;
  t = b & WORD_MASK[(a >> 5) & 1]; r[0] ^= t << 5; r[1] ^= t >> 27;
  t = b & WORD_MASK[(a >> 6) & 1]; r[0] ^= t << 6; r[1] ^= t >> 26;
  t = b & WORD_MASK[(a >> 7) & 1]; r[0] ^= t << 7; r[1] ^= t >> 25;
  t = b & WORD_MASK[(a >> 8) & 1]; r[0] ^= t << 8; r[1] ^= t >> 24;
  t = b & WORD_MASK[(a >> 9) & 1]; r[0] ^= t << 9; r[1] ^= t >> 23;
  t = b & WORD_MASK[(a >> 10) & 1]; r[0] ^= t << 10; r[1] ^= t >> 22;
  t = b & WORD_MASK[(a >> 11) & 1]; r[0] ^= t << 11; r[1] ^= t >> 21;
  t = b & WORD_MASK[(a >> 12) & 1]; r[0] ^= t << 12; r[1] ^= t >> 20;
  t = b & WORD_MASK[(a >> 13) & 1]; r[0] ^= t << 13; r[1] ^= t >> 19;
  t = b & WORD_MASK[(a >> 14) & 1]; r[0] ^= t << 14; r[1] ^= t >> 18;
  t = b & WORD_MASK[(a >> 15) & 1]; r[0] ^= t << 15; r[1] ^= t >> 17;
  t = b & WORD_MASK[(a >> 16) & 1]; r[0] ^= t << 16; r[1] ^= t >> 16;
  t = b & WORD_MASK[(a >> 17) & 1]; r[0] ^= t << 17; r[1] ^= t >> 15;
  t = b & WORD_MASK[(a >> 18) & 1]; r[0] ^= t << 18; r[1] ^= t >> 14;
  t = b & WORD_MASK[(a >> 19) & 1]; r[0] ^= t << 19; r[1] ^= t >> 13;
  t = b & WORD_MASK[(a >> 20) & 1]; r[0] ^= t << 20; r[1] ^= t >> 12;
  t = b & WORD_MASK[(a >> 21) & 1]; r[0] ^= t << 21; r[1] ^= t >> 11;
  t = b & WORD_MASK[(a >> 22) & 1]; r[0] ^= t << 22; r[1] ^= t >> 10;
  t = b & WORD_MASK[(a >> 23) & 1]; r[0] ^= t << 23; r[1] ^= t >> 9;
  t = b & WORD_MASK[(a >> 24) & 1]; r[0] ^= t << 24; r[1] ^= t >> 8;
  t = b & WORD_MASK[(a >> 25) & 1]; r[0] ^= t << 25; r[1] ^= t >> 7;
  t = b & WORD_MASK[(a >> 26) & 1]; r[0] ^= t << 26; r[1] ^= t >> 6;
  t = b & WORD_MASK[(a >> 27) & 1]; r[0] ^= t << 27; r[1] ^= t >> 5;
  t = b & WORD_MASK[(a >> 28) & 1]; r[0] ^= t << 28; r[1] ^= t >> 4;
  t = b & WORD_MASK[(a >> 29) & 1]; r[0] ^= t << 29; r[1] ^= t >> 3;
  t = b & WORD_MASK[(a >> 30) & 1]; r[0] ^= t << 30; r[1] ^= t >> 2;
  t = b & WORD_MASK[(a >> 31) & 1]; r[0] ^= t << 31; r[1] ^= t >> 1;
}

/******************************************************************************
Function `kara_mulX` (X = 2,3,4,5,6):
  Compute r = a * b, where a and b both contain X words,
    using the Karatsuba algorithm.
  In particular, for X = 3, we use a three-way splitting method.
Arguments:
  * `r`:        pointer to the result, with at least 2X words;
  * `a`, `b`:   pointers to the two multipliers, each with X words;
  * `tmp`:      a temporary array with a rough requirement of 4X words;
******************************************************************************/
static inline void kara_mul2(uint32_t *r, const uint32_t *a, const uint32_t *b,
    uint32_t *tmp) {
  poly_mul1(r, a[0], b[0]);
  poly_mul1(r + 2, a[1], b[1]);
  /* The pre-computed value saves one operation, see below. */
  tmp[0] = r[1] ^ r[2];

  poly_mul1(tmp + 1, a[0] ^ a[1], b[0] ^ b[1]);
  r[1] = r[0] ^ tmp[1] ^ tmp[0];
  r[2] = r[3] ^ tmp[2] ^ tmp[0];

  /* Naive approach:
  poly_mul1(tmp, a[0], b[0]);
  poly_mul1(tmp + 2, a[1], b[1]);
  poly_mul1(tmp + 4, a[0] ^ a[1], b[0] ^ b[1]);
  r[0] = tmp[0];
  r[1] = tmp[0] ^ tmp[1] ^ tmp[2] ^ tmp[4];
  r[2] = tmp[1] ^ tmp[2] ^ tmp[3] ^ tmp[5];
  r[3] = tmp[3];
  */
}

static inline void kara_mul3(uint32_t *r, const uint32_t *a, const uint32_t *b,
    uint32_t *tmp) {
  poly_mul1(r, a[0], b[0]);
  poly_mul1(r + 2, a[1], b[1]);
  poly_mul1(r + 4, a[2], b[2]);
  poly_mul1(tmp + 4, a[0] ^ a[1], b[0] ^ b[1]);
  poly_mul1(tmp + 6, a[1] ^ a[2], b[1] ^ b[2]);
  poly_mul1(tmp + 8, a[2] ^ a[0], b[2] ^ b[0]);

  /* The pre-computed values save 6 operations, see below. */
  tmp[0] = r[1] ^ r[2];
  tmp[1] = r[0] ^ tmp[0];
  tmp[2] = r[3] ^ r[4];
  tmp[3] = tmp[2] ^ r[5];

  r[1] = tmp[1] ^ tmp[4];
  r[2] = tmp[1] ^ tmp[2] ^ tmp[5] ^ tmp[8];
  r[3] = tmp[0] ^ tmp[3] ^ tmp[6] ^ tmp[9];
  r[4] = tmp[3] ^ tmp[7];

  /* Naive approach:
  poly_mul1(tmp, a[0], b[0]);
  poly_mul1(tmp + 2, a[1], b[1]);
  poly_mul1(tmp + 4, a[2], b[2]);
  poly_mul1(tmp + 6, a[0] ^ a[1], b[0] ^ b[1]);
  poly_mul1(tmp + 8, a[1] ^ a[2], b[1] ^ b[2]);
  poly_mul1(tmp + 10, a[2] ^ a[0], b[2] ^ b[0]);
  r[0] = tmp[0];
  r[1] = tmp[0] ^ tmp[1] ^ tmp[2] ^ tmp[6];
  r[2] = tmp[0] ^ tmp[1] ^ tmp[2] ^ tmp[3] ^ tmp[4] ^ tmp[7] ^ tmp[10];
  r[3] = tmp[1] ^ tmp[2] ^ tmp[3] ^ tmp[4] ^ tmp[5] ^ tmp[8] ^ tmp[11];
  r[4] = tmp[3] ^ tmp[4] ^ tmp[5] ^ tmp[9];
  r[5] = tmp[5];
  */
}

static void kara_mul4(uint32_t *r, const uint32_t *a, const uint32_t *b,
    uint32_t *tmp) {
  kara_mul2(r, a, b, tmp + 6);
  kara_mul2(r + 4, a + 2, b + 2, tmp + 6);
  tmp[0] = r[2] ^ r[4];
  tmp[1] = r[3] ^ r[5];
  tmp[2] = a[0] ^ a[2];
  tmp[3] = a[1] ^ a[3];
  tmp[4] = b[0] ^ b[2];
  tmp[5] = b[1] ^ b[3];
  kara_mul2(r + 2, tmp + 2, tmp + 4, tmp + 6);
  r[2] ^= tmp[0] ^ r[0];
  r[3] ^= tmp[1] ^ r[1];
  r[4] ^= tmp[0] ^ r[6];
  r[5] ^= tmp[1] ^ r[7];
}

static void kara_mul5(uint32_t *r, const uint32_t *a, const uint32_t *b,
    uint32_t *tmp) {
  kara_mul3(r, a, b, tmp + 9);
  kara_mul2(r + 6, a + 3, b + 3, tmp + 9);
  tmp[0] = r[3] ^ r[6];
  tmp[1] = r[4] ^ r[7];
  tmp[2] = r[5] ^ r[8];
  tmp[3] = a[0] ^ a[3];
  tmp[4] = a[1] ^ a[4];
  tmp[5] = a[2];
  tmp[6] = b[0] ^ b[3];
  tmp[7] = b[1] ^ b[4];
  tmp[8] = b[2];
  kara_mul3(r + 3, tmp + 3, tmp + 6, tmp + 9);
  r[3] ^= tmp[0] ^ r[0];
  r[4] ^= tmp[1] ^ r[1];
  r[5] ^= tmp[2] ^ r[2];
  r[6] ^= tmp[0] ^ r[9];
  r[7] ^= tmp[1];
  r[8] ^= tmp[2];
}

static void kara_mul6(uint32_t *r, const uint32_t *a, const uint32_t *b,
    uint32_t *tmp) {
  kara_mul3(r, a, b, tmp + 9);
  kara_mul3(r + 6, a + 3, b + 3, tmp + 9);
  tmp[0] = r[3] ^ r[6];
  tmp[1] = r[4] ^ r[7];
  tmp[2] = r[5] ^ r[8];
  tmp[3] = a[0] ^ a[3];
  tmp[4] = a[1] ^ a[4];
  tmp[5] = a[2] ^ a[5];
  tmp[6] = b[0] ^ b[3];
  tmp[7] = b[1] ^ b[4];
  tmp[8] = b[2] ^ b[5];
  kara_mul3(r + 3, tmp + 3, tmp + 6, tmp + 9);
  r[3] ^= tmp[0] ^ r[0];
  r[4] ^= tmp[1] ^ r[1];
  r[5] ^= tmp[2] ^ r[2];
  r[6] ^= tmp[0] ^ r[9];
  r[7] ^= tmp[1] ^ r[10];
  r[8] ^= tmp[2] ^ r[11];
}

/******************************************************************************
Function `kara_mul`:
  Compute r = a * b, where a and b contain both n words,
    using the Karatsuba algorithm.
Arguments:
  * `r`:        pointer to the result, with at least 2 * `n` words;
  * `a`, `b`:   pointers to the two multipliers, each with `n` words;
  * `n`:        the length of `a` and `b` (in words);
  * `tmp`:      a temporary array with a rough requirement of 4 * `n` words.
******************************************************************************/
static void kara_mul(uint32_t *r, const uint32_t *a, const uint32_t *b,
    const unsigned int n, uint32_t *tmp) {
  unsigned int i;
  const unsigned int n1 = (n + 1) >> 1;         /* floor(n / 2) */
  const unsigned int n2 = n >> 1;               /* ceil(n / 2) */
  const uint32_t *a2 = a + n1;                  /* higher part of a */
  const uint32_t *b2 = b + n1;                  /* higher part of b */
  uint32_t *r1 = r + n1;
  uint32_t *r2 = r1 + n1;
  uint32_t *r3 = r2 + n1;

  /* Temporary variables for combination. */
  uint32_t *t0 = tmp;
  uint32_t *t1 = tmp + n1;
  uint32_t *t2 = t1 + n1;

  /* Temporary array for recursive Karatsuba calls,
     starting from (tmp + 3 * n1) */
  tmp = t2 + n1;

  /* r_low = karatsuba(a_low, b_low) */
  poly_mul(r, a, b, n1, tmp);
  /* r_high = karatsuba(a_high, b_high) */
  poly_mul(r2, a2, b2, n2, tmp);

  /* r_mid = karatsuba(a_low + a_high, b_low + b_high) */
  for (i = 0; i < n2; i++) {
    t0[i] = r1[i] ^ r2[i];
    t1[i] = a[i] ^ a2[i];
    t2[i] = b[i] ^ b2[i];
  }
  if (n1 != n2) {
    t0[i] = r1[i] ^ r2[i];
    t1[i] = a[i];
    t2[i] = b[i];
  }
  poly_mul(r1, t1, t2, n1, tmp);

  /* combination */
  for (i = 0; i < (n2 << 1) - n1; i++) {
    r1[i] ^= t0[i] ^ r[i];
    r2[i] ^= t0[i] ^ r3[i];
  }
  if (n1 != n2) {
    r1[i] ^= t0[i] ^ r[i];
    r2[i] ^= t0[i];
    i += 1;
    r1[i] ^= t0[i] ^ r[i];
    r2[i] ^= t0[i];
  }
}

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
    const unsigned int n, uint32_t *tmp) {
  if (n <= EXPD_MUL_THRES) {
    switch (n) {
      case 1:
        poly_mul1(r, a[0], b[0]);
        return;
      case 2:
        kara_mul2(r, a, b, tmp);
        return;
      case 3:
        kara_mul3(r, a, b, tmp);
        return;
      case 4:
        kara_mul4(r, a, b, tmp);
        return;
      case 5:
        kara_mul5(r, a, b, tmp);
        return;
      case 6:
        kara_mul6(r, a, b, tmp);
        return;
      default:
        return;
    }
  }
  kara_mul(r, a, b, n, tmp);
}

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
    const unsigned int n, uint32_t *tmp) {
  const uint32_t *a1 = a + n;           /* split a into two parts */
  uint32_t *r1 = r + n;
  uint32_t *t0 = tmp;
  tmp += n;

  poly_mul(r, a, b, n, tmp);
  memcpy(t0, r1, n * sizeof(uint32_t));
  poly_mul(r1, a1, b, n, tmp);

  for (unsigned int i = 0; i < n; i++) r1[i] ^= t0[i];
}


/*============================================================================*\
          Functions for polynomial modular reduction with 32-bit words
\*============================================================================*/

/*******************************************************************************
  The fast reduction algorithm for sparse polynomials is taken from
  the boost library (https://www.boost.org), with some optimizations for
  the MT19937 minimal polynomial.
  The original boost polynomial codes (random/detail/polynomial.hpp)
  are distributed under the Boost Software License:
 
  Copyright Steven Watanabe 2014
  Distributed under the Boost Software License, Version 1.0.
  (See copy at https://www.boost.org/LICENSE_1_0.txt)

*******************************************************************************/

/* Pre-computed properties of the minimal polynomial (phi) for MT19937.
 * ref: https://doi.org/10.1007/978-3-540-85912-3_26
 * These lookup tables are for fast modular reductions with divisor phi. */
#define PHI_NUM_NBITS           134
#define PHI_NUM_BLOCK           34
#define MT19937_POLY_LEN        19937

static const unsigned int phi_bit_pos[PHI_NUM_NBITS] = {
  0,1189,1416,1585,1643,1870,2493,2773,3000,3227,3454,3681,3908,4135,
  4362,4753,5661,6337,6569,7129,7477,7525,7583,7752,7979,8206,
  9505,9901,9969,10128,10693,10761,10920,11089,11147,11157,11215,11321,
  11374,11384,11485,11611,11712,11717,11838,11881,11944,11997,12277,12335,
  12393,12504,12509,12620,12673,12731,12736,12789,12905,12958,12963,13137,
  13185,13190,13243,13301,13412,13528,13533,13639,13697,13760,13813,13866,
  14093,14151,14209,14320,14325,14436,14547,14552,14605,14721,14774,14779,
  14953,15001,15006,15059,15117,15228,15344,15349,15455,15513,15576,15629,
  15682,15909,15967,16025,16136,16141,16252,16363,16368,16421,16537,16590,
  16595,16817,16822,16875,16933,17044,17160,17271,17329,17445,17498,17725,
  17783,17841,17952,18068,18179,18237,18406,18633,18691,18860,19087,19314
};

static const unsigned int phi_block_pos[PHI_NUM_BLOCK] = {
  39875,39252,38629,38006,37383,36760,36137,35514,34891,34268,33645,33022,
  32399,31776,31153,30530,29907,29284,28661,28038,27415,26792,26169,
  25546,24923,24300,23677,23054,22431,21808,21185,20562,19939,19937
};

/******************************************************************************
Function `copy_bits`:
  Copy a segment of bits from a polynomial to the beginning of the other one.
Arguments:
  * `r`:        pointer to the result, with sufficient space;
  * `a`:        pointer to the source, which must not overlap with `r`;
  * `start`:    the position of the first bit in `a` to be copied;
  * `end`:      the position right after the last bit to be copied.
******************************************************************************/
static inline void copy_bits(uint32_t *restrict r, const uint32_t *restrict a,
    const unsigned int start, const unsigned int end) {
  const unsigned int left = MOD_NBIT(start);    /* skipped bits for 1st word */
  const unsigned int right = WORD_SIZE - left;
  const unsigned int len = end - start;
  const unsigned int n = DIV_NBIT(len); /* number of full output words */
  unsigned int i;
  a += DIV_NBIT(start);

  if (left) {
    for (i = 0; i < n; i++)
      r[i] = (a[i] >> left) | (a[i+1] << right);
  }
  else memcpy(r, a, MUL_NBIT(n));

  if ((i = MOD_NBIT(len))) {            /* the last word is used partially */
    r[n] = a[n] >> left;
    if (left && MOD_NBIT(end)) r[n] |= (a[n+1] << right);
    r[n] &= ((UINT32_C(1) << i) - 1);
  }
}

/******************************************************************************
Function `shifted_add`:
  Compute r += a << shift.
Arguments:
  * `r`:        pointer to the result, with at least `n` + 1 words;
  * `a`:        pointer to the addend, with `n` words;
  * `n`:        the length of `a`;
  * `shift`:    the number of bits to be shifted before adding, which must
                be smaller than the word size.
******************************************************************************/
static inline void shifted_add(uint32_t *r, const uint32_t *a,
    const unsigned int n, const unsigned int shift) {
  if (shift == 0) {
    for (unsigned int i = 0; i < n; i++) r[i] ^= a[i];
    return;
  }

  const unsigned int right = WORD_SIZE - shift;
  uint32_t prev = 0;
  for (unsigned int i = 0; i < n; i++) {
    r[i] ^= (a[i] << shift) | (prev >> right);
    prev = a[i];
  }
  r[n] ^= (prev >> right);
}

/******************************************************************************
Function `poly_mod_phi`:
  Compute r %= phi, where phi is the minimal polynomial for MT19937.
Arguments:
  * `r`:        pointer to the result;
  * `tmp`:      a temporary array, with at least the size of one pre-defined
                block.
******************************************************************************/
void poly_mod_phi(uint32_t *r, uint32_t *tmp) {
  for (int i = 0; i < PHI_NUM_BLOCK - 1; i++) {
    unsigned int start = phi_block_pos[i + 1];
    unsigned int end = phi_block_pos[i];
    unsigned int size = DIV_NBIT(end - start + WORD_SIZE - 1);

    copy_bits(tmp, r, start, end);
    for (int j = 0; j < PHI_NUM_NBITS; j++) {
      unsigned int pos = phi_bit_pos[j] + start - MT19937_POLY_LEN;
      shifted_add(r + DIV_NBIT(pos), tmp, size, MOD_NBIT(pos));
    }
    shifted_add(r + DIV_NBIT(start), tmp, size, MOD_NBIT(start));
  }
}

