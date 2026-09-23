
# prand: parallel random number generator

![GitHub](https://img.shields.io/github/license/cheng-zhao/prand.svg)

## Table of Contents

-   [Introduction](#introduction)
-   [Compilation and linking](#compilation-and-linking)
-   [Getting started](#getting-started)
    -   [Random number generator initialisation](#random-number-generator-initialisation)
    -   [Sampling a uniform distribution](#sampling-a-uniform-distribution)
    -   [Sampling a Gaussian distribution](#sampling-a-gaussian-distribution)
    -   [Revising random states](#revising-random-states)
    -   [Releasing memory](#releasing-memory)
    -   [Error handling](#error-handling)
    -   [Examples](#examples)
 -   [Random number generator algorithms](#random-number-generator-algorithms)
 -   [Acknowledgements](#acknowledgements)
 -   [References](#references)

## Introduction

This is an efficient library written in C, for generating _pseudo_ random numbers with multiple streams. In particular, for a given random sequence, multiple streams can be generated from _user-defined_ equally spaced starting points. This is to ensure that for a parallel version of a program, the starting points for different processors/threads can be chosen in a way that the random number sequences are exactly the same as in a serial version of the same program.

For instance, to generate the random  sequence {_r_<sub>1</sub>, _r_<sub>2</sub>, &hellip;, _r_<sub>_n_&times;_m_</sub>},  one can use _n_ threads starting from _r_<sub>_j_&times;_m_+1</sub>, for _j_ = 0, 1, &hellip;, _n_&minus;1, and generate only _m_ numbers for each thread. This is equivalent to generating all (_n_ &times; _m_) numbers in sequence with a single stream.

To this end, efficient jump-ahead algorithms are implemented in this library, to permit arbitrary choices of the starting points of random streams, with relatively low computational costs. And this is hopefully useful for Monte-Carlo simulations that are massively parallelised and require strict reproducibility, especially if the same number of processors/threads is not always available.

This library is compliant with the ISO C99 standard, thread-safe, and does not rely on any external libraries. It is written by Cheng Zhao (&#36213;&#25104;), and is distributed under the [MIT license](LICENSE.txt).

<sub>[\[TOC\]](#table-of-contents)</sub>

## Compilation and linking

On most unix-like systems, the library can be simply compiled with

```console
$ make
$ make install
```

By default a static library `libprand.a` is created in the `lib` subfolder, and a header file `prand.h` is copied to the `include` subfolder, of the current working directory. One can change the `PREFIX` entry in [Makefile](Makefile#L5) to customise the installation path of the library.

To link the library with a program, one has to add the `-lprand` flag for the compilation. And if this library is not installed in the default path for system libraries, the `-I` and `-L` options are also necessary for specifying the path to the header and library files. An example of the [Makefile](example/Makefile) for linking prand is provided in the [example](example) folder.

<sub>[\[TOC\]](#table-of-contents)</sub>

## Getting started

### Random number generator initialisation

Before calling any random number generation functions, the interface of the random number generator has to be initialised by

```c
prand_t *prand_init(const prand_rng_enum type, const uint64_t seed,
    const unsigned int nstream, const uint64_t step, int *err)
```

The arguments are:
-   `type`: a pre-defined enumerator indicating the random number generation algorithm;
-   `seed`: a positive integer for initialising the random number generator;
-   `nstream`: number of streams that can be sampled in parallel;
-   `step`: the length of number sequences between adjacent streams;
-   `err`: a variable for storing the error status, see [Error handling](#error-handling).

If `seed` is set to `0`, then a warning is generated, and the default random seed for the given generator will be used. And if a generator accepts only 32-bit seeds, the supplied `seed` value will be truncated, i.e., reduced modulo 2<sup>32</sup>. Once the interface is initialised with multiple streams, the _i_-th (_i_ = 0, 1, &hellip;, `nstream`&minus;1) stream starts generating numbers from _r_<sub>_i_&times;`step`&plus;1</sub> of the sequence {_r_<sub>1</sub>, _r_<sub>2</sub>, &hellip;}.

The implemented random number generation algorithms, as well as the corresponding seed values and maximum allowed lengths for a single jump ahead operation are listed below.

| Algorithm                                       | Enumerator for `type` | `seed` | Maximum `step`[*](#foot1) |
|-------------------------------------------------|-----------------------|--------|---------------------------|
| MRG32k3a<sup>[\[1\]](#ref1)</sup>               | `PRAND_RNG_MRG32K3A` | 32-bit | 2<sup>63</sup>&minus;1    |
| Mersenne Twister 19937<sup>[\[2\]](#ref2)</sup> | `PRAND_RNG_MT19937`  | 32-bit | 2<sup>63</sup>&minus;1    |

<sub><span id="foot1">*</span> This limitation is only for a single jump ahead operation. The total skipped length can be larger than this value if jumping ahead for multiple times (see [Revising random states](#revising-random-states)).</sub>

All the random number generation routines rely on the interface created by this `prand_init` function. As an example, the interface can be initialised as follows

```c
int err = 0;
prand_t *rng = prand_init(PRAND_RNG_MRG32K3A, 1, 8, 10000, &err);
```

<sub>[\[TOC\]](#table-of-contents)</sub>

### Sampling a uniform distribution

If there is only one stream, one can generate a new random integer by

```c
long x = rng->get(rng->state);                    /* for single stream */
```
Here `rng` is the interface initialised in the previous section. While if there are multiple streams, generating a random integer for the _i_-th stream should be done with

```c
long x = rng->get(rng->state_stream[i]);          /* for multiple streams */
```

The minimum and maximum integer values that can be generated are stored in `rng->min` and `rng->max`, respectively. These values can be used with the `rng->get` function to sample uniform numbers in a custom range. For instance, we have provided the function for generating a double-precision floating-point number uniformly distributed in the range [0, 1):

```c
double z = rng->get_double(rng->state);           /* for single stream */
double z = rng->get_double(rng->state_stream[i]); /* for multiple streams */
```

And these functions evaluate essentially `z = x / (rng->max + 1)`.

Similarly, the following functions are for sampling a uniform double-precision floating-point number in the range (0, 1):

```c
double z = rng->get_double_pos(rng->state);           /* for single stream */
double z = rng->get_double_pos(rng->state_stream[i]); /* for multiple streams */
```

They are equivalent to `z = (x + 1) / (rng->max + 2)`.

<sub>[\[TOC\]](#table-of-contents)</sub>

### Sampling a Gaussian distribution

It is also a common practice to generate random numbers that follows a Gaussian probability distribution. There are several different ways of sampling Gaussian random numbers from a uniform distribution<sup>[\[3\]](#ref3)</sup>. However, to ensure the reproducibility with multiple streams, rejection methods are not appropriate, including the fast and robust Ziggurat method<sup>[\[4\]](#ref4)</sup>. Instead, we recommend the simple [Box&ndash;Muller transform](https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform) method<sup>[\[5\]](#ref5)</sup> for generating random numbers following the Gaussian distribution.

<sub>[\[TOC\]](#table-of-contents)</sub>

### Revising random states

It is sometimes necessary to jump ahead after having sampled some numbers, either for a specific stream or all streams. This can be done with the following functions provided with the random number generator interface:

```c
rng->jump(void *state, const uint64_t step, int *err);
rng->jump_all(rng, const uint64_t step, int *err);
```

With `rng->jump` only the stream associated with `state` is altered. And `state` can be either `rng->state` or `rng->state_stream[i]`, depending on how many streams have been initialised. While for `rng->jump_all`, the first argument must be the interface itself. And once it is called, all the streams will move forward with the length `step`.

The states of the interface can also be reset, with another random seed and/or jump ahead with a new step size from the initial starting point of the sequence:

```c
rng->reset(void *state, const uint64_t seed, const uint64_t step, int *err);
rng->reset_all(rng, const uint64_t seed, const uint64_t step, int *err);
```

Here, `rng->reset` resets the status of a stream indicated by `state`, with a jump ahead `step` counted from the starting point of the random number sequence initialised by `seed`. And `rng->reset_all` is equivalent of re-initialising the random number generator interface without releasing and re-allocating memory, albeit it is not possible to change the generation algorithm and number of streams with this function.

<sub>[\[TOC\]](#table-of-contents)</sub>

### Releasing memory

Once the random number generator is not needed anymore, the interface has to be deconstructed to release the allocated memory, by simply calling

```c
void prand_destroy(prand_t *rng);
```

<sub>[\[TOC\]](#table-of-contents)</sub>

### Error handling

An integer `err` has to be supplied to some of the functions for storing error status. And two macros, `PRAND_IS_ERROR(err)` and `PRAND_IS_WARN(err)`, are defined for checking whether there is an error or warning message, respectively. Furthermore, the following function returns a string describing the problem in more detail:

```c
const char *prand_errmsg(const int err);
``` 

<sub>[\[TOC\]](#table-of-contents)</sub>

### Examples

An example for the usage of this library is provided in the [example](example) folder.

It initialises two versions of the random number generator with the same seed, one with a single stream and the other with multiple streams. And the first floating-point number generated by the multiple streams are compared with the ones generated in sequence by the single stream version.

<sub>[\[TOC\]](#table-of-contents)</sub>

## Random number generator algorithms

| Algorithm              | Reference      | Period                    | Range of integers             | Jump ahead algorithm |
|------------------------|----------------|---------------------------|-------------------------------|-----------------------|
| MRG32k3a               | [\[1\]](#ref1) | 2<sup>191</sup>           | [0, 2<sup>32</sup>&minus;209] | [\[6\]](#ref6)         |
| Mersenne Twister 19937 | [\[2\]](#ref2) | 2<sup>19937</sup>&minus;1 | [0, 2<sup>32</sup>&minus;1]   | [\[7\]](#ref7)         |

<sub>[\[TOC\]](#table-of-contents)</sub>

## Acknowledgements

This library benefits from the following open-source projects:
-   [The BOOST C++ Library](https://www.boost.org/)
-   [https://github.com/vigna/MRG32k3a](https://github.com/vigna/MRG32k3a)
-   [https://github.com/cslarsen/mersenne-twister](https://github.com/cslarsen/mersenne-twister)

<sub>[\[TOC\]](#table-of-contents)</sub>

## References

<span id="ref1">\[1\]</span> L’Ecuyer, 1999, [Good Parameters and Implementations for Combined Multiple Recursive Random Number Generators](https://doi.org/10.1287/opre.47.1.159), _Operations Research_, 47(1):159&ndash;164

<span id="ref2">\[2\]</span> Matsumoto & Nishimura, 1998, [Mersenne Twister: A 623-Dimensionally Equidistributed Uniform Pseudo-Random Number Generator](https://doi.org/10.1145%2F272991.272995), _ACM Trans. Model. Comput. Simul._, 8(1):3&ndash;30 ([home page](http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html))

<span id="ref3">\[3\]</span> Thomas, Luk, Leong & Villasenor, 2007, [Gaussian Random Number Generators](https://doi.org/10.1145/1287620.1287622), _ACM Comput. Surv._, 39(4):11&ndash;es

<span id="ref4">\[4\]</span> Marsaglia & Tsang, 2000, [The Ziggurat Method for Generating Random Variables](https://doi.org/10.18637/jss.v005.i08), _Journal of Statistical Software_, 5(8):1&ndash;7

<span id="ref5">\[5\]</span> Box & Muller, 1958, [A Note on the Generation of Random Normal Deviates](https://doi.org/10.1214/aoms/1177706645), _Ann. Math. Statist._, 29(2):610&ndash;611

<span id="ref6">\[6\]</span> Bradley, du Toit, Tong, Giles & Woodhams, 2011, [Parallelization Techniques for Random Number Generators](https://doi.org/10.1016/B978-0-12-384988-5.00016-4), GPU Computing Gems Emerald Edition, _Morgan Kaufmann_, 231&ndash;246

<span id="ref6">\[7\]</span> Haramoto, Matsumoto & L'Ecuyer, 2008, [A Fast Jump Ahead Algorithm for Linear Recurrences in a Polynomial Space](https://doi.org/10.1007/978-3-540-85912-3_26), Sequences and Their Applications &ndash; SETA 2008, _Springer Berlin Heidelberg_, 290&ndash;298

<sub>[\[TOC\]](#table-of-contents)</sub>

