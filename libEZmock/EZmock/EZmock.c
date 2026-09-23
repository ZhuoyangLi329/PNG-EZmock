/*******************************************************************************
* EZmock.c: this file is part of the EZmock library.

* EZmock: Effective Zel'dovich approximation mock generator.

* Github repository:
        https://github.com/cheng-zhao/EZmock

* Copyright (c) 2023 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include <stdlib.h>
#include <math.h>
#include "EZmock.h"
#include "structs.h"
#include "config.h"
#include "errmsg.h"
#include "linear_pk.h"

#ifdef OMP
#include <omp.h>
#endif

/*============================================================================*\
                       Interfaces of the EZmock generator
\*============================================================================*/

/******************************************************************************
Function `EZmock_init`:
  Initialise the EZmock instance.
Arguments:
  * `Lbox`:     side length of the simulation box;
  * `Ngrid`:    number of grids per side for the density field;
  * `randgen`:  random number generation algorithm;
  * `seed`:     random seed;
  * `nthread`:  number of OpenMP threads, 0 for omp_get_max_threads();
  * `err`:      integer storing the error code.
Return:
  Interface of EZmock generator.
******************************************************************************/
EZMOCK *EZmock_init(const double Lbox, const int Ngrid, const int randgen,
    const uint64_t seed, const int nthread, int *err) {
  /* Validate arguments. */
  if (!err || *err != EZMOCK_SUCCESS) return NULL;
  if (Lbox <= 0) {
    *err = EZMOCK_ERR_ARG_LBOX; return NULL;
  }
  if (Ngrid <= 1 || (Ngrid & 1) || Ngrid > EZMOCK_MAX_GRID_SIZE) {
    *err = EZMOCK_ERR_ARG_NGRID; return NULL;
  }
  switch (randgen) {
    case PRAND_RNG_MRG32K3A:
    case PRAND_RNG_MT19937:
      break;
    default:
      *err = EZMOCK_ERR_ARG_RNG; return NULL;
  }
  if (seed == 0) {
    *err = EZMOCK_ERR_ARG_SEED; return NULL;
  }
  if (nthread < 0) {
    *err = EZMOCK_ERR_ARG_THREAD; return NULL;
  }

  EZMOCK *ez = calloc(1, sizeof(EZMOCK));
  if (!ez) {
    *err = EZMOCK_ERR_MEMORY; return NULL;
  }

  /* Setup anonymous structures. */
  ez->conf = ez->rng = ez->cosmo = ez->pk = ez->mesh = NULL;
  if (!(ez->conf = malloc(sizeof(EZMOCK_CONF))) ||
      !(ez->rng = malloc(sizeof(EZMOCK_RNG)))) {
    *err = EZMOCK_ERR_MEMORY;
    EZmock_destroy(ez); return NULL;
  }

  EZMOCK_CONF *conf = (EZMOCK_CONF *) ez->conf;
#ifdef OMP
  conf->nthread = (nthread == 0) ? omp_get_max_threads() : nthread;
#else
  conf->nthread = 1;
#endif
  EZMOCK_RNG *rng = (EZMOCK_RNG *) ez->rng;
  rng->rng = prand_init(randgen, seed, conf->nthread, 0, err);
  if (PRAND_IS_ERROR(*err) || PRAND_IS_WARN(*err)) {
    *err = EZMOCK_ERR_RNG_SET;
    EZmock_destroy(ez); return NULL;
  }
  rng->ran = randgen;
  rng->seed = seed;

  if (!(ez->cosmo = malloc(sizeof(EZMOCK_COSMO)))) {
    *err = EZMOCK_ERR_MEMORY;
    EZmock_destroy(ez); return NULL;
  }
  EZMOCK_COSMO *cosmo = (EZMOCK_COSMO *) ez->cosmo;
  cosmo->growth2 = cosmo->vfac = HUGE_VAL;
  /* local PNG 生长因子量同样先置为未设置状态（供 EZmock_set_cosmology 填充） */
  cosmo->omega_m = cosmo->png_d0 = cosmo->png_dplus = HUGE_VAL;

  if (!(ez->mesh = calloc(1, sizeof(EZMOCK_MESH)))) {
    *err = EZMOCK_ERR_MEMORY;
    EZmock_destroy(ez); return NULL;
  }

  EZMOCK_MESH *mesh = (EZMOCK_MESH *) ez->mesh;
  mesh->rho_replaced = mesh->psi_ref = false;
  mesh->psi[0] = mesh->psi[1] = mesh->psi[2] = NULL;
  mesh->rhok = mesh->rhok2 = NULL;
  mesh->rho = mesh->rhot = NULL;


  //
  mesh->phik_png = mesh->phik = NULL;

  mesh->phi_png = mesh->phi = NULL;

  mesh->inj[0] = mesh->inj[1] = mesh->inj[2] = NULL;

  conf->Lbox = Lbox;
  conf->Ng = Ngrid;
  conf->fnl = 0.0;    //默认 fNL = 0（高斯初条件）；需要 PNG 时由 EZmock_set_fnl 或配置文件设置
  conf->b_phi = 0.0;  //默认 b_φ = 0（不指示踪物层面注入）；由 EZmock_set_b_phi 设置
  conf->fnl_field = HUGE_VAL;  /* 未设置时 B_PHI=0 退回 fnl，否则关闭场层注入；
                                  由 EZmock_set_fnl_field 或配置文件设置 */

  return ez;
}

/******************************************************************************
Function `EZmock_set_fnl`:
  Set the local PNG parameter used by the enabled PNG mechanisms.
  Field-level injection defaults to this value only when B_PHI is zero.
Arguments:
  * `ez`:       instance of the EZmock generator;
  * `fnl`:      local PNG parameter fNL (0 for Gaussian IC).
******************************************************************************/
void EZmock_set_fnl(EZMOCK *ez, const double fnl) {
  if (!ez) return;
  EZMOCK_CONF *conf = (EZMOCK_CONF *) ez->conf;
  conf->fnl = fnl;
}

/******************************************************************************
Function `EZmock_set_b_phi`:
  Set the calibrated tracer-level PNG input coefficient B_PHI for the
  injection (Ainj = 2*fnl*b_phi)。0 表示关闭该注入（默认）。
Arguments:
  * `ez`:       instance of the EZmock generator;
  * `b_phi`:    coefficient calibrated for this implementation.
******************************************************************************/
void EZmock_set_b_phi(EZMOCK *ez, const double b_phi) {
  if (!ez) return;
  EZMOCK_CONF *conf = (EZMOCK_CONF *) ez->conf;
  conf->b_phi = b_phi;
}

/******************************************************************************
Function `EZmock_set_fnl_field`:
  Set the coefficient of the field-level quadratic term,
  phi_png = phi + fnl_field * phi^2。不调用时 B_PHI=0 使用 fnl，
  B_PHI 非零时使用 0；显式调用本函数可覆盖默认值。
Arguments:
  * `ez`:          instance of the EZmock generator;
  * `fnl_field`:   coefficient of the real-space quadratic term.
******************************************************************************/
void EZmock_set_fnl_field(EZMOCK *ez, const double fnl_field) {
  if (!ez) return;
  EZMOCK_CONF *conf = (EZMOCK_CONF *) ez->conf;
  conf->fnl_field = fnl_field;
}


/******************************************************************************
Function `EZmock_destroy`:
  Release memory allocated for the EZmock generator interface.
Arguments:
  * `ez`:       instance of the EZmock generator.
******************************************************************************/
void EZmock_destroy(EZMOCK *ez) {
  if (!ez) return;
  if (ez->conf) free(ez->conf);
  EZmock_rng_destroy(ez->rng);
  if (ez->cosmo) free(ez->cosmo);
  EZmock_pk_destroy(ez->pk);
  EZmock_mesh_destroy(ez->mesh);
  free(ez);
}
