/*******************************************************************************
* run_mock.c: this file is part of the EZmock program.

* EZmock: Effective Zel'dovich approximation mock generator.

* Github repository:
        https://github.com/cheng-zhao/EZmock

* Copyright (c) 2023 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "define.h"
#include "load_conf.h"
#include "read_pk.h"
#include "EZmock.h"
#include "structs.h"
#include "save_res.h"
#include <stdio.h>
#include <stdlib.h>

#ifdef OMP
#include <omp.h>
#endif

int main(int argc, char *argv[]) {
  CONF *conf;
  if (!(conf = load_conf(argc, argv))) {
    printf(FMT_FAIL);
    P_EXT("failed to load configuration parameters\n");
    return EZMOCK_ERR_CFG;
  }

  int err = 0;
  printf("Initializing the EZmock generator ...");
  fflush(stdout);

  EZMOCK *ez;
  if (!(ez = EZmock_init(conf->Lbox, conf->Ngrid, conf->rng, conf->seed,
      conf->nthread, &err))) {
    P_ERR("%s\n", EZmock_errmsg(err));
    printf(FMT_FAIL);
    P_EXT("failed to initialize the EZmock generator\n");
    conf_destroy(conf);
    return EZMOCK_ERR_INIT;
  }

  /* 把配置里的 local PNG 参数 fNL 传给生成器（0 表示高斯初条件） */
  EZmock_set_fnl(ez, conf->fnl);
  /* 示踪物层面 PNG 注入（B_PHI，Ainj=2*FNL*B_PHI）与场层面二次项系数
     （FNL_FIELD）；B_PHI=0 时不注入，行为与历史版本逐位一致 */
  EZmock_set_b_phi(ez, conf->b_phi);
  EZmock_set_fnl_field(ez, conf->fnl_field);

  printf(FMT_DONE);
  printf("Setting cosmological parameters for EZmock generation ...");
  fflush(stdout);

  if (EZmock_set_cosmology(ez, conf->growth2, conf->vfac, conf->eval_growth,
      conf->redshift, conf->zpk, conf->omega_m, conf->omega_nu, conf->eos_w,
      &err)) {
    P_ERR("%s\n", EZmock_errmsg(err));
    printf(FMT_FAIL);
    P_EXT("failed to set cosmological parameters for EZmock generation\n");
    conf_destroy(conf); EZmock_destroy(ez);
    return EZMOCK_ERR_COSMO;
  }

  /* ---- 中和 cosmo->growth2，并打印生长因子溯源（2026-09-21 改造）----
     本改造版约定 LINEAR_PK 是 z=0 的转移函数 T(k)（如 Tk_0.txt），红移依赖
     完全由 perturb.c 里的 Dplus = D(0)/D(z_out) 承担；而 EZmock_setup_linear_pk()
     还会把输入表乘以 growth2 = (D(z)/D(z_pk))^2（stock 版 P(k) 流程的约定），
     两者叠加即为"生长因子数两遍"：REDSHIFT=1、REDSHIFT_PK=0 时会让位移幅度
     低 (D(1)/D(0))^2 ≈ 0.37 倍。因此这里强制 growth2 = 1。
     （旧版把 REDSHIFT 固定为 0 恰好绕开了这个坑，现在 REDSHIFT 才真正可用；
       REDSHIFT_PK 对本改造版不再有意义。） */
  {
    EZMOCK_COSMO *cosmo = (EZMOCK_COSMO *) ez->cosmo;
    if (conf->eval_growth && conf->zpk != 0)
      printf(FMT_WARN "\n  note: 本改造版忽略 REDSHIFT_PK=%g"
          "（输入表必须是 z=0 的 T(k)，红移换算由 Dplus 负责）\n", conf->zpk);
    if (!conf->eval_growth)
      printf(FMT_WARN "\n  note: 配置里显式给了 GROWTH_PK/VELOCITY_FAC，"
          "其中 GROWTH_PK 会被强制置 1（T(k) 流程不做输入表缩放），"
          "生长因子量仍按 OMEGA_M/DE_EOS_W/REDSHIFT 计算\n");
    cosmo->growth2 = 1.0;
    printf("\n  PNG growth factors: Omega_m = %g, D0 = %.9f, Dplus = %.9f, "
        "Beta = %.12e, vfac = %.6f\n", cosmo->omega_m, cosmo->png_d0,
        cosmo->png_dplus, 1.5 * cosmo->omega_m / (2998. * 2998.) /
        cosmo->png_d0, cosmo->vfac);
  }

  printf(FMT_DONE);

  PK *pk = NULL;
  if (!(pk = read_pk(conf)) ||
      EZmock_setup_linear_pk(ez, pk->k, pk->n, pk->Plin, pk->Pnw, conf->bao_mod,
      conf->logint, &err)) {
    if (err) P_ERR("%s\n", EZmock_errmsg(err));
    printf(FMT_FAIL);
    P_EXT("failed to setup the input power spectrum\n");
    conf_destroy(conf); EZmock_destroy(ez); pk_destroy(pk);
    return EZMOCK_ERR_PK;
  }
  pk_destroy(pk);

  printf("Generating the EZmock density field ...\n");
  
  fflush(stdout);
  

  if (EZmock_create_dens_field(ez, NULL, false, NULL, conf->fixamp, conf->iphase,
      &err)) {
    
    P_ERR("%s\n", EZmock_errmsg(err));
    printf(FMT_FAIL);
    P_EXT("failed to construct the density field\n");
    conf_destroy(conf); EZmock_destroy(ez);
    return EZMOCK_ERR_RHO;
  }

  


  printf(FMT_DONE);
  printf("Populating EZmock tracers ...");
  fflush(stdout);

  const real params[] = {
    conf->rho_c,
    conf->rho_exp,
    conf->pdf_base,
    conf->sigv
  };

  size_t ntracer;
  real *x, *y, *z, *vx, *vy, *vz;

  if (EZmock_populate_tracer(ez, params, conf->Ngal, conf->particle, &ntracer,
      &x, &y, &z, &vx, &vy, &vz, &err)) {
    P_ERR("%s\n", EZmock_errmsg(err));
    printf(FMT_FAIL);
    P_EXT("failed to sample EZmock tracers\n");
    conf_destroy(conf); EZmock_destroy(ez);
    return EZMOCK_ERR_RHO;
  }

  printf(FMT_DONE);
  EZMOCK_COSMO *cosmo = (EZMOCK_COSMO *) ez->cosmo;
  conf->growth2 = cosmo->growth2;
  conf->vfac = cosmo->vfac;
  EZmock_destroy(ez);

  if (save_box(conf, x, y, z, vx, vy, vz, ntracer)) {
    printf(FMT_FAIL);
    P_EXT("failed to save the output tracer catalog\n");
    conf_destroy(conf);
    free(x); free(y); free(z);
    free(vx); free(vy); free(vz);
    return EZMOCK_ERR_SAVE;
  }

  conf_destroy(conf);
  free(x); free(y); free(z);
  free(vx); free(vy); free(vz);
  return 0;
}
