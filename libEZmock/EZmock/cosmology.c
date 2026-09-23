/*******************************************************************************
* cosmology.c: this file is part of the EZmock library.

* EZmock: Effective Zel'dovich approximation mock generator.

* Github repository:
        https://github.com/cheng-zhao/EZmock

* Copyright (c) 2023 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "EZmock.h"
#include "structs.h"
#include "hypergeom.h"
#include "config.h"
#include "errmsg.h"
#include <math.h>

/*============================================================================*\
             Function for computing cosmological growth parameters
\*============================================================================*/

/*******************************************************************************
  Numerical evaluations of exact solutions of linear perturbations,
  in a flat-wCDM cosmology.
  In particular, the following parameters (solutions) are concerned:
     -  delta_m = a * hyp2f1((w-1)/2w, -1/(3w), 1-5/(6w), 1-1/Om(a))
     -  f = 1 + 3 * (w-1) / (6w-5) * (1-1/Om(a))
              * hyp2f1((3w-1)/2w, (3w-1)/(3w), 2-5/(6w), 1-1/Om(a))
              / hyp2f1((w-1)/2w, -1/(3w), 1-5/(6w), 1-1/Om(a))
  The allowed ranges of the relevant cosmological parameters are:
     -  0 < a <= 1
     -  w < -1/3
     -  0 < Om < 1
  Ref: https://doi.org/10.1088/1475-7516/2011/10/010
*******************************************************************************/

/******************************************************************************
Function `cosmo_growth`:
  Compute the linear growth factor and growth rate in a flat wCDM cosmology.
Arguments:
  * `a`:        the scale factor;
  * `ainit`:    reference scale factor for the normalization of growth factor;
  * `Om0`:      matter density parameter at present (z=0);
  * `w`:        dark energy equation of state;
  * `tol`:      the desired tolerance of the numerical evaluation;
  * `maxiter`:  the maximum allowed number of interations;
  * `D`:        the evaluated linear growth factor;
  * `f`:        the evaluated linear growth rate;
  * `err`:      integer storing the error code.
******************************************************************************/
static void cosmo_growth(const double a, const double ainit,
    const double Om0, const double w, const double tol, const int maxiter,
    double *D, double *f, int *err) {
  /* Validate the ranges of arguments. */
  if (a <= 0 || a > 1 || ainit <= 0 || ainit > 1 || Om0 <= 0 || Om0 > 1 ||
      w >= -1.0 / 3.0) {
    *D = *f = HUGE_VAL;
    *err = EZMOCK_ERR_PAR_COSMO; return;
  }

  double tl = (tol < DOUBLE_EPSILON) ? HYPERGEOM_DEFAULT_TOL : tol;
  int mx = (maxiter <= 0) ? HYPERGEOM_DEFAULT_MAXITER : maxiter;

  /* Compute the growth rate. */
  const double iw = 1 / w;
  const double f1 = 0.5 * (1 - iw);                     /* (w - 1) / (2w) */
  const double f2 = -0x1.5555555555555p-2 * iw;         /* -1 / (3w) */
  const double f3 = 1 - 0x1.aaaaaaaaaaaabp-1 * iw;      /* 1 - 5 / (6w) */
  const double f4 = (1 - 1 / Om0) * pow(a, -3 * w);

  double hg1, hg2;
  if (hyp2f1(f1, f2, f3, f4, tl, mx, &hg1) ||
      hyp2f1(f1 + 1, f2 + 1, f3 + 1, f4, tl, mx, &hg2)) {
    *D = *f = HUGE_VAL;
    *err = EZMOCK_ERR_COS_GROWTH; return;
  }

  *f = 1 + 3 * (w - 1) * f4 / (6 * w - 5) * hg2 / hg1;

  /* Compute the growth factor. */
  if (a == ainit) {
    *D = 1;
    return;
  }

  const double f5 = (1 - 1 / Om0) * pow(ainit, -3 * w);
  if (hyp2f1(f1, f2, f3, f5, tl, mx, &hg2)) {
    *D = HUGE_VAL;
    *err = EZMOCK_ERR_COS_GROWTH; return;
  }

  *D = (a * hg1) / (ainit * hg2);
}




double cosmo_growth_f_only(const double a, const double ainit,
    const double Om0, const double w, const double tol, const int maxiter) {
    // 验证参数范围（与之前相同）
    if (a <= 0 || a > 1 || ainit <= 0 || ainit > 1 || Om0 <= 0 || Om0 > 1 ||
        w >= -1.0 / 3.0) {
        return HUGE_VAL; // 返回错误值
    }
 
    double tl = (tol < DOUBLE_EPSILON) ? HYPERGEOM_DEFAULT_TOL : tol;
    int mx = (maxiter <= 0) ? HYPERGEOM_DEFAULT_MAXITER : maxiter;
 
    // 计算增长率 f 所需的参数（与之前相同）
    const double iw = 1 / w;
    const double f1 = 0.5 * (1 - iw);
    const double f2 = -0x1.5555555555555p-2 * iw;
    const double f3 = 1 - 0x1.aaaaaaaaaaaabp-1 * iw;
    const double f4 = (1 - 1 / Om0) * pow(a, -3 * w);
 
    double hg1, hg2;
    // 计算 hypergeometric 函数（假设 hyp2f1 返回 0 表示成功，非 0 表示失败）
    if (hyp2f1(f1, f2, f3, f4, tl, mx, &hg1) != 0 ||
        hyp2f1(f1 + 1, f2 + 1, f3 + 1, f4, tl, mx, &hg2) != 0) {
        return HUGE_VAL; // 返回错误值
    }
 
    // 计算并返回增长率 f（与之前相同，但不需要 D）
    double result_f = 1 + 3 * (w - 1) * f4 / (6 * w - 5) * hg2 / hg1;
    return result_f;
}








/******************************************************************************
Function `cosmo_growth_D_only`:
  只计算生长因子之比 D_A(a)/D_A(ainit)（平直 wCDM 精确解，不含辐射）。
  与 `cosmo_growth` 用的是同一组超几何解：
      D_A(a) ∝ a * hyp2f1((w-1)/(2w), -1/(3w), 1-5/(6w), (1-1/Om0) * a^{-3w})
  归一化约定为“未归一化生长因子”：D_A(a) → a 当 a → 0（早期物质为主阶段）。
  （注意这里的 w 是暗能量状态方程，与 Om0 一样只允许 flat + 常数 w。）
Arguments:
  * `a`:        the scale factor;
  * `ainit`:    归一化参考时刻的标度因子（结果即 D_A(a)/D_A(ainit)）;
  * `Om0`:      今天的物质密度参数（不含中微子，与 cosmo_growth 一致）;
  * `w`:        dark energy equation of state;
  * `tol`:      tolerance of the numerical evaluation (<=0 用默认);
  * `maxiter`:  maximum number of iterations (<=0 用默认).
Return:
  D_A(a)/D_A(ainit)；参数越界或超几何不收敛时返回 HUGE_VAL。

  2026-09-21 修复记录：原实现的分母错用了“导数型”参数组
  (f1+1, f2+1, f3+1) 且与分子共用 a 处的变量（那是算增长率 f 用的公式），
  返回的根本不是生长因子之比（Ωm=0.3175 时给出 2.49 而非 0.789）。
  此函数此前无调用者，故该 bug 没有污染过任何产物。
******************************************************************************/
double cosmo_growth_D_only(const double a, const double ainit,
    const double Om0, const double w, const double tol, const int maxiter) {
    // 验证参数范围（与之前相同）
    if (a <= 0 || a > 1 || ainit <= 0 || ainit > 1 || Om0 <= 0 || Om0 > 1 ||
        w >= -1.0 / 3.0) {
        return HUGE_VAL; // 返回错误值
    }

    double tl = (tol < DOUBLE_EPSILON) ? HYPERGEOM_DEFAULT_TOL : tol;
    int mx = (maxiter <= 0) ? HYPERGEOM_DEFAULT_MAXITER : maxiter;

    // 超几何函数的参数（与 cosmo_growth 相同）
    const double iw = 1 / w;
    const double f1 = 0.5 * (1 - iw);
    const double f2 = -0x1.5555555555555p-2 * iw;
    const double f3 = 1 - 0x1.aaaaaaaaaaaabp-1 * iw;
    // 分子与分母各自在 a / ainit 处求值（z 变量随 a 变化，必须分开算）
    const double z_a = (1 - 1 / Om0) * pow(a, -3 * w);

    if (a == ainit) return 1.;

    double hg1, hg2;
    // 计算 hypergeometric 函数（返回 0 表示成功，非 0 表示失败）
    if (hyp2f1(f1, f2, f3, z_a, tl, mx, &hg1) != 0) {
        return HUGE_VAL; // 返回错误值
    }
    const double z_ainit = (1 - 1 / Om0) * pow(ainit, -3 * w);
    if (hyp2f1(f1, f2, f3, z_ainit, tl, mx, &hg2) != 0) {
        return HUGE_VAL;
    }

    // 计算并返回生长因子之比 D_A(a)/D_A(ainit)
    return (a * hg1) / (ainit * hg2);
}






/*============================================================================*\
          Interface for setting cosmological parameters used by EZmock
\*============================================================================*/

/******************************************************************************
Function `EZmock_set_cosmology`:
  Set or compute structure growth parameters in a flat wCDM cosmology.
  除了 growth2/vfac，本改造版还在这里（动态地）算出 local PNG 注入用的
  omega_m / png_d0 / png_dplus，见下方“PNG 生长因子量”注释。
Arguments:
  * `ez`:       instance of the EZmock generator;
  * `pk_norm`:  linear P(k) renormalization parameter, i.e., (D(z)/D(z_pk))^2;
  * `fHa`:      the factor for computing peculiar velocity, i.e., f*H(a)*a/h;
  * `eval`:     if true, set the above 2 parameters directly, otherwise
                evaluating them in a flat-wCDM cosmology with the following
                parameters;
  * `z`:        redshift of the EZmock snapshot to be produced;
  * `z_pk`:     redshift at which the input linear power spectrum is normalized;
  * `Omega_m`:  matter (without neutrino) density parameter at present (z=0);
  * `Omega_nu`: neutrino density parameter at present (z=0);
  * `w`:        dark energy equation of state;
  * `err`:      integer storing the error code.
Return:
  Zero on success; non-zero on err.
******************************************************************************/
int EZmock_set_cosmology(EZMOCK *ez, const double pk_norm, const double fHa,
    const bool eval, const double z, const double z_pk, const double Omega_m,
    const double Omega_nu, const double w, int *err) {
  if (!err) return EZMOCK_ERR_ARG_ECODE;
  if (*err != EZMOCK_SUCCESS) return *err;
  if (!ez) return EZMOCK_ERR_ARG_EZ;
  EZMOCK_COSMO *cosmo = (EZMOCK_COSMO *) ez->cosmo;

  if (!eval) {  /* set the two parameters for EZmock generation directly */
    if (pk_norm <= 0 || fHa <= 0) return (*err = EZMOCK_ERR_PAR_COSMO);
    /* PNG 生长因子量仍要由宇宙学参数算出（见下），因此这里也做范围检查 */
    if (z < 0 || Omega_m <= 0 || Omega_m > 1 || w >= -1.0 / 3.0)
      return (*err = EZMOCK_ERR_PAR_COSMO);
    cosmo->growth2 = pk_norm;
    cosmo->vfac = fHa;
  }
  else {        /* evaluate the parameters for EZmock generation */
    if (z < 0 || z_pk < 0 || Omega_m <= 0 || Omega_m > 1 ||
        Omega_nu < 0 || Omega_nu >= 1 || w >= -1.0 / 3.0)
      return (*err = EZMOCK_ERR_PAR_COSMO);

    double a = 1 / (1 + z);
    double a_pk = 1 / (1 + z_pk);
    double D, f;
    cosmo_growth(a, a_pk, Omega_m, w, HYPERGEOM_DEFAULT_TOL,
        HYPERGEOM_DEFAULT_MAXITER, &D, &f, err);
    if (*err != EZMOCK_SUCCESS) return *err;

    const double Omn = Omega_m + Omega_nu;
    const double Hubble = 100 * sqrt((1 - Omn) * pow(a, -3 * (1 + w))
        + Omn * pow(a, -3));

    cosmo->growth2 = D * D;
    cosmo->vfac = f * Hubble * a;
  }

  /* ---- local PNG 生长因子量：动态计算（2026-09-21 改造）----
     取代先前写死的 D0=0.789246 / Dplus=1.650（以及 Beta 里的 Ωm=0.315）。
     约定（沿用 2LPT-PNG，见 perturb.c 里的 FNL 注入）：
       D_A(a)  = 未归一化生长因子（精确解 a*hyp2f1(...)，早期 D_A→a，不含辐射）；
       D0      = D_A(1) = g(z=0)，PNG 势能归一化 Beta 的锚点，与输出红移无关；
       Dplus   = D_A(1)/D_A(a_out) = D(z=0)/D(z_out)（D(0)=1 归一化），
                 位移场整体乘 1/Dplus，把 z=0 形状的 T(k) 场缩放到输出红移 z。
     D_A(1) 的绝对归一化用早期锚点技巧：a_ref << 1 时 D_A(a_ref) ≈ a_ref，
     故 D_A(1) = a_ref * [D_A(1)/D_A(a_ref)]（a_ref=1e-3 时误差 ~1e-10）。
     数值交叉验证（2026-09-21）：Ωm=0.3175 时 D0=0.789246093、
     Dplus(z=1)=1.650333193，与 ODE 积分一致到 1e-5，与含辐射的 CAMB 差
     1.4e-4；旧硬编码正是这两个值（分别舍入/截断到 0.789246 与 1.650）。
     注意：这两个量在此处一次性算好，之后改 REDSHIFT/OMEGA_M 需重新调用本函数。 */
  const double a_ref = 1e-3;              /* 早期锚点（<<1 即可，误差 ~a_ref^3） */
  const double a_out = 1 / (1 + z);       /* 输出红移对应的标度因子 */
  const double d_ratio = cosmo_growth_D_only(1.0, a_ref, Omega_m, w,
      HYPERGEOM_DEFAULT_TOL, HYPERGEOM_DEFAULT_MAXITER);
  const double dplus = cosmo_growth_D_only(1.0, a_out, Omega_m, w,
      HYPERGEOM_DEFAULT_TOL, HYPERGEOM_DEFAULT_MAXITER);
  if (!isfinite(d_ratio) || !isfinite(dplus))
    return (*err = EZMOCK_ERR_PAR_COSMO);
  cosmo->omega_m = Omega_m;
  cosmo->png_d0 = a_ref * d_ratio;
  cosmo->png_dplus = dplus;

#ifdef EZMOCK_DEBUG
  printf("\ngrowth2 = %lf    vfac = %lf\n", cosmo->growth2, cosmo->vfac);
  printf("Omega_m = %lf    D0 = %.9f    Dplus = %.9f\n",
      cosmo->omega_m, cosmo->png_d0, cosmo->png_dplus);
#endif
  return EZMOCK_SUCCESS;
}
