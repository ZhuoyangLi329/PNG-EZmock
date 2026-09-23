/*******************************************************************************
* perturb.c: this file is part of the EZmock library.

* EZmock: Effective Zel'dovich approximation mock generator.

* Github repository:
        https://github.com/cheng-zhao/EZmock

* Copyright (c) 2023 Cheng Zhao <zhaocheng03@gmail.com> [GPLv3 license]
 
*******************************************************************************/
/* Macros for the template functions. */
#if !defined(EZMOCK_PT_WHITENOISE) && defined(EZMOCK_PT_LOGPK) && \
  defined(EZMOCK_PT_FIXAMP) && defined(EZMOCK_PT_IPHASE)



/* Macros for generating function names. */
#ifndef CONCAT_FNAME
  #define CONCAT_FNAME(a,b,c,d)         a##b##c##d
#endif

#ifndef EZMOCK_PT_FUNCNAME
  #define EZMOCK_PT_FUNCNAME(a,b,c,d)   CONCAT_FNAME(a,b,c,d)
#endif

/*============================================================================*\
                             Definition validation
\*============================================================================*/

#ifdef EZMOCK_LOGPK_NAME
  #undef EZMOCK_LOGPK_NAME
#endif
#ifdef EZMOCK_FIXAMP_NAME
  #undef EZMOCK_FIXAMP_NAME
#endif
#ifdef EZMOCK_IPHASE_NAME
  #undef EZMOCK_IPHASE_NAME
#endif

#if     EZMOCK_PT_LOGPK == 1
  #define EZMOCK_LOGPK_NAME     _logpk
#elif   EZMOCK_PT_LOGPK == 0
  #define EZMOCK_LOGPK_NAME
#else
  #error "unexpected definition of `EZMOCK_PT_LOGPK`"
#endif

#if     EZMOCK_PT_FIXAMP == 1
  #define EZMOCK_FIXAMP_NAME    _fixamp
#elif   EZMOCK_PT_FIXAMP == 0
  #define EZMOCK_FIXAMP_NAME
#else
  #error "unexpected definition of `EZMOCK_PT_FIXAMP`"
#endif

#if     EZMOCK_PT_IPHASE == 1
  #define EZMOCK_IPHASE_NAME    _iphase
#elif   EZMOCK_PT_IPHASE == 0
  #define EZMOCK_IPHASE_NAME
#else
  #error "unexpected definition of `EZMOCK_PT_IPHASE`"
#endif


/*============================================================================*\
                   Function for displacement field generation
\*============================================================================*/

/******************************************************************************
Function `EZmock_ZA_disp<EZMOCK_LOGPK_NAME><EZMOCK_FIXAMP_NAME>
    <EZMOCK_PT_IPHASE>`:
  Generate the Zel'dovich displacement field given the input power spectrum.
Arguments:
  * `ez`:       instance of the EZmock generator;
  * `plan`:     FFTW plan;
  * `err`:      integer storing the error code.
******************************************************************************/
static void EZMOCK_PT_FUNCNAME(EZmock_ZA_disp, EZMOCK_LOGPK_NAME,
    EZMOCK_FIXAMP_NAME, EZMOCK_IPHASE_NAME)
    (EZMOCK *ez, FFT_PLAN plan, FFT_PLAN plan2,int *err) {
  /* Initialise the random number generator. */
  EZMOCK_RNG *erng = (EZMOCK_RNG *) ez->rng;
  prand_t *rng = erng->rng;

  EZMOCK_CONF *conf = (EZMOCK_CONF *) ez->conf;
  EZMOCK_PK *pk = (EZMOCK_PK *) ez->pk;
  EZMOCK_MESH *mesh = (EZMOCK_MESH *) ez->mesh;

  const int Ngh = conf->Ng >> 1;
  const int Ngk = Ngh + 1;
  const int Ngrid=conf->Ng*conf->Ng*conf->Ng;

  /* ---- 初条件形状参数：仍是写死的（来源待考，与 Tk_0.txt 配套使用）----
     initial_power = Anorm * k^PrimordialIndex，即原初密度场振幅的平方
     （P_prim ∝ k^ns，ns=0.9624），再由 Pfac=(2π/L)^1.5 归一化到格点上。
     注意：若以后换输入表/换原初谱，这两个数需要一起复核。
     （原代码里还有一个未被使用的 k_pivot=0.05 和 h=0.6711，已删除。） */
  const double Anorm = 12275.369233;
  const double PrimordialIndex = 0.9624;

  /* ---- local PNG 生长因子量：2026-09-21 起由宇宙学/红移动态计算 ----
     全部来自 EZmock_set_cosmology() 填好的 cosmo（对应配置项
     OMEGA_M / DE_EOS_W / REDSHIFT，见 cosmology.c 里的注释）：
       Omega_m = Om0（配置 OMEGA_M，不含中微子）；
       D0      = D_A(1) = g(z=0)，未归一化生长因子（早期 D_A→a），
                 固定住 PNG 势能归一化 Beta 的物理含义，与输出红移无关；
       Dplus   = D_A(1)/D_A(a_out) = D(0)/D(z_out)，a_out = 1/(1+REDSHIFT)，
                 位移场整体除以它，把 z=0 形状的输入表缩放到输出红移。
     Beta 沿用 2LPT-PNG 新版约定（与 pyEZmock-master/nested_sampling/model.py 一致）：
       Beta = (3/2) * Om0 * H0^2 / D0 ，其中 H0 用 100 km/s/Mpc、k 用 h/Mpc，
       换算后即 1.5*Om0/2998^2/D0（2998 Mpc = c/(100 km/s)，h 恰好抵消；
       原写法里 UnitLength_in_cm 的因子在分子分母完全抵消，这里化简掉）。
     物理链条：phik = Beta * δ_init / k^2 是 PNG 势能（量纲检查见上面"Ho in units of
     h/Mpc and c=1"）；位移步 twb = Transfer/Dplus/Beta，Beta 在线性部分精确抵消，
     只在 PNG 二次项里以 Fnl*Beta/Dplus 的形式标定 mock 的 fNL 水平。 */
  EZMOCK_COSMO *cosmo = (EZMOCK_COSMO *) ez->cosmo;
  const double Omega_m = cosmo->omega_m;
  const double D0 = cosmo->png_d0;
  const double Dplus = cosmo->png_dplus;
  const double Beta = 1.5 * Omega_m / (2998. * 2998.) / D0;

  fftw_plan myplan1 = fftw_plan_dft_c2r_3d(conf->Ng, conf->Ng, conf->Ng, mesh->phik, mesh->phi, FFTW_ESTIMATE);
  fftw_plan myplan2 = fftw_plan_dft_r2c_3d(conf->Ng, conf->Ng, conf->Ng, mesh->phi, mesh->phik_png, FFTW_ESTIMATE);









#ifdef OMP
  const size_t pnum = ((size_t) conf->Ng * conf->Ng) / conf->nthread;
  const int rem = ((size_t) conf->Ng * conf->Ng) % conf->nthread;
#endif

  const double kfac = M_PI * 2 / conf->Lbox;
  printf("kfac=%f\n",kfac);
#if EZMOCK_PT_LOGPK == 1
  const double logkfac = log(kfac);
#endif
  const double Pfac = pow(2*M_PI/conf->Lbox, 1.5); //目前看来他和2lpt的Pfac不一样，但不需要改？

  /* First dimension: generate initial condition and compute psi[0]. */
#ifdef OMP
#pragma omp parallel num_threads(conf->nthread)
  {
    /* Distribute mesh grids to threads. */
    const int tid = omp_get_thread_num();
    const size_t pcnt = (tid < rem) ? pnum + 1 : pnum;
    const size_t istart = (tid < rem) ? pcnt * tid : pnum * tid + rem;
    const size_t iend = istart + pcnt;
    /* Jump ahead the random states. Note that i=j=k=0 is skipped. */
  #if EZMOCK_PT_FIXAMP == 1
    rng->reset(rng->state_stream[tid], erng->seed,
        (tid == 0) ? 0 : istart * Ngk - 1, err);
  #else
    rng->reset(rng->state_stream[tid], erng->seed,
        (tid == 0) ? 0 : (istart * Ngk - 1) << 1, err);
  #endif
    if (PRAND_IS_ERROR(*err)) *err = EZMOCK_ERR_RNG_JUMP;
#pragma omp barrier     /* synchronize `err` */

    /* Sample the Fourier space density field with OpenMP. */
    if (*err == EZMOCK_SUCCESS) {
      for (size_t idx_ij = istart; idx_ij < iend; idx_ij++) {
        int i = idx_ij / conf->Ng;
        int j = idx_ij % conf->Ng;
        double ki = (i <= Ngh) ? i : i - conf->Ng;
        int ni = (i == 0) ? 0 : conf->Ng - i;           /* index of -i */
        double kni = (ni <= Ngh) ? ni : ni - conf->Ng;
#else
    /* Sample the Fourier space density field sequentially. */
    for (int i = 0; i < conf->Ng; i++) {
      size_t idx_i = (size_t) i * conf->Ng;
      double ki = (i <= Ngh) ? i : i - conf->Ng;
      int ni = (i == 0) ? 0 : conf->Ng - i;             /* index of -i */
      double kni = (ni <= Ngh) ? ni : ni - conf->Ng;
      for (int j = 0; j < conf->Ng; j++) {
        size_t idx_ij = idx_i + j;
#endif
        double kj = (j <= Ngh) ? j : j - conf->Ng;
        int nj = (j == 0) ? 0 : conf->Ng - j;           /* index of -j */
        double kij = ki * ki + kj * kj;
        size_t idx0 = idx_ij * Ngk;

        /* Expand k loop and treat 0 and Nyquist frequencies separately. */
        /* k = 0 */
        if (idx0 == 0) {                /* i = j = k = 0 */
          //原代码这里是rhok，我改成phik
          mesh->phik[0][0] = mesh->phik[0][1] = 0;
          mesh->rhok2[0][0] = mesh->rhok2[0][1] = 0;
          


        }
        else {



          /* Sample randoms for the amplitude and phase first. */
#if EZMOCK_PT_FIXAMP == 0
  #ifdef OMP
          double amp = rng->get_double_pos(rng->state_stream[tid]);
          //amp=fabs(sin(7*i+6*j+5*0+1));


   

  #else
          double amp = rng->get_double_pos(rng->state);
          //amp=fabs(sin(7*i+6*j+5*0+1));





  #endif
#endif
#ifdef OMP
          double phase = rng->get_double(rng->state_stream[tid]) * 2 * M_PI;
          //phase=sin(i+2*j+3*0)*M_PI*2;
          



#else
          double phase = rng->get_double(rng->state) * 2 * M_PI;
          //phase=sin(i+2*j+3*0)*M_PI*2;




#endif

          

          /* Ensure conjugation on the k = 0 plane. */
          /* Fill the field with the lower half plane. */
          if ((i == ni && j <= Ngh) || (i != ni && i <= Ngh)) {
            size_t idx = idx0;
            size_t nidx = ((size_t) ni * conf->Ng + nj) * Ngk;  /* (-i,-j,0) */
            double ksq = kij;
#if EZMOCK_PT_LOGPK == 1
            double kmod = log(ksq) * 0.5 + logkfac;     /* log(sqrt(ksq)) */
#else
            double kmod = kfac * sqrt(ksq);             /* sqrt(ksq) */

            
#endif

            //原代码这里算功率谱，现在这里要算initial_power
            /* initial normalized power. */

            double initial_power = Anorm*(exp( PrimordialIndex * log(kmod) ));//原初扰动


#if EZMOCK_PT_FIXAMP == 0
            amp = -log(amp);
#endif
#if EZMOCK_PT_LOGPK == 1
  #if EZMOCK_PT_FIXAMP == 0
            initial_power = exp(initial_power * 0.5) * sqrt(amp);
  #else
            initial_power = exp(initial_power * 0.5);                           /* sqrt(P) */
  #endif
#else
  #if EZMOCK_PT_FIXAMP == 0
            initial_power = sqrt(initial_power * amp);
  #else
            initial_power = sqrt(initial_power);
  #endif
#endif
            initial_power *= Pfac; //注意这一步，initial_power已经变成了\delta_m
            //initial_power/=(ksq*kfac);



           
            double initial_potential=initial_power * Beta / pow(kmod,2); //原初扰动势能的强度



            /* Generate the Fourier space density. */
#if EZMOCK_PT_IPHASE == 1
            phase += M_PI;
#endif
 
            //将k空间的扰动势能放到格点上,这一步原代码跳步骤了(提前乘了i)，这个版本一步一步来
            mesh->phik[idx][0] = initial_potential*cos(phase);
            mesh->phik[idx][1] = initial_potential*sin(phase);

                                    //学习2LPT，满足下列条件的phik全部至0
            if(i == conf->Ng/2 || j == conf->Ng/2 )
            {
              mesh->phik[idx][0] = 0;
              mesh->phik[idx][1] = 0;
            }


            //源代码nidx也跳步骤了(提前乘了i)，所以共轭关系反掉了，这里一步一步来，就正常的共轭关系
            mesh->phik[nidx][0] = mesh->phik[idx][0];
            mesh->phik[nidx][1] = -mesh->phik[idx][1];

            //

  
          }             /* check of the lower half plane */
        }               /* i = j = k = 0 */

        /* 0 < k < k_ny */
        for (int k = 1; k <= ((conf->Ng - 1) >> 1); k++) {


  




          /* Sample randoms for the amplitude and phase first. */
#if EZMOCK_PT_FIXAMP == 0
  #ifdef OMP
          double amp = rng->get_double_pos(rng->state_stream[tid]);
          //amp=fabs(sin(7*i+6*j+5*k+1));




  #else
          double amp = rng->get_double_pos(rng->state);
          //amp=fabs(sin(7*i+6*j+5*k+1));




  #endif
#endif
#ifdef OMP
          double phase = rng->get_double(rng->state_stream[tid]) * 2 * M_PI;
          //phase=sin(i+2*j+3*k)*M_PI*2;




#else
          double phase = rng->get_double(rng->state) * 2 * M_PI;
          //phase=sin(i+2*j+3*k)*M_PI*2;



#endif
          

          size_t idx = idx0 + k;
          double ksq = kij + k * k;
#if EZMOCK_PT_LOGPK == 1
          double kmod = log(ksq) * 0.5 + logkfac;       /* log(sqrt(ksq)) */
#else
          double kmod = kfac * sqrt(ksq);               /* sqrt(ksq) */
          
#endif
            //原代码这里算功率谱，现在这里要算initial_power
            /* initial normalized power. */
            //double hkmod = kmod*h;
            double initial_power = Anorm*(exp( PrimordialIndex * log(kmod) ));//原初扰动

            //double primo_factor=hkmod;
            //initial_power*=primo_factor;

#if EZMOCK_PT_FIXAMP == 0
          amp = -log(amp);
#endif
#if EZMOCK_PT_LOGPK == 1
  #if EZMOCK_PT_FIXAMP == 0
          initial_power = exp(initial_power * 0.5) * sqrt(amp);
  #else
          initial_power = exp(initial_power * 0.5);                             /* sqrt(P) */
  #endif
#else
  #if EZMOCK_PT_FIXAMP == 0
          initial_power = sqrt(initial_power * amp);
  #else
          initial_power = sqrt(initial_power);
  #endif
#endif
          initial_power *= Pfac;
          //initial_power/=(ksq*kfac);


          /* Ho in units of h/Mpc and c=1, i.e., internal units so far  */

    
          double initial_potential=initial_power * Beta / pow(kmod,2); //原初扰动势能的强度




          /* Generate the Fourier space density. */
#if EZMOCK_PT_IPHASE == 1
          phase += M_PI;
#endif

          //打印i,j,k,phase
          //printf("i=%d,j=%d,k=%d,phase=%f\n",i,j,k,phase);

          //将k空间的扰动势能放到格点上,这一步原代码跳步骤了(提前乘了i)，这个版本一步一步来
          mesh->phik[idx][0] = initial_potential * cos(phase);
          mesh->phik[idx][1] = initial_potential * sin(phase);

                                  //学习2LPT，满足下列条件的phik全部至0
          if(i == conf->Ng/2 || j == conf->Ng/2 || k == conf->Ng/2)
          {
            mesh->phik[idx][0] = 0;
            mesh->phik[idx][1] = 0;
          }

        }               /* for k */

        /* k = k_ny */
        if ((conf->Ng & 1) == 0) {
          /* Sample randoms for the amplitude and phase first. */
#if EZMOCK_PT_FIXAMP == 0
  #ifdef OMP
          double amp = rng->get_double_pos(rng->state_stream[tid]);
          //amp=fabs(sin(7*i+6*j+5*Ngh+1));


  #else
          double amp = rng->get_double_pos(rng->state);
          //amp=fabs(sin(7*i+6*j+5*Ngh+1));


  #endif
#endif
#ifdef OMP
          double phase = rng->get_double(rng->state_stream[tid]) * 2 * M_PI;
          //phase=sin(i+2*j+3*Ngh)*M_PI*2;


#else
          double phase = rng->get_double(rng->state) * 2 * M_PI;
          //phase=sin(i+2*j+3*Ngh)*M_PI*2;

#endif

          /* Ensure conjugation on the k = k_ny plane. */
          /* Fill the field with the lower half plane. */
          if ((i == ni && j <= Ngh) || (i != ni && i <= Ngh)) {
            int k = Ngh;
            size_t idx = idx0 + k;
            size_t nidx = ((size_t) ni * conf->Ng + nj) * Ngk + k;
            double ksq = kij + k * k;
#if EZMOCK_PT_LOGPK == 1
            double kmod = log(ksq) * 0.5 + logkfac;     /* log(sqrt(ksq)) */
#else
            double kmod = kfac * sqrt(ksq);             /* sqrt(ksq) */
            
#endif

            //原代码这里算功率谱，现在这里要算initial_power
            /* initial normalized power. */
            double initial_power = Anorm*(exp( PrimordialIndex * log(kmod) ));//原初扰动


#if EZMOCK_PT_FIXAMP == 0
            amp = -log(amp);
#endif
#if EZMOCK_PT_LOGPK == 1
  #if EZMOCK_PT_FIXAMP == 0
            initial_power = exp(initial_power * 0.5) * sqrt(amp);
  #else
            initial_power = exp(initial_power * 0.5);                           /* sqrt(P) */
  #endif
#else
  #if EZMOCK_PT_FIXAMP == 0
            initial_power = sqrt(initial_power * amp);
  #else
            initial_power = sqrt(initial_power);
  #endif
#endif
            initial_power *= Pfac;
            //initial_power/=(ksq*kfac);



    
            double initial_potential=initial_power * Beta / pow(kmod,2); //原初扰动势能的强度



            /* Generate the Fourier space density. */
#if EZMOCK_PT_IPHASE == 1
            phase += M_PI;
#endif

            //将k空间的扰动势能放到格点上,这一步原代码跳步骤了(提前乘了i)，这个版本一步一步来
            mesh->phik[idx][0] = initial_potential*cos(phase);
            mesh->phik[idx][1] = initial_potential*sin(phase);

                                    //学习2LPT,k=k_ny的时候，phik=0
     
            mesh->phik[idx][0] = 0;
            mesh->phik[idx][1] = 0;
            


            //mesh->rhok2[idx][0] = mesh->phik[idx][0]*ki*DstartFnl/Beta;
            //mesh->rhok2[idx][1] = mesh->phik[idx][1]*ki*DstartFnl/Beta;

            //源代码nidx也跳步骤了(提前乘了i)，所以共轭关系反掉了，这里一步一步来，就正常的共轭关系
            mesh->phik[nidx][0] = mesh->phik[idx][0];
            mesh->phik[nidx][1] = -mesh->phik[idx][1];

            //mesh->rhok2[nidx][0] = mesh->phik[nidx][0] * kni*DstartFnl/Beta;
            //mesh->rhok2[nidx][1] = mesh->phik[nidx][1] * kni*DstartFnl/Beta;
          }     /* check of the lower half plane */
        }       /* Ng is even */
      }         /* for idx_ij (OMP) or for j (no OMP) */
    }           /* *err == EZMOCK_SUCCESS (OMP) or for i (no OMP) */

#ifdef OMP
  }             /* omp parallel */
  if (*err != EZMOCK_SUCCESS) return;
#endif









  


    /* The imaginary parts of white noise field should be 0 at:
     (0,0,0), (0,0,Ngh), (0,Ngh,0), (0,Ngh,Ngh),
     (Ngh,0,0), (Ngh,0,Ngh), (Ngh,Ngh,0), (Ngh,Ngh,Ngh) */
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
     //这个地方也需要改变,原代码注释说是imaginary parts of white noise field should be 0，
     //然而，实际上源代码把real parts设成0,这是因为源代码提前x了i，这里一步一步来，没有提前乘i，所以就把虚部改为0.

  if ((conf->Ng & 1) == 0) {
    size_t idx = Ngh;
    mesh->phik[idx][1] = 0;
    idx = (size_t) Ngh * Ngk;
    mesh->phik[idx][1] =  0;
    idx = (size_t) Ngh * (Ngk + 1);
    mesh->phik[idx][1] = 0;
    idx = (size_t) Ngh * conf->Ng * Ngk;
    mesh->phik[idx][1] = 0;
    idx = Ngh * ((size_t) conf->Ng * Ngk + 1);
    mesh->phik[idx][1] = 0;
    idx = (size_t) Ngh * (conf->Ng + 1) * Ngk;
    mesh->phik[idx][1] =  0;
    idx = ((size_t) Ngh * conf->Ng + Ngh) * Ngk + Ngh;
    mesh->phik[idx][1] = 0;
  }

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////





  fftw_complex *phik_copy = fftw_malloc(sizeof(fftw_complex) * conf->Ng * conf->Ng *Ngk);
  memcpy(phik_copy, mesh->phik, sizeof(fftw_complex) * conf->Ng * conf->Ng *Ngk );

  fftw_execute_dft_c2r(myplan1, phik_copy, mesh->phi);
  
  

  //FFT_EXEC_C2R(plan,mesh->phik,mesh->phi);




  //现在，所有的phi格点都有正确的初始高斯势了
  /* ---- 场层面 PNG 二次项（phi_png = phi + Fnl_field*phi^2）----
     fnl_field 未设置（HUGE_VAL）时退回 conf->fnl，即历史行为；
     B_PHI≠0 的标定模式下配置会把 fnl_field 显式置 0（关掉这个未标定的机制，
     由示踪物层面注入接管，见下面的注入块）。 */
  const double Fnl = conf->fnl;   //fNL 由配置项/接口设置（EZMOCK_CONF.fnl），0 表示纯高斯初条件
  const double Fnl_field = (conf->fnl_field == HUGE_VAL) ? Fnl : conf->fnl_field;
  //对phi进行变换，得到local形式的PNG phi

  printf("开始打印实空间phi\n\n");
  for(int i = 0; i < conf->Ng; i++)
    for(int j = 0; j < conf->Ng; j++)
      for(int k = 0; k < conf->Ng; k++)
            {
             size_t idx = (conf->Ng*i+j)*conf->Ng+k;


             mesh->phi_png[idx] = mesh->phi[idx] + Fnl_field * mesh->phi[idx]*mesh->phi[idx];
             //打印phi_png
             //printf("i=%d,j=%d,k=%d,phi_png=%f\n",i,j,k,1e10*mesh->phi_png[idx]);

   

            }
  //现在，所有的phi_png格点都有正确的local非高斯势了
  
  fftw_complex *phik_png_right = fftw_malloc(sizeof(fftw_complex) * conf->Ng * conf->Ng *Ngk);

  fftw_execute_dft_r2c(myplan2, mesh->phi_png, phik_png_right);



  //再傅里叶逆变换到实空间中
  //FFT_EXEC_R2C(plan2,mesh->phi_png,mesh->phik_png);

  //傅里叶逆变换需要normalize，因此，这个循环要把normailize的因子乘上(1/N^3),同时乘上i
  double normal_factor = 1.0/(conf->Ng*conf->Ng*conf->Ng);


  for (int i = 0; i < conf->Ng; i++) 
  {
    size_t idx_i = (size_t) i * conf->Ng;
    double ki = (i <= Ngh) ? i : i - conf->Ng;
    int ni = (i == 0) ? 0 : conf->Ng - i;             /* index of -i */
    double kni = (ni <= Ngh) ? ni : ni - conf->Ng;
    for (int j = 0; j < conf->Ng; j++) 
    {
      size_t idx_ij = idx_i + j;
      double kj = (j <= Ngh) ? j : j - conf->Ng;
      int nj = (j == 0) ? 0 : conf->Ng - j;           /* index of -j */
      double kij = ki * ki + kj * kj;
      size_t idx0 = idx_ij * Ngk;
      
        /* 0 <= k <= k_ny */
        for (int k = 0; k <Ngk; k++) 
        {
          size_t idx = idx0 + k;
          double ksq = kij + k * k;
          double kmod = kfac * sqrt(ksq);




  
          mesh->phik_png[idx][1] = phik_png_right[idx][0]*normal_factor;
          mesh->phik_png[idx][0] = -phik_png_right[idx][1]*normal_factor;

          
          
        }



    }

    //学习2lpt，phik_png[0].[0]和phik_png[0].[1]都为0
    mesh->phik_png[0][0] = 0;
    mesh->phik_png[0][1] = 0;

  }

  

  //打印phik_png的前几个模式
  //for(int i = 0; i < 3; i++)
    //for(int j = 0; j < 3; j++)
      //for(int k = 0; k < 3; k++)
            //{
             //size_t idx = (conf->Ng*i+j)*Ngk+k;

             //printf("i=%d,j=%d,k=%d,phik_png=%f,%f\n",i,j,k,1e10*mesh->phik_png[idx][0],1e10*mesh->phik_png[idx][1]);

            //}





















  
  //现在，所有的phik_png格点都有正确的local非高斯势了
  //再通过转移函数转换到物质扰动rhok中,注意在傅里叶空间中，第三个指标k的范围和实空间不一样了。



  for (int i = 0; i < conf->Ng; i++) 
  {
    size_t idx_i = (size_t) i * conf->Ng;
    double ki = (i < Ngh) ? i : i - conf->Ng;
    int ni = (i == 0) ? 0 : conf->Ng - i;             /* index of -i */
    for (int j = 0; j < conf->Ng; j++) 
    {
      size_t idx_ij = idx_i + j;
      double kj = (j < Ngh) ? j : j - conf->Ng;
      int nj = (j == 0) ? 0 : conf->Ng - j;           /* index of -j */
      double kij = ki * ki + kj * kj;
      size_t idx0 = idx_ij * Ngk;
      
        /* 0 <= k <= k_ny */
        for (int k = 0; k <Ngk; k++) 
        {
          size_t idx = idx0 + k;
          double ksq = kij + k * k;
          double kmod = kfac * sqrt(ksq);
          double Transfer = pk_interp(pk, kmod);

          // kmod==0,Transfer=0
          if (i==0 && j==0 && k==0)
          {
            Transfer=0;
          }


        
          double twb=Transfer/Dplus/Beta;
 

          mesh->rhok2[idx][0] = mesh->phik_png[idx][0]*ki*kfac*twb;
          mesh->rhok2[idx][1] = mesh->phik_png[idx][1]*ki*kfac*twb;
          mesh->rhok3[idx][0] = mesh->phik_png[idx][0]*kj*kfac*twb;
          mesh->rhok3[idx][1] = mesh->phik_png[idx][1]*kj*kfac*twb;
          mesh->rhok4[idx][0] = mesh->phik_png[idx][0]*kfac*k*twb;
          mesh->rhok4[idx][1] = mesh->phik_png[idx][1]*kfac*k*twb;
        }



    }

  }

  






  //有了rhok，我们就可以计算ZA displacement field了,分3个维度
  FFT_EXEC_C2R(plan, mesh->rhok2, mesh->psi[0]);
  FFT_EXEC_C2R(plan, mesh->rhok3, mesh->psi[1]);
  FFT_EXEC_C2R(plan, mesh->rhok4, mesh->psi[2]);


  /* ==================== 示踪物层面 PNG 注入（2026-09-21 新增）====================
     动机：场层面 phi -> phi + FNL*phi^2 对示踪物的 PNG 响应是"涌现"的——它经过
     密度场 -> PDF 映射后才变成示踪物偏置，强度依赖四旋钮/网格/方差，导致标定时
     需要的 FNL 不是物理 f_NL（固定四旋钮下要对上 Quijote f_NL=100 得开到 ~180，
     且换 ngrid 还要再换）。这里改成显式注入到示踪物上，强度由参数直接控制：

     目标：把示踪物过密度场乘上势的调制因子（局部、乘性响应）
         δ_t(x) -> δ_t(x) * [1 + Ainj * φ_G(x)],   Ainj = 2 * fnl * b_phi
     实现成"把所有示踪物平移一个位移场 Ψ"（乘性 n(1+Aφ) 与平移在 A 的一阶严格
     等价：δ' = δ + Aφ + Aφδ + O(Ψ∇δ 的圈图项)，圈图项相对主项 ~1/(b1*M)~1e-3，
     是真实的非局部修正，本就超出本次标定精度）：
         n'(x) = n(x - Ψ(x))  =>  δ̂'(k) = δ̂(k) - i k·Ψ̂(k) + O(Ψ·∇n)
     要求 -i k·Ψ̂(k) = Ainj*φ̂_G(k)，即
         Ψ̂_i(k) = i * k_i * (Ainj/k²) * φ̂_G(k)     （纯 1/k 核，无转移函数）

     与目标响应的对应：注入给 P(k) 的响应为
         L_inj(k) = [P(+A)-P(-A)]/(2P_G) = 2*Ainj/(b1*M(k)) = 4*fnl*b_phi/(b1*M(k)),
     与 Quijote 的 L_Q = 4*f_NL*b_φ/(b1*M) 同形，取 b_phi = b_φ(Q) 即 1:1（b1 由
     高斯标定本来就对上）。

     代码里用与 ZA 位移完全相同的 i·k 结构：mesh->phik 存的就是高斯势的傅里叶
     系数（c2r 出去即实空间 φ_G），rhok2/3/4 在 psi 算完后已空闲，拿来当复用缓冲。
     注意 k=0 模式 phik=0，但 1/k² 会产生 0*∞=NaN，必须显式跳过。 */
  if (mesh->inj[0] && mesh->inj[1] && mesh->inj[2]) {
    const double Ainj = 2.0 * conf->fnl * conf->b_phi;
    printf("PNG tracer-level injection: Ainj = 2*FNL*B_PHI = %.6g\n", Ainj);
    for (int i = 0; i < conf->Ng; i++) {
      size_t idx_i = (size_t) i * conf->Ng;
      double ki = (i < Ngh) ? i : i - conf->Ng;
      for (int j = 0; j < conf->Ng; j++) {
        size_t idx_ij = idx_i + j;
        double kj = (j < Ngh) ? j : j - conf->Ng;
        double kij = ki * ki + kj * kj;
        size_t idx0 = idx_ij * Ngk;
        for (int k = 0; k < Ngk; k++) {
          size_t idx = idx0 + k;
          if (i == 0 && j == 0 && k == 0) {   /* k=0：φ̂=0 且 1/k² 发散，直接置零 */
            mesh->rhok2[idx][0] = mesh->rhok2[idx][1] = 0;
            mesh->rhok3[idx][0] = mesh->rhok3[idx][1] = 0;
            mesh->rhok4[idx][0] = mesh->rhok4[idx][1] = 0;
            continue;
          }
          double ksq = kij + (double) k * k;
          /* kern = Ainj/k_phys²，再乘 (ki·kfac) 等得到 Ψ̂_i(k) 的系数 */
          double kern = Ainj / (kfac * kfac * ksq);
          /* i*φ̂_G：φ̂_G=(re,im)=mesh->phik -> i*φ̂_G=(-im, re) */
          double re = -mesh->phik[idx][1];
          double im =  mesh->phik[idx][0];
          mesh->rhok2[idx][0] = re * ki * kfac * kern;
          mesh->rhok2[idx][1] = im * ki * kfac * kern;
          mesh->rhok3[idx][0] = re * kj * kfac * kern;
          mesh->rhok3[idx][1] = im * kj * kfac * kern;
          mesh->rhok4[idx][0] = re * kfac * k * kern;
          mesh->rhok4[idx][1] = im * kfac * k * kern;
        }
      }
    }
    FFT_EXEC_C2R(plan, mesh->rhok2, mesh->inj[0]);
    FFT_EXEC_C2R(plan, mesh->rhok3, mesh->inj[1]);
    FFT_EXEC_C2R(plan, mesh->rhok4, mesh->inj[2]);
  }






}


#undef EZMOCK_LOGPK_NAME
#undef EZMOCK_FIXAMP_NAME
#undef EZMOCK_IPHASE_NAME

#undef EZMOCK_PT_LOGPK
#undef EZMOCK_PT_FIXAMP
#undef EZMOCK_PT_IPHASE


/******************************************************************************/

/* Macros for the template functions. */
#elif defined(EZMOCK_PT_WHITENOISE) && defined(EZMOCK_PT_LOGPK) && \
  !defined(EZMOCK_PT_FIXAMP) && !defined(EZMOCK_PT_IPHASE)

/* Macros for generating function names. */
#ifndef CONCAT_FNAME2
  #define CONCAT_FNAME2(a,b)            a##b
#endif

#ifndef EZMOCK_PT_FUNCNAME2
  #define EZMOCK_PT_FUNCNAME2(a,b)      CONCAT_FNAME2(a,b)
#endif

/*============================================================================*\
                             Definition validation
\*============================================================================*/

#ifdef EZMOCK_LOGPK_NAME
  #undef EZMOCK_LOGPK_NAME
#endif

#if     EZMOCK_PT_LOGPK == 1
  #define EZMOCK_LOGPK_NAME     _logpk
#elif   EZMOCK_PT_LOGPK == 0
  #define EZMOCK_LOGPK_NAME
#else
  #error "unexpected definition of `EZMOCK_PT_LOGPK`"
#endif


/*============================================================================*\
                   Function for displacement field generation
\*============================================================================*/

/******************************************************************************
Function `EZmock_ZA_disp_wn<EZMOCK_LOGPK_NAME>`:
  Generate the Zel'dovich displacement field given the white noise field and
  the input power spectrum.
Arguments:
  * `ez`:       instance of the EZmock generator;
  * `plan`:     FFTW plan;
  * `err`:      integer storing the error code.
******************************************************************************/
static void EZMOCK_PT_FUNCNAME2(EZmock_ZA_disp_wn, EZMOCK_LOGPK_NAME)
    (EZMOCK *ez, FFT_PLAN plan, int *err) {
  EZMOCK_CONF *conf = (EZMOCK_CONF *) ez->conf;
  EZMOCK_PK *pk = (EZMOCK_PK *) ez->pk;
  EZMOCK_MESH *mesh = (EZMOCK_MESH *) ez->mesh;

  const int Ngk = (conf->Ng >> 1) + 1;
#ifdef OMP
  const size_t pnum = ((size_t) conf->Ng * conf->Ng) / conf->nthread;
  const int rem = ((size_t) conf->Ng * conf->Ng) % conf->nthread;
#endif

  const double kfac = M_PI * 2 / conf->Lbox;
#if EZMOCK_PT_LOGPK == 1
  const double logkfac = log(kfac);
#endif
  const double Pfac = pow(conf->Lbox, -1.5);

  /* First dimension: compute psi[0]. */
#ifdef OMP
#pragma omp parallel num_threads(conf->nthread)
  {
    /* Distribute mesh grids to threads. */
    const int tid = omp_get_thread_num();
    const size_t pcnt = (tid < rem) ? pnum + 1 : pnum;
    const size_t istart = (tid < rem) ? pcnt * tid : pnum * tid + rem;
    const size_t iend = istart + pcnt;

    /* Traverse the Fourier space density field with OpenMP. */
    for (size_t idx_ij = istart; idx_ij < iend; idx_ij++) {
      int i = idx_ij / conf->Ng;
      int j = idx_ij % conf->Ng;
      double ki = (i <= (conf->Ng >> 1)) ? i : i - conf->Ng;
#else
  /* Traverse the Fourier space density field sequentially. */
  for (int i = 0; i < conf->Ng; i++) {
    size_t idx_i = (size_t) i * conf->Ng;
    double ki = (i <= (conf->Ng >> 1)) ? i : i - conf->Ng;
    for (int j = 0; j < conf->Ng; j++) {
      size_t idx_ij = idx_i + j;
#endif
      double kj = (j <= (conf->Ng >> 1)) ? j : j - conf->Ng;
      double kij = ki * ki + kj * kj;
      size_t idx0 = idx_ij * Ngk;

      /* Expand k loop and treat 0 frequency separately. */
      /* k = 0 */
      size_t idx = idx0;
      if (idx == 0) {
        mesh->rhok[0][0] = mesh->rhok[0][1] = 0;
        mesh->rhok2[0][0] = mesh->rhok2[0][1] = 0;
      }
      else {
        double ksq = kij;
#if EZMOCK_PT_LOGPK == 1
        double kmod = log(ksq) * 0.5 + logkfac;         /* log(sqrt(ksq)) */
#else
        double kmod = kfac * sqrt(ksq);                 /* sqrt(ksq) */
#endif
        /* Obtain the power with interpolation. */
        double P = pk_interp(pk, kmod);

#if EZMOCK_PT_LOGPK == 1
        P = exp(P * 0.5);
#else
        P = sqrt(P);
#endif
        P *= Pfac / (ksq * kfac);
        mesh->rhok[idx][0] = -P * mesh->rhok2[idx][1];
        mesh->rhok[idx][1] = P * mesh->rhok2[idx][0];
        mesh->rhok2[idx][0] = mesh->rhok[idx][0] * ki;
        mesh->rhok2[idx][1] = mesh->rhok[idx][1] * ki;
      }

      /* 0 < k <= k_ny */
      for (int k = 1; k < Ngk; k++) {
        idx = idx0 + k;
        double ksq = kij + k * k;
#if EZMOCK_PT_LOGPK == 1
        double kmod = log(ksq) * 0.5 + logkfac;         /* log(sqrt(ksq)) */
#else
        double kmod = kfac * sqrt(ksq);                 /* sqrt(ksq) */
#endif
        /* Obtain the power with interpolation. */
        double P = pk_interp(pk, kmod);

#if EZMOCK_PT_LOGPK == 1
        P = exp(P * 0.5);
#else
        P = sqrt(P);
#endif
        P *= Pfac / (ksq * kfac);
        mesh->rhok[idx][0] = -P * mesh->rhok2[idx][1];
        mesh->rhok[idx][1] = P * mesh->rhok2[idx][0];
        mesh->rhok2[idx][0] = mesh->rhok[idx][0] * ki;
        mesh->rhok2[idx][1] = mesh->rhok[idx][1] * ki;
      }         /* for k */
    }           /* for idx_ij (OMP) or for j (no OMP) */
  }             /* omp parallel (OMP) or for i (no OMP) */

  FFT_EXEC_C2R(plan, mesh->rhok2, mesh->psi[0]);

  /* Second and third dimension: compute psi[1] and psi[2]. */
#ifdef OMP
#pragma omp parallel num_threads(conf->nthread)
  {
    /* Distribute mesh grids to threads. */
    const int tid = omp_get_thread_num();
    const size_t pcnt = (tid < rem) ? pnum + 1 : pnum;
    const size_t istart = (tid < rem) ? pcnt * tid : pnum * tid + rem;
    const size_t iend = istart + pcnt;

    /* Traverse the Fourier space density field with OpenMP. */
    for (size_t idx_ij = istart; idx_ij < iend; idx_ij++) {
      int j = idx_ij % conf->Ng;
#else
  /* Traverse the Fourier space density field sequentially. */
  for (int i = 0; i < conf->Ng; i++) {
    size_t idx_i = (size_t) i * conf->Ng;
    for (int j = 0; j < conf->Ng; j++) {
      size_t idx_ij = idx_i + j;
#endif
      double kj = (j <= (conf->Ng >> 1)) ? j : j - conf->Ng;
      size_t idx0 = idx_ij * Ngk;

      /* 0 <= k <= k_ny */
      for (int k = 0; k < Ngk; k++) {
        size_t idx = idx0 + k;
        mesh->rhok2[idx][0] = mesh->rhok[idx][0] * kj;
        mesh->rhok2[idx][1] = mesh->rhok[idx][1] * kj;
        mesh->rhok[idx][0] *= k;
        mesh->rhok[idx][1] *= k;
      }
    }           /* for idx_ij (OMP) or for j (no OMP) */
  }             /* omp parallel (OMP) or for i (no OMP) */

  FFT_EXEC_C2R(plan, mesh->rhok2, mesh->psi[1]);
  FFT_EXEC_C2R(plan, mesh->rhok, mesh->psi[2]);

  /* Renormalise the displacement fields. */
  const size_t num = (size_t) conf->Ng * conf->Ng * conf->Ng;
  const real norm = 1 / sqrt(num);
  for (int j = 0; j < 3; j++) {
#ifdef OMP
#pragma omp parallel for num_threads(conf->nthread)
#endif
    for (size_t i = 0; i < num; i++) mesh->psi[j][i] *= norm;
  }
}


#undef EZMOCK_LOGPK_NAME

#undef EZMOCK_PT_WHITENOISE
#undef EZMOCK_PT_LOGPK

#endif
