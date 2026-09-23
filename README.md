# ezmock_png — 带 FNL 参数的改造版 EZmock

> **GitHub 源码快照（2026-09-23）**：本仓库取自 NERSC
> `/pscratch/sd/l/lzy/PNG-EZmock/codes/ezmock_png/`，包含源码、构建文件、
> `Tk_0.txt` 和原理 notebook；不包含已编译的 `EZmock`、运行输出或
> Quijote 模拟数据。下文引用的 `../ezmock_png_binary_verification/` 等目录
> 属于 NERSC 工作项目，未随本源码仓库上传。`options.mk` 中的 FFTW 路径是
> NERSC 环境路径；在其他机器编译时可用 `make -C src FFTW_DIR=/你的/FFTW/路径`
> 覆盖。
>
> 本代码基于 [Cheng Zhao 的 EZmock](https://github.com/cheng-zhao/EZmock)。
> 上游整体以 GPLv3 发布，部分源文件标注 MIT；本仓库保留了上游的
> `LICENSE.txt`、`LICENSE_MIT.txt` 和 `libEZmock/prand/LICENSE.txt`。

**这是什么**：主人自己改造的 EZmock（上游原版：Cheng Zhao, https://github.com/cheng-zhao/EZmock）。
相比 stock 版新增 primordial non-Gaussianity 注入，共**两种机制**：

1. **场层面（历史机制）**：位移前在实空间对势场做 `phi_png = phi + FNL * phi^2`；
2. **示踪物层面（2026-09-21 新增）**：对示踪物做位移 `Psi = i·k·(Ainj/k²)·phi_G`，
   即密度调制 `delta_t -> delta_t·[1 + Ainj·phi_G]`，其中 **`Ainj = 2·FNL·B_PHI`**
   （详见下面"示踪物层面 PNG 注入"一节）。

fNL 通过配置关键字 `FNL`（或库接口 `EZmock_set_fnl()`）自由设置，0 表示纯高斯初条件；
`B_PHI` 控制示踪物层面注入强度（0 = 不注入，逐位回退到历史行为）；
`FNL_FIELD` 控制场层面二次项系数（不设时自动推导，见下）。

**来源**：`/global/u2/l/lzy/pyEZmocktest/softdir/EZmock` 的 rsync 副本（2026-09-21，排除 *.o/*.so/example/），
外加 Makefile、options.mk、Tk_0.txt。原树已复原为改造前状态（改动撤回、二进制换回旧版）；
**本目录收录 FNL 改造版的源码**。

**详细原理文档（带公式 + 可复跑数据自检）**：`png_injection_mechanisms.ipynb`
（生成脚本 `tmp_doc_build/build_png_injection_doc.py`，改脚本后重跑即可再生成）。
两种机制的推导、实现表、验证数据一图流都在这份 notebook 里。

## 与 stock 版的关键差异（最容易踩坑的地方）

1. **LINEAR_PK 必须是转移函数 T(k)**（量级 ~O(1)，如自带的 `Tk_0.txt`），**不是 P(k) 表**。
   位移公式 `twb = Transfer/Dplus/Beta` 直接乘输入曲线；喂 P(k)（~1e4 量级）会让位移放大上万倍
   → CIC 装填越界段错误。
2. **PK_INTERP_LOG 必须为 F（或不设，默认 F）**。log 插值模式在本改造版中是坏的：
   perturb.c 里对 k 先 `log(k)` 后又对结果 `log()` → 全场 NaN → 段错误。
3. **初条件形状仍硬编码，生长因子已动态化**（2026-09-21 改造，详见下节）：
   - 写死的只剩初条件形状 `Anorm=12275.369233`、`PrimordialIndex=0.9624`（与 `Tk_0.txt` 配套；
     输入表只在位移步进 `pk_interp()` 里当 `Transfer` 用）。原代码里**未被使用**的
     `k_pivot=0.05` 和 `h=0.6711` 已删除。
   - `D0 / Dplus / Beta` 不再写死，由配置的 `OMEGA_M / DE_EOS_W / REDSHIFT` 动态算出。
4. **双精度构建**（stock 版为单精度），依赖 kirisame conda env 的 libfftw3；因此两版输出
   无法逐位对比。
5. 输出按 rank/PDF 装填：小场变化会重排示踪点位置，**比较不同 FNL 必须用统计量（P(k) 等）**，
   不能逐点比较。

## 生长因子动态化（2026-09-21）

**动机**：原先 `perturb.c` 里写死 `D0=0.789246`、`Dplus=1.650`、`Omega_m=0.315`（连带 `Beta`），
只用配置的 `OMEGA_M` 算 `vfac`。后果是：改 `OMEGA_M/REDSHIFT` 不会改变密度场幅度，
且**速度永远按 z=0 的 `vfac` 生成**——外部管线再按目标红移做 RSD 换算时，RSD 强度会错。

**约定**（沿用 2LPT-PNG，`D_A(a)` = 未归一化增长因子，早期 `D_A→a`）：

| 量 | 公式 | 说明 |
|---|---|---|
| `D0` | `D_A(1) = g(z=0)` | PNG 势能归一化 `Beta` 的锚点，与输出红移无关 |
| `Dplus` | `D_A(1)/D_A(a_out) = D(0)/D(z_out)` | 位移场整体乘 `1/Dplus`，把 z=0 形状的 T(k) 场缩放到输出红移 |
| `Beta` | `1.5·Omega_m/(2998²·D0)` | h 在 k 用 h/Mpc、H0 用 100 km/s/Mpc 时自动抵消 |
| `vfac` | `f(z)·H(z)·a` | 现在在**输出红移**处求值（原先固定 z=0） |

**语义变化（重要）**：`REDSHIFT` 从"只通过 growth2 缩放输入表的摆设"变成**盒子的真实输出红移**，
同时驱动 `Dplus`（密度幅度）与 `vfac`（RSD 速度）。`REDSHIFT_PK` 对本改造版无意义（被忽略并告警）。

**另一个坑：growth2 双重计入**。`linear_pk.c` 会按 stock 版约定把输入表乘
`growth2 = (D(z)/D(z_pk))²`；本版输入表是 z=0 的 T(k)，红移换算交给 `Dplus`，
两者叠加会让位移幅度再低 `(D(1)/D(0))² ≈ 0.37` 倍。因此 `src/main/run_mock.c` 在
`EZmock_set_cosmology()` 之后**强制 `cosmo->growth2 = 1.0`**，并打印生长因子溯源行：

```
PNG growth factors: Omega_m = 0.3175, D0 = 0.789246093, Dplus = 1.650333193, Beta = 6.713659e-08, vfac = 78.796529
```

配置里显式给了 `GROWTH_PK/VELOCITY_FAC` 时（`eval_growth=false` 分支）：`VELOCITY_FAC` 照用，
`GROWTH_PK` 仍被置 1；`D0/Dplus/Beta` 依旧按 `OMEGA_M/DE_EOS_W/REDSHIFT` 算，所以那几个键不能省。

**迁移对照**（`twb ∝ D0/Dplus`，位移幅度）：

| 配置 | D0 | Dplus | vfac | 与旧写死口径的位移幅度比 |
|---|---|---|---|---|
| 旧代码（任何 OMEGA_M） | 0.789246 | 1.650 | z=0 | 1.000（基准） |
| 新代码 `OMEGA_M=0.3175, REDSHIFT=1` | 0.789246093 | 1.650333193 | 78.7965 | 0.9998 |
| 新代码 `OMEGA_M=0.3089, REDSHIFT=1` | 0.784269602 | 1.642643288 | 77.6892 | 0.9981 |
| 新代码 `OMEGA_M=0.3089, REDSHIFT=0` | 0.784269602 | 1.000000000 | 52.1324 | **1.640** |

也就是说：**老配置（`REDSHIFT=0`）在新代码下位移幅度会变成 1.64 倍**，实空间 P0 在拟合区高出
约 50%（A/B 实测 1.41–1.75，k 依赖，见验证目录）；想复现旧行为请把 `REDSHIFT` 设成目标红移 1。

**改动文件**：`perturb.c`（改用 `cosmo->png_d0/png_dplus/omega_m`；`Beta` 公式化简，与原
`UnitLength_in_cm` 形式等价）、`cosmology.c`（新增计算 + 修 `cosmo_growth_D_only` 的 z_a/z_ainit
混用 bug）、`structs.h`/`EZmock.c`/`EZmock.h`（新增 `omega_m/png_d0/png_dplus` 字段与函数声明）、
`src/main/run_mock.c`（中和 growth2 + 打印）。数值交叉验证见
`../ezmock_png_binary_verification/README.md`。

## 示踪物层面 PNG 注入（B_PHI + FNL_FIELD，2026-09-21）

**动机**：场层面机制 `phi_png = phi + FNL·phi²` 注入的是**初条件**非高斯，经过 EZmock 的
PDF/位移/示踪物采样后响应会被压缩——实测固定四旋钮时，想在大尺度对上 Quijote f_NL=100
的 PNG 信号需要把 EZmock 的 FNL 开到 ~180，且这个倍率依赖 ngrid 与目标样本（不是常数修正）。
**示踪物层面注入**绕过整条标定链，直接给最终示踪物叠加上目标密度的 PNG 响应，
使"配置 FNL 是多少就匹配 f_NL 是多少的 N-body"。

**推导（与代码注释一致）**：希望密度场做如下的线性（对 Ainj）调制

    delta_t(x)  ->  delta_t(x) · [1 + Ainj · phi_G(x)],    Ainj = 2 · FNL · B_PHI

把示踪点整体位移 `n'(x) = n(x - Psi(x))`，一阶展开给出

    delta_hat'(k) = delta_hat(k) - i k·Psi_hat(k) + O(Psi·grad n)

要求 `-i k·Psi_hat(k) = Ainj·phi_hat_G(k)`，即

    Psi_hat_i(k) = i · k_i · (Ainj / k²) · phi_hat_G(k)

——纯 1/k² 核（**不乘任何转移函数**，与 ZA 位移的 `twb` 形状项不同）。
`O(Psi·grad n)` 的平流/圈图项相对大小 ~1/(b₁·M)~1e-3，可忽略。
目标响应（对照 Quijote f_NL=100 的 500 realizations P0）：

    L(k) = 2·Ainj/(b₁·M(k)) = 4·FNL·B_PHI/(b₁·M(k))，  M = k²T/(Dplus·Beta)

因此只要取 **`B_PHI = b_phi(N-body 样本)`** 就实现 1:1；`b_phi` 是示踪物对 `phi_G` 的偏置
（Quijote 测量值在 k≈0.008 处约 2.4，拟合区中位 ~1.45，见验证目录）。

**实现**：`perturb.c` 里沿用 `mesh->phik`（phi_G 的 FFT 系数，经 `phik_copy` 在 `plan1`
之后仍然存活），对每个模构造 `(i·phi_hat_G)·k_i·(Ainj/k²)`，做三次 C2R 得实空间位移场
`mesh->inj[0..2]`；`pop_tracer.c` 在生成示踪点后用三线性插值把 `inj` 场采到每个示踪点上、
加位移并 wrap 回盒子。`FNL=0` 或 `B_PHI=0` 时完全不分配 `inj` 场、不走新分支。

**配置语义**（`load_conf.c`）：

| 键 | 含义 | 不设时 |
|---|---|---|
| `B_PHI` | 示踪物层面响应 `Ainj = 2·FNL·B_PHI` | 0（不注入） |
| `FNL_FIELD` | 场层面二次项系数 `phi + FNL_FIELD·phi²` | `B_PHI!=0` → 0（标定模式，关掉未标定的场层面项）；否则 = `FNL`（历史行为） |

- `FNL_FIELD=0` 表示场层面二次项恰好关闭（`phi_png = phi`）；库层面用 `HUGE_VAL`
  哨兵值表示"未设置"（回退到 `FNL`），`EZmock_set_fnl_field()` 未调用时行为同历史版。
- **逐位一致保证**：`B_PHI=0` 且 `FNL_FIELD` 不设时，配置路径、RNG 流、输出与改造前二进制
  逐字节相同（A/B 实测,FNL=180 与 FNL=0 两档均过,见验证目录）。
- 输出头记录两键的值（ASCII：`# B_PHI=... , FNL_FIELD=...`；FITS 键 `B_PHI`/`FNL_FLD`）。

## 改造点文件表

FNL 接入（9 个文件）：

| 文件 | 改动 |
|---|---|
| `libEZmock/EZmock/structs.h` | EZMOCK_CONF 增加 `double fnl;` |
| `libEZmock/EZmock/EZmock.c` | init 设 `conf->fnl = 0.0`；新增 `EZmock_set_fnl()` |
| `libEZmock/EZmock/EZmock.h` | 声明 `EZmock_set_fnl()` |
| `libEZmock/EZmock/perturb.c` | `double Fnl=180;` → `const double Fnl = conf->fnl;`；phi_png 二次项 |
| `src/main/define.h` | `#define DEFAULT_FNL 0.0` |
| `src/main/load_conf.h` | CONF 结构体加 `double fnl;` |
| `src/main/load_conf.c` | 命令行 `-f/--fnl`、配置模板、参数表键 `FNL`、默认值校验、打印 |
| `src/main/run_mock.c` | 读配置后调用 `EZmock_set_fnl(ez, conf->fnl, &err)` |
| `src/main/save_res.c` | ASCII/FITS 输出头记录 FNL 值 |

示踪物层面注入（在 FNL 基础上追加改动）：

| 文件 | 改动 |
|---|---|
| `libEZmock/EZmock/structs.h` | EZMOCK_CONF 加 `b_phi`/`fnl_field`；EZMOCK_MESH 加 `inj[3]` 位移场 |
| `libEZmock/EZmock/EZmock.c` | init 初始化两字段与 `inj` 指针；新增 `EZmock_set_b_phi()` / `EZmock_set_fnl_field()` |
| `libEZmock/EZmock/EZmock.h` | 声明两个 setter |
| `libEZmock/EZmock/structs.c` | `EZmock_mesh_destroy()` 释放 `inj[3]` |
| `libEZmock/EZmock/dens_field.c` | 需要注入时分配 `inj[3]`（`fnl!=0 && b_phi!=0`） |
| `libEZmock/EZmock/perturb.c` | `Fnl_field` 逻辑 + 1/k² 位移场构造（含推导注释） |
| `libEZmock/EZmock/pop_tracer.c` | `apply_png_shift()`：三线性插值 + wrap（OMP 并行） |
| `src/main/define.h` | `#define DEFAULT_B_PHI 0.0` |
| `src/main/load_conf.h` | CONF 加 `double b_phi; double fnl_field;` |
| `src/main/load_conf.c` | 命令行 `--b-phi`/`--fnl-field`、模板、键 `B_PHI`/`FNL_FIELD`、FNL_FIELD 默认推导、校验、打印 |
| `src/main/run_mock.c` | 调用两个新 setter |
| `src/main/save_res.c` | ASCII/FITS 输出头记录 B_PHI/FNL_FIELD |

## 构建

```bash
make -C src      # 注意：顶层 `make EZmock` 因目标名与二进制同名会被跳过，必须 -C src
```
FFTW 路径在 `options.mk`（`FFTW_DIR` 指向 kirisame env，双精度 `-lfftw3 -lfftw3_omp`）。

运行需要找到 libfftw3_omp（系统没有）：
```bash
export LD_LIBRARY_PATH=/global/homes/l/lzy/anaconda3/envs/kirisame/lib
```
注意该 LD_LIBRARY_PATH 会同时遮蔽系统 libgomp/libncursesw（曾导致 gdb 自身无法启动）；
调试别的东西时不要全局 export，用 `gdb -batch -ex "set env LD_LIBRARY_PATH ..."` 只给被调试进程设。

当前二进制 md5：`EZmock` = **2827e1af53df31df8baac982f4ad0577**（2026-09-21，示踪物层面注入改造后编译）。
历史构建物备份：生长因子动态化后的 `764e0ab80de784ca895d4e06e45a7bfe`
（`../ezmock_png_binary_verification/tmp_growth_dynamic/pre_change_backup/EZmock_pre20260921_build` 是更早的
`38e8bd21df3b5277dbeae90693057ac1`）；注入改造前的二进制备份在
`../ezmock_png_binary_verification/tmp_png_response/pre_injection_backup/EZmock_pre20260921_inj`
（md5 `764e0ab80de784ca895d4e06e45a7bfe`，B_PHI=0 时与新二进制逐字节一致）。

## 用法（历史模板：主人自己那份 ezmock.conf）

> 这一段只是**主人早年那份 ezmock.conf 的复刻**（Tk_0.txt + 线性插值 + Ωm=0.3089 + FNL=180），
> 留作参考；notebook 的 `EZMOCK_MODE="modified"` 会自己生成配置，**不需要用到这个文件**。
> 唯一要记住的是：`REDSHIFT` 必须是盒子的真实输出红移（迁移对照见上一节）。

```
BOX_SIZE = 1000
NUM_GRID = 256
NUM_TRACER = 163629
LINEAR_PK = <本目录>/Tk_0.txt      # 转移函数，不是 P(k)！
PK_INTERP_LOG = F                   # 或不设
REDSHIFT_PK = 0                     # 被忽略（输入表必须是 z=0 的 T(k)）
RAND_SEED = 1995
FIX_AMPLITUDE = True
OMEGA_M = 0.3089
REDSHIFT = 1                        # 盒子的真实输出红移！改造前写死生长因子、此处填 0，
                                    # 现在填 0 会把位移幅度放大 1.64 倍（见上一节迁移对照）
RHO_CRITICAL = 1.64
RHO_EXP = 5
PDF_BASE = 0.17
SIGMA_VELOCITY = 275
ATTACH_PARTICLE = False
FNL = 180                           # 新增：fNL 参数
```

示踪物层面注入的配置示例（1:1 匹配 Quijote f_NL=100 的用法）：

```
FNL = 100                           # 物理 f_NL，直接填 N-body 的值
B_PHI = 2.65                        # 示踪物响应；k<0.01 的 1:1 标定值（见验证目录；
                                    # 换目标样本/重调四旋钮后需重标）
FNL_FIELD = 0                       # 关掉场层面二次项（B_PHI 设了之后这是默认）
```

只需要"关掉旧机制、什么都不注入"时：`FNL = 0`（或 `B_PHI = 0`）即可，输出与改造前逐位一致。

## 验证记录

### 生长因子动态化（2026-09-21）

- **数值交叉验证**：Ωm=0.3175 时新代码给出 `D0=0.789246093`、`Dplus(z=1)=1.650333193`，
  与独立 ODE 积分一致到 1e-5、与含辐射的 CAMB 差 1.4e-4；旧写死值正是它们（0.789246 / 1.650）。
- **统计连续性 A/B/C**（各 3 个 fixed-amplitude realization，同种子；A=旧二进制 z=0、
  B=新二进制 z=0、C=新二进制 z=1/Ωm=0.3175）：C 的实空间 P0 与 A 差 **1.58%**（含 realization
  散射），而 B 高出 **49%**——证明按 `REDSHIFT=1` 迁移能保住原幅度，`REDSHIFT=0` 则不行。
- **RSD 修正的物理后果**：旧口径下 python 侧 RSD 系数只有正确值的 0.671 倍（速度在 z=0 生成），
  四极矩 P2/P0 = 0.2276 vs Quijote 0.3608；新口径（z=1）为 0.3185（比值 0.873，仍偏低 ~13%，
  属四旋钮可调范围）。
- **逐胞元计数检验不可用**（负结果）：改造版的示踪点采样本身是随机的（每格先抽高斯再定标记数），
  密度场差 0.02% 就会让 RNG 流错位，占据格集合的 Jaccard 从"同配置逐位相同"的 1.0 掉到 0.013
  （与密度场差 65% 的 0.011 几乎无差别）。**连续性只能在 P(k) 统计层面验证**。
  详见 `../ezmock_png_binary_verification/README.md`。

### 示踪物层面注入（2026-09-21）

- **A/B 逐字节回归**：`B_PHI=0`（历史模式，含 FNL=180 场层面注入开启与 FNL=0 纯高斯两档）
  下，新二进制与注入改造前二进制输出 md5 完全相同（`c6edcf500125b96b57308a55d7178a81` /
  `67a8b7c1c6974f925d0720310764a332`）；显式写 `B_PHI=2.9, FNL_FIELD=0` + `FNL=0` 的
  数据部分也与纯高斯 FNL=0 一致（Ainj=2·0·2.9=0，符合预期）。
- **冒烟**：`FNL=100, B_PHI=2.9, FNL_FIELD=0` 单跑 3.0 s 通过，输出与 FNL=0 明显不同，
  文件头正确记录 `# B_PHI=2.9 , FNL_FIELD=0`。脚本与配置见
  `../ezmock_png_binary_verification/tmp_png_response/inj_smoke/`。
- **1:1 匹配验证**（2026-09-21，`B_PHI=2.9`、`FNL=±100`、8 realizations、`analyze_injection.py`）：
  - **响应 = 常数 × 2Ainj/(b₁·M(k))**：把测量写成 `L_inj·b₁·M/(2Ainj)`，在 k∈[0.008,0.045]
    上是**平的**（中位 0.893），与推导公式完全一致；
  - **ngrid 无关**：同一常数在 Ng=128 为 0.921、Ng=256 为 0.893（差 3%；**注意口径**：
    这两个数都要用**各自网格自己的 EZ 示踪物 b₁** 归一化——两网格的示踪物样本不同，
    b₁(128)/b₁(256)=0.82，用 Quijote 代理 b₁ 会假性放大成 17%）。对照：旧场层面机制
    L(128)/L(256)=0.233（4.3 倍）；即使扣掉 b₁ 差异仍是 0.201（5.0 倍）——真正的网格依赖，
    与示踪物 b₁ 无关（详见 `png_injection_mechanisms.ipynb` §5.5 与自检单元①）；
  - **标定值（定稿口径：主人关心的超大尺度 k<0.01）**：**B_PHI = 2.65**——两轮独立互证：
    v2 轮（B_PHI=2.9）bin1 需 2.647±0.100；v3 轮（B_PHI=2.6）实测 L_inj/L_Q=0.982±0.037，
    响应精确线性 ⇒ 外推 2.65。整体区域（若关心到 k~0.04）1:1 需要 2.57 (k∈[0.006,0.012))
    → 1.82 (k∈[0.03,0.045))（Ng=256 现旋钮下；标定值 ∝ 示踪物 b₁，重调四旋钮后需重标）；
  - **残留形状偏差**：目标响应 `b_φ^Q(k)` 从 2.41 平滑降到 1.50（40%），而本机制的响应
    严格是 1/(b₁M) 形状 ⇒ **常数 B_PHI 只在单一尺度上 1:1**；该漂移是目标样本自身的性质
    （旧机制同方向、受限于噪声看不出）。若要全局 1:1，需要给注入核乘一个响应形状表
    W(k)（k 表文件，仿照 twb），换样本时需重新标定。

### FNL 接入（2026-09-21 早先一轮）

- `FNL=180`（本版） vs 旧二进制硬编码 180：输出**逐字节一致**
  （md5 `e948d1dadf7811cf56c7d96aa6ae4ad6`）。
- `FNL=0` vs `FNL=180`：输出明显不同（位置差中位 ~156 Mpc/h），符合预期。
- 纯仓库路径端到端 smoke：2.4 s 跑通，与上述验证输出逐字节一致（自包含，只依赖
  Tk_0.txt + kirisame libfftw3）。
- 详细的对比日志与配置见 `../ezmock_png_binary_verification/`。
- `../ezmock_fnl0_calibration/manual_tune.ipynb` 已接入：`EZMOCK_MODE="modified"` 时用本目录的
  EZmock + Tk_0.txt + 线性插值，面板参数 `FNL` 直接写进配置；smoke 结果见验证目录 README。
  生长因子动态化后该模式的面板已同步为 `OMEGA_M=0.3175`、`REDSHIFT=TARGET_REDSHIFT(=1)`，
  并加了"REDSHIFT 必须等于目标红移"的断言。
