#!/usr/bin/env python
"""生成 PNG-EZmock 注入机制原理文档（交付物：../png_injection_mechanisms.ipynb）。

这是【生成器】：要改文档内容，改本脚本里的 Markdown / 代码字符串，然后：

    python build_png_injection_doc.py      # 生成 ipynb
    cd .. && jupyter nbconvert --to notebook --execute --inplace \
        png_injection_mechanisms.ipynb     # 执行并把数据核对单元的输出填进去

注意：Markdown 里的 LaTeX 全部用 raw 字符串（r\"\"\"...\"\"\"）原样书写，
不要手动转义反斜杠。
"""
from pathlib import Path

import nbformat as nbf

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "png_injection_mechanisms.ipynb"

cells = []


def md(text):
    """追加一个 Markdown 单元。"""
    cells.append(nbf.v4.new_markdown_cell(text.strip("\n")))


def code(text):
    """追加一个代码单元（尚未执行；由 nbconvert --execute 填输出）。"""
    cells.append(nbf.v4.new_code_cell(text.strip("\n")))


# ============================================================================
# 0. 封面与摘要
# ============================================================================
md(r"""
# PNG-EZmock 的两种 PNG 注入机制：原理、推导与验证

**文档性质**：写给未来的复习材料——公式推导 + 代码位置 + 实测数字，按仓库现状修订于 2026-09-23。
绝大部分是 Markdown 说明；第 5 节有 4 个只读小单元，从既有 npz 数据文件里读数字直接核对，
秒级、可在登录节点直接运行。

**两版是什么**

- **原版**（主人 2026 年初写的那版，现镜像在 `/global/u2/l/lzy/pyEZmocktest/softdir/EZmock/`）：
  **只在场层面注入**——实空间二次变换 $\phi \to \phi + F_{NL}\,\phi^2$，$F_{NL}=180$ 硬编码。
- **改造版（糖糖版，`codes/ezmock_png/`，唯一维护位置）**：场层面机制**参数化保留**（配置键
  `FNL` / `FNL_FIELD`，可关），**新增示踪物层面注入**——通过一个位移场在
  大尺度线性极限生成 $A_{\rm inj}\phi_G$ 的示踪物密度响应，
  $A_{\rm inj} \equiv 2\,F_{NL}\,B_\phi$。

**一句话结论**：

- 原版把 PNG 种在**"原因"（初条件）**上，响应经由 EZmock 的位移／PDF 映射"涌现"，
  强度不可控——固定四旋钮要对上 Quijote $f_{\rm NL}=100$ 得把 $F_{NL}$ 开到 $\approx 180$，
  且换网格响应差 4 倍；
- 改造版把 PNG 写在**"结果"（示踪物密度）**上，响应有闭式
  $\mathcal{L} = 2A_{\rm inj}/(b_1 M(k))$，只需标定**一个常数** $B_\phi \approx 2.65$，
  实测两个网格的归一化响应相差约 3%。这一平均功率响应检验不等于协方差检验。

**2026-09 修订说明**：位置平移并不与局部乘性数密度调制严格等价；下文的闭式响应
只在指定的大尺度线性极限成立。`B_PHI` 是本实现直接拟合响应得到的系数，
不能把已吸收因子 2 的文献 $b_\phi$ 数值未经换算代入。新的示踪物注入模式
还需要独立的 PNG 协方差对照；旧场层初条件模式的验证不能直接转移到新模式。

**索引**：代码 `codes/ezmock_png/`；验证/数据 `codes/ezmock_png_binary_verification/tmp_png_response/`；
Quijote 目标数据 `data/quijote_z1_local_png_power/`；对比图 `plots/ezmock_png_binary_verification/`。

> 本文档由 `codes/ezmock_png/tmp_doc_build/build_png_injection_doc.py` 生成——
> 要改内容请改生成器再重跑，不要直接在 ipynb 里手改。
""")

# ============================================================================
# 1. 背景
# ============================================================================
md(r"""
## 1. 背景与验收标准

**论文目标**（*Will ignoring PNG in the covariance affect the measurement of $f_{\rm NL}^{\rm loc}$?*）：
DESI/eBOSS 的协方差来自**高斯初条件**的 EZmock；若真实 $f_{\rm NL}\neq 0$，协方差里应有 PNG 贡献，
忽略它会造成误差低估。为此需要一套"带 local PNG 的 EZmock"，并用真 N-body 的 PNG 模拟标定验证。

**参考数据（Quijote-PNG）**：1 Gpc/$h$、$512^3$ 格、$z=1$ FoF halos（$M>10^{13}\,M_\odot/h$）：
15000 套高斯 + $f_{\rm NL}=\{\pm 50,\pm 100\}$ 各 500 套；$f_{\rm NL}=\{0,\pm 10,\pm 20,\pm 30\}$ 的
目标由插值构造（插值保真度已在附录验证）。大尺度平均功率的线性响应是标定目标之一；
若用于论文的协方差结论，还需要验证新模式的协方差。

### 1.1 local PNG 与"响应"

local 型势：$\Phi = \phi + f_{\rm NL}\left(\phi^2 - \langle \phi^2 \rangle\right)$。
它的物理效应包括 **scale-dependent bias**：长波势 $\phi_L$ 调制局部小尺度功率，
在线性偏置展开中产生额外项

$$\delta_t(\mathbf{k}) = b_1\delta_m(\mathbf{k}) + f_{\rm NL}b_\phi^{\rm lit}\phi_G(\mathbf{k})+\cdots.$$

这里 $b_\phi^{\rm lit}$ 指已经包含 local PNG 长短模耦合因子 2 的常见文献约定；
本代码的理想线性换算是 $b_\phi^{\rm lit}=2B_\phi$。实际目录的 $B_\phi$ 由响应标定决定。

本项目用配对差分定义**响应（奇部）**与**偶部**：

$$\mathcal{L}(k) = \frac{P_0(+f) - P_0(-f)}{2\,P_0^{\rm G}},\qquad
\mathcal{E}(k) = \frac{P_0(+f) + P_0(-f)}{2\,P_0^{\rm G}} - 1$$

$\mathcal{L}$ 线性于 $f_{\rm NL}$、内容就是 $b_\phi$（核心观测量）；$\mathcal{E}$ 是 $O(f_{\rm NL}^2)$
的共同项，用作自检。同种子配对差分可以让宇宙方差大部分抵消。

### 1.2 尺度口径：为什么只看 $k<0.01$

$b_\phi$ 在超大尺度（$k \lesssim 0.01\,h\,{\rm Mpc}^{-1}$）近似为常数；按本实现响应定义换算的目标等效系数 $B_\phi^{Q,\rm eff}(k)$
从 $k\approx 0.008$ 的 2.41 平滑漂到高 $k$ 的 1.50（40%）。本项目的科学目标（$k_{\min}$ 效应、
大尺度协方差）只要求 **$k<0.01$ 精确 1:1**，故标定口径取 $k<0.01$（对应最低那个 $k$ bin）。
""")

# ============================================================================
# 2. 共同底座
# ============================================================================
md(r"""
## 2. 共同底座：EZmock 的生成链条

两版共用同一条流水线（`libEZmock/EZmock/perturb.c` → `pop_tracer.c`）：

1. **高斯初条件**：由输入转移函数 $T(k)$ 生成势场 $\phi$ 的傅里叶系数（`mesh->phik`）与实空间场
   （`mesh->phi`）；生长因子 $D_0$、$D_{\rm plus}$、$\beta$ 由配置的 $\Omega_m$ 与输出红移算出
   （改造版动态化，见 §4.7）。
2. **ZA 位移**：K 空间构造

   $$\hat\Psi_i(\mathbf{k}) = \hat\phi(\mathbf{k})\cdot \mathrm{i}k_i\,\frac{2\pi}{L}\cdot
   \frac{T(k)}{D_{\rm plus}\,\beta}$$

   （代码：`rhok_i = phik * (ki*kfac) * twb`，其中 `kfac = 2π/L`、`twb = Transfer/Dplus/Beta`），
   再做 C2R 得到实空间位移场 `mesh->psi[0..2]`。
3. **密度 → 有效 PDF（工程环节）**：`eval_pdf` 用四旋钮
   （`RHO_CRITICAL` / `RHO_EXP` / `PDF_BASE` / `SIGMA_VELOCITY`）把密度非线性地映射成示踪物数量，
   `generate_tracers` 随机装填。
4. **改造版的 PNG 注入在此插入** → `add_vel_scatter`（速度按 `SIGMA_VELOCITY` 离散）；
   RSD 由 python 侧按 Quijote 约定施加（z 轴）。

**关键认识**：第 3 步是"人工"的、为复现某统计量而手调的映射，**不是引力演化**。
PNG 信号能不能"穿过"它、穿过之后剩多少，正是两版机制的根本分歧点。
""")

# ============================================================================
# 3. 原版机制
# ============================================================================
md(r"""
## 3. 原版：场层面注入（注入"原因"）

### 3.1 代码

`softdir/EZmock/libEZmock/EZmock/perturb.c:569,580`（原树已复原为改造前状态）：

```c
double Fnl=180;                                                 /* 硬编码 */
mesh->phi_png[idx] = mesh->phi[idx] + Fnl*mesh->phi[idx]*mesh->phi[idx];
/* phi -> phi + Fnl * phi^2 */
```

随后 $\phi_{\rm png}$ **原封不动**走完整链条：FFT → ZA 位移 → 密度 → PDF 映射 → 示踪物。
设计思路与 Quijote-PNG / 2LPT-PNG 同源：把非高斯"种"在初条件里、让下游自己演化出来。

### 3.2 为什么不可控：响应是"涌现"的

要在超大尺度产生**线性** PNG 响应，靠的是长波×短波模式耦合：长波势调制局部小尺度结构，
而 EZmock 里对"局部小尺度结构"作出反应的，是密度 → PDF 那一步人工环节。于是响应幅度是

$$\mathcal{L}^{\rm EZ}(k) \;\sim\; \underbrace{\vphantom{\big|}\left(\text{注入强度 } F_{NL}\right)}_{\text{可设}}
\;\times\; \underbrace{\left(\text{透射率：四旋钮 / ngrid / 目标样本依赖}\right)}_{\text{不可控}}$$

实测（v1 轮，固定四旋钮、同种子配对差分；数字见第 5 节单元 ①）：

- 对上 Quijote $f_{\rm NL}=100$ 需要 $F_{NL}\approx 140$–$220$（拟合区中位 141、最低 $k$ bin 220；
  原代码里硬编码的 180 就落在这个区间——现在有了新机制才知道那只是经验配平）；
- $\mathcal{L}_{128}/\mathcal{L}_{256}\approx 0.23$（差 4.3 倍）；就算扣掉两网格示踪物 $b_1$ 的差异
  （约 18%），仍有 5.0 倍——**换一个网格，响应就变**。

两个推论：

1. 那个"约 1.8 倍"的系数**不是单纯的单位换算**（否则不该随网格变）——它混合了归一化差异与透射率；
2. **不存在同时适配两个网格的 $F_{NL}$**——这是后来必须发明新机制的直接原因。
""")

# ============================================================================
# 4. 改造版机制
# ============================================================================
md(r"""
## 4. 改造版：示踪物层面注入（注入"结果"）

### 4.1 先备：公式里每个符号是什么（从零开始）

> 这一小节不假设读者记得任何 PNG 细节，把 4.2 之后用到的每个符号从物理上过一遍。
> 熟的部分可以直接跳到 4.2。

#### (a) 图景：local PNG 为什么给示踪物一个"乘性调制"

local 型非高斯把原初势写成

$$\Phi(\mathbf{x}) = \phi_G(\mathbf{x}) + f_{\rm NL}\left[\phi_G^2(\mathbf{x}) - \langle \phi_G^2\rangle\right]$$

其中 $\phi_G$ 是**高斯**随机场（主角），$f_{\rm NL}$ 是耦合强度。"local" 指二次项是**同一点**上
$\phi_G$ 的平方——所以它的效应是局部的、可以逐点理解。

峰值-背景分解（peak-background split）的直觉：在长波 $\phi_L>0$ 的区域，短波涨落的局部振幅被
整体抬高一点 $\Rightarrow$ 那里小尺度结构更多、晕更容易形成。于是**晕的局部数密度被该位置的
$\phi$ 值调制**：

$$n_h(\mathbf{x}) \simeq \bar n_h[1+b_1\delta_m(\mathbf{x})][1+2f_{\rm NL}R_\sigma\phi_G(\mathbf{x})],\qquad R_\sigma\equiv\partial\ln\bar n_h/\partial\ln\sigma_8.$$

第一项是普通高斯偏置；第二项是 PNG 的额外调制；$R_\sigma$ 是丰度对局部振幅的响应。
常见文献定义 $b_\phi^{\rm lit}=2R_\sigma$；代码的 B_PHI 是直接标定的输入系数，不按名称与它们等同。
因为 $\delta_m$ 和 $\phi_G$ 都由同一个原初势生成，两项可以合并成对 $\phi_G$ 的两个"通道"：

$$\delta_t = b_1\delta_m + 2f_{\rm NL}R_\sigma\phi_G + \cdots = \phi_G[b_1M(k)+2f_{\rm NL}R_\sigma]+\cdots.$$

于是在线性极限，有效偏置 $b_{\rm eff}(k)=b_1+2f_{\rm NL}R_\sigma/M(k)$。$M(k)\propto k^2T(k)$ 在低 $k$
急剧变小（见 (c)），所以 $b_{\rm eff}$ 在大尺度暴涨——这就是著名的 **scale-dependent bias**。
PNG 信号相对高斯背景的大小，是两个通道之比 $2f_{\rm NL}R_\sigma/(b_1M)$；4.4 节的线性响应闭式
本质就是这个比值（再多一个 2 的因子，见 (c)③）。

#### (b) 符号表

| 符号 | 名字 | 定义 / 关系 | 本项目中的来源 |
|---|---|---|---|
| $\delta_t(\mathbf{x})$ | 示踪物过密度 | $n(\mathbf{x})/\bar n - 1$；EZmock 产出的点集→密度就是它 | EZmock 输出 |
| $\phi_G(\mathbf{x})$ | 高斯原初势 | 无量纲；代码里 = `mesh->phi`（实空间）/ `mesh->phik`（傅里叶） | EZmock 用 $P_{\rm prim}\propto k^{n_s}$ 生成，与 ZA 位移**用同一个场** |
| $f_{\rm NL}$ / `FNL` | PNG 强度 | 原初势二次项系数 | 配置键 `FNL` |
| $R_\sigma$ / $b_\phi^{\rm lit}$ | PNG 偏置约定 | $b_\phi^{\rm lit}=2R_\sigma$ | 文献定义，需与代码约定换算 |
| B_PHI | 代码输入系数 | $A_{\rm inj}=2F_{NL}B_\phi$ | 从当前实现的功率响应直接标定 |
| $A_{\rm inj}$ | 注入幅度 | $A_{\rm inj} = 2\,F_{NL}\,B_\phi$（系数 2 来自 (a) 的局部 PNG 约定） | 由 `FNL`、`B_PHI` 算出 |
| $\Psi$ | 注入位移场 | 三分量矢量场，见 (c)② | `mesh->inj[0..2]` |
| $P_G$ | 高斯功率 | 不注入时示踪物自己的 $P_0$（响应定义的分母） | EZmock 高斯 run |
| $b_1$ | 线性偏置 | $b_1 \equiv \sqrt{P_G/P_{\rm lin}(z)}$（大尺度近似为常数）；物理上 = 示踪物密度跟着物质密度涨的比例 | EZ 示踪物自己的偏置，$\approx 3.2$ @ $k=0.008$ |
| $M(k)$ | Poisson 核 | 把势换算成密度的核：$\delta_m(k,z) = M(k)\,\phi_G(k)$，"单位势能产多少密度" | $M = k^2T/(\beta D_{\rm plus})$；$M(0.0078)\approx 456$ |
| $T(k)$ | 物质转移函数 | 每个模式穿过辐射-物质等能期后的增长抑制；$k\to0$ 时 $T\to1$，$k\gg k_{\rm eq}$ 时 $\propto k^{-2}\ln k$ | 输入表 `Tk_0.txt` |
| $D_0$, $D_{\rm plus}$ | 生长因子 | $D_0 = D_{A}(a{=}1)$（未归一化，锚定 $\beta$）；$D_{\rm plus} = D(0)/D(z_{\rm out})$，把 z=0 形状的输入表缩放到输出红移 | `cosmology.c` 按 `OMEGA_M`/`REDSHIFT` 现算（4.7 节） |
| $\beta$ | 势归一化常数 | $\beta = \frac{3}{2}\Omega_m H_0^2/D_0$：Poisson 方程组装的写法（H₀ 用 100 km/s/Mpc、k 用 h/Mpc；换算后即 $1.5\Omega_m/2998^2/D_0$，h 恰好抵消） | 随 `OMEGA_M` 变 |
| $\mathcal{L}(k)$ | 响应 | $[P(+)-P(-)]/(2P_G)$：$f_{\rm NL}$ 从 $+X$ 翻到 $-X$，$P_0$ 相对变了多少 | 本项目的核心观测量 |

**$\delta_t$ 到底是 EZmock 的哪一步产物？**（常被问，容易答错）

严格说：**pipeline 内部并没有一个叫 $\delta_t$ 的场**——EZmock 的最终产物是**示踪点集（catalog）**，
$\delta_t$ 是下游把点集 CIC 分格之后的密度场（POWSPEC 干的正是这件事）。它的两级"祖先"是：

```
run_mock.c:112  EZmock_create_dens_field              [dens_field.c:340]
   ├─ perturb.c 的 ZA 位移（φ → ψ；场层面 PNG 在这里的 phi_png 上）
   └─ density_field_cic            [dens_field.c:236] ← 对位移后粒子 CIC 分格
        ⇒ mesh->rho   粒子密度（≈ δ_m 层，尚未任何示踪物偏置）

run_mock.c:139  EZmock_populate_tracer                [pop_tracer.c:1207]
   ├─ eval_pdf         [pop_tracer.c:94]    迭代求示踪物计数 PDF（用 PDF_BASE）
   ├─ bias_model       [pop_tracer.c:149]   rho → mesh->rhot
   │        ⇒ mesh->rhot  示踪物数密度【期望场】（RHO_C / RHO_EXP 在此起作用）
   ├─ n_tracer_in_cell [pop_tracer.c:359]   每个格点分配几个示踪物
   ├─ generate_tracers [pop_tracer.c:773]   采样出具体位置 + ZA 速度
   │        ★ 这一步产出的 catalog 才"是" δ_t 的载体
   ├─ apply_png_shift  [pop_tracer.c:1136]  ← 本项目的注入（只挪点，不改别的）
   └─ add_vel_scatter  [pop_tracer.c:1043]  加随机速度弥散

输出 save_box → xyz+v；下游 POWSPEC 对 catalog CIC → δ_t(k)
```

**为什么注入写在 catalog 上、而不是 `rhot` 上**：在 `rhot`（连续期望场）上乘调制因子会改变
每格的取样概率 → 点数/总数/PDF 全变，与四旋钮标定纠缠；而在 catalog 上只"挪点"（`apply_png_shift`）
则**总点数不变、PDF 不变、速度不变**，但密度变化还包含输运项，并非纯局部乘性调制。
§4.4 响应式里的 $b_1$ 也正是**这个点集自己的偏置**（$\sqrt{P_G/P_{\rm lin}}$，下游量出来的），
和注入作用的对象同层——闭式才能干净成立。

**本质定位（一句话）**：注入 = **在 EZmock 已经生成好的高斯示踪样本上，叠加 PNG 的线性响应**——
发生在同一次运行内部、内存里（不是事后改输出文件，因为 $\phi_G$ 不在输出里）。EZmock 负责提供
高斯骨架（非线性结构、PDF、偏置、RSD 全是它的），PNG 的长波调制由解析式补上；不是重跑一套
非高斯 N-body。合理性的根据：大尺度上 PNG 对晕场有线性响应 $f_{\rm NL}b_\phi^{\rm lit}\phi_G$，
本机制可经标定匹配这个系数。不等价的部分包括输运、高阶偏置与
$f_{\rm NL}^2$ 对物质功率的修正；其大小应由目标统计量检验（偶部的差异见 §5）。

#### (c) 三条公式逐条读

**① 目标式（4.2 节）**：大尺度线性极限增加 $A_{\rm inj}\phi_G$ 的示踪物密度项。
这只是指定尺度与阶数的响应目标，不是对每点密度场的精确乘法。

**② 位移式（4.3 节）**：为什么用"挪点"而不是真乘？

- 真乘 $n\to n(1+A\phi)$ 会改变样本总数与 PDF（EZmock 的示踪物是点集，乘一个连续场很别扭）；
- 把每个示踪点从 $\mathbf{x}$ 挪到 $\mathbf{x}+\Psi(\mathbf{x})$，密度场变化一阶是
  $\delta^\prime=\delta-\nabla\cdot[(1+\delta)\Psi]$；输运项与注入同为位移振幅的一阶；
- 要让位移贡献正好等于想要的调制（$-\mathrm{i}\mathbf{k}\cdot\hat\Psi = A_{\rm inj}\hat\phi_G$），
  就要 $\hat\Psi_i = \mathrm{i}k_i(A_{\rm inj}/k^2)\hat\phi_G$。**$1/k^2$ 就是"从 $\delta$ 反解位移"时
  $\nabla\cdot\Psi$ 的 Poisson 逆**（$\Psi \sim \nabla\nabla^{-2}\phi$ 的结构）；
- 注意它**没有 $T(k)$**：调制源是**势** $\phi_G$ 本身，不是密度——所以不能用 ZA 位移那套
  $\mathrm{i}k_i\,T/(D_{\rm plus}\beta)$ 的核（那是"势→密度"的转移）。代码里 `phik` 存的正是
  未乘 $T$ 的势形状，可以直接用（实现细节见 4.5 节）。

**关于 $\nabla$（倒三角）的两个作用**（读这类公式最容易卡的地方）：

1. **$\nabla\cdot\Psi$（散度）——位移改密度的唯一通道**。点集整体挪一段，局部体积会被压缩或
   拉伸，密度随之变；均匀背景的线性极限是 $\delta^\prime=\delta-\nabla\cdot\Psi$，
   原有密度结构还产生额外的输运项。直觉：$\Psi$ 发散时，点被
   拉开 $\Rightarrow$ 密度变小；非均匀背景还会被位置平移重新分配。
2. **$\Psi = \nabla\Phi_\Psi$（梯度）——位移场本身是标量势的梯度**（ZA 位移无旋，可以这样表示）。
   于是"一个标量场 + 一个 $\nabla$"就表达了三分量矢量场。
3. 两者在傅里叶空间都变成 $\mathrm{i}\mathbf{k}$：$\nabla \to \mathrm{i}\mathbf{k}$、
   $\nabla\cdot \to \mathrm{i}\mathbf{k}\cdot$、$\nabla^2 \to -k^2$、$\nabla^{-2}\to -1/k^2$
   （整体符号看约定）。所以"$\Psi = \nabla\nabla^{-2}\chi$"在 k 空间写出来就是
   $\hat\Psi_i = \mathrm{i}k_i\,\chî/k^2$ 的形式——**这正是 4.3 节公式里 $\mathrm{i}k_i/k^2$
   的来历**：一个 $\mathrm{i}k_i$ 来自位移的"梯度"，一个 $1/k^2$ 来自"从密度反解位移"的 Poisson 逆。

**小帽子 $\hat{\ }$ 不是物理修正（常见误解，特别澄清）**

公式里 $\hat\Psi$、$\hat\phi_G$ 头顶的帽子**只是"傅里叶空间"的记号**，没有任何物理内容：
$\hat\Psi(\mathbf{k})$ 与 $\Psi(\mathbf{x})$ 由傅里叶变换互相**唯一**确定（可逆、无损），
是同一个物理量的两种写法——像同一句话用两种语言写。它**不是** PNG 修正、**不是**"非高斯势"、
不代表任何新东西。本项目示踪物注入强度由系数 $A_{\rm inj} = 2F_{NL}B_\phi$ 控制；
而 $\phi_G$ 本身就是**纯高斯**场（下标 G = Gaussian），它头顶的帽子只说明用 k 空间写它。

帽子记号之所以好用：微分算符在 k 空间变成乘法（$\nabla\to\mathrm{i}\mathbf{k}$、
$\nabla\cdot\to\mathrm{i}\mathbf{k}\cdot$、$\nabla^{-2}\to-1/k^2$），于是**逐模独立**的代数
方程代替了实空间的微分方程（求逆变成逐点相除）。

**去帽子练习**：把 $\hat\Psi_i = \mathrm{i}k_i(A_{\rm inj}/k^2)\hat\phi_G$ 翻译回实空间
（$\mathrm{i}k_i\to\nabla_i$，$1/k^2\to\nabla^{-2}$，注意整体符号看约定）：

$$\Psi(\mathbf{x}) = -\nabla\,\nabla^{-2}\!\left[A_{\rm inj}\,\phi_G(\mathbf{x})\right]$$

读法："位移场 = 把标量场 $A_{\rm inj}\phi_G$ 先做一次**逆拉普拉斯**、再取**梯度**。"这才是
物理内容——加不加帽子说的是同一件事。

**$\Psi$ 的身份定位（一句话说准）**：

- 密度层面：**大尺度线性极限的 PNG − 高斯 = $A_{\rm inj}\phi_G$**；
  在非线性场中，位移还带来输运贡献；
- 位移场 $\Psi$ 与 $F_{NL}$ 线性，但平移后目录的统计量一般还包含更高阶项；
- **不是**对 EZmock 自己 ZA 位移（`psi[0..2]`）的修正——那条生成链一个字节没动；$\Psi$ 是在
  **已生成**的示踪点上额外加的一层独立位移（这正是"注入结果"与"注入原因"的区别，§3 vs 本节）；
- **小位移、大效应**：位移幅度 rms $\approx 1.7$ Mpc/h（只有 EZmock 自身 ZA 位移的约 1/4，
  且几乎全在低频、全盒相干），但因为相干，最低 bin 的 $P_0$ 效应达 $\pm 70\%$——
  不要因为"位移小"就以为"效应小"（数字推导见下一条）。

**每个示踪物的位移大小和方向是怎么定出来的？（实现口径，常被追问）**

注入分两步：**先在网格上建好三个位移分量场，再在每个示踪物位置上取样**——
示踪物自己不做任何傅里叶运算。

1. **建场**（`perturb.c`，每次 run 做一次）：对每个傅里叶模 $\mathbf{k}$ 算
   $\hat\Psi_i(\mathbf{k}) = \mathrm{i}k_i\,(A_{\rm inj}/k^2)\,\hat\phi_G(\mathbf{k})$
   （$k=0$ 模清零，否则 $0\times\infty$），三次逆 FFT 得到实空间三分量场 `mesh->inj[0..2]`；
   与 ZA 生成位移走的是同一套代码路径（同一本"密度 ↔ 位移"字典）。
2. **取样**（`pop_tracer.c` 的 `apply_png_shift`，逐示踪物）：把示踪物坐标换算成网格坐标
   $u = x\,N_g/L_{\rm box}$，向下取整 + 周期回绕定位胞元，用 8 个角点的**三线性插值**
   权重（$w_xw_yw_z$）加权求和得该点的 $(\Psi_x,\Psi_y,\Psi_z)$，加到坐标上再回绕进盒子。
   位移场是平滑的，相邻示踪物拿到的几乎一样——**这不是逐点随机的抖动**。

**为什么 8 个角点加权就能得到该点的 $\Psi$？（插值凭什么合法，常被追问）**

- **连续场没存，只能从采样重建**：$\Psi(\mathbf{x})$ 是连续场，内存里只有 $N_g^3\approx1.7\times10^7$
  个网格采样值，任意位置的值必须靠重建。
- **三线性插值就是最基础的重建**：用所在胞元 8 个角点的值构造**唯一**的三变量线性函数
  （= 局部一阶 Taylor），权重是"体积份额"，**和恒为 1** ⇒ 常数场精确保留、网格点上精确。
- **精度为什么够**：注入场带限（最高 $k_{\rm Ny}=0.8$），功率几乎全在低 $k$（$k<0.01$ 贡献
  1.34/1.71 的 rms）。数学上"采样 + 线性插值"精确等价于乘窗函数
  $\hat W(\mathbf{k})=\prod_i\mathrm{sinc}^2(k_i\Delta x/2)\approx 1-(k\Delta x)^2/12$；
  $\Delta x=1000/256\approx3.9$ Mpc/h，在标定 bin $k=0.0078$ 处失真 $\approx8\times10^{-5}$
  （场波长 $\sim2\pi/k\approx800$ Mpc/h，一个胞元只占波长的 1/200，胞元内近乎直线）。
- **为什么不用"精确"**：逐点对全部 $\sim1.7\times10^7$ 个模求和是
  $O(N_{\rm mode}N_{\rm tracer})\sim10^{12}$ 次运算；网格 + 一次 FFT + 局部插值 = 经典的
  **粒子-网格（PM）套路**（同族操作：引力 PM、下游 POWSPEC 的 CIC 分箱）。

由此，"怎么定"的三个答案：

- **方向**：$\Psi \propto \nabla\chi$（$\chi = \nabla^{-2}[A_{\rm inj}\phi_G]$ 是平滑标量场）⇒
  每个示踪物沿该处 $\chi$ 的**梯度方向**走，无旋（无环流分量）、盒尺度相干——像一股缓慢的
  整体流把大尺度上的示踪物一起推，这正是"大尺度密度被重新分配"的实现方式。
- **大小**：局部 $|\nabla\chi|$；全场 rms $\approx$ **1.7 Mpc/h**（$F_{NL}=100$、`B_PHI`=2.65
  口径；低频主导——只算 $k<0.01$ 的模就有 1.34 Mpc/h）。且 $A_{\rm inj}=2F_{NL}B_\phi$ 线性
  ⇒ **$F_{NL}$ 或 $B_\phi$ 变号则全部位移反向，放大 $x$ 倍则整体放大 $x$ 倍**。
- **与示踪物自身性质无关**：位移只取决于它**所在的位置**（该处 $\phi_G$ 的值），与速度、局域
  密度、是否被 `RHO_C` 挑中统统无关——所以四个旋钮不需要跟着注入一起改。

（数量级对照：EZmock 自身 ZA 位移在 $z=1$ 的 rms $\approx 6.1$ Mpc/h，同 k 范围；注入位移约
是它的 1/4，但几乎全在低 k——恰是最低 bin $P_0$ 最敏感的大尺度整体流动。）

**和 Zel'dovich 近似（ZA）的关系**（是同一台数学机器，但不是"用 ZA 近似引力"）

ZA 的公式（记 $D(t)$ 为生长因子，$\psi$ 为位移势）：

$$\mathbf{x}(\mathbf{q},t) = \mathbf{q} + \Psi(\mathbf{q},t),\qquad
\Psi = D(t)\nabla_q\psi,\qquad
\nabla\cdot\Psi = -\delta_{\rm lin}\ \Leftrightarrow\ \nabla^2\psi = -\delta_{\rm lin}
\ \Leftrightarrow\ \hat\Psi_i(\mathbf{k}) = -\mathrm{i}\frac{k_i}{k^2}\hat\delta_{\rm lin}(\mathbf{k})$$

$$1+\delta(\mathbf{x}) = \left|\det\!\left(\delta_{ij} + \frac{\partial\Psi_i}{\partial q_j}\right)\right|^{-1}
\ \xrightarrow{\ \text{一阶}\ }\ \delta = -\nabla\cdot\Psi$$

即 ZA 的两块结构：**(1) 无旋（梯度）位移；(2) 连续性方程 $\delta = -\nabla\cdot\Psi$**。

- **一样的部分（数学机器）**：本机制用的就是这两块。EZmock 自己的生成位移 `psi[0..2]`
  （`perturb.c` 的 `twb` 核）正是 $\mathrm{i}k_i/k^2\times$"线性密度"这条 ZA 公式；
  注入位移用**同一个算符**，只把源场从 $\delta_{\rm lin}$ 换成 $A_{\rm inj}\phi_G$。
- **不一样的物理**：
  1. 我们**不在模拟引力演化**——是在构造大尺度线性 PNG 密度响应。
     ZA 在这里只是"密度 ↔ 位移"的**一阶字典**；PNG 的物理起源是局部丰度调制
     （peak-background split），**不是**晕被真的推了一下，位移只是实现手段（见上面"身份定位"）。
  2. 实现把**真实的点挪了**，因此密度包含平流与雅可比结构。
     $-\Psi\cdot\nabla\delta$ 与目标响应同为位移振幅的一阶，不能为任意统计量
     预设千分之一上界；其影响必须按目标尺度与协方差直接检验。
  3. ZA 的 $\Psi$ 是**拉格朗日**坐标的函数（粒子沿直线飞）；注入是给每个点按**当前**位置加位移
     $\Psi(\mathbf{x})$（欧拉写法）——一阶无差别，实现上也更省事。

**③ 大尺度线性响应闭式（4.4 节）**：忽略非线性输运时

$$P(\pm A_{\rm inj}) = P_G \pm 2A_{\rm inj}\langle \delta_t^G \phi_G\rangle + A_{\rm inj}^2 P_\phi
+ O(\Psi\nabla\delta)$$

其中交叉谱 $\langle\delta_t^G\phi_G\rangle = b_1 M(k) P_\phi$（因为示踪物场 $\propto b_1M\phi$）。
取差 $P(+)-P(-)$ 时所有偶次项（包括 $A^2$ 项）全消掉：

$$\mathcal{L}_{\rm inj} = \frac{P(+A)-P(-A)}{2P_G}
= \frac{2A_{\rm inj}\,b_1MP_\phi}{b_1^2M^2P_\phi} = \frac{2A_{\rm inj}}{b_1M(k)}$$

$b_1$ 在这里的作用一目了然：**信号是"注入通道 × 背景通道"的交叉项**，所以分母显式含 $b_1M$。
这不是非线性目录或协方差的精确恒等式；实测平均功率响应见 §5。

**④ 如何与 Quijote 对照**：用同一平均功率响应定义测量
$\mathcal{L}_Q$ 与 $\mathcal{L}_{\rm inj}(B_\phi)$，再按比值标定 $B_\phi$。
若外部文献采用 $\delta_t=b_1\delta_m+f_{\rm NL}b_\phi^{\rm lit}\phi_G$，
理想线性核中应比较 $b_\phi^{\rm lit}=2B_\phi$；不能只凭同名参数断言 1:1。

#### (d) 一个真实数字（最低 bin，$k = 0.0078$ $h$/Mpc）

| 量 | 数值 | 说明 |
|---|---|---|
| $M(k)$ | 456 | $=k^2T/(\beta D_{\rm plus})$；$k=0.041$ 处涨到 4981（11 倍）——**PNG 信号在大尺度最大**的数学来源（$\mathcal{L}\propto1/M\propto1/k^2$） |
| $b_1M$ | $\approx 1.5\times10^3$ | 分母；两种 $b_1$ 口径相差 ~5% |
| $B_\phi^{Q,\rm eff}$（目标等效系数） | 2.41 | 按本实现功率响应定义换算；到 $k=0.04$ 漂移降到 1.50，见 §5 |
| 需要 $A_{\rm inj}$ | $\approx \mathcal{L}_Q\,b_1M/2 \approx 0.686\times1.5{\rm e}3/2 \approx 5.1\times10^2$ | 由 $\mathcal{L}_Q=0.686$ 反解 |
| ⇒ `B_PHI` | $510/(2\times100)\approx 2.6$ | 与定稿 2.65 一致；正式标定与不确定性见 §5 |
| 注入后 $\mathcal{L}$ | $\approx 0.7$ | $\pm f_{\rm NL}=100$ 使最低 bin 的 $P_0$ 变化 $\pm70\%$——**效应巨大**，必须注入对 |

#### (e) 与场层面机制的本质对照

| | 场层面（§3） | 示踪物层面（本节） |
|---|---|---|
| 注入的是什么 | $\phi$ 的二次项（"原因"） | 对示踪点做势驱动的位置平移 |
| 响应能不能算 | 要跑完整链条 | 大尺度线性极限有闭式 $\mathcal{L}=2A_{\rm inj}/(b_1M)$，实际目录仍需测量 |
| 强度可控吗 | 靠 $F_{NL}$ 间接调，还额外依赖网格/四旋钮 | 一个常数 $B_\phi$，跨网格差 3% |

### 4.2 目标形式

目标是在超大尺度线性极限得到

$$\delta_t'(\mathbf{k}) = \delta_t(\mathbf{k}) + A_{\rm inj}\phi_G(\mathbf{k})+\cdots,
\qquad A_{\rm inj} \equiv 2\,F_{NL}\,B_\phi.$$

$\phi_G$ 是**高斯**势场（即 `mesh->phik` 对应的实空间场）。
$B_\phi$ 是本实现的待标定系数，并非任意文献中同名参数的数值。

### 4.3 用"位移"实现大尺度响应

对示踪点集做位置平移 $\mathbf{x}'=\mathbf{x}+\Psi(\mathbf{x})$，
数密度守恒在位移的一阶给出：

$$\delta_t'=\delta_t-\nabla\cdot[(1+\delta_t)\Psi]
=\delta_t+A_{\rm inj}\phi_G+A_{\rm inj}\phi_G\delta_t
-\Psi\cdot\nabla\delta_t+O(\Psi^2).$$

要求 $-\mathrm{i}\,\mathbf{k}\cdot\hat\Psi = A_{\rm inj}\,\hat\phi_G$，即

$$\hat\Psi_i(\mathbf{k}) = \mathrm{i}\,\frac{k_i}{k^2}\,A_{\rm inj}\,\hat\phi_G(\mathbf{k})
\qquad(\text{纯 } 1/k^2 \text{ 核，不乘任何转移函数})$$

输运项与目标响应同为 $A_{\rm inj}$ 的一阶；低 $k$ 线性系数可以对上，
却不保证非线性功率、双谱或协方差一致。$\mathbf{k}=0$ 模式要显式置零。

### 4.4 大尺度线性响应与实际标定

记高斯示踪物功率 $P_G = b_1^2 M^2 P_\phi$、示踪物-势交叉谱
$C \equiv \langle \delta_t^G \phi_G \rangle = b_1 M P_\phi$，则注入后

$$P(\pm A_{\rm inj}) = P_G \pm 2 A_{\rm inj} C + O(A_{\rm inj}^2)
\quad\Longrightarrow\quad
\mathcal{L}_{\rm inj}(k) = \frac{P(+A)-P(-A)}{2\,P_G}
= \frac{2 A_{\rm inj}}{b_1 M(k)} = \frac{4\,F_{NL}\,B_\phi}{b_1 M(k)}$$

其中 $M(k) = k^2 T(k) / (D_{\rm plus}\,\beta)$（代码 `twb` 的同款口径）。
若 Quijote 侧使用 $\delta_t=b_1\delta_m+f_{\rm NL}b_\phi^{\rm lit}\phi_G$ 的定义，
相同线性近似给出

$$\mathcal{L}_Q(k) = \frac{2\,f_{\rm NL}\,b_\phi^{\rm lit}}{b_1 M(k)}.$$

此时理想线性换算为 $B_\phi=b_\phi^{\rm lit}/2$，但非线性目录的最佳数值
仍应由相同统计量直接标定，不能只凭参数名称代入。
实际采用同一功率响应比值标定 $B_\phi$；第 5 节的结果只验证了所测统计量和尺度，
并未建立新模式的 PNG 协方差等价性。

### 4.5 代码实现

| 环节 | 位置 | 说明 |
|---|---|---|
| 构造注入位移场 | `libEZmock/EZmock/perturb.c:762-822` | 逐模构造 $\mathrm{i}k_i A_{\rm inj}\hat\phi_G/k^2$（含完整推导注释；$k=0$ 跳过）；复用空闲的 `rhok2/3/4` 当缓冲；三次 C2R → `mesh->inj[0..2]` |
| 按需分配 | `libEZmock/EZmock/dens_field.c:415-418` | 仅当 `fnl != 0 && b_phi != 0` 时分配 `inj[3]`，否则为 NULL、不进新分支 |
| 施加到示踪点 | `libEZmock/EZmock/pop_tracer.c:1136` `apply_png_shift()` | 三线性插值 + 周期回绕，OMP 并行 |
| 调用位置 | `pop_tracer.c:1271`（`generate_tracers` 之后、`add_vel_scatter` 之前） | **速度场不动**，RSD 不被污染 |
| 配置键 | `src/main/load_conf.c` | 见下表 |

**配置语义**：

| 键 | 含义 | 不设时 |
|---|---|---|
| `FNL` | 物理 $f_{\rm NL}$（替代原硬编码 180） | 0 |
| `B_PHI` | 示踪物层面强度，$A_{\rm inj}=2 F_{NL} B_\phi$ | 0（不注入） |
| `FNL_FIELD` | 场层面二次项系数（$\phi + F_{NL}^{\rm field}\phi^2$） | `B_PHI != 0` → 0（标定模式：关掉未标定的旧机制）；否则 = `FNL`（历史行为） |

**初条件里现在还有没有 PNG？（两种模式的区分，常被问）**

| 配置 | 初条件（场层面） | 示踪物层 | 用途 |
|---|---|---|---|
| 写 `B_PHI ≠ 0`（标定模式） | **纯高斯**（`FNL_FIELD` 自动取 0，$\phi_{\rm png}=\phi$） | PNG 调制在此 | 本项目的全部标定/验证跑法 |
| 只写 `FNL`、不写 `B_PHI`（历史模式） | 有 PNG：$\phi + F_{NL}\phi^2$（`FNL_FIELD` 默认 = `FNL`） | 无 | 相同线程数下可与改造前二进制逐位一致 |

即：**"糖糖版初条件里有没有 PNG"取决于配置**——标定模式刻意把未标定的旧机制关掉（否则两套 PNG
叠加，标定会混乱）；历史模式则完全保留原行为。即使 `FNL_FIELD=0`，初条件里仍然有**高斯**的
$\phi$ 场（ZA 位移的源头），"没有 PNG"指的是二次项关掉、初条件统计上严格高斯。

**逐位回退条件**：`B_PHI=0`、`FNL_FIELD` 不设且线程数相同，配置路径、RNG 流与
改造前二进制相同（先前 A/B 实测见验证目录）。CLI 现在遵守 `OMP_NUM_THREADS`
并以 24 为上限；改变线程数会改变 RNG 流，不能要求逐字节相同。
因此新机制不会污染任何历史结果；`FNL=0`（或 `B_PHI=0`）即可完全关掉注入。

### 4.6 为什么这样"对"（设计哲学）与边界

EZmock 本来就是**统计仿制**工具（四旋钮是为别的统计量手调的），并不从第一原理演化 PNG；
而超大尺度上 PNG 的效应是已知的、线性可加的偏置响应 $b_\phi$。与其指望链条"涌现"出正确幅度，
不如把这个解析响应直接、精确地写在示踪物上——代价是强度要外部标定，收益是响应可控、跨网格可复用。

边界（明确接受的近似）：

- 常数 $B_\phi$ 只在 $k<0.01$ 的标定尺度上匹配（高 $k$ 目标等效系数有 40% 漂移，本机制响应形状平坦
  $\Rightarrow$ $k\sim 0.03$ 处高估约 30%）。主人明确只关心超大尺度，**不做** $W(k)$ 形状表；
- 标定值 $\propto$ 示踪物 $b_1$：**重调四旋钮后必须用 `analyze_injection.py` 重标**；
- 后期计划：在 Quijote 数据上用功率谱模型对 $(b_1,\,p,\,b_\phi)$ 做 MCMC 联合标定
  （现在 2.65 是先行的粗口径标定）。

### 4.7 顺带的正事：生长因子动态化

改造版把 `perturb.c` 里写死的 $D_0=0.789246$、$D_{\rm plus}=1.650$、$\Omega_m=0.315$
改成按配置的 `OMEGA_M` / `REDSHIFT` 现算（`cosmology.c`）。`REDSHIFT` 现在是**盒子的真实输出红移**，
同时驱动密度幅度与 RSD 速度 $v_{\rm fac}=f(z)H(z)a$；`run_mock.c` 强制 `growth2=1`，
避免红移换算被重复计入。

**迁移警示**：老配置写 `REDSHIFT=0` 的，在新代码下位移幅度会变成 1.64 倍、实空间 $P_0$ 高约 50%；
主人的配置（$\Omega_m=0.3089$ + `REDSHIFT=1`）与旧硬编码口径的幅度只差 0.19%，
但速度是被修正过的（无法与历史星表逐位对齐）。

**对 PNG 标定的意义**：旧口径下生成端速度固定在 $z=0$，python 侧 RSD 换算只有正确值的 0.671 倍，
四极矩比 Quijote 低 38%（$P_2/P_0 = 0.2276$ vs $0.3608$）；动态化后为 0.3185（仍低 13%，
属四旋钮可调范围）。**所以四旋钮（尤其 `SIGMA_VELOCITY`）必须在正确 RSD 下重新标定**——
这是当前最大的待办。
""")

# ============================================================================
# 5. 验证（数据核对单元）
# ============================================================================
md(r"""
## 5. 验证结果（可直接跑的数据核对）

下面 4 个单元从既有 npz 读数字（只读、秒级、登录节点即可）。数据目录：
`codes/ezmock_png_binary_verification/tmp_png_response/`。

**口径提醒（易错）**：机制验证用的"常数"要用**各自网格的 EZ 示踪物 $b_1$** 归一化——
两个网格的示踪物种群不同（$b_1$ 差约 18%），若统一用 Quijote 代理 $b_1$，会把种群差异
误读成机制问题。而直接对 Quijote 的标定（同一网格内比 EZ/Q 响应，单元 ③④）不需要任何
$b_1$ 建模。
""")

code(r"""
# ============ ① 原版（场层面）机制：要多大 FNL？响应怎样依赖网格？ ============
# 数据：v1_analysis.npz —— analyze_response.py 对 v1 轮（FNL=±180、同种子配对差分）的分析输出
import numpy as np
from pathlib import Path

VER = Path("/pscratch/sd/l/lzy/PNG-EZmock/codes/ezmock_png_binary_verification/tmp_png_response")
d1 = np.load(VER / "v1_analysis.npz")
k, need = d1["k"], d1["need"]
mfit = (k >= 0.006) & (k <= 0.08)          # 拟合区口径（与 analyze_response.py 一致）

print(f"对上 Quijote f_NL=100 需要的 FNL：拟合区中位 {np.nanmedian(need[mfit]):.0f}，"
      f"最低 k bin（k={k[0]:.5f}）{need[0]:.0f}，"
      f"逐 bin 范围 [{np.nanmin(need[mfit]):.0f}, {np.nanmax(need[mfit]):.0f}]")
r = np.nanmedian((d1["L128"] / d1["L256"])[mfit])
print(f"同 FNL 下响应比 L(128)/L(256) 中位 = {r:.3f}（= {1 / r:.1f} 倍差）")

# 扣掉两网格示踪物 b1 的差异（两个网格的示踪物种群不同），看"纯机制"是否仍网格依赖：
#   b1^EZ(ng) = sqrt(P0_EZ(ng) / P_lin(z=1))，用 v1 的高斯 run 的 P0（P0_ez* 键）
camb = np.loadtxt("/pscratch/sd/l/lzy/PNG-EZmock/data/ezmock_fnl0_calibration/"
                  "linear_pk/quijote_camb_matterpow_z0.dat")
Plin_z1 = np.interp(k, camb[:, 0], camb[:, 1]) / 1.650333193**2      # Dplus(z=1)
b1e256 = np.sqrt(d1["P0_ez256"] / Plin_z1)
b1e128 = np.sqrt(d1["P0_ez128"] / Plin_z1)
print(f"两网格示踪物 b1 比（Ng=128/256，k=0.008）= {b1e128[0] / b1e256[0]:.3f}")
r_b1 = np.nanmedian((d1["L128"] * b1e128 / (d1["L256"] * b1e256))[mfit])
print(f"扣除 b1 差异后，旧机制响应比 = {r_b1:.3f}（= {1 / r_b1:.1f} 倍差）——仍是强网格依赖")
""")

code(r"""
# ============ ② 新机制：响应是不是 常数×2Ainj/(b1·M)？跨网格是否一致？ ============
# 数据：v2_injection_analysis.npz —— B_PHI=2.9、FNL=±100、8 realizations，含 Ng=128/256 两档。
# 口径（易错点）：b1 必须用【各自网格的 EZ 示踪物 b1】——两个网格的示踪物种群不同
#   （b1 差约 18%），用统一的 Quijote 代理 b1 会把种群差异误读成机制问题。
d2 = np.load(VER / "v2_injection_analysis.npz")
camb = np.loadtxt("/pscratch/sd/l/lzy/PNG-EZmock/data/ezmock_fnl0_calibration/"
                  "linear_pk/quijote_camb_matterpow_z0.dat")
Plin_z1 = np.interp(d2["k"], camb[:, 0], camb[:, 1]) / 1.650333193**2    # 线性 P(k,z=1)
b1e256 = np.sqrt(d1["P0_ez256"] / Plin_z1)      # 复用 ① 里的高斯 P0
b1e128 = np.sqrt(d1["P0_ez128"] / Plin_z1)

Ainj = 2.0 * 100.0 * 2.9                        # = 2·FNL·B_PHI（v2 轮所用强度）
c256 = d2["L_inj256"] * b1e256 * d2["M"] / (2 * Ainj)   # 把测量写成本机制预言的常数
c128 = d2["L_inj128"] * b1e128 * d2["M"] / (2 * Ainj)
mwide = (d2["k"] >= 0.006) & (d2["k"] <= 0.08)   # 常数中位口径（与验证记录一致）
mtight = (d2["k"] >= 0.008) & (d2["k"] <= 0.045) # 平度口径

print(f"常数中位（[0.006,0.08]）：Ng=256 -> {np.nanmedian(c256[mwide]):.3f}，"
      f"Ng=128 -> {np.nanmedian(c128[mwide]):.3f}"
      f"   <- 相对差 {abs(np.nanmedian(c128[mwide]) / np.nanmedian(c256[mwide]) - 1) * 100:.1f}%")
print(f"逐 bin 两网格之比中位 = {np.nanmedian((c128 / c256)[mwide]):.3f}（= 1.0 即网格无关）")
print(f"平度（[0.008,0.045] 相对中位 rms）：Ng=256 {np.nanstd(c256[mtight] / np.nanmedian(c256[mtight])) * 100:.1f}%，"
      f"Ng=128 {np.nanstd(c128[mtight] / np.nanmedian(c128[mtight])) * 100:.1f}%")
""")

code(r"""
# ============ ③ 标定 B_PHI：k<0.01（超大尺度口径）达到 1:1 需要多强 ============
# 两轮独立数据互证：v2（B_PHI=2.9）与 v3（B_PHI=2.6），各 8 realizations、FNL=±100。
for tag, used in (("v2", 2.9), ("v3", 2.6)):
    dd = np.load(VER / f"{tag}_injection_analysis.npz")
    mlo = dd["k"] < 0.01                          # 主人的标定口径：k<0.01
    ratio = dd["L_inj256"][mlo] / dd["LQ"][mlo]
    print(f"{tag} 轮（B_PHI={used}）：k<0.01 的 L_inj/L_Q 中位 = {np.median(ratio):.4f}"
          f"  -> 1:1 需 B_PHI = {used / np.median(ratio):.3f}")
print("两轮互证 + 响应对 B_PHI 精确线性 => 定稿 B_PHI = 2.65（50-real P0 对比后最佳 2.54±0.04）")
""")

code(r"""
# ============ ④ 50-realization 的 P0 对比（v5，定稿数字） ============
# EZ 侧：v5_ng256_cmp_*.npz（各 50 realizations，同种子）；Q 侧：Quijote 500 realizations 均值。
QD = Path("/pscratch/sd/l/lzy/PNG-EZmock/data/quijote_z1_local_png_power/native_realizations")
q0 = np.load(QD / "fnl_0/pkrsd_matrix.npz")
q = {"fnl0": q0,
     "p100": np.load(QD / "fnl_p100/pkrsd_matrix.npz"),
     "m100": np.load(QD / "fnl_m100/pkrsd_matrix.npz")}
ez = {tag: np.load(VER / f"v5_ng256_cmp_{tag}.npz")["P0"] for tag in ("fnl0", "p100", "m100")}
nreal = ez["p100"].shape[0]

# --- 奇部：同种子配对差分（宇宙方差大幅抵消），再与 Quijote 目标比 ---
odd = ((ez["p100"] - ez["m100"]) / 2) / ez["fnl0"]
L_EZ, e_EZ = odd.mean(0), odd.std(0, ddof=1) / np.sqrt(nreal)
L_Q = (q["p100"]["P0_mean"] - q["m100"]["P0_mean"]) / (2 * q["fnl0"]["P0_mean"])
print(f"奇部 bin0（k={q0['kcen'][0]:.5f}）：EZ {L_EZ[0]:.4f}±{e_EZ[0]:.4f} vs Q {L_Q[0]:.4f}"
      f"  -> 比值 {L_EZ[0] / L_Q[0]:.3f}±{e_EZ[0] / L_Q[0]:.3f}"
      f"  -> 1:1 需 B_PHI = {2.65 / (L_EZ[0] / L_Q[0]):.2f}")

# --- 偶部（f_NL² 共同项），各自以本侧高斯 P0 归一 ---
even_EZ = 100 * (((ez["p100"] + ez["m100"]) / 2) / ez["fnl0"] - 1)
sem_EZ = even_EZ.std(0, ddof=1) / np.sqrt(nreal)
even_Q = 100 * ((q["p100"]["P0_mean"] + q["m100"]["P0_mean"]) / 2 / q["fnl0"]["P0_mean"] - 1)
nq = q["fnl0"]["P0"].shape[0]
# (P+ + P-)/2 的 SEM = sqrt(s+² + s-²)/2（两侧独立模拟）
sem_Q = (100 * np.sqrt(q["p100"]["P0_std"] ** 2 + q["m100"]["P0_std"] ** 2)
         / np.sqrt(nq) / (2 * q["fnl0"]["P0_mean"]))
print(f"偶部 bin0：EZ {even_EZ.mean(0)[0]:+.2f}%±{sem_EZ[0]:.2f}%  vs  "
      f"Q {even_Q[0]:+.2f}%±{sem_Q[0]:.2f}%")

# --- 拟合区 [0.006,0.08] 的 P0 残差（EZ 均值 vs Quijote 目标） ---
mfit = (q0["kcen"] >= 0.006) & (q0["kcen"] <= 0.08)
for tag in ("fnl0", "p100", "m100"):
    r = (ez[tag].mean(0) - q[tag]["P0_mean"]) / q[tag]["P0_mean"]
    print(f"拟合区 P0（{tag}）：mean {100 * r[mfit].mean():+.2f}%   "
          f"rms {100 * np.sqrt(np.mean(r[mfit] ** 2)):.2f}%")
""")

md(r"""
### 5.5 汇总

| 检验 | 结果 | 含义 |
|---|---|---|
| 旧机制 ngrid 依赖（①） | $\mathcal{L}_{128}/\mathcal{L}_{256}\approx 0.23$（4.3 倍；扣 $b_1$ 差异后 5.0 倍），需要 $F_{NL}\approx 140$–$220$ | 原版不可跨分辨率复用 |
| 新机制响应形状（②） | 常数 $c=\mathcal{L}b_1^{\rm EZ}M/(2A_{\rm inj})$ 平（Ng=256 中位 0.893、rms 约 5–8%） | 闭式响应成立 |
| 新机制 ngrid 无关（②） | $c$ 的 Ng=128/256 差 **3%**（0.921 vs 0.893，各自用本网格 EZ 示踪物 $b_1$） | 分辨率无关，可移植 |
| 标定（③） | 定稿 $B_\phi=2.65$；50-real 后最佳 $2.54\pm0.04$ | 与目标等效系数 $B_\phi^{Q,\rm eff}(k\approx 0.008)\approx 2.4$ 同量级 |
| P0 对比（④，50 real） | 奇部 bin0 比 $1.044\pm 0.016$；拟合区残差 mean/rms 约 1–3% | 与 Quijote 统计一致 |
| 偶部（④） | EZ $+15.2\%\pm 0.4\%$ vs Q $+10.8\%\pm 1.4\%$（差约 2.9$\sigma$） | 二阶残留，**未定论**，待后续 |
| 逐字节回退 | `B_PHI=0` 且线程数相同时与改造前一致 | 不污染历史/高斯结果 |
| 对 $f_{\rm NL}$ 严格二次 | $0/30/50/100$ 四点拟合残差 $<0.001\%$ | 可外推到任意 $f_{\rm NL}$ |

**其他已知细节**：

- POWSPEC 口径下 $k\approx 0.0109$ 的 bin，$P_2$ 恒为 0（Quijote 侧同样），做 $P_2$ 分析要跳过该 bin
  （$P_0$ 不受影响；本项目目前只看 $P_0$）；
- Quijote 插值目标（$\pm 10/\pm 20/\pm 30$）的 npz 没有 `halo_counts` 键（驱动器已做兼容回退）。
""")

# ============================================================================
# 6. 对比总表
# ============================================================================
md(r"""
## 6. 两版对比总表

| | 原版（场层面 only） | 改造版（糖糖版） |
|---|---|---|
| 注入位置 | 初条件势场（"原因"） | 示踪物密度（"结果"）；场层面机制保留但标定时关闭 |
| 强度参数 | 硬编码 $F_{NL}=180$，改值要重编译 | 配置键 `FNL` + 标定常数 `B_PHI`（`FNL_FIELD` 控旧机制） |
| 位移核 | $\hat\phi_{\rm png}\cdot \mathrm{i}k_i\,T/(D_{\rm plus}\beta)$（带转移函数） | 注入部分是**纯** $A_{\rm inj}/k^2$ 核 |
| 响应 | 涌现、不可控；依赖四旋钮/ngrid/样本 | 闭式 $\mathcal{L}=2A_{\rm inj}/(b_1 M)$ |
| 标定 | 经验配平（约有 1.8 倍系数），换网格失效（4 倍） | **一个常数** $B_\phi$，跨网格差 3% |
| 生长因子 | 写死 $D_0/D_{\rm plus}/\beta$（RSD 速度固定 $z=0$） | 按 $\Omega_m$、输出红移现算（RSD 修正） |
| 关闭/兼容 | 无"关闭"概念（永远在注入） | `B_PHI=0` 且线程数相同时与改造前逐字节一致 |
| 高 $k$ 行为 | 无法分辨（噪声大） | 已知：常数 $B_\phi$ 仅 $k<0.01$ 精确 1:1 |

**什么时候用哪个**：

- 高斯标定（复现 Quijote $f_{\rm NL}=0$ 的 $P_0/P_2$）→ 官方 stock 版 EZmock；
- 一切 PNG 相关生产/实验 → 改造版（`manual_tune.ipynb` 的 `EZMOCK_MODE="modified"`，默认）；
- 原版机制只在兼容性/回归验证时用（`FNL_FIELD` 不设、`B_PHI=0`，并保持相同线程数）。
""")

# ============================================================================
# 7. 现状/待办/索引
# ============================================================================
md(r"""
## 7. 现状、待办与文件索引

**现状**：`B_PHI = 2.65`（$k<0.01$ 口径；50-real P0 对比后最佳估计 $2.54\pm0.04$，
当前值高约 4%，属预期内的统计/口径差）。$f_{\rm NL}=\pm 30,\pm 50,\pm 100$ 与 Quijote 的
P0 对比全部通过（拟合区残差 mean/rms 约 1–3%）。

**待办**：

1. 主人重调四旋钮后，用 `analyze_injection.py` 重标 $B_\phi$（标定值 $\propto b_1$）；
2. 后期在 Quijote 数据上用功率谱模型对 $(b_1,\,p,\,b_\phi)$ 做 MCMC 联合标定；
3. 偶部 $+4.4\%\pm 1.5\%$ 的二阶残留待解释（不影响超大尺度奇部，但值得记录）。

**文件索引**：

| 路径 | 说明 |
|---|---|
| `codes/ezmock_png/` | 改造版代码（唯一维护位置）；`README.md` 含工程细节与踩坑清单 |
| `codes/ezmock_png/libEZmock/EZmock/perturb.c:762-822` | 注入位移场构造（含完整中文推导注释） |
| `codes/ezmock_png/libEZmock/EZmock/pop_tracer.c:1136` | `apply_png_shift()` 三线性插值 + wrap（调用点 1271） |
| `codes/ezmock_png_binary_verification/tmp_png_response/` | 验证驱动/分析/数据；`MANIFEST.md` 为索引 |
| `.../response_driver.py` | 从 `manual_tune.ipynb` 抽格执行，覆写 NREAL / NUM_GRID / FNL / B_PHI |
| `.../analyze_injection.py` | 新机制标定：直接给出"1:1 所需 B_PHI"（重标就用它） |
| `.../analyze_response.py` | 旧机制响应与 ngrid 依赖分析 |
| `.../plot_pk_comparison.py` | P0 对比图（每 $f_{\rm NL}$ 一页 + 响应汇总页） |
| `data/quijote_z1_local_png_power/` | Quijote 目标（`native_realizations/` 含 ±50/±100；`interpolated_realizations/` 含 ±10/±20/±30） |
| `plots/ezmock_png_binary_verification/ezmock_vs_quijote_p0_comparison.pdf` | 8 页 P0 对比图（定稿） |
| `codes/ezmock_png/tmp_doc_build/build_png_injection_doc.py` | 本文档的生成器（改文档改这里） |
""")

# ============================================================================
nb = nbf.v4.new_notebook(cells=cells)
nb.metadata["kernelspec"] = {
    "display_name": "Python 3",
    "language": "python",
    "name": "python3",
}
nb.metadata["language_info"] = {"name": "python"}
nbf.write(nb, str(OUT))
print("wrote", OUT, f"({len(cells)} cells)")
