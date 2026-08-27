# Optical Conductivity Note 设计

日期：2026-08-27

## 目标

在 DMRG branch 的 living note
`dmrg/report/dmrg_summary.tex` 中新增一个完整的线性光学响应章节。章节统一说明旧
ED 脚本中的显式 Kubo 本征态求和、当前 response-Lanczos/resolvent 方法，以及
regular、Drude 和 total conductivity 的实部与虚部。章节同时收录旧的线性响应
PDF 和当前经过验证的 `4x6, Np=4, m=2` FCI benchmark 图与数值。

本章只讨论线性光电导，完全不包含 shift current 或其他二阶响应。

## 修改位置与文档风格

- 修改 `dmrg/report/dmrg_summary.tex`。
- 新章节放在 `Raw Branch-Jump Note` 之后、`Open Checks` 之前。
- 由于文档类是 `article`，用户所说的“新开一章”落实为一个新的顶层
  `\section{Optical Conductivity: Kubo Sum and Response Lanczos}`。
- 正文沿用现有 note 的英文风格；公式、文件名和数值约定与当前代码一致。
- 新增 `subcaption` 宏包，用于历史图和 benchmark 图的紧凑分组。

## 统一电磁耦合与算符定义

章节首先固定均匀笛卡尔矢势约定

\[
  H(A_x)=H(0)+A_xJ_x+\frac{A_x^2}{2}K_{xx}+O(A_x^3),
\]

\[
  J_x=\left.\frac{\partial H}{\partial A_x}\right|_{A_x=0},
  \qquad
  K_{xx}=\left.\frac{\partial^2 H}{\partial A_x^2}\right|_{A_x=0}.
\]

当前有限环面实现以

\[
  H_A(k)=H(k+A_x\hat x)
\]

作为 `H(Ax)`、`Jx` 和 `Kxx` 的唯一 source of truth。`Ax` 是每单位笛卡尔长度
的均匀矢势；对 `RectLat4x6()`，总边界扭角是 `Phi_x=4Ax`。代码使用
`Jx=dH/dAx` 的整体符号约定，并设置 `e=hbar=1`。

## 谱表示与旧显式 Kubo 求和

对基态 `|g>` 定义

\[
  Q=1-|g\rangle\langle g|,
  \qquad |f\rangle=QJ_x|g\rangle,
\]

以及

\[
  \Delta_n=E_n-E_g,
  \qquad w_n=|\langle n|J_x|g\rangle|^2,
  \qquad z=\omega+i\eta.
\]

完整本征态求和给出

\[
  G(z)=\langle f|[z-(H-E_g)]^{-1}|f\rangle
      =\sum_{n\ne g}\frac{w_n}{z-\Delta_n},
\]

\[
  \chi(z)=G(z)+G(-z)
  =\sum_{n\ne g}w_n
   \left(\frac{1}{z-\Delta_n}-\frac{1}{z+\Delta_n}\right),
\]

\[
  \chi_0=\chi(0)=2\operatorname{Re}G(0)
  =-2\sum_{n\ne g}\frac{w_n}{\Delta_n}.
\]

旧 `ED/shift_current.jl::compute_optical_conductivity` 中启用的纵向 Kubo 表达式
可以按激发态逐项整理为

\[
  \sigma_{xx}^{\mathrm{old}}(z)
  =\frac{2\pi i}{A_{\mathrm{old}}}
   \frac{\chi(z)-\chi_0}{z}.
\]

因此旧公式的结构对应当前定义中的 regular conductivity，而不是包含独立
diamagnetic term 的 total conductivity。章节将展示这一代数关系，同时明确：旧图
使用的电流算符和面积约定与当前 benchmark 不同，所以旧图只能作为历史结果，
不能直接比较绝对峰高。

## 当前 regular、Drude 与 total conductivity

定义

\[
  D_{xx}=\langle g|K_{xx}|g\rangle+\chi_0,
\]

\[
  \sigma_{xx}^{\mathrm{Drude}}(z)
  =\frac{2\pi i}{A}\frac{D_{xx}}{z},
\]

\[
  \sigma_{xx}^{\mathrm{reg}}(z)
  =\frac{2\pi i}{A}\frac{\chi(z)-\chi_0}{z},
\]

\[
  \sigma_{xx}^{\mathrm{total}}(z)
  =\sigma_{xx}^{\mathrm{Drude}}(z)
   +\sigma_{xx}^{\mathrm{reg}}(z)
  =\frac{2\pi i}{A}
   \frac{\langle K_{xx}\rangle+\chi(z)}{z}.
\]

章节用表格明确以下来源和物理含义：

- `Re sigma_reg`：有限频率吸收谱；有限 `eta` 下为展宽峰。
- `Im sigma_reg`：regular response 的反应性/色散部分，是实部的
  Kramers--Kronig partner。
- `Re sigma_Drude`：有限 `eta` 下展宽的零频 Drude 峰；`eta -> 0+` 时收缩为
  零频 delta contribution。
- `Im sigma_Drude`：Drude 的反应性 `1/omega` pole（有限 `eta` 时平滑）。
- `Re/Im sigma_total`：相应 regular 与 Drude 分量之和。

章节还必须说明：有限环面、固定 twist 和固定 momentum branch 上计算的 `Dxx`
是该有限尺寸能量分支的曲率组合。当前 benchmark 的负 `Dxx` 不应直接解释成
热力学极限中的“负耗散 Drude weight”；后者需要跟踪全局基态、扭角依赖并做
尺寸外推。

当前 DAT 的列映射也写入正文：exact total、Lanczos total、exact regular 和
Lanczos regular 各自都有实部和虚部。现有 PDF 只选择绘制
`Re sigma_reg` 与 `Im sigma_total`，不代表其余分量没有计算。

## Exact 与 response Lanczos 的关系

章节明确说明两者使用同一个 `H,Jx,Kxx` 和同一个 Kubo 定义：

- `exact` 在有限 `m=2` sector 中完整对角化 `H`，保存所有 `Delta_n,w_n`，直接
  计算谱求和；它是有限尺寸、固定 sector、有限 broadening 的 exact reference，
  不是热力学极限 exact。
- response Lanczos 从 `|f>=QJx|g>` 出发，把 `H-Eg` 三对角化，只保存
  `alpha_j,beta_j`，用 continued fraction 计算同一个 `G(z)`。
- 当 Krylov 空间覆盖 response subspace 时，Lanczos 与 exact 只差浮点误差；
  当前 `M=892` benchmark 的 scaled maximum error 是
  `5.5927404437300286e-14`。

## 历史线性响应图片

从 `/Users/shajianyu/CMP_manybody/ED/组会` 复制以下六张 PDF 到新的、由 note
自包含管理的目录 `dmrg/report/figures/optical_response/`：

1. `optical_conductivity.pdf`
2. `optical_conductivity_0.001_-0.359.pdf`
3. `optical_conductivity_0.065.pdf`
4. `optical_conductivity_0.065_-0.373.pdf`
5. `long_x.pdf`
6. `longitudinal_y.pdf`

目标文件名保留原始 stem，只把不利于 LaTeX 路径的字符机械规范化。图片按自身
title/axis 分为 historical transverse/general 与 historical longitudinal 两组。

图片 caption 只写从图片或文件名能够确认的信息。文件名中的 `0.001`、`0.065`
标为 recorded broadening；`-0.359`、`-0.373` 在没有可靠 provenance 时只作为
原文件标签保留，不猜测其物理含义。

历史图前必须有醒目的解释框或段落，列出已确认的限制：

- 旧代码启用 `v_eigen_x=v_eigen_1`，而笛卡尔组合被注释；
- 旧 velocity operator construction 复用了会加入 interaction 的 Hamiltonian
  apply path；
- 旧面积因子使用 `L1*L2`，当前物理面积是 `12sqrt(3)`；
- 当前修正前后的零场 Hamiltonian 相同，但折叠 hopping 的一、二阶 `Ax`
  导数不同。

因此历史图用于记录旧计算形状和开发过程，不作为当前物理算符的 validation。

## 当前 `4x6` FCI benchmark 图片与数值

从 optical-response worktree 复制以下三张 PDF 到同一 figures 目录：

1. `optical_conductivity_regular.pdf`
2. `optical_conductivity_total.pdf`
3. `lanczos_convergence.pdf`

正文记录：

```text
geometry/sector          4x6 torus, Np=4, m=2
t1,t3,V1,V2,V3           1.0,0.2,1.0,0.0,0.0
omega grid               0:0.001:10
eta                      0.065
E0                       -9.313575887962
<Kxx>                    4.47923524291
||QJxg||^2               10.8443081433
Dxx                      -0.885036768331
regular peak             omega=4.107, Re sigma_reg=9.34982602717
final response dimension 892
exact/Lanczos error       5.5927404437300286e-14
```

图片 caption 解释 exact 蓝线被最后绘制的 Lanczos 曲线覆盖，而不是 exact 未画。

## 章节结构

新增 section 包含以下 subsections：

1. `Electromagnetic convention and operators`
2. `Explicit sum-over-states Kubo formula`
3. `Regular, Drude, and total conductivity`
4. `Exact diagonalization versus response Lanczos`
5. `Historical linear-response plots`
6. `Corrected 4x6 FCI benchmark`
7. `What can and cannot be compared`

最后一个 subsection 用表格对照旧结果与当前结果的 current operator、area、
response content、numerical method 和 validation status。

## 验证与交付

- 检查 copied PDFs 全部存在且页数为一。
- 编译 `dmrg_summary.tex`，解决所有 fatal LaTeX errors。
- 至少运行两遍 LaTeX 或使用 `latexmk`，保证引用和目录稳定。
- 使用 Poppler 把最终 PDF 的所有页面渲染为 PNG，逐页检查：公式不越界、表格
  不溢出、图片和 caption 不重叠、字体可读、没有空白或裁切页面。
- 检查 Git diff 只包含本章、必要宏包和受控 figures；保留 DMRG worktree 中现有
  unrelated untracked files，不修改或删除它们。

## 非目标

- 不实现或讨论 shift current。
- 不重新运行旧 ED 脚本。
- 不把旧图重新解释为当前正确算符的结果。
- 不修改 DMRG、ED 或 Lanczos production code。
- 不合并 DMRG branch 与 optical-response branch。
