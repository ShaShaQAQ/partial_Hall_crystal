# 有限环面上物理正确的光学算符设计

日期：2026-08-26

## 目标

修正 `4×6, Np=4, m=2` 光电导 benchmark 中的电磁耦合层。保留已经验证的
response-Lanczos、continued-fraction resolvent 和 Kubo 后处理，只替换错误的
`H(Ax)`、`Jx`、`Kxx` 构造，使它们严格来自同一个 Bloch 哈密顿量族

\[
H_A(k)=H(k+A_x\hat x).
\]

所有 Julia 测试、完整对角化、benchmark 和画图继续只在 `W003` 运行；Mac
只进行源代码编辑、Git、rsync 和轻量静态检查。

## 已确认的根因

当前实现先在有限离散动量网格上把 `H(k)` Fourier 变换为折叠后的 `t(R)`，
然后为每个折叠后的有向格点对保留一个位移，并写成

\[
t_{ij}(A_x)=t_{ij}(0)e^{iA_xd^x_{ij}}.
\]

这个次序在 hopping range 小于环面所有半周长时可行，但不适用于当前
`Lx=4`、含第三近邻 `dx=±2=±Lx/2` 的模型。两条未折叠路径在零场下已经被
有限 Fourier 变换合并到同一个矩阵元。正确的有限场矩阵元是

\[
t_+e^{2iA_x}+t_-e^{-2iA_x},
\]

而当前实现得到

\[
(t_++t_-)e^{2iA_x}.
\]

因此当前代码虽然在 `A_x=0` 给出正确的哈密顿量和能谱，却给出错误的
一、二阶导数。现有有限差分测试只对错误的 `build_hops_Ax` 自身求差分，
不能检验物理耦合。

在 `W003` 上以 `H(k+A_x)` 为独立参考，已经复现：

```text
E0                              -9.3135758879616
relative error of current Jx     1.100677899174417
relative error of Kxx            2.1589512290624735
old ||Q Jx g||^2                 3.3968963658100124
physical ||Q Jx g||^2           10.844308142783673
old regular peak                 omega=4.323, Re sigma=1.2604123159
physical regular peak            omega=4.107, Re sigma=9.3498260266
```

## 物理约定

固定 Bloch/embedding 规范不变：

- `get_Hk(k, t1, t3)` 仍是模型的 Bloch 哈密顿量；
- `SUBLAT_POS` 仍定义轨道 embedding；
- 均匀笛卡尔矢势通过 `k -> k + Ax * xhat` 引入；
- `Ax` 的单位是每单位笛卡尔长度的均匀矢势，不是总边界扭角。对
  `RectLat4x6()` 的水平超胞周期 `Lx=4`，两者的关系是
  `Phi_x = 4Ax`；
- 本代码以 `Jx = dH/dAx` 作为算符符号约定。若把电子电荷另外写成
  `q=-e`，物理电流的整体负号应统一放在外部电荷约定中，不在此处
  混入 Fourier 变换；
- 单位保持 `e = hbar = 1`；
- 系统面积仍由真实超胞矢量计算。

定义

\[
H^{(0)}(k)=H(k),\qquad
H^{(1)}_x(k)=\partial_{k_x}H(k),\qquad
H^{(2)}_{xx}(k)=\partial_{k_x}^2H(k).
\]

在有限离散 `k` 网格上分别进行同一个带 embedding 的逆 Fourier 变换：

\[
t_R^{(n)}=\frac1{N_k}\sum_k
e^{-ik\cdot(R+\delta_\alpha-\delta_\beta)}H_x^{(n)}(k).
\]

最后才把每个 `t_R^(n)` 折叠到有限环面并构造 many-body hopping list：

\[
H=\sum t_R^{(0)}c^\dagger c+H_{int},\qquad
J_x=\sum t_R^{(1)}c^\dagger c,\qquad
K_{xx}=\sum t_R^{(2)}c^\dagger c.
\]

相互作用不依赖 `A_x`，所以只进入 `H`，不进入 `J_x` 或 `K_{xx}`。

## 解析导数

对当前模型的三个向量 `a_j`，

\[
g_j(k)=2t_1\cos(k\cdot a_j),\qquad
g_0(k)=2t_3\sum_j\cos(2k\cdot a_j).
\]

使用解析导数

\[
\partial_{k_x}g_j=-2t_1a_{j,x}\sin(k\cdot a_j),
\]

\[
\partial_{k_x}^2g_j=-2t_1a_{j,x}^2\cos(k\cdot a_j),
\]

\[
\partial_{k_x}g_0=-4t_3\sum_j a_{j,x}\sin(2k\cdot a_j),
\]

\[
\partial_{k_x}^2g_0=-8t_3\sum_j a_{j,x}^2\cos(2k\cdot a_j).
\]

解析导数避免用数值差分生成生产算符；数值差分只作为独立测试 oracle。

## 软件边界

### `shared/hoppings.jl`

新增或重构以下职责：

1. `get_Hk_x_derivatives(k, t1, t3)`：返回 Bloch 层的 `Hk`、`dHdkx`、
   `d2Hdkx2`。
2. `fourier_to_real_Ax(lat, t1, t3, Ax)`：直接 Fourier 变换
   `get_Hk(k + [Ax,0])`，作为物理 `H(Ax)`。
3. `fourier_to_real_x_derivatives(lat, t1, t3)`：分别 Fourier 变换解析
   `Hk`、`dHdkx`、`d2Hdkx2`。
4. 一个只负责把给定 `tR` 折叠成 hopping list 的 helper。它对折叠到
   同一 `(target, source)` 的贡献求和，而不是保留遇到的第一条；也不猜测
   单一未折叠位移，或用手工添加的共轭反向键替代 Fourier 变换本身给出的
   有向矩阵元。不同未折叠路径的信息必须在 Fourier 变换前由
   `H(k+Ax)` 保留。
5. `build_hops_Ax` 改为调用 `fourier_to_real_Ax`。
6. `build_hops_x_derivatives` 改为调用解析 Bloch 导数，并继续返回
   `(hops, jx_hops, kxx_hops)`，从而不改变 benchmark 调用接口。

现有 `build_hops(lat, ..., phi_y)` 暂不扩展为完整的笛卡尔 `Ay` 光学接口；
本次只保证其 `phi_y=0` 行为和旧能谱不回归。`RealSpaceBond` 若仍被谱流代码
使用可以保留，但不得再作为 `Ax` 光学算符的 source of truth。

### `shared/optical_response.jl`

不修改 response source、Lanczos 递推、continued fraction、exact spectral
reference 或 Kubo 分解。它们已经对任意给定的 Hermitian `H,J,K` 通过直接
矩阵验证。

### `case4_4x6_Np4/benchmark_optical_response.jl`

继续使用 `build_hops_x_derivatives`，但输出元数据增加算符构造约定，避免旧
JLD2 与新结果混淆。旧 `optical_response_output` 必须在 `W003` 上重新生成并
同步回本地。

## TDD 与独立物理验证

### RED 1：Bloch 解析导数

对多个非对称 `k` 点比较解析 `dH/dkx`、`d2H/dkx2` 与
`get_Hk(k±delta*xhat)` 的中心差分。该测试在新 API 尚不存在时先失败。

### RED 2：有限环面单粒子导数

把 hopping list 累加为 `Ns×Ns` 单粒子矩阵，检查

\[
J_x=\frac{H(+\delta A)-H(-\delta A)}{2\delta A},
\]

\[
K_{xx}=\frac{H(+\delta A)-2H(0)+H(-\delta A)}{\delta A^2}.
\]

测试必须使用 `RectLat4x6()`、`t3=0.2`，从而覆盖 `dx=±Lx/2` 的回归情形。
该测试针对 `H(k+A_x)`，不能再以单一位移 Peierls list 作为 oracle。

### RED 3：many-body sector 导数

在 `Np=4, m=2` 中物化 `H(Ax)`，直接比较 many-body `Jx,Kxx` 与有限差分。
同时验证：

```text
E0                         -9.3135758879616
||Q Jx g||^2               10.8443081428  (target tolerance 1e-8)
<Kxx>                       4.47922837      (finite-difference tolerance 1e-5)
```

### GREEN：完整响应 benchmark

修正算符后重新运行 exact sum-over-states 和 response-Lanczos。验收条件：

- `H,Jx,Kxx` Hermiticity error `< 1e-12`；
- `H(0)` 与修正前零场 Hamiltonian relative error `< 1e-12`；
- many-body 一阶有限差分 relative error `< 1e-8`；
- many-body 二阶有限差分 relative error `< 1e-5`；
- `abs(<g|QJx|g>) < 1e-11`；
- exact 与最终 Lanczos regular curve scaled error `< 1e-8`；
- 使用同一个物理 `Jx` 时，旧显式求和公式与新 `regular` 公式逐点一致；
- `omega=0:0.001:10` 的新 JLD2、DAT 和 PDF 全部重新生成。

`Drude weight` 不预设符号或强制为零；只要求 exact 与 Lanczos 一致，并明确
输出物理 `Kxx` 所得到的值。

## 非目标

本次不实现：

- `Jy`、`Kyy`、`sigma_xy`；
- many-body shift current；
- 30-site 或多节点生产版本；
- 旧 `/Users/shajianyu/CMP_manybody/ED` 目录的修复；
- 一般自动微分框架。

这些工作只有在本次 `4×6` 物理算符 benchmark 通过后再分别设计。
