# 一般光电导的 response-Lanczos benchmark 设计

日期：2026-08-26

## 目标

在现有的 `4×6`、`Np=4` 小系统 ED case 上，把纵向光电导从显式激发态求和改写成从

\[
|f_x\rangle=QJ_x|g\rangle
\]

出发的 Lanczos/resolvent 计算，并用完整本征谱逐点验证。第一阶段只证明这个数值方法和一般 Kubo 公式正确，不实现 shift current，也不改多节点程序。

实现放在 Git 仓库 `partial_Hall_crystal` 中，复用 `case4_4x6_Np4` 与 `shared` 里的模块化 ED 代码。旧目录 `/Users/shajianyu/CMP_manybody/ED/shift_current.jl` 仅用于核对历史参数和公式，不作修改。

## Benchmark case

- 晶格：`RectLat4x6()`，24 sites、12 unit cells
- 粒子数：`Np=4`
- 参数：`t1=1.0`、`t3=0.2`、`V1=1.0`、`V2=V3=0`
- 旧代码动量：`kpoint=(2,0)`；在新代码中对应 sector label `m=2`
- 默认展宽：`eta=0.065`
- 默认频率范围：`omega=0:0.001:10`
- 单位：`e=hbar=1`
- 系统面积：由真实超胞矢量计算，不使用裸的 `L1*L2`

## 电磁耦合和算符定义

统一把外参量定义为均匀矢势 `A_x`：

\[
t_{ij}(A_x)=t_{ij}\exp(iA_xd^x_{ij}).
\]

由同一个 Peierls 耦合生成

\[
J_x=\left.\frac{\partial H}{\partial A_x}\right|_{A_x=0}
=\sum_{ij}i d^x_{ij}t_{ij}c_i^\dagger c_j,
\]

\[
K_{xx}=\left.\frac{\partial^2 H}{\partial A_x^2}\right|_{A_x=0}
=-\sum_{ij}(d^x_{ij})^2t_{ij}c_i^\dagger c_j.
\]

`d_ij` 必须来自 hopping 构造时的未折叠物理位移，不能从周期边界条件规约后的 site index 推回。密度-密度相互作用不依赖 `A_x`，因此不进入 `J_x` 和 `K_xx`。

旧代码把 boundary flux `Phi_x` 和均匀矢势 `A_x=Phi_x/L_x` 混用了；新 benchmark 不继承这个归一化。显式求和和 Lanczos 两条路径必须使用完全相同的新 `J_x`、`K_xx`，从而把“算法是否正确”和“历史公式是否正确”分开。

## 一般 Kubo 结果

对选定 sector 的归一化基态 `|g>` 和基态能量 `E_g`，定义

\[
Q=1-|g\rangle\langle g|,
\qquad
|f_x\rangle=QJ_x|g\rangle,
\qquad
\bar H=H-E_g.
\]

响应 resolvent 为

\[
G_{xx}(z)=\langle f_x|(z-\bar H)^{-1}|f_x\rangle,
\qquad z=\omega+i\eta.
\]

一般纵向电导定义为

\[
\sigma_{xx}(z)=\frac{2\pi i}{A_{\rm sys}z}
\left[\langle g|K_{xx}|g\rangle+G_{xx}(z)+G_{xx}(-z)\right].
\]

不再用激发态求和强制 Drude weight 为零。单独输出

\[
D_{xx}=\langle K_{xx}\rangle+G_{xx}(0)+G_{xx}(-0)
=\langle K_{xx}\rangle
-2\sum_{n\ne g}\frac{|\langle n|J_x|g\rangle|^2}{E_n-E_g}.
\]

并将总响应拆成

\[
\sigma_{xx}^{\rm Drude}(z)=\frac{2\pi iD_{xx}}{A_{\rm sys}z},
\]

\[
\sigma_{xx}^{\rm reg}(z)=\frac{2\pi i}{A_{\rm sys}z}
\left[G_{xx}(z)+G_{xx}(-z)-G_{xx}(0)-G_{xx}(-0)\right].
\]

## 软件结构

新增两个入口，不改动旧 `shift_current.jl`：

1. `shared/optical_response.jl`
   - 从 hopping 的真实位移构造 `H`、`J_x`、`K_xx` 所需的 hopping lists；
   - 计算并投影 `J_x|g>`；
   - 实现三项 response Lanczos；
   - 从三对角系数计算 continued-fraction resolvent；
   - 提供一般 Kubo、Drude 和 regular 部分的纯数值函数。
2. `case4_4x6_Np4/benchmark_optical_response.jl`
   - 构造已确认的 benchmark case；
   - 在 `m=2` sector 中运行完整对角化和 response Lanczos；
   - 扫描 Lanczos 步数并输出误差、数据和对比图。

如需扩展现有 hopping 模块，只新增保留物理位移的构造接口；现有 `build_hops` 的返回值和已有 case 行为保持不变。

## 两条独立计算路径

### 完整本征态参考

对小 sector 完整对角化，直接计算

\[
G_{xx}^{\rm exact}(z)=
\sum_{n\ne g}
\frac{|\langle n|J_x|g\rangle|^2}
{z-(E_n-E_g)}.
\]

同一组本征态还用于计算 exact `D_xx` 和一般 Kubo 曲线。这条路径仅是小系统测试参考，不进入未来 30-site production 路径。

### Response Lanczos

以 `q1=f_x/norm(f_x)` 为起点，对 `H-E_g` 执行三项递推，保存

\[
\{\alpha_j\}_{j=1}^{M},\qquad
\{\beta_j\}_{j=1}^{M-1}.
\]

构造三对角矩阵 `T_M` 后计算

\[
G_{xx}^{(M)}(z)=\|f_x\|^2[(zI-T_M)^{-1}]_{11}.
\]

频率和展宽只进入小矩阵后处理，不触发新的 many-body Lanczos。benchmark 扫描 `M=50,100,200,400`，并在需要时继续到 sector 维数或 Krylov breakdown。

## 数值检查和错误处理

运行响应前必须检查：

- `H`、`J_x`、`K_xx` 的 Hermiticity；
- 基态归一化与 `norm(H*g-E_g*g)`；
- `abs(imag(dot(g,K_xx*g)))` 足够小；
- `dot(g,f_x)` 在投影后接近零；
- `norm(f_x)` 非零，否则报告该 sector 没有 `x` 方向光学权重；
- Lanczos 的 `alpha` 实、`beta` 非负；
- Krylov breakdown 被当作正常终止并记录，不产生除零。

`omega=0` 不直接用总电导的 `1/z` 形式解释；代码分别报告 `D_xx` 和有限 `eta` 的频谱。

## 测试和验收

### 单元测试

使用一个小型随机 Hermitian 矩阵和 Hermitian `J`、`K`：

- 在 Krylov 空间闭合时，continued fraction 与直接矩阵逆一致；
- exact sum-over-states 与 resolvent 一致；
- Drude/regular 分解之和等于 total；
- ground-state projection 去掉零能极点。

### 物理 case 集成测试

在 `4×6, Np=4, m=2` 上比较：

- `E_g`；
- `norm(f_x)^2`；
- `<K_xx>`；
- `D_xx`；
- `Re/Im G_xx(omega+i eta)`；
- `Re/Im sigma_xx^reg` 和 total curve。

使用尺度化最大误差

\[
\epsilon_X=
\frac{\max_\omega|X_{\rm Lanczos}(\omega)-X_{\rm exact}(\omega)|}
{\max(1,\max_\omega|X_{\rm exact}(\omega)|)}.
\]

必须观察到误差随 `M` 增加而下降；最终 `M` 的目标误差为 `1e-8`。若三项递推的有限精度使该阈值不可达，应先定位正交性问题，不能通过放宽阈值掩盖。

## 输出

benchmark 输出到独立子目录，避免覆盖现有能谱和结构因子文件：

- JLD2：模型参数、`E_g`、`K_xx`、`D_xx`、`norm(f_x)^2`、`alpha`、`beta`、exact/Lanczos 曲线和误差；
- 文本数据：`omega`、exact/Lanczos 的 real/imag regular/total conductivity；
- PDF：exact 与不同 `M` 的叠图、误差随 `M` 的收敛图；
- 终端摘要：算符检查、残差、Drude weight、每个 `M` 的误差和耗时。

## 非目标

本阶段不包含：

- `sigma_yy` 或 `sigma_xy`；
- many-body shift current；
- ground-state multiplet 平均；
- 30-site 计算；
- Slurm、多节点调度或 production 内存优化。

只有小系统 benchmark 通过后，才为这些内容分别设计后续阶段。
