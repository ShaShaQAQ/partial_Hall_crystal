# 30-site FQAHC/PHC 能谱前置展示设计

## 目标

在 `dmrg/report/dmrg_summary.tex` 的“30-site FQAHC/PHC 正式响应计算”之前，新增当前光学参数点的大体系 ED 能谱及数值说明，使读者先看到基态流形与同 sector 激发，再阅读光电导结果。

## 数据与参数

能谱数据来自：

`/Users/shajianyu/CMP_manybody/PHC_Multinode/PHC_Multinode/case3_30sites_Np12/output/spectrum_Np12.dat`

对应参数为

\[
N_s=30,\quad N_{\rm uc}=15,\quad N_p=12,\quad
t_1=1,\quad t_3=0.2,\quad V_1=100,\quad V_2=V_3=0.
\]

数据包含动量 sector $k=0,\ldots,14$，每个 sector 保存最低 8 个本征值，纵轴使用相对能量 $E-E_0$。

## 图片设计

重新生成一张双面板图片，不直接复用旧 PDF：

- 左面板绘制 15 个 sector 中全部 120 个已保存低能本征值，用于显示整体低能谱结构。
- 右面板放大 $E-E_0\lesssim0.30$ 区域，突出每个 sector 的最低态以及第 16 个能级。
- 最低 15 态使用醒目颜色；其余态使用中性颜色。
- 图内或图注明确列出系统尺寸、粒子数、hopping 与相互作用参数。

图片保存到 `dmrg/report/figures/optical_response/phc30_manybody_spectrum.png`，绘图脚本保存到 `dmrg/report/plot_phc30_spectrum.py`，输入数据复制到 `dmrg/report/data/optical_response/phc30_spectrum_V1_100.dat`，保证 note 可以独立复现。

## Note 内容与位置

在现有修正后的 $4\times6$ FCI benchmark 之后、30-site 光学响应之前新增小节“当前参数下的 30-site ED 能谱”。正文说明：

- 最低流形实际包含 15 个而不是 12 个近简并态；$N_p=12$ 是粒子数。
- 最低 15 态覆盖全部 $k=0,\ldots,14$，每个 sector 恰好一个。
- 流形内部散布为
  \[
  E_{15}-E_1=0.0287122.
  \]
- 第 16 态位于
  \[
  E_{16}-E_1=0.215794,
  \]
  因而流形到更高态的有限尺寸分离为
  \[
  E_{16}-E_{15}=0.187082.
  \]
- $m=5$ 与 $m=10$ 是最低的两个 sector；当前光电导分别从这两个基态计算。
- 均匀 $J_x$ 保持总动量，因此后面的光谱只直接探测同 sector 且具有非零电流矩阵元的激发。

## 验证标准

- 绘图脚本从保存的数据文件独立生成双面板图片。
- 图片包含 15 个 sector，每个 sector 8 个点；低能放大图能清楚识别 15 态流形。
- LaTeX 小节位于 30-site 光学响应之前，参数与数据文件 header 一致。
- XeLaTeX 编译成功，最终 PDF 无缺图、乱码、截断或重叠。
- 逐页渲染并检查新增能谱页及其相邻光学响应页。
