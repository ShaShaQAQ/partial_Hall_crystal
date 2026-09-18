# DMRG 结果汇总 note 中文化与光学结果补充设计

## 目标

将 `dmrg/report/dmrg_summary.tex` 整体改写为中文，并把最新的 30-site FQAHC/PHC 光电导结果补入现有光学响应章节。本文档仍定位为计算结果汇总 note，而不是以完整证据链组织的论文草稿。

## 组织原则

- 保留现有章节顺序、结果出现顺序和图片分组，不重新组织物理论证结构。
- 标题、正文、图注、表格与总结均改为中文；公式、变量名、sector 标签、代码路径和数据文件名保持原样。
- 保留已有结果，不删除历史图或尚未完成的检查项。
- 不把不同尺寸、填充或相互作用参数下的结果当作连续参数扫描来比较。

## 逐图说明格式

每组结果说明以下信息：

1. 系统尺寸、粒子数、填充和相互作用参数；
2. 使用的数值方法和输入数据；
3. 横纵轴、谱函数及 regular/Drude/total 等约定；
4. 图中可直接读出的峰位、曲线重合关系和数值结论；
5. 必要时注明该图不能单独支持的物理结论。

## 光电导算法表述

- 传统 exact Kubo 曲线来自 ED 的 sum-over-states 表达，需要取得相关激发本征态并显式求和。
- response-Lanczos 同样属于 ED。它从投影后的源向量 $QJ_x|g\rangle$ 出发，在 Krylov 空间中计算 resolvent 或 continued fraction，因而无需求出并保存全部激发态。
- response-Lanczos 不是 DMRG 动力学响应方法。多节点程序用于大 Hilbert 空间中的基态 ED 以及后续 Krylov 迭代。
- 旧 $4\times6$ FCI benchmark 中 exact 与 Lanczos 曲线在绘图精度内重合，因此后画出的曲线会遮住先画出的曲线；图注与正文都应说明这一点。

## 新增 30-site FQAHC/PHC 结果

在现有光学响应章节后部加入四幅结果：

- $\operatorname{Re}\sigma_{xx}^{\rm reg}$；
- $\operatorname{Im}\sigma_{xx}^{\rm reg}$；
- $\operatorname{Re}\sigma_{xx}^{\rm total}$；
- $\operatorname{Im}\sigma_{xx}^{\rm total}$。

正文记录计算参数 $N_s=30$、$N_p=12$、$t_1=1$、$t_3=0.2$、$V_1=100$、$V_2=V_3=0$、$\eta=0.065$、$M=600$，并说明 $m=5$ 与 $m=10$ 曲线最大差约为 $3.2\times10^{-7}$，代表同一对称相关物理光谱。

标出主要 regular absorption peaks：$\omega\simeq1.524$ 与 $2.248$；较弱峰只作为数值结果列出。文字将这些峰称为 $J_x$ 可耦合的 current-active neutral excitations，不将其直接等同于完整多体能隙。

## 物理解读边界

- 旧 FCI 光谱用于读取 $q=0$ bright neutral excitation、optical gap、振子强度，并验证 exact Kubo 与 response-Lanczos 的一致性。
- 当前 FQAHC/PHC 光谱用于读取晶体/折叠相中哪些同 sector 激发具有非零电流矩阵元，以及 bright/dark mode 和谱权重分布。
- 光电导不能单独确定拓扑基态简并、many-body Chern number、分数 Hall 量子化、任意子性质、CDW 波矢或 phason。
- 当前结构因子数据对应 $V_1=10,V_2=V_3=2$，而光电导对应 $V_1=100,V_2=V_3=0$；note 中不得用前者直接给后者的光谱峰作归属。

## 验证标准

- LaTeX 能以支持中文的引擎完整编译，无未解析引用或缺图错误。
- 输出 PDF 中中文字体正常，公式、表格、图片与图注不存在截断、重叠或乱码。
- 四幅 30-site 光电导图片全部出现，并各自具有计算方法和结果说明。
- 文中所有光学结果均明确标记为 ED；response-Lanczos 的作用和 exact sum-over-states 的区别准确。
- 逐页渲染 PDF 并完成视觉检查。
