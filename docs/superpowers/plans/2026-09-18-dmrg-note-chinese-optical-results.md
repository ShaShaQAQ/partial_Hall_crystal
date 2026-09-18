# DMRG 结果汇总 Note 中文化实施计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 在不改变现有章节和结果顺序的前提下，将 DMRG 结果汇总 note 完整改写为中文，并补全旧 FCI 与 30-site FQAHC/PHC 光电导图的计算方法和直接结果说明。

**Architecture:** 只修改单一 LaTeX 源文件中的语言与说明文字，保留公式、数值、标签和图片资源。使用 XeLaTeX/latexmk 生成 PDF，再用 Poppler 将所有页面渲染为 PNG，完成文本检查与视觉检查。

**Tech Stack:** LaTeX、ctex、latexmk、XeLaTeX、Poppler (`pdftoppm`)、ripgrep。

---

## 文件结构

- 修改：`dmrg/report/dmrg_summary.tex`——中文结果汇总、公式、表格和全部图注。
- 生成：`dmrg/report/dmrg_summary.pdf`——最终可阅读文档。
- 读取但不修改：`dmrg/report/figures/optical_response/*`——历史 FCI、修正 benchmark 与 30-site PHC 图片。
- 临时生成：`tmp/pdfs/dmrg_summary-*.png`——逐页视觉检查文件，检查后可删除。

### Task 1：启用中文排版并翻译 DMRG/flux 章节

**Files:**
- Modify: `dmrg/report/dmrg_summary.tex:1-309`

- [ ] **Step 1：检查中文编译环境**

Run:

```bash
command -v latexmk
command -v xelatex
kpsewhich ctex.sty
```

Expected: 三条命令均返回可执行文件或 `ctex.sty` 的路径。

- [ ] **Step 2：启用中文并翻译文档元信息**

将导言区和标题改为：

```tex
\documentclass[11pt]{article}
\usepackage[UTF8,scheme=plain]{ctex}
...
\title{DMRG 圆柱基准与电荷泵浦结果汇总}
\author{Partial Hall Crystal Project}
\date{\today}
```

- [ ] **Step 3：逐段翻译前六个章节**

保持所有数值和公式不变，将章节标题固定为：

```tex
\section{文档目的}
\section{模型与几何}
\section{圆柱 ED 与 DMRG 基准}
\section{Lx=15 warm-start 谱流计算}
\section{电荷泵浦观测量}
\section{通量插入结果}
\section{原始分支跳变说明}
```

表头、图例、图注和正文全部翻译为中文；路径、运行目录、参数名和 `converged=true` 保持原样。

- [ ] **Step 4：检查前半部分是否仍残留英文叙述**

Run:

```bash
sed -n '1,309p' dmrg/report/dmrg_summary.tex | rg -n '^(This|The|For|Thus|No |Supporting|Generated|Only )|\\caption\{[A-Za-z]'
```

Expected: 除专有名词、代码路径和必要英文术语外，不再出现完整英文叙述或英文图注。

### Task 2：翻译光电导公式和算法说明

**Files:**
- Modify: `dmrg/report/dmrg_summary.tex:310-487`

- [ ] **Step 1：保留公式并翻译定义说明**

保留 Eqs. `optical-operators` 至 `real-imag-source`，将章节层级改为：

```tex
\section{光电导：Kubo 求和与 response-Lanczos}
\subsection{电磁场约定与算符}
\subsection{显式态求和的 Kubo 公式}
\subsection{regular、Drude 与 total 电导}
\subsection{完全对角化与 response-Lanczos}
```

逐项说明 regular 实部是有限频率吸收，虚部是色散响应；total 是 regular 与 Drude 的复数和。

- [ ] **Step 2：明确两种方法都属于 ED**

正文必须包含以下含义明确的表述：

```tex
这里的 exact sum-over-states 与 response-Lanczos 都是 ED 方法，而不是
DMRG 动力学响应。前者显式获得相关激发本征态并求和；后者从
$QJ_x|g\rangle$ 出发，在 Krylov 空间内计算同一个 resolvent，因此无需求出
或保存全部激发态，适合 Hilbert 空间很大的 ED sector。
```

- [ ] **Step 3：检查方法归属和关键术语**

Run:

```bash
rg -n '都.*ED|不是.*DMRG|QJ_x|Krylov|全部激发态|sum-over-states' dmrg/report/dmrg_summary.tex
```

Expected: 能定位到同一小节中的 ED/DMRG 区分、源向量和无需全部激发态的说明。

### Task 3：翻译并补全每组光学图片的结果说明

**Files:**
- Modify: `dmrg/report/dmrg_summary.tex:488-692`

- [ ] **Step 1：翻译历史光学结果说明**

保留历史图，并在图前说明它们是旧算符与旧面积约定下的 ED Kubo 结果，只作历史记录；不得把未解析文件标签赋予物理意义。

- [ ] **Step 2：补全修正后的 4x6 FCI benchmark**

正文和图注必须写清：

```tex
蓝色 exact 曲线由有限动量 sector 的完整 ED 本征态求和得到；红色虚线是
从 $QJ_x|g\rangle$ 出发的 response-Lanczos 结果。两者最大缩放误差为
$5.59\times10^{-14}$，所以在主图分辨率下重合；蓝线并非没有绘制，而是
被重合曲线遮住。
```

同时记录主吸收峰 `\(\omega=4.107\)` 和 `\(\operatorname{Re}\sigma_{xx}^{\rm reg}=9.349826\)`。

- [ ] **Step 3：补全 30-site FQAHC/PHC 四幅图**

保留四幅 regular/total 实虚部图片，说明基态由多节点 ED 获得、响应由 response-Lanczos ED 计算。记录：

```tex
N_s=30,\quad N_p=12,\quad V_1=100,\quad \eta=0.065,\quad M=600,
```

并写清 `m=5` 与 `m=10` 曲线最大差约 `\(3.22\times10^{-7}\)`，它们是对称相关 sector 的同一物理响应，而不是 exact/Lanczos 对照。

- [ ] **Step 4：记录直接光谱结果及解读边界**

记录主要峰：

```tex
(\omega,\operatorname{Re}\sigma_{xx}^{\rm reg})
\simeq(1.524,4.55784),\qquad(2.248,1.33201),
```

将它们称为 `\(J_x\)` 可耦合的 current-active 中性激发。说明其他动量 sector 的低能态不会被均匀 `\(J_x\)` 直接看到，因此 optical gap 不等于完整 many-body spectral gap。

- [ ] **Step 5：翻译比较表与未完成检查**

翻译原有“可比较/不可比较”和 open checks；明确历史图与修正图的峰高、低频行为、tensor 标签不能直接定量比较。

### Task 4：编译、自动检查和 PDF 视觉检查

**Files:**
- Modify if needed: `dmrg/report/dmrg_summary.tex`
- Generate: `dmrg/report/dmrg_summary.pdf`
- Generate temporarily: `tmp/pdfs/dmrg_summary-*.png`

- [ ] **Step 1：以 XeLaTeX 编译两遍**

Run:

```bash
cd dmrg/report
latexmk -xelatex -interaction=nonstopmode -halt-on-error dmrg_summary.tex
```

Expected: exit code 0，生成 `dmrg_summary.pdf`。

- [ ] **Step 2：检查日志、引用与四幅新图**

Run:

```bash
rg -n 'LaTeX Warning|Undefined control sequence|Missing character|not found|multiply defined' dmrg/report/dmrg_summary.log
rg -n 'phc30_sigma_xx_(regular|total)_(real|imag)\.png' dmrg/report/dmrg_summary.tex
```

Expected: 日志无未解析引用、缺字或缺图错误；第二条命令返回四个唯一图片引用。

- [ ] **Step 3：检查 PDF 文本中的关键表述**

Run:

```bash
pdftotext dmrg/report/dmrg_summary.pdf - | rg '都属于 ED|不是 DMRG|4.107|1.524|2.248|3.22'
```

Expected: 六类关键内容全部可从最终 PDF 提取。

- [ ] **Step 4：渲染全部页面**

Run:

```bash
mkdir -p tmp/pdfs
pdftoppm -png -r 130 dmrg/report/dmrg_summary.pdf tmp/pdfs/dmrg_summary
```

Expected: 每一页都生成对应 PNG，且页数与 `pdfinfo` 输出一致。

- [ ] **Step 5：逐页检查并修复排版问题**

检查所有页面的中文字体、页码、表格宽度、公式断行、图像锐度、图注、浮动位置和章节转换；若发现截断、重叠、乱码或孤立标题，修复 `dmrg_summary.tex` 后重新执行 Task 4。

- [ ] **Step 6：最终差异检查**

Run:

```bash
git diff --check -- dmrg/report/dmrg_summary.tex
git diff --stat -- dmrg/report/dmrg_summary.tex dmrg/report/dmrg_summary.pdf
```

Expected: `git diff --check` 无输出；差异只涉及已批准的 note 源文件和重新生成的 PDF。
