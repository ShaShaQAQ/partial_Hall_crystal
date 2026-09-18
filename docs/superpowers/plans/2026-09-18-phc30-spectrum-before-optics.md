# 30-site FQAHC/PHC 能谱前置展示实施计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 从当前 $V_1=100,V_2=V_3=0$ 的 30-site 多节点 ED 数据生成可复现的双面板能谱图，并在对应光学响应之前写清参数与 15 态近简并流形。

**Architecture:** 将外部计算结果复制为 note 自包含的数据文件；用一个小型 Python/Matplotlib 脚本负责解析、校验、统计和绘图；以测试固定 120 个数据点及三个低能流形指标；最后将图片与数值说明插入现有 LaTeX。

**Tech Stack:** Python 3、NumPy、Matplotlib、pytest、XeLaTeX、Poppler。

---

### Task 1：加入能谱数据和红灯测试

**Files:**
- Create: `dmrg/report/data/optical_response/phc30_spectrum_V1_100.dat`
- Create: `dmrg/report/test_plot_phc30_spectrum.py`

- [ ] **Step 1：复制当前参数点的完整能谱数据**

数据文件保留原 header：

```text
# k  E-E0  [30sites Np=12 V1=100.0 V2=0.0 V3=0.0 t'=0.2 multi-node]
```

并保留 $15\times8=120$ 行本征值。

- [ ] **Step 2：编写解析和低能统计测试**

测试导入 `load_spectrum` 与 `summarize_low_energy`，断言：

```python
assert len(rows) == 120
assert set(rows[:, 0].astype(int)) == set(range(15))
assert all(np.count_nonzero(rows[:, 0] == k) == 8 for k in range(15))
assert summary["manifold_size"] == 15
assert summary["manifold_width"] == pytest.approx(0.0287121891)
assert summary["sixteenth_energy"] == pytest.approx(0.2157939408)
assert summary["separation"] == pytest.approx(0.1870817517)
```

- [ ] **Step 3：运行测试并确认因绘图模块尚不存在而失败**

Run:

```bash
pytest -q dmrg/report/test_plot_phc30_spectrum.py
```

Expected: collection error，提示 `plot_phc30_spectrum` 无法导入。

### Task 2：实现绘图脚本并生成图片

**Files:**
- Create: `dmrg/report/plot_phc30_spectrum.py`
- Generate: `dmrg/report/figures/optical_response/phc30_manybody_spectrum.png`

- [ ] **Step 1：实现数据解析和统计函数**

脚本提供：

```python
def load_spectrum(path: Path) -> np.ndarray: ...
def summarize_low_energy(rows: np.ndarray, manifold_size: int = 15) -> dict[str, float]: ...
def make_figure(rows: np.ndarray, output: Path) -> None: ...
```

`load_spectrum` 检查两列数据、sector 为 $0\ldots14$ 且每个 sector 恰有 8 个能级；`summarize_low_energy` 对全局能量排序并返回 $E_{15}-E_1$、$E_{16}-E_1$ 和 $E_{16}-E_{15}$。

- [ ] **Step 2：实现双面板图**

左面板画全部 120 个点；右面板放大 $-0.005\le E-E_0\le0.30$。最低 15 态画为红色圆点，其余态画为蓝灰色小圆点；右图加虚线标记 $E_{15}$ 和 $E_{16}$。标题或图内文字列出

```text
Ns=30, Nuc=15, Np=12, t1=1, t3=0.2, V1=100, V2=V3=0
```

- [ ] **Step 3：运行测试并确认通过**

Run:

```bash
pytest -q dmrg/report/test_plot_phc30_spectrum.py
```

Expected: `2 passed` 或更多，无失败。

- [ ] **Step 4：生成最终图片**

Run:

```bash
python3 dmrg/report/plot_phc30_spectrum.py
```

Expected: 输出 `phc30_manybody_spectrum.png`，脚本打印 120 个点、15 个 sector 及三个低能指标。

### Task 3：插入 LaTeX 并验证 PDF

**Files:**
- Modify: `dmrg/report/dmrg_summary.tex`
- Generate: `dmrg/report/dmrg_summary.pdf`

- [ ] **Step 1：在 30-site 光学响应之前加入能谱小节**

小节标题为：

```tex
\subsection{当前参数下的 30-site ED 能谱}
```

正文列出完整参数、每 sector 8 个本征值、15 态流形三个指标，并明确 $N_p=12$ 是粒子数而非基态简并数。

- [ ] **Step 2：加入双面板能谱图片与光学选择定则连接**

图片引用：

```tex
\includegraphics[width=0.96\linewidth]{figures/optical_response/phc30_manybody_spectrum.png}
```

图后说明 $m=5,10$ 是最低 sector，均匀 $J_x$ 保持总动量，后续光电导只看到同 sector 的 current-active 激发。

- [ ] **Step 3：强制编译并检查日志**

Run:

```bash
cd dmrg/report
latexmk -g -xelatex -interaction=nonstopmode -halt-on-error dmrg_summary.tex
```

Expected: exit code 0，无缺图、缺字、未定义引用或未定义控制序列。

- [ ] **Step 4：渲染并检查新增页面**

Run:

```bash
mkdir -p tmp/pdfs/phc30_spectrum_note
pdftoppm -png -r 140 dmrg/report/dmrg_summary.pdf tmp/pdfs/phc30_spectrum_note/page
```

检查新增能谱图参数可读、120 个点未被裁切、右图能识别 15 态流形，并确认它位于 30-site 光学响应之前。

- [ ] **Step 5：最终自动检查**

Run:

```bash
pytest -q dmrg/report/test_plot_phc30_spectrum.py
git diff --check -- dmrg/report/plot_phc30_spectrum.py dmrg/report/test_plot_phc30_spectrum.py dmrg/report/dmrg_summary.tex
```

Expected: 测试全过且 `git diff --check` 无输出。
