# ED 激发计算整合实施计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 将多节点能谱、谱流、静态结构因子、response-Lanczos 和跨方法结果 note 整合到一个可在 W003 与 GitHub 同步的仓库结构中。

**Architecture:** 以 `feature/spectral-flow` 为集成基线，合并 optical-response 和 DMRG 历史；用 `ed/`、`dmrg/`、`reports/` 和 `results/` 表达长期职责。大波函数只存 W003，Git 中用 TOML 清单和 SHA256 建立可追溯引用。

**Tech Stack:** Git worktree、Julia、JLD2、TOML、Python/pytest、XeLaTeX、W003 SSH。

---

### Task 1：建立集成历史

**Files:** Git history only.

- [ ] 将 `feature/optical-response-lanczos` 合入 `integration/ed-excitations`。
- [ ] 解决公共 ED 文件冲突，并运行受影响的小体系测试。
- [ ] 将 `DMRG` 已提交历史合入，保留 `dmrg/` 模块和报告资产。
- [ ] 提交合并结果，不纳入其他工作树未提交文件。

### Task 2：先写目录与数据策略测试

**Files:**
- Create: `ed/tests/test_repository_layout.py`
- Modify: `.gitignore`

- [ ] 编写失败测试，要求 `ed/shared`、`ed/cases`、`ed/projected`、
  `reports/fqahc_results` 和 `results` 存在。
- [ ] 编写失败测试，扫描 Julia 文件并拒绝指向旧根目录 `../shared` 的 include。
- [ ] 编写失败测试，要求 Git ignore `*.jld2`、集群日志和 runtime output。
- [ ] 运行 pytest 并确认因目录尚未迁移而失败。

### Task 3：迁移 ED 与报告目录

**Files:**
- Move: `shared/` -> `ed/shared/`
- Move: `case*/` -> `ed/cases/case*/`
- Move: `projected_ed/` -> `ed/projected/`
- Move: `test_csr.jl` -> `ed/tests/test_csr.jl`
- Move: `dmrg/report/` -> `reports/fqahc_results/`

- [ ] 使用 `git mv` 完成机械迁移。
- [ ] 更新所有 Julia include、集群提交脚本、README 和 LaTeX 构建路径。
- [ ] 扩展 `.gitignore`，禁止波函数和 runtime output 进入 Git。
- [ ] 运行目录测试并确认通过。

### Task 4：导入 W003 验证过的 30-site 响应优化

**Files:**
- Modify: `ed/shared/basis.jl`
- Modify: `ed/shared/ksector.jl`
- Modify: `ed/shared/hamiltonian.jl`
- Modify: `ed/shared/optical_response.jl`
- Modify: `ed/cases/case3_30sites_Np12/run_optical_response.jl`
- Create/modify: `ed/tests/test_optical_response.jl`

- [ ] 先加入 translation cache 和 threaded response 的失败回归测试。
- [ ] 只导入 `/Users/shajianyu/CMP_manybody/PHC_Multinode` 中已在 W003 验证的差异。
- [ ] 在 W003 运行小体系 exact/response-Lanczos benchmark。

### Task 5：建立 CDW 结果包和静态结构因子输出

**Files:**
- Create: `results/phc30_np12_v1_100_cdw/parameters.toml`
- Create: `results/phc30_np12_v1_100_cdw/manifest.toml`
- Create: `results/phc30_np12_v1_100_cdw/data/*.dat`
- Create: `results/phc30_np12_v1_100_cdw/figures/*`
- Modify/Create: `ed/cases/case3_30sites_Np12/compute_cdw_structure_factor.jl`

- [ ] 写失败测试，要求能谱分类为三个基态，且第四态 gap 为 `0.0212691488`。
- [ ] 写失败测试，要求结构因子文件包含 `k=0,5,10`、三态平均和一致性诊断。
- [ ] 在 W003 从保存波函数重新计算 unit-cell normalized `N(q)`。
- [ ] 同步小型数据和图片回 Git；在清单中记录大波函数路径、大小和 SHA256。

### Task 6：修正能谱图和中文结果 note

**Files:**
- Modify: `reports/fqahc_results/plot_phc30_spectrum.py`
- Modify: `reports/fqahc_results/test_plot_phc30_spectrum.py`
- Modify: `reports/fqahc_results/dmrg_summary.tex`
- Generate: `reports/fqahc_results/figures/optical_response/phc30_manybody_spectrum.png`

- [ ] 先修改测试，要求 `manifold_size=3`、正确三态劈裂和三态到第四态 gap。
- [ ] 重画能谱并加入静态结构因子图。
- [ ] 将 `V1=100` 相关表述统一为强耦合三周期 CDW，并说明现有光学只算 `k=5,10`。
- [ ] 编译并渲染 PDF，检查新增页面和光学章节衔接。

### Task 7：W003、GitHub 和本地三方同步

- [ ] 在 W003 创建独立 checkout 并运行 scoped Julia/Python/LaTeX 验证。
- [ ] 提交所有小型代码、数据、图和清单。
- [ ] 推送 `integration/ed-excitations` 到 `origin`。
- [ ] 在 W003 fetch/checkout 远端分支，确认本地、GitHub 与 W003 commit 相同。
- [ ] 记录仍留在 W003 的大文件及其哈希，不上传 Git。
