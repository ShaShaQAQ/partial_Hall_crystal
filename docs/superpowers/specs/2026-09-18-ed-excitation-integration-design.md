# ED 激发计算整合设计

## 目标

将当前分散在 `feature/spectral-flow`、`feature/optical-response-lanczos`、
`DMRG` 和未受 Git 管理的 `PHC_Multinode` 中的 ED 能谱、谱流、静态结构因子与
光学响应工作整合到 `partial_Hall_crystal` 的单一开发线上，同时保持 DMRG/iDMRG
代码的独立模块边界。

## 仓库边界

仓库按计算方法和长期职责组织，而不是用长期分支区分方法：

```text
partial_Hall_crystal/
  ed/
    shared/                 # full two-band ED 公共基矢、晶格和 Hamiltonian
    projected/              # projected-band ED
    cases/                  # 具体有限尺寸几何与填充
    tests/                  # ED 回归测试和仓库结构测试
  dmrg/                     # finite/infinite DMRG
  reports/fqahc_results/    # 跨 ED/DMRG 的中文结果汇总
  results/                  # 可复现的小型派生数据、图片和运行清单
  docs/
```

`integration/ed-excitations` 从 `feature/spectral-flow` 建立，因为该分支已包含 CSR、
多节点 ED 和 flux insertion。将 `feature/optical-response-lanczos` 合入该分支，解决
公共 ED 文件冲突时同时保留多节点稀疏实现和经过小体系 benchmark 的物理电流算符。
DMRG 分支的已提交内容随后合入，使 `dmrg/` 和现有中文 note 进入同一历史。

## 目录迁移

- 根目录的 `shared/` 移至 `ed/shared/`。
- `case*` 移至 `ed/cases/`，保留原 case 名称以减少运行脚本变化。
- `projected_ed/` 移至 `ed/projected/`。
- 根目录 ED 测试移至 `ed/tests/`。
- `dmrg/report/` 整体移至 `reports/fqahc_results/`；它汇总 ED 和 DMRG，不能继续被
  命名为纯 DMRG 报告。
- 所有 Julia `include`、PBS/Slurm 路径和报告生成路径随迁移更新。

## 运行数据策略

Git 跟踪：

- 参数文件与运行清单；
- 能谱、结构因子和光学响应等小型文本数据；
- 最终图及其绘图脚本；
- 测试、日志摘要、数据哈希和生成代码的 Git commit。

Git 不跟踪：

- `*.jld2` 基态波函数；
- 稀疏矩阵、checkpoint、完整集群日志和临时输出。

每个正式参数点使用 `results/<run_id>/manifest.toml`。清单记录模型参数、基态
sector、代码 commit、W003 数据绝对路径、文件大小和 SHA256。当前强耦合 CDW
参数点命名为 `phc30_np12_v1_100_cdw`，不再标记为已确认 FQAHC。

## 当前 CDW 结果修正

当前 `V1=100,V2=V3=0` 能谱应分为：

- `k=0,5,10` 三个 CDW 基态，三态内部劈裂约 `1.245629e-4`；
- 第四态位于 `0.0213937117`，与三态流形相隔约 `0.0212691488`；
- 随后的 12 个态是低能激发，不属于 15 重基态简并。

静态结构因子使用已保存的三个基态波函数分别计算，并输出三态平均与 sector 间
最大差异。现有实现使用 `1/Nuc` 归一化，因此正文、代码注释和数据 header 必须统一
写为 unit-cell normalized，不能再误写成 `1/Ns`。

## W003 与 Git 同步

- GitHub 保存 `integration/ed-excitations` 分支的小型、可复现内容。
- W003 使用独立 checkout `/home/public/shajy/codex/partial_Hall_crystal-ed-excitations`，
  不修改已有生产目录。
- 已有大波函数保留在 `/home/public/shajy/codex_runs/phc30_optical_legacy/`；正式运行
  清单引用并校验这些文件，不复制到 Git。
- 所有 Julia 测试、结构因子重算和大体系验证在 W003 执行，本机只做 Git、文本与
  轻量绘图操作。

## 验收条件

1. 集成分支同时包含 multi-node/spectral-flow 与 response-Lanczos 历史。
2. 目录契约测试证明 ED、DMRG、报告和结果路径分离，大波函数模式被忽略。
3. ED 小体系测试和 response-Lanczos benchmark 在 W003 通过。
4. 当前 CDW 能谱图和 note 不再出现“15 重近简并基态”的错误判断。
5. 三个 CDW 基态的静态结构因子在 W003 从波函数重算，`q=5,10` 峰及一致性被记录。
6. GitHub 分支与 W003 checkout 指向同一 commit。
