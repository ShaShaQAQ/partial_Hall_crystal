# 30-site、Np=13 激发性质计算设计

## 目标

在与现有 `Np=12` 结果完全相同的 30-site 两带 ED 模型中，把粒子数改为
`Np=13`，分别计算以下两个参数点：

1. 强耦合点：`t1=1, t3=0.2, V1=100, V2=V3=0`；
2. 温和点：`t1=1, t3=0.2, V1=10, V2=V3=2`。

每个参数点产出全动量能谱、可复用的基态波函数、静态结构因子和
response-Lanczos 线性光电导，并将参数、算法、图和主要数值结果补入现有中文
结果汇总 `reports/fqahc_results/dmrg_summary.tex`。

## 物理口径

- 系统大小固定为 `Ns=30`、`Nuc=15`、`Np=13`，即 site filling 为
  `13/30`，每原胞粒子数为 `13/15`。
- 使用修正后的 `TiltedLat30` 晶格和 Fourier/平移约定。旧 `Np=13` 投影 ED
  数据只作为历史参照，不代替本次 full-ED 计算。
- 不预设强耦合点一定有三重基态，也不预设温和点一定有 15 态流形。每个动量
  sector 保存最低 8 个能级，汇总后依据排序能谱中的低能内部宽度和通往下一态的
  间隔识别实际低能流形。
- 结构因子和光学响应对识别出的低能流形逐态计算，并给出等权平均；同时保留逐
  sector/逐态数据，以免平均掩盖简并破缺或有限尺寸差异。
- note 中只根据本次能谱和结构因子描述晶体序与低能简并，不以这些量单独宣称
  拓扑序。若没有 flux insertion 或多体 Chern 数证据，温和参数点只能标为候选态。

## 代码边界

现有 `ed/cases/case3_30sites_Np12` 驱动中有多处把 `Np=12`、结果目录和流形大小
写死。本次把可复用计算入口参数化为任意合法 `Np`，但不重写共享 ED 内核：

- `run_spectrum.jl`：接受 `--Np 13`，每个 PBS array 元素负责一个动量 sector，
  保存 8 个能量和该 sector 最低态波函数、残差及完整模型元数据。
- 基态流形分析：从 15 个 partial JLD2 文件读取能谱和波函数，检查粒子数、模型
  参数、晶格指纹、残差和 sector 完整性；识别流形后输出能谱表、结构因子表、
  diagnostics、parameters 和带 SHA256 的远端波函数 manifest。
- `run_optical_response.jl`：从保存态中的 `Np` 重建 basis 和 sector，计算
  `J_x|g>` 出发的 response Lanczos；保存 total/regular/Drude 的实部与虚部、
  Lanczos 系数、残差和全部响应元数据。
- PBS 驱动：为两组 `Np=13` 参数使用独立 run directory，避免覆盖 `Np=12`
  数据。每个重型作业请求 W003 的 24 核/90 GB 单节点资源。
- 结果打包与作图：接受 result ID、`Np`、参数和实际流形 sector 列表，不再假设
  温和点必有 sectors `0:14` 或写死 `Np=12`。

## 数据流与存储

每个参数点使用独立远端目录：

```text
/home/public/shajy/codex_runs/<result_id>/
├── spectrum/partial_0.jld2 ... partial_14.jld2
└── optical/sector_<k>_optical_response.{jld2,dat}
```

大型 JLD2 波函数和完整响应检查点只保存在 W003。Git 跟踪：

```text
results/<result_id>/
├── parameters.toml
├── manifest.toml
├── optical_manifest.toml
└── data/
    ├── spectrum.dat
    ├── structure_factor_ground_manifold.dat
    ├── manifold_diagnostics.txt
    ├── optical_diagnostics.txt
    └── optical/*.dat
```

manifest 必须记录 W003 绝对路径、文件大小和 SHA256，使每幅图都能追溯到原始态。

## 调度顺序

1. 在 W003 上对两个参数点分别提交 15 个 sector 的能谱 array 作业。
2. 验证 30 个 partial 文件均存在、元数据一致、波函数归一且残差小于 `1e-7`。
3. 分别汇总能谱并识别实际低能流形。
4. 对流形中的每个态计算静态结构因子。
5. 对流形中的每个态提交 response-Lanczos 光学作业，默认沿用现有
   `M=600`、`eta=0.065`、`omega=0:0.002:10`。
6. 打包小型数据、生成能谱图、结构因子图和四幅 `sigma_xx` 图：regular/total
   各自的实部与虚部。
7. 将两组 `Np=13` 结果作为新小节补入现有中文 note，并重新编译 PDF。

若实测流形包含全部 15 个 sector，则计算全部 15 条光学曲线；若只有更小的隔离
流形，则只计算该流形并在 note 和 manifest 中明确 sector 列表。若能谱没有清晰
隔离流形，则先保留全能谱和全部 sector 基态，结构因子与光学平均使用最低能的
有限尺寸准简并组，并明确写出选择阈值和不确定性。

## 测试与验收

代码改动遵循先失败测试、再最小实现：

- 小粒子数/伪造元数据测试验证 `Np` 参数可传播到 spectrum、structure-factor 和
  optical 驱动，且错误粒子数会被拒绝。
- 打包与绘图测试验证任意流形 sector 列表、`Np=13` 表头和结果 ID，不允许悄悄
  回退到 `Np=12` 或固定 15 态。
- W003 正式结果要求：每个 sector 的基态残差 `<1e-7`，波函数范数接近 1，全部
  partial 文件参数/晶格指纹一致，光学输入态再次通过残差检查，数值表无 NaN/Inf。
- 最终运行 ED 测试、报告测试、绘图和 LaTeX 编译；PDF 无缺图、无 overfull box。
- 本地、W003 checkout 与 GitHub `integration/ed-excitations` 最终指向同一 commit；
  W003 现有未跟踪 PBS 日志不删除、不覆盖。

