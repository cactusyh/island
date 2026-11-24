# J103 · 从聚合物配方到周期分子动力学

同一份 JSON 配方指定聚合物架构、重复单元及连接方式、链数、初始密度、力场和 MD 协议。运行器依次完成建图、PACKMOL 周期装填、完整图参数分配、NVT、NPT 和生产段。交联配方调用已有的周期固化流程，然后进入相同的三个 MD 阶段。输出包含实际力场来源、原子映射、原始 LAMMPS 输入、轨迹与续算文件。

```bash
python -m pip install -e '.[epoxy]' radonpy-pypi==0.2.11
python -m island.workflows.recipe examples/recipes/pe-gaff2.json work/pe \
  --lammps /path/to/lmp --packmol /path/to/packmol
```

PACKMOL 需支持周期盒，验证使用 21.1.1。LAMMPS 需要支持所选势函数，包括 GAFF 的 Fourier 二面角、PCFF 的 Class-II 交叉项、OPLS 的相应势函数和长程静电。本轮使用 LAMMPS 29 Aug 2024、EMC 9.4.4。

## 配方

```json
{
  "architecture": "block",
  "polymer": {
    "repeat_units": {"A": "[*:1]CC[*:2]", "B": "[*:1]CCO[*:2]"},
    "blocks": [["A", 3], ["B", 2], ["A", 1]]
  },
  "chains": 4,
  "density_g_cm3": 0.5,
  "seed": 2026,
  "forcefield": {"name": "pcff"},
  "md": {
    "temperature_k": 300,
    "pressure_atm": 1,
    "timestep_fs": 0.5,
    "nvt_steps": 1000,
    "npt_steps": 4000,
    "production_steps": 4000,
    "production_ensemble": "npt",
    "sample_every": 100,
    "segment_steps": 2000
  }
}
```

`density_g_cm3` 控制初始装填密度。NPT 中体积随动力学变化。示例三个阶段共 4.5 ps，结果将热力学平衡标为 `unmeasured`；研究体相性质时需要增加平衡和采样长度，并实际检查收敛。

| 架构 | 结构参数 | 可运行示例 |
|---|---|---|
| `homopolymer` | 一个 `repeat_units` 定义和 `dp` | [GAFF2 聚乙烯](../examples/recipes/pe-gaff2.json) |
| `block` | 至少两种重复单元，`blocks` 给出顺序和每段长度，也接受逐单元 `sequence` | [PCFF 醚三嵌段](../examples/recipes/block-pcff.json) |
| `branched` | 主链的 `dp`、`blocks` 或 `sequence`，加 `branches` | [OPLS-AA 双支链嵌段](../examples/recipes/branched-opls.json) |
| `crosslinked` | `cure` 指定 DGEBA 与 PACM 或 PETMP 的数量和反应协议 | [PCFF 环氧–胺](../examples/recipes/epoxy-pcff.json) |

每个主链配方只选择 `dp`、`blocks`、`sequence` 中的一项。有限链默认氢封端，也接受 `end_groups` 的头、尾单连接点 SMILES。支链使用例如 `{"anchor_site": 2, "psmiles": "*CC"}`，其中 `anchor_site` 是修改前主链图的一基位点编号，支链只有一个连接星号。`packing/chain.json` 和 `packing/site-map.json` 保存最终图与每条链的原子映射。多条链由同一指定拓扑的独立刚体放置得到。

有限链构建复用 `PolymerSpec` 的立体化学规则。受支持的单中心均聚物可选择 `tacticity`。该规范化器目前会拒绝含指定双键 E/Z 的配方和修改后无法保留声明构型的接点，具体原因写入诊断。已保有立体化学的完整 RDKit 图仍可使用直接参数分配接口。[现有图的用法](phase_4j94.md)

## 力场选择

GAFF 系列使用 `{"name":"gaff2","charge":"gasteiger"}` 等显式电荷设置，沿用 RadonPy 的原生类型、系数、LJ 混合与 1–4 缩放。已有完整图也可以直接调用 `assign_polymer_parameters(graph, "gaff2", charge="gasteiger")`。本轮数值验证采用 Gasteiger 电荷，其他已接入的量化电荷方法需要对应的量化计算依赖。

PCFF 和 OPLS 使用参数库的原生电荷。`opls` 指向 2024 OPLS-AA；也可以写完整的年份和 AA/UA 库名。扩展通过 `profile` 显式选择，例如 `pcff-carboxamide-v1` 或 `opls-epoxy-amine-v1`。程序将扩展安装到输出目录，并记录源文件指纹。自备 EMC 库通过 `field_dir` 传入，CLI 中相对路径以配方文件所在目录为基准。请求的力场或具体化学环境分配失败时，会保存该次原生诊断。

环氧–胺支持 PCFF 和带环氧扩展的 2024 OPLS-AA。环氧–硫醇支持 PCFF。`cure` 沿用 `CureConfig` 的反应设置，包括 `dgeba`、`pacm` 或 `petmp`、`target`、`cutoff`、`pre_steps`、`batch_steps`、`post_steps` 等。随机种子和装填密度统一使用顶层配置。结果同时记录实际转化率、目标是否达到和停止原因，固化步骤采用现有的 0.5 fs 协议，`md.timestep_fs` 用于固化后的三个阶段。[环氧反应与验证范围](phase_4j102.md)

## 续算和结果

```bash
python -m island.workflows.recipe examples/recipes/pe-gaff2.json work/pe --stop-after nvt
python -m island.workflows.recipe examples/recipes/pe-gaff2.json work/pe
```

两次调用需使用相同配方和相同可执行文件。程序校验已完成阶段中的文件 SHA-256，NPT 和生产段从上一段的 LAMMPS 二进制 restart 继承坐标、速度与盒尺寸。每隔 `segment_steps` 重新初始化长程静电网格，适应体积变化。生产段可选择 NVT 或 NPT。各段保存自己的局部步数和时长，`result.json` 给出三个 MD 阶段的总时长，固化记录另保存在 `cure` 中。

`--stop-after` 支持有限链的 `packing`、所有配方的 `assignment`、`nvt` 和 `npt`。分段已完成后可直接继续，失败的 MD 阶段从上一完整阶段重新运行。固化过程继续使用原有逐批检查点。失败目录保留为 `failed-阶段-编号`，`status.json` 给出出错阶段和诊断，同一输出目录由进程锁保护。

| 文件 | 内容 |
|---|---|
| `recipe.json`、`runtime.json`、`field/hashes.json` | 完整配方、模拟程序指纹和参数源指纹 |
| `assignment/` | 初始数据、完整系数和参数分配说明 |
| `nvt/`、`npt/`、`production/` | 各段 `md.in`、系数、轨迹、末态、热力学日志和 restart |
| `result.json` | 完成后的实际体系、各阶段结果及固化目标达成情况 |
| `status.json` | 运行、暂停、完成或失败状态 |

在任一完整 MD 阶段目录执行 `lmp -in md.in` 可独立复跑该段。恢复整份配方时使用原始运行目录及其校验文件。

## 实际验证

通过公开 CLI 独立运行了 15 个组合。均聚物为四条 PE 八聚体，嵌段为四条乙烯/氧乙烯 A₃B₂A₁，支链体系为四条乙烯/苯乙烯 A₃B₂A₁，每条链在 2、4 号位点分别接乙基和丙基。三种有限链架构都分别使用 GAFF、GAFF2、PCFF 和 2024 OPLS-AA。醚嵌段的 OPLS-AA 显式选用 J100 分型扩展。12/12 组合完成装填、完整参数分配及 1,000 步 NVT、4,000 步 NPT、4,000 步 NPT 生产段。每条轨迹实际解析全部粒子 ID、有限坐标和速度及最后一个采样步。

| 交联配方 | 力场 | 固化后转化率 | 固化端点绕盒秩 | 后续 MD |
|---|---|---:|---:|---|
| 4 DGEBA + 2 PACM，种子 41 | PCFF | 75% | 1 | 三段均完成 |
| 4 DGEBA + 2 PACM，种子 41 | OPLS-AA 环氧扩展 | 75% | 1 | 三段均完成 |
| 4 DGEBA + 2 PETMP，种子 173 | PCFF | 100% | 3 | 三段均完成 |

力场分配、已有整图与环氧流程及 J103 回归测试共 93 项通过，随后增加一项 OPLS 缺项诊断测试并复跑全部 13 项 J103 测试。GAFF2 在同一苯乙醚构象上的 LAMMPS 能量分项与逐原子力，对照 RadonPy 原生数据文件，差值均在 10⁻¹⁰ 以内。短配方在 NVT 后暂停再恢复，与直接运行的生产段末态逐字节相同；检查点内容变化和配方变更会被识别。15 个矩阵案例也全部再次通过完成阶段的校验复用。

这些结果验证短时周期模拟的执行流程。当前例子的链长、链数和 每个案例 4.5 ps 的 MD 采样，尚不足以建立材料密度、Tg 或力学性质的收敛结论。[逐案例阶段结果](evidence/phase_4j103/matrix.json) · [验证范围和归档指纹](evidence/phase_4j103/summary.json) · [原始输入、热力学日志和末态](evidence/phase_4j103/native-runs.zip)

归档包含所有 15 个组合、45 个 MD 阶段的独立复跑文件及末态。完整轨迹保留在均聚物 GAFF2、嵌段 PCFF、支链 OPLS-AA 和环氧–胺 PCFF 四个代表组合中。每个阶段可在解压后的对应目录执行 `lmp -in md.in`。要重新生成完整矩阵及工作流检查点，可运行下列命令。

```bash
python scripts/phase_4j103_recipes.py work/j103-matrix --workers 3   --lammps /path/to/lmp --packmol /path/to/packmol
```
