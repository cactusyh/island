# J105 · 多架构周期模拟覆盖率与装填构象重试

J105 用同一套公开配方入口测量从 PSMILES 到周期 MD 的各阶段覆盖情况。固定面板包含 PE、PS、PEG、PCL、PA6、PMMA，每类各有均聚物、嵌段和乙基支链体系，并分别使用 GAFF2、PCFF 和 2024 OPLS-AA，共 54 个组合。PCFF 显式选择 `pcff-carboxamide-v1`，OPLS 显式选择 `opls-polymer-v1`，GAFF2 使用 Gasteiger 电荷。

## 装填改进

首次装填沿用 J103 的 ETKDG 构象。PACKMOL 明确报告空间约束不可满足或未达到装填容差时，再用随机坐标初始化的 ETKDG 生成构象，最多尝试 4 次（含首次）。后续构象采用确定性种子，仍通过原有图适配器的坐标和立体化学校验。PACKMOL 种子、原子顺序、化学拓扑、链数、装填容差和请求的密度保持不变。

`packing/attempt-000/` 等目录保留每次的构象、PACKMOL 输入和原始日志；`packing/packing-attempts.json` 记录构象生成方式、种子和结果。最后成功的一次同时发布到原有 `packing/` 文件路径，保持后续分配及续算接口兼容。输入或运行程序错误不触发构象重试；达到上限仍失败时，整个阶段进入已有失败检查点，保留诊断。

重试不保证长链或高密度体系一定可装填。当前仍以同一构象的刚体副本进行装填，尚未引入多链独立构象、软势压实或长链逐段生长。

## 可复跑的覆盖率面板

```bash
python scripts/benchmark_polymer_architectures.py work/j105 \
  --workers 3 --lammps /path/to/lmp --packmol /path/to/packmol
```

每个组合使用 4 条六重复单元链、氢封端、初始密度 0.6 g/cm³、种子 173。嵌段为 A₃B₃，B 通常为 PE；PE 的 B 为 PEG。支链体系在主链可替换的 C–H 位点接一条乙基支链。PA6 的氢封端碳酰基端是醛端，这里按实际生成的图验证，不能代表所有实验尼龙端基。

MD 协议为 300 K、1 atm、0.5 fs 步长，依次进行 1,000 步 NVT、8,000 步 NPT 和 4,000 步 NPT 生产段，共 6.5 ps。成功标准包括每段末步、轨迹粒子 ID、有限坐标与速度、热力学记录及非空 restart 文件；`single_point` 记录最小化后第一次 `run 0` 的有限势能。它与原始构型的独立能量精度验证不是同一项检查。

脚本逐案例写出 `results.json`，并按力场、化学类型和架构更新 `summary.json`。失败输入保留在分母中；未完成面板标记 `complete=false`。复跑相同输出目录时复用已通过文件校验的阶段，失败案例从最后一个完整阶段继续。

可直接运行的配方：[PCL / OPLS](../examples/recipes/pcl-opls.json)、[PA6 / PCFF](../examples/recipes/pa6-pcff.json)、[PEG 支链 / OPLS](../examples/recipes/peg-branched-opls.json)。例如：

```bash
python -m island.workflows.recipe examples/recipes/pcl-opls.json work/pcl \
  --lammps /path/to/lmp --packmol /path/to/packmol
```

## 实际结果

原来的单构象流程完成 42/54 个组合，其余 12 个都停在 PACKMOL 装填阶段。加入有界构象重试后，完成 **54/54** 个组合，其中恢复 12 个，原成功案例退化 0 个。恢复案例保留相同密度、链数、力场和 MD 协议。

| 化学类型 | GAFF2 完成 MD | PCFF 完成 MD | OPLS 完成 MD |
|---|---:|---:|---:|
| PE | 3/3 | 3/3 | 3/3 |
| PS | 3/3 | 3/3 | 3/3 |
| PEG | 3/3 | 3/3 | 3/3 |
| PCL | 3/3 | 3/3 | 3/3 |
| PA6 | 3/3 | 3/3 | 3/3 |
| PMMA | 3/3 | 3/3 | 3/3 |

每格分母 3 为均聚、嵌段、支链三个配方。42 项相关回归测试通过。

[原始流程结果](evidence/phase_4j105/baseline.json) · [逐案例最终结果](evidence/phase_4j105/results.json) · [分阶段统计](evidence/phase_4j105/summary.json) · [验证和文件指纹](evidence/phase_4j105/verification.json) · [原始模拟文件](evidence/phase_4j105/native-runs.zip)

归档保留各 MD 阶段的 `md.in`、数据、系数、起始 restart（适用时）、末态和原始日志，可在各阶段目录独立执行 `lmp -in md.in`。另外保留 PE/GAFF2 均聚、PEG/OPLS 支链、PCL/OPLS 均聚和 PA6/PCFF 均聚的三段完整轨迹。归档内 `manifest.json` 列出每个文件的 SHA-256。它是阶段独立复跑归档，不包含恢复整份配方所需的全部检查点文件。使用面板命令的 `--archive work/native-runs.zip` 可重新生成同类归档。

## 验证边界与后续开发

这组面板衡量小体系、短链和短时模拟的执行能力。交联体系仍使用 J100/J102/J103 已验证的环氧–胺、环氧–硫醇流程，本面板没有扩大交联反应类型。54 个配方也不构成对任意嵌段组成、支化程度、端基或分子量分布的保证。

下一步应优先处理三个问题：

1. **更多化学环境的力场分配。** 根据 J104 的 982 条诊断逐类补齐 OPLS 分型与电荷增量，并继续独立核对原始参数来源；PCFF 也要保留其剩余缺项的真实失败结果。
2. **长链和密集体系的建模。** 引入多链独立构象、可配置的压实/退火协议，并验证拓扑、立体化学与周期边界；将刚体 PACKMOL 装填的限制与力场覆盖率分开记录。
3. **平衡及物性可信度。** 扩展链长、链数、种子与模拟时长，检查密度、链尺寸和应力的收敛，再建立 Tg、弹性模量等物性的实验或独立模拟对照。目前 `equilibration` 仍为 `unmeasured`。
