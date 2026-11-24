# J94 · 在已建聚合物图上分配 PCFF / OPLS 并运行动力学

已有嵌段或支链时，把完整图交给 `assign_polymer_parameters(graph, field)`。嵌段之间的键、侧链接点和端基都参与 EMC 分型及成键参数查找。可以直接使用构建器返回的 `MolecularSystem`、`normalize(PolymerSpec(...))` 返回的 `ChemicalGraph`、`FinalChemicalGraph` 或 RDKit `Mol`。图已确定链长和端基，因此调用时省略 `dp`。

## 从现有结构开始

在仓库根目录安装 Python 依赖。EMC 的 setup 脚本使用 Perl。动力学需要 LAMMPS 可执行文件，编译时需包含所选势函数，例如本例用到的 MOLECULE、CLASS2 和 KSPACE。

```bash
python -m pip install -e '.[chemistry]' emc-pypi==2026.11.19
```

```python
from island.builders import build_block_copolymer
from island.forcefields.assignment import assign_polymer_parameters, run_emc_md
from island.forcefields.assignment_core import PolymerSpec, normalize
from island.forcefields.assignment_core.schema import Branch

# 已经有这些图时，从下面的分配调用开始即可。
ether_block = build_block_copolymer([("*CC*", 3), ("*CCO*", 2), ("*CC*", 1)])
branch = normalize(
    PolymerSpec((("A", "*CC*"),), 4, branches=(Branch(4, "*CC"),))
)

for graph, field, directory in [
    (ether_block, "pcff", "ether-pcff"),
    (ether_block, "opls/2012/opls-aa", "ether-opls-aa"),
    (branch, "opls/2024/opls-ua", "branch-opls-ua"),
]:
    assignment = assign_polymer_parameters(graph, field, seed=2026)
    run = run_emc_md(
        assignment, directory, lammps="lmp", temperature_k=300,
        timestep_fs=0.5, steps=20_000, sample_every=1000,
    )
    print(assignment.field, run.time_fs, run.directory / "trajectory.lammpstrj")
```

`run_emc_md` 在新目录中写出完整 LAMMPS 输入，先进行共轭梯度最小化，再初始化速度并执行固定体积的 Nosé–Hoover NVT 积分。初始帧取自最小化后的结构。温度是热浴设定值，短链的瞬时温度会波动。

| 文件 | 用途 |
| --- | --- |
| `chain.data`、`chain.params`、`md.in` | 拓扑、初始坐标、质量、电荷、势函数与全部系数，以及积分设置。在这三个文件的新副本所在目录执行 `lmp -in md.in` 即可复跑。 |
| `assignment.json` | 输入结构、SMILES 遍历对应的位点、LAMMPS 粒子与位点的映射、补加的氢、坐标来源、EMC 版本和参数文件哈希。 |
| `trajectory.lammpstrj` | 含粒子 ID、分子 ID、类型、展开坐标（Å）和速度（Å/fs）。可在 OVITO 载入。 |
| `final.restart`、`log.lammps`、`result.json` | 末态续算文件、热力学日志、通过步数与逐帧检查后的结果摘要。 |

## 结构、坐标和参数覆盖

- 全原子计算保留每个输入位点的编号；隐式氢得到新的位点编号。联合原子模型用一个粒子表示重原子及被折叠的氢。`assignment.site_to_atom_id[site_id]` 可以追踪指定接点在轨迹里的粒子编号，`assignment.site_groups[atom_id]` 给出它代表的全部位点。分配时核对分子归属、元素、映射后的键端点、位点覆盖和各分子的总电荷。
- 已有三维坐标时沿用输入构象。只有图或二维图时，以 RDKit ETKDG 生成三维构象。已有全周期正交盒时保留其尺寸，按照最小镜像键向量展开有限分子。指定盒中的多分子体系需要提供已装填坐标。没有盒时，各方向在坐标范围两侧各留 12 Å；这些示例是单链稀盒。跨边界绕盒成环的共价网络需要显式周期成键描述，超出这个有限图入口的范围。
- 带 E/Z 的输入可用保留方向信息的 RDKit 分子。入口会核对指定构型与用于模拟的三维坐标。已经在上游图格式转换中丢失的 E/Z 需从保有构型的源结构传入。有同位素标记时使用全原子模型。
- `opls` 默认为 EMC 发行版 `opls/2024/opls-aa`。在这里的醚嵌段上，2024 发行库和本次核对的 PolyJarvis 2024 修订库都缺少 `{c4, o2e}` 电荷增量。接口保留 EMC 诊断并终止该次分配。表中的全原子醚结果显式使用完整覆盖该图的 2012 库。其他自备 EMC 库可通过 `field_dir` 指定。[参数库选择与 PolyJarvis 修订库安装](reference_assignment_zh.md)

## 实际运行的体系

用 EMC 9.4.4 和 LAMMPS 29 Aug 2024，对下列每个列出的组合均执行 20,000 步，步长 0.5 fs，目标温度 300 K。每条 10 ps 轨迹有 21 帧，包含第 0 步和最后一步。

| 已建分子图 | 完成积分的参数库 |
| --- | --- |
| 乙烯 / 氧乙烯 / 乙烯三嵌段 A₃B₂A₁，40 个原子 | PCFF、2012 OPLS-AA、2024 OPLS-UA |
| 乙烯 / 苯乙烯 / 乙烯三嵌段 A₃B₂A₁，58 个原子 | PCFF、2024 OPLS-AA |
| 乙烯四聚体的 4 号位点接乙基，32 个原子 | PCFF、2024 OPLS-AA、2024 OPLS-UA |
| 在乙烯 / 苯乙烯三嵌段的 2、4 号位点分别接乙基和丙基，73 个原子 | PCFF、2024 OPLS-AA |
| 乙烯 / 苯乙烯 A₃₀B₁₀，DP 40，342 个原子 | PCFF、2024 OPLS-AA |

这些单链运行验证了所列图的接点分型、参数导出及有限温度积分。研究熔体密度、玻璃化转变或输运性质，还需要对应的多链初始装填、平衡协议和采样长度。[逐组结果](evidence/phase_4j94/summary.json) · [12 组输入、轨迹、日志和续算文件](evidence/phase_4j94/graph_md_runs.zip)

将双支链三嵌段的 PCFF 输入从压缩包单独解出后，用同一 LAMMPS 版本复跑，末步势能同为 39.12334665 kcal/mol。J94 与现有参考分配测试共 72 项通过，其中包括聚乙烯整图与原重复单元入口的原生系数对照，以及对嵌段接点、支化点、跨盒分子、同位素和双键构型的检查。[整图回归测试](../tests/test_phase_4j94.py)

从仓库根目录可用以下命令重新生成整批数据，输出路径应为新目录。

```bash
python scripts/phase_4j94_graph_md.py ./j94-runs --lammps /path/to/lmp
```
