# 从聚合物重复单元分配力场参数

`island.forcefields.assignment.assign_polymer_parameters` 接受线性双连接点 PSMILES、力场名称和聚合度。GAFF 族通过 RadonPy 返回带类型、参数和可选部分电荷的 RDKit 分子。PCFF 和 OPLS 族通过 EMC 返回可配套交给 LAMMPS 的拓扑与系数文件。

## 安装与使用

从仓库根目录安装所需的后端。以下命令安装本次短链验证使用的版本。EMC 的 setup 程序还需要 Perl。

```bash
python -m pip install -e '.[chemistry]' radonpy-pypi==0.2.11 emc-pypi==2026.11.19
```

```python
from pathlib import Path
from island.forcefields.assignment import assign_polymer_parameters

peg = assign_polymer_parameters("*CCO*", "gaff2", dp=3, charge="gasteiger", atom_order="radonpy")
for atom in peg.mol.GetAtoms():
    site = peg.site_ids[atom.GetIdx()]
    print(site, atom.GetProp("ff_type"), atom.GetDoubleProp("AtomicCharge"))

pcff = assign_polymer_parameters("*CCO*", "pcff", dp=3)
Path("peg.data").write_text(pcff.data)
Path("peg.params").write_text(pcff.parameters)
```

### GAFF 与 GAFF2

- `gaff`、`gaff2` 和 `gaff2_mod` 分别加载 RadonPy 的对应参数库。原子属性包含 `ff_type`、`ff_sigma`、`ff_epsilon`，键属性包含 `ff_type`、`ff_k`、`ff_r0`。`mol.angles`、`mol.dihedrals` 和 `mol.impropers` 保留原子索引和完整系数。结果还包含 1–4 缩放、RadonPy 版本和数据库 SHA-256。
- 对照 RadonPy 的氢封端链时使用 `atom_order="radonpy"`。RadonPy 遍历原子的顺序会影响少数共轭结构的交替类型，例如 PI547 的 `cc/cd` 和 `nc/nd`。该选项采用 RadonPy 构建器的原子和键顺序，同时把 ISLAND 链的立体化学、坐标及位点编号映射过去。DP 1 采用封端单体的规范 SMILES 顺序。默认沿用 ISLAND 的原子顺序。两种模式均可通过 `site_ids[atom.GetIdx()]` 找回 ISLAND 位点。
- 默认用氢封端。ISLAND 顺序支持 `head_end_group` 和 `tail_end_group`，例如 `"*C"`。重复单元可用两个 `*`、定向的 `[*:1]` / `[*:2]`，也可用 RadonPy 的 `[3H]` 连接点。若已有两个星号，额外的氚原子保留同位素含义。链内和端基的 E/Z 与四面体标记在连接后重建。氢封端产生的末端 `CH₂=` 不再有 E/Z。
- `assign_radonpy(mol, "gaff2", ...)` 可直接接收 RDKit 分子并在副本上分配，沿用输入图、原子顺序、构象和已有电荷。要复现给定 RadonPy 拓扑，可传入同一份显式氢分子。对于 ISLAND 的 `MolecularSystem`，使用其已有位点编号。需保留输入的双键构型时传 RDKit 分子或原始 PSMILES。
- `charge=None` 仅分配力场参数，保留 RDKit 输入中已有的 `AtomicCharge` 属性。可选 `gasteiger`、`zero`，以及安装量化依赖后的 `RESP`、`ESP`、`Mulliken`、`Lowdin`。PSMILES 入口为量化方法生成三维构象。`charge_options` 透传 RadonPy 的计算选项。四种量化方法都在氢封端的 `*C*`、DP 1 上用 Psi4、HF/STO-3G 和 `opt=False` 做过实际求值。`db_file`、`use_mdl`、`retry_mdl` 可指定参数文件及芳香性选项。

### PCFF 与 OPLS

- `pcff` 和 `pcff_ore` 载入对应的 EMC Class-II 库，并保留成键交叉项。`opls`、`oplsaa` 都指向 `opls/2024/opls-aa`。也支持精确指定年份和分辨率，例如 `opls/2012/opls-ua`。可选年份为 2012、2017 和 2024。全原子模型用氢封端，联合原子模型用甲基封端。
- `dp` 是每条链的重复单元数，`nchains` 是链数，`density` 是初始建盒密度（g/cm³）。默认单链、0.1 g/cm³。返回的 `field`、`field_sha256`、`emc_version` 和 `emc_smiles` 可用于追溯。EMC 安装目录优先取 `emc_root`，其次取 `EMC_ROOT`，最后取 Python 包自带版本。`field_dir` 可单独指定参数目录。
- 默认 `stereo_policy="strict"` 拒绝 EMC 无法保持的双键方向。复现 PolyJarvis 时显式选 `stereo_policy="polyjarvis"`，将方向键转成未指定构型，保留四面体标记。结果的 `bond_stereo_omitted` 标记这一处理。若研究需要固定顺反构型，请保留严格模式。
- 实际构建会校验所载入的参数文件、链数、原子与键数、系数有限性、逐链净电荷，并检测缺系数警告。当 EMC 报告大环识别深度不足时，用原参数库和自动环搜索深度重试，结果记录 `automatic_ring_depth`。对稠环意外多加氢的情况，用同一结构的 Kekulé 表示重试。全原子同位素输入按连通图定位并分出新的质量类型，复制原有 LJ 系数与类型对覆盖。PI240 和 PI854 已验证氘数和同坐标 LAMMPS 势能。联合原子同位素或无法唯一判定化学位置的标记会报错。

## 与参考实现的实测结果

### RadonPy 的 1,077 条聚合物记录

对 PI1070.csv 中的 1,077 个 PSMILES，分别用 ISLAND 和独立调用的 RadonPy 0.2.11 构建氢封端三聚体，并分配 Gasteiger 电荷。GAFF 和 GAFF2 均为 **1,077/1,077**。使用 `atom_order="radonpy"`，逐原子类型完全一致，LJ 参数和电荷的绝对差均小于 `1e-10`，逐键、逐角、逐项 proper / improper 的索引和系数一致。41 个含方向键的重复单元在 ISLAND 输出中都保有 E/Z，其中 25 个与 RadonPy 默认拼链器生成的立体 SMILES 有差异。逐条文件保留了双方的结构。这批结果使用 RDKit 2024.09.6。[固定版本的数据集](https://github.com/RadonPy/RadonPy/blob/c9d52bcafc05c1083acb593895e9573868538961/data/PI1070.csv)

### PolyJarvis 的公开成功记录

作者公开的两个试建批次中，有 586 个重复单元在至少一种目标力场下成功，共 611 个不同的“单元 × 力场”组合。ISLAND 在相同的 DP 4、每盒 2 条链、初始密度 0.9 g/cm³ 下复现了 **611/611**。ISLAND 的构建种子为 2026。这里使用作者版 2024 OPLS-AA，其他力场使用 EMC 发行库。[作者的逐条试建记录](https://github.com/arz-2/PolyJarvis/tree/c9fbfcccf617e59a89bda75bc4b166c6426ed4e8/docs/ff_coverage_sweep)

| 力场 | 成功复现 / 作者报告成功 |
| --- | ---: |
| PCFF | 575/575 |
| PCFF-ORE | 12/12 |
| OPLS-AA 2024，作者版 | 9/9 |
| OPLS-AA 2012 | 1/1 |
| OPLS-UA 2012 | 13/13 |
| OPLS-UA 2024 | 1/1 |

### 同一组 982 条 PSMILES 的逐库覆盖率

下表对每一个力场都使用相同的 982 个 PSMILES、DP 4、2 条链和种子 2026。严格模式要求保留指定双键构型，PolyJarvis 模式允许把它交给 EMC 作为未指定构型。严格列按同一批输出中成功且未省略双键方向的记录统计。[逐体系结果与失败诊断（CSV.gz）](evidence/reference_assignment.csv.gz)

| 固定参数库 | PolyJarvis 模式成功 | 其中省略了 E/Z | 严格模式成功 |
| --- | ---: | ---: | ---: |
| PCFF，EMC 发行库 | 678/982（69.0%） | 32 | 646/982 |
| OPLS-AA 2024，EMC 发行库 | 270/982（27.5%） | 8 | 262/982 |
| OPLS-AA 2024，作者修订库 | 205/982（20.9%） | 0 | 205/982 |

发行版 OPLS-AA 覆盖了修订版在这批数据上的全部 205 条，并多覆盖 65 条。另一项独立的 PDMS 三聚体测试需要作者补丁。PI508 的十二元环在自动加深环搜索后完成了 PCFF 分配和 LAMMPS 单点求值。其余输入的诊断包括缺分型规则、电荷增量或扭转参数。这里的比例适用于该数据集、指定端基和参数库。

逐条文件共 5,127 行。`field=opls` 对应 2024 发行库，`field=opls/2024/opls-aa` 对应本次安装的作者修订库。`published_positive` 表示该“单元 × 库”是否出现在作者的成功记录中，GAFF 列组另记录直接调用 RadonPy 的成功与逐项对照。`diagnostic` 保留失败时的关键诊断，复验脚本输出更完整的后端消息。

### 论文的八个代表体系

比较对象是 PolyJarvis v3 补充材料中的 8 个试构建案例，其中包括 PTFE。表内 EMC 结果使用 emc-pypi 2026.11.19 自带的 EMC 9.4.4（2026-10-07 构建）和随包参数。每个成功组合分别在三聚体和表列聚合度上做过单链分配。两种长度均做了 LAMMPS `run 0` 单点求值。[PolyJarvis 论文与补充材料](https://arxiv.org/html/2604.02537v3) · [作者的构建输入](https://github.com/arz-2/PolyJarvis/tree/c9fbfcccf617e59a89bda75bc4b166c6426ed4e8/data)

| 重复单元 | 论文 DP | ISLAND PCFF | ISLAND OPLS-AA | 论文 OPLS-AA |
| --- | ---: | --- | --- | --- |
| PE | 179 | 成功 | 成功 | 成功 |
| PEG | 114 | 成功 | 缺电荷增量 | 缺参数 |
| PLLA | 70 | 成功 | 缺电荷增量 | 缺参数 |
| aPS | 49 | 成功 | 成功 | 成功 |
| sPVC | 41 | 成功 | 缺电荷增量 | 缺参数 |
| PEEK | 18 | 成功 | 缺扭转系数 | 缺分型规则 |
| PSU | 12 | 成功 | 成功 | 缺分型规则 |
| PTFE | 50 | 成功 | 成功 | 成功 |

论文的 PCFF 试构建也是 8/8。表中 OPLS-AA 缺项来自三聚体真实分配时的报错。PSU 的新增覆盖依赖当前发行版的参数库，作者发布的本地修改版本具有不同的文件哈希。单点检查衡量参数文件的可执行性，密度、玻璃化转变温度及其他凝聚态性质需要另行验证。

在同一组 8 条显式氢封端链上，GAFF 与 GAFF2 都与独立调用 RadonPy 的原子类型、LJ、电荷、键、角以及每一项 proper / improper 系数逐项一致。测试覆盖 RadonPy 0.2.9 的三聚体与表列聚合度，以及 RadonPy 0.2.11 的三聚体与表列聚合度。[RadonPy 力场实现](https://github.com/RadonPy/RadonPy/tree/develop/radonpy/ff)

在 14 个三聚体组合上，还直接运行了作者同一修订的 `smiles_to_emc.py`（固定种子 2026）。比较了全部 `*_coeff` 行、质量、逐原子类型和电荷，以及带原子编号的键、角、二面角和 improper 表。8 个标准 PCFF 案例、1 个 PCFF-ORE、5 个 OPLS-AA/UA 案例均一致，其中包括作者扩展库的 PDMS。[PolyJarvis 试建脚本](https://github.com/arz-2/PolyJarvis/blob/c9fbfcccf617e59a89bda75bc4b166c6426ed4e8/mcp-servers/mcp-emc-server/smiles_to_emc.py)

### 更换 OPLS 参数来源

作者提供了扩展 OPLS-AA 的三个参数补丁。以下命令安装到独立目录，依次核对基线、补丁和最终文件的 SHA-256，需要系统命令 `patch`。离线安装可追加 `--checkout /path/to/PolyJarvis`。本次在作者版本上为氢封端 PDMS 三聚体分配了 32 个全原子位点并完成 LAMMPS 单点求值。随包 OPLS-AA 对同一输入报告缺少分型规则。[作者的参数清单与补丁](https://github.com/arz-2/PolyJarvis/tree/c9fbfcccf617e59a89bda75bc4b166c6426ed4e8/emc_fields)

```bash
python -m island.forcefields.polyjarvis_fields ./polyjarvis-field
```

```python
pdms = assign_polymer_parameters(
    "*[Si](C)(C)O*", "opls/2024/opls-aa", dp=3, field_dir="./polyjarvis-field"
)
```

### 复验

```bash
python -m pytest -q tests/test_reference_assignment.py
ISLAND_ASSIGNMENT_PAPER_DP=1 ISLAND_ASSIGNMENT_LAMMPS=/path/to/lmp \
    python -m pytest -q tests/test_reference_assignment.py
```

测试会按可用依赖运行。设置 `ISLAND_ASSIGNMENT_LAMMPS` 会对成功输出执行单点检查，所选 LAMMPS 需要支持 MOLECULE、CLASS2 和 KSPACE 等对应势函数。设置 `ISLAND_POLYJARVIS_CHECKOUT` 指向作者仓库可再测试参数安装器和直接脚本对照。下方大样本逐条检查参数分配，LAMMPS 单点另在上述代表体系上运行。

### 公开数据集复验命令

下列命令从仓库根目录运行。RadonPy 的数据文件名为 PI1070，此固定版本实际含 1,077 条记录。PolyJarvis 对其中 982 条走 EMC 的记录做了按类别和失败情况选择力场的试建。`--complete-emc-fields` 额外把所选力场逐一用于全部 982 条。可省略该选项来只回放作者报告成功的组合。CSV 记录输入、实际输出结构、每项比较和错误信息。[RadonPy 数据文件](https://github.com/RadonPy/RadonPy/blob/c9d52bcafc05c1083acb593895e9573868538961/data/PI1070.csv) · [PolyJarvis 公开试建记录](https://github.com/arz-2/PolyJarvis/tree/c9fbfcccf617e59a89bda75bc4b166c6426ed4e8/docs/ff_coverage_sweep)

```bash
source_url=https://raw.githubusercontent.com/arz-2/PolyJarvis/c9fbfcccf617e59a89bda75bc4b166c6426ed4e8
for arm in arm0 arm1; do
    curl -fL "$source_url/docs/ff_coverage_sweep/$arm.jsonl" -o "$arm.jsonl"
done
curl -fL https://raw.githubusercontent.com/RadonPy/RadonPy/c9d52bcafc05c1083acb593895e9573868538961/data/PI1070.csv -o PI1070.csv
python scripts/benchmark_reference_assignment.py --radonpy-panel PI1070.csv \
    --polyjarvis-archives arm0.jsonl arm1.jsonl --patched-fields ./polyjarvis-field \
    --complete-emc-fields pcff opls/2024/opls-aa opls --workers 4 --output results.csv
```
