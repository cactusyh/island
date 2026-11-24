# 从聚合物重复单元分配力场参数

`island.forcefields.assignment.assign_polymer_parameters` 接受线性双连接点 PSMILES、力场名称和聚合度，构造一条两端显式氢封端的链。GAFF、GAFF2 调用 RadonPy；OPLS-AA、PCFF 调用 EMC。调用方可以取出逐原子电荷和实际的成键、非键系数。对于已有 ISLAND 分子或已选择其他端基的 GAFF 链，可将其直接传给 `assign_radonpy(system, ...)`。

## 安装与使用

从仓库根目录安装所需的后端。以下命令安装本次短链验证使用的版本。EMC 的 setup 程序还需要 Perl。

```bash
python -m pip install -e '.[chemistry]' radonpy-pypi==0.2.11 emc-pypi==2026.11.19
```

```python
from pathlib import Path
from island.forcefields.assignment import assign_polymer_parameters

peg = assign_polymer_parameters("*CCO*", "gaff2", dp=3, charge="gasteiger")
for atom in peg.mol.GetAtoms():
    site = peg.site_ids[atom.GetIdx()]
    print(site, atom.GetProp("ff_type"), atom.GetDoubleProp("AtomicCharge"))

pcff = assign_polymer_parameters("*CCO*", "pcff", dp=3)
Path("peg.data").write_text(pcff.data)
Path("peg.params").write_text(pcff.parameters)
```

- **GAFF / GAFF2**　`mol` 是保留 RadonPy 属性的 RDKit 分子。原子上有 `ff_type`、`ff_sigma`、`ff_epsilon`，键上有 `ff_type`、`ff_k`、`ff_r0`。`mol.angles`、`mol.dihedrals` 和 `mol.impropers` 保存原子索引及完整的谐角、多项 Fourier 扭转、CVFF improper 系数。长度用 Å，能量用 kcal/mol，平衡角和相位用度。结果还包含实际的 1–4 缩放、RadonPy 版本和参数库 SHA-256。
- **电荷与坐标**　GAFF 族默认 `charge=None`，此时仅分配 LJ 和成键参数，结果中没有部分电荷。可显式选 `gasteiger`、`zero`，或使用安装了相应量化依赖的 RadonPy 计算 `RESP`、`ESP`、`Mulliken`、`Lowdin`。统一入口为前三种情形使用二维示意坐标，为量化方法生成三维构象；需要保留已有三维构象时使用 `assign_radonpy`。下表的 GAFF 验收使用 Gasteiger 电荷。
- **OPLS-AA / PCFF**　选择 `opls` 或 `oplsaa` 会载入 EMC 的 `opls/2024/opls-aa`，选择 `pcff` 会载入其 Class-II 数据库。返回的 `data` 和 `parameters` 是同一次 EMC 构建的 LAMMPS 文件内容，前者含 EMC 原子编号、电荷与相互作用实例，后者含系数、势函数、混合与特殊对设置。PCFF 保留各类交叉项。两份文件需配套使用。EMC 在 0.1 g/cm³ 下构建单链初始坐标。
- **来源与错误**　安装目录依次取显式 `emc_root` 参数、环境变量 `EMC_ROOT`、Python 包自带的 EMC。结果记录 EMC 构建版本和所选参数文件的 SHA-256。缺类型或系数时，异常包含 EMC 的诊断。EMC 路径支持四面体手性标记，对其无法保留的双键 E/Z 标记会报错。

## 与参考实现的实测结果

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

### 更换 OPLS 参数来源

运行作者的硅氧烷扩展时，可将 `emc_root` 指向含其已打补丁的 `field/opls/2024` 的 EMC 安装。本次按作者清单校验了三个补丁后文件的 SHA-256，在该版本上为 `*[Si](C)(C)O*` 的氢封端三聚体分配了 32 个原子的 OPLS-AA 参数，并通过 LAMMPS 单点求值。随包库对这个输入报告缺分型规则。此例使用作者选择的增量、电荷及端基扭转参数，具有独立的参数来源。[PolyJarvis 参数清单与补丁](https://github.com/arz-2/PolyJarvis/tree/c9fbfcccf617e59a89bda75bc4b166c6426ed4e8/emc_fields)

### 复验

```bash
python -m pytest -q tests/test_reference_assignment.py
ISLAND_ASSIGNMENT_PAPER_DP=1 ISLAND_ASSIGNMENT_LAMMPS=/path/to/lmp \
    python -m pytest -q tests/test_reference_assignment.py
```

测试会按可用依赖运行。设置 `ISLAND_ASSIGNMENT_LAMMPS` 会对成功输出执行单点检查，所选 LAMMPS 需要支持 MOLECULE、CLASS2 和 KSPACE 等对应势函数。完整覆盖数字来自安装两个分配后端并启用 LAMMPS 的运行。
