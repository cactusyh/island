# J107：完整图覆盖率、预检、缺项定位与来源绑定缓存

本阶段从 J106 提交 `eec859a` 创建 `codex/j107-final-graph-coverage`。J106 已推送至 GitHub。J93 的独立 WIP 分支、J104/J105 的阶段边界及历史审查收据保持独立。

J107 按审查规划补齐完整图的覆盖率测量和诊断接口。参数分配、数值正确性、实际 MD 与平衡/物性仍分层报告。本阶段不添加零参数、通用替代项或未经独立来源支持的数值参数，也不把预检或缓存命中当作新的科学验证。

## 使用方式

### 先检查最终图

```python
from island.forcefields.assignment_core import PolymerSpec, normalize
from island.forcefields.assignment_diagnostics import preflight_assignment

graph = normalize(PolymerSpec((("A", "*CCO*"),), 4))
report = preflight_assignment(graph, "opls", field_dir="work/opls-polymer")
```

预检不生成三维构象、不运行 EMC。它检查封端闭壳层图、参数库元素声明、UA 同位素限制，并报告组分形式电荷、已指定立体中心、位点编号和半径二局部化学环境。离子/反离子即使净电荷抵消也会提示其需要来源支持的离子参数。`ready_for_assignment=true` 仅表示基本结构/域检查通过；`parameter_completeness` 和 `md` 保持 `unmeasured`。

```bash
python -m island.forcefields.assignment_diagnostics '*CCO*' \
  --field opls --dp 4 --field-dir work/opls-polymer --output preflight.json
```

### 完整分配、能力查询与缓存

```python
from island.forcefields.assignment import assign_polymer_parameters

assignment = assign_polymer_parameters(
    "*CCO*", "opls", dp=4, field_dir="work/opls-polymer",
    cache_dir="work/assignment-cache",
)
```

缓存默认关闭。启用后，仅复用成功的完整图 EMC 分配；键绑定图的实际原子顺序与立体信息、位点编号、全部坐标、盒边界、随机种子、整个指定参数目录、EMC 可执行文件/脚本/模块/模板、Perl、相关适配代码及依赖版本。复用前仍验证输入图和几何。缓存不能放在其绑定的参数或运行源码目录内部；文件校验失败会显式报错，分配期间来源改变也不会写入缓存。JSON 原子写入避免半文件读取，不使用 pickle。缓存不绕过三维构象验证，也不适用于 GAFF 或旧 EMC 重复单元模式。

```bash
python -m island.forcefields.capabilities probe opls --repeat '*CCO*' --dp 4 \
  --field-dir work/opls-polymer --cache-dir work/assignment-cache
```

能力查询接口已修正 J106 后的兼容问题：默认检查并分配完整图，不再无条件传入旧模式的 `stereo_policy`。如需旧行为，显式指定 `--emc-mode repeat`。支持自定义 `--head-end-group` / `--tail-end-group`。查询结果保留 `single_point`、`dynamics`、`properties` 的 `unmeasured` 状态。

### 查看缺项

`assignment_failure(error, graph)` 保留完整异常文本、多类别标签、缺少的类型组合，以及 EMC 明确报告的组/位点。能够映射至最终图时返回真实稳定编号和对应局部环境；原生错误只报告类型组合时会明确标记 `backend_did_not_identify_sites`，不捏造原子定位。一次失败可能只暴露首批缺项，不意味着已穷举该输入的全部缺参数。

## 大面板发现并修复的端部芳香氮问题

部分聚酰亚胺重复单元经 RDKit 芳香性识别后，连接点邻居是带取代基的芳香 N，其 `NoImplicit` 并不为真。移除末端 dummy 后，需要显式的 H 计数才能成为 `[nH]`；仅靠隐式价态补齐会导致无法 Kekulize。J107 对这类默认端基补 H，DP 1/2/4 的吡咯与聚酰亚胺构图均与显式 `*[H]` 封端的独立结果一致。其他默认封端和自定义端基保持原语义。

重测脚本也保证两条入口独立执行：旧 EMC 路径不会因 ISLAND 最终图预检失败而被提前跳过。早期被这项构图检查阻断的记录单独保留，在修复后复测；最终比较使用独立执行后的结果。

## 完整图环扫描策略

原完整图入口无条件使用 EMC 自动环扫描，在较大多环分子上造成显著额外计算。J107 与原生重复单元入口采用一致策略：先用发行版标准深度，仅在后端明确提示 `increase scanning depth?` 时自动重试。结果的 `automatic_ring_depth` 字段记录是否触发重试；未改变力场数值参数、最终位点映射或返回的输入坐标。小模型的两种扫描方式逐项系数/拓扑一致；单独聚酰亚胺探针从超过原时限改善为约 15 秒完成。

大型环图的 `Inconsistent connectivity in rings` 不是已解决的问题。独立探针将 PI736 / PCFF 强制恢复为自动环扫描，仍于约 18.5 秒报相同错误；至少这个案例不能靠恢复自动深度解决。探针与脚本单独保存，不计入主面板分母，也不把结论外推至所有环图失败。

最终面板对所有完整图任务重新执行。旧模式中已完成、来源哈希相同且确实进入原生分配的结果直接保留；被早期构图预检阻断的旧模式记录重新执行。初次探索运行作为未完成的历史记录单独保存，不与最终分母混淆。

## 固定 982 条成对重测

使用原有 `reference_assignment.csv.gz` 的 982 个唯一输入；PCFF 采用 `pcff-carboxamide-v1`，OPLS 采用 `opls-polymer-v1`。每个输入在每种力场分别走完整图和旧重复单元入口，共 3,928 个任务。均 DP 4、单链、种子 2026；旧 EMC 建盒密度设为 0.03 g/cm³，与完整图后端内部临时构建设置一致。完整图返回的是经过验证的输入构象，旧模式返回 EMC 生成的构象，两条路径的坐标不作为覆盖率比较条件。

旧模式显式使用 `stereo_policy="polyjarvis"`，另统计其未省略双键方向的成功数；完整图保留指定的立体化学。J104 的表格采用两链、0.9 g/cm³，本轮不将不同构建设置的数字混为同一实验。

额外加入 34 个输入，各使用两种力场：Si/P/S/Cl/Br/F、明确手性与 E/Z、带电重复单元，氢/甲基/TMS 端基，以及三个显式离子组合。补充面板和主面板分开统计，共计 3,996 个任务。每个任务运行在独立进程组，默认 120 s 时限；超时保留在分母中，归为运行/几何限制而不是参数库缺失。新运行器会保留超时前已完成的预检。

本次最终记录 **3,996 / 3,996**，任务唯一性和分母检查通过。

| 主面板入口 | 分配成功 / 总数 | 成功率 | 未省略指定双键方向的成功数 |
|---|---:|---:|---:|
| PCFF 旧重复单元 | 899 / 982 | 91.5% | 865 |
| PCFF 完整图 | 880 / 982 | 89.6% | 880 |
| OPLS 旧重复单元 | 406 / 982 | 41.3% | 398 |
| OPLS 完整图 | 440 / 982 | 44.8% | 440 |

成对比较使用相同输入、来源、DP 和种子：

| 力场 | 两者成功 | 仅完整图成功 | 仅旧入口成功 | 两者失败 |
|---|---:|---:|---:|---:|
| PCFF | 875 | 5 | 24 | 78 |
| OPLS | 405 | 35 | 1 | 541 |

**旧入口成功、完整图失败的案例仍是明确缺口，不能用净增成功数抵消：**

| 力场 / 原因 | 个数 | 输入 ID |
|---|---:|---|
| OPLS / atom_typing | 1 | PI837 |
| PCFF / EMC 环连接检查拒绝 | 20 | PI1012, PI1013, PI1014, PI1021, PI1023, PI1025, PI1026, PI1028, PI1038, PI1042, PI1043, PI1044, PI1048, PI1078, PI1100, PI1101, PI1114, PI736, PI913, PI969 |
| PCFF / timeout | 3 | PI1017, PI1039, PI1047 |
| PCFF / torsion | 1 | PI936 |

环图失败的 RDKit 构图及预检均已通过，不能视作输入结构无效；PI936 缺少 `[cp:c1:o_2:c_1]` 扭转项；PI837 为含萘环重复单元的完整图分型失败，旧入口分型成功，具体原因仍需独立定位。逐例原始诊断见 [入口差异证据](evidence/phase_4j107/repeat-only-cases.json)。

完整图在 PI1017、PI1039、PI1047 上的两种力场均达到 120 秒时限，共 6 次，全部计入失败分母；没有 worker_error。只有完成参数分配的任务计为成功，超时不能据此断言参数缺失。

完整图主面板失败类别（类别可重叠）：

| 类别 | PCFF | OPLS |
|---|---:|---:|
| `angle` | 0 | 9 |
| `atom_typing` | 11 | 396 |
| `charge_balance` | 0 | 4 |
| `charge_increment` | 39 | 65 |
| `graph` | 21 | 20 |
| `timeout` | 3 | 3 |
| `torsion` | 28 | 74 |

主面板含相应元素的输入子集（**不是该元素所有类型的覆盖率**）：

| 元素 | PCFF 成功 / 含该元素输入 | OPLS 成功 / 含该元素输入 |
|---|---:|---:|
| N | 312 / 383 | 17 / 383 |
| O | 658 / 760 | 292 / 760 |
| Si | 未包含 | 未包含 |
| P | 1 / 16 | 0 / 16 |
| S | 65 / 77 | 10 / 77 |
| F | 112 / 134 | 4 / 134 |
| Cl | 29 / 41 | 26 / 41 |
| Br | 15 / 16 | 0 / 16 |
| I | 2 / 2 | 0 / 2 |

无法完成元素分类的输入数：PCFF 0，OPLS 0，仍保留在主分母中。

补充的 34 个完整图输入：

| 力场 | 分配成功 / 总数 |
|---|---:|
| PCFF | 29 / 34 |
| OPLS | 21 / 34 |

| 补充输入（H 为默认氢端基） | PCFF | OPLS |
|---|---|---|
| PA6-H | 成功 | 成功 |
| PA6-methyl | 成功 | 成功 |
| PDMS-H | 失败：charge_increment | 失败：atom_typing |
| PDMS-TMS | 成功 | 失败：atom_typing |
| PDMS-methyl | 失败：charge_increment | 失败：atom_typing |
| PEG-H | 成功 | 成功 |
| PEG-methyl | 成功 | 成功 |
| PET-H | 成功 | 成功 |
| PET-methyl | 成功 | 成功 |
| PP_chiral-H | 成功 | 成功 |
| PP_chiral-methyl | 成功 | 成功 |
| PTFE-H | 成功 | 成功 |
| PTFE-methyl | 成功 | 成功 |
| PU-H | 成功 | 成功 |
| PU-methyl | 成功 | 成功 |
| PVBr-H | 成功 | 失败：charge_increment |
| PVBr-methyl | 成功 | 失败：charge_increment |
| PVC-H | 成功 | 失败：charge_increment |
| PVC-methyl | 成功 | 失败：charge_increment |
| ammonium-H | 成功 | 成功 |
| ammonium-chloride | 成功 | 成功 |
| ammonium-methyl | 成功 | 成功 |
| cis_diene-H | 成功 | 成功 |
| cis_diene-methyl | 成功 | 成功 |
| lithium-ether | 失败：charge_balance | 成功 |
| phosphoester-H | 失败：charge_increment | 失败：atom_typing |
| phosphoester-methyl | 失败：charge_increment | 失败：atom_typing |
| polysulfone-H | 成功 | 成功 |
| polysulfone-methyl | 成功 | 成功 |
| polythioether-H | 成功 | 失败：torsion |
| polythioether-methyl | 成功 | 失败：torsion |
| sodium-acetate | 成功 | 成功 |
| sulfonated_PS-H | 成功 | 失败：atom_typing |
| sulfonated_PS-methyl | 成功 | 失败：atom_typing |

补充面板显示端基会改变可用性：PCFF 的 PDMS 只有 TMS 端基模型成功，默认 H 和甲基端基分别缺少 `{osi, si}`、`{o, sio}` 电荷增量。两种力场的磷酸酯模型均失败；OPLS 的 PVC/PVBr 缺 `{c4, cl1}` / `{c4, br1}` 电荷增量，聚硫醚缺扭转项。显式离子组合的分配成功仅说明这几个固定模型通过分配与电荷检查，没有进行离子体系 MD 或热力学验证。

原参数库未新增数值项。本轮提升主要来自完整端部/接点处理及运行路径修复；当前结果不能证明已经覆盖任意 PSMILES、任意交联图或全面达到 RadonPy/PolyJarvis 的能力。GAFF/GAFF2 未参加本轮 982 条成对覆盖试验。

成对结果、按元素分组、缺项类型组合和局部化学环境频率见证据目录。元素分组可重叠；无法完成结构预检的输入单列，不从总分母消失。计数只表示在固定 DP/端基/参数来源上的成功，不能外推为任意链长或所有交联环境的覆盖率。

复现：

```bash
python -m island.forcefields.emc_profiles pcff-carboxamide-v1 work/pcff
python -m island.forcefields.emc_profiles opls-polymer-v1 work/opls-polymer
python scripts/benchmark_final_graphs.py \
  --pcff-field work/pcff --opls-field work/opls-polymer \
  --out work/j107-panel --workers 4 --timeout 120
```

原始输入、来源和时限不同会拒绝复用原运行目录。失败和成功均保存逐例 JSON；`summary.json` 只有在全部唯一任务齐全时才标记 `complete=true`。缺项类别可以重叠，不应把类别之和当作失败体系数。

## 独立化学、能量/力与 MD 检查

- 复验七个 OPLS 源残基的所有类型/电荷及四个外推接点模型；PCFF 复验伯/仲/叔酰胺和双酰胺的来源电荷及重排不变性。
- 扩展检查 DMES、EMES、MSME、ETBR、ETCL 五个含硫/卤素模型，所有原子类型和电荷与 `.define` 独立残基定义一致；元素通过源 `.prm` 声明读取，避免把 Br/Cl 错当作单字符元素。
- DMSO 残基虽存在，当前库仍缺少 `[c4s:s3:c4s:h1]`、`[h1:c4s:s3:o1s]` 扭转项，明确保留为失败对照。
- 醚接点、酰胺接点、支化酯在 PCFF/OPLS 上共六组扰动构型进行独立有限差分和原子重排检查。直接 Ewald；重排能量分项阈值为 `1e-6 kcal/mol`，力阈值为 `1e-5 kcal/(mol Å)`；每个模型选第一个非 C/H 原子的三个坐标方向做中心差分，步长 `1e-4 Å`，力导数误差阈值为 `5e-3 kcal/(mol Å)`，观测最大值为 `3.011e-3`。重排比较覆盖全部原子的力，最大误差 `1.954e-14 kcal/(mol Å)`，能量分项最大误差 `1.066e-13 kcal/mol`。这些是数值一致性检查，不是实验精度验证，也不是对全部原子逐一有限差分。
- PCFF 的 TMS 封端 PDMS、OPLS 的 PEG/PE/PEG 嵌段链、OPLS 的甲基封端 PA6，各用两条链、初始 0.3 g/cm³、种子 173，完成 200 步 NVT + 200 步 NPT + 200 步 production，步长 0.5 fs，共 300 fs。保留三阶段轨迹、能量、restart、参数来源与检查点；不标记为已平衡。

```bash
PYTHONPATH=.:src ISLAND_ASSIGNMENT_LAMMPS=/path/to/lmp ISLAND_PACKMOL=/path/to/packmol \
  python scripts/validate_j107.py work/j107-validation \
  --pcff-field work/pcff --opls-field work/opls-polymer
PYTHONPATH=.:src python scripts/check_opls_extended_sources.py \
  work/opls-polymer work/extended-sources.json
PYTHONPATH=.:src python -m pytest -q tests/test_phase_4j107.py
```

## 验证记录与证据边界

实现提交为 `dd5866b`；本报告与复现证据另行提交。

- 最终实现上的 J107/J106/J104 回归共 **90 passed**，其中 J107 为 31 项，覆盖缓存污染与来源变化、超时分母、独立入口、芳香氮端基、环扫描重试条件及结果一致性。
- 芳香氮修复后的构图/立体化学相关检查 **111 passed**；参考分配及 J94/J98/J103/J105 回归 **93 passed**。这两组在最后一次环扫描策略修改前执行，不能当作最终提交上重跑的结果。
- 测试组之间有重叠，不把各日志数字相加当作独立测试数量，也不声称整个仓库测试集全部通过。既有 J82 历史 README 字节断言的基线问题未由本阶段修复。
- 最终实现上实际执行的来源、能量/力与 MD 验证共 **11 个顶层检查组，全部通过**。短程 MD 仅证明指定模型在这些设置下完成流程，不证明长时间稳定、平衡或物性精度。

逐例结果保留原始错误、来源哈希、端基与运行设置；归档中的绝对路径记录当次运行位置，复现时使用上方命令指定本机路径。单独保存的测试文本日志仅去除行尾空白；压缩包内科学计算原始产物保持不变。证据校验命令见 [证据目录](evidence/phase_4j107/README.md)。

## 限制及后续

本轮覆盖率针对 PCFF 和 OPLS-AA 的指定配置，不包含 OPLS-UA，也没有扩大 GAFF/GAFF2 的网络参数化范围。来源绑定缓存是本地加速机制，不能为参数质量背书。预检有元素声明也不代表该元素的全部价态、化学环境或相互作用已覆盖。离子、高极性体系需要相容的电荷及长程相互作用策略；没有自动添加反离子或改变用户电荷态。

下一步应先对本轮“仅旧入口成功”的案例建立独立的小模型和环处理复现，区分完整接点/端基新增环境、分型遍历差异及 EMC 大型环图限制，再补来源支持的缺项。OPLS 的原子分型和电荷增量缺口应与扭转项分开处理；不得把大型图后端限制计作化学参数缺失。

J108 的通用反应模板、周期贯穿网络及 GAFF 网络接入仍是后续工作；J109 的长链规模化和 J110 的平衡/物性验收未由本阶段短程 MD 替代。应按本轮完整图失败的来源类型和频率继续补充经过独立验证的化学环境，避免用无来源通用参数提高表面成功率。
