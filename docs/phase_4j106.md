# J106：端基、立体化学和参数分配入口一致性

J106 从 `codex/j105-architecture-coverage` 的 `66cf1a2` 创建，分支为 `codex/j106-unified-polymer-graphs`。该父提交仅归档独立审查报告和复现证据；J105 实现仍为 `2bfd33f`，J104 为 `9a63a5a`。

## 修复内容

- 默认终端 dummy 被氢替换时，显式处理方括号原子的 `NoImplicit` 属性。`*C[C@H](C)*` 和 `*O[Si](C)(C)*` 不再因缺一个终端氢产生意外自由基。`add_hydrogens=False` 仍不额外添加默认封端 H 位点，重原子及分支锚点编号不变。默认氢封端并不表示每种聚合物的实验端基：例如酰胺重复单元的羰基端可能成为醛，用户应按实际化学指定端基。
- 共用构建器按连接点替换映射转移四面体邻居顺序和双键控制原子。重复单元、显式端基、支链的 E/Z 不再在入口转换时被丢弃。核心 `Bond` 保存可选 `stereo`、`stereo_atoms`，RDKit 往返、系统序列化、最终图身份、PCFF 图记录、多链打包均传递这两项。控制原子与键端点不相邻时拒绝输入。旧的无双键立体信息记录仍可读取。
- ETKDG 和局部模板构象都校验双键的三维构型。局部模板同时保留对应的双键控制关系。
- tacticity 仍采用 `final_graph_absolute_cip` 约定。端基导致首/末重复单元的指定中心成为非手性中心时，`stereochemical_sequence` 中该位置为 `null`，并记录 `terminal_achiral_site_ids`；完整请求单独保存在 `requested_stereochemical_sequence`。其余中心继续严格分配和验证，不会把链内配置失败当作端基例外。
- GAFF 字符串入口直接使用共用构建器，移除原先在力场层单独修补立体化学的链构建路径。保留 `atom_order="radonpy"` 用于原生遍历顺序对照。
- PCFF/OPLS 的字符串入口默认先构建完整图，再走与已有图相同的 `assign_emc_graph` 路径；支持自定义端基，返回 `EMCGraphAssignment`。本次未增加或修改力场数值参数。

## EMC API 迁移

```python
from island.forcefields.assignment import assign_polymer_parameters

# 默认：完整图分型，支持实际端基及保留的立体化学。
a = assign_polymer_parameters("*CCO*", "pcff", dp=3)

# 历史 EMC 重复单元建模：明确用于旧基准或 EMC 自己的多链打包。
b = assign_polymer_parameters(
    "*CCO*", "pcff", dp=3, emc_mode="repeat",
    nchains=2, density=0.9, stereo_policy="polyjarvis",
)
```

默认模式的 `site_groups` / `input_smiles` 替代旧结果特有的 `emc_smiles` 等字段；AA 与 UA 均按输入的氢封端完整图处理。旧模式保留 UA 甲基封端等原有约定。`nchains`、`density`、`stereo_policy` 属于旧模式，不能在默认模式中静默触发另一种构图行为。recipe 多链工作流仍使用自己的打包流程。

三个历史大样本 benchmark 脚本已显式使用 `emc_mode="repeat"`，保持 J98/J100/J104 的基准口径和历史收据不变；没有把旧模式的 982 条结果宣称为新的完整图覆盖率。

## 验证与证据

[证据目录](evidence/phase_4j106/)保存日志、引擎结果、完整 MD 文件压缩包及逐文件 SHA-256。

| 检查 | 结果与范围 |
| --- | --- |
| J106 新回归 | 51 项通过；DP 1/2/3/6、默认/显式 H 与自定义端基、非单调位点 ID、反向连接点标签、端部非手性中心、支链 E/Z、序列化、图身份、局部模板、三维反例和多链编号映射 |
| 构建器与原生参数对照 | 166 项通过，其中包含当时的 49 项 J106 测试；后增加的 2 项由上述 51 项日志覆盖。包含 GAFF 原子类型/参数与 RadonPy 对照、EMC 原生参数和旧模式测试 |
| 核心图、序列化、打包、局部模板 | 119 项通过，其中包含当时的 45 项 J106 测试；与上一行有重叠，不能相加当作独立测试数 |
| J94/J98/J101–J105 | 46 项通过 |
| J82 隔离测试 | 103 项通过，1 项历史 README 收据校验失败；同一项在未修改的 J105 `2bfd33f` 上同样失败，基线日志一并保存 |
| 真实 PCFF/OPLS 入口对照 | PEG、PA6、PET、PU × 2 个力场，共 8 组：原子类型、电荷、参数文件、完整输入 SMILES、位点映射一致；使用 `pcff-carboxamide-v1` 和 `opls-polymer-v1` |
| 真实 GAFF/GAFF2 入口对照 | 聚丙烯手性重复单元和含 E/Z 的二烯链 × 2 个力场，共 4 组：规范异构 SMILES 和原子类型多重集一致 |
| 实际周期性 MD | GAFF2 聚丙烯、PCFF cis 双键链、OPLS cis 双键链，均 DP 3、2 条链、初始 0.3 g/cm³；各完成 40 步 NVT、40 步 NPT、40 步 production，0.5 fs/步，共 60 fs，保留轨迹、能量、restart、参数哈希和检查点 |

MD 仅验证构图→分型→打包→实际积分的连接正确性，不是平衡或物性验证。首次验证脚本把完成状态拼成 `complete` 而接口为 `completed`；修正后在新目录重跑 15 组真实引擎检查，全部通过。归档的是修正后完整执行的结果。

重现需要 RDKit、RadonPy、EMC、支持所需势函数的 LAMMPS，以及 PACKMOL。验证使用的 Python/RDKit 版本见 `engine-results.json`，模拟程序和参数文件的哈希在压缩包的运行记录中。

```bash
PYTHONPATH=src python -m pytest -q tests/test_phase_4j106.py
PYTHONPATH=src ISLAND_ASSIGNMENT_LAMMPS=/path/to/lmp ISLAND_PACKMOL=/path/to/packmol \
  python scripts/validate_j106.py /tmp/new-j106-validation
```

J82 需要原有隔离启动器，不可直接与普通 pytest 混跑：

```bash
PYTHONPATH=.:src python -m scripts.run_phase_4j72_c1 --output /tmp/new-j106-j82 \
  -- -m pytest -q -p scripts.phase_4j73_scope tests/test_phase_4j82.py
```

## 当前边界

J106 修复入口和化学图保真问题，不证明已全面达到 RadonPy/PolyJarvis 的体系覆盖率。参数库缺失、价态限制、依赖 QM 的电荷计算仍可阻断分配；没有用零参数或通用替代项掩盖失败。未指定手性中心的 `*CC(C)*` 仍不能凭 `tacticity` 自动获得可靠的立构解释；多中心、多种重复单元的普适 meso/racemic 立构控制仍需后续实现。周期盒中有绕盒闭合路径的无限网络适配仍不是本阶段解决范围。旧版写出时已经丢失的 E/Z 信息也无法凭空恢复。

J93 原有工作保存在独立 worktree `../island-worktrees/j93-isolated-wip`、分支 `codex/j93-isolated-wip`，WIP 快照 `e5b6dc5`；四个原有已跟踪修改文件均与隔离前 SHA-256 一致。保留 stash `41a3c73e8dc99ff5cd461eb398ef66e5636041ef` 作为备份。该快照仅表示保存原工作，不宣称本轮完成或验证了 J93；未混入 J104/J105/J106。
