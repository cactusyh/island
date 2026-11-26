# J109：逐链构象、长度/序列分布与多组分

有限链 recipe 现在为每条链分别生成确定性随机种子的构象，再交给 PACKMOL 装箱。相同化学结构的链不再共用一个内部构象。已有均聚物、嵌段和支链入口继续使用同一最终化学图与参数分配路径。

新增 `mixture` 架构，`components` 中每个组分声明名称、有限链架构、重复单元和链数。可指定逐链 `lengths`，或给出离散 `length_distribution` 的 `values`/`probabilities`。两者不能同时给出，也不能与同一组分的 `dp`、`blocks`、`sequence` 同时指定。概率必须为正且和为一。

`copolymer` 组分可通过 `sequence_probabilities` 指定各重复单元概率，也可继续显式提供序列。随机选择结果会具体化为逐链 DP、序列和种子，保存到 `packing/chains.json`，便于复现与检查。每条链、每个原子的组分身份和链内编号保存在 `packing/site-map.json`。

示例见 `examples/recipes/j109_mixture.json`；所有组分使用 recipe 顶层指定的同一个力场。某一组分缺参数时整体分配失败，不静默替换力场。

```bash
python -m island.workflows.recipe examples/recipes/j109_mixture.json work/j109 \
  --lammps "$ISLAND_ASSIGNMENT_LAMMPS" --packmol "$ISLAND_PACKMOL"
PYTHONPATH=.:src python scripts/validate_j109_scale.py work/j109-scale
```

## 实测

- 新增测试连同 recipe/端基立体化学回归：**71 passed**；另有随机长度/随机序列混合物的真实构建、MD 与恢复测试 **1 passed**。
- 从 J105 三力场矩阵抽取完全相同的 18 个独立几何输入，保持原来的密度、链数与种子：**18/18 装箱成功**。该复测只验证几何，不冒充重新完成 54 个 MD 配方。
- 含不同链长 PE、PE/PEG 嵌段和 PE 支链的混合体系完成实际参数分配和周期 MD；映射、独立构象和恢复检查通过。
- DP 10、50、100 分别使用 seed 173/719、1/2 条链，共 **6/6** 组完成构建和 GAFF2 短时 MD。最大模型为 1,204 个显式原子。
- 每组为 20+20+20 步、0.5 fs，总 30 fs；密度 0.1 g/cm³。记录包含实际耗时、累计进程峰值 RSS、逐链回转半径、种子和输入。不是高密度熔体、长时间链松弛或性能扩展性的验收。

本阶段验证了 DP100 的这些 PE 实例，不保证任何 PSMILES/DP100 都能成功。长刚性链、高支化链和高密度装箱仍可能失败；保留四次确定性构象/装箱尝试及失败日志。未引入通用近似参数。

证据目录：`docs/evidence/phase_4j109/`。
