# J108：可配置反应单体与周期 GAFF 网络

本阶段把 GAFF/GAFF2 接入已有的周期反应、动力学、检查点与 recipe 流程。每次成键后对当前完整显式氢图重新分配参数；保留原子顺序、速度、力场数据库哈希、RadonPy 版本与 1–4 缩放。环氧–胺和环氧–硫醇共用原有映射反应，不删除原子，氢从给体转移至环氧氧。

GAFF 网络电荷目前明确限定为 **Gasteiger**。分配不要求把绕行周期网络展开或嵌入为一个有限分子。每个连通分量的总电荷必须与形式电荷相符；源库、事件记录、阶段参数或状态文件变动会阻止续跑。

## 配置

```json
{
  "architecture": "crosslinked",
  "forcefield": {"name": "gaff2", "charge": "gasteiger"},
  "density_g_cm3": 0.3,
  "seed": 2026,
  "cure": {
    "dgeba": 2, "pacm": 1,
    "reaction_template": "epoxy_amine",
    "epoxy_smiles": "CC1CO1",
    "hardener_smiles": "NCCN",
    "cutoff": 4.5, "target": 0.25,
    "pre_steps": 20, "batch_steps": 20, "post_steps": 20, "max_idle": 1
  },
  "md": {"nvt_steps": 20, "npt_steps": 20, "production_steps": 20, "sample_every": 10, "segment_steps": 20}
}
```

这是小型功能测试输入，不是已经平衡的交联材料。`dgeba` 保留历史名称，表示所选环氧单体的分子数；`pacm` 表示所选胺固化剂数量。`reaction_template` 可显式选 `epoxy_amine` 或 `epoxy_thiol`，并与固化剂计数字段一致。选择硫醇反应时用 `petmp` 替代 `pacm`，并提供相应 `hardener_smiles`。不提供 SMILES 时保持原来的 DGEBA/PACM 或 DGEBA/PETMP。

自定义单体必须是单个已封端 SMILES 图，环氧需匹配现有末端环氧模板，胺需为支持的中性脂肪胺。酰胺氮不会当作胺参与反应。转化率与目标成键数按实际环氧和给体氢数计算，不再假定每个自定义环氧单体有两个反应位点。

自定义反应单体目前使用 GAFF/GAFF2；PCFF/OPLS 保持已有已编译参考环境库的范围。`cutoff` 必须小于当前最短盒长的一半。反应由几何搜索和规则驱动，不代表真实反应势垒或反应动力学速率。

## 验证

- J108 新增 10 项测试通过：两种反应 × 两种 GAFF 的原生拓扑/参数一致性，原子与形式电荷守恒，自定义官能团数量，来源变动拒绝，实际固化、续跑、事件/参数损坏拒绝。
- 人工构造的绕行周期环被正确分配；整体平移一个周期盒后 LAMMPS 力与能量不变。该模型验证周期处理，不作为合理网络构象的证明。
- 实际环氧–胺和环氧–硫醇短时固化通过；相关旧流程的 30 项集成与回归通过。

证据目录：`docs/evidence/phase_4j108/`。复现命令：

```bash
PYTHONPATH=.:src python -m pytest tests/test_phase_4j108.py -q
python -m island.workflows.recipe examples/recipes/j108_gaff2_custom_epoxy.json work/j108 \
  --lammps "$ISLAND_ASSIGNMENT_LAMMPS" --packmol "$ISLAND_PACKMOL"
```
