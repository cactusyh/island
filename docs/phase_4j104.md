# J104 · OPLS-AA 聚合物分型扩展

`opls-polymer-v1` 基于固定版本 EMC 2024 OPLS-AA，在 J100 的环氧/胺规则基础上扩展醚、酯与 N-烷基酰胺环境。安装器核验原始 `.define`、`.prm`、`.top` 的 SHA-256。数值参数文件沿用原始字节，修改集中在分型环境。

源库 EACE、MACE 的酯烷氧碳使用 `c4oe`，NMAM、DMAM 的 N-烷基酰胺碳使用 `c4`。新增模板覆盖相应的较长或有支化的烷基取代。MACE 的原始氢模板 `H(C(O)(H)(H))` 同时匹配醚甲基；本扩展补入该残基的羰基环境，让 DMET、EMET 的醚氢使用原始残基指定的 `h1`。`h1` 与 `h1co` 的 LJ 参数有差异，因此这项检查同时验证类型和电荷。

七个源残基 DMET、EMET、MACE、EACE、ACAM、NMAM、DMAM 的全部原子类型及电荷均与 `.define` 中的独立定义匹配。另用较长的仲酰胺、叔酰胺、有支链的酯和醚作为独立验证分子，在扰动坐标上重排原子并重新分型，比较直接 Ewald 下的能量分项和逐原子力。[源文件指纹](evidence/phase_4j104/profile.json) · [化学与数值对照](evidence/phase_4j104/chemical-checks.json)

```bash
python -m island.forcefields.emc_profiles opls-polymer-v1 work/opls-polymer
python scripts/check_opls_polymer_profile.py work/opls-polymer work/checks --lammps lmp
python scripts/benchmark_opls_polymers.py docs/evidence/reference_assignment.csv.gz \
  work/opls-polymer work/replay.jsonl.gz --workers 4
```

J103 配方中的 `forcefield` 可以写成 `{"name":"opls","profile":"opls-polymer-v1"}`。该选择会在运行目录中安装并记录参数库。完整分子图与 EMC 重复单元接口分别衡量覆盖情况，分配通过后仍需要对周期装填和动力学进行验证。

## 固定 982 条 PSMILES 的完整复测

全部记录采用相同的氢封端、DP 4、每盒 2 条链、初始密度 0.9 g/cm³ 与种子 2026。下表的严格列只统计成功且未省略双键方向的输入。

| 参数库 / 显式分型扩展 | PolyJarvis 模式分配成功 | 严格立体化学成功 |
|---|---:|---:|
| EMC 原始 2024 OPLS-AA | 270/982（27.5%） | 262/982（26.7%） |
| J100 `opls-epoxy-amine-v1` | 335/982（34.1%） | 327/982（33.3%） |
| J104 `opls-polymer-v1` | **406/982（41.3%）** | **398/982（40.5%）** |

相对原始库增加 136 条成功记录，相对 J100 增加 71 条；两个基线的原有成功记录均无退化。数值系数库 `.prm` 和源残基 `.define` 保持原始 SHA-256。这个结果衡量分配覆盖率，不能证明新增的全部聚合物已通过动力学或物性验证。

剩余 576 条失败按诊断优先级归入：分型规则 413 条、电荷增量 83 条、扭转参数 71 条、角参数 5 条、总电荷检查 4 条。一条输入可能存在多个缺项，完整原始诊断保存在逐条结果中。当前结果不能外推为任意 PSMILES、任意端基或任意交联反应的覆盖率。

[统计与源文件指纹](evidence/phase_4j104/summary.json) · [J104 全部 982 条结果](evidence/phase_4j104/replay.jsonl.gz) · [J100 同条件复测](evidence/phase_4j104/epoxy-replay.jsonl.gz)

J104 自动化测试共 8 项通过，包括源残基全部原子的类型与电荷核对。独立能量/力检查验证同一参数化在原子重编号前后的一致性；它没有提供与实验或独立拟合参数的精度比较。
