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
