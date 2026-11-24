# J99 羧酰胺 PCFF 分型

EMC 随包 PCFF 的 FRC 明确区分酰胺共轭氮 `n` 和脂肪胺氮 `na`。旧面板中 232 条失败记录包含 `c_1/na` 的电荷增量错误。独立安装的 `pcff-carboxamide-v1` 对中性、非芳香且直接连到 C–C(=O) 的三配位氮增加分型模板，沿用全部原生 PCFF 系数。原生 `c_1/n` 增量为 0/0，`c/n` 为 +0.2108/−0.2108，`h*/n` 为 +0.4395/−0.4395。安装前校验两个源文件，保留原始 FRC 的 SHA-256。

```bash
python -m island.forcefields.emc_profiles pcff-carboxamide-v1 work/amide-field
python -m island.forcefields.capabilities probe pcff --repeat '*CC(C(=O)NC)*' --dp 4 --field-dir work/amide-field
python scripts/check_pcff_amides.py work/amide-field work/amide-check
python scripts/benchmark_pcff_amides.py docs/evidence/reference_assignment.csv.gz work/amide-field work/amide.jsonl.gz --workers 4
```

独立验证使用丙酰胺及取代物和含两个酰胺的分子，核对每个酰胺氮的电荷，并在同坐标、重排编号后重新分型，逐分项比较能量与力。[验证结果](evidence/phase_4j99/controls.json)
