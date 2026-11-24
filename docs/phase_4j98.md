# J98 · 查询单个结构的力场覆盖情况

已有 PSMILES 时，用 `--repeat` 和 `--dp` 进行实际参数分配。已有嵌段或支链聚合物的完整连接图时，用 `--mol` 读入保留键序、立体化学和显式氢的 MDL Mol 文件，或在 Python 中调用 `probe(mol, field)`。成功结果记录所用参数库的 SHA-256；失败结果保留原生诊断。EMC 能定位分型失败原子时，整图结果会给出对应的源位点 ID。GAFF 系列探针使用 Gasteiger 电荷。探针通过后，单点能、动力学和材料性质仍分别标为 `unmeasured`。

```bash
python -m island.forcefields.capabilities probe pcff --repeat '*CCO*' --dp 4
python -m island.forcefields.capabilities probe gaff2 --repeat '*CCO*' --dp 4
python -m island.forcefields.capabilities probe opls --mol branched.mol
python -m island.forcefields.capabilities baseline docs/evidence
```

`baseline` 汇总仓库中固定版本的记录。GAFF 和 GAFF2 的 RadonPy 对照各为 1,077/1,077。另一组 982 个 PSMILES 上，随包 PCFF、随包 2024 OPLS-AA 和 PolyJarvis 修订的 2024 OPLS-AA 在允许省略双键方向时分别有 678、270、205 个成功，要求保留指定双键方向时分别为 646、262、205。PolyJarvis 公开的 611 个成功“输入 × 力场”组合均可复现。已做过 MD 的 12 个固定图案例包括 7 个嵌段、3 个支链、2 个支链嵌段。聚合度、力场版本、立体化学策略和验证层次均影响这些数字的适用范围。[数据来源与逐项对照](reference_assignment_zh.md) · [固定图动力学](phase_4j94.md)
