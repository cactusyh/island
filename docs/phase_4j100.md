# J100 · DGEBA/PACM 的 OPLS-AA 周期固化

`opls-epoxy-amine-v1` 是显式选择的 EMC 2024 OPLS-AA 分型扩展。对取代环氧、芳基醚及脂肪伯/仲/叔胺扩展拓扑匹配，使用原生 EOXI、ANIS、EMET 与脂肪胺系列的类型和增量。参数文件、几何 LJ 混合及 1–4 的 LJ/库仑 0.5 缩放均来自同一发行库。下列证据适用于这个扩展的整图和周期固化入口。普通重复单元入口中的 PEG 仍可能因 EMC 的重复单元边界缺少环境而失败。

```bash
python -m island.forcefields.emc_profiles opls-epoxy-amine-v1 work/opls
python -m island.epoxy compile work/opls-library --field opls/2024/opls-aa --field-dir work/opls
python -m island.epoxy verify work/parity --library work/opls-library
python -m island.epoxy cure work/network --library work/opls-library --config recipe.json
```

化学计量用 DGEBA 4 个、PACM 2 个，沿用 J97 的种子、半径和热压协议。两套力场由相同的装填配方独立积分。同一构型上对 DGEBA、PACM、两种加成物逐个比较整图与局部参数后的 LAMMPS 能量分项和力；12 个构型通过。另按原始残基的预期值检查了每个 N 和 O 的电荷，伯/仲/叔胺氮分别为 −0.90/−0.78/−0.63 e。[参数文件指纹](evidence/phase_4j100/profile.json) · [原子电荷](evidence/phase_4j100/source-charges.json) · [逐项对照](evidence/phase_4j100/parity.json) · [周期不变性](evidence/phase_4j100/periodic.json)

三个相同配方和协议的种子 41、173、2026，OPLS 最终开环率分别为 75%、87.5%、87.5%，绕盒秩分别为 1、2、1。PCFF 在这三个种子上的终态转化率相同，两套模型的反应事件和最终结构各自独立。该比较没有检验体相性质相等。[逐批记录与 PCFF 对照](evidence/phase_4j100/cures.json) · [可复跑输入与轨迹](evidence/phase_4j100/opls-example.zip)

额外测试的连续 REACTER 路线在 OPLS 下出现了局部失稳。60,000 步协议在第 4,642 步报告缺角原子，缩短步长至 0.1 fs、起始距离上限至 3 Å 后在第 18,602 步仍出现同类错误。该入口目前限于经过验证的 PCFF。OPLS 使用上面的分批路线，显式转移氢坐标并做最小化松弛，再进入下一段 MD。

另外检查了相同醚环境在完整有限分子图中的迁移。二甲氧基乙烷 `COCCOC` 及羟基封端短醚 `CCOCCOCCOCCO` 在发行库报缺电荷增量，在扩展后均完成分配。乙氧基乙烷与苯甲醚在两个版本均成功，同批酰胺对照仍失败。PSMILES 的 `*CCO*` 走 EMC 拼链入口仍有缺项，表明完整分子图的新增覆盖与重复单元边界处理有区别。[逐个输入和原始类型对](evidence/phase_4j100/ether-challenge.json)
