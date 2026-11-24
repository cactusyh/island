# J102 · 硫醇–环氧周期交联

除了 DGEBA/PACM 的两次胺加成，分批固化入口现在也支持 DGEBA 与四官能硫醇 PETMP（季戊四醇四(3-巯基丙酸酯)）。使用纯 DGEBA n=0，每 2 个 DGEBA 对 1 个 PETMP 提供相等的 4 个环氧和 4 个巯基。文献表 1 的 Neat (PETMP) 实验配方为每克 1.64 mmol DGEBA、0.82 mmol PETMP，实验还含催化剂。这里模拟净拓扑变化，未显式模拟催化剂、去质子化中间体或绝对反应时间。[公开配方及反应研究](https://doi.org/10.3390/polym7040680)

两套模板共用周期搜索、映射开环、氢转移、参数更新、松弛、绕盒分析及续算。模板指定硬化剂结构、供氢原子、官能度、事件字段和参考产物。PETMP 每个 S–H 可反应一次：切断环氧末端 C–O 与 S–H，生成 C–S 与 O–H，同一个氢迁移到氧。计数项为剩余巯基、已生成硫醚和巯基转化率。没有离去分子。验证单体和硫醚加成物时用原生整图 EMC PCFF 做三种构型上的逐项能量与逐原子力比较。[9 组单点对照](evidence/phase_4j102/parity.json)

在预试验中，5.5 Å 距离规则为种子 41 接受一个预期 C–S–C 角 169.7°、接触距离 5.43 Å 的事件，最小化结束时最大力分量仍为约 8.0 × 10⁵ kcal/(mol·Å)，该批次的 MD 失败。PETMP 的分批规则因此增加成键几何检查。中心取原生 PCFF 的 `c,sc,c` 平衡角 97.5°，默认半宽 45° 是本流程的选事件参数。事件记录保留预期角度，拒绝统计保留 `donor_angle`。可用 `donor_angle_tolerance` 调整半宽；比较不同规则时应分别报告转化率和终止条件。该规则只用于 PETMP。

```bash
python -m island.epoxy compile work/thiol-lib --chemistry petmp
python -m island.epoxy verify work/thiol-parity --library work/thiol-lib
echo '{"dgeba":4,"petmp":2,"seed":173,"cutoff":5.5,"max_idle":6}' > work/thiol.json
python -m island.epoxy cure work/thiol --library work/thiol-lib --config work/thiol.json --stop-after 2
python -m island.epoxy cure work/thiol --library work/thiol-lib --config work/thiol.json
```

PACM 的原有配置和检查点格式继续可用。PETMP 必须选对应的 PCFF 参数库，程序会拒绝不同化学之间的检查点或参数库混用。

两个 310 原子样品均完成固化后 NPT。种子 41 的环氧及巯基转化率均为 75%，绕盒秩为 1；种子 173 的两类转化率均为 100%，绕盒秩为 3。种子 173 在第 2 个完整阶段暂停再恢复，与单次运行的全部 20 个阶段记录及状态文件逐一相同。两组端点足以演示新化学的周期固化与续算，尚不足以推导热固性材料的体相性质。[事件与网络统计](evidence/phase_4j102/networks.json) · [重启对照](evidence/phase_4j102/replay.json) · [周期平移复核](evidence/phase_4j102/periodic.json) · [参数库、输入和轨迹](evidence/phase_4j102/thiol-example.zip)
