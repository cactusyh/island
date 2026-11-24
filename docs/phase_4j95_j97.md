# J95–J97 · DGEBA/PACM 周期交联与性质验证

ISLAND 可以从未反应的 DGEBA（双酚 A 二缩水甘油醚，n=0）和 PACM（二氨基二环己基甲烷）出发，在正交周期盒内用 PCFF 固化并继续运行 LAMMPS。实测 7 个体系都形成至少沿一个方向贯穿周期盒的共价网络。下面的参数覆盖结论限于这组化学，贯穿结果限于给定盒尺寸、搜索半径和采样时长。[开发规划](crosslink_modeling_review_zh.md) · [环境与力场指纹](evidence/phase_4j95_j97/run-environment.json)

## J95 · 每次加成后的化学和参数

伯胺攻击环氧末端 CH₂ 后变成仲胺，可再攻击一个环氧形成叔胺。每次反应断开 C–O 和 N–H，形成 C–N 和 O–H。氧和转移的氢保留原子身份，各元素数量及形式总电荷守恒。记录中的原子编号从 0 开始，对应 LAMMPS 编号加 1。参考结构包括两个单体，以及在同一氮上分别完成一次和两次加成的完整产物。[化学定义](../src/island/epoxy/chemistry.py) · [化学与重启测试](../tests/test_phase_4j95_j97.py)

| EMC 9.4.4 力场 | 四个参考结构的实测结果 |
| --- | --- |
| PCFF | 4/4 可完整分配。使用 J94 原生 EMC 整图结果建立局部原子环境与相互作用对应关系，保留部分电荷、非键参数、键、角、二面角、improper 及 Class-II 交叉项。遇到未知环境会报错。 |
| 发行版 OPLS-AA 2024 | 0/4 可完整分配。DGEBA 缺少醚键电荷增量 `{c4, o2e}`，PACM 有连在氮上的氢缺分型，一次、二次加成物还分别涉及 `{n3h, c4}` 和 `{n3, c4}`。诊断包含源结构、EMC 位点和原子编号映射。[逐结构诊断](evidence/phase_4j95_j97/opls-coverage.json) |

四个模型各取松弛、扰动及扰动后随机重排原子编号三种情况。对同一组坐标分别用 EMC 原生输入和 ISLAND 输出计算键、角、二面角、improper、范德华、短程库仑及长程能量，并比较所有原子的力。12/12 组通过。最大能量分项差为 1.51 × 10⁻¹¹ kcal/mol 以内，最大力分量差为 1.53 × 10⁻¹⁰ kcal/(mol·Å) 以内。[逐构型结果](evidence/phase_4j95_j97/native-parity.json) · [参数映射实现](../src/island/epoxy/parameters.py)

## J96 · 固化、绕盒与续算

PACKMOL 周期装填后，先在 450 K、1 atm 做 NPT。每批按最小镜像 C···N 距离选取互不冲突的反应，可在已连接的团簇内成环。本轮半径上限 5.5 Å、每批最多两次加成、最小环长 6。更新拓扑、新接点参数及转移氢的位置，随后做能量最小化和 NPT 松弛。在新坐标上重新搜索，连续 6 批没有合格事件即停止，最后再做 NPT。批间保留速度并重新建立热浴。批次数和模拟时长不能用于推算实验反应速率。[固化配置](../src/island/epoxy/cure.py) · [周期搜索和统计](../src/island/epoxy/network.py)

每个阶段保存坐标、速度、盒尺寸、随机状态、候选/拒绝数、事件和热力学量。事件含 C/O/N/H 编号、初始接触距离、成键镜像及成键时的团簇关系。另存本帧非零键镜像（其余为零），沿共价环累加镜像位移得到绕盒向量及其秩。对一个已贯穿且有可反应候选的实际快照整体平移、重新包裹后，候选和绕盒秩一致；用直接 Ewald 比较，最大能量分项差 2.20 × 10⁻¹¹ kcal/mol，最大力分量差 1.27 × 10⁻¹⁰ kcal/(mol·Å) 以内。中断恢复测试在完成加成后重启，最终坐标、盒、速度逐元素完全相等，事件和热力学记录也一致。[周期不变性结果](evidence/phase_4j95_j97/periodic-invariance.json) · [续算测试](../tests/test_phase_4j95_j97.py)

部分 LAMMPS 日志报告 `Inconsistent image flags`。输入使用周期坐标，存在非零绕盒环时，单个原子的一组镜像编号无法同时展开环上所有键。分析局部键向量和贯穿性时应使用对应阶段的 `bond-images.json` 和盒尺寸；整网的有限展开坐标不适用于这类环。

同一套化学还生成 LAMMPS REACTER 的两组模板与原子映射。由相同的未固化压实快照另起 `fix bond/react` 连续轨迹，在 4 Å 起始距离上限下积分 60,000 步（12 ps），观察到 3 次伯胺加成和 3 次仲胺加成，环氧转化率 75%，绕盒秩 1。用已发生的事件在 ISLAND 重建拓扑，与 LAMMPS 原生反应产物在同一终态坐标比较，最大能量分项差 2.11 × 10⁻¹¹ kcal/mol 以内，最大力分量差 1.12 × 10⁻¹⁰ kcal/(mol·Å) 以内。[连续固化结果](evidence/phase_4j95_j97/continuous.json) · [两种后端的输入和轨迹](evidence/phase_4j95_j97/pcff-cure-example.zip)

## J97 · 转化率与网络结构

化学计量为每 2 个 DGEBA 配 1 个 PACM，分别提供 4 个环氧和 4 个胺活泼氢。化学计量组两类转化率相等。表中首次贯穿是第一个绕盒秩大于零的已保存帧，单批可以跨过多个转化率。每个体系都完成固化后 NPT。

| DGEBA + PACM | 原子数 | 种子 | 首次观测贯穿时环氧转化率 | 终态环氧 / 胺氢转化率 | 终态绕盒秩 | 最大团簇原子占比 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 4 + 2 | 278 | 41 | 50% | 75% / 75% | 1 | 100% |
| 4 + 2 | 278 | 173 | 75% | 87.5% / 87.5% | 2 | 100% |
| 4 + 2 | 278 | 2026 | 87.5% | 87.5% / 87.5% | 1 | 100% |
| 8 + 4 | 556 | 41 | 68.75% | 75% / 75% | 1 | 67.6% |
| 8 + 4 | 556 | 173 | 62.5% | 93.75% / 93.75% | 2 | 100% |
| 8 + 4 | 556 | 2026 | 81.25% | 81.25% / 81.25% | 1 | 100% |
| 4 + 3（胺过量） | 319 | 41 | 75% | 100% / 66.67% | 2 | 100% |

6 个化学计量样品均因在固定半径和等待时间内找不到合格接触而终止。胺过量组消耗了全部 8 个环氧及 12 个活泼氢中的 8 个，终态含 1 个伯胺、2 个仲胺和 3 个叔胺。共保存 109 个阶段，其中 95 个为固化批次。统计文件逐批列出单体连接图的环秩、度为 1 的悬挂单体数及链段长度（以单体间键数计），并列出按转化率统计贯穿概率的观测数、实际采样转化率及提前停滞数。[逐体系和逐阶段数据](evidence/phase_4j95_j97/cure-study.json)

两个尺寸在 75% 环氧转化率时都观察到 2/3 的种子贯穿。HTPolyNet 公开算例在 51.0% 帧未见周期偏移冲突、69.0% 帧已见冲突；理想等反应性无环 A₂+B₄ 模型给出约 57.7%。本次 6 个化学计量样品首次观测贯穿的范围为 50%–87.5%。盒尺寸、批次间隔、反应准则和力场均会影响这个观测值，三个种子不足以外推宏观凝胶点。[HTPolyNet 周期网络示例](https://htpolynet.readthedocs.io/en/latest/user-guide/periodic-networks.html)

## J97 · 密度、体积收缩、玻璃化与拉伸

将同一装填样品固化前、固化后的快照分别在 300 K、1 atm 做 6 ps NPT，取后半段体积均值，计算 `100 × (1 − 固化后平均体积 / 固化前平均体积)`。下表为每种尺寸 3 个种子的均值 ± 种子间样本标准差。

| 配方 | 固化后密度 (g/cm³) | 体积收缩率 |
| --- | ---: | ---: |
| 4 + 2 | 0.967 ± 0.017 | −0.08 ± 4.33% |
| 8 + 4 | 0.982 ± 0.043 | 1.70 ± 2.90% |

商业 Epon 828/PACM 致密样品的文献密度约为 1.2 g/cm³，高于两组模拟均值。含 0.1–1.0 wt% 颜料的商业 YD-128/PACM 样品报告 0.66–1.12% 的 TMA 表观收缩。两者在原料、热历史或量测定义上与模拟有差异。目前缺少与这里纯 n=0 单体、转化率和温压历史都匹配的净树脂收缩对照；模拟中种子间离散度也超过平均收缩量。[密度参考](https://doi.org/10.1021/acsomega.0c04035) · [收缩参考，表 4](https://doi.org/10.3390/polym14051028)


8 + 4 配方的种子 173 达到 93.75% 环氧转化率。把固化前后各自的 NPT 延长到 30 ps，300 K 下对应平均密度为 1.016、1.021 g/cm³，配对体积收缩为 0.485%。这是同一种子的采样时长对照。PACM 商业异构体比例和商业树脂的低聚物分布也会影响实验对比。

| 性质 | 实际协议 | 结果 |
| --- | --- | --- |
| 表观 Tg | 对种子 173 从 550 K 冷至 250 K，用比容–温度曲线拟合连续双斜率。要求高温斜率大于正的低温斜率、折点位于搜索区间内部，且相对直线的 BIC 改善大于 6。 | 4 + 2 在 10、20 K/ps 下均未得到受支持的折点。8 + 4 在 4 K/ps 下得到 442.75 K（BIC 改善 7.94）；在 2 K/ps 下最佳折点顶到搜索上界 490 K，按同一准则拒绝。 |
| 方向拉伸模量 | 对上述两个样品在 300 K 分别沿 x、y、z 拉伸至 4%，控制横向压力。取 0.1%–3% 应变区间分块回归，要求斜率为正且 R² ≥ 0.8。 | 4 + 2 的两种应变速率 0.002、0.005 ps⁻¹ 下，六个拟合均未通过。8 + 4 在 0.0005 ps⁻¹ 下仅 z 向通过（10.41 GPa，R² = 0.857）；在 0.001 ps⁻¹ 下仅 x 向通过（4.25 GPa，R² = 0.803）。 |

文献中商业 DGEBA/PACM 的 DMA Tg 为 175 °C（448.15 K），25 °C 储能模量为 2.28 GPa，拉伸杨氏模量为 2.60 ± 0.30 GPa。储能模量与拉伸模量是不同观测量。文献树脂的环氧当量约为 188 g/eq，本算例采用纯 n=0 单体。一个冷却速度下的数值接近，其中另一个速度的拟合被拒绝；方向拉伸结果又随速率和方向变化。这些轨迹验证了测量流程，目前不足以报告材料的收敛 Tg、各向同性杨氏模量或实验精度的收缩率。[文献配方与力学数据，表 2](https://pmc.ncbi.nlm.nih.gov/articles/PMC7563135/) · [全部拟合及体积结果](evidence/phase_4j95_j97/properties.json) · [34 组测量的原始数据和可运行输入](evidence/phase_4j95_j97/property-samples.zip)

## 复现

需要 Python 3.11 或更高、Perl、支持周期盒的 PACKMOL，以及包含 MOLECULE、CLASS2 和 KSPACE 的 LAMMPS。连续反应还需要带 REACTION 的 LAMMPS。下例假定可执行文件分别叫作 `packmol`、`lmp` 和 `lmp-reaction`，可通过命令行选项改成实际路径。验证所用版本分别为 PACKMOL 21.1.1、LAMMPS 29 Aug 2024 和 22 Jul 2025 Update 4。

```bash
python -m pip install -e '.[epoxy]'
python -m island.epoxy compile work/library
python -m island.epoxy coverage work/opls
python -m island.epoxy verify work/parity --library work/library
python -m island.epoxy study work/study --library work/library
python -m island.epoxy verify-periodic work/pbc --library work/library --cure work/study/d4-p2-s173
python -m island.epoxy reacter work/live --library work/library --cure work/study/d4-p2-s173 --lammps lmp-reaction --steps 60000
python -m island.epoxy characterize work/properties --study work/study --extended
```

只运行单个体系时，下面的命令在两个完整 MD 阶段后暂停。第二次运行会从最近的完整阶段恢复。输出目录会校验配方和参数库指纹。

```bash
mkdir -p work
echo '{"dgeba":4,"pacm":2,"seed":173,"cutoff":5.5,"max_idle":6}' > work/recipe.json
python -m island.epoxy cure work/network --library work/library --config work/recipe.json --stop-after 2
python -m island.epoxy cure work/network --library work/library --config work/recipe.json
```

代表性归档中 `library/` 包含四个原生参考结构及编译参数，`cure/` 包含 4 + 2、种子 173 的全部批次、轨迹及重启快照，`reacter/` 包含对应的连续轨迹、反应模板及两套终态单点输入。性质归档中每个目录都可用 `lmp -in property.in` 独立重跑，`samples.dat` 包含每 100 步的温度、密度、体积、三个盒边长及三个对角压力分量。[固化归档](evidence/phase_4j95_j97/pcff-cure-example.zip) · [性质归档](evidence/phase_4j95_j97/property-samples.zip)

本轮执行 J95–J97、J94、参考力场分配、J20 和 J31 的相关回归，共 143 项通过。两次执行固化矩阵的 109 份坐标、盒与速度检查点逐份相同。两份 ZIP 均通过 CRC 检查，SHA-256 记录可用于下载后验证。[测试范围与归档校验](evidence/phase_4j95_j97/verification.json)
