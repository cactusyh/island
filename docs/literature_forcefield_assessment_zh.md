# island_ffref 逐篇正文、SI 与力场归属评估

这批资料中有多篇实质性参数来源。最适合 OPLS 扩展的是电解质、离子液体、特定 API 和磷酸二甲酯等窄域；Class-II 来源主要是 CFF93/COMPASS 的特定聚合物/官能团扩展，包括补读的 Sun 1995 聚硅烷；Smith 2004 和 PolCA 2021 则应保持为各自独立的有机硅模型。**适合作为参数来源，不表示已成为可执行 ISLAND profile。**

本报告按 29 个来源组记录，覆盖当前目录的 **55 个文件：50 个 PDF、4 个 ZIP、1 个 XLSX**；其中有 27 个论文正文/SI 来源组、1 组 LigParGen 工作簿和 1 个 PCFF 工程归档。29 组均已完成本次阅读评估；组 29 的完整 PDF 补读见下文，同组截断 checkpoint 单独标记，不用于参数结论。相同论文的正文、SI、附件和重复副本合并分析，每个文件仍单独记录哈希。原 J44 的 52 个文件全部保持原哈希；新增的 Smith PDF checkpoint 与正文完全相同。27 个论文组包含只有原始 OPLS SI 的组，并不表示 27 篇都有本地正文。

这里同时检查 **SI（Supporting Information，补充材料）** 与 **Si（硅）**。阅读重点是正文势函数/方法/参数表及所供 SI、ZIP 中的实际参数与 convention。只使用本地文件，未下载论文、未网络补参。缺少的 SI、引用基础库、DOI 或许可证据均明确保留缺口；不能把下载水印视为再分发授权。

本次是文献归属评估，未做参数导入、引擎实现或数值资格验证。没有修改历史 PCFF/OPLS/GAFF 数据、J17–J45 证据或已注册 profile，也没有运行 MD。下文是新增阅读结论，不重写 J44/J45 的历史收据。`production_validated=False`，`simulation_readiness="not_established"` 保持不变。

## 最值得推进的路径

| 目标 | 优先来源 | 应加入的路径 | 尚不能省略的工作 |
| --- | --- | --- | --- |
| 全原子 PDMS 专用模型 | Smith 2004，组 23 | 独立 external PDMS Buckingham 模型 | 特殊弯曲–扭转、柔性/约束版、端基与空白系数语义 |
| Class-II 硅氧烷 | Sun–Rigby 1997，组 17 | external CFF/PCFF-related Class-II | 同版基础行、耦合/省略规则、类型与非键约定 |
| Class-II 聚硅烷 / Si–H、Si–C、Si–Si | Sun 1995，组 29 | external PCFF/CFF93 polysilane Class-II | 同版烷烃/苯依赖、空白 k4、混合/1–4、逐图耦合；不提供 Si–O–Si |
| OPLS-compatible 甲基 PDMS | Xiang 2025 SI S4，组 24 | 外部 OPLS-compatible PDMS | 保留 J45 的单位、端基、电荷、特殊对等 blocked 缺口 |
| 有机硅液相联合原子 | Jorge 2021 PolCA，组 25 | 独立 external PolCA-UA | UA 映射、约束、RB、无 1–4、Si 电荷算法与极化后处理 |
| PEO/LiTFSI 电解质 | Liang 2024 SI S4–S7，组 04 | 外部 OPLS-compatible 电解质 | 精确端基、pair override、基础 convention 和独立数值复现 |
| 机器可读 OPLS 窄域 | Lopes 2025 ZIP/SI，组 12 | 外部 OPLS-based 特定 API | 库/电荷模型映射；有 X-site 的域另验 |
| 有机磷 | Murzyn 2013 表 4–5，组 15 | 外部 OPLS-compatible DMP− 或 DMPH | 分开离子/中性状态、基础行和特殊对 |
| 磷腈 | Sun–Ren–Fried 1998 附录，组 16 | 外部 COMPASS Class-II 磷腈 | 状态特异分型、最终 increment、基础行与全部耦合 |
| Class-II 简单控制域 | Sun 1998 附录 B，组 18 | external published COMPASS 烷烃/苯 | 逐项耦合与独立能量/力验证 |

没有一篇可以不经逐图参数完备与数值审核就直接并入通用 native PCFF/OPLS-AA。上述优先级按参数内容和归属确定，没有声明新的“全目录最完整可执行域”。尤其不能因更适合 OPLS 的简单函数就把别篇 Class-II 降阶，也不能因同用 LJ12–6 就把 PolCA、额外站点模型、12–6–4 修正当成普通 OPLS。

## 逐组阅读记录

页码均指本地 PDF 的一基页码，必要时补印刷页码；SI 缺失指当前目录未供，不表示论文从未发表 SI。完整路径/哈希/书目和已有许可证据见 [文件与评估清单](evidence/literature_review_2026_10_27/source-review.json)。

### 01. Revised Atomic Charges for OPLS Force Field Model of Poly(ethylene oxide): Benchmarks and Applications in Polymer Electrolyte

Chan-En Fang; Yi-Chen Tsai; Christoph Scheurer; Chi-Cheng Chiu；2021；Polymers；DOI `10.3390/polym13071131`

本地文件：

- [elctrolytes/.ipynb_checkpoints/polymers-13-01131-checkpoint.pdf](../island_ffref/elctrolytes/.ipynb_checkpoints/polymers-13-01131-checkpoint.pdf)
- [elctrolytes/polymers-13-01131-s001.pdf](../island_ffref/elctrolytes/polymers-13-01131-s001.pdf)
- [elctrolytes/polymers-13-01131.pdf](../island_ffref/elctrolytes/polymers-13-01131.pdf)
- [opls/polymers-13-01131.pdf](../island_ffref/opls/polymers-13-01131.pdf)

**适合归属：外部 OPLS-compatible 的 PEO/PEODME 电荷扩展。** 正文 PDF p3 式 1–3 给出谐键、谐角、OPLS torsion 和 LJ/Coulomb；p4 图 1 直接包含原始 OPLS、Barbosa 和 OPLSR 的逐原子电荷。5 页 SI 的表 S1–S3 是体系组成和时长，图 S1–S8 是 MSD，不是完整成键参数库。主要贡献是 B3LYP/6-311G* 多构象 RESP 电荷与验证，其他项保留 OPLS。作者另定义聚合物/离子分别为 0.8/0.55 的缩放模型，必须另立版本，不能静默更改现有电荷。

**Si 与导入边界：** 无硅氧烷参数。还须固定所依赖的 OPLS、TFSI 参数版本、端基电荷映射及特殊对约定。图中电荷确实存在，不能因为不在文本表中就判为缺失。优先级较高，适用于聚醚电解质。

### 02. Differentiable OPLS Force Field Parameterization for Ionic Electrolytes and High-Throughput Application to Lithium-ion Batteries

Haichao Huang; Zilin Chen; Qi Liu; Tianqi Zhao; Yunpei Liu; Guotao Qiu; Jianhui Chen; Zhen Li; Wenshuo Liang; Minsung Cho; Manxue Zhang; Feiyu Kang; Xiaolong Zou; Yidan Cao; Xushan Zhao; Ziqi Cheng; Ye Mei；2026；arXiv preprint；DOI `10.48550/arXiv.2607.04633`

本地文件：

- [elctrolytes/2607.04633v2.pdf](../island_ffref/elctrolytes/2607.04633v2.pdf)

**适合归属：外部优化 OPLS 电解液模型及可微参数化方法。** 正文 p8–10 描述 68 类、136 个 LJ 优化变量、LigParGen 1.14×CM1A 起点和盐电荷因子 0.7。p24/p28 说明 BF4−/ClO4− 优化存在局限，部分生产体系保留原始参数。该 46 页 PDF 在 p30 起内嵌 SI：p31–32 是类型/优化细节，p33–39 是收敛、平台对比、组分图。**136 个变量不等于提供了 136 个最终系数**；本地文件未给出可重建完整模型的最终系数总表/发布件。

**Si 与导入边界：** p8 与内嵌 SI p32 提到含硅/埋藏 Si 沿用原 OPLS；这不提供独立的 Si–O–Si/PDMS 完备证据。其筛选保留 PF6−，不支持普遍有机磷。需作者最终参数文件、完整类型规则、逐分子电荷；不能把作者的类比转移算法当作已有来源行。当前更适合方法参考。

### 03. Force Field Development and MD Simulations of Poly(ethylene oxide)/LiBF4 Polymer Electrolytes

Oleg Borodin; Grant D. Smith; Richard Douglas；2003；Journal of Physical Chemistry B；DOI `10.1021/jp027539z`

本地文件：

- [elctrolytes/J-Phys-Chem-B_2003_107_6824_Borodin_PEOLiBF4ff.pdf](../island_ffref/elctrolytes/J-Phys-Chem-B_2003_107_6824_Borodin_PEOLiBF4ff.pdf)

**适合归属：独立外部 PEO/LiBF4 多体极化模型，或与其分开的有效两体模型。** 正文 p4 表 2 给电荷/极化率，表 3 给排斥和色散参数，同页正文给 B–F 键及 F–B–F 角；p6 表 5 给 POM 电荷/极化率。MB/TB 与 FF-1/FF-3 的区别见 p2–4。没有本地单独 SI。

**Si 与导入边界：** 无 Si。含显式极化、指数排斥及有效两体修正，不能替换为普通 LJ12–6。Li 的零极化率/色散是作者明确规定的值，区别于缺行补零。PEO 基础库和完整势函数依赖前两篇系列论文；需独立极化后端和同版基础来源，不宜直接补到固定电荷 OPLS/PCFF。

### 04. Charge Scaling in Classical Force Fields for Lithium Ions in Polymers

Dongyue Liang; Yuxi Chen; Chuting Deng; Juan J. de Pablo；2024；ACS Macro Letters；DOI `10.1021/acsmacrolett.4c00368`

本地文件：

- [elctrolytes/mz4c00368.pdf](../island_ffref/elctrolytes/mz4c00368.pdf)
- [elctrolytes/mz4c00368_si_001.pdf](../island_ffref/elctrolytes/mz4c00368_si_001.pdf)

**适合归属：优先整理外部 OPLS-compatible PEO/LiTFSI 模型。** 正文 p3 表 2 给盐电荷因子 **0.79**，Li–O(PEO) 的 epsilon/sigma 为 **0.06 kcal/mol、2.39 Å**，Li–O(TFSI) 为 **0.07 kcal/mol、2.31 Å**。13 页 SI 的 p7 表 S4/S5 给非键/电荷和键，p8 表 S6/S7 给角及 K1–K4 扭转；p2 给 CP2K、12 Å cutoff 与 SPME。这是实质性参数来源，不能归为只有缩放说明。

**Si 与导入边界：** 无 Si，含 Li、TFSI 的 S/F。修订只缩放盐电荷，保留 PEO 电荷；S4 离子电荷为缩放前值。两个 Li–O 是显式 pair override，不能被普通混合规则覆盖。仍须复核基础 OPLS 的 1–4、混合/排除、谐势 1/2 因子、端基与原子映射。不得与 OPLSR 的 0.8/0.55 混用。这是本地 OPLS 电解质的高优先候选。

### 05. Accelerating amorphous polymer electrolyte screening by learning to reduce errors in molecular dynamics simulated properties

Tian Xie; Arthur France-Lanord; Yanming Wang; Jeffrey Lopez; Michael A. Stolberg; Megan Hill; Graham Michael Leverick; Rafael Gomez-Bombarelli; Jeremiah A. Johnson; Yang Shao-Horn; Jeffrey C. Grossman；2022；Nature Communications；DOI `10.1038/s41467-022-30994-1`

本地文件：

- [elctrolytes/s41467-022-30994-1.pdf](../island_ffref/elctrolytes/s41467-022-30994-1.pdf)

**适合归属：PCFF+ 聚合物电解质筛选/误差学习工作流与基准。** 正文 Methods（PDF p8）采用 PCFF+，TFSI 电荷因子 0.7，并引用 Supplementary Table 7；该 SI 未在本地提供。重点是短/长轨迹误差学习及大量聚合物筛选，没有独立可重建的 PCFF+ 数值库。

**Si 与导入边界：** 未给可导入的 Si 行。PCFF+ 名称不能证明历史 pinned PCFF 完整，也不能由模拟成功反推所有源参数。适合借鉴筛选/基准流程；补参需要其引用的 SI/实际 FF 发布件。

### 06. Development of PolyParGen Software to Facilitate the Determination of Molecular Dynamics Simulation Parameters for Polymers

Makoto Yabe; Kazuki Mori; Kazuyoshi Ueda; Minoru Takeda；2019；Journal of Computer Chemistry, Japan - International Edition；DOI `10.2477/jccjie.2018-0034`

本地文件：

- [engineering/5_2018-0034.pdf](../island_ffref/engineering/5_2018-0034.pdf)
- [opls/5_2018-0034.pdf](../island_ffref/opls/5_2018-0034.pdf)

**适合归属：OPLS/LigParGen 聚合物分片参数化的工程参考。** 正文 p2 实际实现使用本地 LigParGen，把大图分片、加 H 封端、重叠合并；p3 明确将残余电荷平均分配到原子。两份 PDF 是同文，未供额外 SI。论文不是一套独立 PCFF/OPLS 系数库。

**Si 与导入边界：** 没有证明 Si–O–Si 覆盖。分片思想可参考，但均摊残余电荷不能进入当前禁止电荷重分配的来源路径。需实际生成件、基础库/QM/电荷来源，且在最终图上审查拼接边界。不能仅凭自动工具能运行便宣称参数完备。

### 07. LAMMPS-pcff-main.zip：PCFF 工程归档

工程资产，无论文书目。

本地文件：

- [engineering/LAMMPS-pcff-main.zip](../island_ffref/engineering/LAMMPS-pcff-main.zip)

**适合归属：PCFF Class-II 导出/等价映射工程对照；这不是论文。** ZIP 有 README.md、pcff.frc，没有 LICENSE。FRC 的 SHA-256 是 `e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c`，与既有 pinned PCFF 一致，没有新增来源版本。

**Si 与导入边界：** README 的 PIB 类型、电荷增量和 coefficient 分区示例可供工程核对；其中缺 improper 置零不能采纳。相同 FRC 重新打包不补全混合 Si/P 耦合，README 也不能取代化学/参数来源证据。

### 08. MolMod – an open access database of force fields for molecular simulations of fluids

Simon Stephan; Martin T. Horsch; Jadran Vrabec; Hans Hasse；2019；Molecular Simulation；DOI `10.1080/08927022.2019.1601191`

本地文件：

- [engineering/MolMod   an open access database of force fields for molecular simulations of fluids.pdf](../island_ffref/engineering/MolMod%20%20%20an%20open%20access%20database%20of%20force%20fields%20for%20molecular%20simulations%20of%20fluids.pdf)

**适合归属：将来独立的刚体/多位点流体模型接口。** 正文 p3 的模型是刚性单/多位点 LJ12–6、点电荷、偶极与四极；p6 表 1 是物质/文献目录，不是所有模型的系数表，p4–5 给形式/跨引擎表示。没有本地具体模型导出件或独立 SI。

**Si 与导入边界：** 本地这篇不提供可直接重建的 PDMS 参数集。不能把数据库概览变成通用 OPLS/PCFF 库；选定模型必须连同刚体几何、多极矩、混合规则一并取得。适合数据库/溯源设计参考。

### 09. Type Label Framework for Bonded Force Fields in LAMMPS

Jacob R. Gissinger; Ilia Nikiforov; Yaser Afshar; Brendon Waters; Moon-ki Choi; Daniel S. Karls; Alexander Stukowski; Wonpil Im; Hendrik Heinz; Axel Kohlmeyer; Ellad B. Tadmor；2024；Journal of Physical Chemistry B；DOI `10.1021/acs.jpcb.3c08419`

本地文件：

- [engineering/jp3c08419.pdf](../island_ffref/engineering/jp3c08419.pdf)
- [engineering/jp3c08419_si_001.pdf](../island_ffref/engineering/jp3c08419_si_001.pdf)

**适合归属：类型标签、拓扑映射与 LAMMPS 导出的工程改进。** 正文 p2–5 介绍 bonded 形式和类型语义；7 页 SI 给丙烷 VMD/TopoTools 脚本、质量、电荷与拓扑。SI p1 明确说力场系数要由输入文件或 OpenKIM 另供。它有数值示例，但不是新开发的参数库。

**Si 与导入边界：** 正文讨论 SiO2 等基准不构成 PDMS 参数发布。type label 不能代替最终图分型，也不能赋予执行权限。工程价值高，新增参数价值有限。

### 10. Atomistic Polymer Modeling: Recent Advances and Challenges in Building and Parametrization Workflows

Hannah N. Turney; Micaela Matta；2025；Macromolecules；DOI `10.1021/acs.macromol.5c01166`

本地文件：

- [engineering/ma5c01166.pdf](../island_ffref/engineering/ma5c01166.pdf)
- [engineering/ma5c01166_si_001.pdf](../island_ffref/engineering/ma5c01166_si_001.pdf)

**适合归属：PSMILES、建链、参数化与 FAIR 溯源的方法参考。** 正文 p3 表 1 比较工具，p5–7 讨论电荷、torsion、ML 和直接化学感知。9 页 SI 是检索准则、聚合物类别和 82 篇调查文献，非系数表。

**Si 与导入边界：** SI 中的 silicone 检索词不是含硅参数。不能把综述列举的 OPLS、GAFF、ML 等路径拼成一套力场。适合作为设计/文献索引，不作为补参源。

### 11. Development of OPLS-AA Force Field Parameters for 68 Unique Ionic Liquids

Somisetti V. Sambasivarao; Orlando Acevedo；2009；Journal of Chemical Theory and Computation；DOI `10.1021/ct900009a`

本地文件：

- [halogen/ct900009a.pdf](../island_ffref/halogen/ct900009a.pdf)
- [halogen/ct900009a_si_001.pdf](../island_ffref/halogen/ct900009a_si_001.pdf)
- [opls/ct900009a.pdf](../island_ffref/opls/ct900009a.pdf)
- [p/ct900009a.pdf](../island_ffref/p/ct900009a.pdf)

**适合归属：外部 OPLS-compatible 离子液体扩展。** 正文 p2 给几何混合和 1–4 非键减半；p4–6 表 1–5 给电荷/LJ/torsion，图 3–4 给类型/电荷映射。15 页 SI 的 p3–6 表 S1–S5 补键、角与扭转，后面为构象能量和盒尺寸。确实包含可提取参数。

**Si 与导入边界：** 无 PDMS。PF6− 是明确的 P 来源，但只对应六氟磷酸根，不能转用于磷酸酯、磷腈或氧化膦。应选择一致的通用或 EMIM/BMIM、EPyr/BPyr 专用 charge/torsion 版本，核对基础 OPLS 与阴离子刚性/柔性处理；非键表存在不证明柔性拓扑全部齐全。离子液体窄域优先级高。

### 12. An all-atom force field for MD simulations on organosulfur and organohalogen active pharmaceutical ingredients developed from experimental sublimation enthalpies and single crystal X-ray diffraction data

Cátia S. D. Lopes; Manuel E. Minas da Piedade; Carlos E. S. Bernardes；2025；Physical Chemistry Chemical Physics；DOI `10.1039/D5CP01216C`

本地文件：

- [halogen/d5cp01216c.pdf](../island_ffref/halogen/d5cp01216c.pdf)
- [halogen/d5cp01216c1_suppl.zip](../island_ffref/halogen/d5cp01216c1_suppl.zip)
- [halogen/d5cp01216c2_suppl.pdf](../island_ffref/halogen/d5cp01216c2_suppl.pdf)
- [s/d5cp01216c.pdf](../island_ffref/s/d5cp01216c.pdf)
- [s/d5cp01216c1_suppl.zip](../island_ffref/s/d5cp01216c1_suppl.zip)
- [s/d5cp01216c2_suppl.pdf](../island_ffref/s/d5cp01216c2_suppl.pdf)

**适合归属：外部 OPLS-based 有机硫/卤素 API 特定模型。** 正文 p5 给 OPLS 起点/新扭转；44 页 SI 的 p33–36 表 S24–S28 给五种 API 的多套电荷，p37–40 给晶体/升华焓验证。ZIP 有 `opls.lib`、DLPGEN `.gen`、PCF 和优化结构；库头第 7–17 行明确 kJ/mol、sigma 用 Å、几何混合、1–4 LJ/Q 均 0.5。

**Si 与导入边界：** 无 Si/PDMS 或通用 P。它是很有价值的机器可读 OPLS 资产，适合从一个确切 API 起步。README 明说不同模型须用 SI 替换电荷，默认库不是所有模型的最终版本；碘 X-site 需额外站点、几何与力传播。必须绑定 `.gen` 原子名、库行、选定电荷列与授权；不需 X-site 的分子可优先审查。

### 13. Molecular Modeling of Energetic Materials: The Parameterization and Validation of Nitrate Esters in the COMPASS Force Field

S. W. Bunte; H. Sun；2000；Journal of Physical Chemistry B；DOI `10.1021/jp991786u`

本地文件：

- [halogen/jp991786u.pdf](../island_ffref/halogen/jp991786u.pdf)
- [pcff/jp991786u.pdf](../island_ffref/pcff/jp991786u.pdf)

**适合归属：外部 COMPASS Class-II 硝基/硝酸酯扩展。** 正文 p4 给 Class-II、LJ9–6 与 bond increment 形式；PDF p11–12 附录表 12 有 n3o/o12/o2n、equivalence、电荷增量、四阶键角、torsion、Wilson oop、LJ、BB、BA、angle-torsion 和 angle-angle-torsion。没有独立 SI，参数确实在正文附录。

**Si 与导入边界：** 无 Si，硝酸酯也不是磷酸酯。文件位于 halogen 不改变其化学归属。o2n 到 o2 等映射要保留作用类别，不能跨 PCFF 挪行；这是增量发布，需同版 COMPASS 基础和未列耦合的明确规则。较适合 Class-II 硝酸酯窄域。

### 14. LigParGen gkx312_supp.xlsx：跨引擎比较工作簿

2017；Nucleic Acids Research；DOI `10.1093/nar/gkx312`

本地文件：

- [opls/gkx312_supp.xlsx](../island_ffref/opls/gkx312_supp.xlsx)

**适合归属：OPLS/LigParGen 跨引擎能量基准；本地只有工作簿。** `Comparison` 工作表 A/B 列是分子/SMILES，C–V 列为 GROMACS、NAMD、OpenMM、BOSS 的总能、键、角、torsion、非键能，单位 kcal/mol。对应 DOI `10.1093/nar/gkx312`，正文未供。

**Si 与导入边界：** 这些是单点能，不是 LJ、torsion 或电荷系数，也不能推断 PDMS 覆盖。缺对应坐标、BOSS/OPLS版本和参数文件时，SMILES 不能重现同一个单点能。适合作验证参考，不作系数源。

### 15. Refined OPLS All-Atom Force Field Parameters for n-Pentadecane, Methyl Acetate, and Dimethyl Phosphate

Krzysztof Murzyn; Maciej Bratek; Marta Pasenkiewicz-Gierula；2013；Journal of Physical Chemistry B；DOI `10.1021/jp408162d`

本地文件：

- [opls/jp408162d_si_002.pdf](../island_ffref/opls/jp408162d_si_002.pdf)
- [p/jp408162d.pdf](../island_ffref/p/jp408162d.pdf)
- [p/jp408162d_si_001.zip](../island_ffref/p/jp408162d_si_001.zip)
- [p/jp408162d_si_002.pdf](../island_ffref/p/jp408162d_si_002.pdf)

**适合归属：三个分别版本化的外部 OPLS-compatible 修订域。** 正文 p5 表 1 给 n-十五烷非键、电荷和新 torsion；p7 表 3/4 给乙酸甲酯及 DMP−/DMPH 的非键/电荷；p8 表 5 给 P 相关 bonded 项。6 页 PDF SI 主要为模拟/熔点对照；ZIP 的 README 和 PDB/GRO 是构型/速度，没有隐藏参数库。

**Si 与导入边界：** 无 Si；这是很具体的有机磷 OPLS 来源。作者明确把 RESP 电荷与调整后的 OPLS 非键组成新模型，不是 GAFF/RESP 工作流产物。原始/新参数、DMP−/中性 DMPH 不可混用，也不能补 PCFF 的 P 类型。通用基础项、两个状态的最终图、特殊对和水模型要单独复核。先独立做 DMP− 或 DMPH 的优先级较高。

### 16. The COMPASS force field: parameterization and validation for phosphazenes

H. Sun; P. Ren; J. R. Fried；1998；Computational and Theoretical Polymer Science

本地文件：

- [p/main.pdf](../island_ffref/p/main.pdf)
- [pcff/compass_ff.pdf](../island_ffref/pcff/compass_ff.pdf)

**适合归属：外部 COMPASS Class-II 中性磷腈及明确覆盖的取代基。** `p/main.pdf` 和 `pcff/compass_ff.pdf` 是同篇不同字节扫描。PDF p16–18 附录确有 n2=/p4=、非键、bond increments、四阶键角、torsion、oop 和多种耦合；p4–5/p14 说明 P–N 电荷分离约束为 0.35e 的处理。无独立 SI；DOI 在本地未可靠核定。

**Si 与导入边界：** 不补 Si。不能只看 ESP/验证表而忽略附录参数；但最终 constrained increment 与原始 QM 电荷不同，不可直接使用原始 ESP。P–N、P–F/P–Cl 只适用于明确磷腈连接，不转给 PF5、磷酸酯、膦氧或带电 P。基础 c4/h1/o2 行、扫描标签/列、全部 cross-term 和 end-group 需同版独立重建，不能用宽泛 >P 模板取代。磷腈窄域价值高，审查成本也高。

### 17. Polysiloxanes: ab initio force field and structural, conformational and thermophysical properties

Huai Sun; David Rigby；1997；Spectrochimica Acta Part A

本地文件：

- [pcff/1-s2.0-S1386142597000139-main.pdf](../island_ffref/pcff/1-s2.0-S1386142597000139-main.pdf)

**适合归属：外部 CFF/PCFF-related Class-II 硅氧烷模型。** 正文 p3 给四阶键角、torsion、oop、cross-term 与 LJ9–6，并说明按需要纳入耦合；PDF p21–22 附录 A 有 si4/si4c/o2 定义、O–Si increment、Si–O 键、六类角、多类 torsion（含显式 wildcard）、LJ、BB、BA。没有另附 SI。卷页为 53 (1997) 1301–1323，DOI 在本地未可靠核定。

**Si 与导入边界：** 是直接相关的 Si–O–Si 来源，区分含 H 和全非 H 取代 Si 的非键类型；不能把这里的 si4 与 pinned FRC 的 si 自动等同。基础 Si–C/C–H 行、特殊对/混合、角色/对称性和正式省略规则仍需重建。“按需拟合”不授权对任意新图缺失耦合补零。若目标是保留 Class-II 的硅氧烷，这是高优先来源；不意味着 Si–N/Si–S/卤素全覆盖。

### 18. COMPASS: An ab Initio Force-Field Optimized for Condensed-Phase Applications—Overview with Details on Alkane and Benzene Compounds

H. Sun；1998；Journal of Physical Chemistry B

本地文件：

- [pcff/COMPASS, An ab Initio Force-Field Optimized forCondensed-Phase Applications.pdf](../island_ffref/pcff/COMPASS%2C%20An%20ab%20Initio%20Force-Field%20Optimized%20forCondensed-Phase%20Applications.pdf)

**适合归属：已发表 COMPASS Class-II 烷烃/芳烃的独立外部控制域。** 正文 p4 表 1 给类型/等价；PDF p24–25 附录 B 给 alkane/benzene 的增量、四阶键角、torsion、oop、LJ9–6、BB、BB1–3、BA、AA、end/middle-bond-torsion、angle-torsion、angle-angle-torsion。p23 附录 A 是训练/验证模型分子清单。无独立 SI。JPCB 102 (1998), 7338–7364；本地提取文字未直接核出 DOI，未联网补录。

**Si 与导入边界：** 附录 A 列出 silane/disiloxane 等，不等于附录 B 公布了这些分子的完整系数。烷烃/苯参数真实且丰富，适合 Class-II 数值控制；硅域另审。需逐域解析耦合/等价、明确零/空白与版本，保持既有 published COMPASS 比较源及 native PCFF 身份分离。

### 19. Force field for computation of conformational energies, structures, and vibrational frequencies of aromatic polyesters

H. Sun；1994；Journal of Computational Chemistry；DOI `10.1002/jcc.540150708`

本地文件：

- [pcff/J Comput Chem - July 1994 - Sun - Force field for computation of conformational energies  structures  and vibrational.pdf](../island_ffref/pcff/J%20Comput%20Chem%20-%20July%201994%20-%20Sun%20-%20Force%20field%20for%20computation%20of%20conformational%20energies%20%20structures%20%20and%20vibrational.pdf)

**适合归属：外部 CFF93 芳香聚酯 PHBA/PET 及苯甲酸酯模型。** PDF p5 表 I 给类型/LJ，p7 表 II 给增量，p9–10 表 IV 给 valence 和 cross-coupling；p8 明确不重发烷基/苯基基础参数。未供独立 SI。

**Si 与导入边界：** 无 Si。可以提取酯基 Class-II 行，须区分原始 HF/6-31G* 数据与作者缩放后的最终参数，并取得同代基础烃、全部耦合和非键约定。不能从 OPLS/GAFF 填缺，适合作聚酯外部扩展。

### 20. Computer simulations of poly(ethylene oxide): force field, PVT diagram and cyclization behaviour

David Rigby; Huai Sun; B. E. Eichinger；1997；Polymer International；DOI `10.1002/(SICI)1097-0126(199711)44:3<311::AID-PI880>3.0.CO;2-H`

本地文件：

- [pcff/Polymer International - 1999 - Rigby - Computer simulations of poly ethylene oxide   force field  pvt diagram and.pdf](../island_ffref/pcff/Polymer%20International%20-%201999%20-%20Rigby%20-%20Computer%20simulations%20of%20poly%20ethylene%20oxide%20%20%20force%20field%20%20pvt%20diagram%20and.pdf)

**适合归属：外部 CFF/COMPASS-related Class-II 醇、醚与聚醚。** 题首页是 Polymer International 44 **(1997)** 311–330；文件名中的 1999 不是出版年。PDF p17–20 附录有类型/等价、增量、四阶键角、torsion、oop、LJ9–6、BB、BB1–3、BA、AA 和 torsion cross terms。正文 p3–4 给形式/参数化；并非只有 PVT 性质。无独立 SI。

**Si 与导入边界：** 无 Si，不能把醚 O–C 转给 Si–O。o2/o2e/o2h、c4o/h1o 的实际/形式类型要保留角色。需同代基础烃、特殊对和完整图交互审核，同称 CFF 不意味着任意版本可混用。Class-II 聚醚扩展优先级较高。

### 21. An ab Initio CFF93 All-Atom Force Field for Polycarbonates

Huai Sun; Stephen J. Mumby; Jon R. Maple; Arnold T. Hagler；1994；Journal of the American Chemical Society；DOI `10.1021/ja00086a030`

本地文件：

- [pcff/ja00086a030.pdf](../island_ffref/pcff/ja00086a030.pdf)

**适合归属：外部 CFF93 碳酸酯/聚碳酸酯。** 正文 p2 给能量形式，p4 给参数化，PDF p9–10 表 10 给类型、valence、oop、cross terms、非键与电荷相关参数；p4 明确烷基/苯基基础另有文献。未供独立 SI。

**Si 与导入边界：** 无 Si/P 专用行。碳酸酯与磷酸酯的名称相近不构成可转移性；扫描原子标号、同代基础库、全部 Class-II 依赖和端基仍需审核。适合碳酸酯溶剂/聚合物窄域，不作硅/磷补参。

### 22. Development and Testing of the OPLS All-Atom Force Field on Conformational Energetics and Properties of Organic Liquids

William L. Jorgensen; David S. Maxwell; Julian Tirado-Rives；1996；Journal of the American Chemical Society；DOI `10.1021/ja9621760`

本地文件：

- [pcff/ja11225.pdf](../island_ffref/pcff/ja11225.pdf)

**适合归属：原始 OPLS-AA 数值/convention 的参考来源。** `ja11225.pdf` 实际是 **13 页 SI，不是正文**：p2 给形式及 1–4 因子 0.5，p3–8 表 1–5 给非键/电荷，p9 表 6 给键角，p10–13 表 7–9 给 torsion。p1 有个人使用、再分发须获许可的明确声明。

**Si 与导入边界：** 放在 pcff 文件夹不改变 OPLS 归属；这份 SI 不建立 PDMS 的 Si–O 库。显式零可据原表存证，缺项不能置零。本地可阅读不等于可无条件重发整表。适合核对现有 OPLS 与外部派生模型，图特异覆盖、版本和发布许可另判。

### 23. A Quantum Chemistry Based Force Field for Poly(dimethylsiloxane)

James S. Smith; Oleg Borodin; Grant D. Smith；2004；Journal of Physical Chemistry B；DOI `10.1021/jp047434r`

本地文件：

- [pcff/jp047434r.pdf](../island_ffref/pcff/jp047434r.pdf)
- [pcff/.ipynb_checkpoints/jp047434r-checkpoint.pdf](../island_ffref/pcff/.ipynb_checkpoints/jp047434r-checkpoint.pdf)

**适合归属：独立外部 PDMS Smith2004 模型。** 建议名 `island_external_pdms_smith2004_v1`，这里只是命名建议，未注册。正文 PDF p3（印刷 20342）表 3 是不同硅氧烷/端基电荷；p5（20344）表 4 有 9 个 Buckingham pair rows、3 类键、6 类角、4 类 torsion、单位与组合规则；p6 式 6 是 Si–O–Si–O torsion 的 `sin²(theta)/sin²(theta0)` 调制。p3 III.A 指定 1–2/1–3 排除、1–4 非键乘 0.5。无本地独立 SI。

**Si 与导入边界：** **确实有 PDMS 参数，且参数量可观。** Buckingham exp-6/Coulomb、四阶成键项与特殊 bend–torsion 不等同于标准 OPLS12–6 或 PCFF9–6。表 3 区分内部 Si、trimethyl 端 Si、SiH3 等电荷；表 4 没有 Si–H bonded 行，因此 `[SiH3]O[SiH3]` 不自动可执行。也不能把 Si–H 缺口推广为所有甲基封端 PDMS 的相同缺口。需确定甲基封端图域、组合规则、柔性/约束版及空白高阶格的正式含义；完整保留其 Hamiltonian，优先级很高。新 checkpoint 与正文哈希一致。

### 24. Force-Field Benchmark for Polydimethylsiloxane: Density, Heat Capacity, Isothermal Compressibility, Viscosity and Thermal Conductivity

Zhirui Xiang; Chao Gao; Teng Long; Lifeng Ding; Tianhang Zhou; Zhenghao Wu；2025；Journal of Physical Chemistry B；DOI `10.1021/acs.jpcb.4c08471`

本地文件：

- [pcff/jp4c08471.pdf](../island_ffref/pcff/jp4c08471.pdf)
- [pcff/jp4c08471_si_001.pdf](../island_ffref/pcff/jp4c08471_si_001.pdf)

**适合归属：分别审查五个外部 PDMS 模型，不能合并。** 14 页 SI 的 S1(p2)、S2(p3)、S3(p4–5)、S4(p6)、S5(p7) 分别为 UA、HUANG、COMPASS、OPLS、Dreiding；p8 起给 NPA/Qeq 和验证。J44 的 **136 个系数候选只来自这份 SI 的五张表**，不是全目录所有数值行的数量。

**Si 与导入边界：** 直接覆盖甲基 PDMS/硅氧烷。S4 是全原子 OPLS-compatible 候选；J45 已记录其末端类型/电荷、单位、1–4/混合、improper 适用性等缺口。S3 有 Si–O 键长与 BA/BB 平衡依赖不一致、improper 分类疑点，不能自行更正后冒充原值；HUANG 也须按实际混合 Hamiltonian 审查，不能因有 class2 项就假定其与完整 PCFF 等同。未声明 SiH3/OH/O–C 端基不能推广。参数候选和基准价值高，当前已审域仍 blocked。

### 25. New Force-Field for Organosilicon Molecules in the Liquid Phase

Miguel Jorge; Andrew W. Milne; Maria Cecilia Barrera; José R. B. Gomes；2021；ACS Physical Chemistry Au；DOI `10.1021/acsphyschemau.1c00014`

本地文件：

- [pcff/pg1c00014.pdf](../island_ffref/pcff/pg1c00014.pdf)
- [pcff/pg1c00014_si_001.pdf](../island_ffref/pcff/pg1c00014_si_001.pdf)

**适合归属：独立外部 PolCA 联合原子有机硅模型。** 正文 p5 式 9–15 给约束键、谐角、RB torsion、LJ12–6/Coulomb、Lorentz–Berthelot，并明确取消 1–4；p6 表 1–2 给 DDEC/类型电荷，p7 表 3–4 给键长/角，p8 表 5 给 torsion，p9 表 6 给 LJ。76 页 SI：S1–S70 给实验/极化依据，p47 S71 给 QM 溶剂，p48–65 给键角来源/扫描及 S72 torsion，p68–72 给验证。

**Si 与导入边界：** 覆盖烷基硅烷、烷氧基硅烷、硅醇、桥氧硅氧烷，包括六甲基二硅氧烷，但 CHx 是联合原子、键为约束，部分性质含极化后处理。作者有明示等同假设和 Si 电荷中性约束算法，不能在当前禁止重分配的路径中静默套用，也不能拆成全原子 OPLS/PCFF。需同版 PolCA 烷基基础、RB/约束与单位重建；表 3 的 OH–H 数值/单位还需来源核验，不能自行修正疑似小数点。适合独立 UA 小分子分支，不自动支持长 PDMS/Si–H；优先级高。

### 26. Optimizing cation–π force fields for molecular dynamics studies of competitive solvation in conjugated organosulfur polymers for lithium–sulfur batteries

Diptesh Gayen; Yannik Schütze; Sébastien Groh; Joachim Dzubiella；2025；Physical Chemistry Chemical Physics；DOI `10.1039/D4CP04484C`

本地文件：

- [s/d4cp04484c.pdf](../island_ffref/s/d4cp04484c.pdf)
- [s/d4cp04484c1_suppl.pdf](../island_ffref/s/d4cp04484c1_suppl.pdf)

**适合归属：独立外部 OPLS-AA/corr. 的 Li–π 交互修正。** 正文 p3 式 1–2 为 LJ12–6/Coulomb **加 `−C4/r4`**；p7 表 4 给 C4 和特定 pair sigma。17 页 SI 的 S1–S5 给拟合设置/修正参数，图 S1 给 CA/CW/CS/S/SH 类型，后有能量面与 PMF 对比。

**Si 与导入边界：** 无 Si，适合含硫共轭电极/电解液。不能只拷 sigma 而丢掉 ion-induced-dipole 项；需保留 pair override 和盐/溶剂基础，重建 PTBT/端基电荷、C4 白名单及引擎力导数/单位。是专项模型，不是普通 OPLS 或 GAFF 的通用更新。

### 27. Large Scale Molecular Model Construction of Xishan Bituminous Coal

Zhiqiang Zhang; Qiannan Kang; Shuai Wei; Tao Yun; Guochao Yan; Kefeng Yan；2017；Energy & Fuels；DOI `10.1021/acs.energyfuels.6b02623`

本地文件：

- [s/ef6b02623.pdf](../island_ffref/s/ef6b02623.pdf)

**适合归属：含 S/N/O 的结构/建模测试参考。** 正文表 1–4 是元素、NMR、XPS 和结构数据；p6 说明 SI 为三维煤模型 PDB，但本地未供。没有完整 force-constant 库。

**Si 与导入边界：** 研究去矿物质有机煤 C/H/O/N/S，不是硅氧烷。XPS binding energy、原子含量、模拟密度都不是键/角力常数；采用 MM/MD 不证明底层参数来源可重建。不建议作为 PCFF/OPLS 新增数值源。

### 28. Improved Description of Sulfur Charge Anisotropy in OPLS Force Fields: Model Development and Parameterization

Xin Cindy Yan; Michael J. Robertson; Julian Tirado-Rives; William L. Jorgensen；2017；Journal of Physical Chemistry B；DOI `10.1021/acs.jpcb.7b04233`

本地文件：

- [s/jp7b04233.pdf](../island_ffref/s/jp7b04233.pdf)
- [s/jp7b04233_si_001.pdf](../island_ffref/s/jp7b04233_si_001.pdf)

**适合归属：分别版本化的外部 OPLS-AAx 与 OPLS/CM5x 硫模型。** 正文 p3–4/表 1 给新电荷/LJ和新 torsion 讨论；19 页 SI 从 p2 讲拟合，p9 表 S2 给 XS–S 距离、几何、力常数及涉及 XS 的显式零 torsion，其余表/图给性质和 torsion 曲线验证。

**Si 与导入边界：** 无 Si，适合硫醇/硫醚/二硫键等明确环境。离原子站点不能折回 S 后仍称保留各向异性；额外粒子或虚拟站点的几何、约束、力传播都属于模型定义。OPLS-AAx 与 CM5x 电荷须分开。须重建数值 torsion 和基础库，不从曲线估值；SI 的角力常数 degree/radian 约定须特别核验。硫扩展价值高，但需额外站点后端支持。

### 29. Ab Initio Calculations and Force Field Development for Computer Simulation of Polysilanes

H. Sun；1995；Macromolecules 28(3), 701–712；DOI `10.1021/ma00107a006`

本地文件：

- [pcff/ma00107a006.pdf](../island_ffref/pcff/ma00107a006.pdf)
- [pcff/.ipynb_checkpoints/ma00107a006-checkpoint.pdf](../island_ffref/pcff/.ipynb_checkpoints/ma00107a006-checkpoint.pdf)（截断快照，见下文）

**适合归属：高优先级的外部 PCFF/CFF93 Class-II 聚硅烷与烷基/苯基硅烷扩展。** 已重新读取用户补全的本地 PDF：12 页，4,132,865 字节，SHA-256 为 `9c6d6315b91e0c15a3cba5f90f011a862da79ae193bfe66900bfbdc9c5c8fe48`，存在完整 EOF。正文 PDF p2（印刷 p702）明确采用 CFF93-type all-atom 力场，p6（印刷 p706）将最终修订后的模型称为 **PCFF**，p11（印刷 p711）直接给出参数附录。因此这不是只有方法描述的论文。原先“副本不可解析”的结论已由本次完整阅读替换；原 J44/J45 历史证据未变。

**先辨认 Si 化学域：本篇的 PDMS 是 poly(dimethylsilane)，重复结构为 `–Si(CH3)2–`，主链是 Si–Si；不是 poly(dimethylsiloxane) 的 `–O–Si(CH3)2–`。** PDF p1 和 p8 均明确这一点。p3 表 2 的 14 个模型覆盖硅烷、二/三/四硅烷、甲基/乙基/二甲基硅烷，以及苯基取代衍生物；不包含氧桥。它能推进 Si–H、Si–C、Si–Si 和 Si–芳香碳窄域，没有给出 Si–O、Si–O–Si 或磷参数，不能解决当前硅氧烷 PDMS 的来源缺口。

| 证据位置（一基 PDF 页码） | 实际内容 | 对 ISLAND 的意义 |
| --- | --- | --- |
| p3，§III 势函数与电荷定义 | 键/角二至四阶展开、三项 Fourier torsion、谐和 out-of-plane；BB、BA、AA、BT、AT、AAT；LJ9–6 与 Coulomb | 保留 Class-II 全形式；不能删耦合后称原模型 |
| p4，表 5 | Si–H、Si–C、Si–Cp、C–H、Cp–H 五类 bond increment；Si/H 的初值与最终 vdW 参数 | 是可追溯的电荷及非键数值来源；应区分初值/最终值 |
| p11，附录 intrinsic terms | 4 条键、10 条角、7 条 torsion、1 条芳香 out-of-plane 记录 | Si–Si、Si–H、Si–C、Si–Cp 有显式行；`Cp` 指 phenyl carbon |
| p11，附录 cross-coupling | 8 条 BB、9 条 BA、15 条 AA、7 条 torsion–bond、3 条 angle–torsion、6 条 angle–torsion–angle 记录 | 有实际耦合数值；记录数不等于任意最终图的完整交互覆盖 |
| p5–6，表 9–12；p7 表 13；p9–10 表 17 | 作者报告的液态内聚能/密度、二硅烷频率/几何/转动势垒、构象能和晶体比较 | 可作为后续独立重建基准；不代表 ISLAND 已复现 |

例如，p11 的 Si–C 键给 `b0=1.899 Å`，`k2=189.65`、`k3=−279.42`、`k4=307.51`，对应单位分别为 kcal·mol⁻¹·Å⁻²/⁻³/⁻⁴；H–Si–C 角给 `θ0=112.10°`、`k2=36.48 kcal·mol⁻¹·rad⁻²`、`k3=−12.81 kcal·mol⁻¹·rad⁻³`。表 5 的 Si–H/Si–C 增量分别为 `0.126/0.135 e`；最终 Si 的 `ε=0.190 kcal/mol`、`D=4.450 Å`，H 的 `ε=0.023 kcal/mol`、`D=2.940 Å`。这些只是说明数值行存在的例子，未形成可执行参数包。附录表格已用 PDF 内嵌页图核对，不能只依赖 OCR 列对齐。

**电荷、非键与依赖边界：** p3 定义电荷为有方向的键增量之和，并明确同类原子间增量为零；p4 说明表 4 的 ESP 拟合结果经平均及作者的 0.9 缩放得到表 5 的最终增量。不能把表 4 电荷当最终图电荷，也不能对表 5 再乘 0.9 或在运行时归一化。p3 的 LJ 形式为 `ε[2(r*/r)^9−3(r*/r)^6]`，不是 OPLS 的 LJ12–6；D/r* 与后端距离参数的映射须独立校验。正文规定至少相隔三个共价键的分子内原子对参与非键作用，包含 1–4；但本篇未单独列出异类混合公式和显式 1–4 缩放表，不能套用 OPLS 的 0.5。p5 明确部分固定参数来自烷烃/苯基础力场；表 5 未独立列出 C/Cp 的 vdW，附录也未重列完整 C–H/C–C/芳环基础成键库，需沿作者明示的同版来源核定，不能跨家族补齐。

**尚不能直接导入的原因：** 附录多条角的 k4 栏为空白，没有在本篇给出统一的“空白即零”规则；这与 torsion 等行明写的 `0.00` 不同。未列出的端键/其他 Class-II 子项也须确定是原模型不适用、定义冗余或确实缺失，不能由转换器初始化补零。应按附录脚注保留中心原子、角顺序和 BT 的中心键角色；通配 torsion 与精确行并存，须验证匹配优先级和 Cp/H 类型映射。作者在正文声明的参数转移，只授权沿其指定模型依赖核验，不意味着可将本篇某行直接拼入现有 pinned PCFF。

**SI、许可与建议：** 当前仅供这份正文，参数附录内嵌于 p11，未有独立配套 SI/机器参数文件；首页有 ©1995 ACS 版权标记及机构下载水印，未据此建立参数或原文的公开再分发许可。建议将甲基硅烷/二硅烷列为后续 external PCFF/CFF93 聚硅烷窄域候选，与组 17 的硅氧烷模型分开。可用工作名 `island_external_sun1995_polysilane_class2_v1`，**本次未创建或注册该 profile**。不能直接加入 OPLS-AA；改成 OPLS 需要另行拟合与验证，而不是单位换算。它也不能证明 native PCFF、COMPASS 或 msi2lmp 的全面覆盖。

**附带快照说明：** 新增的 `pcff/.ipynb_checkpoints/ma00107a006-checkpoint.pdf` 为 2,097,152 字节，无 EOF，pypdf 报 `PdfStreamError: Stream has ended unexpectedly`；其全部字节与完整正文的同长度前缀一致。它是不可独立解析的截断快照，不是额外论文、SI 或可执行参数库，未用于补参。完整正文的阅读结论不受影响。逐页定位、简短原文及新旧文件身份见 [本次补读证据](evidence/literature_review_2026_10_27/ma00107a006-reading.json)。

## 如何理解“适合加入”

应分别记录：数值行存在、化学分型明确、电荷完整、所需交互完整、后端形式支持、独立能量/力通过、工作流验证和科学验证。本文只推进来源辨认与归属；没有宣布任何新 profile 已通过后面的执行/科学门槛。

建议的实施次序取决于目标：若坚持全原子 PDMS，先完整重建 Smith2004 或继续澄清 Xiang S4，二者分开；若坚持 Class-II 硅氧烷，先 Sun–Rigby 及其同版依赖；若目标是 Si–H/Si–C/Si–Si，则优先核验新增的 Sun1995 聚硅烷附录及基础依赖，不能将其中的 PDMS 误作硅氧烷；若先做易审计的 OPLS 窄域，则优先 Liang 的 PEO/LiTFSI、Lopes 的单一 API 或 Murzyn 的单一磷酸酯状态。PolCA 适合另开 UA 支持，不能作为全原子接口的隐式替代。

每次只取一个作者明确规定的模型组合。基础库借用只有作者明确声明且版本可绑定时才可形成依赖链；不能由开发者跨家族填缺。作者给出的真实零、约束、拟合/缩放算法须保留出处；缺行、未标单位、未解释空白和图特异电荷歧义继续阻断。文献模型使用的电荷拟合/中性约束也不自动授权在 ISLAND 当前路径进行电荷重分配。

本地可读取、论文开放获取和参数再分发权限是不同事实。原始 OPLS SI 的明确限制、API ZIP 的来源/电荷模型、无 LICENSE 的工程归档均需单独核定；对未知权限不能直接断言“未授权”或“可公开重发”。本报告不复制整套原文系数表，也不重新发布本地 PDF/ZIP。

原有 Ewald/PME、后端压力、cutoff recovery 与 long-NPT 限制保持不变；这份阅读评估不构成实验一致性、平衡或生产可用性声明。
