# `pcff/main.pdf` 对 PDMS 参数缺口的核对

已阅读 PDF 全部 7 页、4 张表及 6 幅图。**这篇论文不能补齐目前 PDMS 的 Class-II 力场参数；本次可用于填补缺口的新增系数为 0 条。** 它提供的是已有力场的应用、比较和构象计算结果，没有公开所用完整参数数据库，也没有列出目前缺失的耦合系数。

本次仅增加阅读报告与证据，不改参数、不创建执行 profile、不运行 MD。

## 论文身份

- 文件：`island_ffref/pcff/main.pdf`，263138 bytes，7 页。
- 标题：*Evaluations of forcefields for aromatic polysiloxanes, and some applications to poly(diphenylsiloxane)*。
- 作者：R. D. Patil、J. E. Mark。
- 期刊：*Computational and Theoretical Polymer Science* **10** (2000), 189–195。
- 文件印出的 PII：`S1089-3156(99)00079-3`。未从该文件核实 DOI，不将 PII 推测转换为 DOI。
- SHA-256：`c518da36d84d4826bd62d7fd4a9ecf1b784dfd8a0667b507a196637861a663d3`。
- 版权标注：© 2000 Elsevier Science Ltd., All rights reserved。此处不据此推定参数数据库的使用或再分发授权。

论文研究的主要对象是 **PDPS，聚二苯基硅氧烷**，重复单元为 `[Si(C6H5)2–O–]`。它与目标 PDMS `[Si(CH3)2–O–]` 同属硅氧烷，但苯基和甲基取代环境不同；其晶体比较对象是六苯基环三硅氧烷 HPT，也不是六甲基二硅氧烷 HMDS。

## 全文实际提供的数值

下表 PDF 页码从 1 开始，括号内为印刷页码。

| 位置 | 内容 | 能否作为当前缺失系数 |
| --- | --- | --- |
| 第 2 页（190），表 1 | HPT 晶格、空间群、体积、实验密度 | 不能；是特定芳香族环状分子的晶体参考数据 |
| 第 3 页（191），表 2 | O–Si–O / Si–O–Si 两种连续键对的两个转动角及 RIS 统计权重 | 不能；不是三原子角势系数，也不是 Fourier 扭转系数或 Class-II 耦合系数 |
| 第 4 页（192），表 3 | 两套力场的模拟能量及密度 | 不能；是模拟输出，不是单个相互作用参数 |
| 第 5 页（193），表 4 | PDPS 在 300 K 下的链尺寸、特征比及回转半径相关量 | 不能；是构象统计结果 |
| 第 7 页（195），讨论此前 PDMS 的 RMMC 工作 | Si–O–Si 与 O–Si–O 的角度约 147°、108° | 不能；没有配套力常数、原子类型、版本及完整势函数定义，不能作为本地 FRC 的平衡角替换值 |

图 1 是 PDPS 的结构/转动角示意，图 2–3 是密度与能量随时间的结果，图 4–6 是链尺寸及温度依赖。已核对 PDF 内的图像资源，没有隐藏在图中的数值参数表；PDF 也没有嵌入附件或参数附录。

表 2 的“Si–O–Si bond pair”描述围绕两条连续键的两个转动自由度，不应因标签相同而解释为 `Si–O–Si` 三原子角势。即使在已知温度下由统计权重求相对能量，也不能唯一分解出键、角、扭转、非键和多种耦合系数。

## 与当前缺口逐项比较

当前基线为 [PDMS 源参数补充 v2](pdms_source_completion_v2.md)，详细键见 [missing-source-keys.json](evidence/pdms_source_completion_v2/missing-source-keys.json)。本次只读取这些记录，没有重新定义其模型。

| 缺失项 | 当前缺失的定向源查询键数 | 本论文是否给出可用数值行 |
| --- | ---: | --- |
| angle–angle | 12 | 无 |
| bond–bond 1–3 | 7 | 无 |
| end-bond–torsion | 7 | 无 |
| middle-bond–torsion | 7 | 无 |
| angle–torsion | 7 | 无 |
| angle–angle–torsion | 7 | 无 |

因此 47 个保留中心/方向角色的查询键、5 个图合计 2166 个缺失站点请求均未减少。这些数字不是独立待拟合常数个数。本论文也没有补上当前模型需要核实的扭转相位、angle–angle–torsion 公式/单位、同版混合规则、1–4 与排除约定。

## 为什么不能把 PCFF2 参数拼接到本地 PCFF

第 2 页（190）§1.3 将其使用的 `pcff2` 称为 `compass 1.0`，并明确说明：

> Both pcff and pcff2 use the same equations to calculate the potential energy of a system, but the parameters used in the equation are different [30].

这句话是模型版本边界的证据，不能解释为两者数值可互换。论文没有给出其实际数据库的文件哈希、完整类型表、参数行或可复现的版本快照，不能据此识别为当前钉住的 `pcff.frc`。
此前审计也已确认 Sun–Rigby 1997 论文与本地 FRC 的 Si–O、Si–O–Si、LJ 系数存在数值差异。共同的力场家族和函数形式不足以授权跨版本拼接。

## 这篇论文仍有哪些用途

1. **模型选择依据。** 它说明已有 PCFF/PCFF2 实现曾用于含芳香取代基的硅氧烷；不能据此证明当前本地分发完整。
2. **有限范围的参考结果。** 表 1、3 的 HPT 实验密度为 1.220 g/cm³，两模型平均模拟密度分别为 1.184、1.217 g/cm³。这些可作为以后重建同一芳香族体系的比较资料，不能作为当前 PDMS profile 已验证的结论。
3. **追溯真实参数来源。** 第 7 页参考文献 [29] 指向已读的 Sun 1995；[31] 指向已读的 Sun–Rigby 1997；[30] 指向 MSI 1996 年 *Discover, Forcefield simulation user guides*；[32]、[33] 分别列出 Rigby–Sun–Eichinger 1997 和 Sun–Ren–Fried 1998。引用本身没有提供这些来源的缺失参数或使用授权。

要补齐当前模型，仍需获得**同一明确版本、允许使用且逐行可核实**的上述六类耦合表和配套约定。若改用另一套完整文献模型，应建立独立外部身份并重新验证，不能覆盖原 PCFF 或把缺项填零。

本次结论为 `parameter_incomplete`，未新增可执行参数或科学验收结论。保持 `production_validated=False`、`simulation_readiness="not_established"`。

## 证据与复查

来源、逐页/逐表位置、短原文证据和缺口基线校验值见 [review.json](evidence/pcff_main_pdms_review/review.json)，报告及证据的校验值见 [manifest.json](evidence/pcff_main_pdms_review/manifest.json)。

可使用 `pypdf` 读取全部页面复查；该过程只依赖本地 PDF，不需要网络：

```bash
../miniforge3/bin/python - <<'PY'
import hashlib
from pathlib import Path
from pypdf import PdfReader
p = Path('island_ffref/pcff/main.pdf')
assert hashlib.sha256(p.read_bytes()).hexdigest() == 'c518da36d84d4826bd62d7fd4a9ecf1b784dfd8a0667b507a196637861a663d3'
for n, page in enumerate(PdfReader(p).pages, 1):
    print(f'\n=== PDF page {n} ===\n{page.extract_text()}')
PY
```
