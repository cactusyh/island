# PDMS 有依据的参数补充与剩余缺口

本次已将可核实的 PDMS 参数整理为独立数据记录
`island_external_pdms_pcff_source_supplement_v2`，并补充了此前 8 条角势空白的完整 FRC 源行证据。
**当前资料仍不足以组成完整、可执行的 PDMS Class-II 力场。** 未创建执行 profile，未执行 PDMS MD，也未修改原生 PCFF、OPLS 或此前 Sun 文献参数库。

## 依据和版本选择

本次使用两个互相隔离的证据层：

1. 此前已提取的 Sun 1995 / Sun–Rigby 1997 共 105 条文献记录，仅用于核对。
2. 本地 LAMMPS 分发的完整 `pcff.frc`，用于本次最终图的源行、类型、质量和电荷增量审计；所有指派候选来自这一文件，不从不同论文或 COMPASS、OPLS、GAFF 拼接系数。

1995 年文献为 *Ab Initio Calculations and Force Field Development for Computer Simulation of Polysilanes*，DOI `10.1021/ma00107a006`，其中 PDMS 指 poly(dimethylsilane)，主链为 Si–Si。
1997 年文献为 *Polysiloxanes: ab initio force field and structural, conformational and thermophysical properties*，Appendix 在 PDF 第 21–22 页；本地 DOI 未核实，仍保留为空。

本地分发文件的 SHA-256：

```text
pcff.frc
e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c
pcff.rlb
82087ba7e9117d3689ee93870b35cd21ddbc2db6b16142cdcf7696b5af7dcc36
pcff_templates.dat
dbde1d51aae2ce3a38d6bcb3e70c2c4ca2116a69904ea8c1b8ede58f06a60093
```

论文与本地 FRC 并非同一数值版本，例如：

| 参数 | 1997 论文 | 本地 FRC |
| --- | --- | --- |
| Si–O：r0, K2, K3, K4 | 1.640, 350.1232, −517.3424, 673.7067 | 1.6562, 306.1232, −517.3424, 673.7067（行 1775） |
| Si–O–Si：theta0, K2, K3, K4 | 159.0, 8.5000, −13.4188, −4.1785 | 157.0260, 9.0740, −19.5576, 8.5000（行 2365） |
| 甲基取代 Si 的 LJ9–6：rmin, epsilon | 4.290, 0.131 | 4.2840, 0.07000（行 3357） |
| O 的 LJ9–6：rmin, epsilon | 3.300, 0.080 | 3.3500, 0.24000（行 3347） |

本次没有把论文中的 Si–O 项覆盖进 FRC，也没有把 FRC 的 C/H 项宣称为 1997 论文的同版基础参数。
源记录按原始单位保存：kcal/mol、Å、平衡角 degree、角位移 radian、质量 dalton、电荷 e；数值转换为恒等转换，保留全部原始十进制列。

## 8 条角势的补充证据

这些原论文行的 `k4=null` 保持原样。下表是本地 FRC 明确打印的**整行**，另存于新证据目录；不能把论文舍入后的前三项与 FRC 的第四项拼成新系数。
前 3 列与论文显示精度一致，不等于证明论文未印出的项必然为零。
8 行均为版本 2.2、reference 9；FRC reference 9 标注 Huai Sun、1-November-94、added parameters for silanes。

| 原论文标签 | FRC 行 / 标签 | theta0 | K2 | K3 | K4 |
| --- | --- | ---: | ---: | ---: | ---: |
| Si–Si–H | 2410 / h si si | 112.0893 | 22.5062 | −11.5926 | 0.0000 |
| H–Si–H | 2409 / h si h | 108.6051 | 32.5415 | −8.3164 | 0.0000 |
| Si–Si–Si | 2411 / si si si | 114.2676 | 24.9501 | −19.5949 | 0.0000 |
| H–Si–C | 2406 / c si h | 112.0977 | 36.4832 | −12.8094 | 0.0000 |
| Si–C–H | 2196 / h c si | 112.0355 | 28.7721 | −13.9523 | 0.0000 |
| Si–Si–C | 2407 / c si si | 113.0000 | 19.4692 | −34.3471 | 0.0000 |
| Cp–Cp–Si | 2282 / cp cp si | 120.0000 | 30.4689 | −23.5439 | 0.0000 |
| Si–C–C | 2158 / c c si | 112.6700 | 39.5160 | −7.4430 | 0.0000 |

**PDMS 必须使用硅氧烷类型 `sio`，不能直接使用以上硅烷 `si` 行。** 例如实际 PDMS 的 Si–C–H 对应行 2197：`h c sio 111.5360 30.2481 -15.5255 0.0000`。
这与行 2196 的 `h c si` 有明确数值差别。

## 实际补入的数据

完整源行位于 [source-rows.json](evidence/pdms_source_completion_v2/source-rows.json)，每行记录文件哈希、行号、原文、源版本及 reference、类型、字段、逐列单位、原值/转换值、明确印出的零及证据校验值。
其中 39 行用于有效 PDMS/二硅氧烷图的审计，8 行仅用于上表的硅烷文献核对；不是 47 个新拟合参数。

| 类别 | 独立完整源行数（包括硅烷核对） |
| --- | ---: |
| 四次键势 | 4 |
| 四次角势 | 16 |
| 三项扭转 | 5 |
| LJ9–6 | 5 |
| bond–bond | 8 |
| bond–angle | 8 |
| angle–angle | 1 |

另在每个图绑定记录中保留了类型/质量源行、普通等价映射、电荷增量源行、组分/分子归属、原子顺序、全部相互作用和平衡依赖。审计按精确标签、再按源定义的普通等价关系选择，不把 `sio` 变成 `si`，不以通配符或低阶势替代缺少的 Class-II 项。

所有有效图的本地源电荷增量均已覆盖；独立 Decimal 重建证实总电荷严格为零，没有归一化或电荷再分配。例如 HMDS：Si `+0.520 e`、O `−0.440 e`、甲基 C `−0.259 e`、H `+0.053 e`；DP4 中内部 Si 为 `+0.640 e`，端 Si 为 `+0.520 e`。
这些电荷只属于本地 FRC 诊断记录，不能标成 1997 论文电荷或 msi2lmp 自动生成的电荷。

电荷重建的零原子基准与缺少系数时填零是不同操作：本次每条键都有显式源增量行，源零增量亦有原文证据；任一增量缺行即阻断。

## 仍不能补造的耦合项

最终图来自 `[SiH3]O[SiH3]`、`C[Si](C)(C)O[Si](C)(C)C` 和 PSMILES `[*:1]O[Si](C)(C)[*:2]`。
聚合物使用 head `[*][Si](C)(C)C`、tail `[*]C` 并显式补氢。DP 计重复单元；头端基额外含一个 Si。
DP1 与 HMDS 同构，但稳定 ID 和原子顺序不同，分别绑定。

| 最终图 | 原子数 | 键数 | 类型/电荷覆盖 | 缺失的具体耦合请求 |
| --- | ---: | ---: | --- | ---: |
| disiloxane | 9 | 8 | 完整 | 54 |
| HMDS | 27 | 26 | 完整 | 378 |
| PDMS DP1 | 27 | 26 | 完整 | 378 |
| PDMS DP2 | 37 | 36 | 完整 | 528 |
| PDMS DP4 | 57 | 56 | 完整 | 828 |

以 HMDS 为例，49 个 angle、49 个 bond–bond、49 个 bond–angle 均有完整源行，60 个 proper torsion 主项也有源行。
但 angle–angle 缺 78 个请求；bond–bond 1–3、end-bond–torsion、middle-bond–torsion、angle–torsion、angle–angle–torsion 各缺 60 个请求，总计 378。
因此阻碍完整性的并非只有之前提到的 8 个角势空白。

5 个图合计 2166 个缺失的站点请求，归为 47 个**保留方向和中心角色的源查询键**，不是 2166 个独立待拟合常数。
见 [missing-source-keys.json](evidence/pdms_source_completion_v2/missing-source-keys.json)。`study.tar.xz` 中的 `missing-interactions.json` 保存每个站点、标签、普通等价标签和平衡依赖。
其中 AA 的中心和共享臂不能随意反转；不同方向的查询未被简单合并成一个可转用系数。

例如缺少 `end_bond-torsion_3: hsi sio osi sio` 与 `hc c3 sio c3` 等所对应的源行。
这些耦合系数不能从键长、平衡角或单独扭转势推导，也不能因为某个转换器把内存初始化为零就认定物理系数为零。
`CO[Si](C)(C)C` 保持 `typing_blocked`，不把 Si–O–C 自动视为 Si–O–Si。

## 独立重建与执行边界

- 独立 ASCII 解析器重建全部 47 个完整数值行、明确零列、相互作用数、类型/质量和 Decimal 电荷，不调用生产解析器或解析决策器。
- 当原值转为角度浮点比较时，明确统一 `degree * (pi/180)` 的运算顺序，没有放宽容差；源 Decimal 保持完全相同。
- 本地重新编译 `msi2lmp`，对同一最终图内部生成 CAR/MDF，仅作转换器对照。5 个图均因缺少 end-bond–torsion 退出 17，没有使用 `-ignore`。
- 转换器可执行文件 SHA-256：`59f3718d4e537d584011315b76d0bf80ca105b43468f87b384cfa4ce19d76123`。编译器、完整源码哈希和 revision 见 `converter-build.json`；这不代表对其他构建的二进制哈希作相同承诺。
- 记录不可变；加载需独立 SHA-256，且对当前最终图和钉住的源文件重新计算，拒绝重新签名后的系数、单位、电荷、源、图、ID/归属、依赖及执行权限篡改。
- `require_executable()` 抛出带完整结果与逐项缺口的 `PDMSCompletionBlockedError(status="parameter_incomplete")`。数据记录不是 C5 执行 bundle，也不会绕过现有 adapter。

尚未运行本次 PDMS 的能量/力、最小化或 NVT/NPT，故没有数值/工作流验收结论。
原生有限硅烷及历史 J17–J45 身份不受影响；原有 105 条 Sun 记录也未改写。

## 下一步如何真正完成

优先获取同一版本、可授权使用的缺失耦合参数表及对应函数/排除/1–4 约定，逐行补上上述源键；新的来源需建立独立模型身份，不能无说明地覆盖本地源。
如果无此资料，只能针对明确定义的模型做新的量化计算与联合参数拟合：显式构象/扭转扫描、Hessian 与多项 Class-II 耦合的联合拟合及独立验证，而不能声称这些新值来自论文。
还需解决本地 FRC 与参考后端的扭转/angle–angle–torsion 公式、单位和完整适用性，再通过能量分项、力、有限差分、最小化及续跑验证。
目前没有足够证据给这些缺项提供唯一可靠数值。

## 复现

在仓库根目录、支持 RDKit 的 Python 环境执行：

```bash
export PYTHONPATH=.:src
python scripts/complete_pdms_source_parameters.py \
  --distribution ../lammps/lammps/tools/msi2lmp/frc_files \
  --catalog docs/evidence/sun_siloxane_parameters/parameters.json \
  --references island_ffref --output /tmp/pdms-source-v2
python scripts/audit_pdms_source_completion.py \
  --study /tmp/pdms-source-v2 \
  --distribution ../lammps/lammps/tools/msi2lmp/frc_files \
  --output /tmp/pdms-source-v2-independent.json
python scripts/complete_pdms_source_parameters.py \
  --distribution ../lammps/lammps/tools/msi2lmp/frc_files \
  --catalog docs/evidence/sun_siloxane_parameters/parameters.json \
  --references island_ffref --output /tmp/pdms-source-v2 \
  --reconstruct --require-executable
# 上一命令必须退出 1：有效的 parameter_incomplete 结果，不是异常通过。
```

本地转换器对照（`--output` 应为尚未存在的目录）：

```bash
cc -O2 -o /tmp/pdms-msi2lmp \
  ../lammps/lammps/tools/msi2lmp/src/*.c -lm
python scripts/compare_pdms_source_converter.py \
  --study /tmp/pdms-source-v2 \
  --distribution ../lammps/lammps/tools/msi2lmp/frc_files \
  --executable /tmp/pdms-msi2lmp \
  --executable-sha256 59f3718d4e537d584011315b76d0bf80ca105b43468f87b384cfa4ce19d76123 \
  --output /tmp/pdms-converter-v2
```

若编译器/构建环境不同，需独立记录并核对新二进制哈希，再传入其真实值。不得直接信任待验证 bundle 内自报的哈希。

```bash
python -m pytest -q tests/test_pdms_source_completion.py
python scripts/validate_pdms_source_completion.py \
  --directory docs/evidence/pdms_source_completion_v2 \
  --distribution ../lammps/lammps/tools/msi2lmp/frc_files \
  --references island_ffref
```

完整验证计数、运行命令、环境和哈希记录于新目录的 `receipt.json` / `manifest.json`；测试日志与逐图重建材料分别存于 `validation.tar.xz`、`study.tar.xz`。

本次新增测试 27 项全部通过；支持的完整套件 3618 项通过、10 项跳过，J22–J45 回归包含于其中。Ruff、compileall、5 个环境的 pip check、差异空白检查、独立源行/PDF 重建、迁移和新进程重建，以及 J43/J44/J45/原 Sun 收据验证均通过。3279 个已有跟踪文件、16 个上一轮 Sun 补充文件和原 J44 的 52 个本地源文件字节哈希均保持不变。

保持 `production_validated=False`、`simulation_readiness="not_established"`，聚合物 Ewald/PME 未验证、后端压力不可用、cutoff recovery 和 long NPT 仍阻断；无平衡态或实验符合性结论。
