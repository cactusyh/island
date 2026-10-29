# Sun 文献硅氧烷参数补充

已将两篇本地论文中的 **105 条参数记录**加入独立、版本化的文献参数库：
`island_external_sun_siloxane_class2_data_v1`。这是可供逐行查询和最终图缺项审计的数值补充，尚非可执行力场配置。原有 PCFF、OPLS-AA 和 J17–J45 证据均未修改。

## 来源与允许的联合范围

| 来源 | 本地文件 | 提取位置 | 记录数 |
| --- | --- | --- | ---: |
| H. Sun, *Ab Initio Calculations and Force Field Development for Computer Simulation of Polysilanes*, Macromolecules 28 (1995), 701–712; DOI `10.1021/ma00107a006` | `island_ffref/pcff/ma00107a006.pdf` | PDF 第 4 页 Table 5、第 11 页 Appendix | 77 |
| H. Sun, D. Rigby, *Polysiloxanes: ab initio force field and structural, conformational and thermophysical properties*, Spectrochimica Acta Part A 53 (1997), 1301–1323 | `island_ffref/pcff/1-s2.0-S1386142597000139-main.pdf` | PDF 第 21–22 页 Appendix | 28 |

1997 年论文的 DOI 未由本地证据核实，保留为空。两篇 PDF 的 SHA-256 分别为：

```text
1995  9c6d6315b91e0c15a3cba5f90f011a862da79ae193bfe66900bfbdc9c5c8fe48
1997  65bda3a5b14cfd05bc778f2306aa6b7221dda8124d66b4d7c7fd817853d01fe7
```

1995 年论文中的 PDMS 指 **poly(dimethylsilane)，Si–Si 主链**；不能直接当成聚二甲基硅氧烷的完整模型。1997 年论文 §3.2 明确说明，连接甲基/苯基的 Si–C 键及相关角、扭转参数沿用先前硅烷模型，并在参考文献 [28c] 引用 1995 年论文。这允许有依据地逐项核对继承关系，但不授权任意合并两套参数。

1997 年 §3.3 重新拟合了 Si/O 的 LJ9–6 参数，因此 1995 年的 Si 非键参数仅作为比较记录。1997 年 C/H 基础参数引用的 [38] 是 submitted work；本次没有用仓库原生 PCFF、OPLS 或其他论文的参数替代它。原文摘录保存在参数库的 `inheritance_evidence` 中。版权及再分发权限状态也保留在库中；未复制原始 PDF。

## 实际补充的数值

| 参数类别 | 1995 | 1997 |
| --- | ---: | ---: |
| 键电荷增量 | 5 | 1 |
| LJ9–6 | 2 | 3 |
| 四次键势 | 4 | 1 |
| 四次角势 | 10 | 6 |
| 三项扭转 | 7 | 7 |
| 面外势 | 1 | 0 |
| bond–bond | 8 | 5 |
| bond–angle | 9 | 5 |
| angle–angle | 15 | 0 |
| middle-bond–torsion | 7 | 0 |
| angle–torsion | 3 | 0 |
| angle–angle–torsion | 6 | 0 |

例如，1997 年原表包含：

| 源标签 / 项 | 原始十进制数值 |
| --- | --- |
| `o2-si4` 四次键势：r0, k2, k3, k4 | `1.640, 350.1232, -517.3424, 673.7067` |
| `si4-o2-si4` 四次角势：theta0, k2, k3, k4 | `159.0, 8.5000, -13.4188, -4.1785` |
| `c4-si4-o2` 四次角势 | `114.9, 23.0218, -31.3993, 24.9814` |
| `si4-o2-si4` bond–bond | `41.1143` |
| `si4-o2-si4` bond–angle 左、右系数 | `28.6686, 28.6686` |
| `o2-si4` 电荷增量 | `-0.2225 e` |
| `si4c` LJ9–6：rmin, epsilon | `4.290 Å, 0.131 kcal/mol` |

键势采用 `Σ(k_n·Δr^n)`，角势采用 `Σ(k_n·Δtheta^n)`。平衡角以度记录，角位移及相应系数的量纲按弧度表示。每一列均单独保存单位；LJ 的距离参数为势阱位置 `rmin`，不能直接当成 LJ12–6 的 sigma。尚未执行能量或距离单位变换，转换策略为 `identity_decimal_v1`，保留原表十进制精度。源表的单位排版需结合正文能量公式理解；记录明确注明这一依据。

1995 年的 **8 条角势 k4 原表为空**，JSON 中为 `null`，不会变成零。原文明确印出的零另存于 `explicit_zero_fields`。提取器记录了 3 处经页面图像核对的 OCR 更正，以及 1995 年 Table 5 只选 refined LJ 列、调整字段顺序的操作。

## 查询、缺项审计与执行边界

参数位于 [parameters.json](evidence/sun_siloxane_parameters/parameters.json)，逐行数值证据位于 [source-rows.tsv](evidence/sun_siloxane_parameters/source-rows.tsv)。每条记录含源论文、源文件哈希、PDF 页码/表位置/坐标框、原始 OCR、显式 OCR 更正、标签、原值/转换值、逐列单位和证据校验值。数值重建复用 J43 合约。

```python
from island.forcefields.siloxane_literature import load_siloxane_parameters

supplement = load_siloxane_parameters(
    "docs/evidence/sun_siloxane_parameters/parameters.json",
    "d23244f11d6be8fb683ab01b1f586ee466f60c3b65d10b8da9379eda111961ca",
    references="island_ffref",
)
row, = supplement.lookup("sun1997", "quartic_bond", "o2-si4")
assert row["original_values"][0] == "1.640"
```

`lookup` 只进行原始标签精确查询，不自动反转、套用等价规则或通配符。参数库外部 SHA-256 与代码中的独立内容身份共同约束记录；重新签名也不能授权修改参数、单位、零项或执行标志。

`assess_siloxane_parameters(system, supplement, references=...)` 从显式氢最终图枚举所有 Class-II 相互作用、平衡依赖、电荷增量、质量及策略需求。图包括稳定原子 ID、顺序、组分和分子归属。候选行按来源隔离，分为原文存在但依赖未解决、原表留空、需核实继承的 1995 行、未经验证的通配候选、原文缺行。图审计的候选匹配包含可能的反序，但**不选择系数**，`selected_row` 始终为空。

检查实例为 `[SiH3]O[SiH3]`、hexamethyldisiloxane、PDMS DP1/DP2/DP4 和 `CO[Si](C)(C)C` 阴性对照。聚合物从 `[*:1]O[Si](C)(C)[*:2]` 构建，显式 head cap 为 `[*][Si](C)(C)C`、tail cap 为 `[*]C`，再补全氢。DP 计重复单元数，头端基另含一个 Si；DP1 与 HMDS 化学结构相同但图顺序不同。不会将 methoxy O–C 端基默认为 Si–O–Si 类型。

尚缺或未获验证的依赖包括：

- 同版本 C/H 基础键、角、非键参数及质量/类型定义；例如 HMDS 中 C–H、H–C–H 的完整行。
- 1995 年空白 k4 的含义、继承角势的完整性以及电荷与剩余耦合项是否可继承。
- 最终图所需端键扭转、1–3 bond–bond 等耦合项的参数或明确不适用依据；正文“when required”不足以把所有未列项视为零。
- `si4/si4c` 正文与 Appendix 宽泛类型描述之间的适用范围。
- 完整的异类对组合规则、明确的 1–4 权重，以及全部相互作用的角色/相位映射。
- 独立后端、数值及工作流验证。

因此没有精确电荷向量、参数指派、`PreparedForceField` 或新执行 profile。`require_executable()` 始终抛出携带完整结果的分类阻断异常。数值缺项返回 `parameter_incomplete`；O–C、P、异常价态、带电或自由基等未建立标签的环境返回 `typing_ambiguous`。缺失资产、哈希或证据错误分别返回 `asset_unavailable` / `source_incompatible`。

报告中的 `missing_dependencies` 包含每个尚未完成图绑定指派的请求，其中可能已有可用源行；该计数**不等于论文缺少数值行的计数**。每个请求保留状态、标签、候选来源、空白字段及平衡依赖，便于逐项继续补齐。

## 复现

从仓库根目录执行。数值提取及独立 PDF 重建使用安装 `pypdf` 的 Python；图审计和测试使用支持 RDKit 的项目环境。全部输入均为显式本地路径，不联网。

```bash
export PYTHONPATH=.:src
python scripts/extract_sun_siloxane_parameters.py \
  --references island_ffref --output /tmp/sun-parameters-regenerated
python scripts/audit_sun_siloxane_parameters.py \
  --catalog docs/evidence/sun_siloxane_parameters/parameters.json \
  --sha256 d23244f11d6be8fb683ab01b1f586ee466f60c3b65d10b8da9379eda111961ca \
  --references island_ffref
python scripts/assess_sun_siloxane_parameters.py \
  --catalog docs/evidence/sun_siloxane_parameters/parameters.json \
  --sha256 d23244f11d6be8fb683ab01b1f586ee466f60c3b65d10b8da9379eda111961ca \
  --references island_ffref --output /tmp/sun-siloxane-study
python scripts/assess_sun_siloxane_parameters.py \
  --catalog docs/evidence/sun_siloxane_parameters/parameters.json \
  --sha256 d23244f11d6be8fb683ab01b1f586ee466f60c3b65d10b8da9379eda111961ca \
  --references island_ffref --output /tmp/sun-siloxane-study \
  --reconstruct --require-executable
```

最后一条是执行验收，预期退出码 **1**；纯数据审计和重建成功返回 0。

```bash
python -m pytest -q tests/test_siloxane_literature.py
python scripts/validate_sun_siloxane_parameters.py \
  docs/evidence/sun_siloxane_parameters --references island_ffref
```

验证结果与文件哈希见同目录 `receipt.json` / `manifest.json`，完整图报告在 `study.tar.xz`。没有进行新模型的能量、力、最小化或 MD，也不声称 msi2lmp 一致性。

专项测试 32 项通过。支持环境的全套测试收集 3,601 项；一次重建测试的父/子进程分别导入了测试期间修改前后的候选分类，首次失败，在固定最终代码后单项复跑通过。初次运行和复跑日志均保留在 `validation.tar.xz`，全套汇总保留这一说明。最终版本的完整专项测试也重新通过。Ruff、compileall、五个环境的 pip check 及 J43/J44/J45 历史收据验证通过。工作区补充尚未提交；实现文件与历史基线的精确哈希记录在 manifest 中。

保持 `production_validated=False`、`simulation_readiness="not_established"`。聚合物 Ewald/PME 未验证、后端压力不可用、cutoff recovery 与 long NPT 仍阻断，不增加平衡态或实验吻合声明。
