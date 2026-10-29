# OPLS-compatible PDMS 文献参数补充

已在 `island.forcefields.oplsaa.pdms_literature` 增加 Xiang 2025 表 S4 的独立文献数据接口，可以查询、导出、按独立 SHA-256 加载，并对最终图运行原 J45 完备性审计。
**目前补充的是可追溯的原始参数记录，尚不能自动给 PDMS 分配可执行的 OPLS 力场。**

来源：Xiang et al., *Force-Field Benchmark for Polydimethylsiloxane: Density, Heat Capacity, Isothermal Compressibility, Viscosity and Thermal Conductivity*, J. Phys. Chem. B **129** (2025), 1864–1873，DOI `10.1021/acs.jpcb.4c08471`。
表 S4 位于本地 `pcff/jp4c08471_si_001.pdf` 第 6 页；正文讨论电荷方法位于 PDF 第 6–7 页。参数置于 `pcff` 文件夹不改变其 OPLS-compatible 归属。

本次复用 J44/J45 的不可变源记录及提取结果，不重写历史证据。

## 已补充内容

| 类别 | 行数 | 来源内容 |
| --- | ---: | --- |
| 非键 | 4 | O、Si、C、H 的同类原子 LJ 行 |
| 键 | 3 | Si–O、Si–C、C–H，原文 `bond_style harmonic` |
| 角 | 6 | H–C–H、Si–C–H、Si–O–Si、C–Si–C、C–Si–O、O–Si–O |
| proper torsion | 4 | 原文 `dihedral_style opls`，保留全部四列 |
| 部分电荷 | 4 | O、Si、C、H，各保留 `inside monomer` / `end monomer` 两列 |

总计 17 行系数和 4 行双列电荷记录；电荷列不指派给具体原子，不构造或归一化电荷向量。
每行保留原始十进制字符串、标签、原文、PDF 页码/提取行号、DOI、源文件 SHA-256、历史证据校验值及新记录校验值。
明确印出的零值保留其列索引；缺少的相互作用仍返回空查询结果，不能补零。

数据身份：`island_external_opls_compatible_pdms_xiang2025_supplement_v1`。
它是 data identity；`profile_identity=None`、`execution_profiles=[]`、`execution_authorized=False`。
没有注册原生 OPLS-AA 或外部执行 profile，也没有覆盖 Foyer XML、原生 PCFF、GAFF 或历史类型/电荷数据。

## 为何尚不能用于 MD

原文的系数行并未建立完整的可执行来源约定。本次保留 J45 的全部未决项：

- 数值表的完整单位、力常数定义及转换，不能由数值大小或 OPLS 名称猜测。
- 内部与末端单元的精确逐原子类型、电荷列映射；未提供的 Si–H、O–H、O–C 环境不能继承甲基 PDMS 参数。
- 不同类型间的 LJ 混合、1–2/1–3 排除、1–4 LJ/Coulomb 缩放及特殊对。
- improper/inversion 的适用性、完整质量/同位素策略和同版参数来源。
- 扭转方程/转换的独立重建、完整非键/长程后端配置、配套来源及使用依据。

`source_units`、`converted_values`、`conversion_code_version` 均保持 `None`；没有进行隐式 kcal/kJ、Å/nm 或 Fourier/RB 换算。
论文明确使用的 OPLS cutoff 1.1 nm 和 PPPM 仍保留于原 J45 证据，但不据此推导其他缺失策略。
不能从原生 OPLS 的默认设置、另一论文或 COMPASS 表借入这些约定。

六个历史最终图的审计继续阻断：HMDS 为 `parameter_incomplete`；SiH3 二硅氧烷、原 J45 使用含 O–C/O–H 端基的显式 H repeat/DP2/DP4，以及混合 Si–O–C 图为 `typing_ambiguous`。
这些结果针对相应的最终图与端基，不是对所有端基方案都作相同结论。

## Python 使用

```python
from pathlib import Path
from island.forcefields.literature_closure import DomainSelection
from island.forcefields.oplsaa.pdms_literature import load_opls_pdms_literature

selection = DomainSelection(
    Path("docs/evidence/phase_4j45/selection.json").read_text(),
    Path("docs/evidence/phase_4j44/source-audit.json").read_text(),
)
supplement = load_opls_pdms_literature(
    "docs/evidence/opls_pdms_literature_v1/parameters.json",
    "a8c1e757fb2558820dd2f1f9d2c8d0c9db73dae4fe7069c2709c5fcc4f42e016",
    selection=selection,
    references="island_ffref",
)
row = supplement.lookup("bond", ["o", "si"])[0]
print(row["original_values"])  # ['374', '1.64']；未赋予可执行单位
print(row["source_evidence"]["location"])
```

查询严格使用原文标签、大小写和顺序；不会自动反转、匹配通配符或推导原子类型。
已加载的 supplement 仍是数据记录，不能作为原生 `ForceFieldRequest("oplsaa", ...)` 的 options。

对 ISLAND 最终图可调用：

```python
result = supplement.assess(system, references="island_ffref")
print(result.payload["status"])
print(result.payload["missing_rows"])
result.require_executable()  # 抛出带完整结果的 DomainClosureBlockedError
```

审计每次重新验证所引用的两份本地 PDF，不扫描或下载参数。缺少文件为 `asset_unavailable`，源被更改为 `source_incompatible`。
结果沿用 J45 的最终图、稳定 ID、顺序、组分/分子归属、相互作用请求及校验机制。

## 命令行复现

在支持 RDKit 的环境中执行：

```bash
export PYTHONPATH=.:src
python scripts/supplement_opls_pdms_literature.py \
  --selection docs/evidence/phase_4j45/selection.json \
  --catalog docs/evidence/phase_4j44/source-audit.json \
  --references island_ffref --output /tmp/opls-pdms-literature
```

导出成功退出 0，仅表示数据导出完成。添加 `--require-executable` 必须退出 1；它不会启动数值计算。
独立固定的输出文件校验值为上例 `a8c1e757...42e016`；验证时须从此独立记录取得完整值，不能信任待验证文件自报的哈希。

```bash
python scripts/supplement_opls_pdms_literature.py \
  --selection docs/evidence/phase_4j45/selection.json \
  --catalog docs/evidence/phase_4j44/source-audit.json \
  --references island_ffref \
  --reconstruct /tmp/opls-pdms-literature/parameters.json \
  --sha256 a8c1e757fb2558820dd2f1f9d2c8d0c9db73dae4fe7069c2709c5fcc4f42e016
```

独立 PDF 重建沿用不导入生产参数解析器的原审计器，需要本地安装的 `pypdf`：

```bash
python scripts/audit_phase_4j45_sources.py --references island_ffref \
  --catalog docs/evidence/phase_4j44/source-audit.json \
  --selection docs/evidence/phase_4j45/selection.json
python -m pytest -q tests/test_opls_pdms_literature.py
```

源哈希、导出结果和验证信息见 [evidence/opls_pdms_literature_v1](evidence/opls_pdms_literature_v1/manifest.json)。
本次只验证数据与阻断边界；没有 PDMS 能量/力、MD、检查点或科学验收结果。
保持 `production_validated=False`、`simulation_readiness="not_established"`；Ewald/PME、压力、cutoff recovery 和 long NPT 限制不变。
