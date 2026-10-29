# Xiang 2025 外部 OPLS-compatible PDMS：端基与非键约定核查

本次为原文献补充增加独立的 v2 审计层：可根据正文图 7，对一个显式 H、三甲基硅基封端的线性 PDMS 最终图记录 E/M 分区和原始 Qeq 电荷。HMDS 的唯一氧可直接判定；较长链必须显式提供一个端氧的稳定 ID。

这推进了“末端类型与电荷映射”，但仍未证明作者原始脚本中的全部逐原子映射。完整系数单位、LJ 混合、1–4／排除规则没有从本地论文中查明，因此没有生成执行 profile、准备数值后端或运行 MD。原生 PCFF、原生 OPLS-AA 和先前 v1 数据补充均未修改。

来源是 Xiang et al., *Force-Field Benchmark for Polydimethylsiloxane: Density, Heat Capacity, Isothermal Compressibility, Viscosity and Thermal Conductivity*, J. Phys. Chem. B **129** (2025), 1864–1873，DOI `10.1021/acs.jpcb.4c08471`。

| 本地来源 | SHA-256 |
| --- | --- |
| `pcff/jp4c08471.pdf` | `04ecee67243157ada68ef96e6f67d5bf31743e8b207efb9b3e4263a823dbc1fc` |
| `pcff/jp4c08471_si_001.pdf` | `6e86aa64085466e286b5c59878c2288f5f16da1b7a5298d37ba1833a88cd602f` |

这里只核查上述显式本地来源，没有下载作者仓库，也没有借用其他力场的数据。文件所在的 `pcff/` 目录不决定模型归属。

## 图 7 补充了什么

正文 PDF 第 7 页（期刊页 1870）图 7 是嵌入图像，普通文本提取漏掉了其中的表和端基定义。本次直接查看该图，并固定其图像资源及 SHA-256；独立审计器验证图像字节，再从 SI 表 S4 重建同一组电荷。图中端基的文字属于人工视觉转录，不能把图像校验称为自动 OCR 或第二位独立人员的化学审阅。

图 7 将 E 片段写为 `–OSi(CH₃)₃`、`–Si(CH₃)₃`，M 片段写为 `–Si(CH₃)₂O–`。Qeq 列为：

| 元素 | M 电荷 / e | E 电荷 / e |
| --- | ---: | ---: |
| O | −0.6 | −0.578 |
| Si | 0.9 | 0.7 |
| C | −0.342 | −0.386 |
| H | 0.064 | 0.083 |

SI 第 5 页指明 Qeq 来自 Materials Studio 7.0；正文第 7 页明确后续模拟采用 Qeq。图中的 DFT/NPA 列是另一套电荷，不能与 Qeq 混合。

分区规则限定为一条 `Me₃Si–O–(SiMe₂–O)n–SiMe₃` 链：

- 两个端 Si 及其全部甲基 C/H 使用 E 列。
- **一个**明确指定的端氧使用 E 列，其余 O 使用 M 列。
- 内部 Si 及其甲基 C/H 使用 M 列。
- HMDS 只有一个氧，无取向歧义；较长链若未提供 `e_oxygen_site_id`，返回 `typing_ambiguous / explicit_E_oxygen_required`。
- 指定的氧必须与端 Si 相连。内部氧、多个氧、布尔值、字符串和不存在的 ID 均被拒绝。
- 不按最小原子 ID 选择端氧。示例命令显式传入封端构建器的 head ID，并把最终图和该选择一起记录；库本身不隐式读取 head 元数据替用户选择。

按原始十进制数：

```text
M 单元 = 0.9 − 0.6 + 2×(−0.342) + 6×0.064 = 0
两个 E 端基与一个 E 氧 = 2×[0.7 + 3×(−0.386) + 9×0.083] − 0.578 = 0
```

给两个端部附近的氧都套用 E 电荷，会使较长链产生 `+0.022 e` 偏差。这不是可以重新分配电荷修复的问题；实现保留原数值，仅接受显式的分区。两种相反取向会得到不同的图绑定记录，不能宣称它们是同一作者脚本的映射。记录中 `author_atom_mapping_verified=False`。

Si–H、Si–OH、Si–O–C、P、离子、自由基、非显式 H、错误价态、分支或环状硅氧烷均不从本表外推。旧 J45 含 O–C/O–H 端基的 fixture 结果保持原样。

## 哪些约定仍缺来源

| 问题 | 本地证据 | 当前结论及还需要的内容 |
| --- | --- | --- |
| 单位 | 正文 OPLS 截断 1.1 nm；表 S4 的截断参数为 11；SI 图 S1 的总能量轴为 kcal/mol | 与 Å/kcal 的解释一致，但没有完整声明所有力常数单位和方程，不能自动执行换算。需要原始 `units`、键角定义、扭转方程与输出单位。 |
| 端基／电荷 | 图 7 和表 S4 的 Qeq E/M 列 | 可验证上述源一致分区；较长链需显式取向，作者脚本的实际分区仍待核实。 |
| LJ 混合 | 表 S4 只列 4 个同类原子的 LJ 行 | 需要作者的 `pair_modify mix` 或所有异类 `pair_coeff` 及覆盖顺序。不能套用 ISLAND 原生 OPLS 的默认混合规则。 |
| 1–4／排除 | 本地未找到明确的 LJ/Coulomb 缩放与特殊对设置 | 需要作者的 `special_bonds`、排除指令和实际生成拓扑；不能猜测为 0.5。 |
| 其他执行依赖 | 表 S4 列 harmonic 键角、opls 二面角 | 还需质量表、improper/inversion 的适用性或来源参数、Moltemplate 与原始输入版本、长程配置。 |

即使将来确认 `units real`，仍应按作者实际方程核实 harmonic 力常数的二倍因子、Å/nm 转换及 OPLS 四项的 1/2 系数，不能只给数值换单位标签。本次没有输出转换后的可执行系数。

补齐剩余证据的明确来源是正文参考文献 22 的作者仓库 `https://github.com/COMO-lab/PDMS_benchmark`，或用户提供的同版原始脚本。在获得网络读取授权前，本次没有访问该仓库。核实时应固定提交、全部 include 文件及生成文件的哈希，沿实际覆盖顺序读取，检查其系数是否与表 S4 相同；当前仓库版本不能自动视为论文发表时版本。

## 接口与复现

新增数据身份为 `island_external_opls_compatible_pdms_xiang2025_conventions_v2`，不属于执行 profile。`profile_identity=None`、`execution_authorized=False`、`acceptance_gate=1`。

```python
from island.forcefields.oplsaa.pdms_conventions import assess_opls_pdms_conventions

result = assess_opls_pdms_conventions(
    system, supplement, references="island_ffref",
    e_oxygen_site_id=explicit_terminal_oxygen_id,
)
print(result.payload["exact_charge_vector"])  # 原始十进制字符串；歧义图为 None
result.validate_integrity(system, references="island_ffref")
result.require_executable()  # 始终抛出携带完整结果/缺失依赖的 classified exception
```

`supplement` 的加载方法见 [v1 文献补充](opls_pdms_literature_supplement.md)。本层重放验证 v1 的固定身份、图、稳定 ID、原子顺序、组分／分子身份、化学状态、端氧选择和证据。独立 SHA-256 加载失败、源文件改变和重新签名的伪造单位／规则／执行权限均被拒绝。坐标移动不改变化学身份。

在支持 RDKit 的环境中：

```bash
export PYTHONPATH=.:src
../island-validation/phase4e2-env/bin/python scripts/assess_opls_pdms_conventions.py \
  --selection docs/evidence/phase_4j45/selection.json \
  --catalog docs/evidence/phase_4j44/source-audit.json \
  --references island_ffref --output /tmp/opls-pdms-conventions-v2
```

该命令仅生成数据审计；添加 `--require-executable` 必须返回 1。输出包含 HMDS、正确 TMS 封端 DP2/DP4、未指定取向的 DP2/DP4、SiH₃ 二硅氧烷负对照。这里 DP 指 repeat 数，显式头端基额外贡献一个 Si，不能混同于原 J45 的端基或 Si 总数。

```bash
python scripts/audit_opls_pdms_convention_sources.py \
  --references island_ffref \
  --evidence docs/evidence/opls_pdms_conventions_v2/source-evidence.json
../island-validation/phase4e2-env/bin/python -m pytest -q tests/test_opls_pdms_conventions.py
../island-validation/phase4e2-env/bin/python scripts/validate_opls_pdms_conventions.py \
  docs/evidence/opls_pdms_conventions_v2 --repo . --references island_ffref
```

独立审计器只导入标准库与 pypdf，不调用生产参数解析器。迁移加载使用 `load_opls_pdms_conventions(path, independent_sha256, supplement=..., references=...)`；具体文件校验值列在新证据目录的 `report.json` 和 `manifest.json` 中。

本次最终代码的针对性与相关回归共 232 passed、1 skipped；其中本层 focused 43 passed。运行 Ruff、compileall、五个环境 pip check、J44/J45 receipt/manifest 验证及历史字节校验。该小范围变更未重新运行全项目 supported suite，不将历史全套结果当作本次结果。详细命令和结果见 [验证记录](evidence/opls_pdms_conventions_v2/validation.json)。

保留 `production_validated=False`、`simulation_readiness="not_established"`。没有 PDMS 能量／力、MD 或检查点验收；Ewald/PME、压力、cutoff recovery、long NPT 限制不变。
