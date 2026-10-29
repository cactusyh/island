OPLS 现在有一个独立的、可继续扩展的参数入口。新元素和新化学环境由显式外部参数源的 SMARTS 规则和参数行定义，无须逐个加入原生元素白名单。原生 Foyer XML、原生 OPLS 身份、J17–J45 以及已有 PDMS 文献审计保持原样。

这次完成的是外部参数的准备、不可变分配、OpenMM 求值、最小化及有限非周期动力学连接。没有新增 PDMS 物理参数，也没有把作者的周期 PDMS 模型改成有限分子模型。当前没有可执行的 PDMS 参数包。

使用方式：

```python
from island.chemistry import from_smiles
from island.forcefields.extensible import (
    ForceFieldRequest, ExternalOPLSOptions, ExternalOPLSError,
    prepare_forcefield, create_evaluator,
)
from island.minimization import minimize_geometry

system = from_smiles("CC", random_seed=2026)  # 也可传入已经封端、补氢的 PSMILES 最终图
request = ForceFieldRequest(
    "external_opls",
    ExternalOPLSOptions("/explicit/path/source.json", "<independent SHA-256>"),
)
prepared = prepare_forcefield(system, request)
# prepared.payload 保留完整诊断；未准备完成时 profile_identity 为 None。
try:
    evaluator = create_evaluator(system, prepared)
except ExternalOPLSError as error:
    print(error.status, error.diagnostics)
    # error.result 是完整准备记录。此时不能继续数值计算。
else:
    energy_and_forces = evaluator.evaluate()
    minimum = minimize_geometry(system, evaluator)
```

显式扩展命名空间是 `island.forcefields.extensible`。其中的原生 `oplsaa`、`pcff`、`gaff` 等请求继续交给原有实现。现有工作流配置不会自动选择外部参数；外部请求需要显式选择 `external_opls`。本入口返回独立的 `ExternalOPLSPreparation`，不伪装成历史 `PreparedForceField`。

后续加入参数时，准备一份完整的外部 XML 源和配套 JSON 描述。XML 可以使用 Foyer 的 `AtomTypes` SMARTS、覆盖关系、谐键、谐角、六系数 RB 或四系数 OPLS proper torsion，以及显式电荷和 LJ 行。四系数约定为 `1/2 Σ k_n [1+(-1)^(n+1) cos(n φ)]`；RB 保留全部六个系数，包括常数项。使用几何 LJ 混合，1–2/1–3 排除，1–4 LJ/Coulomb 缩放必须在源中分别写出。缺行不会填零；源中明确存在的零会连同原文证据保存。

描述文件的身份为 `island_external_opls_source_v1`；每套参数必须有独立的 `island_external_opls_<name>_vN` 身份。描述包含 XML 相对路径和哈希、引用/DOI、类型电荷方法、模型范围、用途授权的文件哈希与逐行原文、转换说明和科学限制。`write_source_manifest()` 需要调用者显式提供 `model`，不会默认选择截止、长程或 improper 约定。完整结构见 [source.py](../src/island/forcefields/opls_external/source.py)，实际构造例子见 [probe_extensible_opls.py](../scripts/probe_extensible_opls.py) 的有限 CH 对照来源。

数值输入单位固定为 OpenMM XML 的 nm、kJ/mol、rad、e；内部分配记录显式保存 Å 转换、原始十进制文本、所选源行、候选行、通配符选择、稳定原子 ID、顺序、组件、分子身份、源质量和完整电荷。逐组件电荷必须精确中性，不修正电荷余差。图质量与源质量不一致时拒绝，不替换图或源数值。类型/参数匹配采用正反角色对称、精确行优先、冲突候选拒绝；这是一套保守的版本化外部解析规则，不宣称对任意 XML 都与 Foyer 首行选择完全相同。

v1 可执行范围是显式氢、无同位素、逐原子中性、闭壳层、有限非周期 `NoCutoff`、几何 LJ、受支持的谐键/角与 RB/OPLS proper torsion，且源必须明确声明该化学域不需要 impropers。带显式 improper、其他势函数、Class-II 耦合、虚拟位点、成对覆盖、额外排除规则、带电/自由基、切换势和周期 PPPM/PME 不能由本适配器执行。源码可以保存被阻止的模型声明，但不得把它转换成受支持的模型。未来增加同一受支持势函数下的参数主要增加来源数据；增加新势函数或长程方法仍需要新的实现和独立验证。

不自动合并原生表与文献表，不做条目覆盖或来源发现。若原文明确授权基底与扩展组合，应在完整外部源及其文献证据中明确列出，再形成新的整体哈希；不能用缺项时回退来实现组合。本层验证来源/记录的一致性，不能仅凭哈希、许可证声明或成功求值授予科学有效性。文献提取仍须通过已有 J43/J44/J45 合同和独立审计。

准备记录使用 `save_preparation()` 保存；`load_preparation()` 需要独立的记录哈希、最终图、显式源路径及源描述的独立哈希，并重做源、图、类型、电荷和每个相互作用的核对。每次求值前再次核对实际源文件，不能只验证曾经缓存的字节。搬移目录保持相对文件名时身份稳定。坐标更新和项目产生的坐标历史不属于化学身份，修改原子 ID、插入顺序、拓扑、组件或化学构建元数据则失效。

乙烷和丁烷使用未改动的真实 Foyer 源完成了以下检查：

| 检查 | 结果 |
| --- | --- |
| 原生类型/电荷/质量/顺序、相互作用数量及选中源行 | 一致；独立 DOM 重建 9/14 条不同源行 |
| 原生能量差 / 最大力差 | ≤ 3.6×10⁻¹⁴ kJ/mol / 6.5×10⁻¹³ kJ/(mol·Å) |
| 独立 NumPy 分量能量与 OpenMM | 最大差 ≤ 3.2×10⁻¹⁴ kJ/mol |
| 全原子、全笛卡尔分量有限差分 | 两种步长；最大误差 ≤ 1.21×10⁻⁵ kJ/(mol·Å) |
| 最小化 | 均收敛，返回能量重新核验 |
| 有限 Langevin 控制 | 4 步，2+2 步重启的最终状态和完整 RNG 精确相同 |
| checkpoint 序列化往返 | 字节精确；分段与连续运行的执行历史/计数不同，不声称整个 checkpoint 相同 |

纯解析测试另用明确标注的人工软件数据检查 Si 和 F 的来源驱动类型规则。这些测试数据没有物理授权，不作为 PDMS/含氟力场参数交付。测试还检查缺少键/角/扭转、冲突通配符、电荷/单位/缩放/混合变更、失效图、重签名伪造记录、来源变更、搬移和所有阻止状态的新进程重建。

PDMS 的实际状态：五个最终图（SiH₃ 二硅氧烷、HMDS、封端 DP1/DP2/DP4）都能通过新的图检查。没有外部来源时返回 `asset_unavailable`；使用仅声明 CH 的对照来源时，逐项列出缺少的 O/Si 范围并返回 `typing_blocked`。不会创建执行配置或数值引擎。

现有作者来源审计仍为 `scientifically_unqualified`：DFT/NPA 输入不能充当论文 Qeq 模型；作者实际使用 `lj/charmm/coul/long 9 11`、几何混合、0/0/0.5 特殊键及 PPPM。授权/再分发证据、对新图的原始类型规则、端基电荷映射与独立后端资格仍需补齐。本适配器没有用 `NoCutoff` 替代该模型，也没有把旧审计变成执行许可。

复现命令：

```bash
# 实际来源正对照；--work 和 --output 使用新路径
PYTHONPATH=.:src OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 \
../island-validation/phase4g1-env/bin/python scripts/probe_extensible_opls.py \
  --xml ../island-validation/phase4g1-env/lib/python3.11/site-packages/foyer/forcefields/xml/oplsaa.xml \
  --license ../island-validation/foyer-4g1-source/LICENSE.rst \
  --work /tmp/opls-extension-probe --output /tmp/opls-extension-probe.json

# 不提供来源，PDMS 必须返回非零并保存分类诊断；不启动求值
PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python scripts/external_opls.py \
  --psmiles '[*:1]O[Si](C)(C)[*:2]' --dp 2 \
  --head-end-group '[*][Si](C)(C)C' --tail-end-group '[*]C' \
  --output /tmp/opls-pdms-blocked --evaluate

# 完整专门测试，在固定的 Foyer 环境执行所有数值正对照
PYTHONPATH=.:src OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 \
ISLAND_OPLSAA_XML=../island-validation/phase4g1-env/lib/python3.11/site-packages/foyer/forcefields/xml/oplsaa.xml \
ISLAND_FOYER_LICENSE=../island-validation/foyer-4g1-source/LICENSE.rst \
../island-validation/phase4g1-env/bin/python -m pytest -q tests/test_extensible_opls.py
```

证据和检验摘要在 [opls_extensible_v1](evidence/opls_extensible_v1/receipt.json)。保留 `production_validated=False`、`simulation_readiness="not_established"`。聚合物 Ewald/PME、压力、cutoff recovery、长 NPT 均未因此获得验证；本阶段没有 LAMMPS 导出/执行、实验符合性或平衡态声明。

组件电荷核算使用精确有理数，原始十进制字符串仍完整保存。补充的回归例子使用人工的 `-0.3 + 3×0.1 + 1e-40`：它必须被拒绝，不能因为默认 Decimal 精度将结果舍入为 `0E-28` 就获准执行。这个修正仅作用于新的外部适配器。

常规检查命令：

```bash
PYTHONPATH=.:src OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 \
../island-validation/phase4e2-env/bin/python -m pytest -q
../island-validation/phase4e2-env/bin/python -m ruff check .
../island-validation/phase4e2-env/bin/python -m compileall -q src scripts tests
for py in \
  ../island-validation/phase4e2-env/bin/python \
  ../island-validation/phase4g1-env/bin/python \
  ../island-validation/phase4j32-analysis-env/bin/python \
  ../miniforge3/envs/md_rg/bin/python \
  ../island-validation/phase4j36-native-env/bin/python; do
  "$py" -m pip check
done
git diff --check
PYTHONPATH=.:src ../island-validation/phase4e2-env/bin/python \
  scripts/validate_phase_4j44_receipt.py docs/evidence/phase_4j44
PYTHONPATH=.:src ../island-validation/phase4e2-env/bin/python \
  scripts/validate_phase_4j45_receipt.py docs/evidence/phase_4j45
# 新回执：使用独立保存的 manifest SHA-256；本次值保存在 validation.json
PYTHONPATH=.:src ../island-validation/phase4e2-env/bin/python \
  scripts/validate_extensible_opls.py docs/evidence/opls_extensible_v1 \
  --repo . --sha256 '<independent manifest SHA-256>'
```

nm→Å 和相应键常数的十进制换算也使用精确有理数后再转为后端浮点值，避免进程中的 Decimal 精度设置改变参数。低精度上下文和清空缓存后的独立重建已纳入回归。
