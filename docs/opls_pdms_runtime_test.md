# 当前 ISLAND OPLS 能否模拟 PDMS：实际入口测试

结论：**当前不能通过 ISLAND 的 OPLS 公共入口模拟 PDMS。** 本次实际调用 `prepare_forcefield`，不是仅根据文献审计推断。固定的原生 Foyer OPLS XML 没有 Si atom type，原生 OPLS 的元素／价态检查也没有 Si。因此准备阶段终止，未创建 PDMS 的 `PreparedForceField`、求值器或 MD 系统。

先前文献提取和作者源审计仍有效，但它们是外部数据记录，没有接成一个经过数值验收的执行 adapter。本次不修改原生 OPLS、PCFF、历史数据或已封存收据。

## 真实入口结果

统一使用固定 XML `c78ccb763cda33e3456a3f8b4ed8a8f0361b10c92c33aa02f2a9abf618c163e4`，Foyer 源 revision `dd2f6eaa0ec271432ccd0c5c729f17a3ea3364bd`，验证过安装代码的全部既有 pins。环境为 Foyer 1.2.0、OpenMM 8.6.1、Reference 平台、RDKit 2026.3.6。

| 输入 | 原子数 | OPLS 准备结果 |
| --- | ---: | --- |
| `[SiH3]O[SiH3]` | 9 | `typing_blocked` |
| 六甲基二硅氧烷 HMDS | 27 | `typing_blocked` |
| PSMILES、TMS 封端 DP1 | 27 | `typing_blocked` |
| PSMILES、TMS 封端 DP2 | 37 | `typing_blocked` |
| PSMILES、TMS 封端 DP4 | 57 | `typing_blocked` |

聚合物使用 `[*:1]O[Si](C)(C)[*:2]`，显式 `head_end_group="[*][Si](C)(C)C"`、`tail_end_group="[*]C"` 并补全 H。DP 指 repeat 数；头端基额外贡献一个 Si，故 DP2 有 3 个 Si，DP4 有 5 个 Si。DP1 和 HMDS 化学组成相同，分别验证 SMILES 与 PSMILES 两条构建路径，不能计为两个不同化学域。

这些中性显式 H 图都通过现有 `island_external_pdms_final_graph_v1` 的 Si/C/O/H 化学检查。另行测试通用 `final_graph` 时，它因历史合同没有 Si 而报告 `ValidationError: Invalid valence at site ...: 4.0`；这不意味着上述分子实际存在错误 Si 价态。

测试保留该错误，再直接测试允许不显式传入通用图记录的现有 `prepare_forcefield(system, ForceFieldRequest("oplsaa", OPLSOptions(xml)))` 路径。它仍在原生 OPLS 元素检查处抛出 `OPLSAssignmentError: Unsupported or inconsistent element/atomic number`。`typing_blocked` 是测试报告对该原生异常的分类；没有声称原生异常自带此 status 字段，也没有把 Si 专用图伪装成通用最终图。

## 环境与数值引擎对照

在同一安装和 XML 下，对照真实通过 `prepare_forcefield → create_evaluator → evaluate → minimize_geometry`，没有替换参数或 monkeypatch 源 pins：

| 对照 | 初始能量 / kJ mol⁻¹ | 最小化后 / kJ mol⁻¹ | 最小化结果 |
| --- | ---: | ---: | --- |
| 乙烷 | 30.924085 | 8.165926 | force_converged，最终能量独立复核通过 |
| 正丁烷 | 55.114967 | 8.671195 | force_converged，最终能量独立复核通过 |
| 苯 | 43.978486 | 27.951823 | force_converged，最终能量独立复核通过 |

全部初始力为有限数值。最小化使用现有默认选项，包括 fmax 阈值 0.1 kJ mol⁻¹ Å⁻¹；未扩大容差。三种对照结果证明 OPLS/OpenMM 环境能工作，不能替代 PDMS 的能量或力验证。对照和 PDMS 均未运行 MD。

另行重放作者源 v3 审计的执行门，仍返回 `scientifically_unqualified`、退出码 1。原因保持为：作者输入的 DFT/NPA 电荷＋CHARMM 切换 LJ 与论文 Qeq＋截断 LJ 分支不同，许可与新图类型规则来源未确认，也未完成独立数值和工作流验收。此次不选择分支、不用默认 OPLS 非键实现替换作者模型。

## 复现

本次只增加测试工具 [probe_opls_pdms_runtime.py](../scripts/probe_opls_pdms_runtime.py) 和新证据。显式使用本地固定安装，不联网、不发现或借用其他参数：

```bash
export PYTHONPATH=.:src
OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 \
../island-validation/phase4g1-env/bin/python scripts/probe_opls_pdms_runtime.py \
  --xml ../island-validation/phase4g1-env/lib/python3.11/site-packages/foyer/forcefields/xml/oplsaa.xml \
  --output /tmp/opls-pdms-runtime-probe.json --require-pdms-executable
```

目标文件必须不存在。预期退出码 **1**：`boundary_probe_passed=true` 表示预期阻断和有机对照符合实测，`native_PDMS_executable=false` 才是 PDMS 的执行结论。去掉最后一个参数时，0 只表示诊断探针通过。探针若发现未知异常或意外的 Si 准备成功，不会继续构造 PDMS 数值引擎。

有真实 XML 的原生 OPLS 测试：

```bash
ISLAND_OPLSAA_XML=../island-validation/phase4g1-env/lib/python3.11/site-packages/foyer/forcefields/xml/oplsaa.xml \
OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 \
../island-validation/phase4g1-env/bin/python -m pytest -q \
  tests/test_oplsaa_parameters.py::test_real_source_noncontiguous_permutation
```

作者源／PDMS 边界测试使用已下载并固定的作者文件；不会下载参数：

```bash
ISLAND_PDMS_AUTHOR_SOURCES=/tmp/island-pdms-author-fetch-verification \
ISLAND_LAMMPS_DOCS=/inspire/ssd/project/sais-suiren-foundation-model-prediction/public/yh/lammps/lammps/doc/src \
OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 \
../island-validation/phase4e2-env/bin/python -m pytest -q \
  tests/test_opls_pdms_author.py tests/test_opls_pdms_conventions.py
```

本次 **62＋1＝63 passed，0 skipped**。包含作者源与记录篡改拒绝、新进程迁移、以及真实 Foyer XML 的非连续 atom ID／能量／力映射测试。原生对照的能量与收敛最小化另见 [runtime-report.json](evidence/opls_pdms_runtime_probe_v1/runtime-report.json)。作者执行门的精确命令见 [validation.json](evidence/opls_pdms_runtime_probe_v1/validation.json)。

Ruff、compileall、diff／新文件空白检查及新收据／manifest 的字节校验通过，3,438 个既有文件的字节不变。本次是实际入口诊断，没有重跑全项目 full suite，也没有把历史验证计数算作本次结果。

当前需要的是一个独立的、明确选定文献分支的外部 OPLS-compatible PDMS adapter，以及源绑定类型／电荷／全部势函数的独立数值验收；修改原生 OPLS 的允许元素列表并不能补齐这些内容。保持 `production_validated=False`、`simulation_readiness="not_established"`，Ewald/PME、压力、cutoff recovery、long NPT 限制不变。
