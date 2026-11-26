# J107B：来源约束的覆盖扩展与完整重测

本阶段增加显式 `opls-polymer-v2`，只修改分型规则，`.prm` 数值参数和 `.define` 来源 residue 保持字节一致。新增卤代烷、饱和中性胺和支化芳醚环境；不使用通用零参数或自动近似参数补齐 PCFF/OPLS 缺项。PCFF 继续使用已有 carboxamide-v1 来源。

## 固定面板

982 输入 × 2 力场 × 完整图/旧入口，加 34 挑战输入 × 2 力场，共 **3,996/3,996** 项记录齐全。

| 力场 | 入口 | J107 历史成功 | 本次成功 | 本次严格立体化学 |
|---|---|---|---|---|
| PCFF | graph | 880/982 | 904/982 | 904/982 |
| PCFF | repeat | 899/982 | 899/982 | 865/982 |
| OPLS | graph | 440/982 | 448/982 | 448/982 |
| OPLS | repeat | 406/982 | 408/982 | 400/982 |

挑战图：PCFF **29/34**，OPLS **25/34**。这些都是参数分配结果，不能视作全部完成 MD 或平衡。

相对 J107，共 38 个“案例/力场/入口”由失败变成功，**0 个旧成功案例回退**，4 项成功案例的类型计数发生变化（完整记录见 comparison.json）。J107 原任务时限为 120 秒；本次为 180 秒，图路径每个后端命令 30 秒，因此总增益包含 J107A 图处理修复、执行预算和 v2 profile 的共同作用，不能全归为新规则。

对全部新增成功的 OPLS 案例，以相同当前代码和时限重新使用 v1 来源：**8 项**在 v1 下仍失败、v2 成功；**6 项**在 v1 下也成功。逐案例源哈希、结果和对照日志均归档；这不是额外完整重测所有 v1 输入后的全覆盖数字。

## 来源与数值验证

8 个独立来源 residue（MEAN、ETAN、DMAN、TMAN、ETCL、ETBR、ETIO、ANIS）的原子类型与电荷逐项一致。源元素与键型解析复用可靠解析器，避免从类型名称猜测芳香性。5 个开发探针单独标记，不再误称为留出集。

4 个未用于规则创建的较长模型完成原子置换与有限差分验证，最大力梯度误差为 **0.00253179 kcal/(mol Å)**。PVC、PVBr、聚胺和芳醚各完成一个 **300 fs** 的周期 NVT/NPT/生产配方。它们验证可运行性，不证明平衡或参数的实验精度。

## GAFF / GAFF2

固定 982 个 DP4、氢封端完整图，GAFF 与 GAFF2 均 **982/982** 成功，分别与独立调用的原生 RadonPy 比较全部导出类型、参数与 Gasteiger 电荷：**各 982/982 一致**。该面板直接测试无坐标完整图；不把图分配当作三维构建、装箱或 MD 覆盖。

另有每个力场 12 个实际构象和原生拓扑导出比较，包括 PE、手性 PP、PEG、PET、PA6、cis 二烯、卤代链、芳醚及明确配置的嵌段/支链，全部通过。采用的是 RadonPy 自身参数与经验补项行为，没有声称 RESP/ESP、反应力场或任意化学空间均已验证。

相关回归 **125 passed**。OPLS 的 P/Si、部分硫/芳香氮、离子及缺失增量/扭转环境仍需来源支持；PCFF 也保留所有真实缺项。不能据此宣称 PCFF/OPLS 已完整覆盖任意 PSMILES。

## 复现

```bash
python -m island.forcefields.emc_profiles opls-polymer-v2 work/opls-v2
PYTHONPATH=.:src python scripts/validate_j107b_coverage.py work/source-checks --field-dir work/opls-v2
PYTHONPATH=.:src python scripts/benchmark_final_graphs.py \
  --pcff-field work/pcff-carboxamide-v1 --opls-field work/opls-v2 \
  --out work/paired --workers 6 --timeout 180 --backend-timeout 30
PYTHONPATH=.:src python scripts/benchmark_gaff_native_graphs.py work/gaff-native
PYTHONPATH=.:src python scripts/validate_gaff_graph_coverage.py work/gaff-architectures.json
PYTHONPATH=.:src python scripts/validate_j107b_dynamics.py work/dynamics work/opls-v2
```

证据与各归档成员校验和见 `docs/evidence/phase_4j107b/`；J107/J107A 历史证据未改写。
