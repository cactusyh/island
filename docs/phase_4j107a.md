# J107A：完整图兼容性修复

从 J107 `11e0651` 新建 `codex/j107a-graph-compatibility`。实现提交 `bef9b63`，本报告与证据另行提交。J107 已推送 GitHub，其原始覆盖率证据未改写。

## 结果与归因

J107 的 25 个仅旧入口成功案例在本次最终重测中 **25/25 分配成功**。包括 PCFF 的 20 个环连接拒绝、3 个任务超时、1 个扭转缺项，以及 OPLS 的 1 个分型失败。参数文件未修改，没有补零、替代数值或退回重复单元入口。

在隔离的 J106 `eec859a` 工作树中，以同样输入、DP4、种子 2026、参数目录重放这 25 例，均未成功。这说明这批入口差异早于 J107；J106/J107 的超时与后端报错可能不同，不能把失败类别变化等同于化学变化。基线每任务限时 120 秒。

## 实现

- **等价图遍历重试**：默认先使用原输入遍历；仅在芳香图的已识别环处理/分型/缺参错误后，尝试规范遍历、反序及固定种子 1–5 的重排，最多八种表示。每种表示从相同完整图生成，保留真实原子映射；成功仍须通过元素、全部键连接、组分、电荷及位点完整性检查。返回的是经验证的输入坐标，指定立体化学不靠 EMC 的临时构象表达。
- **按实际环尺寸限制搜索**：以 RDKit 环基中的最大环尺寸与 8 的较大者设置 EMC 深度；本批十二元环从无界自动搜索改为深度 12。后端明确要求进一步深度时仍保留自动重试。`ring_scan_depth` 和 `automatic_ring_depth` 区分实际策略。环基深度不是任意复杂图必定成功的证明。
- **记录尝试与阶段**：`serialization`、`assignment_attempts` 保存选中的表示、每次失败原文及耗时；`timings_seconds` 分开记录输入几何与后端耗时。异常携带 `assignment_stage`；无法确定阶段的外层进程超时不冒充参数缺失。
- **可选后端时限**：`backend_timeout` 限制每次脚本生成/EMC 执行，超时终止对应进程组。明确提供时限时，芳香图可在超时后尝试下一等价表示；这是单次外部命令上限，不是整个请求的总预算。复测运行器另设整个任务 180 秒上限。
- **缓存隔离**：显式表示策略属于缓存键，不能用另一种策略的缓存替代对照实验。

```python
from island.forcefields.assignment import assign_polymer_parameters

assignment = assign_polymer_parameters(
    psmiles, "pcff", dp=4, field_dir="work/pcff",
    backend_timeout=30,
)
print(assignment.serialization, assignment.ring_scan_depth)
print(assignment.assignment_attempts, assignment.timings_seconds)
```

显式 `serialization="input"` 禁用遍历重试；也可指定 `canonical`、`reverse`、`permuted-1` 至 `permuted-5`，用于独立对照。该选项不改变用户请求的化学结构或端基。

## 验证

- 最终版本 J107A/J107/J106/J104/J94 共 **120 项测试通过**。覆盖序列化与组分/键映射、超时、缓存策略隔离、环尺寸选择及既有图分配行为。
- PI736、PI936、PI837、PI1047 四组结构分别以原原子顺序及独立种子 41 重排运行：全部类型、电荷一致；原子映射还原后的逐原子力最大差约 `1.14e-13 kcal/(mol Å)`，能量分项最大差约 `1.03e-12 kcal/mol`。对照使用实际 LAMMPS、直接 Ewald、同一输入构象。
- 最终 25 例全部重新执行，使用 `pcff-carboxamide-v1` / `opls-polymer-v1`，四个工作进程、每个后端命令 30 秒、整个任务 180 秒上限；保留每例参数、数据、原始错误及来源绑定清单。
- 初始仅 20/25 成功的探索记录另存于压缩包 `final/`，最终 25/25 记录位于 `release/`，不得混作同一次运行。

```bash
python scripts/validate_graph_compatibility.py \
  --pcff-field work/pcff --opls-field work/opls-polymer \
  --out work/j107a --workers 4 --timeout 180 --backend-timeout 30
ISLAND_ASSIGNMENT_LAMMPS=/path/to/lmp python scripts/check_graph_retry_parity.py \
  work/j107a-parity --pcff-field work/pcff --opls-field work/opls-polymer
python -m pytest -q tests/test_phase_4j107a.py tests/test_phase_4j107.py \
  tests/test_phase_4j106.py tests/test_phase_4j104.py tests/test_phase_4j94.py
```

[证据目录](evidence/phase_4j107a/)含汇总、J106 基线、重排数值检查、120 项回归日志及带逐文件哈希的原始产物。

本阶段验收针对固定 25 个差异案例，**未将 880+24 或 440+1 当作重新测得的整个 982 面板覆盖率**。完整面板及化学扩展在 J107B 另测。等价表示重试解决的是后端处理问题，不保证所有分型顺序都产生相同参数；四组独立数值对照也不能替代所有化学环境的来源验证。本阶段未新增长时间 MD 或物性可信度结论。
