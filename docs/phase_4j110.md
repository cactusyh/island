# J110：有界自动平衡诊断与恢复

recipe 可通过顶层 `equilibration` 在固定 NVT/NPT 之后、生产阶段之前追加 NPT 段。控制器检查累计热力学和结构时间序列；达到连续通过次数时停止，否则最多执行预先声明的段数。

```json
"equilibration": {
  "max_segments": 10,
  "steps_per_segment": 20000,
  "consecutive_passes": 2,
  "production_on_exhaustion": false
}
```

示例见 `examples/recipes/j110_adaptive_pe.json`。段长须与 `md.sample_every` 对齐。`policy` 可显式提供现有 `EquilibrationPolicy` 参数，展开后的默认门槛和用户指定值随 recipe 完整保存；本次验证未放宽默认 512 个保留样本、ESS 50 和漂移/分块判据。

检测温度、势能、密度、平均键长和有限链平均平方回转半径，使用固定 burn-in、ACF、ESS、分块独立性及漂移。压力使用采样/独立性、绝对漂移和目标置信范围；不对可能接近零的压力均值使用相对漂移门槛。周期键路径展开检测绕行；交联网络以键长、组分数量和绕行秩检查结构，不强行为无限网络定义有限分子回转半径。

## 输出与恢复

- `status=budget_exhausted`、`equilibrated=false`：预算用尽且未通过，默认没有生产阶段。
- 显式 `production_on_exhaustion=true` 可运行探索性生产，仍保留 `equilibrated=false` 和停止原因。
- 通过时 `equilibrated="diagnostics_passed"`；证据等级是 `single_run_stationarity`，`physical_validation=false`。通过短期平稳性诊断并不证明长链松弛完成、尺寸收敛或物性准确。
- 未配置 `equilibration` 的旧配方继续固定时长流程，`equilibrated="unmeasured"`。
- 每段保存真实 LAMMPS restart、热力学采样、轨迹、诊断及父检查点哈希；事件、参数、可执行程序、源库、样本或策略改变会触发相应校验失败。
- `--stop-after equilibration` 可停在首段检查点；用同一配方、目录再次运行恢复。调整预算或诊断策略应使用新目录，避免改变已经启动的实验条件。

## 验证与科学边界

最终 J107A–J110 联合回归 **154 passed**；随后针对有效策略默认值固化、恢复和旧构造器接口兼容性再验证 **12 passed**。包括人为平稳/漂移/温度偏离数据的统计检查、连续通过与样本下限、中断续跑、损坏拒绝、未收敛默认跳过生产以及显式探索性生产。

真实 MD 对照为三类体系各四组：基准、换种子、增大尺寸、延长诊断段，合计 **12/12** 流程通过验收：

| 体系 | 力场 | 对照结果 |
|---|---|---|
| PE | OPLS polymer-v2 | 40 或 80 个诊断样本，未收敛 |
| PA6 | PCFF carboxamide-v1 | 40 或 80 个诊断样本，未收敛 |
| DGEBA/PACM 反应体系 | GAFF2 / Gasteiger | 形成 2–6 条交联键，40 或 80 个诊断样本，未收敛 |

三个材料组都严格保留默认统计门槛，所以小预算结束时均返回 `budget_exhausted`。这些数据证明控制器能识别证据不足并可靠停止；**没有得到任何已平衡材料或可靠物性结果**。网络案例绕行秩为 0，尚未验证贯通凝胶的收敛。

```bash
PYTHONPATH=.:src python scripts/validate_j110_equilibration.py work/j110-validation
python -m island.workflows.recipe examples/recipes/j110_adaptive_pe.json work/j110 \
  --lammps "$ISLAND_ASSIGNMENT_LAMMPS" --packmol "$ISLAND_PACKMOL"
```

复现输入、原始日志/轨迹/重启文件、逐组判定与校验清单见 `docs/evidence/phase_4j110/`。后续科学验收需要更长时长、独立重复与尺寸比较，不能用本阶段通过的流程测试代替。
