# J107 复现证据

对应实现 `dd5866b`，解释与复现命令见 [阶段报告](../../phase_4j107.md)。

在仓库根目录运行（仅需 Python 标准库，不会重新进行参数分配或 MD）：

```bash
python docs/evidence/phase_4j107/verify.py
```

校验全部文件及压缩包内每个成员的 SHA256，重新计算覆盖率汇总，检查任务唯一性、输入与来源、实际复用数量及科学检查结果。

| 文件 | 含义 |
|---|---|
| `inputs.json` | 3,996 个任务、固定 982 输入文件哈希、参数库/EMC 来源、120 秒时限 |
| `results.jsonl.gz` / `summary.json` | 最终逐例结果及统计；失败和超时仍在分母内 |
| `final-panel-artifacts.tar.xz` | 最终每例输入、进度、结果、后端日志及运行清单 |
| `reused-native-records.json` | 从本次初始运行保留的 1,316 个旧入口结果及来源；所有完整图结果重新执行 |
| `repeat-only-cases.json` | 旧入口成功、完整图失败的全部主面板案例，保留完整诊断 |
| `ring-auto-probe.json` / `ring-auto-probe-script.txt` | PI736 强制自动环扫描仍失败的独立探针，不计入主面板 |
| `challenge-outcomes.json` | 34 个补充输入在 PCFF/OPLS 上的结果和失败类别 |
| `exploratory-incomplete-results.jsonl.gz` | 修复前探索运行的 2,702 个结果，不能当作完整面板统计 |
| `exploratory-interrupted-study.json` / `exploratory-incomplete-artifacts.tar.xz` | 探索运行中止原因及原始记录，与最终结果分开保留 |
| `numerical-and-md.json` / `numerical-md-artifacts.tar.xz` | 最终版本的来源、有限差分、重排及三种短程周期 MD 结果和原始产物 |
| `extended-source-checks.json` | 五个硫/卤素源模型的类型/电荷匹配，以及 DMSO 缺扭转项的失败对照 |
| `cli-preflight.json` / `cli-probe.json` | CLI 预检和完整图分配示例 |
| `*-tests.log` / `*-regressions.log` | 各次测试原始结果，单独文本仅去除行尾空白；测试组有重叠 |
| `environment.json` | Python/RDKit/NumPy 版本、外部程序哈希与基准工作线程配置 |
| `*-manifest.json` / `manifest.json` | 压缩包成员及整个证据集合的校验清单 |

`final-code-tests.log` 是最终实现上的 90 项回归。其他测试日志记录开发过程中的检查，时间边界见阶段报告；不将各日志通过数相加。`numerical-and-md.json` 的 11 个顶层组均在最终实现上执行。

最终 JSONL 统一刷新失败类别；每例压缩包保留当时后端原始诊断。数值/MD 压缩包包含 LAMMPS 生成的输入、能量、轨迹、restart 等文件。路径是当次运行记录，复现时需按本机安装路径调整命令。覆盖率成功不等于精度验证；本次每种 MD 只有 300 fs，不代表平衡。
