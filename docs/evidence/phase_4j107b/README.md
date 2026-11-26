# J107B 复现证据

- `source-checks.json`：8 个来源 residue 的逐原子类型/电荷比较；5 个开发探针。开发探针不计作留出模型。
- `dynamics-checks.json` / `dynamics.tar.gz`：4 个独立留出模型的原子置换/有限差分检查；4 个 300 fs 周期 MD 配方与原始轨迹、系数、日志。
- `gaff-native-*`：固定 982 输入 × GAFF/GAFF2，各自对照独立调用的原生 RadonPy 类型、参数和 Gasteiger 电荷；这是无坐标完整图测试，不计作 1,964 个 MD 模型。
- `gaff-architectures.json`：每个力场 12 个实际构象模型，包含真正的嵌段/支链接点和端基、立体化学输入，逐项比较原生拓扑/系数导出。
- `paired-*`：982 × 2 力场 × 2 入口及 34 × 2 挑战图的完整分母，共 3,996 项。
- `comparison.json`：与保留不变的 J107 结果逐案例比较；OPLS 新增成功案例以当前代码、旧 v1 来源另行复测，区分代码/时限和 profile 的作用。
- `tests.log`：125 项相关回归通过。
- `profile.json`：来源和安装后文件哈希；未改变数值参数文件。
- `superseded-draft.tar.gz`：接续开发前的初步报告和小面板，仅供审计。其“留出模型”等标签已更正，不能用该旧报告代替本目录正式结果。

复现脚本位于 `scripts/benchmark_final_graphs.py`、`scripts/benchmark_gaff_native_graphs.py`、`scripts/validate_gaff_graph_coverage.py`、`scripts/validate_j107b_coverage.py`、`scripts/validate_j107b_dynamics.py`。`manifest.json` 绑定实现、文件与归档成员校验和。
