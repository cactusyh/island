# J109 证据

- `scale-checks.json` / `scale.log`：六组 PE DP10/50/100 构建与 MD、精确配方、种子、链数、耗时、逐链回转半径、累计父进程峰值 RSS（KiB，不包含外部程序独立峰值）。
- `scale-runs.tar.gz`：六组的构象、装箱、映射、原生参数、MD 轨迹、日志和 checkpoint。
- `mixture-runs.tar.gz`：混合链长、嵌段、支链以及独立构象/旧入口回归的原始输出。
- `island-j109-tests.log`：71 项相关测试通过。
- `random-mixture.tar.gz` / `island-j109-distribution-test.log`：随机长度/序列混合物的真实 MD 和恢复测试。
- `packing-checks.json` / `packing-replay.tar.gz` / `packing.log`：原 J105 矩阵的 18 个独立几何输入，全部装箱成功。只计算几何覆盖，不计算为新增力场或 MD 证据。

采样仅用于功能/规模验证，不构成密集熔体或平衡物性结果。实现、文件与归档成员由 `manifest.json` 校验。
