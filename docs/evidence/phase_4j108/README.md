# J108 证据

`final-network-runs.tar.gz` 和对应 `island-j108-final-tests.log` 保留最终 10 项测试的原始参数、MD、周期环、续跑及故意损坏测试文件。故意损坏用于验证拒绝行为，不能把它们当作未修改的科学输入。

`release-runs.tar.gz` / `island-j108-release-tests.log` 是加入闭壳层拒绝检查之前的 30 项 J108/旧流程回归。新增拒绝检查不改变这些非自由基模型的参数或动力学；最终 10 项再次覆盖网络正例和拒绝行为。

运行命令见 `docs/phase_4j108.md`；环境设置需要真实 EMC、RadonPy、RDKit、LAMMPS 和 PACKMOL。`manifest.json` 绑定最终实现、文件及每个归档成员的 SHA-256。
