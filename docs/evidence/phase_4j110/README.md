# J110 证据

`checks.json` / `real-runs.tar.gz` 保存 12 组独立配方和完整运行文件：OPLS PE、PCFF PA6、GAFF2 DGEBA/PACM，分别做种子、尺寸和时长对照。全部严格报告预算耗尽、未收敛；交联模型绕行秩为 0。

`run-source/` 保存产生上述科学运行的实际脚本/模块快照。此后 J108 额外收紧了自由基拒绝、显式模板选择和事件记录校验，J110 将有效统计默认值完整写入 recipe；没有改变这些模型的数值参数或 MD 协议。快照避免把先前运行误记为用后续源文件生成。

`island-j110-finaltests.log` 是 41 项控制器/混合物/网络/recipe 集成测试。`integration.log` 保留较早的 153 项联合回归。`release-integration.log` 为包含随机序列 MD 的最终 154 项联合回归；`policy-tests.log` 保留策略持久化后的 11 项检查；`compatibility-tests.log` / `policy-checkpoint-runs.tar.gz` 是最终 12 项策略与接口兼容性检查，包括故意损坏样本以确认拒绝恢复。

诊断通过仅为单次轨迹统计平稳性，既不等同于物理平衡，也不证明物性精度。本目录没有已平衡的材料结果。`manifest.json` 校验最终实现、运行源快照、输出文件和每个压缩包成员。
