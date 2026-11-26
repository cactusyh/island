# J111 复现证据

- `results.jsonl.gz` / `summary.json`：最终 2,032 项完整图记录，包含失败。PCFF 920/982、OPLS 449/982；挑战输入分别 29/34、25/34。
- `initial-panel.tar.gz`：首次完整 2,032 项；其中 PCFF 结果用于最终报告，OPLS 是修复环状胺规则前的版本。
- `opls-final-panel.tar.gz`：最终 OPLS profile 的全部 1,016 项重新运行。
- `initial-inputs.json` / `opls-final-inputs.json`：各自的输入、源库、EMC 程序及脚本哈希。
- `combine.py`：选择上述完整 PCFF / OPLS 面板并比较历史 J107B；不替换超时或构象失败。
- `comparison.json`：19 项历史失败变成功，2 项成功变失败，17 项成功结果类型计数变化。
- `profile-controls.*`：全部 21 项结果变化的当前代码旧 profile 对照，并复核新 profile 的失败项；19 项仍为旧失败/新成功，两个失败新旧均存在。
- `type-controls.json` / `type-controls.log`：17 项类型计数变化的当前代码旧 profile 对照。4 项甲酰端基类型差异来自 profile，其余 13 项当前新旧一致。
- `geometry-controls.*`：PI250 在 RDKit 2026.03.6 下新旧 profile 都成功；原面板使用 RDKit 2024.09.6，仍保留其失败。
- `numerical-md-checks.json` / `numerical-md.tar.gz`：6 个数值模型和 5 个 150 fs 周期配方全部通过。关闭库仑查表仅用于有限差分检查，不改变生产 MD 默认设置。
- `development-validation.tar.gz`：有限差分超差、加严 Ewald 后仍超差、关闭查表的对照、未解决的芳香酰胺/亚砜，以及环状胺修复检查。
- `regression-initial.log`（152 passed）、`profile-final-tests.log`（36 passed）：相关回归与最终来源模型、接点、负例测试。

原始第三方源库按哈希引用；运行归档保留生成的系数、数据、轨迹、重启、输入和日志。`environment.json` 区分批量/数值运行与 pytest 环境；历史 manifest 未显式绑定 RDKit。不能把历史计数净差直接称为单一分型修改的贡献。

`manifest.json` 绑定数值验证实现快照、最终控制脚本和全部文件/归档成员哈希。归档证明运行与来源，不证明材料平衡或实验物性精度。
