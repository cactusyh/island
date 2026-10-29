# 用户提供的 PDMS_benchmark ZIP 核查

已检查 `island_ffref/pcff/PDMS_benchmark-main.zip`。压缩包完整，32 个文件全部与此前审计的作者仓库提交 `0053fda8449797e7e2b0790426cfc52dcd023e3e` 相同。它提供了本地、可核验的作者输入副本，没有增加新的参数、许可或构建模板。本次没有运行作者脚本或 PDMS MD，没有修改历史证据或原生力场。

对应论文是 Xiang et al., *Force-Field Benchmark for Polydimethylsiloxane: Density, Heat Capacity, Isothermal Compressibility, Viscosity and Thermal Conductivity*, *J. Phys. Chem. B* **129** (2025), 1864–1873，DOI `10.1021/acs.jpcb.4c08471`。

| 核对项 | 结果 |
| --- | --- |
| ZIP 大小 | 28,850,437 字节 |
| ZIP SHA-256 | `8b0e8b74303b39b1b0650be3a15b8dd201b20027465f95ffc3f2989da123fea4` |
| ZIP comment | `0053fda8449797e7e2b0790426cfc52dcd023e3e` |
| 条目 | 43 个，含 32 个文件；解压后共 110,154,891 字节 |
| 完整性 | 全部文件 CRC 校验通过 |
| 固定版本比较 | 全部 32 个文件的长度和 Git blob SHA-1 一致；另记录每个文件的 SHA-256 |
| 目录 | OPLS、COMPASS、Dreiding、HUANG、UA 参数与通用模拟输入 |
| 许可证文件 | 没有 LICENSE、COPYING、COPYRIGHT 或 NOTICE |
| 原始构建模板/代码 | 没有 Moltemplate `.lt`、Python、C/C++ 或 shell 文件 |
| 独立 `.charges` 文件 | 只有 `param/COMPASS/system.in.charges` |

ZIP comment 本身不能证明版本；这里同时核对了此前固定、带独立 SHA-256 的完整 Git tree 响应。逐文件结果见 [archive-audit.json](evidence/pdms_user_archive_v1/archive-audit.json)。

## OPLS 输入的实际含义

`param/OPLS/opls_system.in` 引用的 init、data、settings 均在同一目录存在。输入包含真实参数和带电荷的原子数据，**并非只有论文说明**。这只是静态依赖检查，不表示已经通过运行和科学验证。

本次将 ZIP 中五份已固定的 README/OPLS 文件复制到临时目录，重新运行独立原文审计：51 处行引用、17 条相关系数和单位转换、32,448 个原子、64 条链全部与旧审计一致。保留原始十进制电荷字面值，总和为 `-2.8800E-13 e`，没有归一化。结果见 [independent-opls-reconstruction.json](evidence/pdms_user_archive_v1/independent-opls-reconstruction.json)。

实际有效指令仍是：

```lammps
units real
bond_style harmonic
angle_style harmonic
dihedral_style opls
improper_style harmonic
pair_style lj/charmm/coul/long 9.0 11.0
pair_modify mix geometric
special_bonds lj/coul 0.0 0.0 0.5
kspace_style pppm 0.0001
```

电荷对应已审计的图 7 DFT/NPA 分支。第 7 行 `lj/cut/coul/long` 被注释；不能把它当作当前运行设置。它与论文主要 Results 的 Qeq/表 S4 分支有区别，不能宣称已经取得完整 Qeq 模型。OPLSAA/M (2015) 注释也不能单独证明完整模型版本。

另外检查了 ZIP 中全部 10 份 `.data`。通用 `input/` 数据使用类型 6、13、39、44，对应电荷约为 C −0.588、H +0.106、O −0.89、Si +1.43/+1.255；其类型与电荷对应 COMPASS 数据分支，与 OPLS 的类型 98、122、868、871 和 DFT/NPA 电荷不同。因此不能将通用输入中的结构或 COMPASS `.charges` 当成 OPLS Qeq 的缺失输入。原始字面值及计数保存在 JSON 中，没有修约后导入。

静态检查还发现 HUANG、UA 的入口引用 `.in.settings`，实际文件名为 `.in.setting`。本次仅记录，不修改作者文件，也不把其他模型的依赖问题当成 OPLS 的问题。

## 对继续接入的影响

这份 ZIP 没有改变此前的结论：可准确追溯的是作者固定提交中的外部 DFT/NPA + CHARMM 切换 LJ 分支。论文正文的 CC-BY 4.0 已有证据；独立仓库的许可范围没有新增声明。README 表达了模拟复现、参考和修改参数的用途意图，但未明确整体再分发范围。没有据此推断本地研究使用一概禁止。

新图的端基/电荷候选映射和长程人工对照仍见 [执行边界审计](opls_pdms_execution_boundary.md)。本次没有新增执行 profile，没有用下载成功代替真实 PDMS 的能量、力、最小化、MD 或续跑验收。原有作者精度对照未通过的结果保持不变。

保持 `production_validated=False`、`simulation_readiness="not_established"`；Ewald/PME、压力、cutoff recovery 和 long-NPT 限制不变。

## 复现

下列命令只读取显式本地资产，不联网、不解压执行脚本，不覆盖已存在的输出文件。退出码 0 表示源副本核对通过，不表示 PDMS 模拟获得执行资格。

```bash
python scripts/audit_pdms_author_archive.py \
  --archive island_ffref/pcff/PDMS_benchmark-main.zip \
  --sha256 8b0e8b74303b39b1b0650be3a15b8dd201b20027465f95ffc3f2989da123fea4 \
  --tree docs/evidence/pdms_license_model_confirmation_v1/github_tree.body \
  --output /tmp/pdms-user-archive-audit.json
```

验证包含：迁移 ZIP/固定 tree 后的独立进程逐字节重建、错误 ZIP 哈希拒绝、参数源内容改变但重新计算 ZIP 哈希后的拒绝、固定 tree 篡改拒绝，以及历史 manifest 校验。既有 3476 个受保护文件和许可核对 manifest 的 21 个文件校验通过。Ruff、format、compileall、diff check 通过；本次没有更改生产代码，没有重跑全套测试，不沿用旧测试数量作为本次结果。验证详情和校验和见 [receipt.json](evidence/pdms_user_archive_v1/receipt.json)。
