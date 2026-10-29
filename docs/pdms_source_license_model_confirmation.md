# PDMS 作者仓库：许可范围与模型版本的公开证据核对

本次核查确认了 `param/OPLS` 的技术分支。没有找到将论文 CC-BY 4.0 明确扩展至独立 GitHub 文件的授权；这表示授权范围尚未确认，不能推导为所有本地研究使用均被禁止。作者没有收到本次核查发出的邮件或消息。

## 许可范围

论文 DOI 为 [10.1021/acs.jpcb.4c08471](https://doi.org/10.1021/acs.jpcb.4c08471)。本地 PDF 首页与 [Crossref 的出版方登记记录](https://api.crossref.org/works/10.1021/acs.jpcb.4c08471) 一致：正式发表版本（`content-version=vor`）为 CC-BY 4.0，起始日期 2025-02-03。这是正文许可的直接证据。

补充材料有独立 DOI [10.1021/acs.jpcb.4c08471.s001](https://doi.org/10.1021/acs.jpcb.4c08471.s001)，Crossref 明确登记它是论文的组成部分，但未提供单独 `license` 字段。缺少字段本身不能证明没有许可。ACS 正文、条款页及 ACS Figshare 的相关公开请求返回 HTTP 403；旧版 author-guidelines/licensing URL 则跳到新首页。因此本次没有取得足够证据确认 SI 的完整许可范围。

作者仓库现在位于 [WuGroup-XJTLU/PDMS_benchmark](https://github.com/WuGroup-XJTLU/PDMS_benchmark)，repository ID `902160948`。核对最新提交的完整 tree（42 个条目，未截断），没有 LICENSE、COPYING、COPYRIGHT、NOTICE 文件；GitHub repository `license=null`，许可证接口返回 404。没有 tag 或 release。README 指出这些是论文模拟输入，并说明可以参考或修改参数以运行模拟，但没有写明整套文件的 CC-BY 声明、再分发权或第三方参数的权利范围。

所以可确认作者公开这些文件是为了模拟复现；**尚不能将整套独立 GitHub 输入作为已获 CC-BY 授权的参数包集成和分发**。本地研究运行、复制/修改整个文件及公开再分发应分别判断。GitHub 的公开查看/fork 权限也不等于全面的开源再分发许可，参见 [GitHub 官方许可说明](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/licensing-a-repository)。

## 精确版本与模型分支

固定仓库最新提交为 `0053fda8449797e7e2b0790426cfc52dcd023e3e`（2025-01-26）；OPLS 目录最后变更来自初始提交 `8e39aeb13835a3810f9848aa8d5d157c32cb8b6c`（2024-12-12）。2025-01-26 的提交修改的是 HUANG 数据，不表示 OPLS 有新版本。

实际 include 顺序仍为 init → read_data → settings。[固定版本 init](https://github.com/WuGroup-XJTLU/PDMS_benchmark/blob/0053fda8449797e7e2b0790426cfc52dcd023e3e/param/OPLS/opls_system.in.init) 的有效指令如下：

| 项目 | 仓库 `param/OPLS` 实际输入 | 论文主要结果 / 表 S4 |
| --- | --- | --- |
| 电荷 | 与图 7 DFT/NPA 列一致；已有 v3 审计逐一核对全部 32448 个原子 | Qeq 列；正文明确主要 Results 使用 Qeq |
| LJ 函数 | 第 8 行 `lj/charmm/coul/long 9.0 11.0` | `lj/cut/coul/long 11 11` |
| 混合 | 第 12 行 `pair_modify mix geometric` | 不能从表 S4 单独确认完整指令 |
| 排除和 1–4 | 第 13 行 `special_bonds lj/coul 0.0 0.0 0.5` | 不能从表 S4 单独确认完整指令 |
| 长程 | 第 14 行 `kspace_style pppm 0.0001` | 需要该分支的完整输入核对 |

第 7 行的截断 LJ 指令是注释，不能作为运行配置。第 9–11 行提到 OPLSAA/M (2015) 的文字也是注释，不能据此确认当前参数就是完整的 OPLS-AA/M 2015。

数据文件首行明确标识 `LAMMPS data file via write_data, version 22 Aug 2018, timestep = 2100000`。这确认了该数据的写出版本；它不能证明论文全部生产模拟都使用相同版本。未找到原始 Moltemplate `.lt` 定义、精确生成版本或正式的模型 release 身份。

本次还核对固定提交下五份通用 generation/equilibration/production/viscosity/thermal-conductivity 输入。它们引用 `system.in.init`、`system.in.settings`、`system.in.charges`；完整 tree 只有 COMPASS 提交了 `.charges` 文件，没有对应 OPLS Qeq 的闭合输入。不能从这些通用模板推断一套完整的 OPLS Qeq 运行版本，也不能借用 COMPASS 电荷。

## 对 ISLAND 的含义

可准确标识的是“固定作者提交中的 DFT/NPA + CHARMM 切换 LJ 外部 PDMS 输入分支”。不能将其标为论文 Qeq 主结果的完整复现，也不能称为原生 OPLS-AA 完成。

仍待作者或出版方明确的是：独立仓库文件的使用/修改/再分发许可、第三方参数来源范围、正式 Qeq 分支的完整输入，以及各分支的生成工具版本。已有[询问草稿](pdms_author_permission_request.md)，尚未发送。没有以公开页面访问失败为由宣称资产违法或禁止所有本地运行。

本次没有修改参数、代码、准备阶段门槛或历史证据。查询地址、HTTP 状态、元数据原文和精确输入引用保存在 [receipt.json](evidence/pdms_license_model_confirmation_v1/receipt.json)。保持 `production_validated=False`、`simulation_readiness="not_established"`。
