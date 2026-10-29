# PDMS 作者 OPLS 输入：单位、末端映射与非键约定审计

经用户授权，已读取论文引用的作者仓库，固定提交并独立重建其 OPLS 输入。原先缺少的单位、LJ 混合、1–4／排除约定，以及作者已有拓扑中的 E/M 电荷分区，现有明确来源。**这些结论适用于本次固定的作者输入：它使用 DFT/NPA 电荷和 CHARMM 切换 LJ，与论文表 S4 的 Qeq／截断 LJ 模型不同。** 因此不能直接宣布论文模型已可执行。

这是新增 v3 数据审计，不是执行 profile。原生 OPLS-AA、PCFF 和先前 v1/v2 文献补充均保持原样；没有准备数值系统或运行 PDMS MD。

## 来源与版本

论文为 Xiang et al., *Force-Field Benchmark for Polydimethylsiloxane: Density, Heat Capacity, Isothermal Compressibility, Viscosity and Thermal Conductivity*, J. Phys. Chem. B **129** (2025), 1864–1873，DOI `10.1021/acs.jpcb.4c08471`。本地论文及 SI 的两个哈希仍为：

| 文件（相对 `island_ffref`） | SHA-256 |
| --- | --- |
| `pcff/jp4c08471.pdf` | `04ecee67243157ada68ef96e6f67d5bf31743e8b207efb9b3e4263a823dbc1fc` |
| `pcff/jp4c08471_si_001.pdf` | `6e86aa64085466e286b5c59878c2288f5f16da1b7a5298d37ba1833a88cd602f` |

正文引用的 `https://github.com/COMO-lab/PDMS_benchmark` 现转向 `https://github.com/WuGroup-XJTLU/PDMS_benchmark`，GitHub repository ID 为 `902160948`。固定提交为 `0053fda8449797e7e2b0790426cfc52dcd023e3e`，日期 `2025-01-26T07:03:32Z`；对应 tree 为 `5b3bcc406226a9f5af9fcaee4abf5a465e9abc26`。该版本的 OPLS 路径历史只列出初始提交 `8e39aeb13835a3810f9848aa8d5d157c32cb8b6c`（2024-12-12）；这不能证明它就是论文最终生产运行版本。

| 下载并审计的文件 | SHA-256 |
| --- | --- |
| `README.md` | `79f2c4c56fe32f44b344b1d0c6c13c628fd1ef9b64941cda80e8d62d8866191a` |
| `param/OPLS/opls_system.in` | `dfd22dac270ad11bf50001573e7d581ea924e287036e34aa85cded67913d6cdd` |
| `param/OPLS/opls_system.in.init` | `4ae3754ef4f428f60eff9f1a33f80a040ec17b45e41ac1e48bba1be3d3ee2592` |
| `param/OPLS/opls_system.in.settings` | `d22a19757398bdf9fb04536f4acf74a7ad3dabce686aaba65d58675bba0016e5` |
| `param/OPLS/opls_system.data` | `8b03e6e8f3c51525f011d9d64e02d69d81faab78cd17ba66676b208d1c4b9c24` |

下载仅使用该仓库的 GitHub API；同时核对大小、Git blob SHA-1 和 SHA-256。作者脚本未执行。原始下载保存在仓库外，未打包重新分发完整参数库。API 的 `license` 字段为 null，完整 tree 未列出许可证文件；用途及再分发许可仍记为 `not_established`。README 将内容说明为论文输入，但不能据此补写一个许可证。

## 四类缺口的核查结果

实际覆盖顺序为 `opls_system.in` 第 5 行 include init → 第 10 行 read_data → 第 14 行 include settings。审计同时对比 data 内的系数和后续 settings 的覆盖值。

| 项目 | 作者实际证据 | 本次结论 |
| --- | --- | --- |
| 单位 | init 第 1 行 `units real` | Å、kcal/mol、fs、e；各势函数的方程另由本地 LAMMPS 文档核实。 |
| LJ 混合 | init 第 12 行 `pair_modify mix geometric` | sigma、epsilon 都使用几何混合；检查 settings 未发现异类 pair override。 |
| 1–4 与排除 | init 第 13 行 `special_bonds lj/coul 0.0 0.0 0.5` | LJ 与 Coulomb 的 1–2／1–3 权重为 0，1–4 为 0.5；检查本次 include 闭包没有额外 `neigh_modify exclude`。 |
| 末端／电荷 | data 的实际原子、电荷和键图 | 每条链有两个 E 端 Si(CH₃)₃，且只有一个 E 氧；与正文图 7 的 DFT/NPA 列逐原子一致。此处验证的是作者已有拓扑，不是任意新图的通用自动类型规则。 |

init 第 7 行的 `#pair_style lj/cut/coul/long 11.0 11.0` 是注释；**实际启用的是第 8 行 `pair_style lj/charmm/coul/long 9.0 11.0`**。不能为了与 SI 对齐而选取注释行。`kspace_style pppm 0.0001` 是作者设置记录，不代表 ISLAND 的 polymer Ewald/PPPM 验证已经完成。

作者 data 的 `Pair Coeffs # zero` 中 `# zero` 是段标题注释；行内包含实际 epsilon/sigma 及相同的 epsilon14/sigma14。该注释不是可据以填零的物理参数。

## 已重建的系数与拓扑

从实际使用的 type ID 出发选取 **17 条**系数：4 条同类 LJ、3 条键、6 条角和 4 条 proper torsion。这些数值均与表 S4 的对应行相同；没有引入 Moltemplate 大库中未使用的其他条目。每条保存 data 行、最终 settings 行、SI 行、原始十进制值、单位换算和证据 checksum。

| 相互作用 | 原始系数（LAMMPS real） |
| --- | --- |
| O–Si、C–Si、C–H 键 | `(374, 1.64)`、`(187, 1.86)`、`(340, 1.09)`；分别为 K 和 r₀ |
| Si–O–Si、C–Si–C、O–Si–C、O–Si–O 角 | `(20,145)`、`(60,110)`、`(60,100)`、`(60,110)`；分别为 K 和 θ₀ |
| H–C–H、H–C–Si 角 | `(33,107.8)`、`(35,109.5)` |
| Si–O–Si–C、O–Si–O–Si、H–C–Si–O 扭转 | 作者与 SI 明确列出 `(0,0,0,0)`；不是补入缺失项。 |
| H–C–Si–C 扭转 | `(0,0,0.18,0)` |
| H/O/Si/C 同类 LJ | epsilon 分别 `0.03/0.14/0.1/0.066`，sigma 分别 `2.5/2.9/4.0/3.5`。 |

对原始图重新枚举角和二面角，核查恰好覆盖 32,448 个原子、64 条分子链、32,384 个键、61,504 个角和 77,568 个二面角，且没有重复／缺失项。此拓扑有 **0 个 improper interaction**，不能外推为其他分子也无需 improper。质量来自作者 data：H 1.008、O 15.999、Si 28.086、C 12.011。

每条链共 507 个原子（含显式 H），两个端 Si 各连接三个甲基，内部 Si 各连接两个甲基和两个氧。原始电荷对应正文图 7 的 DFT 列：

| 元素 | M / e | E / e |
| --- | ---: | ---: |
| O | −1.261 | −1.255 |
| Si | 2.121 | 1.874 |
| C | −1.102 | −1.083 |
| H | 0.224 | 0.2225 |

审计逐项比对作者 `write_data` 的 binary64 序列化字面值，不用容差归类、不编辑电荷。原始序列化十进制数每条链之和为 `−4.50E−15 e`，整个体系为 `−2.8800E−13 e`；两者原样保存，没有归一化。这与原始印刷十进制列的精确中性不矛盾。

单位换算保留势函数方程：

- LAMMPS harmonic 键 `E=K(r−r₀)²`；若映射到 OpenMM 的 `E=k(r−r₀)²/2`，则 `k=836.8 K`，单位 kJ mol⁻¹ nm⁻²；`r₀,nm=0.1 r₀,Å`。
- harmonic 角同理 `k=8.368 K`，单位 kJ mol⁻¹ rad⁻²；θ₀ 保存为精确表达式 `θ₀,deg*pi/180`。
- OPLS 扭转保留 `Σ Kₙ[1+(−1)ⁿ⁺¹cos(nφ)]/2`，Kₙ 乘 4.184 转为 kJ/mol。
- LJ epsilon 乘 4.184，sigma 乘 0.1。

上述为数据换算证据，尚未构建 OpenMM 系统。7 份本地 LAMMPS 方程／命令文档的哈希与行引用也记录在 audit 中，不能把这些本地文档版本视为原始作者引擎版本。

## 还不能执行的原因

论文正文和 SI 表 S4 使用 **Qeq 电荷、`lj/cut/coul/long 11 11`**；作者当前 OPLS 输入使用 **DFT/NPA 电荷、`lj/charmm/coul/long 9 11`**。相同的系数并不使这两套总势能相同。本次没有把作者输入的混合／排除规则默认为论文 Qeq 模型的已验证规则，也没有切换电荷或 pair style 来消除差别。

还缺：模型分支及论文运行版本的确认、明确的资产使用许可、用于新最终图的 Moltemplate 类型规则与版本来源，以及独立引擎能量／力和完整工作流验收。Si–H、Si–OH、Si–O–C、P、带电／自由基等环境均未从该作者拓扑外推。

新增身份为 `island_external_opls_pdms_author_input_audit_v3`，payload identity 为 `9f28383c1bd37522bb2309aaec31d1ad709a0913b7f67ba150f45b30d323651d`。状态为 `scientifically_unqualified`，`profile_identity=None`、`execution_authorized=False`、`acceptance_gate=1`。`require_executable()` 始终给出携带完整记录与缺失依赖的分类异常；原生 `ForceFieldRequest` 不接受该数据审计对象。

## 复现与验收

显式下载是独立的开发命令；运行时和离线审计不会联网。该命令只读取上述固定作者提交的元数据及 5 份文件，目标目录必须不存在：

```bash
export PYTHONPATH=.:src
../island-validation/phase4e2-env/bin/python scripts/opls_pdms_author.py \
  --fetch-author-source-to /tmp/pdms-author-source
```

随后离线生成数据审计：

```bash
../island-validation/phase4e2-env/bin/python scripts/opls_pdms_author.py \
  --sources /tmp/pdms-author-source \
  --engine-docs /inspire/ssd/project/sais-suiren-foundation-model-prediction/public/yh/lammps/lammps/doc/src \
  --references island_ffref \
  --selection docs/evidence/phase_4j45/selection.json \
  --catalog docs/evidence/phase_4j44/source-audit.json \
  --output /tmp/opls-pdms-author-audit --require-executable
```

预期退出码 **1**，分类状态 `scientifically_unqualified`；去掉 `--require-executable` 才是数据导出的成功码 0。下载工具也实际验证过：先核对 commit 对应的 tree，再按 tree SHA 读取；避免把 GitHub 以 commit 查询 tree 时回显的 SHA 误作 tree SHA。

导出记录重建：将上述 `--output ...` 替换为以下两个参数，其余输入保持显式提供：

```text
--reconstruct docs/evidence/opls_pdms_author_v3/record.json
--sha256 33fb0a3f5be078b29532123c9563b47068fd54ca827fc0df0f0160454e1d0e49
```

独立解析器只导入标准库，不调用生产源解析器。传入审计 JSON 的独立哈希后，重建 5 个文件、51 条原始行引用、17 条系数、全部原子电荷和各分子的电荷总和：

```bash
python scripts/audit_opls_pdms_author_sources.py \
  --sources /tmp/pdms-author-source \
  --audit docs/evidence/opls_pdms_author_v3/audit.json \
  --sha256 9c753877262dc679c8ba6f1d9a2a72f21eef66beaa94ed44ca64d2a947cfd2af
```

完整新证据清单、收据和 3,417 个既有文件的字节校验：

```bash
../island-validation/phase4e2-env/bin/python scripts/validate_opls_pdms_author.py \
  docs/evidence/opls_pdms_author_v3 --repo . \
  --sha256 <independent-manifest-sha256> \
  --sources /tmp/pdms-author-source --references island_ffref \
  --engine-docs /inspire/ssd/project/sais-suiren-foundation-model-prediction/public/yh/lammps/lammps/doc/src
```

manifest 校验成功的 0 只表示数据与历史字节验证通过，不能替代执行验收的非零结果。独立哈希应从可信传递的记录取得，不能用被测文件自己报告的哈希代替。

本次最终相关回归为 **251 passed、1 skipped**，其中新增作者源测试 **19 passed**；跳过项是需要单独提供 pinned Foyer XML 的既有原生 OPLS 测试。测试命令和实际结果保存在 [validation.json](evidence/opls_pdms_author_v3/validation.json)。测试显式设置 `ISLAND_PDMS_AUTHOR_SOURCES` 和 `ISLAND_LAMMPS_DOCS`；没有这些输入就跳过作者源测试，不自动发现或下载资产。覆盖篡改后重新签名的电荷／模型／单位／权限记录、源字节改变、明确下载边界和新进程迁移重建。此前 v2 的图映射负对照及 OPLS 相关回归继续保留。

这次小范围审计未重新运行项目 full supported suite；不把历史全套计数当作本次结果。没有 PDMS 能量、力、最小化或 MD 验收。保持 `production_validated=False`、`simulation_readiness="not_established"`；polymer Ewald/PME、压力、cutoff recovery、long NPT 限制不变。
