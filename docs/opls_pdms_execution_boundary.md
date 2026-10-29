# PDMS：源授权、最终图映射与作者长程模型的验证边界

这次完成了新最终图的来源绑定映射和长程势函数的独立数值对照。**实际 PDMS 体系尚未获得执行资格，没有运行 PDMS MD。** 资产许可仍有缺口；高精度的人工数值对照不能替代作者 PDMS 原设置的验证。既有 J17–J45、PDMS v1/v2/v3 和外部 OPLS v1 的文件保持不变。

## 授权证据

Xiang et al., *J. Phys. Chem. B* **129** (2025), 1864–1873，DOI `10.1021/acs.jpcb.4c08471`，本地 `island_ffref/pcff/jp4c08471.pdf` 首页明确写有 `This article is licensed under CC-BY 4.0`。本次独立读取 PDF 正文和元数据，并核对图 7 原始图像哈希。论文正文中的使用授权现在有明确证据，需要保留作者、引用、许可和修改说明。

这不自动许可独立的 GitHub 文件。再次查询固定的论文作者仓库后，最新提交仍为 `0053fda8449797e7e2b0790426cfc52dcd023e3e`，API `license=null`，许可证接口返回 404。README 的运行说明没有被改写成开源许可。支持信息本地文件中也未找到单独的明确授权条款，未擅自将正文许可扩展至独立资产。用户允许我们读取与分析，不能替作者授予这些资产的权利。

需要作者确认的具体范围已写成[授权与模型确认请求草稿](pdms_author_permission_request.md)，尚未发送。已有真实授权时，应保留原文、授权人身份、来源、适用文件/提交、允许用途和再分发范围；文件 checksum 本身不能证明授权真实。

## 新图的端基、类型和电荷

新增 [pdms_mapping.py](../src/island/forcefields/opls_external/pdms_mapping.py) 产生不可变的 `island_external_pdms_author_graph_mapping_v4` 候选记录。它不是执行 profile，也没有注册到原生 OPLS 或 PCFF。

图 7 明确定义 E 端片段 `–OSi(CH3)3`、`–Si(CH3)3` 和 M 片段 `–Si(CH3)2O–`。已有 v3 对作者原始 64 条链的重建提供逐原子的独立来源核对；这次将相同的有向片段定义应用到经过化学验证的 HMDS、三甲基硅端基 DP2 和 DP4 新图。对存在两个候选端氧的链，调用者必须显式指定 `e_oxygen_site_id`；PSMILES builder 的 head oxygen 可作为显式输入。没有按最小原子 ID、容器顺序或位置猜测方向。

每条候选记录绑定：最终图、原子顺序、组分/分子身份、端氧选择、精确来源质量、类型 ID、E/M 标签、源哈希、旧审计身份和完整电荷向量。元素、价态、带电/自由基、质量、源字节、端氧、图和记录的篡改会拒绝。未指定方向为 `typing_ambiguous`，Si–H、Si–OH、Si–O–C、P 等没有从甲基 PDMS 外推。

| 新图 | 显式原子数 | E 原子数 | 作者原始序列化电荷的精确和 / e |
| --- | ---: | ---: | ---: |
| HMDS | 27 | 27 | `3/10000000000000000` |
| TMS-capped DP2 | 37 | 27 | `1/5000000000000000` |
| TMS-capped DP4 | 57 | 27 | `0` |

保留作者 DFT/NPA 分支的 binary64 `write_data` 十进制字面值，没有用印刷值替换它们，也没有归一化。DP4 的和本来就是零，不是电荷修复的结果。这条路径没有放宽既有有限体系 external OPLS v1 的精确中性规则。原始 Moltemplate 规则/版本仍未取得；这里明确标为从图 7 重建的规则，不能声称恢复了原脚本或得到作者对短链适用性的确认。

独立脚本 [reconstruct_pdms_graph_mapping.py](../scripts/reconstruct_pdms_graph_mapping.py) 不导入 ISLAND 生产解析器；从源 data、电荷字面值、最终图和显式端氧独立核对三个新图的每个类型与电荷。

## 长程后端验证

新增 [qualify_pdms_longrange_backend.py](../scripts/qualify_pdms_longrange_backend.py) 只运行固定的八粒子人工数学对照，不接受 PDMS 参数，不产生可用执行 profile，也不能绕过准备阶段的阻断。

验证函数为 harmonic bonds/angles、四项 OPLS torsions、几何混合 LJ、1–2/1–3 排除和 1–4 的 LJ/Coulomb 各 0.5。LJ 保留作者使用的 CHARMM 9–11 Å 势能切换：

```text
S(r) = (rc²-r²)² (rc²+2r²-3ri²) / (rc²-ri²)³
ri=9 Å; rc=11 Å
```

该式是 r² 的多项式，不能换成 OpenMM 默认的五次切换函数。没有 LJ tail correction。电静力比较 LAMMPS PPPM、OpenMM PME 和独立直接 Ewald 求和；独立求和保留实空间、倒空间、自能、净电荷背景和特殊对修正。控制体系电荷中性；原作者数据的微小净电荷对实际 PDMS 后端的处理尚未验收。

| 对照 | LAMMPS / OpenMM 精度参数 | 最大分项能量差 kJ/mol | 最大力差 kJ/(mol·Å) | 预设阈值结果 |
| --- | --- | ---: | ---: | --- |
| 作者精度设置 | `1e-4` / `1e-4` | `6.4533e-2` | `4.5329e-2` | `scientifically_unqualified` |
| 中间精度 | `1e-7` / `1e-7` | 约 `1.76e-3` | 约 `9.88e-4` | 未通过 |
| 独立高精度方程对照 | `1e-10` / `1e-7` | 约 `1.33e-6` | 约 `3.28e-6` | 通过 |

作者精度对照保留 LAMMPS 默认 Coulomb table；中间/高精度方程对照显式 `table 0`。较粗对照的失败全部保留，没有加宽阈值，没有将高精度设置冒充作者原输入。PPPM 和 PME 的精度参数也不意味着两者采用相同网格。LAMMPS 和 OpenMM 的 Coulomb 单位常数细微差异在记录中保留。

在两个位移步长 `1e-4 Å`、`5e-5 Å` 上，对全部 24 个力分量进行独立方程有限差分；高精度对照最大误差小于 `5e-6 kJ/(mol·Å)`。这验证了人工对照下的函数与单位，不是实际 PDMS 能量、力、最小化、NVT/NPT、续跑或 RNG 验收。

## 复现

运行时不联网，不发现参数文件。新图审计显式提供全部来源：

```bash
export PYTHONPATH=.:src
PY=../island-validation/phase4g1-env/bin/python
$PY scripts/audit_pdms_execution_boundary.py \
  --sources /tmp/island-pdms-author-fetch-verification \
  --references island_ffref \
  --engine-docs /inspire/ssd/project/sais-suiren-foundation-model-prediction/public/yh/lammps/lammps/doc/src \
  --selection docs/evidence/phase_4j45/selection.json \
  --catalog docs/evidence/phase_4j44/source-audit.json \
  --author-record docs/evidence/opls_pdms_author_v3/record.json \
  --author-record-sha256 33fb0a3f5be078b29532123c9563b47068fd54ca827fc0df0f0160454e1d0e49 \
  --output /tmp/pdms-graph-candidates --require-executable
```

预期为正常分类报告，退出码 **1**。数据审计成功不代表可以运行 MD。

长程数学控制：

```bash
OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 $PY \
  scripts/qualify_pdms_longrange_backend.py \
  --executable ../island-validation/phase4j22/lmp-kspace \
  --executable-sha256 6659313ef24e7a8477008c98d19c9410cd0a04e47cff0a9d91216f1040bf9191 \
  --output /tmp/longrange-equation-controls
```

输出目录必须尚不存在。成功码 0 仅表示固定数学对照完成，并保留作者精度对照未通过的状态；没有物理 PDMS 资格。

测试、独立重建、迁移和历史哈希的命令与结果见新增 [receipt.json](evidence/opls_pdms_execution_boundary_v4/receipt.json)。本次新增测试 **30 passed**，相关回归 **271 passed**；最后对长程脚本细化后的三个相关测试另行重跑并通过。全套 supported suite 本次未重新运行；历史 full-suite 计数不能作为本次结果。Ruff、compileall、五个 pip check 和 diff check 均通过，3476 个既有文件未改变。

新证据的 manifest 需要由独立哈希校验：

```bash
$PY scripts/validate_pdms_execution_boundary.py \
  docs/evidence/opls_pdms_execution_boundary_v4 \
  --repo . --sha256 <independent-manifest-sha256>
```

## 后续执行门槛

仍需确认源资产许可及短链适用范围，并在明确选择的作者 DFT 分支上完成真实 PDMS 数据导出、独立能量/力、最小化、允许的短 NVT、迁移及精确续跑验收。论文 Qeq 主结果模型需取得与其一致的完整输入，不能借用 DFT 分支的不同 pair-style 约定来宣称论文模型复现成功。

保持 `production_validated=False`、`simulation_readiness="not_established"`。polymer Ewald/PME、backend pressure、cutoff recovery、long NPT 的既有限制不变；没有生产、平衡或实验吻合结论。

针对本次修改的测试可用以下显式来源复现：

```bash
ISLAND_PDMS_AUTHOR_SOURCES=/tmp/island-pdms-author-fetch-verification \
ISLAND_LAMMPS_DOCS=/inspire/ssd/project/sais-suiren-foundation-model-prediction/public/yh/lammps/lammps/doc/src \
ISLAND_LAMMPS_EXECUTABLE=/inspire/ssd/project/sais-suiren-foundation-model-prediction/public/yh/island-validation/phase4j22/lmp-kspace \
ISLAND_LAMMPS_EXECUTABLE_SHA256=6659313ef24e7a8477008c98d19c9410cd0a04e47cff0a9d91216f1040bf9191 \
OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
$PY -m pytest tests/test_pdms_execution_boundary.py -q
```
