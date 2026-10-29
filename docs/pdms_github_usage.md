# 使用 GitHub 分支中的 PDMS 研究入口

代码位于 `codex/pdms-public-research-workflow` 分支。安装与首次执行：

```bash
git clone --branch codex/pdms-public-research-workflow https://github.com/cactusyh/island.git
cd island
python3.11 -m venv .venv
.venv/bin/python -m pip install -e '.[chemistry,evaluation,minimization,dev]' 'openmm==8.6.1'

# 用户自行提供已获准使用的作者 ZIP；文件路径必须显式传入。
OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
.venv/bin/python scripts/pdms_md.py \
  --archive /absolute/path/to/PDMS_benchmark-main.zip \
  --archive-sha256 8b0e8b74303b39b1b0650be3a15b8dd201b20027465f95ffc3f2989da123fea4 \
  --dp 4 --box-angstrom 40 --temperature 300 --timestep-fs 0.5 \
  --steps 1000 --report-interval 100 --output /tmp/island-pdms-first-run
```

输出目录必须尚不存在。`--prepare-only` 可先检查最终图、来源及参数映射，并导出 LAMMPS 文件。实际 MD 使用 OpenMM Reference；LAMMPS 输入导出不会自动启动 LAMMPS MD。完整 Python API、续跑方式和实际数值检查见 [PDMS 用户指南](pdms_user_workflow.md)。

本分支包含 ISLAND 代码、测试、审计记录和生成的验证结果，没有上传本地论文、补充材料、作者参数 ZIP、数值引擎或虚拟环境。运行此入口需要显式 ZIP 和独立 SHA-256，不需要整个 `island_ffref` 目录，也不需要 CAR/MDF。没有提供来源时返回分类阻止结果，不会启动 MD。

完整历史回执校验 `scripts/validate_pdms_public.py` 还核对本地论文和早期审计快照；仅克隆仓库无法完成这些本地资产的核验。历史清单与回执保留原始字节，缺少资产不能视为校验通过。发布目录中的 `docs/evidence/pdms_publication_v1/` 单独记录提交内容、未发布资产及独立检出检查，不替代原来的科学验收。

仅支持已检查的三甲基硅封端、中性、显式氢线性 PDMS。参数来自 Xiang et al. (2025), DOI `10.1021/acs.jpcb.4c08471` 作者输入的 DFT/NPA 分支。profile 为 `island_external_opls_pdms_xiang2025_dft_local_v1`，保持外部 OPLS 兼容身份。论文与独立代码资产的许可范围分别记录；没有将论文许可自动扩大到作者 ZIP 的再分发。

当前结果是短链局部研究 NVT：不表示体相平衡或材料性质验证；没有 NPT、压力、长时间生产模拟或通用 LAMMPS 能量一致性结论。`production_validated=False`、`simulation_readiness="not_established"` 保持不变。
