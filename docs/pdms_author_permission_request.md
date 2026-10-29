# 待发送的作者授权与模型确认请求

这是一份供用户审核使用的草稿，尚未发送给作者。不能把本文件当作许可或执行授权。

**Subject:** Permission and model clarification for the PDMS OPLS input files accompanying DOI 10.1021/acs.jpcb.4c08471

Dear Prof. Wu and Dr. Xiang,

We are developing an explicitly external PDMS model adapter in ISLAND. Your article, *Force-Field Benchmark for Polydimethylsiloxane: Density, Heat Capacity, Isothermal Compressibility, Viscosity and Thermal Conductivity*, is marked CC-BY 4.0. We have also inspected the paper-cited PDMS_benchmark repository, now at https://github.com/WuGroup-XJTLU/PDMS_benchmark, commit `0053fda8449797e7e2b0790426cfc52dcd023e3e`. We found no explicit repository license.

Could you please confirm the following?

1. May we use and adapt the OPLS input files for local research simulations? Separately, may we redistribute an attributed, derived parameter bundle and bounded graph-typing implementation? Please specify the license or permitted scope, including any restrictions on third-party Moltemplate parameters. We will not presume that the article's license licenses separate repository assets or third-party content.
2. Does the article's CC-BY grant cover the numerical tables in the accompanying `jp4c08471_si_001.pdf`, or is there a separate grant?
3. The active repository input uses DFT/NPA charges, `lj/charmm/coul/long 9.0 11.0`, geometric sigma/epsilon mixing, `special_bonds lj/coul 0 0 0.5` and `pppm 0.0001`. Table S4 instead uses Qeq charges and `lj/cut/coul/long 11 11`. Which exact inputs and engine/Moltemplate versions correspond to each published branch? Please provide the complete Qeq branch if available.
4. Figure 7 defines E fragments `–OSi(CH3)3` and `–Si(CH3)3`, and M fragments `–Si(CH3)2O–`. The repository assigns exactly one E oxygen to each trimethylsilyl-capped chain. May this oriented fragment rule be applied to other explicitly capped chain lengths? Are HMDS and short DP2/DP4 oligomers within the intended domain? Original typing/build templates would help establish this scope.

For identification, the inspected `param/OPLS/opls_system.data` SHA-256 is `8b03e6e8f3c51525f011d9d64e02d69d81faab78cd17ba66676b208d1c4b9c24`; the active init SHA-256 is `4ae3754ef4f428f60eff9f1a33f80a040ec17b45e41ac1e48bba1be3d3ee2592`.

We will preserve attribution and original charge values, distinguish your external extension from native OPLS-AA, and avoid claims of experimental agreement or production qualification without independent validation.

Thank you.
