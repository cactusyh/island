"""Offline transcriptions of pinned engine definitions; never author authority."""

CITATIONS = {
    "coulomb_prefactor": [
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 177,
            "raw_text": "    force->qqr2e = 332.06371;    // see also "
            "force->qqr2d_lammps_real",
            "source_file": "src/update.cpp",
            "source_sha256": "362978558ca39ccb12d78ba6a59b1cbac339af545732c2bfc4af75ba01319f59",
        },
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 178,
            "raw_text": "  qqrd2e = qqr2e / dielectric;",
            "source_file": "src/force.cpp",
            "source_sha256": "412ec56b1ed4212ef6f91a0c106d24de8bd005a5e72620ce92b02f8abf5a0041",
        },
    ],
    "dielectric": [
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 64,
            "raw_text": "  dielectric = 1.0;",
            "source_file": "src/force.cpp",
            "source_sha256": "412ec56b1ed4212ef6f91a0c106d24de8bd005a5e72620ce92b02f8abf5a0041",
        },
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 178,
            "raw_text": "  qqrd2e = qqr2e / dielectric;",
            "source_file": "src/force.cpp",
            "source_sha256": "412ec56b1ed4212ef6f91a0c106d24de8bd005a5e72620ce92b02f8abf5a0041",
        },
    ],
    "erfc_kernel": [
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 140,
            "raw_text": "            expm2 = exp(-grij*grij);",
            "source_file": "src/KSPACE/pair_lj_charmm_coul_long.cpp",
            "source_sha256": "8efd83f9c3569b21598177adbf253c03a62fad0767d2d30147d0fe8f686b42e5",
        },
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 141,
            "raw_text": "            t = 1.0 / (1.0 + EWALD_P*grij);",
            "source_file": "src/KSPACE/pair_lj_charmm_coul_long.cpp",
            "source_sha256": "8efd83f9c3569b21598177adbf253c03a62fad0767d2d30147d0fe8f686b42e5",
        },
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 142,
            "raw_text": "            erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * "
            "expm2;",
            "source_file": "src/KSPACE/pair_lj_charmm_coul_long.cpp",
            "source_sha256": "8efd83f9c3569b21598177adbf253c03a62fad0767d2d30147d0fe8f686b42e5",
        },
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 143,
            "raw_text": "            prefactor = qqrd2e * qtmp*q[j]/r;",
            "source_file": "src/KSPACE/pair_lj_charmm_coul_long.cpp",
            "source_sha256": "8efd83f9c3569b21598177adbf253c03a62fad0767d2d30147d0fe8f686b42e5",
        },
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 144,
            "raw_text": "            forcecoul = prefactor * (erfc + "
            "EWALD_F*grij*expm2);",
            "source_file": "src/KSPACE/pair_lj_charmm_coul_long.cpp",
            "source_sha256": "8efd83f9c3569b21598177adbf253c03a62fad0767d2d30147d0fe8f686b42e5",
        },
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 19,
            "raw_text": "  static constexpr double EWALD_F = 1.12837917;",
            "source_file": "src/ewald_const.h",
            "source_sha256": "3b3abb6ada91a5f39cf507c4c121b0df6a662ed60c4aabfd04b2fe7861abb07c",
        },
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 20,
            "raw_text": "  static constexpr double EWALD_P = 0.3275911;",
            "source_file": "src/ewald_const.h",
            "source_sha256": "3b3abb6ada91a5f39cf507c4c121b0df6a662ed60c4aabfd04b2fe7861abb07c",
        },
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 21,
            "raw_text": "  static constexpr double A1 = 0.254829592;",
            "source_file": "src/ewald_const.h",
            "source_sha256": "3b3abb6ada91a5f39cf507c4c121b0df6a662ed60c4aabfd04b2fe7861abb07c",
        },
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 22,
            "raw_text": "  static constexpr double A2 = -0.284496736;",
            "source_file": "src/ewald_const.h",
            "source_sha256": "3b3abb6ada91a5f39cf507c4c121b0df6a662ed60c4aabfd04b2fe7861abb07c",
        },
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 23,
            "raw_text": "  static constexpr double A3 = 1.421413741;",
            "source_file": "src/ewald_const.h",
            "source_sha256": "3b3abb6ada91a5f39cf507c4c121b0df6a662ed60c4aabfd04b2fe7861abb07c",
        },
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 24,
            "raw_text": "  static constexpr double A4 = -1.453152027;",
            "source_file": "src/ewald_const.h",
            "source_sha256": "3b3abb6ada91a5f39cf507c4c121b0df6a662ed60c4aabfd04b2fe7861abb07c",
        },
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 25,
            "raw_text": "  static constexpr double A5 = 1.061405429;",
            "source_file": "src/ewald_const.h",
            "source_sha256": "3b3abb6ada91a5f39cf507c4c121b0df6a662ed60c4aabfd04b2fe7861abb07c",
        },
    ],
    "lj_tail_correction": [
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 102,
            "raw_text": "  tail_flag = 0;",
            "source_file": "src/pair.cpp",
            "source_sha256": "4bb8edec0e5f98971d50fd7ca069ffa3201a934348f5aec005b180d8fa775004",
        },
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 230,
            "raw_text": "",
            "source_file": "doc/src/pair_charmm.rst",
            "source_sha256": "2273f0fc5c03fdccbbc54e42f132d92526ca8d237866c4948d9a935f7e918109",
        },
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 231,
            "raw_text": "None of the *lj/charmm* or *lj/charmmfsw* pair "
            "styles support the",
            "source_file": "doc/src/pair_charmm.rst",
            "source_sha256": "2273f0fc5c03fdccbbc54e42f132d92526ca8d237866c4948d9a935f7e918109",
        },
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 232,
            "raw_text": ":doc:`pair_modify <pair_modify>` tail option for "
            "adding long-range tail",
            "source_file": "doc/src/pair_charmm.rst",
            "source_sha256": "2273f0fc5c03fdccbbc54e42f132d92526ca8d237866c4948d9a935f7e918109",
        },
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 233,
            "raw_text": "corrections to energy and pressure, since the "
            "Lennard-Jones portion of",
            "source_file": "doc/src/pair_charmm.rst",
            "source_sha256": "2273f0fc5c03fdccbbc54e42f132d92526ca8d237866c4948d9a935f7e918109",
        },
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 234,
            "raw_text": "the pair interaction is smoothed to 0.0 at the cutoff.",
            "source_file": "doc/src/pair_charmm.rst",
            "source_sha256": "2273f0fc5c03fdccbbc54e42f132d92526ca8d237866c4948d9a935f7e918109",
        },
    ],
    "neighbor_list_policy": [
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 19,
            "raw_text": "         N = delay building neighbor lists "
            "until this many steps since last build",
            "source_file": "doc/src/neigh_modify.rst",
            "source_sha256": "9c993eed51d00385852dce31c42755afe83e57a91217d809962fb938cf3457c4",
        },
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 20,
            "raw_text": "       *every* value = M",
            "source_file": "doc/src/neigh_modify.rst",
            "source_sha256": "9c993eed51d00385852dce31c42755afe83e57a91217d809962fb938cf3457c4",
        },
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 21,
            "raw_text": "         M = consider building neighbor lists "
            "every this many steps",
            "source_file": "doc/src/neigh_modify.rst",
            "source_sha256": "9c993eed51d00385852dce31c42755afe83e57a91217d809962fb938cf3457c4",
        },
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 22,
            "raw_text": "       *check* value = *yes* or *no*",
            "source_file": "doc/src/neigh_modify.rst",
            "source_sha256": "9c993eed51d00385852dce31c42755afe83e57a91217d809962fb938cf3457c4",
        },
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 24,
            "raw_text": "         *no* = always build on 1st step where "
            "*every* and *delay* are conditions are "
            "satisfied",
            "source_file": "doc/src/neigh_modify.rst",
            "source_sha256": "9c993eed51d00385852dce31c42755afe83e57a91217d809962fb938cf3457c4",
        },
    ],
    "pair_table_bits": [
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 104,
            "raw_text": "  ncoultablebits = 12;",
            "source_file": "src/pair.cpp",
            "source_sha256": "4bb8edec0e5f98971d50fd7ca069ffa3201a934348f5aec005b180d8fa775004",
        },
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 392,
            "raw_text": "  // linear lookup tables of length N = 2^ncoultablebits",
            "source_file": "src/pair.cpp",
            "source_sha256": "4bb8edec0e5f98971d50fd7ca069ffa3201a934348f5aec005b180d8fa775004",
        },
    ],
    "periodic_image_policy": [
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 203,
            "raw_text": "  if (slabflag == 0 && domain->nonperiodic > 0)",
            "source_file": "src/KSPACE/pppm.cpp",
            "source_sha256": "8281e42208df957c5193731c746fd667c3fa0981e49c6744d625d85fc29b307b",
        },
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 391,
            "raw_text": "  // z dimension for 3d PPPM is zprd since "
            "slab_volfactor = 1.0",
            "source_file": "src/KSPACE/pppm.cpp",
            "source_sha256": "8281e42208df957c5193731c746fd667c3fa0981e49c6744d625d85fc29b307b",
        },
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 72,
            "raw_text": "  slabflag = wireflag = 0;",
            "source_file": "src/kspace.cpp",
            "source_sha256": "f73ab2c2c7bad11ed9211e144894c584fdebac884fa33c5c2d48eed5d6eed07e",
        },
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 74,
            "raw_text": "  slab_volfactor = 1;",
            "source_file": "src/kspace.cpp",
            "source_sha256": "f73ab2c2c7bad11ed9211e144894c584fdebac884fa33c5c2d48eed5d6eed07e",
        },
    ],
    "precision": [
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 358,
            "raw_text": '    mesg += "  using " LMP_FFT_PREC " precision " '
            'LMP_FFT_LIB "\\n";',
            "source_file": "src/KSPACE/pppm.cpp",
            "source_sha256": "8281e42208df957c5193731c746fd667c3fa0981e49c6744d625d85fc29b307b",
        }
    ],
    "pressure_convention": [
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 41,
            "raw_text": "   P = \\frac{N k_B T}{V} + \\frac{1}{V "
            "d}\\sum_{i=1}^{N'} \\vec r_i \\cdot \\vec f_i",
            "source_file": "doc/src/compute_pressure.rst",
            "source_sha256": "6278e7f80ad6e34de4b1f4912171bdec3214f4f868b4567a11aaafdd243daedf",
        },
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 79,
            "raw_text": "compute temperature or ke and/or the virial.  "
            "The *virial* keyword means",
            "source_file": "doc/src/compute_pressure.rst",
            "source_sha256": "6278e7f80ad6e34de4b1f4912171bdec3214f4f868b4567a11aaafdd243daedf",
        },
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 80,
            "raw_text": "include all terms except the kinetic energy *ke*\\ .",
            "source_file": "doc/src/compute_pressure.rst",
            "source_sha256": "6278e7f80ad6e34de4b1f4912171bdec3214f4f868b4567a11aaafdd243daedf",
        },
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 81,
            "raw_text": "",
            "source_file": "doc/src/compute_pressure.rst",
            "source_sha256": "6278e7f80ad6e34de4b1f4912171bdec3214f4f868b4567a11aaafdd243daedf",
        },
    ],
    "virial_convention": [
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 41,
            "raw_text": "   P = \\frac{N k_B T}{V} + \\frac{1}{V "
            "d}\\sum_{i=1}^{N'} \\vec r_i \\cdot \\vec f_i",
            "source_file": "doc/src/compute_pressure.rst",
            "source_sha256": "6278e7f80ad6e34de4b1f4912171bdec3214f4f868b4567a11aaafdd243daedf",
        },
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 60,
            "raw_text": "A symmetric pressure tensor, stored as a "
            "6-element vector, is also",
            "source_file": "doc/src/compute_pressure.rst",
            "source_sha256": "6278e7f80ad6e34de4b1f4912171bdec3214f4f868b4567a11aaafdd243daedf",
        },
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 61,
            "raw_text": "calculated by this compute.  The six components "
            "of the vector are",
            "source_file": "doc/src/compute_pressure.rst",
            "source_sha256": "6278e7f80ad6e34de4b1f4912171bdec3214f4f868b4567a11aaafdd243daedf",
        },
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 62,
            "raw_text": "ordered :math:`xx,` :math:`yy,` :math:`zz,` "
            ":math:`xy,` :math:`xz,`",
            "source_file": "doc/src/compute_pressure.rst",
            "source_sha256": "6278e7f80ad6e34de4b1f4912171bdec3214f4f868b4567a11aaafdd243daedf",
        },
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 63,
            "raw_text": ":math:`yz.` The equation for the :math:`(I,J)` "
            "components (where",
            "source_file": "doc/src/compute_pressure.rst",
            "source_sha256": "6278e7f80ad6e34de4b1f4912171bdec3214f4f868b4567a11aaafdd243daedf",
        },
        {
            "authority": "backend_only",
            "citation": "LAMMPS 10 Sep 2025 development, revision "
            "36f7c81ced8dd456430b13647e48f0c7d7bd8faa",
            "doi": None,
            "line": 176,
            "raw_text": "    force->nktv2p = 68568.415;",
            "source_file": "src/update.cpp",
            "source_sha256": "362978558ca39ccb12d78ba6a59b1cbac339af545732c2bfc4af75ba01319f59",
        },
    ],
}
