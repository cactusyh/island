"""Pinned local engine definitions, transcribed offline; no runtime discovery."""

PAPER_PBC = {
    "source_file": "jp4c08471.pdf",
    "source_medium": "PDF",
    "source_sha256": "04ecee67243157ada68ef96e6f67d5bf31743e8b207efb9b3e4263a823dbc1fc",
    "doi": "10.1021/acs.jpcb.4c08471",
    "page": 3,
    "line": 3,
    "raw_text": "All simulations are performedusing theLAMMPS (Large-\n"
    "ScaleAtomic/MolecularMassivelyParallelSimulator)package\n"
    "with periodic boundary conditions applied in the x, y, and z\n"
    "directions.27 ",
    "scope": "paper Methods: all simulations; boundary only, no Qeq validation transfer",
}

ENGINE_EVIDENCE = {
    "doc/src/angle_harmonic.rst": {
        "authority": "pinned_engine_definition; not author authority",
        "citations": [{"line": 33, "raw_text": "   E = K (\\theta - \\theta_0)^2"}],
        "sha256": "4c03ae7bc4409ce827035cee2dc14b632a34f3af662c3aafd1edbe41d40e2ec0",
    },
    "doc/src/bond_harmonic.rst": {
        "authority": "pinned_engine_definition; not author authority",
        "citations": [{"line": 33, "raw_text": "   E = K (r - r_0)^2"}],
        "sha256": "a092eec0207f0079c005c92c0c0c930c0b20d1c2015c031b22876d65451d7ef6",
    },
    "doc/src/compute_pressure.rst": {
        "authority": "pinned_engine_definition; not author authority",
        "citations": [
            {
                "line": 41,
                "raw_text": "   P = \\frac{N k_B "
                "T}{V} + \\frac{1}{V "
                "d}\\sum_{i=1}^{N'} "
                "\\vec r_i \\cdot \\vec "
                "f_i",
            },
            {
                "line": 79,
                "raw_text": "compute temperature or "
                "ke and/or the virial.  "
                "The *virial* keyword "
                "means",
            },
            {
                "line": 66,
                "raw_text": "components related to "
                "the kinetic energy "
                "tensor and the second "
                "term",
            },
        ],
        "sha256": "6278e7f80ad6e34de4b1f4912171bdec3214f4f868b4567a11aaafdd243daedf",
    },
    "doc/src/dihedral_opls.rst": {
        "authority": "pinned_engine_definition; not author authority",
        "citations": [
            {
                "line": 35,
                "raw_text": "   E = & \\frac{1}{2} K_1 "
                "[1 + \\cos(\\phi)] + "
                "\\frac{1}{2} K_2 [1 - "
                "\\cos(2 \\phi)] + \\\\",
            },
            {
                "line": 38,
                "raw_text": "Note that the usual 1/2 "
                "factor is not included in "
                "the K values.",
            },
        ],
        "sha256": "221ebe20ead28e92658245511c644e9cf255f3058a4e799158df0b628ded6d6a",
    },
    "doc/src/kspace_style.rst": {
        "authority": "pinned_engine_definition; not author authority",
        "citations": [
            {
                "line": 473,
                "raw_text": "Note that the long-range "
                "electrostatic solvers in "
                "LAMMPS assume conducting",
            },
            {
                "line": 252,
                "raw_text": "   must be used if energy "
                "and/or pressure are "
                "quantities of interest,",
            },
        ],
        "sha256": "96058c166ccbe8547d13a80e6bf94a72cc28276b65334187388a761aeced5d93",
    },
    "doc/src/pair_charmm.rst": {
        "authority": "pinned_engine_definition; not author authority",
        "citations": [
            {
                "line": 197,
                "raw_text": "Note that :math:`\\sigma` "
                "is defined in the LJ "
                "formula as the "
                "zero-crossing",
            },
            {
                "line": 179,
                "raw_text": "specified for these styles; "
                "if only 2 arguments are "
                "used in the",
            },
            {
                "line": 121,
                "raw_text": "LJ and Coulombic "
                "interactions with an energy "
                "switching function which",
            },
        ],
        "sha256": "2273f0fc5c03fdccbbc54e42f132d92526ca8d237866c4948d9a935f7e918109",
    },
    "doc/src/units.rst": {
        "authority": "pinned_engine_definition; not author authority",
        "citations": [
            {"line": 100, "raw_text": "For style *real*, these are the units:"},
            {"line": 103, "raw_text": "* distance = Angstroms"},
            {"line": 105, "raw_text": "* energy = kcal/mol"},
            {"line": 107, "raw_text": "* force = (kcal/mol)/Angstrom"},
            {
                "line": 112,
                "raw_text": "* charge = multiple of electron charge (1.0 is a proton)",
            },
        ],
        "sha256": "442a188db85bacf28bbbe47541382c09aa6ba6f2a8271bc2aadbf35d16987834",
    },
    "src/KSPACE/pair_lj_charmm_coul_long.cpp": {
        "authority": "pinned_engine_definition; not author authority",
        "citations": [
            {
                "line": 167,
                "raw_text": "            switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *",
            },
            {
                "line": 169,
                "raw_text": "            switch2 = 12.0*rsq * (cut_ljsq-rsq) *",
            },
        ],
        "sha256": "8efd83f9c3569b21598177adbf253c03a62fad0767d2d30147d0fe8f686b42e5",
    },
    "src/force.cpp": {
        "authority": "pinned_engine_definition; not author authority",
        "citations": [{"line": 64, "raw_text": "  dielectric = 1.0;"}],
        "sha256": "412ec56b1ed4212ef6f91a0c106d24de8bd005a5e72620ce92b02f8abf5a0041",
    },
    "src/pair.cpp": {
        "authority": "pinned_engine_definition; not author authority",
        "citations": [{"line": 104, "raw_text": "  ncoultablebits = 12;"}],
        "sha256": "4bb8edec0e5f98971d50fd7ca069ffa3201a934348f5aec005b180d8fa775004",
    },
    "src/update.cpp": {
        "authority": "pinned_engine_definition; not author authority",
        "citations": [
            {
                "line": 177,
                "raw_text": "    force->qqr2e = 332.06371;    // "
                "see also force->qqr2d_lammps_real",
            },
            {"line": 176, "raw_text": "    force->nktv2p = 68568.415;"},
        ],
        "sha256": "362978558ca39ccb12d78ba6a59b1cbac339af545732c2bfc4af75ba01319f59",
    },
}

CONVERTER_PINS = {
    "forcefields/oplsaa/pdms_author.py": "1ed30d4abd73b62b59d299b06adb9b2d32c28d70d2f52cd17c36d6e4a0ba9a4e",
    "forcefields/pdms_opls.py": "b36f24bbc7cc943db4a976e6db1ec23da49ae99be829926e782edc9b22cec848",
    "pdms.py": "24db78b6bbb438c80d1d3497ff865305db103f1ca7a8bc661e29efb7b71fcb5e",
}

ENGINE_EVIDENCE.update(
    {
        "doc/src/pair_modify.rst": {
            "authority": "pinned_engine_definition; not author authority",
            "citations": [
                {
                    "line": 23,
                    "raw_text": "       *mix* value = "
                    "*geometric* or *arithmetic* "
                    "or *sixthpower*",
                }
            ],
            "sha256": "dcd3f84bf95784da917a87e16d154ed5f978b7491948db7c4b835fa9037f13cc",
        },
        "doc/src/special_bonds.rst": {
            "authority": "pinned_engine_definition; not author authority",
            "citations": [
                {
                    "line": 115,
                    "raw_text": "factor on 1-2 atom pairs, "
                    "which are pairs of atoms "
                    "directly bonded to",
                },
                {
                    "line": 116,
                    "raw_text": "each other.  The second "
                    "coefficient is the "
                    "weighting factor on 1-3",
                },
                {
                    "line": 119,
                    "raw_text": "on 1-4 atom pairs which "
                    "are those separated by 3 "
                    "bonds (e.g. the first",
                },
            ],
            "sha256": "4a166e664256cdaef62ffd2c76eec0020542bb9a83e1a995d8edc01a1960aba6",
        },
    }
)
