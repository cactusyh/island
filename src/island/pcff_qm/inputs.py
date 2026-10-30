"""Five explicit engine input adapters. This module never starts a QM engine."""

import json
from pathlib import Path

from island.pcff_completion import _verified, inventory, require, signed, write
from island.pcff_completion.qm import ADAPTERS

from .geometry import geometries, graph_record
from .targets import PROTOCOL, validate_target


def _body(elements, xyz):
    return "\n".join(
        e + " " + " ".join(format(v, ".17g") for v in x)
        for e, x in zip(elements, xyz, strict=True)
    )


def deck(engine, body, *, hessian=False):
    """Proposed PBE0/def2-SVP protocol; CP2K basis libraries remain explicit assets.

    Templates request gradients (forces are minus gradients) and central Hessians.
    Native outputs need an externally pinned importer; no output is synthesized.
    """
    require(engine in ADAPTERS, "unsupported_qm_adapter")
    if engine == "psi4":
        return (
            "# Input preparation only; PBE0, no dispersion; no optimization\n"
            f"molecule mol {{\n0 1\n{body}\nunits angstrom\nno_reorient\nno_com\n}}\n"
            "set {\n basis def2-svp\n scf_type pk\n e_convergence 1e-10\n d_convergence 1e-8\n maxiter 200\n}\n"
            "gradient('pbe0')\n" + ("hessian('pbe0')\n" if hessian else "")
        )
    if engine == "pyscf":
        return (
            "# Input only. Keep stdout and package version; no optimization/dispersion.\n"
            "from pyscf import gto, dft, hessian\n"
            f"mol = gto.M(atom={body!r}, unit='Angstrom', basis='def2-svp', charge=0, spin=0, symmetry=False)\n"
            "mf = dft.RKS(mol)\nmf.xc = 'pbe0'\nmf.conv_tol = 1e-10\n"
            "mf.conv_tol_grad = 1e-8\nmf.max_cycle = 200\n"
            "mf.kernel()\nassert mf.converged\nprint('ENERGY_HARTREE', mf.e_tot)\n"
            "print('GRADIENT_HARTREE_PER_BOHR', mf.nuc_grad_method().kernel().tolist())\n"
            + (
                "print('HESSIAN_HARTREE_PER_BOHR2', hessian.RKS(mf).kernel().tolist())\n"
                if hessian
                else ""
            )
        )
    if engine == "orca":
        return (
            "# No dispersion or optimization. Input only; preserve .out and .engrad/.hess.\n"
            "! PBE0 def2-SVP EnGrad" + (" NumFreq" if hessian else "") + "\n"
            "%scf\n MaxIter 200\n TolE 1e-10\n TolRMSP 1e-8\nend\n"
            f"* xyz 0 1\n{body}\n*\n"
        )
    if engine == "gaussian":
        return (
            "#p PBE1PBE/Def2SVP SCF=(Conver=10,MaxCycle=200) NoSymm "
            + ("Freq" if hessian else "Force")
            + "\n\nJ57 input only; isolated singlet, no dispersion or optimization\n\n0 1\n"
            + body
            + "\n\n"
        )
    # PBE0 all-electron GAPW. The basis file and Poisson implementation must be
    # supplied and pinned before execution; these placeholders are visible assets.
    return (
        "# TEMPLATE: supply/pin DEF2_BASIS_FILE and CP2K_KIND_DECLARATIONS before execution.\n"
        "# PBE0/def2-SVP, no dispersion. Libraries may not be silently substituted.\n"
        "&GLOBAL\n PROJECT j57\n RUN_TYPE "
        + ("VIBRATIONAL_ANALYSIS" if hessian else "ENERGY_FORCE")
        + "\n&END GLOBAL\n"
        "&FORCE_EVAL\n METHOD Quickstep\n &DFT\n"
        "  BASIS_SET_FILE_NAME <DEF2_BASIS_FILE>\n  CHARGE 0\n  MULTIPLICITY 1\n"
        "  &QS\n   METHOD GAPW\n  &END QS\n"
        "  &POISSON\n   PERIODIC NONE\n   POISSON_SOLVER WAVELET\n  &END POISSON\n"
        "  &SCF\n   EPS_SCF 1.0E-10\n   MAX_SCF 200\n  &END SCF\n"
        "  &XC\n   &XC_FUNCTIONAL PBE\n    &PBE\n     SCALE_X 0.75\n     SCALE_C 1.0\n    &END PBE\n   &END XC_FUNCTIONAL\n"
        "   &HF\n    FRACTION 0.25\n    &INTERACTION_POTENTIAL\n     POTENTIAL_TYPE COULOMB\n    &END INTERACTION_POTENTIAL\n   &END HF\n  &END XC\n &END DFT\n"
        " &SUBSYS\n  &CELL\n   ABC 40 40 40\n   PERIODIC NONE\n  &END CELL\n"
        f"  &COORD\n{body}\n  &END COORD\n"
        "  <CP2K_KIND_DECLARATIONS: exact def2-SVP all-electron BASIS_SET and POTENTIAL ALL per element>\n"
        " &END SUBSYS\n &PRINT\n  &FORCES ON\n  &END FORCES\n &END PRINT\n&END FORCE_EVAL\n"
    )


def emit(system, assessment, metadata, output, *, engines=tuple(ADAPTERS)):
    p = _verified(system, assessment)
    require(
        engines and len(set(engines)) == len(engines) and set(engines) <= set(ADAPTERS),
        "unsupported_qm_adapter",
    )
    out = Path(output)
    require(not out.exists(), "output_directory_already_exists")
    samples = geometries(system, assessment, metadata["scan"])
    out.mkdir(parents=True)
    elements = [system.topology.sites[i].element for i in p["atom_order"]]
    labels = [e + str(i) for e, i in zip(elements, p["atom_order"], strict=True)]
    xyz0 = [system.coordinates.get(i).tolist() for i in p["atom_order"]]
    (out / "input-geometry.xyz").write_text(
        f"{len(elements)}\nexact J54 input; not optimized\n{_body(elements, xyz0)}\n"
    )
    write(out / "metadata.json", metadata)
    write(out / "graph.json", graph_record(system, p))
    write(
        out / "atom-labels.json",
        {"atom_order": p["atom_order"], "labels": labels, "elements": elements},
    )
    write(out / "geometries.json", samples)
    for sample in samples:
        body = _body(elements, sample["coordinates"])
        sid = sample["id"]
        (out / (sid + ".xyz")).write_text(
            f"{len(elements)}\n{p['graph_identity']} {sid} input_only\n{body}\n"
        )
        central = all(v == 0 for v in sample["offsets"])
        for engine in engines:
            ext = "py" if engine == "pyscf" else "inp"
            (out / f"{sid}.{engine}.{ext}").write_text(
                deck(engine, body, hessian=central)
            )
    result = signed(
        {
            "schema": "island.pcff-qm-input-bundle.v2",
            "metadata_identity": metadata["identity"],
            "graph_identity": p["graph_identity"],
            "atom_order": p["atom_order"],
            "qm_protocol": PROTOCOL,
            "sample_count": len(samples),
            "files": inventory(out),
            "engines": list(engines),
            "engine_assets": {e: "not_probed_in_reproducible_bundle" for e in engines},
            "cp2k_template_assets": [
                "exact all-electron def2-SVP basis library",
                "element KIND declarations",
                "cell/vacuum and grid convergence study",
            ],
            "adapter_validation": "input syntax prepared; actual engine compatibility not executed or claimed",
            "expected_outputs": [s["expected_output_filename"] for s in samples],
            "output_schema": "island.pcff-qm-reference.v2",
            "native_output_units": {
                "energy": "hartree",
                "gradient": "hartree/bohr",
                "force_sign": "force = minus gradient; explicit importer conversion required",
            },
            "QM_executed": False,
            "reference_data_status": "not_supplied",
            "parameter_completeness": "incomplete",
            "execution_authority": False,
        }
    )
    write(out / "manifest.json", result)
    return result


def generate(system, assessment, target, output, *, engines=tuple(ADAPTERS)):
    return emit(
        system,
        assessment,
        validate_target(system, assessment, target),
        output,
        engines=engines,
    )


def validate_bundle(root, expected_manifest_sha256):
    from island.pcff_completion import digest, pin

    root = Path(root)
    require(
        pin(root / "manifest.json") == expected_manifest_sha256,
        "bundle_manifest_hash_mutation",
    )
    m = json.loads((root / "manifest.json").read_text())
    require(
        m.get("identity") == digest({k: v for k, v in m.items() if k != "identity"}),
        "bundle_identity_mutation",
    )
    files = inventory(root)
    files.pop("manifest.json")
    require(files == m["files"], "bundle_inventory_mutation")
    return m
