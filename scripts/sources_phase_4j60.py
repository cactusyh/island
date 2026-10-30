import lzma
import subprocess
from pathlib import Path

from island.pcff_completion import pin, signed, write

E = Path("docs/evidence/phase_4j60")
U = Path("../island-validation/phase4h4-upstream/lammps")
files = {
    "FRC": "tools/msi2lmp/frc_files/pcff.frc",
    "converter": "tools/msi2lmp/src/GetParameters.c",
    "topology": "tools/msi2lmp/src/MakeLists.c",
    "dihedral_doc": "doc/src/dihedral_class2.rst",
    "improper_doc": "doc/src/improper_class2.rst",
    "pair_doc": "doc/src/pair_class2.rst",
    "special_doc": "doc/src/special_bonds.rst",
    "units_doc": "doc/src/units.rst",
    "dielectric_doc": "doc/src/dielectric.rst",
    "boundary_doc": "doc/src/boundary.rst",
    "dihedral_code": "src/CLASS2/dihedral_class2.cpp",
    "pair_code": "src/CLASS2/pair_lj_class2_coul_cut.cpp",
    "real_units_code": "src/update.cpp",
    "mixing_doc": "doc/src/pair_modify.rst",
}
records = {}
for k, f in files.items():
    p = U / f
    dest = E / "sources" / f.replace("/", "__")
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.with_suffix(dest.suffix + ".xz").write_bytes(lzma.compress(p.read_bytes()))
    records[k] = {
        "file": f,
        "sha256": pin(p),
        "source_commit": subprocess.check_output(
            ["git", "-C", str(U), "rev-parse", "HEAD"], text=True
        ).strip(),
        "citation": "LAMMPS pinned distribution and official documentation; Sun JPCB 102 7338 (1998) cited by engine docs; FRC internal references retained",
        "doi": "10.1021/jp980939v"
        if k in ("dihedral_doc", "improper_doc", "pair_doc")
        else None,
        "doi_scope": "engine documentation citation, not identity of the FRC distribution",
    }


def ev(k, lo, hi):
    p = U / files[k]
    return {
        **records[k],
        "line_range": [lo, hi],
        "raw_text": "\n".join(p.read_text().splitlines()[lo - 1 : hi]),
    }


def convention(name, value, unit, evidence, authority, derivation):
    return signed(
        {
            "name": name,
            "interpreted_value": value,
            "units": unit,
            "evidence": evidence,
            "authority": authority,
            "implementation": "island_research_class2_sih_disiloxane_j60_v1",
            "reconstruction": derivation,
        }
    )


items = [
    convention(
        "bond_bond13_applicability",
        "Retain all six selected pilot BB13 paths as one type-symmetric fitted unknown; do not copy the converter cp-only omission rule",
        "kcal/mol/angstrom^2",
        [
            ev("FRC", 3616, 3626),
            ev("converter", 506, 529),
            ev("converter", 1222, 1240),
            ev("dihedral_doc", 43, 57),
        ],
        "declared_research_extension_from_source_equation_and_engine_form",
        "The pinned converter performs BB13 lookup only if a cp type appears. No cp occurs here. This does not establish a physical zero for sio/osi/hsi. J54-J60 explicitly target the missing Si-H coefficient, and the selected LAMMPS class2 Hamiltonian supports BB13 for every declared dihedral type. All six sites remain explicit unknown contributions.",
    ),
    convention(
        "ordering_and_phi",
        "I-J-K-L path; cis=0; signed Phi from atan2(axis cross projected I-J dot projected L-K, projection dot); full path reversal invariant",
        "radian",
        [ev("dihedral_code", 147, 236)],
        "documented_implementation_convention",
        "Projection dot/(norms) equals (c0+costh12*costh23)/(sin(theta1)*sin(theta2)); signed triple product agrees with LAMMPS branch.",
    ),
    convention(
        "main_torsion",
        "Use FRC plus-cos with declared LAMMPS geometry; export K=V, executable phase=(source phase+pi) modulo 2pi",
        "kcal/mol,radian",
        [ev("FRC", 2655, 2664), ev("converter", 342, 350), ev("dihedral_doc", 39, 57)],
        "source_equation_plus_declared_mapping_research_extension",
        "1-cos(n*phi-(delta+pi))=1+cos(n*phi-delta), including constant. Converter copies phases unchanged and does not establish this mapping. Alternative source geometry is not proven; this choice is research-owned.",
    ),
    convention(
        "angle_angle_torsion",
        "LAMMPS cosine AAT selected for research variant; literal FRC linear-Phi retained and not converted",
        "kcal/mol/radian^2",
        [ev("FRC", 5183, 5189), ev("dihedral_doc", 43, 52), ev("converter", 489, 504)],
        "declared_research_variant_conflicts_with_literal_FRC",
        "For nonzero angle product, linear Phi gains 2pi*K*dt1*dt2 per turn; cosine returns unchanged. No nonzero constant rescaling or phase shift makes them globally equal. New AAT coefficient is fitted, never copied from the linear-Phi section.",
    ),
    convention(
        "torsion_torsion",
        "Outside the selected LAMMPS class2 sum; no numeric zero row",
        "model_selection",
        [ev("FRC", 5519, 5526), ev("dihedral_doc", 43, 57)],
        "declared_research_model_scope",
        "An empty FRC section does not establish native absence. The exact chosen finite sum contains no five-site TT. All nine graph paths explicitly marked model-inapplicable.",
    ),
    convention(
        "Wilson",
        "degree-three centers only under pinned topology converter; none in selected graph; degree-four angle-angle terms retained",
        "topology",
        [ev("topology", 409, 438), ev("improper_doc", 28, 57)],
        "converter_defined_research_topology",
        "Degrees: Si4,O2,H1. No degree3; eight neighbor triples at degree4 supply AA, not invented Wilson zero coefficients.",
    ),
    convention(
        "special_pairs",
        "special_bonds lj/coul 0 0 1 angle no dihedral no; full pairs beyond graph distance3",
        "dimensionless",
        [ev("special_doc", 114, 143)],
        "declared_implementation_convention",
        "Graph shortest paths bind weights: 8 direct and 13 two-bond pairs excluded; 6 three-bond and 9 four-bond pairs included. No per-pair fit weights.",
    ),
    convention(
        "LJ96",
        "2*epsilon*(rmin/r)^9-3*epsilon*(rmin/r)^6; sixth-power mixing",
        "kcal/mol,angstrom",
        [ev("FRC", 3254, 3264), ev("pair_doc", 79, 94), ev("pair_doc", 129, 140)],
        "source_equation_and_backend_defined",
        "rmin_ij=((ri^6+rj^6)/2)^(1/6); epsilon_ij=2*sqrt(ei*ej)*ri^3*rj^3/(ri^6+rj^6). Source converted epsilon kJ divided by4.184; rmin retained, not 12-6 sigma.",
    ),
    convention(
        "Coulomb",
        "332.06371*qi*qj/(dielectric*r); dielectric1; source charge vector unchanged",
        "kcal/mol,angstrom,e",
        [ev("real_units_code", 145, 163), ev("dielectric_doc", 10, 45)],
        "backend_defined_and_explicit_research_choice",
        "LAMMPS real prefactor pinned; dielectric1 explicitly declared, not author-confirmed. No charge normalization or fitting.",
    ),
    convention(
        "isolation",
        "boundary f f f; all pairs direct within100A cutoff; no switching, shifting, tail, periodic images, KSPACE or pressure claim; containing box[-60,60]^3",
        "angstrom",
        [ev("pair_doc", 43, 58), ev("boundary_doc", 35, 65)],
        "declared_isolated_molecule_implementation_policy",
        "All pilot pairs are within100A; direct finite molecule gives full included pair sums. Cutoff crossing rejected. Box has no energetic volume role and pressure not assessed.",
    ),
]
write(
    E / "hamiltonian-definition.json",
    signed(
        {
            "internal_identity": "island_research_class2_sih_disiloxane_j60_v1",
            "registered_profile": False,
            "native_PCFF_claim": False,
            "conventions": items,
            "sources": records,
            "author_confirmations": {
                "path": "docs/evidence/phase_4j57/audit/author-attested-confirmations.json",
                "sha256": pin(
                    "docs/evidence/phase_4j57/audit/author-attested-confirmations.json"
                ),
                "kind": "user_attested_author_confirmation",
                "original_author_artifact_status": "not_yet_supplied",
            },
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
    ),
)
