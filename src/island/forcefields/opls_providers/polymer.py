"""Explicit finite-oligomer requests using the unchanged J82 graph builder."""

from island.forcefields.assignment_core.schema import (
    PolymerSpec,
    data,
    identity,
    require,
)
from island.forcefields.assignment_core.topology import molecule, normalize
from island.forcefields.bundle_adapters.schema import keys, text

SETTINGS = (
    "psmiles",
    "repeat_count",
    "end_groups",
    "protonation_state",
    "charge_state",
    "stereochemistry",
    "tacticity",
    "stereo_seed",
    "intended_export_engine",
)


def finite_graph(settings):
    keys(settings, SETTINGS, "explicit_finite_chain_settings_required")
    for name in ("psmiles", "protonation_state", "stereochemistry", "tacticity"):
        text(settings[name], "explicit_setting_required:" + name)
    require(
        isinstance(settings["end_groups"], list)
        and len(settings["end_groups"]) == 2
        and all(isinstance(s, str) and s.strip() for s in settings["end_groups"]),
        "two_explicit_caps_required",
    )
    require(
        settings["protonation_state"] == "as_encoded",
        "protonation_must_be_encoded_in_graph",
    )
    require(type(settings["charge_state"]) is int, "explicit_integer_charge_required")
    require(type(settings["stereo_seed"]) is int, "explicit_stereo_seed_required")
    require(
        settings["tacticity"]
        in ("not_applicable", "isotactic", "syndiotactic", "atactic"),
        "explicit_tacticity_required",
    )
    require(
        settings["stereochemistry"] in ("achiral", "explicit", "tacticity_generated"),
        "explicit_stereochemistry_required",
    )
    require(
        settings["intended_export_engine"] in ("gromacs", "lammps", "openmm"),
        "unsupported_export_engine",
    )
    graph = normalize(
        PolymerSpec(
            (("repeat", settings["psmiles"]),),
            settings["repeat_count"],
            end_groups=tuple(settings["end_groups"]),
            tacticity=None
            if settings["tacticity"] == "not_applicable"
            else settings["tacticity"],
            stereo_seed=settings["stereo_seed"],
        )
    )
    require(not graph.blockers, "finite_graph_blocked:" + ",".join(graph.blockers))
    from rdkit import Chem

    mol = molecule(graph)
    require(len(Chem.GetMolFrags(mol)) == 1, "one_finite_chain_required")
    require(
        sum(a.GetFormalCharge() for a in mol.GetAtoms()) == settings["charge_state"],
        "formal_charge_state_mismatch",
    )
    require(
        not any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()),
        "radical_scope_unsupported",
    )
    centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if settings["stereochemistry"] == "achiral":
        require(
            not centers and settings["tacticity"] == "not_applicable",
            "achiral_claim_disagrees_with_final_graph",
        )
    else:
        require(
            centers and all(label != "?" for _, label in centers),
            "unassigned_stereochemistry",
        )
        require(
            (settings["stereochemistry"] == "tacticity_generated")
            == (settings["tacticity"] != "not_applicable"),
            "stereochemistry_tacticity_mismatch",
        )
    return graph


def graph_metadata(settings):
    graph = finite_graph(settings)
    return {
        "settings": settings,
        "settings_identity": identity(settings),
        "graph_identity": graph.identity,
        "atom_order": list(graph.atom_order),
        "graph_sha256": identity(data(graph)),
        "architecture": graph.architecture,
        "scope": "this_finite_capped_oligomer_only",
    }


def benchmarks():
    definitions = {
        "PE": ("*CC*", ["*C", "*C"], False),
        "PEO_PEG": ("*CCO*", ["*O", "*C"], False),
        "PS": ("*C[C@H](c1ccccc1)*", ["*C", "*C"], True),
        "PLA": ("*O[C@H](C)C(=O)*", ["*C", "*O"], True),
        "PAN": ("*C[C@H](C#N)*", ["*C", "*C"], True),
    }
    return {
        name: {
            "psmiles": psmiles,
            "repeat_count": 2,
            "end_groups": caps,
            "protonation_state": "as_encoded",
            "charge_state": 0,
            "stereochemistry": "tacticity_generated" if chiral else "achiral",
            "tacticity": "isotactic" if chiral else "not_applicable",
            "stereo_seed": 2026,
            "intended_export_engine": "gromacs",
        }
        for name, (psmiles, caps, chiral) in definitions.items()
    }
