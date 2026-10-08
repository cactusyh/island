"""Retained real-source DP3 assignments packed through J21 for J22 verification."""

from pathlib import Path

from island.forcefields import PreparedForceFieldSources, load_prepared_forcefield
from island.graph import final_graph
from island.packing import PeriodicPackingConfig, pack_multichain_periodic
from island.periodic import PeriodicForceFieldConfig, prepare_periodic_forcefield

ROOT = Path(__file__).resolve().parents[2]
RETAINED = ROOT / "island-validation/phase4i2"
FRC = ROOT / "lammps/lammps/tools/msi2lmp/frc_files/pcff.frc"
XML = ROOT / "island-validation/foyer-4g1-upstream/foyer/forcefields/xml/oplsaa.xml"


def retained_fixture(family, *, method="pme", tolerance=1e-10):
    sources = (
        PreparedForceFieldSources(pcff_frc=FRC)
        if family == "pcff"
        else PreparedForceFieldSources(opls_xml=XML)
        if family == "oplsaa"
        else None
    )
    loaded = load_prepared_forcefield(RETAINED / family / "relocated", sources=sources)
    native = loaded.prepared.native_result
    packed = pack_multichain_periodic(
        [("polyethylene-DP3", loaded.system, final_graph(loaded.system))],
        config=PeriodicPackingConfig(box_lengths=(30.0, 30.0, 30.0), seed=2026),
        provenance="J22 single retained polymer rigid placement",
        evidence=("docs/evidence/phase_4i2.json",),
    )
    if family == "pcff":
        source = native.payload["source"]
        mix = "sixthpower_9_6"
        scales = (
            native.payload["special_pairs"]["lj"][2],
            native.payload["special_pairs"]["coulomb"][2],
        )
    elif family == "oplsaa":
        source = {
            **loaded.prepared.source.identity,
            "sha256": loaded.prepared.source.identity["xml_sha256"],
        }
        mix = "geometric_12_6"
        scales = (0.5, 0.5)
    else:
        r = native.imported_result
        source = {
            "family": family.upper(),
            "sha256": r.source_sha256,
            "source": r.source,
        }
        mix = "lorentz_berthelot_12_6"
        scales = (r.nonbonded_policy.lj_scale_14, r.nonbonded_policy.coulomb_scale_14)
    config = PeriodicForceFieldConfig(
        family="OPLS-AA" if family == "oplsaa" else family.upper(),
        source=source,
        box_lengths=(30.0, 30.0, 30.0),
        nonbonded_cutoff=10.0,
        lj_mixing_rule=mix,
        one_four_scaling=scales,
        electrostatics_method=method,
        pme_tolerance=tolerance,
        provenance="J22 retained real-source DP3 periodic verification",
        evidence=("docs/evidence/phase_4i2.json", packed.plan.identity),
    )
    result = prepare_periodic_forcefield(
        packed.system, packed.graph, loaded.prepared, config
    )
    return result, sources
