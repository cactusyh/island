"""Authentic-source frozen J17 negative controls and unchanged guanidinium record."""

import argparse
from pathlib import Path

from island.charge_references.records import pack
from island.chemistry import from_smiles
from island.forcefields import (
    ForceFieldRequest,
    PCFFOptions,
    PreparedForceFieldError,
    prepare_forcefield,
)
from island.forcefields.pcff import (
    PCFFAutomaticChargeResult,
    assign_pcff_parameters,
    bind_pcff_types,
    define_pcff_model,
    load_pcff_source,
    provide_pcff_charges,
    special_pair_policy,
)
from island.forcefields.pcff.fallbacks import MSI_POLICY
from island.workflows import storage


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", required=True, type=Path)
    p.add_argument("--historical", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    a = p.parse_args()
    source = load_pcff_source(a.source)
    historical = storage.decode(storage.read_json(a.historical))
    old = historical[MSI_POLICY]["native_charge_record"]
    charge = PCFFAutomaticChargeResult(pack(old), source)
    charge.validate_integrity()
    native = charge.payload["native_charge_record"]
    assert not charge.complete and abs(native["total_charge"] - 0.9999) < 1e-12
    assert native["formal_charge"] == 1
    system = from_smiles("COC", random_seed=2026)
    types = bind_pcff_types(
        system,
        source,
        {
            i: {"C": "c", "H": "h", "O": "oc"}[s.element]
            for i, s in system.topology.sites.items()
        },
        provenance="FRC generic SP3 carbon and ether/carbon-hydrogen descriptions",
        evidence_references=["pcff.frc atom_types c/h/oc"],
    )
    # Explicitly synthetic zero charges isolate interaction completeness. No
    # chemical or scientific acceptance is attributed to this vector.
    charges = provide_pcff_charges(
        system,
        types,
        dict.fromkeys(system.topology.sites, 0.0),
        unit="elementary_charge",
        provenance="software completeness control; zero vector",
        evidence_references=["docs/phase_4j17.md frozen incomplete ether control"],
        component_totals={min(system.topology.sites): 0},
        total_charge=0,
        resolution_policy=MSI_POLICY,
    )
    assignment = assign_pcff_parameters(
        system, types, charges, resolution_policy=MSI_POLICY
    )
    model = define_pcff_model(
        assignment, special_pairs=special_pair_policy(lj=(0, 0, 1), coulomb=(0, 0, 1))
    )
    assert charges.complete and not model.payload["model_definition_complete"]
    try:
        prepare_forcefield(
            system,
            ForceFieldRequest(
                "pcff",
                PCFFOptions(
                    a.source,
                    (0, 0, 1),
                    (0, 0, 1),
                    resolution_policy=MSI_POLICY,
                    typed_graph=types,
                    graph_charges=charges,
                ),
            ),
        )
    except PreparedForceFieldError as error:
        rejection = str(error)
        assert rejection.startswith("Incomplete PCFF model:")
    else:
        raise AssertionError("Incomplete interactions were published")
    result = {
        "passed": True,
        "guanidinium": {
            "historical_charge_identity": charge.identity,
            "total_charge": native["total_charge"],
            "formal_charge": native["formal_charge"],
            "complete": charge.complete,
            "historical_sha256": storage.checksum(a.historical.read_bytes()),
        },
        "incomplete_ether": {
            "types_valid": True,
            "provided_charges_valid": charges.complete,
            "charge_origin": "software_test_zero_vector",
            "scientific_acceptance": False,
            "model_complete": False,
            "diagnostics": model.payload["diagnostics"],
            "rejection": rejection,
        },
    }
    storage.publish(a.output, storage.json_bytes(result))
    print(result)


if __name__ == "__main__":
    main()
