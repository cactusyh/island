"""Execute the predeclared J3 matrix; incomplete source coverage never exits zero.

No automatic seed/source substitutions. Native records are retained separately
from independently checked labels/charges and later numerical/workflow receipts.
"""

import argparse
import json
import time
from pathlib import Path

from validate_pcff_expanded import CASES, independent_charges, labels_for
from validate_pcff_fallbacks import EXTRA

from island.chemistry import from_smiles
from island.forcefields.pcff import (
    assign_automatic_pcff_charges,
    assign_pcff_parameters,
    assign_pcff_source_types,
    define_pcff_model,
    load_pcff_source,
    special_pair_policy,
    type_pcff_atoms,
)
from island.forcefields.pcff.domains import PROFILE_NAME
from island.forcefields.pcff.fallbacks import DOMAIN_POLICY
from island.workflows.bundle import system_data
from island.workflows.storage import encode, json_bytes, publish


def expected_labels(system, case):
    expected = case["expected_non_generated_labels"]
    if expected is None:
        return None
    ids = sorted(system.topology.sites)
    if case["name"] in ("heavy_water", "hydrogen", "hydrogen_sulfide"):
        return dict(zip(ids, expected, strict=True))
    heavy = [i for i in ids if system.topology.sites[i].element != "H"]
    labels = dict(zip(heavy, expected, strict=True))
    hlabels = {
        "s": "hs",
        "n": "hn",
        "npc": "hn",
        "nn": "hn2",
        "n=": "hn",
        "n=1": "hn",
        "n=2": "hn",
        "p": "hp",
        "o": "ho2",
        "oh": "ho2",
        "sio": "hsi",
        "sh": "hs",
        "n4": "hn",
    }
    for i in set(ids) - set(heavy):
        parents = [
            b.site2 if b.site1 == i else b.site1
            for b in system.topology.bonds.values()
            if i in b.key
        ]
        if len(parents) != 1:
            raise ValueError("Expected explicit terminal H")
        labels[i] = hlabels.get(labels[parents[0]], "hc")
    return labels


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--declaration",
        type=Path,
        default=Path("docs/evidence/phase_4j3_declaration.json"),
    )
    args = parser.parse_args()
    root = args.output
    root.mkdir(parents=True, exist_ok=False)
    declaration = json.loads(args.declaration.read_text())
    source = load_pcff_source(args.source)
    assert source.identity == declaration["source"]
    publish(root / "declaration.json", json_bytes(declaration))
    cases = [
        {
            "name": n,
            "smiles": sm,
            "domain": "original_j2",
            "expected_non_generated_labels": expected,
        }
        for n, sm, expected, _ in CASES
    ]
    cases += [
        {
            "name": n,
            "smiles": sm,
            "domain": "original_j2",
            "expected_non_generated_labels": None,
        }
        for n, sm in EXTRA
    ]
    cases += declaration["new_cases"]
    report = {
        "source": source.identity,
        "cases": [],
        "production_validated": False,
        "simulation_readiness": "not_established",
        "bounded_gate": False,
        "full_source_complete": False,
    }
    for case in cases:
        name = case["name"]
        out = {
            **case,
            "typed": False,
            "charges_complete": False,
            "model_complete": False,
        }
        directory = root / name
        directory.mkdir()
        start = time.monotonic()
        try:
            system = from_smiles(case["smiles"], random_seed=declaration["seed"])
            publish(directory / "system.json", json_bytes(encode(system_data(system))))
            t = type_pcff_atoms(system, source, profile=PROFILE_NAME)
            publish(directory / "typing.json", t.json_text.encode())
            out.update(
                typed=t.complete,
                typing_identity=t.identity,
                typing_diagnostics=t.payload["diagnostics"],
            )
            if not t.complete:
                continue
            out["types"] = sorted(set(t.assignments.values()))
            explicit = assign_pcff_source_types(
                system,
                source,
                t.assignments,
                provenance="checked J3 automatic assignment, not an independent typing oracle",
                profile=PROFILE_NAME,
            )
            assert explicit.assignments == t.assignments
            if (
                case["domain"] == "original_j2"
                and case["expected_non_generated_labels"] is not None
            ):
                expected = labels_for(system, case["expected_non_generated_labels"])
            else:
                expected = expected_labels(system, case)
            out["expected_labels_match"] = (
                expected == t.assignments if expected is not None else None
            )
            if expected is not None and expected != t.assignments:
                out["expected_label_discrepancy"] = {
                    "expected": expected,
                    "actual": t.assignments,
                }
            q = assign_automatic_pcff_charges(
                system, t, resolution_policy=DOMAIN_POLICY
            )
            publish(directory / "charges.json", q.json_text.encode())
            out.update(
                charges_complete=q.complete,
                charge_identity=q.identity,
                charge_diagnostics=q.payload["native_charge_record"]["diagnostics"],
            )
            if not q.complete:
                continue
            try:
                reference, rows = independent_charges(source.raw, system, t.assignments)
                error = max(abs(reference[i] - q.charges[i]) for i in reference)
                assert error <= 1e-12
                out["independent_charge"] = {
                    "max_error": error,
                    "rows": rows,
                    "scope": "independent raw direct/ordinary oracle; same validated explicit labels",
                }
            except ValueError as error:
                out["independent_charge_unavailable"] = str(error)
            a = assign_pcff_parameters(system, t, q, resolution_policy=DOMAIN_POLICY)
            publish(directory / "assignment.json", a.json_text.encode())
            m = define_pcff_model(
                a, special_pairs=special_pair_policy(lj=(0, 0, 1), coulomb=(0, 0, 1))
            )
            publish(directory / "model.json", m.json_text.encode())
            out.update(
                model_complete=m.payload["model_definition_complete"],
                model_identity=m.identity,
                diagnostics=m.payload["diagnostics"],
                coverage=a.payload["coverage"],
            )
        except Exception as error:  # noqa: BLE001 -- preserve every declared failed case
            out["error"] = type(error).__name__ + ": " + str(error)
        finally:
            out["wall_seconds"] = time.monotonic() - start
            publish(directory / "outcome.json", json_bytes(out))
            report["cases"].append(out)
            print(
                name,
                out["typed"],
                out["charges_complete"],
                out["model_complete"],
                out.get("error", ""),
                flush=True,
            )
    publish(root / "report.json", json_bytes(report))
    return (
        1  # Full orchestration/numerical gates are evaluated by the separate summary.
    )


if __name__ == "__main__":
    raise SystemExit(main())
