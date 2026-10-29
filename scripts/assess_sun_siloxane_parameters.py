"""Graph-bound Sun siloxane parameter audit; --require-executable stays nonzero."""

import argparse
import json
from pathlib import Path

from island.builders import build_linear_polymer
from island.charge_references.records import pack, unpack
from island.chemistry import from_smiles
from island.forcefields.pdms_literature import INPUT_ERRORS
from island.forcefields.siloxane_literature import (
    SiloxaneParameterAssessment,
    SiloxaneParameterError,
    assess_siloxane_parameters,
    load_siloxane_parameters,
)
from island.workflows.bundle import system_data, system_from


def fixtures():
    for name, smiles in [
        ("disiloxane", "[SiH3]O[SiH3]"),
        ("hexamethyldisiloxane", "C[Si](C)(C)O[Si](C)(C)C"),
        ("si_o_c_negative", "CO[Si](C)(C)C"),
    ]:
        yield name, from_smiles(smiles, random_seed=2026)
    for dp in (1, 2, 4):
        yield (
            f"pdms_dp{dp}",
            build_linear_polymer(
                "[*:1]O[Si](C)(C)[*:2]",
                dp=dp,
                random_seed=2026,
                add_hydrogens=True,
                # Trimethylsilyl head and methyl tail yield Si-O-Si backbones
                # without introducing a methoxy O-C end group. dp counts repeats;
                # the explicit head cap contributes one additional silicon atom.
                head_end_group="[*][Si](C)(C)C",
                tail_end_group="[*]C",
            ),
        )


def write(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def run(catalog, pin, references, output, *, reconstruct=False):
    supplement = load_siloxane_parameters(catalog, pin, references=references)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    summaries = []
    if reconstruct:
        stored_report = json.loads((output / "report.json").read_text())
    for name, system in fixtures():
        if reconstruct:
            rebuilt_system = system_from(
                unpack((output / f"{name}-graph.json").read_text())
            )
            # Also compare against the canonical fixture to reject changed
            # end groups, molecule identities or atom IDs in saved records.
            expected = assess_siloxane_parameters(
                system, supplement, references=references
            )
            stored = SiloxaneParameterAssessment(
                (output / f"{name}-assessment.json").read_text()
            )
            stored.validate_integrity(rebuilt_system, supplement, references=references)
            if stored != expected:
                raise SiloxaneParameterError(
                    "source_incompatible", "fixture_graph_changed"
                )
            result = stored
        else:
            result = assess_siloxane_parameters(
                system, supplement, references=references
            )
            (output / f"{name}-graph.json").write_text(pack(system_data(system)))
            (output / f"{name}-assessment.json").write_text(result.json_text)
        p = result.payload
        summaries.append(
            {
                "fixture": name,
                "status": p["status"],
                "graph_identity": p["requirements"]["graph_identity"],
                "interaction_requests": len(p["row_associations"]),
                "requests_with_primary_source_rows": sum(
                    bool(r["primary_source_rows"])
                    for r in p["row_associations"].values()
                ),
                "requests_with_prior_silane_candidates": sum(
                    bool(r["prior_silane_candidates"])
                    for r in p["row_associations"].values()
                ),
                "requests_with_partial_rows": sum(
                    bool(r["partial_source_rows"])
                    for r in p["row_associations"].values()
                ),
                "missing_dependencies": len(p["missing_dependencies"]),
            }
        )
    report = {
        "supplement_identity": supplement.identity,
        "fixtures": summaries,
        "execution_profiles": [],
        "backend_executions": [],
        "acceptance_gate": 1,
        "native_pcff_acceptance_gate": 1,
        "limitations": supplement.payload["limitations"],
    }
    if reconstruct:
        if stored_report != report:
            raise SiloxaneParameterError("source_incompatible", "study_report_changed")
    else:
        write(output / "report.json", report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--references", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--reconstruct", action="store_true")
    parser.add_argument("--require-executable", action="store_true")
    args = parser.parse_args()
    try:
        result = run(
            args.catalog,
            args.sha256,
            args.references,
            args.output,
            reconstruct=args.reconstruct,
        )
    except SiloxaneParameterError as error:
        print(json.dumps({"status": error.status, "diagnostics": error.diagnostics}))
        return 1
    except INPUT_ERRORS:
        print(
            json.dumps(
                {
                    "status": "source_incompatible",
                    "diagnostics": [{"stage": "stored_graph_or_evidence_invalid"}],
                }
            )
        )
        return 1
    print(json.dumps(result, indent=2))
    return 1 if args.require_executable else 0


if __name__ == "__main__":
    raise SystemExit(main())
