"""One declared, retained-template GAFF2 provided-charge interoperability run."""

import argparse
from math import fsum, isfinite
from pathlib import Path

from island.fragment_charges import load_fragment_record
from island.workflows import storage
from island.workflows.bundle import system_from


def compare_charges(expected, actual):
    """Case-specific gates; never change charges or engine defaults."""
    if set(expected) != set(actual) or not expected:
        raise ValueError("Charge coverage mismatch")
    if not all(isfinite(v) for v in (*expected.values(), *actual.values())):
        raise ValueError("Nonfinite charges")
    maximum = max(abs(expected[s] - actual[s]) for s in expected)
    total = fsum(actual.values())
    if maximum > 1e-8 or abs(total) > 1e-6:
        raise ValueError(
            "Serialization acceptance failed: per-site or molecular charge"
        )
    return {"total_e": total, "max_per_site_difference_e": maximum}


def inputs(source):
    outcomes = storage.read_json(source / "outcomes.json")
    row = next(r for r in outcomes["cases"] if r["case"] == "pe")
    for relative, checksum in row["file_sha256"].items():
        if storage.checksum((source / "pe" / relative).read_bytes()) != checksum:
            raise ValueError(f"Source checksum mismatch: {relative}")
    template = load_fragment_record(
        source / "pe/template.json", expected_cache_key=row["cache_key"]
    )
    assignment = load_fragment_record(source / "pe/integration/assignment.json")
    if assignment.payload["template"] != template.payload:
        raise ValueError("Assignment uses a different retained template")
    return template, assignment, row


def declaration(source):
    template, assignment, row = inputs(source)
    q = assignment.payload["data"]["assignments"]
    # This specific supplied vector lies on the six-decimal MOL2 grid. For an
    # arbitrary 20-atom vector six-decimal rounding alone could reach 1e-5 e.
    if len(q) != 20 or max(abs(v) for v in q.values()) >= 1:
        raise ValueError("Not the declared bounded PE DP3 vector")
    quantization = fsum(abs(v - float(f"{v:.6f}")) for v in q.values())
    # E16.8: |18.2223*q| < 18.2223; half ULP <= 5e-8 scaled units.
    bound = quantization + 20 * (5e-11 + 5e-8 / 18.2223)
    if bound >= 1e-6:
        raise ValueError("Unsupported serialization bound")
    return {
        "schema": "island_fragment_interoperability_v1",
        "template_identity": template.identity,
        "assignment_identity": assignment.identity,
        "source_outcomes_sha256": storage.checksum(
            (source / "outcomes.json").read_bytes()
        ),
        "source_files": row["file_sha256"],
        "force_field": row["force_field"],
        "method": "provided",
        "timeout_seconds": 600,
        "outer_retries": 0,
        "fragment_qm_calls": 0,
        "molecular_tolerance_e": 1e-6,
        "per_site_tolerance_e": 1e-8,
        "algebraic_tolerance_e": 1e-12,
        "mol2_vector_quantization_bound_e": quantization,
        "conservative_molecular_serialization_bound_e": bound,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }


def execute(source, output, amberhome):
    import parmed

    from island.forcefields import AmberToolsOptions, AmberToolsParameterizationEngine

    plan = storage.read_json(output / "declaration.json")
    if declaration(source) != plan:
        raise ValueError("Declaration or retained inputs changed")
    if {p.name for p in output.iterdir()} != {"declaration.json"}:
        raise ValueError("Refuse retry or existing execution output")
    template, assignment, _ = inputs(source)
    p = assignment.payload
    system = system_from(p["target_system"])
    q = p["data"]["assignments"]
    result = {"status": "failed", "declaration": plan}
    try:
        prep = AmberToolsParameterizationEngine().parameterize(
            system,
            AmberToolsOptions(
                "gaff2",
                "provided",
                provided_charges=q,
                charge_source="fragment-derived; assignment=" + assignment.identity,
                charge_tolerance=1e-6,
                timeout_seconds=600,
                amberhome=amberhome,
                work_root=output,
                retain_success_artifacts=True,
            ),
        )
        prep.validate_integrity(system)
        snapshot = prep.to_parameterized_system(system)
        if (
            prep.record["force_field_data"]["sha256"]
            != plan["force_field"]["data_sha256"]
        ):
            raise ValueError("GAFF2 data identity mismatch")
        directory = next(output.glob("island-ambertools-*"))
        if list(directory.glob("sqm*")):
            raise ValueError("Unexpected SQM artifact in provided path")
        parm = parmed.load_file(str(directory / "result.prmtop"))
        mapping = prep.record["lineage"]["prmtop_index_to_site_id"]
        actual = {mapping[i]: a.charge for i, a in enumerate(parm.atoms)}
        result.update(compare_charges(q, actual))
        result.update(
            status="passed",
            record_signature=prep.record_signature,
            import_signature=prep.imported_result.result_signature,
            snapshot_created=True,
            snapshot_sites=snapshot.system.number_of_sites,
            provided_charge_provenance={
                "assignment_identity": assignment.identity,
                "template_identity": template.identity,
                "transformations": {
                    k: template.payload[k] for k in ("policy", "hydrogen_policy")
                },
            },
        )
        storage.publish(
            output / "preparation.json",
            storage.json_bytes(storage.encode(dict(prep.record))),
        )
    except Exception as error:  # noqa: BLE001 -- durable failed gate, no retry
        result["failure"] = f"{type(error).__name__}: {error}"
    result["artifact_sha256"] = {
        str(p.relative_to(output)): storage.checksum(p.read_bytes())
        for p in output.rglob("*")
        if p.is_file()
    }
    inputs(source)  # original failed experiment remains byte-for-byte intact
    storage.publish(output / "outcome.json", storage.json_bytes(storage.encode(result)))
    return result["status"] == "passed"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--amberhome")
    parser.add_argument("--execute-declared", action="store_true")
    args = parser.parse_args()
    if args.execute_declared:
        return 0 if execute(args.source, args.output, args.amberhome) else 1
    plan = declaration(args.source)
    args.output.mkdir(parents=True, exist_ok=False)
    storage.publish(args.output / "declaration.json", storage.json_bytes(plan))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
