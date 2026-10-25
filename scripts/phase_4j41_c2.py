"""Execute C2 source-completion decisions and fresh C1/silane reference replay."""

import argparse
import shutil
import subprocess
import sys
from hashlib import sha256
from pathlib import Path

from island.charge_references.records import unpack
from island.forcefields.pcff.mixed_sip import load_pcff_sip_distribution
from island.forcefields.pcff.sip_completion import (
    SILANE_PROFILE,
    assess_pcff_sip_completion,
    prepare_pcff_sip_completion,
)
from island.workflows import storage
from island.workflows.bundle import system_from
from scripts.pcff_j41_c2_reference import compiled_searches, source_audit
from scripts.pcff_j41_reference import require


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    path.write_bytes(storage.json_bytes(storage.encode(value)))


def read(path):
    return storage.decode(storage.read_json(path))


def execute(a):
    root = a.output.resolve()
    root.mkdir(exist_ok=False)
    cmd = [
        sys.executable,
        str(Path(__file__).with_name("phase_4j41_c1.py").resolve()),
        "--distribution",
        str(a.distribution.resolve()),
        "--converter",
        str(a.converter.resolve()),
        "--matcher",
        str(a.matcher.resolve()),
        "--lammps",
        str(a.lammps.resolve()),
        "--output",
        str(root / "c1-replay"),
    ]
    run = subprocess.run(cmd, capture_output=True, text=True, check=False)
    (root / "replay.log").write_text(run.stdout + run.stderr)
    write(root / "replay-command.json", {"argv": cmd, "exit_code": run.returncode})
    require(run.returncode == 1, "Expected source-incomplete replay")
    d = load_pcff_sip_distribution(a.distribution)
    for name in d.hashes:
        shutil.copyfile(a.distribution / name, root / name)
    shutil.copyfile(a.matcher, root / "matcher-leaf.so")
    report = {
        "schema": "island.j41-c2-execution.v1",
        "distribution_hashes": d.hashes,
        "matcher_sha256": digest(a.matcher),
        "converter_sha256": digest(a.converter),
        "lammps_sha256": digest(a.lammps),
        "cases": [],
        "acceptance_gate": 1,
        "production_validated": False,
        "simulation_readiness": "not_established",
        "ewald_status": "ewald_unvalidated",
        "pressure_status": "pressure_unavailable",
        "long_npt_status": "blocked",
    }
    queries = {}
    old = read(root / "c1-replay/report.json")
    for prior in old["cases"]:
        name = prior["name"]
        folder = root / name
        folder.mkdir()
        system = system_from(read(root / "c1-replay" / name / "system.json"))
        c1 = unpack((root / "c1-replay" / name / "assessment.json").read_text())
        decision = assess_pcff_sip_completion(system, d, domain=c1["domain"])
        require(
            decision.payload["c1_identity"] == prior["assessment_identity"],
            "Changed C1 identity",
        )
        (folder / "decision.json").write_text(decision.json_text)
        audit = source_audit(system, decision)
        write(folder / "raw-source-audit.json", audit)
        for query in audit.get("unique_searches", []):
            queries[query["family"], query["namespace"], tuple(query["types"])] = query
        row = {
            "name": name,
            "decision_identity": decision.identity,
            "c1_identity": decision.payload["c1_identity"],
            "graph_identity": c1["graph_identity"],
            "source_complete": decision.payload["source_complete"],
            "execution_profile": decision.payload["execution_profile"],
            "blocking_interaction_count": len(
                decision.payload["blocking_interactions"]
            ),
            "typing_status": audit["status"],
            "model_published": False,
            "new_domain_dynamics_executed": False,
            "converter": prior["converter"],
            "workflow_status": "blocked_before_model_publication",
        }
        if decision.payload["execution_authorized"]:
            model = prepare_pcff_sip_completion(
                system, decision, profile=SILANE_PROFILE
            )
            require(
                name == "silane_control" and model.identity == prior["model_identity"],
                "Changed silane model",
            )
            row.update(
                model_identity=model.identity,
                workflow_status="unchanged_silane_fresh_reference_replay_passed",
            )
        report["cases"].append(row)
        print(name, row["workflow_status"], flush=True)
    compiled = compiled_searches(
        d.source.raw, a.matcher.resolve(), [queries[k] for k in sorted(queries)]
    )
    write(root / "compiled-completion-searches.json", compiled)
    report["compiled_search_count"] = len(compiled)
    write(root / "report.json", report)
    return 1


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in (
        "distribution",
        "converter",
        "matcher",
        "lammps",
        "output",
        "reconstruct",
    ):
        p.add_argument("--" + name, type=Path)
    a = p.parse_args()
    if a.reconstruct:
        from validate_pcff_dynamics import reconstruct

        from scripts.validate_phase_4j41_c2_receipt import validate_study

        print(validate_study(a.reconstruct))
        for kind in ("nve", "baoab"):
            reconstruct(a.reconstruct / "c1-replay/baseline" / ("silane-" + kind))
        print("fresh-process C2 reconstruction passed")
        return 0
    require(
        all((a.distribution, a.converter, a.matcher, a.lammps, a.output)),
        "Execution paths required",
    )
    return execute(a)


if __name__ == "__main__":
    raise SystemExit(main())
