"""Predeclared J9 final-graph matrix: construction/perception/source diagnostics.

Only the separately declared P_DP2 vertical slice earns numerical/workflow gates.
This matrix does not promote every observed type or repeat to general coverage.
"""

import argparse
from collections import deque
from pathlib import Path

from island.builders import (
    build_linear_polymer,
    build_polymer_from_sequence,
    build_random_copolymer,
)
from island.forcefields import PCFFOptions
from island.forcefields.pcff import (
    edit_polymer_topology,
    prepare_pcff_polymer,
    select_pcff_profile,
)
from island.forcefields.pcff.registry import LINKED_BENZENOID
from island.forcefields.pcff.source import PIN
from island.workflows import storage
from island.workflows.bundle import system_data


def endpoints(system):
    carbons = sorted(i for i, a in system.topology.sites.items() if a.element == "C")
    neighbors = {i: [] for i in system.topology.sites}
    for b in system.topology.bonds.values():
        neighbors[b.site1].append(b.site2)
        neighbors[b.site2].append(b.site1)
    for i in carbons:
        if system.topology.sites[i].metadata["repeat_unit_index"] != 0:
            continue
        h = [j for j in neighbors[i] if system.topology.sites[j].element == "H"]
        if len(h) != 1:
            continue
        distances, todo = {i: 0}, deque([i])
        while todo:
            k = todo.popleft()
            for j in neighbors[k]:
                if j not in distances:
                    distances[j] = distances[k] + 1
                    todo.append(j)
        for j in carbons:
            if (
                system.topology.sites[j].metadata["repeat_unit_index"] != 2
                or distances[j] < 6
            ):
                continue
            hj = [n for n in neighbors[j] if system.topology.sites[n].element == "H"]
            if len(hj) == 1:
                return (h[0], hj[0]), (i, j)
    raise ValueError("No declared crosslink endpoints")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument(
        "--declaration", default="docs/evidence/phase_4j9_polymer_declaration.json"
    )
    a = p.parse_args()
    d = storage.read_json(a.declaration)
    if storage.checksum(Path(a.source).read_bytes()) != d["source_sha256"]:
        raise ValueError("Wrong source identity")
    a.output.mkdir(parents=True, exist_ok=False)
    storage.publish(a.output / "declaration.json", storage.json_bytes(d))
    options = PCFFOptions(
        a.source,
        (0, 0, 1),
        (0, 0, 1),
        typing_profile="island_pcff_source_graph_v1",
        source_profile=select_pcff_profile(LINKED_BENZENOID, sha256=PIN["sha256"]),
    )
    report = []
    for case in d["matrix"]:
        work = a.output / case["name"]
        work.mkdir()
        outcome = {"case": case["name"], "numerical_verified": False}
        try:
            if case["name"] == "random_PM":
                system = build_random_copolymer(
                    {k: d["psmiles"][k] for k in ("P", "M")},
                    dp=4,
                    fractions={"P": 0.5, "M": 0.5},
                    sequence_seed=2026,
                    random_seed=2026,
                )
            elif case["name"].startswith("P_DP") or case["name"] == "crosslinked_P_DP3":
                system = build_linear_polymer(
                    d["psmiles"]["P"], dp=len(case["sequence"]), random_seed=2026
                )
            else:
                system = build_polymer_from_sequence(
                    d["psmiles"],
                    case["sequence"],
                    random_seed=2026,
                    head_end_group=case.get("head_end_group"),
                    tail_end_group=case.get("tail_end_group"),
                )
            if case["name"] == "crosslinked_P_DP3":
                h, bond = endpoints(system)
                system = edit_polymer_topology(
                    system, remove_hydrogens=h, add_crosslinks=[bond]
                )
                outcome["fixed_edit"] = {"removed_hydrogens": h, "crosslink": bond}
            storage.publish(
                work / "system.json",
                storage.json_bytes(storage.encode(system_data(system))),
            )
            before = system.to_dict()
            diagnostic = prepare_pcff_polymer(system, options, mode="diagnostic")
            storage.publish(
                work / "diagnostic.json", storage.json_bytes(storage.encode(diagnostic))
            )
            outcome.update(
                model_complete=diagnostic.get("model_complete", False),
                sites=system.number_of_sites,
                graph_unchanged=system.to_dict() == before,
                source_sha256=PIN["sha256"],
                system_sha256=storage.checksum((work / "system.json").read_bytes()),
            )
            if outcome["model_complete"]:
                result = prepare_pcff_polymer(system, options)
                storage.publish(work / "final-graph.json", result.json_text.encode())
                outcome.update(
                    strict_prepared=True,
                    native=result.prepared.native_result.identity,
                    prepared=result.prepared.identity,
                    evidence=result.identity,
                    roles={
                        k: sum(
                            bool(r["roles"][k]) for r in result.payload["interactions"]
                        )
                        for k in ("crosslink", "inter_repeat", "end_group")
                    },
                )
        except Exception as exc:  # noqa: BLE001 -- retain diagnostic evidence, never replace cases
            outcome["error"] = type(exc).__name__ + ": " + str(exc)
        report.append(outcome)
        storage.publish(work / "outcome.json", storage.json_bytes(outcome))
    storage.publish(
        a.output / "report.json",
        storage.json_bytes({"cases": report, "full_source_complete": False}),
    )
    print(
        [(r["case"], r.get("strict_prepared", False), r.get("error")) for r in report]
    )
    return (
        0
        if all(r.get("strict_prepared") and r["graph_unchanged"] for r in report)
        else 1
    )


if __name__ == "__main__":
    raise SystemExit(main())
