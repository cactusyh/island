"""Bounded supporting scans; these never add or merge numerical fit targets."""

from island.pcff_completion import _verified

from .geometry import coordinate, graph_record
from .records import record
from .targets import PROTOCOL


def supporting_scans(system, assessment, fixture):
    p = _verified(system, assessment)
    elements = {i: system.topology.sites[i].element for i in p["atom_order"]}
    neighbors = {i: set() for i in p["atom_order"]}
    for b in system.topology.bonds.values():
        neighbors[b.site1].add(b.site2)
        neighbors[b.site2].add(b.site1)
    paths = {2: [], 3: [], 4: []}

    def visit(path):
        if len(path) >= 2:
            paths[len(path)].append(path)
        if len(path) == 4:
            return
        for i in sorted(neighbors[path[-1]] - set(path)):
            visit(path + [i])

    for i in sorted(neighbors):
        visit([i])
    scans = []

    def choose(label, pattern):
        match = [
            path
            for path in paths[len(pattern)]
            if [elements[i] for i in path] == pattern
        ]
        if match:
            return min(match)
        return None

    specifications = []
    for label, kind, pattern in (
        ("Si-O-bond", "bond", ["Si", "O"]),
        ("Si-C-bond", "bond", ["Si", "C"]),
        ("Si-O-Si-angle", "angle", ["Si", "O", "Si"]),
        ("Si-O-Si-O-torsion", "torsion", ["Si", "O", "Si", "O"]),
        ("terminal-C-Si-O-Si-torsion", "torsion", ["C", "Si", "O", "Si"]),
    ):
        path = choose(label, pattern)
        if path:
            specifications.append((label, [coordinate(kind, path)]))
    path = choose("Si-O-Si", ["Si", "O", "Si"])
    if path:
        specifications.append(
            (
                "joint-bond-bond",
                [coordinate("bond", path[:2]), coordinate("bond", path[1:])],
            )
        )
        specifications.append(
            (
                "joint-bond-angle",
                [coordinate("bond", path[:2]), coordinate("angle", path)],
            )
        )
    path = choose("Si-O-Si-O", ["Si", "O", "Si", "O"])
    if path is None:
        path = choose("C-Si-O-Si", ["C", "Si", "O", "Si"])
    if path is not None:
        specifications.extend(
            [
                (
                    "joint-angle-torsion",
                    [coordinate("angle", path[:3]), coordinate("torsion", path)],
                ),
                (
                    "joint-end-bond-torsion",
                    [coordinate("bond", path[:2]), coordinate("torsion", path)],
                ),
                (
                    "joint-middle-bond-torsion",
                    [coordinate("bond", path[1:3]), coordinate("torsion", path)],
                ),
                (
                    "joint-angle-angle-torsion",
                    [
                        coordinate("angle", path[:3]),
                        coordinate("angle", path[1:]),
                        coordinate("torsion", path),
                    ],
                ),
            ]
        )
    # AA supporting target uses an exact source-missing site, never a path angle
    # with an invented silicon environment.
    row = next((r for r in p["missing_terms"] if r["family"] == "angle-angle"), None)
    if row:
        from .geometry import target_scan

        specifications.append(("joint-angle-angle", target_scan(row)["coordinates"]))
    for label, coords in specifications:
        # 3-point auxiliary grid; all numerical targets retain their full 5x5 grid.
        for c in coords:
            c["displacements"] = c["displacements"][::2]
            c["scan_spacing"] *= 2
        scans.append(
            record(
                {
                    "schema": "island.pcff-qm-supporting-scan.v1",
                    "fixture": fixture,
                    "label": label,
                    "graph_identity": p["graph_identity"],
                    "graph": graph_record(system, p),
                    "atom_order": p["atom_order"],
                    "charge_identity": p["charge_identity"],
                    "charge_state": 0,
                    "spin_state": PROTOCOL["spin_state"],
                    "qm_protocol": PROTOCOL,
                    "scan": {
                        "coordinates": coords,
                        "purpose": "supporting diagnostic only",
                        "joint_fit": False,
                    },
                    "creates_missing_target": False,
                    "QM_executed": False,
                    "status": "input_prepared",
                    "execution_authority": False,
                }
            )
        )
    return tuple(scans)
