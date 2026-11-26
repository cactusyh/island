"""Periodic candidate search, integer cycle closure, and functional-group statistics."""

from collections import Counter
from dataclasses import asdict

import numpy as np
from rdkit import Chem
from scipy.spatial import cKDTree

from .chemistry import AMINE, epoxides


def image(a, b, box):
    return -np.rint((b - a) / box).astype(int)


def select(
    mol,
    xyz,
    box,
    cutoff,
    rng,
    limit,
    min_ring=6,
    probability=1.0,
    template=AMINE,
    donor_angle_tolerance=None,
):
    if not 0 < cutoff < min(box) / 2 or not 0 <= probability <= 1 or limit < 0:
        raise ValueError("Invalid periodic reaction search settings")
    rings = epoxides(mol)
    ns = {n: h for n, h in template.donors(mol).items() if h}
    candidates = []
    rejected = Counter()
    if rings and ns:
        nids = list(ns)
        cells = cKDTree(xyz[nids] % box, boxsize=box)
        for c, o in rings.items():
            for j in sorted(cells.query_ball_point(xyz[c] % box, cutoff)):
                n = nids[j]
                shift = image(xyz[c], xyz[n], box)
                dx = xyz[n] + shift * box - xyz[c]
                candidates.append((np.linalg.norm(dx), c, o, n, shift))
    # Random tie breaker makes the policy insensitive to periodic roundoff.
    ranks = rng.random(len(candidates))
    order = sorted(
        range(len(candidates)), key=lambda j: (round(candidates[j][0], 8), ranks[j])
    )
    events = []
    records = []
    used = set()
    provisional = Chem.RWMol(mol)
    for j in order:
        distance, c, o, n, shift = candidates[j]
        if c in used or n in used:
            rejected["site_conflict"] += 1
            continue
        if len(events) >= limit:
            rejected["batch_limit"] += 1
            continue
        geometry = {}
        if (
            template.product_donor_angle_deg is not None
            and donor_angle_tolerance is not None
        ):
            # The existing C-S and proposed S-C vectors must use the same image of S.
            (neighbor,) = tuple(
                a.GetIdx()
                for a in mol.GetAtomWithIdx(n).GetNeighbors()
                if a.GetAtomicNum() != 1
            )
            u = xyz[neighbor] + image(xyz[n], xyz[neighbor], box) * box - xyz[n]
            v = xyz[c] - shift * box - xyz[n]
            cosine = (
                np.dot(u, v) / (np.linalg.norm(u) * distance)
                if distance and np.linalg.norm(u)
                else np.nan
            )
            angle = float(np.degrees(np.arccos(np.clip(cosine, -1, 1))))
            if (
                not np.isfinite(angle)
                or abs(angle - template.product_donor_angle_deg) > donor_angle_tolerance
            ):
                rejected["donor_angle"] += 1
                continue
            geometry["donor_angle_deg"] = angle
        provisional.RemoveBond(c, o)
        path = Chem.GetShortestPath(provisional, c, n)
        provisional.AddBond(c, o, Chem.BondType.SINGLE)
        if path and len(path) < min_ring:
            rejected["short_ring"] += 1
            continue
        if rng.random() >= probability:
            rejected["probability"] += 1
            continue
        e = template.event_type(c, o, n, ns[n][0])
        template.apply(provisional, e)
        events.append(e)
        records.append(
            dict(
                **asdict(e),
                **geometry,
                image=shift.tolist(),
                distance_A=float(distance),
                stage=template.max_additions + 1 - len(ns[n]),
                same_component=bool(path),
            )
        )
        used.update([c, n])
    return events, {
        "candidate_count": len(candidates),
        "rejected": dict(rejected),
        "events": records,
        "box_A": list(box),
    }


def statistics(
    mol,
    xyz,
    box,
    nd,
    npac,
    template=AMINE,
    *,
    initial_epoxides=None,
    initial_donor_h=None,
):
    adj = [[n.GetIdx() for n in a.GetNeighbors()] for a in mol.GetAtoms()]
    wraps = []
    components = Chem.GetMolFrags(mol)
    for group in components:
        offsets = {group[0]: np.zeros(3, dtype=int)}
        stack = [group[0]]
        while stack:
            i = stack.pop()
            for j in adj[i]:
                proposed = offsets[i] + image(xyz[i], xyz[j], box)
                if j not in offsets:
                    offsets[j] = proposed
                    stack.append(j)
                elif i < j and np.any(delta := proposed - offsets[j]):
                    wraps.append(delta.tolist())
    # Count network strands between branch/end sites on the heavy skeleton.
    heavy = {
        a.GetIdx(): [
            j for j in adj[a.GetIdx()] if mol.GetAtomWithIdx(j).GetAtomicNum() > 1
        ]
        for a in mol.GetAtoms()
        if a.GetAtomicNum() > 1
    }
    from island.graph.network import components as graph_components

    # Monomer-level multigraph distinguishes chemical aromatic/cyclohexyl rings
    # from cycles created by crosslinking. Preserve parallel edges.
    owner = [a.GetIntProp("monomer") for a in mol.GetAtoms()]
    edges = []
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if owner[i] != owner[j] and i in heavy and j in heavy:
            edges.append((owner[i], owner[j]))
    ma = {i: set() for i in set(owner)}
    incident = {i: [] for i in ma}
    for k, (a, b) in enumerate(edges):
        ma[a].add(b)
        ma[b].add(a)
        incident[a].append((k, b))
        incident[b].append((k, a))
    deg = {i: len(v) for i, v in incident.items()}
    seen = set()
    lengths = []
    for root in sorted(ma, key=lambda i: deg[i] == 2):
        for k, j in incident[root]:
            length = 0
            while k not in seen:
                seen.add(k)
                length += 1
                if deg[j] != 2:
                    break
                k, j = next((e, n) for e, n in incident[j] if e != k)
            if length:
                lengths.append(length)
    nh = Counter(map(len, template.donors(mol).values()))
    remaining = len(epoxides(mol))
    nh0 = (
        initial_donor_h
        if initial_donor_h is not None
        else template.sites_per_hardener * template.max_additions * npac
    )
    conversion = 1 - sum(k * v for k, v in nh.items()) / nh0
    groups = (
        {
            "amine_h_conversion": conversion,
            "primary": nh[2],
            "secondary": nh[1],
            "tertiary": nh[0],
        }
        if template == AMINE
        else {
            "thiol_h_conversion": conversion,
            "remaining_thiols": nh[1],
            "formed_thioethers": nh[0],
        }
    )
    return {
        "epoxy_conversion": 1
        - remaining / (initial_epoxides if initial_epoxides is not None else 2 * nd),
        **groups,
        "remaining_epoxides": remaining,
        "crosslinks": len(edges),
        "largest_component_atom_fraction": max(map(len, components))
        / mol.GetNumAtoms(),
        "monomer_cycle_rank": len(edges) - len(ma) + len(graph_components(ma)),
        "dangling_monomers": sum(d == 1 for d in deg.values()),
        "strand_lengths_monomer_bonds": sorted(lengths),
        "heavy_terminal_sites": sum(len(v) == 1 for v in heavy.values()),
        "closure_vectors": wraps,
        "winding_rank": int(np.linalg.matrix_rank(wraps)) if wraps else 0,
        "volume_A3": float(np.prod(box)),
        "density_g_cm3": sum(a.GetMass() for a in mol.GetAtoms())
        / (0.602214076 * np.prod(box)),
    }


def bond_images(mol, xyz, box):
    """Nonzero minimum-image shifts for wrapped sites i<j. Zero is implicit."""
    positions = np.asarray(xyz) % box
    if not np.isfinite(positions).all() or not np.isfinite(box).all() or min(box) <= 0:
        raise ValueError("Invalid periodic coordinates")
    result = []
    for b in mol.GetBonds():
        i, j = sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))
        delta = (positions[j] - positions[i]) / box
        if np.any(np.isclose(np.abs(delta), 0.5, rtol=0, atol=1e-10)):
            raise ValueError("Ambiguous half-box bond image")
        shift = -np.rint(delta).astype(int)
        if np.any(shift):
            result.append([i, j, *shift.tolist()])
    return sorted(result)
