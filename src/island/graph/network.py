"""Force-field-neutral topology facts for bounded periodic network analysis."""

from collections import Counter

from island.exceptions import AnalysisError


def components(adjacency):
    """Canonical components, independent of insertion order."""
    remaining = set(adjacency)
    result = []
    while remaining:
        pending = [min(remaining)]
        found = set()
        while pending:
            current = pending.pop()
            if current in found:
                continue
            found.add(current)
            pending.extend(adjacency[current] - found)
        remaining -= found
        result.append(tuple(sorted(found)))
    return tuple(result)


def network_topology(graph):
    """Validate authoritative bonds and keep molecular ancestry separate from clusters.

    Degree distribution uses all atoms. Branches, strands and dangling ends use
    the induced heavy-atom skeleton; hydrogens are never polymer strand ends.
    A strand is a maximal edge-disjoint degree-two path, with every crosslink
    endpoint also a junction. Pure degree-two rings are counted once.
    """
    graph.validate_integrity()
    p = graph.payload
    periodic = p["periodic"]
    if periodic is None or periodic["periodic"] != [True, True, True]:
        raise AnalysisError("Network analysis requires three periodic boundaries")
    adjacency = {site["id"]: set() for site in p["sites"]}
    heavy_ids = {site["id"] for site in p["sites"] if site["element"] != "H"}
    for row in p["bonds"]:
        a, b = row["sites"]
        adjacency[a].add(b)
        adjacency[b].add(a)
    heavy = {i: adjacency[i] & heavy_ids for i in sorted(heavy_ids)}
    crosslinks = []
    for row in p["crosslink_provenance"]:
        edge = tuple(sorted(row["sites"]))
        if edge in crosslinks or edge[1] not in adjacency[edge[0]]:
            raise AnalysisError("Duplicate or missing provenance crosslink bond")
        crosslinks.append(edge)
    crosslink_ends = {i for edge in crosslinks for i in edge}
    junctions = {i for i in heavy if len(heavy[i]) != 2} | (crosslink_ends & heavy_ids)
    visited = set()
    strands = []

    def walk(start, neighbour):
        length = 0
        previous, current = start, neighbour
        while True:
            edge = tuple(sorted((previous, current)))
            if edge in visited:
                return length
            visited.add(edge)
            length += 1
            if current in junctions:
                return length
            following = sorted(heavy[current] - {previous})
            if not following:
                return length
            previous, current = current, following[0]

    for start in sorted(junctions):
        for neighbour in sorted(heavy[start]):
            if tuple(sorted((start, neighbour))) not in visited:
                strands.append(walk(start, neighbour))
    for start, neighbours in heavy.items():
        for neighbour in sorted(neighbours):
            if tuple(sorted((start, neighbour))) not in visited:
                strands.append(walk(start, neighbour))
    clusters = components(adjacency)
    return {
        "components": [list(c) for c in clusters],
        "component_membership": {str(i): n for n, c in enumerate(clusters) for i in c},
        "degree_distribution": dict(
            sorted(Counter(map(len, adjacency.values())).items())
        ),
        "heavy_degree_distribution": dict(
            sorted(Counter(map(len, heavy.values())).items())
        ),
        "branch_points": [i for i in heavy if len(heavy[i]) > 2],
        "cluster_size_distribution": dict(sorted(Counter(map(len, clusters)).items())),
        "cycle_rank": len(p["bonds"]) - len(adjacency) + len(clusters),
        "crosslink_bonds": [list(edge) for edge in sorted(crosslinks)],
        "strand_lengths_bonds": sorted(strands),
        "dangling_end_ids": [i for i in heavy if len(heavy[i]) == 1],
        "heavy_atom_count": len(heavy),
    }
