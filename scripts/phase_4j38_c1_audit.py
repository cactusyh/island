"""Independent raw-topology/Foyer/XML audit. No ISLAND charge/graph importer."""

import hashlib
import json
import math
import sys
import xml.etree.ElementTree as ET
from collections import Counter
from decimal import Decimal
from pathlib import Path


def require(condition, diagnostic):
    if not condition:
        raise ValueError(diagnostic)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(
        json.dumps(
            value, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode()
    ).hexdigest()


def components(graph):
    ids = [a["id"] for a in graph["sites"]]
    require(ids == sorted(set(ids)), "Duplicate, omitted or reordered atom IDs")
    adj = {i: set() for i in ids}
    edges = set()
    for b in graph["bonds"]:
        a, c = b["sites"]
        edge = tuple(sorted((a, c)))
        require(a != c and edge not in edges, "Invalid or duplicate bond")
        edges.add(edge)
        adj[a].add(c)
        adj[c].add(a)
    remaining, groups = set(ids), []
    while remaining:
        seen, pending = set(), [min(remaining)]
        while pending:
            i = pending.pop()
            if i not in seen:
                seen.add(i)
                pending.extend(adj[i] - seen)
        groups.append(sorted(seen))
        remaining -= seen
    require(groups == graph["components"], "Wrong component membership")
    require(
        graph["component_membership"]
        == {str(i): k for k, group in enumerate(groups) for i in group},
        "Wrong component index mapping",
    )
    return adj, groups


def aromatic(evidence, xml, pin_path="src/island/forcefields/oplsaa/pin.json"):
    import foyer
    import parmed
    from rdkit import Chem

    pin = json.loads(Path(pin_path).read_text())
    require(sha(xml) == pin["xml_sha256"], "Changed authoritative XML source")
    upstream = Path(foyer.__file__).parent
    require(foyer.__version__ == pin["foyer_version"], "Wrong Foyer version")
    for file, expected in pin["implementation_files"].items():
        require(sha(upstream / file) == expected, "Changed upstream typer: " + file)
    graph, plan = evidence["final_graph"], evidence["packing_plan"]
    chemical = {
        k: v
        for k, v in graph.items()
        if k not in ("graph_identity", "history_identity")
    }
    chemical["transformations"] = []
    require(
        graph["graph_identity"]
        == digest({"payload": chemical, "sha256": digest(chemical)}),
        "Changed final graph identity",
    )
    adj, groups = components(graph)
    sites = {a["id"]: a for a in graph["sites"]}
    require(len(groups) == 2 and len(sites) == 44, "Wrong frozen aromatic inventory")
    require(
        plan["output_graph_identity"] == graph["graph_identity"], "Stale packed graph"
    )
    # Independently rebuild the external typing topology directly from raw graph.
    p = parmed.Structure()
    atoms = {}
    for i, a in sites.items():
        atom = parmed.Atom(
            name=a["element"], atomic_number=a["atomic_number"], mass=a["mass"]
        )
        p.add_atom(atom, "M", 1)
        atoms[i] = atom
    for b in graph["bonds"]:
        a, c = b["sites"]
        p.bonds.append(parmed.Bond(atoms[a], atoms[c], order=b["order"]))
    typemap = foyer.Forcefield(forcefield_files=str(xml)).run_atomtyping(
        p, use_residue_map=False
    )
    root = ET.parse(xml)
    table = {a.attrib["type"]: a.attrib for a in root.findall("./NonbondedForce/Atom")}
    definitions = {a.attrib["name"]: a.attrib for a in root.findall("./AtomTypes/Type")}
    rows, summaries = [], []
    for index, ids in enumerate(groups):
        chain = evidence["complete_chains"][index]
        mapping = plan["site_id_remapping"][index]
        local = {a["id"]: a for a in chain["sites"]}
        require(
            sorted(mapping.values()) == ids and len(mapping) == len(local),
            "Packing map omitted/duplicated atoms",
        )
        require(set(mapping) == {str(i) for i in local}, "Wrong local atom mapping")
        for old, new in mapping.items():
            a, b = local[int(old)], sites[new]
            for field in ("element", "formal_charge", "metadata"):
                require(
                    a[field] == b[field], "Repeat/end-group remapping changed: " + field
                )
            require(
                graph["molecule_membership"][str(new)]
                == plan["molecule_membership_remapping"][index][old],
                "Wrong molecule membership",
            )
        chains = {sites[i]["metadata"]["chain_id"] for i in ids}
        require(chains == {f"J38-aromatic-{index}"}, "Wrong chain membership")
        # An independent known complete biphenyl chemical graph, not a second
        # invocation of the polymer builder. Explicit H and formal charge checked.
        rw = Chem.RWMol()
        rdmap = {}
        for i in ids:
            a = sites[i]
            atom = Chem.Atom(a["element"])
            atom.SetFormalCharge(a["formal_charge"])
            atom.SetNoImplicit(True)
            rdmap[i] = rw.AddAtom(atom)
        for b in graph["bonds"]:
            a, c = b["sites"]
            if a in rdmap:
                rw.AddBond(
                    rdmap[a],
                    rdmap[c],
                    Chem.BondType.AROMATIC if b["aromatic"] else Chem.BondType.SINGLE,
                )
        mol = rw.GetMol()
        Chem.SanitizeMol(mol)
        expected = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1-c1ccccc1"))
        require(
            Chem.MolToSmiles(mol) == Chem.MolToSmiles(expected),
            "Incomplete or wrong end-group graph",
        )
        require(
            Counter(sites[i]["element"] for i in ids) == {"C": 12, "H": 10},
            "Missing/extra atoms",
        )
        require(
            len({graph["molecule_membership"][str(i)] for i in ids}) == 1,
            "Split molecule membership",
        )
        for repeat in (0, 1):
            members = [
                i for i in ids if sites[i]["metadata"]["repeat_unit_index"] == repeat
            ]
            require(
                Counter(sites[i]["element"] for i in members) == {"C": 6, "H": 5},
                "Wrong repeat/end H mapping",
            )
        endpoints = [
            mapping[str(chain["metadata"]["polymer"][k])]
            for k in ("head_site_id", "tail_site_id")
        ]
        require(endpoints == evidence["endpoints"][index], "Changed endpoint mapping")
        require(
            all(
                sum(sites[j]["element"] == "H" for j in adj[i]) == 1 for i in endpoints
            ),
            "Missing end-group hydrogen",
        )
        for i in ids:
            a = sites[i]
            source_type = typemap[list(sites).index(i)]["atomtype"]
            text = table[source_type]["charge"]
            rows.append(
                {
                    "atom_id": i,
                    "element": a["element"],
                    "component": index,
                    "molecule": graph["molecule_membership"][str(i)],
                    "formal_charge": a["formal_charge"],
                    "metadata": a["metadata"],
                    "neighbors": sorted(adj[i]),
                    "source_type": source_type,
                    "source_charge_decimal": text,
                    "charge_e": float(text),
                    "source_location": f'NonbondedForce/Atom[@type="{source_type}"]',
                    "source_definition": definitions[source_type],
                    "attached_hydrogens": sum(
                        sites[j]["element"] == "H" for j in adj[i]
                    ),
                }
            )
        subset = [r for r in rows if r["component"] == index]
        total = sum((Decimal(r["source_charge_decimal"]) for r in subset), Decimal())
        formal = sum(sites[i]["formal_charge"] for i in ids)
        require(
            formal == 0 and total == Decimal("-0.230"),
            "Frozen source residual changed; requires separate review",
        )
        summaries.append(
            {
                "atom_ids": ids,
                "chain": next(iter(chains)),
                "formal_total_e": formal,
                "source_total_decimal_e": str(total),
                "residual_decimal_e": str(total - formal),
                "fsum_e": math.fsum(r["charge_e"] for r in subset),
                "carbon_without_hydrogen": [
                    r["atom_id"]
                    for r in subset
                    if r["element"] == "C" and r["attached_hydrogens"] == 0
                ],
                "formula": "C12H10",
                "endpoints": endpoints,
            }
        )
    return {
        "schema": "island.j38-c1-aromatic-source-audit.v1",
        "final_graph_identity": graph["graph_identity"],
        "packing_identity": plan["identity"],
        "evidence_identity": digest(evidence),
        "source_sha256": sha(xml),
        "pin_sha256": sha(pin_path),
        "units": "elementary_charge",
        "typing": "independent raw graph -> pinned Foyer.run_atomtyping",
        "rows": rows,
        "components": summaries,
        "total_decimal_e": "-0.460",
        "classification": "source_incomplete",
        "stage": "charge_completeness",
        "cause": "Pinned aromatic carbon rule selects -0.115 for all 12 carbons, including two inter-ring ipso carbons without H; ten +0.115 H leave -0.230 per chain. No exact-graph charge-complete source set is available.",
        "normalized": False,
        "redistributed": False,
        "dynamics_executed": False,
        "runtime": {
            "python": sys.version,
            "executable": str(Path(sys.executable).resolve()),
            "executable_sha256": sha(Path(sys.executable).resolve()),
            "foyer_version": foyer.__version__,
        },
    }


def resp_inventory(requested):
    """Audit declared retained raw archives; never transfer a charge vector."""
    import tarfile
    import tempfile

    from scripts.phase_4j37_reference import audit

    anchors_path = Path("scripts/phase_4j37_anchors.json")
    anchors = json.loads(anchors_path.read_text())
    runtime = json.loads(Path("docs/evidence/phase_4j36/environment.json").read_text())
    inventory = []
    with tempfile.TemporaryDirectory(prefix="j38-c1-raw-resp-") as directory:
        for archive, expected in anchors["archives"].items():
            require(sha(archive) == expected, "Changed raw quantum archive")
            with tarfile.open(archive) as handle:
                handle.extractall(directory, filter="data")
        for kind in ("linear", "network"):
            charges_path = Path(
                f"docs/evidence/phase_4j36_c2/memory-{kind}-charges.json"
            )
            require(
                sha(charges_path) == anchors["cases"][kind]["charges_sha256"],
                "Changed reviewed charge record",
            )
            charges = json.loads(charges_path.read_text())
            report = audit(directory, kind, charges, runtime)
            data = json.loads(
                (Path(directory) / f"quantum-{kind}/input.json").read_text()
            )
            inventory.append(
                {
                    "kind": kind,
                    "final_graph_identity": data["final_graph_identity"],
                    "input_sha256": sha(Path(directory) / f"quantum-{kind}/input.json"),
                    "charge_record_sha256": sha(charges_path),
                    "independent_review_identity": report["identity"],
                    "atom_count": len(data["graph"]["sites"]),
                }
            )
    results = {}
    for case, graph_id in requested.items():
        matches = [r for r in inventory if r["final_graph_identity"] == graph_id]
        require(
            not matches,
            "New exact graph-specific asset requires separate review before this frozen correction can accept it",
        )
        results[case] = {
            "status": "asset_unavailable",
            "stage": "charge_provenance",
            "final_graph_identity": graph_id,
            "charge_origin": "RESP/Psi4",
            "exact_graph_match": False,
            "dynamics_executed": False,
            "reason": "No independently reviewed raw Psi4/RESP geometry, ESP, charge and executable provenance for this exact final graph in the declared retained asset inventory. PEO artifacts cannot be transferred.",
        }
    return {
        "schema": "island.j38-c1-exact-resp-inventory.v1",
        "anchors_sha256": sha(anchors_path),
        "archives": anchors["archives"],
        "runtime_sha256": sha("docs/evidence/phase_4j36/environment.json"),
        "inventory": inventory,
        "requested": results,
        "search_scope": "Declared reviewed J36/J37 raw archive inventory only; no claim of an exhaustive filesystem search",
    }
