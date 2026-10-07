"""J18 graph acceptance: deterministic chains, explicit crosslinks, periodic bundle.

Force-field assignment is attempted only on the supported nonperiodic linear
slice. Periodic and crosslinked graph persistence is neutral and conservative.
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path

from island.graph import (
    build_psmiles_graph,
    construct_periodic_box,
    create_crosslinks,
    load_final_graph_bundle,
    save_final_graph_bundle,
)
from island.workflows import storage


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=False)
    out = {
        "schema": "island_j18_graph_acceptance_v1",
        "passed": False,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    try:
        s3, g3, generation3 = build_psmiles_graph(
            "[*:1]CC[*:2]", dp=3, random_seed=2026
        )
        _s10, g10, generation10 = build_psmiles_graph(
            "[*:1]CC[*:2]", dp=10, random_seed=2026
        )
        assert (
            g3.identity
            == build_psmiles_graph("[*:1]CC[*:2]", dp=3, random_seed=2026)[1].identity
        )
        assert len(g3.payload["sites"]) == 20 and len(g10.payload["sites"]) == 62
        carbons = [i for i, atom in s3.topology.sites.items() if atom.element == "C"]
        targets = (carbons[0], carbons[-1])
        for carbon in targets:
            hydrogen = next(
                i
                for i in s3.topology.neighbors(carbon)
                if s3.topology.sites[i].element == "H"
            )
            s3.topology.remove_site(hydrogen)
            s3.coordinates._positions.pop(hydrogen)
        crosslinked, _gc, crosslink = create_crosslinks(s3, [targets], seed=2026)
        periodic, gp, _box = construct_periodic_box(
            crosslinked, (40, 40, 40), seed=2026
        )
        bundle = a.output / "periodic-bundle"
        save_final_graph_bundle(periodic, gp, bundle)
        moved = a.output / "relocated-bundle"
        bundle.rename(moved)
        loaded_system, loaded_graph, _ = load_final_graph_bundle(moved)
        assert loaded_graph.identity == gp.identity and loaded_system.box.lengths == (
            40,
            40,
            40,
        )
        # Child reconstruction forbids RDKit/force-field/scientific imports.
        child = a.output / "child.json"
        code = """
import builtins, json, sys
orig=builtins.__import__
def blocked(name,*a,**k):
    if name.split('.')[0] in {'rdkit','openmm','scipy','foyer','parmed'}: raise AssertionError(name)
    return orig(name,*a,**k)
builtins.__import__=blocked
from island.graph import load_final_graph_bundle
s,g,q=load_final_graph_bundle(sys.argv[1])
print(json.dumps({'graph':g.identity,'periodic':g.periodic,'sites':len(g.payload['sites'])}))
"""
        result = subprocess.run(
            [sys.executable, "-c", code, str(moved)],
            capture_output=True,
            text=True,
            check=False,
        )
        assert result.returncode == 0, result.stderr
        child.write_text(result.stdout)
        out.update(
            passed=True,
            dp3_sites=len(g3.payload["sites"]),
            dp10_sites=len(g10.payload["sites"]),
            generation_identity=generation3.identity,
            dp10_generation_identity=generation10.identity,
            crosslink_identity=crosslink.identity,
            periodic_identity=gp.identity,
            bundle_identity=gp.identity,
            child_process=True,
            scientific_imports_blocked=True,
            forcefield_scope={
                "pcff_linear_external": "covered by J17 and graph binding regression",
                "periodic": "precise unsupported diagnostic in current backends",
                "crosslinked": "neutral persistence; parameter coverage remains source-dependent",
            },
        )
    except Exception as error:  # noqa: BLE001 -- retain failed acceptance receipt
        out.update(error_type=type(error).__name__, error=str(error))
    storage.publish(a.output / "outcome.json", storage.json_bytes(out))
    print(json.dumps(out, indent=2))
    return int(not out["passed"])


if __name__ == "__main__":
    raise SystemExit(main())
