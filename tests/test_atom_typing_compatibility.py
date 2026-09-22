import os
import subprocess
import sys

import pytest

pytest.importorskip("rdkit")

from island.chemistry import from_smiles
from island.forcefields.typing import (
    RDKitSmartsAtomTypingEngine,
    island_demo_v1_ruleset,
)


def test_compatibility_ignores_coordinates_and_unrelated_metadata() -> None:
    system = from_smiles("CCO", random_seed=7)
    ruleset = island_demo_v1_ruleset()
    result = RDKitSmartsAtomTypingEngine().type_system(system, ruleset)
    changed = system.copy()
    changed.coordinates.translate([100.0, -3.0, 8.0])
    changed.metadata["coordinate_generation"] = {"method": "replacement"}
    for site in changed.topology.sites.values():
        site.metadata["repeat_unit_index"] = 99
    assert result.is_compatible_with(changed, ruleset)


def test_compatibility_rejects_typing_relevant_chemistry_changes() -> None:
    system = from_smiles("CCO", random_seed=7)
    ruleset = island_demo_v1_ruleset()
    result = RDKitSmartsAtomTypingEngine().type_system(system, ruleset)

    charge_changed = system.copy()
    charge_changed.topology.sites[1].formal_charge = 1
    assert not result.is_compatible_with(charge_changed, ruleset)

    bond_changed = system.copy()
    bond = next(iter(bond_changed.topology.bonds.values()))
    bond_changed.topology.remove_bond(bond.site1, bond.site2)
    bond_changed.topology.add_bond(bond.site1, bond.site2, order=2)
    assert not result.is_compatible_with(bond_changed, ruleset)

    id_changed = system.copy()
    site = id_changed.topology.remove_site(1)
    site.id = 1001
    id_changed.topology.add_site(site)
    assert not result.is_compatible_with(id_changed, ruleset)


def test_result_records_versioned_reproducibility_data() -> None:
    result = RDKitSmartsAtomTypingEngine().type_system(
        from_smiles("CCO", random_seed=7), island_demo_v1_ruleset()
    )
    assert result.ruleset_name == "island_demo_v1"
    assert result.ruleset_version == "1.0.0"
    assert result.engine_name == "rdkit_smarts"
    assert result.engine_version
    assert result.dependency_versions["rdkit"]
    assert len(result.graph_signature) == 64
    assert len(result.ruleset_signature) == 64
    assert len(result.typing_signature) == 64


def test_typing_interfaces_import_without_rdkit() -> None:
    code = r"""
import importlib.abc
import sys

class BlockRDKit(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == "rdkit" or fullname.startswith("rdkit."):
            raise RuntimeError("typing interface eagerly imported RDKit")
        return None

sys.meta_path.insert(0, BlockRDKit())
from island.forcefields import AtomTypingRule, AtomTypingRuleSet
from island.forcefields.typing import AtomTypingEngine, AtomTypingResult
assert "rdkit" not in sys.modules
assert AtomTypingRule.__name__ == "AtomTypingRule"
assert AtomTypingRuleSet.__name__ == "AtomTypingRuleSet"
assert AtomTypingEngine.__name__ == "AtomTypingEngine"
assert AtomTypingResult.__name__ == "AtomTypingResult"
"""
    environment = os.environ.copy()
    environment.pop("PYTHONPATH", None)
    completed = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        check=False,
        env=environment,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr
