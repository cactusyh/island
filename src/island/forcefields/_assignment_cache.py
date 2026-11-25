"""Opt-in, content-bound cache of successful EMC graph assignments.

Coordinates, site IDs, source files and execution code are part of the key.
A cache is a local acceleration artifact, never an independent validation result.
"""

import json
import os
import platform
import shutil
from dataclasses import asdict
from hashlib import sha256
from pathlib import Path
from tempfile import NamedTemporaryFile

SCHEMA = "island_emc_graph_cache_v1"


def encode(value):
    value = json.loads(json.dumps(value, allow_nan=False))
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def digest(value):
    return sha256(encode(value).encode()).hexdigest()


def file_hash(path):
    with Path(path).open("rb") as handle:
        from hashlib import file_digest

        return file_digest(handle, "sha256").hexdigest()


def tree_hash(path):
    return {
        p.relative_to(path).as_posix(): file_hash(p)
        for p in sorted(Path(path).rglob("*"))
        if p.is_file()
    }


def _code_identity(code):
    # marshal output can change after CPython quickening/interning even when
    # the effective function is unchanged. Hash stable public code attributes.
    from types import CodeType

    def constant(value):
        if isinstance(value, CodeType):
            return _code_identity(value)
        if isinstance(value, (tuple, frozenset)):
            members = [constant(v) for v in value]
            return {
                type(value).__name__: sorted(members, key=repr)
                if isinstance(value, frozenset)
                else members
            }
        return {type(value).__name__: repr(value)}

    return {
        "bytecode": code.co_code.hex(),
        "constants": [constant(c) for c in code.co_consts],
        "names": code.co_names,
        "variables": code.co_varnames,
        "freevars": code.co_freevars,
        "cellvars": code.co_cellvars,
        "args": (code.co_argcount, code.co_posonlyargcount, code.co_kwonlyargcount),
        "flags": code.co_flags,
        "exceptions": code.co_exceptiontable.hex(),
    }


def source_binding(root, binary, directory):
    """Hash current bytes on every call, including custom includes and EMC modules."""

    import numpy
    import rdkit

    from island.forcefields import _emc_graph_input, _emc_isotopes, emc, emc_graph

    here = Path(__file__).parent
    code = [
        here / name
        for name in (
            "_assignment_cache.py",
            "emc.py",
            "emc_graph.py",
            "_emc_graph_input.py",
            "_emc_isotopes.py",
        )
    ]
    return {
        "field_tree": tree_hash(directory),
        "emc_executable": file_hash(binary),
        "emc_scripts": tree_hash(root / "scripts"),
        "emc_lib": tree_hash(root / "lib"),
        "emc_templates": tree_hash(root / "templates"),
        "perl": file_hash(shutil.which("perl")),
        "perl5lib": {
            p: tree_hash(Path(p))
            for p in os.environ.get("PERL5LIB", "").split(os.pathsep)
            if p
        },
        "implementation": {p.name: file_hash(p) for p in code},
        "loaded_implementation": {
            name: digest(_code_identity(function.__code__))
            for name, function in {
                "build": emc_graph._build_script,
                "run": emc._run,
                "map": emc_graph._map_output,
                "assign": emc_graph.assign_emc_graph,
                "input": _emc_graph_input.graph_input,
                "groups": _emc_graph_input.emc_groups,
                "isotopes": _emc_isotopes.split_isotope_types,
            }.items()
        },
        "python": platform.python_version(),
        "rdkit": rdkit.__version__,
        "numpy": numpy.__version__,
    }


def request(source, field, seed, binding):
    from rdkit import Chem

    mol = Chem.Mol(source.mol)
    mol.RemoveAllConformers()
    return {
        "schema": SCHEMA,
        "field": field,
        "seed": seed,
        "binding": binding,
        "molecule_json": json.loads(Chem.MolsToJSON([mol])),
        "sites": source.sites,
        "generated": list(source.generated),
        "xyz": source.xyz.tolist(),
        "bounds": source.bounds.tolist(),
        "coordinate_source": source.coordinate_source,
    }


def read_cache(directory, payload):
    from island.exceptions import ForceFieldError
    from island.forcefields.emc_graph import EMCGraphAssignment

    path = Path(directory) / (digest(payload) + ".json")
    if not path.exists():
        return None
    try:
        envelope = json.loads(path.read_text())
        if set(envelope) != {"request", "result", "sha256"}:
            raise ValueError("invalid envelope")
        if envelope["request"] != payload or envelope["sha256"] != digest(
            {"request": payload, "result": envelope["result"]}
        ):
            raise ValueError("checksum or request mismatch")
        data = envelope["result"]
        data["site_groups"] = {int(k): tuple(v) for k, v in data["site_groups"].items()}
        data["added_hydrogen_ids"] = tuple(data["added_hydrogen_ids"])
        data["input_smiles_site_ids"] = tuple(data["input_smiles_site_ids"])
        result = EMCGraphAssignment(**data)
        mapped = [s for group in result.site_groups.values() for s in group]
        if (
            sorted(mapped) != sorted(payload["sites"])
            or result.field != payload["field"]
        ):
            raise ValueError("mapping or field mismatch")
        return result
    except (ValueError, TypeError, KeyError) as exc:
        raise ForceFieldError(
            f"Invalid assignment cache entry {path.name}: {exc}"
        ) from exc


def write_cache(directory, payload, result):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    record = {"request": payload, "result": asdict(result)}
    envelope = {**record, "sha256": digest(record)}
    path = directory / (digest(payload) + ".json")
    temporary = None
    try:
        with NamedTemporaryFile(
            mode="w", dir=directory, prefix=".cache-", suffix=".tmp", delete=False
        ) as handle:
            temporary = Path(handle.name)
            handle.write(encode(envelope) + "\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
