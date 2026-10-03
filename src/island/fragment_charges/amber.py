"""Charge-only AmberTools adapter using the existing audited runner and lineage.

No prmtop, parameter import or fabricated ImportedAmberResult is constructed.
"""

import os
import tempfile
from math import isfinite
from pathlib import Path

from island.charge_references.records import pack
from island.workflows import storage
from island.workflows.bundle import system_from

from .core import FragmentChargeResult, boundary, charge_data, need


def _mol2(text, system, names):
    # Offline evidence check. Live generation additionally uses the established
    # RDKit-based full aromatic/bond semantic verifier.
    sections = {}
    current = None
    for line in text.splitlines():
        if line.startswith("@<TRIPOS>"):
            current = line.strip()
            need(current not in sections, "Duplicate MOL2 section")
            sections[current] = []
        elif current and line.strip() and not line.startswith("#"):
            sections[current].append(line.split())
    index = {}
    q = {}
    positions = {}
    for row in sections["@<TRIPOS>ATOM"]:
        need(len(row) >= 9, "Malformed MOL2 atom")
        n = int(row[0])
        name = row[1]
        need(
            name in names and n not in index and names[name] not in q,
            "MOL2 lineage mismatch",
        )
        s = names[name]
        token = row[5].lower()
        element = (
            "Cl"
            if token.startswith("cl")
            else "Br"
            if token.startswith("br")
            else token[0].upper()
        )
        need(element == system.topology.sites[s].element, "MOL2 element mismatch")
        index[n] = s
        q[s] = float(row[8])
        positions[s] = tuple(map(float, row[2:5]))
        need(all(isfinite(v) for v in (*positions[s], q[s])), "Nonfinite MOL2")
    need(set(q) == set(system.topology.sites), "Incomplete MOL2 coverage")
    seen = set()
    aromatic_tokens = {}
    for row in sections["@<TRIPOS>BOND"]:
        key = tuple(sorted((index[int(row[1])], index[int(row[2])])))
        need(
            key in system.topology.bonds and key not in seen,
            "Changed MOL2 connectivity",
        )
        seen.add(key)
        bond = system.topology.bonds[key]
        if bond.aromatic:
            aromatic_tokens[key] = row[3]
        need(
            row[3] in ("1", "2", "ar")
            if bond.aromatic
            else row[3] == str(int(bond.order)),
            "Changed bond semantics",
        )
    need(seen == set(system.topology.bonds), "Incomplete MOL2 bonds")
    # For source aromatic components, accept explicit aromatic bonds or a
    # neutral Kekule representation with the same per-atom pi valence. Do not
    # accept an all-single ring as equivalent aromatic evidence offline.
    remaining = {s for edge in aromatic_tokens for s in edge}
    while remaining:
        component = {min(remaining)}
        pending = list(component)
        while pending:
            site = pending.pop()
            for edge in aromatic_tokens:
                if site in edge:
                    for other in edge:
                        if other not in component:
                            component.add(other)
                            pending.append(other)
        remaining -= component
        tokens = {
            edge: token
            for edge, token in aromatic_tokens.items()
            if edge[0] in component
        }
        if set(tokens.values()) == {"ar"}:
            continue
        need("ar" not in tokens.values(), "Mixed aromatic/Kekule component unsupported")
        for site in component:
            atom = system.topology.sites[site]
            donor = atom.element in {"O", "S"} or (
                atom.element == "N" and len(system.topology.neighbors(site)) == 3
            )
            need(
                sum(token == "2" for edge, token in tokens.items() if site in edge)
                == (0 if donor else 1),
                "Invalid aromatic pi-valence evidence",
            )
    return q, positions


def validate_evidence(fragment, q, backend, evidence):
    from island.forcefields.ambertools.engine import FATAL, SQM_FAILURE, SQM_SUCCESS
    from island.forcefields.ambertools.lineage import generated_atom_name

    system = system_from(fragment["system"])
    need(system.number_of_sites <= 100, "AM1-BCC limit is 100")
    need(
        set(backend) == {"method", "force_field", "settings", "tools", "outcome"}
        and backend["force_field"]["family"] in {"gaff", "gaff2"},
        "AM1-BCC fragment backend is GAFF/GAFF2 only",
    )
    need(
        backend["settings"]["max_atoms"] == 100
        and backend["settings"]["spin_multiplicity"] == 1
        and type(backend["settings"]["timeout_seconds"]) in (float, int)
        and backend["settings"]["timeout_seconds"] > 0
        and isfinite(backend["settings"]["timeout_seconds"]),
        "Invalid charge settings",
    )
    need(
        backend["outcome"] == "sqm_completed_and_typed_lineage_validated",
        "Unsuccessful charge calculation",
    )
    need(
        set(backend["settings"])
        == {"max_atoms", "spin_multiplicity", "timeout_seconds"},
        "Unsupported charge settings",
    )
    need(
        set(backend["force_field"]) == {"family", "data_sha256"},
        "Invalid force field identity",
    )
    need(
        set(evidence)
        == {
            "files",
            "sha256",
            "stage",
            "artifact_dir",
            "names",
            "chemical_graph_identity",
        },
        "Invalid calculation evidence",
    )
    need(type(evidence["artifact_dir"]) is str, "Invalid historical artifact path")
    need(
        type(backend["settings"]["max_atoms"]) is int
        and type(backend["settings"]["spin_multiplicity"]) is int,
        "Integer calculation settings required",
    )
    files = evidence["files"]
    need(set(files) == set(evidence["sha256"]), "Artifact inventory mismatch")
    for name, value in files.items():
        need(
            storage.checksum(value.encode()) == evidence["sha256"][name],
            "Artifact checksum mismatch",
        )
    names = {
        generated_atom_name(system.topology.sites[s].element, k): s
        for k, s in enumerate(sorted(system.topology.sites))
    }
    need(evidence["names"] == names, "Generated-name mapping changed")
    stage = evidence["stage"]
    need(
        set(stage)
        == {
            "stage",
            "command",
            "returncode",
            "elapsed_seconds",
            "stdout_sha256",
            "stderr_sha256",
            "stdout_tail",
            "stderr_tail",
        },
        "Incomplete stage evidence",
    )
    need(
        type(stage["returncode"]) is int
        and type(stage["elapsed_seconds"]) in (int, float)
        and isfinite(stage["elapsed_seconds"])
        and stage["elapsed_seconds"] >= 0,
        "Invalid stage outcome",
    )
    command = stage["command"]
    expected = command_for(
        backend["tools"]["executables"]["antechamber"], backend["force_field"]["family"]
    )
    need(
        stage["stage"] == "antechamber"
        and stage["returncode"] == 0
        and command == expected,
        "Command/method contradiction",
    )
    for stream in ("stdout", "stderr"):
        need(
            stage[stream + "_sha256"]
            == evidence["sha256"]["antechamber." + stream + ".log"],
            "Stage log digest mismatch",
        )
        need(
            not FATAL.search(files["antechamber." + stream + ".log"]), "Fatal stage log"
        )
    text = files["sqm.out"]
    need(
        SQM_SUCCESS.search(text)
        and not SQM_FAILURE.search(text)
        and not FATAL.search(text),
        "SQM completion/convergence unverified",
    )
    actual, _ = _mol2(files["typed.mol2"], system, names)
    need(actual == q, "Raw charges differ from typed artifact")
    _, xyz = _mol2(files["input.mol2"], system, names)
    need(
        all(
            abs(xyz[s][j] - system.coordinates.get(s)[j]) <= 5.000001e-7
            for s in xyz
            for j in range(3)
        ),
        "Input geometry mismatch",
    )
    tools = backend["tools"]
    need(
        set(tools)
        == {
            "executables",
            "executable_sha256",
            "tool_versions",
            "package",
            "force_field_data",
        }
        and set(tools["executables"])
        == set(tools["executable_sha256"])
        == {"antechamber", "parmchk2", "tleap"},
        "Incomplete executable provenance",
    )
    need(
        Path(tools["force_field_data"]["path"]).name
        == backend["force_field"]["family"] + ".dat",
        "Selected force-field path/family contradiction",
    )
    need(
        isinstance(tools["tool_versions"], dict)
        and all(
            type(k) is str and type(v) is str for k, v in tools["tool_versions"].items()
        ),
        "Malformed tool versions",
    )
    need(
        tools["force_field_data"]["sha256"] == backend["force_field"]["data_sha256"],
        "Force field data identity mismatch",
    )
    for h in [
        *tools["executable_sha256"].values(),
        tools["force_field_data"]["sha256"],
    ]:
        need(
            type(h) is str and len(h) == 64 and set(h) <= set("0123456789abcdef"),
            "Invalid tool/data digest",
        )
    need(
        evidence["chemical_graph_identity"] == fragment["data"]["graph_identity"],
        "Charge/chemical identity mismatch",
    )


def command_for(executable, family):
    return [
        executable,
        "-i",
        "input.mol2",
        "-fi",
        "mol2",
        "-o",
        "typed.mol2",
        "-fo",
        "mol2",
        "-at",
        family,
        "-c",
        "bcc",
        "-nc",
        "0",
        "-m",
        "1",
        "-s",
        "2",
        "-j",
        "4",
        "-du",
        "yes",
        "-pf",
        "no",
    ]


class AM1BCCFragmentChargeBackend:
    """Explicit GAFF or GAFF2 charge calculation; no other family fallback."""

    @boundary
    def __init__(self, *, force_field, work_root, amberhome=None, timeout_seconds=600):
        from island.forcefields.ambertools.models import AmberToolsOptions

        need(
            force_field in {"gaff", "gaff2"},
            "AM1-BCC fragment backend requires GAFF/GAFF2",
        )
        self.options = AmberToolsOptions(
            force_field,
            "am1bcc",
            work_root=Path(work_root),
            amberhome=Path(amberhome) if amberhome else None,
            timeout_seconds=timeout_seconds,
            max_atoms=100,
            retain_success_artifacts=True,
        )

    @boundary
    def calculate(self, fragment):
        from island.forcefields.ambertools.engine import (
            _discover,
            _require_artifact,
            _run_stage,
        )
        from island.forcefields.ambertools.lineage import (
            parse_and_validate_mol2,
            prepare_input,
        )

        f = fragment.payload
        system = system_from(f["system"])
        prepared = prepare_input(system, None, max_atoms=100)
        options = self.options
        chain = _discover(options)
        options.work_root.mkdir(parents=True, exist_ok=True)
        directory = Path(
            tempfile.mkdtemp(prefix="island-fragment-", dir=options.work_root)
        )
        (directory / "input.mol2").write_text(prepared.mol2_text)
        (directory / "lineage.json").write_bytes(storage.json_bytes(prepared.names))
        stage = _run_stage(
            "antechamber",
            command_for(chain.executables["antechamber"], options.force_field),
            directory,
            options.timeout_seconds,
            environment={**os.environ, "AMBERHOME": str(chain.amberhome)},
        )
        for name in ("typed.mol2", "sqm.out"):
            _require_artifact("antechamber", directory / name, directory)
        _, q, _ = parse_and_validate_mol2(
            directory / "typed.mol2",
            system,
            prepared.names,
            prepared.expected_cip,
            stage="antechamber",
        )
        files = {
            p.name: p.read_text() for p in sorted(directory.iterdir()) if p.is_file()
        }
        tools = {
            "executables": chain.executables,
            "executable_sha256": {
                k: storage.checksum(Path(v).read_bytes())
                for k, v in chain.executables.items()
            },
            "tool_versions": chain.versions,
            "package": chain.package,
            "force_field_data": {
                "path": str(chain.data_file),
                "sha256": storage.checksum(chain.data_file.read_bytes()),
                "header": chain.data_file.read_text().splitlines()[0],
            },
        }
        p = {
            "schema": "island_fragment_charges_v1",
            "fragment": f,
            "raw_charges": q,
            "backend": {
                "method": "am1bcc",
                "force_field": {
                    "family": options.force_field,
                    "data_sha256": tools["force_field_data"]["sha256"],
                },
                "settings": {
                    "timeout_seconds": options.timeout_seconds,
                    "max_atoms": 100,
                    "spin_multiplicity": 1,
                },
                "tools": tools,
                "outcome": "sqm_completed_and_typed_lineage_validated",
            },
            "evidence": {
                "files": files,
                "sha256": {k: storage.checksum(v.encode()) for k, v in files.items()},
                "stage": stage,
                "artifact_dir": str(directory),
                "names": prepared.names,
                "chemical_graph_identity": f["data"]["graph_identity"],
            },
            "data": {},
        }
        p["data"] = charge_data(p)
        result = FragmentChargeResult(pack(p))
        result.validate_integrity()
        return result
