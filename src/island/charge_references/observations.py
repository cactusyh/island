"""Data-only post-BCC observations, including historically failed imports.

An observation is deliberately not a ChargeReference or parameterized system.
Artifact-backed creation reads exact manifest paths. Offline validation rechecks
embedded evidence, without optional scientific libraries or executable calls.
"""

import re
from dataclasses import dataclass
from math import fsum, isfinite
from pathlib import Path, PurePosixPath

from island.core.coordinate_provenance import coordinate_hash
from island.dynamics._checkpoint_data import strict_load
from island.exceptions import ChargeReferenceError
from island.forcefields.ambertools.lineage import generated_atom_name
from island.forcefields.ambertools.models import SQM_SUCCESS
from island.workflows import storage
from island.workflows.bundle import system_from

from .correspondence import repeat_correspondence, require
from .records import ChargeReference, pack, unpack

OBSERVATION_SCHEMA = "island_post_bcc_observation_v1"


def restricted_mol2(text, system, names):
    """Strict H/C/O single-bond reader for the already verified PE/PEO scope."""
    sections = {}
    current = None
    for line in text.splitlines():
        if line.startswith("@<TRIPOS>"):
            current = line[9:].strip()
            require(current not in sections, "Duplicate MOL2 section")
            sections[current] = []
        elif current and line.strip():
            sections[current].append(line.split())
    indices = {}
    charges = {}
    positions = {}
    for row in sections["ATOM"]:
        index = int(row[0])
        name = row[1]
        require(
            index > 0 and index not in indices and name in names,
            "Invalid MOL2 index/name",
        )
        site = names[name]
        require(site not in charges, "Duplicate MOL2 atom")
        require(
            row[5][0].upper() == system.topology.sites[site].element,
            "MOL2 element mismatch",
        )
        q = float(row[8])
        xyz = tuple(map(float, row[2:5]))
        require(
            len(xyz) == 3 and all(isfinite(v) for v in (*xyz, q)), "Nonfinite MOL2 data"
        )
        indices[index] = site
        charges[site] = q
        positions[site] = xyz
    require(set(charges) == set(system.topology.sites), "Incomplete MOL2 coverage")
    edges = set()
    bond_ids = set()
    for row in sections["BOND"]:
        i = int(row[0])
        edge = tuple(sorted((indices[int(row[1])], indices[int(row[2])])))
        require(
            i > 0 and i not in bond_ids and edge not in edges and row[3] == "1",
            "Invalid/duplicate MOL2 bond",
        )
        bond_ids.add(i)
        edges.add(edge)
    require(edges == set(system.topology.bonds), "MOL2 connectivity mismatch")
    counts = sections["MOLECULE"][1]
    require(
        int(counts[0]) == len(charges) and int(counts[1]) == len(edges),
        "MOL2 declared inventory mismatch",
    )
    return indices, charges, positions


def prmtop_charge_data(text, system, names):
    """Read charge evidence only, not a second parameter/import implementation."""
    from decimal import Decimal

    def field(name, width):
        sections = re.split(r"%FLAG " + name + r"\s*\n", text)
        require(len(sections) == 2, "Missing/duplicate prmtop field")
        body = sections[1].split("\n", 1)[1].split("%FLAG")[0]
        return [
            line[i : i + width].strip()
            for line in body.splitlines()
            for i in range(0, len(line), width)
            if line[i : i + width].strip()
        ]

    atom_names = field("ATOM_NAME", 4)
    numbers = field("ATOMIC_NUMBER", 8)
    values = field("CHARGE", 16)
    require(
        len(atom_names) == len(names)
        and len(set(atom_names)) == len(names)
        and set(atom_names) == set(names),
        "prmtop atom inventory mismatch",
    )
    charges = {}
    for name, number, token in zip(atom_names, numbers, values, strict=True):
        site = names[name]
        require(
            int(number) == system.topology.sites[site].atomic_number,
            "prmtop element mismatch",
        )
        q = float(Decimal(token) / Decimal("18.2223"))
        require(isfinite(q), "Nonfinite prmtop charge")
        charges[site] = q
    return charges


def _observation_data(evidence):
    require(
        set(evidence) == {"manifest_json", "case", "files", "historical_reference"},
        "Unexpected observation evidence",
    )
    manifest = strict_load(evidence["manifest_json"])
    rows = [r for r in manifest["cases"] if r["case"] == evidence["case"]]
    require(len(rows) == 1, "Missing/ambiguous manifest case")
    row = rows[0]
    checks = row["source_checksums"]
    files = evidence["files"]
    require(
        set(files) == set(checks), "Evidence differs from exact manifest file inventory"
    )
    paths = {}
    for path, expected in checks.items():
        parts = PurePosixPath(path)
        require(
            not parts.is_absolute() and ".." not in parts.parts and str(parts) == path,
            "Unsafe manifest path",
        )
        require(
            storage.checksum(files[path].encode()) == expected,
            "Artifact checksum mismatch",
        )
        require(parts.name not in paths, "Ambiguous artifact basename/directory")
        paths[parts.name] = path
    backend = {
        str(PurePosixPath(p).parent)
        for n, p in paths.items()
        if n not in {"input-system.json", "preparation.json"}
    }
    require(
        len(backend) == 1 and "." not in backend, "Ambiguous backend artifact directory"
    )

    def text(name):
        return files[paths[name]]

    source = storage.decode(strict_load(text("input-system.json")))
    system = system_from(source)
    correspondence = repeat_correspondence(system)
    require(
        correspondence["compatible"] and 0 < system.number_of_sites <= 100,
        "Unsupported reference scope",
    )
    require(
        row["psmiles"] == correspondence["definition"]
        and type(row["dp"]) is int
        and row["dp"] == correspondence["dp"],
        "Case chemical identity mismatch",
    )
    require(
        type(row["sites"]) is int and row["sites"] == system.number_of_sites,
        "Manifest site count mismatch",
    )
    generation = system.metadata["polymer"]["coordinate_generation"]
    seeds = {k: generation[k] for k in ("template_seed", "assembly_seed")}
    require(
        type(row["seed"]) is int and all(v == row["seed"] for v in seeds.values()),
        "Case seed mismatch",
    )
    names = {
        generated_atom_name(system.topology.sites[s].element, i): s
        for i, s in enumerate(sorted(system.topology.sites))
    }
    lineage = strict_load(text("lineage.json"))
    require(
        lineage.get("schema") == "island_ambertools_lineage_v1"
        and set(lineage)
        == {"schema", "input_name_to_site_id", "input_index_to_site_id"},
        "Unsupported lineage schema",
    )
    require(
        all(
            type(v) is int
            for key in ("input_name_to_site_id", "input_index_to_site_id")
            for v in lineage[key].values()
        ),
        "Noninteger lineage site ID",
    )
    require(lineage["input_name_to_site_id"] == names, "Source name lineage mismatch")
    indices, _, input_xyz = restricted_mol2(text("input.mol2"), system, names)
    require(
        lineage["input_index_to_site_id"] == {str(i): s for i, s in indices.items()},
        "Source index lineage mismatch",
    )
    require(
        all(
            abs(input_xyz[s][j] - system.coordinates.get(s)[j]) <= 5.000001e-7
            for s in input_xyz
            for j in range(3)
        ),
        "Input coordinates disagree with source serialization",
    )
    _, charges, xyz = restricted_mol2(text("typed.mol2"), system, names)
    bcc = {}
    for line in text("ANTECHAMBER_AM1BCC.AC").splitlines():
        if not line.startswith("ATOM "):
            continue
        fields = line.split()
        name = fields[2][:4]
        require(
            fields[2] == name + "MOL" and name in names and names[name] not in bcc,
            "Invalid BCC atom lineage",
        )
        q = float(fields[-2])
        require(isfinite(q), "Nonfinite BCC charge")
        bcc[names[name]] = q
    require(bcc == charges, "Typed MOL2 differs from post-BCC charges")
    # SQM input carries named identities; its final Cartesian table carries indices.
    require(
        re.search(r"qm_theory\s*=\s*'AM1'", text("sqm.in"))
        and re.search(r"qmcharge\s*=\s*0\s*[,/]", text("sqm.in")),
        "Unsupported SQM method/charge",
    )
    require(
        " -p gaff2" in text("antechamber.stdout.log")
        and "/sqm -O -i sqm.in -o sqm.out" in text("antechamber.stdout.log"),
        "Missing GAFF2/SQM execution evidence",
    )
    require(
        re.search(
            r'^source\s+"[^"\n]*/leaprc\.gaff2"\s*$', text("leap.in"), re.MULTILINE
        ),
        "Unsupported force field evidence",
    )
    # The retained command path writes/reloads three-decimal AC coordinates
    # before writing sqm.in. Compare the actual intermediate, not input MOL2
    # directly to a presumed four-decimal rounding of the original frame.
    ac_xyz = {}
    for line in text("ANTECHAMBER_AC.AC").splitlines():
        if not line.startswith("ATOM "):
            continue
        fields = line.split()
        name = fields[2][:4]
        require(
            fields[2] == name + "MOL" and name in names and names[name] not in ac_xyz,
            "Invalid pre-QM AC lineage",
        )
        xyz_ac = tuple(map(float, fields[4:7]))
        require(
            len(xyz_ac) == 3 and all(isfinite(v) for v in xyz_ac),
            "Invalid pre-QM AC coordinates",
        )
        site = names[name]
        require(
            all(abs(xyz_ac[j] - input_xyz[site][j]) <= 5.000001e-4 for j in range(3)),
            "Pre-QM AC differs from three-decimal input serialization",
        )
        ac_xyz[site] = xyz_ac
    require(set(ac_xyz) == set(charges), "Incomplete pre-QM AC geometry")
    sqm_input = [
        r.split()
        for r in text("sqm.in").splitlines()
        if re.match(r"^\s*\d+\s+\w+\s", r)
    ]
    require(
        len(sqm_input) == len(names) and {r[1] for r in sqm_input} == set(names),
        "SQM input coverage mismatch",
    )
    for r in sqm_input:
        require(
            int(r[0]) == system.topology.sites[names[r[1]]].atomic_number,
            "SQM element mismatch",
        )
        sqm_start = tuple(map(float, r[2:5]))
        require(
            len(sqm_start) == 3
            and all(isfinite(v) for v in sqm_start)
            and all(sqm_start[j] == ac_xyz[names[r[1]]][j] for j in range(3)),
            "SQM input coordinates disagree with retained pre-QM AC",
        )
    output = text("sqm.out")
    from island.forcefields.ambertools.engine import FATAL, NONZERO_ERRORS, SQM_FAILURE

    require(
        SQM_SUCCESS.search(output)
        and not FATAL.search(output)
        and not SQM_FAILURE.search(output),
        "SQM completion/convergence unverified",
    )
    # Last full Cartesian table; indices and elements must resolve via sqm.in.
    sqm_xyz = {}
    for line in output.rsplit("Final Structure", 1)[-1].splitlines():
        r = line.split()
        if len(r) == 7 and r[0] == "QMMM:" and r[1].isdigit() and r[2].isdigit():
            i = int(r[1])
            require(
                1 <= i <= len(sqm_input) and int(r[2]) == i, "SQM output index mismatch"
            )
            s = names[sqm_input[i - 1][1]]
            require(
                s not in sqm_xyz and r[3] == system.topology.sites[s].element,
                "SQM output lineage mismatch",
            )
            value = tuple(map(float, r[4:]))
            require(all(isfinite(v) for v in value), "Invalid SQM geometry")
            sqm_xyz[s] = value
    require(set(sqm_xyz) == set(system.topology.sites), "Missing final SQM geometry")
    for stage in ("antechamber", "parmchk2", "tleap"):
        logs = text(stage + ".stdout.log") + "\n" + text(stage + ".stderr.log")
        require(
            not FATAL.search(logs) and not NONZERO_ERRORS.search(logs),
            "Failed stage diagnostic",
        )
    require(
        re.search(r"Exiting LEaP: Errors = 0;", text("leap.log")) is not None,
        "Unverified tleap completion",
    )
    require(
        "am1bcc -i ANTECHAMBER_AM1BCC_PRE.AC -o ANTECHAMBER_AM1BCC.AC"
        in text("antechamber.stdout.log"),
        "Missing recorded BCC execution",
    )
    status = row["status"]
    require(status in {"passed", "failed"}, "Invalid historical status")
    historical = evidence["historical_reference"]
    if status == "passed":
        ref = ChargeReference(pack(historical))
        rp = ref.payload
        require(
            ref.identity == row["reference_identity"],
            "Original reference identity mismatch",
        )
        require(rp["system"] == source, "Original reference source mismatch")
        prep = storage.decode(strict_load(text("preparation.json")))
        require(
            {k: prep[k] for k in ("record", "record_signature")} == rp["preparation"]
            and prep["import_provenance"] == rp["import_content"]["provenance"],
            "Original preparation mismatch",
        )
        require(
            rp["preparation"]["record_signature"] == row["record_signature"]
            and rp["preparation"]["record"]["imported_result_signature"]
            == row["import_signature"],
            "Historical signature mismatch",
        )
        for name, expected in rp["preparation"]["record"]["artifact_sha256"].items():
            require(checks[paths[name]] == expected, "Reference artifact mismatch")
        require(
            rp["import_content"]["source_sha256"] == checks[paths["result.prmtop"]],
            "Reference prmtop source mismatch",
        )
    else:
        require(
            historical is None and "preparation.json" not in paths,
            "Failed observation cannot claim successful preparation",
        )
        require(
            type(row.get("failure")) is str
            and "import: Charge assignment is incomplete" in row["failure"],
            "Failure is not the retained charge-import rejection",
        )
    # prmtop is retained separately; it is not the raw projection input.
    prmtop = prmtop_charge_data(text("result.prmtop"), system, names)
    if status == "passed":
        require(
            all(
                abs(prmtop[s] - rp["charge_result"]["assignments"][s]["charge"]) < 1e-14
                for s in prmtop
            ),
            "Original imported charges disagree with prmtop",
        )
    return {
        "system": source,
        "system_sha256": storage.checksum(storage.json_bytes(storage.encode(source))),
        "correspondence": correspondence,
        "seeds": seeds,
        "case": row["case"],
        "historical_status": status,
        "historical_failure": row.get("failure"),
        "historical_reference_identity": row.get("reference_identity"),
        "preparation_signature": row.get("record_signature"),
        "import_signature": row.get("import_signature"),
        "source_manifest_sha256": storage.checksum(evidence["manifest_json"].encode()),
        "artifact_paths": paths,
        "artifact_sha256": checks,
        "raw_stage": "typed.mol2 post-BCC",
        "prmtop_charges": prmtop,
        "prmtop_total_e": fsum(prmtop.values()),
        "prmtop_minus_raw_e": {s: prmtop[s] - charges[s] for s in sorted(charges)},
        "raw_charges": charges,
        "raw_total_e": fsum(charges[s] for s in sorted(charges)),
        "formal_charge_e": 0,
        "raw_residual_e": fsum(charges[s] for s in sorted(charges)),
        "atom_count": len(charges),
        "input_coordinate_fingerprint": coordinate_hash(
            {s: tuple(system.coordinates.get(s)) for s in sorted(charges)}
        ),
        "post_sqm_coordinates_angstrom": sqm_xyz,
        "post_sqm_coordinate_fingerprint": coordinate_hash(sqm_xyz),
        "typed_coordinates_angstrom": xyz,
        "typed_coordinate_fingerprint": coordinate_hash(xyz),
        "stage_outcome_evidence": "SQM completion, fatal-diagnostic checks, BCC command and LEaP zero errors; failed imports have no signed stage return codes",
        "unit": "elementary_charge",
        "production_validated": False,
        "simulation_readiness": "not_established",
    }


@dataclass(frozen=True)
class RawChargeObservation:
    json_text: str

    def validate_integrity(self):
        try:
            p = unpack(self.json_text)
            require(
                set(p) == {"schema", "evidence", "data"}
                and p["schema"] == OBSERVATION_SCHEMA,
                "Invalid observation schema",
            )
            require(
                pack(p["data"]) == pack(_observation_data(p["evidence"])),
                "Observation differs from source evidence",
            )
        except ChargeReferenceError:
            raise
        except Exception as error:
            raise ChargeReferenceError(f"Invalid raw observation: {error}") from error

    @property
    def payload(self):
        self.validate_integrity()
        return unpack(self.json_text)

    @property
    def identity(self):
        return storage.checksum(storage.json_bytes(storage.encode(self.payload)))


def observe_retained_case(root, manifest_path, case, *, historical_reference=None):
    """Verify all exact listed artifacts before parsing scientific contents."""
    try:
        root = Path(root)
        manifest_json = Path(manifest_path).read_text()
        manifest = strict_load(manifest_json)
        rows = [r for r in manifest["cases"] if r["case"] == case]
        require(len(rows) == 1, "Missing/ambiguous manifest case")
        directory = storage.child(root, case)
        checks = rows[0]["source_checksums"]
        files = {}
        for name, expected in checks.items():
            raw = storage.child(directory, name).read_bytes()
            require(storage.checksum(raw) == expected, "Artifact checksum mismatch")
            files[name] = raw.decode("utf-8")
        listed = {
            PurePosixPath(p).parts[0] for p in checks if len(PurePosixPath(p).parts) > 1
        }
        actual = {
            p.name
            for p in directory.iterdir()
            if p.is_dir() and p.name.startswith("island-ambertools-")
        }
        require(actual == listed, "Ambiguous/unlisted artifact directory")
        if rows[0]["status"] == "passed" and historical_reference is None:
            # Reconstruct only historical successes through the existing importer.
            # Every input path is the unique exact manifest entry, never a glob.
            by_name = {PurePosixPath(p).name: p for p in checks}
            require(len(by_name) == len(checks), "Ambiguous artifact paths")
            from island.workflows.bundle import preparation_from

            from .records import create_charge_reference

            system = system_from(
                storage.decode(strict_load(files[by_name["input-system.json"]]))
            )
            prep_data = storage.decode(strict_load(files[by_name["preparation.json"]]))
            prep = preparation_from(
                system, prep_data, storage.child(directory, by_name["result.prmtop"])
            )
            historical_reference = create_charge_reference(system, prep)
        evidence = {
            "manifest_json": manifest_json,
            "case": case,
            "files": files,
            "historical_reference": historical_reference.payload
            if historical_reference is not None
            else None,
        }
        result = RawChargeObservation(
            pack(
                {
                    "schema": OBSERVATION_SCHEMA,
                    "evidence": evidence,
                    "data": _observation_data(evidence),
                }
            )
        )
        result.validate_integrity()
        return result
    except ChargeReferenceError:
        raise
    except Exception as error:
        raise ChargeReferenceError(f"Cannot observe retained case: {error}") from error


def load_raw_observation(path):
    try:
        result = RawChargeObservation(Path(path).read_text())
        result.validate_integrity()
        return result
    except Exception as error:
        raise ChargeReferenceError(f"Cannot load raw observation: {error}") from error
