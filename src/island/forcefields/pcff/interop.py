"""Optional CAR/MDF interoperability only; never used by native PCFF runtime.

CAR/MDF do not encode all ISLAND semantics. A required, checksummed sidecar
preserves the full original graph and provenance. Loading checks both physical
files against that graph, rather than trusting a matching envelope checksum.
"""

import math
import os
import shutil
import tempfile
from pathlib import Path

from island.workflows import storage
from island.workflows.bundle import system_data, system_from

from .source import FLAGS, boundary, require


def _bytes(system, types, charges):
    system.validate()
    require(
        system.representation == "atomistic" and system.box is None,
        "Finite atomistic CAR/MDF export only",
    )
    ids = sorted(system.topology.sites)
    require(
        set(types) == set(charges) == set(ids),
        "Exact export type/charge coverage required",
    )
    require(
        all(
            type(i) is int
            and type(types[i]) is str
            and types[i]
            and not any(c.isspace() for c in types[i])
            and type(charges[i]) in (float, int)
            and math.isfinite(charges[i])
            for i in ids
        ),
        "Malformed export labels/charges",
    )
    names = {sid: "A" + str(n + 1) for n, sid in enumerate(ids)}
    neighbors = {sid: [] for sid in ids}
    for b in system.topology.bonds.values():
        require(
            b.order in (1, 1.5, 2, 3) and b.aromatic == (b.order == 1.5),
            "Lossy/unsupported bond semantics",
        )
        neighbors[b.site1].append((b.site2, b.order))
        neighbors[b.site2].append((b.site1, b.order))
    car = [
        "!BIOSYM archive 3",
        "PBC=OFF",
        "ISLAND optional final-graph interoperability",
        "!DATE fixed",
    ]
    columns = [
        "element",
        "atom_type",
        "charge_group",
        "isotope",
        "formal_charge",
        "charge",
        "switching_atom",
        "oop_flag",
        "chirality_flag",
        "occupancy",
        "xray_temp_factor",
        "connections",
    ]
    mdf = (
        ["!BIOSYM molecular_data 4", "#topology"]
        + [f"@column {i} {k}" for i, k in enumerate(columns, 1)]
        + ["@molecule chain", ""]
    )
    for i in ids:
        a = system.topology.sites[i]
        require(
            a.element and not any(c.isspace() for c in a.element), "Invalid element"
        )
        xyz = system.coordinates.get(i)
        car.append(
            names[i]
            + " "
            + " ".join(format(v, ".17g") for v in xyz)
            + f" XXXX 1 {types[i]} {a.element} {charges[i]:.17g}"
        )
        isotope = a.metadata.get("isotope", 0)
        require(type(isotope) is int and isotope >= 0, "Malformed isotope")
        mdf.append(
            f"XXXX_1:{names[i]} {a.element} {types[i]} ? {isotope} {a.formal_charge} {charges[i]:.17g} 0 0 8 1.0 0.0 "
            + " ".join(
                names[j] + "/" + format(order, ".17g")
                for j, order in sorted(neighbors[i])
            )
        )
    car += ["end", "end"]
    mdf += ["", "!", "#end"]
    return {
        "system.car": ("\n".join(car) + "\n").encode(),
        "system.mdf": ("\n".join(mdf) + "\n").encode(),
    }, names


@boundary
def export_msi2lmp_car_mdf(system, output_dir, *, types, charges, provenance=None):
    """Exclusively publish checked CAR/MDF + lossless sidecar. No converter invoked.

    Supplied types/charges remain explicitly external unless their native record
    identity is included in provenance. Export does not certify their chemistry.
    """
    files, names = _bytes(system, types, charges)
    payload = {
        "schema": "island_msi2lmp_interop_v1",
        "system": system_data(system),
        "types": dict(types),
        "charges": dict(charges),
        "atom_names": names,
        "provenance": provenance,
        "optional_interoperability_only": True,
        **FLAGS,
    }
    sidecar = storage.json_bytes(storage.encode(payload))
    target = Path(output_dir)
    require(not target.exists(), "Export target exists")
    require(target.parent.is_dir(), "Export parent must exist")
    candidate = Path(tempfile.mkdtemp(prefix=".interop-", dir=target.parent))
    reserved = False
    published = False
    try:
        files["system.json"] = sidecar
        for name, data in files.items():
            storage.publish(candidate / name, data)
        manifest = {
            "schema": "island_msi2lmp_interop_manifest_v1",
            "files": {n: storage.checksum(b) for n, b in files.items()},
        }
        storage.publish(candidate / "manifest.json", storage.json_bytes(manifest))
        load_msi2lmp_car_mdf(candidate)
        target.mkdir()  # Atomic exclusive reservation, never overwrite an old export.
        reserved = True
        os.replace(candidate, target)
        published = True
        storage.sync_directory(target.parent)
    finally:
        if candidate.exists():
            shutil.rmtree(candidate)
        if reserved and not published:
            target.rmdir()  # Remove only our empty reservation, not unexpected contents.
    return load_msi2lmp_car_mdf(target)


@boundary
def load_msi2lmp_car_mdf(directory):
    """Return an owned original system after reversible physical-file checks."""
    root = Path(directory)
    manifest = storage.read_json(root / "manifest.json")
    require(
        set(manifest) == {"schema", "files"}
        and manifest["schema"] == "island_msi2lmp_interop_manifest_v1"
        and set(manifest["files"]) == {"system.car", "system.mdf", "system.json"},
        "Malformed interop manifest",
    )
    for name, sha in manifest["files"].items():
        require(
            storage.checksum(storage.child(root, name).read_bytes()) == sha,
            "Interop checksum mismatch",
        )
    payload = storage.decode(storage.read_json(root / "system.json"))
    require(
        payload["schema"] == "island_msi2lmp_interop_v1"
        and payload["optional_interoperability_only"] is True
        and all(payload[k] == v for k, v in FLAGS.items()),
        "Invalid interop sidecar",
    )
    system = system_from(payload["system"])
    physical, names = _bytes(system, payload["types"], payload["charges"])
    require(names == payload["atom_names"], "Atom ID/name mapping mismatch")
    for name, raw in physical.items():
        require(
            (root / name).read_bytes() == raw,
            "CAR/MDF graph/coordinate/type/charge contradiction",
        )
    return system
