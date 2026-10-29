"""Independent standard-library graph/charge reconstruction; no ISLAND parser."""

import argparse
import json
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from pathlib import Path

DATA_SHA = "8b03e6e8f3c51525f011d9d64e02d69d81faab78cd17ba66676b208d1c4b9c24"
TYPES = {"H": 98, "O": 122, "Si": 868, "C": 871}
DFT = {
    "O": ["-1.261", "-1.255"],
    "Si": ["2.121", "1.874"],
    "C": ["-1.102", "-1.083"],
    "H": ["0.224", "0.2225"],
}


def pinned(path, pin):
    raw = Path(path).read_bytes()
    if sha256(raw).hexdigest() != pin:
        raise ValueError("independent_source_pin_mismatch:" + str(path))
    return raw


def check(mapping, system, source):
    # Independently decode the data-only envelope; never import production storage.
    encoded = mapping["payload"]
    raw = (
        json.dumps(encoded, sort_keys=True, separators=(",", ":"), allow_nan=False)
        + "\n"
    ).encode()
    if (
        set(mapping) != {"payload", "sha256"}
        or sha256(raw).hexdigest() != mapping["sha256"]
    ):
        raise ValueError("independent_envelope_checksum_mismatch")

    def unpack(value):
        if value is None or type(value) in (str, int, float, bool):
            return value
        if type(value) is not dict or len(value) != 1:
            raise ValueError("independent_malformed_container")
        tag, rows = next(iter(value.items()))
        if tag == "map":
            result = {}
            for key, item in rows:
                if key in result:
                    raise ValueError("independent_duplicate_key")
                result[key] = unpack(item)
            return result
        if tag in ("list", "tuple"):
            return [unpack(v) for v in rows]
        raise ValueError("independent_unknown_container")

    p = unpack(encoded)
    if "payload" in system:
        system = unpack(system["payload"])
    atoms = {r["id"]: r for r in system["sites"]}
    adjacency = {i: set() for i in atoms}
    for bond in system["bonds"]:
        a, b = bond["site1"], bond["site2"]
        adjacency[a].add(b)
        adjacency[b].add(a)
    elements = {i: r["element"] for i, r in atoms.items()}
    ends = {
        i
        for i in atoms
        if elements[i] == "Si"
        and Counter(elements[j] for j in adjacency[i]) == {"C": 3, "O": 1}
    }
    if len(ends) != 2:
        raise ValueError("independent_terminal_connectivity_mismatch")
    eo = p["E_oxygen_site_id"]
    if elements[eo] != "O" or not (adjacency[eo] & ends):
        raise ValueError("independent_E_oxygen_mismatch")
    e = set(ends) | {eo}
    for si in ends:
        for c in adjacency[si]:
            if elements[c] == "C":
                e.add(c)
                e.update(adjacency[c] - {si})
    if len(e) != 27:
        raise ValueError("independent_E_fragment_size_mismatch")
    raw_charges = {t: set() for t in TYPES.values()}
    active = False
    for line in source.decode().splitlines():
        text = line.split("#", 1)[0].strip()
        if not text:
            continue
        if text == "Atoms":
            active = True
            continue
        if active and text[0].isalpha():
            break
        if active:
            cols = text.split()
            raw_charges[int(cols[2])].add(cols[3])
    total = Fraction()
    if len(p["mapping"]) != len(atoms) or {r["site_id"] for r in p["mapping"]} != set(
        atoms
    ):
        raise ValueError("independent_mapping_coverage_mismatch")
    if p["atom_order"] != [r["id"] for r in system["sites"]]:
        raise ValueError("independent_atom_order_mismatch")
    for row in p["mapping"]:
        i = row["site_id"]
        element = elements[i]
        fragment = "E" if i in e else "M"
        printed = DFT[element][int(i in e)]
        literal = format(float(printed), ".16e")
        if (
            row["fragment"] != fragment
            or row["element"] != element
            or row["author_type_id"] != TYPES[element]
            or row["printed_Figure7_DFT_charge"] != printed
            or row["exact_author_serialized_charge"] != literal
            or literal not in raw_charges[TYPES[element]]
        ):
            raise ValueError("independent_charge_or_type_mismatch")
        total += Fraction(literal)
    if total != Fraction(p["charge_policy"]["exact_serialized_sum_fraction"]):
        raise ValueError("independent_charge_sum_mismatch")
    if p["execution_profile"] is not None or p["execution_authorized"]:
        raise ValueError("independent_forged_execution_authority")
    return {
        "atoms_checked": len(atoms),
        "E_atoms_checked": len(e),
        "exact_charge_sum_fraction": str(total),
        "author_source_sha256": DATA_SHA,
        "independent_row_mapping_equal": True,
        "numerical_PDMS_execution": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--mapping-sha256", required=True)
    parser.add_argument("--system", type=Path, required=True)
    parser.add_argument("--system-sha256", required=True)
    parser.add_argument("--author-data", type=Path, required=True)
    args = parser.parse_args()
    result = check(
        json.loads(pinned(args.mapping, args.mapping_sha256)),
        json.loads(pinned(args.system, args.system_sha256)),
        pinned(args.author_data, DATA_SHA),
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
