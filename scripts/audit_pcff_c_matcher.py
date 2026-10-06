"""Execute pinned msi2lmp leaf matchers as a validation control, never assignment.

The extracted upstream code and compiled library stay in the external output
folder. Applying this matcher to automatic rows is a control, not a capability
claim for the full converter (which does not supplement auto-equivalences).
"""

import argparse
import ctypes
import json
import subprocess
from hashlib import sha256
from pathlib import Path

from pcff_j5_reference import match, read_source

PINS = {
    "GetParameters.c": "add196c115ed9066aa5eba6da5bd8036bd29a3619395385140276e108c769591",
    "Forcefield.h": "4d6c8fa7772d235a3e0fa5ed7398da11b7ced28794170214d754f79d23c3c139",
}
FRC_PIN = "e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c"


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--msi2lmp-source", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    raw = a.source.read_bytes()
    if sha256(raw).hexdigest() != FRC_PIN:
        raise ValueError("Wrong FRC hash")
    for name, expected in PINS.items():
        if sha256((a.msi2lmp_source / name).read_bytes()).hexdigest() != expected:
            raise ValueError(f"Pinned implementation mismatch: {name}")
    a.output.mkdir(parents=True, exist_ok=False)
    original = (a.msi2lmp_source / "GetParameters.c").read_text()
    start, end = (
        original.index("int find_match(int n,"),
        original.index("double get_r0(int typei,"),
    )
    leaf = original[start:end]
    prefix = '#include <string.h>\n#include "Forcefield.h"\nint match_types(int,int,char[][5],char[][5],int*);\n'
    wrapper = """
int select_rows(int n, int count, char query[][5], struct FrcFieldData *data, int *reverse) {
  struct FrcFieldItem item = {0};
  item.entries = count;
  item.data = data;
  return find_match(n, query, item, reverse);
}
"""
    (a.output / "leaf.c").write_text(prefix + leaf + wrapper)
    library = (a.output / "leaf.so").resolve()
    command = [
        "cc",
        "-shared",
        "-fPIC",
        "-Wall",
        "-Werror",
        "-I",
        str(a.msi2lmp_source.resolve()),
        str(a.output / "leaf.c"),
        "-o",
        str(library),
    ]
    compiled = subprocess.run(
        command, text=True, capture_output=True, timeout=60, check=False
    )
    (a.output / "compile.json").write_text(
        json.dumps(
            {
                "command": command,
                "returncode": compiled.returncode,
                "stdout": compiled.stdout,
                "stderr": compiled.stderr,
            },
            indent=2,
        )
    )
    compiled.check_returncode()
    atom = ctypes.c_char * 5

    class Row(ctypes.Structure):
        _fields_ = [
            ("version", ctypes.c_float),
            ("ref", ctypes.c_int),
            ("types", atom * 6),
            ("parameters", ctypes.c_double * 8),
        ]

    lib = ctypes.CDLL(str(library))
    lib.select_rows.argtypes = [
        ctypes.c_int,
        ctypes.c_int,
        ctypes.POINTER(atom),
        ctypes.POINTER(Row),
        ctypes.POINTER(ctypes.c_int),
    ]
    lib.select_rows.restype = ctypes.c_int
    rows = read_source(raw)
    table = rows["quadratic_angle", "cff91_auto"]
    cases = []
    for query in (
        ("h_", "c=_", "s'_"),
        ("h_", "c_", "h_"),
        ("h_", "abs_", "s'_"),
        ("h_", "c'_", "h_"),
    ):
        observations = []
        for backwards in (False, True):
            for reversed_file in (False, True):
                labels = query[::-1] if backwards else query
                selected_table = table[::-1] if reversed_file else table
                data = (Row * len(table))()
                for row, values in zip(data, selected_table, strict=True):
                    for i, label in enumerate(values["types"]):
                        if len(label.encode()) > 4:
                            raise ValueError("Pinned matcher label truncation")
                        row.types[i].value = label.encode()
                    row.version = float(values["version"])
                q = (atom * 3)(*(atom(*label.encode()) for label in labels))
                orientation = ctypes.c_int(0)
                selected = lib.select_rows(
                    3, len(table), q, data, ctypes.byref(orientation)
                )
                record = None if selected == -1 else selected_table[selected]
                observations.append(
                    {
                        "reverse_query": backwards,
                        "reverse_table_control": reversed_file,
                        "source_row": record["id"] if record else None,
                        "values": list(map(str, record["values"])) if record else None,
                        "backwards": orientation.value,
                    }
                )
        independent = match(rows, "quadratic_angle", "cff91_auto", query)
        cases.append(
            {
                "query": query,
                "pinned_c_control": observations,
                "independent_ambiguity_check": independent,
            }
        )
    # These gates establish the matcher semantics only, not charge/model coverage.
    wildcard, _generic, absent, exact = cases
    assert wildcard["independent_ambiguity_check"]["status"] == "ambiguous"
    assert (
        wildcard["pinned_c_control"][0]["values"]
        != wildcard["pinned_c_control"][1]["values"]
    )
    assert exact["independent_ambiguity_check"]["status"] == "assigned"
    assert {o["source_row"] for o in exact["pinned_c_control"]} == {
        "quadratic_angle:cff91_auto:1807"
    }
    assert all(o["source_row"] is None for o in absent["pinned_c_control"])
    receipt = {
        "revision": "e891a3e10973c1a729e391a0aefaa02fd70f8c0f",
        "source_hashes": PINS,
        "frc_sha256": FRC_PIN,
        "line_start": original[:start].count("\n") + 1,
        "line_end": original[:end].count("\n"),
        "library_sha256": sha256(library.read_bytes()).hexdigest(),
        "cases": cases,
        "matcher_control_passed": True,
        "runtime_rule_authorized": False,
        "scope": "Leaf matcher, supplied patterns; not a full converter automatic-equivalence or charge calculation",
    }
    (a.output / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(
        json.dumps({"matcher_control_passed": True, "runtime_rule_authorized": False})
    )


if __name__ == "__main__":
    main()
