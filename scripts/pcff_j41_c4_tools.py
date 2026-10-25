"""Build the authorized converter and a diagnostic-only trace copy locally."""

import platform
import shutil
import subprocess
import sys
from pathlib import Path

from scripts.phase_4j41_c2 import digest, write


def command(argv, cwd=None):
    r = subprocess.run(
        argv, cwd=cwd, capture_output=True, text=True, timeout=180, check=False
    )
    return {
        "argv": list(map(str, argv)),
        "exit_code": r.returncode,
        "stdout": r.stdout,
        "stderr": r.stderr,
    }


def build(source, output, lammps):
    output.mkdir(parents=True, exist_ok=False)
    hashes = {}
    for p in sorted(source.iterdir()):
        if p.suffix in (".c", ".h"):
            hashes[p.name] = digest(p)
    records = {}
    for kind in ("original", "trace"):
        dest = output / kind
        dest.mkdir()
        for n in hashes:
            shutil.copyfile(source / n, dest / n)
        if kind == "trace":
            path = dest / "GetParameters.c"
            text = path.read_text()
            start = text.index("int find_match(int n,")
            end = text.index("int match_types(int n,", start)
            part = text[start:end]
            part = part.replace(
                "  int k,match;", "  int k,match;\n  int initial = *backwards_ptr;"
            )
            old = "  if (match) return k;\n  else return -1;"
            trace = """  fprintf(stderr,"C4_MATCH %s %d %d %d %d",item.keyword,n,
          match ? k : -1,*backwards_ptr,initial);
  for (int t=0; t<n; ++t) fprintf(stderr," %s",types1[t]);
  fprintf(stderr,"\\n");
"""
            if old not in part:
                raise ValueError("Unrecognized converter lookup implementation")
            part = part.replace(old, trace + old)
            path.write_text(text[:start] + part + text[end:])
        exe = dest / "msi2lmp"
        argv = ["cc", "-O2", "-o", str(exe), *map(str, sorted(dest.glob("*.c"))), "-lm"]
        record = command(argv)
        records[kind] = {
            **record,
            "source_hashes": {n: digest(dest / n) for n in hashes},
            "executable_sha256": digest(exe) if exe.exists() else None,
        }
        write(output / (kind + "-build.json"), records[kind])
        if record["exit_code"]:
            return {"status": "converter_unavailable", "builds": records}
    git = command(["git", "-C", str(source), "rev-parse", "HEAD"])
    help_run = command([str(lammps), "-h"])
    import numpy
    import openmm
    import rdkit

    metadata = {
        "status": "built",
        "source_hashes": hashes,
        "builds": records,
        "compiler": command(["cc", "--version"]),
        "revision": git,
        "local_source_diff": command(["git", "-C", str(source), "diff", "--", "."]),
        "converter_help": command([str(output / "original/msi2lmp"), "-h"]),
        "lammps": {"path": str(lammps), "sha256": digest(lammps), "help": help_run},
        "python": sys.version,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "numpy": numpy.__version__,
        "rdkit": rdkit.__version__,
        "openmm": openmm.__version__,
        "cpuinfo": Path("/proc/cpuinfo").read_text(),
        "trace_scope": "diagnostic writes only; original converter is independent execution authority",
        "build_flags": ["-O2", "-lm"],
    }
    write(output / "build.json", metadata)
    return metadata


def source_contract(source):
    """Complete raw algorithm files with hashes; no inferred runtime profile."""
    names = (
        "ReadCarFile.c",
        "ReadMdfFile.c",
        "SearchAndFill.c",
        "ReadFrcFile.c",
        "GetParameters.c",
        "InitializeItems.c",
        "MakeLists.c",
        "WriteDataFile.c",
    )
    return {
        "files": {
            n: {"sha256": digest(source / n), "raw": (source / n).read_text()}
            for n in names
        },
        "atom_typing": "ReadCarFile consumes supplied labels; ReadMdfFile consumes connectivity and overwrites charge with MDF precision; no chemical perception",
        "source_precedence": "first matching section; identical ordered key replaces in original table slot only for strictly larger C float version; equal version retains first",
        "ordinary": "direct exact/reversed pass then file-order wildcard pass; on failure one ordinary family map then same search; no recursive mapping",
        "automatic": "ReadFrcFile/InitializeItems do not load auto_equivalence, cff91_auto or bond_increments",
        "role_permutations": "Wilson uses six outer-arm permutations at fixed second-site center without reversal; angle-angle uses three shared-angle slots with two explicit role-preserving permutations; its shared backwards flag changes from -1 to 0 on success, allowing unintended full reversal on later slots; C4 records and rejects role-changing selections",
        "angle_angle_limit": "find_angleangle_data succeeds if any of three kloc entries exists; absent entries retain initialization zero; C4 requires every applicable source term independently",
        "bond13_limit": "GetParameters quo_cp gates bond-bond_1_3 lookup on cp-prefixed types; skipped lookup retains zero and cannot supply a physical missing row",
        "derived": "get_r0/get_t0 reuse slot 0 of already resolved graph bond/angle coefficients; symmetric abbreviated BA/EBT/AT vectors are copied cyclically by SearchAndFill; no missing force constant can be derived from these operations",
        "charges": "converter copies CAR then MDF charge values; ISLAND independent source increment construction remains separate with component conservation, no normalization",
        "style": "Class II; lj/class2/coul/cut finite oracle; sixthpower; special_bonds lj/coul 0 0 1",
        "ignore_authorized": False,
        "runtime_profile_authorized": False,
    }
