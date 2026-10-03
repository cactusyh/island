"""Pinned upstream selection witness; never used as ISLAND's evaluated model."""

import inspect
import json
from hashlib import sha256
from math import isclose, sqrt
from pathlib import Path
from tempfile import TemporaryDirectory

from .adapter import check_foyer_installation
from .source import require


def upstream_system(system, source, types):
    check_foyer_installation()
    import openmm
    from foyer import Forcefield
    from openmm import app

    pin = json.loads(Path(__file__).with_name("resolution_pin.json").read_text())
    methods = {
        "bond": app.forcefield.HarmonicBondGenerator.createForce,
        "angle": app.forcefield.HarmonicAngleGenerator.postprocessSystem,
        "rb": app.forcefield.RBTorsionGenerator.createForce,
    }
    require(
        openmm.__version__ == pin["openmm"]
        and all(
            sha256(inspect.getsource(fn).encode()).hexdigest()
            == pin["matching_methods"][key]
            for key, fn in methods.items()
        ),
        "Unverified OpenMM matching implementation",
    )

    topology = app.Topology()
    residue = topology.addResidue("M", topology.addChain())
    atoms = {}
    for sid in sorted(system.topology.sites):
        a = system.topology.sites[sid]
        atoms[sid] = topology.addAtom(
            a.element,
            app.element.Element.getByAtomicNumber(a.atomic_number),
            residue,
            id=types[sid],
        )
    for a, b in sorted(system.topology.bonds):
        topology.addBond(atoms[a], atoms[b])
    with TemporaryDirectory(prefix="island-opls-source-") as directory:
        file = Path(directory) / "oplsaa.xml"
        file.write_bytes(source.xml)
        ff = Forcefield(forcefield_files=str(file))
        model = ff.createSystem(
            topology,
            nonbondedMethod=app.NoCutoff,
            constraints=None,
            rigidWater=False,
            removeCMMotion=False,
        )
    return model


def verify_resolved(system, source, types, resolved):
    import openmm as mm
    from openmm import unit

    model = upstream_system(system, source, types)
    ids = sorted(system.topology.sites)
    require(
        model.getNumConstraints() == 0
        and not any(model.isVirtualSite(i) for i in range(len(ids))),
        "Unexpected constraints/virtual sites",
    )
    units = {
        "bonds": (unit.nanometer, unit.kilojoule_per_mole / unit.nanometer**2),
        "angles": (unit.radian, unit.kilojoule_per_mole / unit.radian**2),
        "rb": (unit.kilojoule_per_mole,) * 6,
    }
    expected_names = {
        mm.HarmonicBondForce: ("bonds", "getNumBonds", "getBondParameters", 2),
        mm.HarmonicAngleForce: ("angles", "getNumAngles", "getAngleParameters", 3),
        mm.RBTorsionForce: ("rb", "getNumTorsions", "getTorsionParameters", 4),
    }
    seen = set()
    for force in model.getForces():
        if type(force) in expected_names:
            kind, count, get, n = expected_names[type(force)]
            require(kind not in seen, "Duplicate upstream force family")
            seen.add(kind)
            actual = {}
            for i in range(getattr(force, count)()):
                row = getattr(force, get)(i)
                key = tuple(ids[int(j)] for j in row[:n])
                key = min(key, key[::-1])
                require(key not in actual, "Duplicated upstream interaction")
                actual[key] = [
                    x.value_in_unit(u)
                    for x, u in zip(row[n:], units[kind], strict=True)
                ]
            expected = {}
            for row in resolved[kind]:
                raw = row["original"]
                if kind in ("bonds", "angles") and raw["k"] == 0:
                    continue  # Explicit zero source record retained; upstream omits force.
                values = [
                    raw[k]
                    for k in (
                        ("length", "k")
                        if kind == "bonds"
                        else ("angle", "k")
                        if kind == "angles"
                        else tuple(f"c{i}" for i in range(6))
                    )
                ]
                expected[min(row["sites"], row["sites"][::-1])] = values
            require(set(actual) == set(expected), f"Upstream {kind} coverage mismatch")
            require(
                all(
                    all(
                        isclose(a, b, rel_tol=1e-13, abs_tol=1e-13)
                        for a, b in zip(actual[k], expected[k], strict=True)
                    )
                    for k in expected
                ),
                f"Upstream {kind} selection mismatch",
            )
        elif type(force) is mm.NonbondedForce:
            require("nonbonded" not in seen, "Duplicate nonbonded force")
            seen.add("nonbonded")
            require(
                force.getNonbondedMethod() == mm.NonbondedForce.NoCutoff,
                "Unexpected periodic/cutoff model",
            )
            for i, sid in enumerate(ids):
                q, sigma, epsilon = force.getParticleParameters(i)
                row = resolved["sites"][sid]
                actual = (
                    q.value_in_unit(unit.elementary_charge),
                    sigma.value_in_unit(unit.nanometer),
                    epsilon.value_in_unit(unit.kilojoule_per_mole),
                )
                expected = (
                    row["charge_e"],
                    row["sigma_angstrom"] / 10,
                    row["epsilon_kj_mol"],
                )
                require(
                    all(
                        isclose(a, b, rel_tol=1e-13, abs_tol=1e-13)
                        for a, b in zip(actual, expected, strict=True)
                    ),
                    "Upstream LJ/charge mismatch",
                )
            exceptions = {}
            for i in range(force.getNumExceptions()):
                a, b, q, sigma, epsilon = force.getExceptionParameters(i)
                key = tuple(sorted((ids[int(a)], ids[int(b)])))
                require(key not in exceptions, "Duplicate upstream pair")
                exceptions[key] = (
                    q.value_in_unit(unit.elementary_charge**2),
                    epsilon.value_in_unit(unit.kilojoule_per_mole),
                )
            require(
                set(exceptions) == {p["sites"] for p in resolved["pairs"]},
                "Upstream pair inventory mismatch",
            )
            for pair in resolved["pairs"]:
                a, b = (resolved["sites"][i] for i in pair["sites"])
                q, eps = exceptions[pair["sites"]]
                require(
                    isclose(
                        q,
                        a["charge_e"] * b["charge_e"] * pair["coulomb_scale"],
                        abs_tol=1e-13,
                    )
                    and isclose(
                        eps,
                        sqrt(a["epsilon_kj_mol"] * b["epsilon_kj_mol"])
                        * pair["lj_scale"],
                        abs_tol=1e-13,
                    ),
                    "Upstream pair scaling mismatch",
                )
        else:
            raise ValueError(f"Unsupported upstream force {type(force).__name__}")
    require(
        seen == {"bonds", "angles", "rb", "nonbonded"},
        "Incomplete upstream force inventory",
    )
