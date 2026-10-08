"""Periodic OpenMM forces built only from validated numerical records."""

from math import sqrt


def build(data, config, mm):
    model = mm.System()
    sites = {s["id"]: s for s in data["sites"]}
    ids = sorted(sites)
    index = {sid: i for i, sid in enumerate(ids)}
    for sid in ids:
        model.addParticle(sites[sid]["mass"])
    lengths = config["box_lengths"]
    model.setDefaultPeriodicBoxVectors(
        mm.Vec3(lengths[0] / 10, 0, 0),
        mm.Vec3(0, lengths[1] / 10, 0),
        mm.Vec3(0, 0, lengths[2] / 10),
    )
    pcff = data["family"] == "PCFF"
    if pcff:
        from island.evaluation.pcff import expression

        families = sorted({r["family"] for r in data["pcff_terms"]})
        for group, family in enumerate(families):
            rows = [r for r in data["pcff_terms"] if r["family"] == family]
            row = rows[0]
            force = mm.CustomCompoundBondForce(len(row["sites"]), expression(family))
            for i in range(len(row["coefficients"])):
                force.addPerBondParameter(f"c{i}")
            for i in range(len(row["equilibria"])):
                force.addPerBondParameter(f"e{i}")
            for row in rows:
                force.addBond(
                    [index[i] for i in row["sites"]],
                    row["coefficients"] + row["equilibria"],
                )
            force.setUsesPeriodicBoundaryConditions(True)
            force.setName(family)
            force.setForceGroup(group)
            model.addForce(force)
    else:
        bond, angle = mm.HarmonicBondForce(), mm.HarmonicAngleForce()
        for row in data["bonds"]:
            bond.addBond(
                *(index[i] for i in row["sites"]), row["length"] / 10, row["k"] * 100
            )
        for row in data["angles"]:
            angle.addAngle(*(index[i] for i in row["sites"]), row["theta"], row["k"])
        forces = [("bond", bond), ("angle", angle)]
        if data["family"] == "OPLS-AA":
            proper = mm.RBTorsionForce()
            for row in data["torsions"]:
                proper.addTorsion(*(index[i] for i in row["sites"]), *row["rb"])
            forces.append(("rb_proper", proper))
        else:
            for name, key in (("proper", "torsions"), ("improper", "impropers")):
                force = mm.PeriodicTorsionForce()
                for row in data[key]:
                    for k, n, phase in row["fourier"]:
                        force.addTorsion(*(index[i] for i in row["sites"]), n, phase, k)
                forces.append((name, force))
        for group, (name, force) in enumerate(forces):
            force.setName(name)
            force.setForceGroup(group)
            force.setUsesPeriodicBoundaryConditions(True)
            model.addForce(force)
    coul = mm.NonbondedForce()
    coul.setName("coulomb")
    coul.setForceGroup(30)
    coul.setNonbondedMethod(
        mm.NonbondedForce.PME
        if config["electrostatics_method"] == "pme"
        else mm.NonbondedForce.Ewald
    )
    coul.setCutoffDistance(config["nonbonded_cutoff"] / 10)
    coul.setEwaldErrorTolerance(config["pme_tolerance"])
    coul.setUseDispersionCorrection(False)
    coul.setUseSwitchingFunction(False)
    coul.setExceptionsUsePeriodicBoundaryConditions(True)
    if pcff:
        expression = "eps*(2*(rmin/r)^9-3*(rmin/r)^6);eps=2*sqrt(epsilon1*epsilon2)*sigma1^3*sigma2^3/(sigma1^6+sigma2^6);rmin=((sigma1^6+sigma2^6)/2)^(1/6)"
        scaled_expression = "epsilon*(2*(sigma/r)^9-3*(sigma/r)^6)"
    else:
        mixing = (
            "sqrt(sigma1*sigma2)"
            if data["family"] == "OPLS-AA"
            else "(sigma1+sigma2)/2"
        )
        expression = f"4*sqrt(epsilon1*epsilon2)*(x^12-x^6);x=({mixing})/r"
        scaled_expression = "4*epsilon*((sigma/r)^12-(sigma/r)^6)"
    lj = mm.CustomNonbondedForce(expression)
    for name in ("sigma", "epsilon"):
        lj.addPerParticleParameter(name)
    lj.setName("lj")
    lj.setForceGroup(31)
    lj.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
    lj.setCutoffDistance(config["nonbonded_cutoff"] / 10)
    lj.setUseSwitchingFunction(False)
    lj.setUseLongRangeCorrection(False)
    scaled = mm.CustomBondForce(scaled_expression)
    scaled.addPerBondParameter("sigma")
    scaled.addPerBondParameter("epsilon")
    scaled.setName("lj14")
    scaled.setForceGroup(31)
    scaled.setUsesPeriodicBoundaryConditions(True)
    for sid in ids:
        row = sites[sid]
        coul.addParticle(row["charge"], 1, 0)
        lj.addParticle([row["sigma"] / 10, row["epsilon"]])
    for pair in data["pairs"]:
        a, b = pair["sites"]
        x, y = sites[a], sites[b]
        coul.addException(
            index[a], index[b], x["charge"] * y["charge"] * pair["coulomb"], 1, 0
        )
        lj.addExclusion(index[a], index[b])
        if pair["lj"]:
            s1, s2 = x["sigma"], y["sigma"]
            eps = sqrt(x["epsilon"] * y["epsilon"])
            if pcff:
                mean = (s1**6 + s2**6) / 2
                sigma, eps = mean ** (1 / 6), eps * s1**3 * s2**3 / mean
            else:
                sigma = sqrt(s1 * s2) if data["family"] == "OPLS-AA" else (s1 + s2) / 2
            scaled.addBond(index[a], index[b], [sigma / 10, eps * pair["lj"]])
    for force in (coul, lj, scaled):
        model.addForce(force)
    return model
