"""J61 constrained baseline/coupling fit, with held-out evaluation kept separate."""

from itertools import permutations, product

import numpy as np

from island.pcff_completion import digest, require

from .stretch_baseline import Configuration, source_baseline


def aligned_rms(left, right):
    a, b = np.array(left, float), np.array(right, float)
    a -= a.mean(axis=0)
    b -= b.mean(axis=0)
    values = []
    for first, second in product(permutations((3, 4, 5)), permutations((6, 7, 8))):
        for permutation in ([0, 1, 2, *first, *second], [2, 1, 0, *second, *first]):
            q = b[permutation]
            u, _, v = np.linalg.svd(q.T @ a)
            rotation = u @ np.diag([1, 1, np.linalg.det(u @ v)]) @ v
            values.append(
                float(np.sqrt(np.mean(np.sum((q @ rotation - a) ** 2, axis=1))))
            )
    return min(values)


def leakage(training, validation, historical=()):
    require(len(validation) >= 2, "fresh_validation_required")
    result = []
    pool = [*training, *historical]
    for index, sample in enumerate(validation):
        comparisons = [
            (
                s.get("key", s.get("id", "geometry")),
                aligned_rms(sample["coordinates"], s["coordinates"]),
            )
            for s in [*pool, *validation[:index]]
        ]
        name, rms = min(comparisons, key=lambda row: row[1])
        require(rms > 1e-7, "fresh_validation_symmetry_leakage")
        result.append(
            {
                "key": sample.get("key", str(index)),
                "nearest_seen_configuration": name,
                "symmetry_rigid_RMS_A": rms,
            }
        )
    return {"automorphisms": 72, "rows": result, "independently_sampled_basins": False}


def metrics(model, samples, parameters, offset):
    records = []
    for sample in samples:
        r = Configuration(model, sample["coordinates"]).evaluate(parameters)
        de = float(r["energy"] + offset - sample["energy"])
        df = r["forces"] - np.asarray(sample["forces"])
        records.append(
            {
                "key": sample.get("key", sample.get("id")),
                "geometry_identity": digest(sample["coordinates"]),
                "energy_error_kcal_mol": de,
                "force_errors_kcal_mol_A": df.tolist(),
                "force_RMS_kcal_mol_A": float(np.sqrt(np.mean(df**2))),
                "force_max_kcal_mol_A": float(np.max(abs(df))),
            }
        )
    return {
        "configurations": records,
        "energy_RMSE": float(
            np.sqrt(np.mean([s["energy_error_kcal_mol"] ** 2 for s in records]))
        ),
        "force_RMS": float(
            np.sqrt(
                np.mean(
                    [np.mean(np.square(s["force_errors_kcal_mol_A"])) for s in records]
                )
            )
        ),
        "force_max": max(s["force_max_kcal_mol_A"] for s in records),
    }


class Problem:
    """One stationary/curvature block, plus explicitly selected E/F training data."""

    def __init__(self, model, training, center, q, curvature, *, step=0.002):
        require(
            len(training) >= 4 and step > 0, "insufficient_training_or_invalid_step"
        )
        for sample in [*training, center]:
            require(
                sample["graph_identity"] == model["graph_identity"]
                and sample["atom_order"] == model["atom_order"],
                "baseline_reference_mapping_changed",
            )
            force = np.asarray(sample["forces"])
            require(
                force.shape == (9, 3)
                and np.isfinite(force).all()
                and np.isfinite(sample["energy"]),
                "invalid_baseline_reference",
            )
        self.model = model
        self.training = training
        self.configurations = [Configuration(model, s["coordinates"]) for s in training]
        self.stationary = Configuration(model, center["coordinates"])
        self.center = center
        self.q = np.asarray(q, float)
        require(
            self.q.shape == (27, 8)
            and np.max(abs(self.q.T @ self.q - np.eye(8))) < 1e-10,
            "invalid_stretch_basis",
        )
        self.curvature = np.asarray(curvature, float)
        require(
            self.curvature.shape == (8, 8) and np.isfinite(self.curvature).all(),
            "invalid_curvature_target",
        )
        self.step = step
        x = np.asarray(center["coordinates"])
        self.plus = [
            Configuration(model, x + (step * self.q[:, i]).reshape(9, 3))
            for i in range(8)
        ]
        self.minus = [
            Configuration(model, x - (step * self.q[:, i]).reshape(9, 3))
            for i in range(8)
        ]
        self.upper = np.triu_indices(8)
        self.names = self.stationary.names
        baseline = source_baseline(model)
        self.lower = np.r_[baseline[:2] - 0.08, [50, 50], np.full(20, -1000)]
        self.upper_bounds = np.r_[baseline[:2] + 0.08, [1000, 1000], np.full(20, 1000)]
        self._last = None

    def mm_curvature(self, parameters):
        a, b = (
            [c.evaluate(parameters) for c in self.plus],
            [c.evaluate(parameters) for c in self.minus],
        )
        response = np.column_stack(
            [
                -(x["forces"] - y["forces"]).ravel() / (2 * self.step)
                for x, y in zip(a, b, strict=True)
            ]
        )
        response_jac = np.stack(
            [
                -(
                    x["force_parameter_jacobian"] - y["force_parameter_jacobian"]
                ).reshape(27, 24)
                / (2 * self.step)
                for x, y in zip(a, b, strict=True)
            ],
            axis=1,
        )
        h = self.q.T @ response
        j = np.einsum("ia,ibp->abp", self.q, response_jac)
        return (h + h.T) / 2, (j + j.transpose(1, 0, 2)) / 2

    def system(self, parameters):
        parameters = np.asarray(parameters)
        if self._last is not None and np.array_equal(self._last[0], parameters):
            return self._last[1:]
        calculated = [c.evaluate(parameters) for c in self.configurations]
        ref_s, ref_p = self.training[0], calculated[0]
        residuals, jacobians = [], []
        for sample, p in zip(self.training, calculated, strict=True):
            residuals.extend(
                [
                    np.array(
                        [
                            (
                                p["energy"]
                                - ref_p["energy"]
                                - sample["energy"]
                                + ref_s["energy"]
                            )
                            / 0.05
                        ]
                    ),
                    (p["forces"] - np.asarray(sample["forces"])).ravel() / 0.5,
                ]
            )
            jacobians.extend(
                [
                    (
                        p["energy_parameter_jacobian"]
                        - ref_p["energy_parameter_jacobian"]
                    )[None, :]
                    / 0.05,
                    p["force_parameter_jacobian"].reshape(27, 24) / 0.5,
                ]
            )
        stationary = self.stationary.evaluate(parameters)
        h, hj = self.mm_curvature(parameters)
        residuals.extend(
            [
                (stationary["forces"] - np.asarray(self.center["forces"])).ravel()
                / 0.05
                / np.sqrt(27),
                (h - self.curvature)[self.upper] / 5 / 6,
            ]
        )
        jacobians.extend(
            [
                stationary["force_parameter_jacobian"].reshape(27, 24)
                / 0.05
                / np.sqrt(27),
                hj[self.upper] / 5 / 6,
            ]
        )
        r, j = np.concatenate(residuals), np.vstack(jacobians)
        self._last = (parameters.copy(), r, j)
        return r, j

    def fit(self, initial, *, fixed_baseline=False, maximum_evaluations=400):
        from scipy.optimize import least_squares

        initial = np.array(initial, float)
        require(
            initial.shape == (24,) and np.isfinite(initial).all(),
            "all_parameters_required",
        )
        require(
            np.all(initial >= self.lower) and np.all(initial <= self.upper_bounds),
            "initial_parameter_bounds",
        )
        active = np.arange(4, 24) if fixed_baseline else np.arange(24)

        def full(x):
            p = initial.copy()
            p[active] = x
            return p

        result = least_squares(
            lambda x: self.system(full(x))[0],
            initial[active],
            jac=lambda x: self.system(full(x))[1][:, active],
            bounds=(self.lower[active], self.upper_bounds[active]),
            x_scale="jac",
            ftol=1e-10,
            xtol=1e-10,
            gtol=1e-10,
            max_nfev=maximum_evaluations,
        )
        parameters = full(result.x)
        residual, jac = self.system(parameters)
        a = jac[:, active]
        norms = np.linalg.norm(a, axis=0)
        require(np.all(norms > 1e-12), "baseline_unidentifiable_zero_column")
        singular = np.linalg.svd(a / norms, compute_uv=False)
        rank = int(sum(singular > singular[0] * 1e-10))
        condition = (
            float(singular[0] / singular[-1]) if singular[-1] > 0 else float("inf")
        )
        require(
            rank == len(active) and condition < 1000,
            "baseline_joint_unidentifiable",
            rank=rank,
            condition=condition,
        )
        require(result.success, "baseline_optimizer_not_converged")
        offset = (
            self.training[0]["energy"]
            - self.configurations[0].evaluate(parameters)["energy"]
        )
        stationary = self.stationary.evaluate(parameters)
        h, _ = self.mm_curvature(parameters)
        covariance = (
            np.linalg.pinv(a.T @ a)
            * (residual @ residual)
            / max(len(residual) - len(active), 1)
        )
        sd = np.sqrt(np.maximum(np.diag(covariance), 0))
        denominators = np.outer(sd, sd)
        correlation = np.divide(
            covariance,
            denominators,
            out=np.zeros_like(covariance),
            where=denominators > 0,
        )
        return {
            "model": "J60_source_baseline_coupling_only"
            if fixed_baseline
            else self.stationary.evaluate(parameters)["model"],
            "origin": "fitted_external research-extension candidates",
            "parameters": dict(zip(self.names, parameters.tolist(), strict=True)),
            "active_parameter_names": [self.names[i] for i in active],
            "rank": rank,
            "scaled_condition": condition,
            "scaled_singular_values": singular.tolist(),
            "training_reference_offset": float(offset),
            "weighted_residual_SS": float(residual @ residual),
            "training": metrics(self.model, self.training, parameters, offset),
            "stationary_force_RMS": float(np.sqrt(np.mean(stationary["forces"] ** 2))),
            "stationary_force_errors": (
                stationary["forces"] - np.asarray(self.center["forces"])
            ).tolist(),
            "projected_curvature": h.tolist(),
            "projected_curvature_errors": (h - self.curvature).tolist(),
            "projected_curvature_relative_Frobenius": float(
                np.linalg.norm(h - self.curvature) / np.linalg.norm(self.curvature)
            ),
            "optimizer": {
                "implementation": "scipy.optimize.least_squares",
                "method": "trf",
                "evaluations": result.nfev,
                "status": result.status,
                "message": result.message,
                "optimality": result.optimality,
            },
            "active_bounds": [
                self.names[active[i]]
                for i, v in enumerate(result.active_mask)
                if v != 0
            ],
            "residual_only_sigma_proxy": sd.tolist(),
            "parameter_correlation": correlation.tolist(),
            "uncertainty_limitations": "Cartesian residual covariance is a proxy, not independent-measurement confidence. Hessian responses are one correlated block; full source/method/protocol systematics omitted.",
            "production_validated": False,
            "simulation_readiness": "not_established",
        }


def remaining_discrepancy(problem, parameters):
    """Training/stationary/curvature unit probes, with no new fitted coefficient.

    Fractions overlap. They test competing fixed-baseline/selected-form directions;
    they do not prove that a source row is wrong or supply a missing row.
    """
    from .hamiltonian import Jet, harmonic
    from .stretch_baseline import bond_kind

    def probes(configuration):
        out = {}

        def add(name, value):
            if name not in out:
                out[name] = Jet(0, np.zeros(27))
            out[name] += value

        for row in problem.model["interaction_inventory"]:
            f = row["family"]
            if f == "quartic_bond":
                k = bond_kind(problem.model, row["sites"])
                delta = configuration.bond_geometry[row["id"]] - parameters[k]
                for power in (3, 4):
                    add(
                        f"stretch_{'SiO' if k == 0 else 'SiH'}_power{power}",
                        delta**power,
                    )
            elif f == "quartic_angle":
                labels = tuple(problem.model["typing"][str(i)] for i in row["sites"])
                labels = min(labels, labels[::-1])
                delta = (
                    configuration.angle_geometry[row["id"]]
                    - row["assignment"]["values"][0]
                )
                for power in (1, 2):
                    add("angle_" + ",".join(labels) + f"_power{power}", delta**power)
            elif f == "torsion_3":
                add(
                    "fixed_main_torsion_n3_amplitude",
                    1 + harmonic(configuration.cosines[tuple(row["sites"])], 3),
                )
        for a in (4, 5, 6):
            for b in (7, 8, 9):
                add(
                    "outside_model_torsion_torsion_probe",
                    configuration.cosines[(a, 1, 2, 3)]
                    * configuration.cosines[(1, 2, 3, b)],
                )
        return out

    training = [probes(c) for c in problem.configurations]
    center = probes(problem.stationary)
    plus, minus = [probes(c) for c in problem.plus], [probes(c) for c in problem.minus]
    residual, jacobian = problem.system(parameters)
    total = float(residual @ residual)
    orthogonal = (
        residual - jacobian @ np.linalg.lstsq(jacobian, residual, rcond=1e-10)[0]
    )
    rows = []
    for name in sorted(center):
        feature = []
        for row in training:
            feature.extend(
                [
                    np.array([(row[name].v - training[0][name].v) / 0.05]),
                    -row[name].d / 0.5,
                ]
            )
        response = np.column_stack(
            [
                (a[name].d - b[name].d) / (2 * problem.step)
                for a, b in zip(plus, minus, strict=True)
            ]
        )
        h = problem.q.T @ response
        h = (h + h.T) / 2
        feature.extend([-center[name].d / 0.05 / np.sqrt(27), h[problem.upper] / 30])
        z = np.concatenate(feature)
        outside = z - jacobian @ np.linalg.lstsq(jacobian, z, rcond=1e-10)[0]
        cosine = (
            float(
                outside
                @ orthogonal
                / np.linalg.norm(outside)
                / np.linalg.norm(orthogonal)
            )
            if np.linalg.norm(outside) > 1e-10 and np.linalg.norm(orthogonal) > 1e-10
            else 0.0
        )
        rows.append(
            {
                "probe": name,
                "fraction_weighted_residual_SS_projected": cosine**2,
                "outside_current_parameter_Jacobian_norm": float(
                    np.linalg.norm(outside)
                ),
            }
        )
    return {
        "weighted_training_constraint_residual_SS": total,
        "fraction_outside_full24_Jacobian": float(orthogonal @ orthogonal / total)
        if total
        else 0.0,
        "probes": sorted(
            rows,
            key=lambda r: r["fraction_weighted_residual_SS_projected"],
            reverse=True,
        ),
        "interpretation": "Overlapping local unit-energy polynomial/periodic probes after removing the full24 Jacobian; fractions must not be added. Angular linear probes diagnose equilibrium-force bias, not a separately assigned linear PCFF source row. TT is a counterfactual outside this selected Hamiltonian, not an automatic applicable or zero term.",
        "held_out_used": False,
        "additional_coefficients_fitted_or_created": [],
    }
