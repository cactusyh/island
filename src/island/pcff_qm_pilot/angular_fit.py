"""Preregistered bounded nested models; training-only selection and uncertainty."""

import numpy as np

from island.pcff_completion import require

from .angular_baseline import Configuration, source_extensions
from .stretch_baseline import source_baseline


def metrics(configurations, samples, values, offset):
    rows = []
    for c, s in zip(configurations, samples, strict=True):
        r = c.evaluate(values)
        df = r["forces"] - np.array(s["forces"])
        rows.append(
            {
                "key": s.get("key", s.get("id")),
                "energy_error_kcal_mol": float(r["energy"] + offset - s["energy"]),
                "force_errors_kcal_mol_A": df.tolist(),
                "force_RMS_kcal_mol_A": float(np.sqrt(np.mean(df * df))),
                "force_max_kcal_mol_A": float(np.max(abs(df))),
            }
        )
    return {
        "configurations": rows,
        "energy_RMSE": float(
            np.sqrt(np.mean([r["energy_error_kcal_mol"] ** 2 for r in rows]))
        ),
        "force_RMS": float(
            np.sqrt(np.mean([r["force_RMS_kcal_mol_A"] ** 2 for r in rows]))
        ),
        "force_max": max(r["force_max_kcal_mol_A"] for r in rows),
    }


class Problem:
    def __init__(self, base, q, target, full_nonstretch_response):
        self.model = base.model
        self.training = base.training
        self.center = base.center
        self.q = np.array(q)
        self.target = np.array(target)
        self.full_response = np.array(full_nonstretch_response)
        self.configurations = [
            Configuration(self.model, s["coordinates"]) for s in self.training
        ]
        self.stationary = Configuration(self.model, self.center["coordinates"])
        self.names = self.stationary.names
        self.upper = np.triu_indices(self.q.shape[1])
        self.step = 0.002
        self._last = None
        require(
            self.q.shape == (27, 12)
            and np.max(abs(self.q.T @ self.q - np.eye(12))) < 1e-10,
            "angular_combined_projection_mapping",
        )
        x = np.array(self.center["coordinates"])
        self.plus = [
            Configuration(self.model, x + (self.step * self.q[:, i]).reshape(9, 3))
            for i in range(12)
        ]
        self.minus = [
            Configuration(self.model, x - (self.step * self.q[:, i]).reshape(9, 3))
            for i in range(12)
        ]
        _, k2 = source_extensions(self.model)
        b = source_baseline(self.model)
        self.lower = np.r_[
            b[:2] - 0.08, [50, 50], np.full(20, -1000), -5, -0.08, -0.08, 1 - k2
        ]
        self.upper_bounds = np.r_[
            b[:2] + 0.08, [1000, 1000], np.full(20, 1000), 5, 0.08, 0.08, 300 - k2
        ]

    def curvature(self, p):
        a, b = [c.evaluate(p) for c in self.plus], [c.evaluate(p) for c in self.minus]
        response = -np.column_stack(
            [(v["forces"] - w["forces"]).ravel() for v, w in zip(a, b, strict=True)]
        ) / (2 * self.step)
        response_j = -np.stack(
            [
                (v["force_parameter_jacobian"] - w["force_parameter_jacobian"]).reshape(
                    27, 28
                )
                for v, w in zip(a, b, strict=True)
            ],
            axis=1,
        ) / (2 * self.step)
        h = self.q.T @ response
        j = np.einsum("ia,ibp->abp", self.q, response_j)
        return (h + h.T) / 2, (j + j.transpose(1, 0, 2)) / 2, response

    def system(self, p):
        p = np.asarray(p)
        if self._last is not None and np.array_equal(p, self._last[0]):
            return self._last[1:]
        calculated = [c.evaluate(p) for c in self.configurations]
        anchor = calculated[0]
        r, j = [], []
        for sample, c in zip(self.training, calculated, strict=True):
            r.extend(
                [
                    np.array(
                        [
                            (
                                c["energy"]
                                - anchor["energy"]
                                - sample["energy"]
                                + self.training[0]["energy"]
                            )
                            / 0.05
                        ]
                    ),
                    (c["forces"] - np.array(sample["forces"])).ravel() / 0.5,
                ]
            )
            j.extend(
                [
                    (
                        c["energy_parameter_jacobian"]
                        - anchor["energy_parameter_jacobian"]
                    )[None, :]
                    / 0.05,
                    c["force_parameter_jacobian"].reshape(27, 28) / 0.5,
                ]
            )
        stationary = self.stationary.evaluate(p)
        h, hj, _ = self.curvature(p)
        r.extend(
            [
                (stationary["forces"] - np.array(self.center["forces"])).ravel()
                / 0.05
                / np.sqrt(27)
                * 10,
                (h - self.target)[self.upper] / 5 / np.sqrt(len(self.upper[0])),
            ]
        )
        j.extend(
            [
                stationary["force_parameter_jacobian"].reshape(27, 28)
                / 0.05
                / np.sqrt(27)
                * 10,
                hj[self.upper] / 5 / np.sqrt(len(self.upper[0])),
            ]
        )
        r, j = np.concatenate(r), np.vstack(j)
        self._last = (p.copy(), r, j)
        return r, j

    def fit(self, initial, count, maximum_evaluations=160):
        from scipy.optimize import least_squares

        require(count in (24, 25, 27, 28), "unregistered_nested_parameter_set")
        initial = np.array(initial, float)
        require(
            initial.shape == (28,)
            and np.all(initial >= self.lower)
            and np.all(initial <= self.upper_bounds),
            "nested_initial_bounds",
        )

        def full(p):
            a = initial.copy()
            a[:count] = p
            return a

        result = least_squares(
            lambda p: self.system(full(p))[0],
            initial[:count],
            jac=lambda p: self.system(full(p))[1][:, :count],
            bounds=(self.lower[:count], self.upper_bounds[:count]),
            x_scale="jac",
            ftol=1e-10,
            xtol=1e-10,
            gtol=1e-9,
            max_nfev=maximum_evaluations,
        )
        p = full(result.x)
        r, j = self.system(p)
        a = j[:, :count]
        norms = np.linalg.norm(a, axis=0)
        scaled = a / np.maximum(norms, 1e-30)
        _, sv, vt = np.linalg.svd(scaled, full_matrices=False)
        rank = int(sum(sv > sv[0] * 1e-10))
        condition = float(sv[0] / max(sv[-1], 1e-300))
        cov = (
            vt.T * np.where(sv > sv[0] * 1e-10, 1 / np.maximum(sv, 1e-300) ** 2, 0)
        ) @ vt
        cov = (
            cov
            / np.outer(np.maximum(norms, 1e-30), np.maximum(norms, 1e-30))
            * (r @ r)
            / max(len(r) - count, 1)
        )
        sigma = np.sqrt(np.maximum(np.diag(cov), 0))
        correlation = np.divide(
            cov,
            np.outer(sigma, sigma),
            out=np.zeros_like(cov),
            where=np.outer(sigma, sigma) > 0,
        )
        offset = (
            self.training[0]["energy"] - self.configurations[0].evaluate(p)["energy"]
        )
        h, _, response = self.curvature(p)
        stationary = self.stationary.evaluate(p)
        delta = response[:, 8:] - self.full_response
        colerr = np.sqrt(np.mean(delta * delta, axis=0))
        refcol = np.sqrt(np.mean(self.full_response**2, axis=0))
        bounds = [self.names[i] for i, v in enumerate(result.active_mask) if v]
        training = metrics(self.configurations, self.training, p, offset)
        gates = {
            "optimizer_converged": bool(result.success),
            "full_rank_and_condition": rank == count and condition < 1000,
            "no_active_bounds": not bounds,
            "training_energy": training["energy_RMSE"] <= 0.1,
            "training_forces": training["force_RMS"] <= 0.5,
            "stationarity": float(np.sqrt(np.mean(stationary["forces"] ** 2))) <= 0.05,
            "combined_curvature": float(
                np.linalg.norm(h - self.target) / np.linalg.norm(self.target)
            )
            <= 0.05,
            "nonstretch_full_response": bool(np.all(colerr <= 0.5 + 0.05 * refcol)),
        }
        return {
            "parameters": dict(zip(self.names, p.tolist(), strict=True)),
            "active_parameter_names": self.names[:count],
            "parameter_count": count,
            "rank": rank,
            "column_scaled_condition": condition,
            "scaled_singular_values": sv.tolist(),
            "active_bounds": bounds,
            "optimizer": {
                "success": bool(result.success),
                "status": result.status,
                "message": result.message,
                "evaluations": result.nfev,
                "optimality": result.optimality,
            },
            "weighted_SS": float(r @ r),
            "training_reference_offset": float(offset),
            "training": training,
            "stationary_force_RMS": float(np.sqrt(np.mean(stationary["forces"] ** 2))),
            "projected_curvature": h.tolist(),
            "projected_curvature_relative_error": float(
                np.linalg.norm(h - self.target) / np.linalg.norm(self.target)
            ),
            "nonstretch_Cartesian_response": response[:, 8:].tolist(),
            "nonstretch_response_column_RMS_errors": colerr.tolist(),
            "nonstretch_reference_column_RMS": refcol.tolist(),
            "residual_covariance_sigma_proxy": sigma.tolist(),
            "residual_covariance_correlation_proxy": correlation.tolist(),
            "training_gates": gates,
            "training_qualifies": all(gates.values()),
            "uncertainty_scope": "Local residual covariance proxy only. Correlated Hessian entries, Cartesian forces and method/source/model systematics are not independent measurements. No scientific acceptance from optimizer convergence.",
        }
