"""Explicit experimental whole-molecule L2 projection; never an assignment engine."""

from collections.abc import Mapping
from dataclasses import dataclass
from math import fsum, isfinite, sqrt
from pathlib import Path

from island.exceptions import ChargeReferenceError
from island.workflows import storage
from island.workflows.bundle import system_from

from .correspondence import repeat_correspondence, require
from .observations import RawChargeObservation
from .records import pack, unpack

POLICY = "uniform_molecular_l2_v1"
CONSERVATION_TOLERANCE_E = 1e-12
PROJECTION_SCHEMA = "island_experimental_charge_projection_v1"


def project_molecular_charges(system, raw_charges, *, policy):
    """Return owned numerical diagnostics, not a parameterized molecule.

    Minimize sum((q'-q)**2) subject to sum(q')=Q_formal. All explicit
    atoms get the same offset, including H; no last-atom residual repair.
    """
    try:
        require(policy == POLICY, "Explicit supported projection policy required")
        corr = repeat_correspondence(system)
        require(
            corr["compatible"] and 0 < system.number_of_sites <= 100,
            "Unsupported projection scope",
        )
        require(isinstance(raw_charges, Mapping), "Charge mapping required")
        ids = sorted(system.topology.sites)
        require(
            all(type(s) is int for s in raw_charges) and set(raw_charges) == set(ids),
            "Exact integer site coverage required",
        )
        require(
            all(type(q) in (int, float) and isfinite(q) for q in raw_charges.values()),
            "Finite nonboolean charges required",
        )
        raw = {s: float(raw_charges[s]) for s in ids}
        total = fsum(raw.values())
        target = sum(system.topology.sites[s].formal_charge for s in ids)
        offset = (target - total) / len(ids)
        projected = {s: raw[s] + offset for s in ids}
        residual = fsum(projected.values()) - target
        require(
            isfinite(offset)
            and all(isfinite(v) for v in projected.values())
            and abs(residual) <= CONSERVATION_TOLERANCE_E,
            "Projection cannot satisfy fixed numerical conservation tolerance",
        )
        changes = [projected[s] - raw[s] for s in ids]
        repeats = []
        for r in corr["repeats"]:
            members = [
                s for g in r["groups"] for s in (g["site_id"], *g["hydrogen_ids"])
            ]
            a = fsum(raw[s] for s in members)
            b = fsum(projected[s] for s in members)
            repeats.append(
                {
                    "repeat_index": r["repeat_index"],
                    "role": r["role"],
                    "central": r["central"],
                    "atom_count": len(members),
                    "raw_total_e": a,
                    "projected_total_e": b,
                    "change_e": b - a,
                    "expected_change_e": len(members) * offset,
                }
            )
        return {
            "policy": POLICY,
            "numerical_tolerance_e": CONSERVATION_TOLERANCE_E,
            "raw_charges": raw,
            "projected_charges": projected,
            "raw_total_e": total,
            "raw_residual_e": total - target,
            "projected_total_e": fsum(projected.values()),
            "projected_residual_e": residual,
            "target_formal_charge_e": target,
            "atom_count": len(ids),
            "uniform_offset_e": offset,
            "maximum_correction_e": max(map(abs, changes)),
            "rms_correction_e": sqrt(fsum(x * x for x in changes) / len(changes)),
            "repeat_changes": repeats,
            "unit": "elementary_charge",
            "description": "post-BCC charges with explicit conservation projection; not recovered QM precision",
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
    except ChargeReferenceError:
        raise
    except Exception as error:
        raise ChargeReferenceError(f"Cannot project charges: {error}") from error


@dataclass(frozen=True)
class ChargeProjection:
    json_text: str

    def validate_integrity(self):
        try:
            p = unpack(self.json_text)
            require(
                set(p)
                == {"schema", "observation", "observation_identity", "projection"}
                and p["schema"] == PROJECTION_SCHEMA,
                "Invalid projection schema",
            )
            observation = RawChargeObservation(pack(p["observation"]))
            data = observation.payload["data"]
            require(
                observation.identity == p["observation_identity"],
                "Observation identity mismatch",
            )
            calculated = project_molecular_charges(
                system_from(data["system"]),
                data["raw_charges"],
                policy=p["projection"]["policy"],
            )
            require(
                pack(calculated) == pack(p["projection"]),
                "Projection/metrics differ from independently recomputed values",
            )
        except ChargeReferenceError:
            raise
        except Exception as error:
            raise ChargeReferenceError(f"Invalid projection: {error}") from error

    @property
    def payload(self):
        self.validate_integrity()
        return unpack(self.json_text)

    @property
    def identity(self):
        return storage.checksum(storage.json_bytes(storage.encode(self.payload)))


def project_charge_observation(observation, *, policy):
    try:
        require(
            isinstance(observation, RawChargeObservation), "Expected raw observation"
        )
        p = observation.payload
        data = p["data"]
        result = ChargeProjection(
            pack(
                {
                    "schema": PROJECTION_SCHEMA,
                    "observation": p,
                    "observation_identity": observation.identity,
                    "projection": project_molecular_charges(
                        system_from(data["system"]), data["raw_charges"], policy=policy
                    ),
                }
            )
        )
        result.validate_integrity()
        return result
    except ChargeReferenceError:
        raise
    except Exception as error:
        raise ChargeReferenceError(f"Cannot project observation: {error}") from error


def load_charge_projection(path):
    try:
        result = ChargeProjection(Path(path).read_text())
        result.validate_integrity()
        return result
    except Exception as error:
        raise ChargeReferenceError(f"Cannot load projection: {error}") from error
