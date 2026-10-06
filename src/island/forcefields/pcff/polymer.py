"""Final-graph native preparation and explicit topology-update diagnostics.

No interop files or external converter participate in these operations. Repeat
and crosslink annotations describe physical roles; they never select atom types.
"""

from copy import deepcopy
from dataclasses import dataclass
from itertools import combinations, pairwise

from island.charge_references.records import pack, unpack
from island.core import Coordinates
from island.workflows.bundle import system_data, system_from

from .automatic import chemical_graph, type_pcff_atoms
from .charges import identity
from .source import FLAGS, boundary, require


def _roles(
    system, native, *, schema="island_pcff_final_graph_evidence_v2", prepared=None
):
    require(
        schema
        in (
            "island_pcff_final_graph_evidence_v1",
            "island_pcff_final_graph_evidence_v2",
        ),
        "Unsupported final-graph evidence schema",
    )
    extended = schema.endswith("v2")
    assignment = unpack(native.assignment.json_text)
    crosslinks = {
        tuple(sorted(p))
        for p in system.metadata.get("topology_edits", {}).get("crosslinks", [])
    }
    bonds = {tuple(sorted((b.site1, b.site2))) for b in system.topology.bonds.values()}
    ends = {
        system.metadata.get("polymer", {}).get(k)
        for k in ("head_site_id", "tail_site_id")
    }
    ends.discard(None)
    end_hydrogens = {
        j for i, j in bonds if i in ends and system.topology.sites[j].element == "H"
    } | {i for i, j in bonds if j in ends and system.topology.sites[i].element == "H"}
    rows = []
    for row in assignment["assignments"]:
        sites = row["sites"]
        meta = [deepcopy(system.topology.sites[i].metadata) for i in sites]
        edges = [
            tuple(sorted((a, b)))
            for a, b in (combinations(sites, 2) if extended else pairwise(sites))
        ]
        rows.append(
            {
                "assignment": deepcopy(row),
                "ordered_site_provenance": [
                    {"site": i, "metadata": m} for i, m in zip(sites, meta, strict=True)
                ],
                "roles": {
                    "crosslink": any(e in crosslinks for e in edges),
                    "end_group": any(m.get("end_group_role") for m in meta)
                    or (extended and bool(set(sites) & (ends | end_hydrogens))),
                    "inter_repeat": len(
                        {(m.get("chain_id"), m.get("repeat_unit_index")) for m in meta}
                    )
                    > 1,
                    **(
                        {"chain_end_sites": sorted(set(sites) & ends)}
                        if extended
                        else {}
                    ),
                    "authoritative_bond_edges": [list(e) for e in edges if e in bonds],
                },
            }
        )
    return {
        "schema": schema,
        "system": system_data(system),
        "chemical_graph_identity": identity(chemical_graph(system)),
        "model_identity": native.identity,
        "source": native.assignment.source.identity,
        "interactions": rows,
        "strict_source_rows": True,
        **(
            {
                "prepared_identity": prepared.identity,
                "operational_profile_identity": prepared.operational_profile.identity,
            }
            if extended
            else {}
        ),
        **FLAGS,
    }


@dataclass(frozen=True)
class PCFFFinalGraphPreparation:
    """Owned evidence plus the existing owned facade; native fingerprints stay native."""

    json_text: str
    prepared: object

    @boundary
    def validate_integrity(self, system):
        from island.forcefields import PreparedForceField

        require(
            type(self.prepared) is PreparedForceField, "Owned native facade required"
        )
        self.prepared.validate_integrity(system)
        require(
            self.prepared.operational_profile is not None,
            "Strict evidence requires an operational profile",
        )
        native = self.prepared.native_result
        require(
            unpack(native.assignment.json_text)["parameter_coverage_complete"],
            "Strict evidence has missing source assignments",
        )
        require(
            all(
                t["origin"] == "source_row" and t["source_rows"]
                for t in unpack(native.json_text)["terms"]
            ),
            "Strict evidence cannot claim policy zeros",
        )
        require(
            pack(
                _roles(
                    system,
                    self.prepared.native_result,
                    schema=unpack(self.json_text)["schema"],
                    prepared=self.prepared,
                )
            )
            == self.json_text,
            "Final-graph provenance/preparation contradiction",
        )

    @property
    def payload(self):
        data = unpack(self.json_text)
        self.validate_integrity(system_from(data["system"]))
        return data

    @property
    def identity(self):
        return identity(self.payload)


@boundary
def prepare_pcff_polymer(system, options, *, mode="strict"):
    """Prepare a final construction/edit result, or return a non-executable ledger.

    Strict mode requires an explicit operational profile and every active source
    row. Diagnostic mode intentionally has no prepared-result field.
    """
    from island.forcefields import ForceFieldRequest, PCFFOptions, prepare_forcefield

    require(type(options) is PCFFOptions, "Validated PCFFOptions required")
    options.__post_init__()
    require(mode in ("strict", "diagnostic"), "Unknown PCFF preparation mode")
    if mode == "diagnostic":
        from .model import special_pair_policy
        from .operational import inspect_pcff_operational_support
        from .source import load_pcff_source

        source = load_pcff_source(options.source_path)
        typing = type_pcff_atoms(system, source, profile=options.typing_profile)
        if not typing.complete:
            return {
                "diagnostic_only": True,
                "typing": typing.payload,
                "model_complete": False,
                **FLAGS,
            }
        if options.resolution_policy is None:
            from .automatic import assign_automatic_pcff_charges
            from .class2 import assign_pcff_parameters
            from .model import define_pcff_model

            charges = assign_automatic_pcff_charges(system, typing)
            assignment = (
                assign_pcff_parameters(system, typing, charges)
                if charges.complete
                else None
            )
            model = (
                define_pcff_model(
                    assignment,
                    special_pairs=special_pair_policy(
                        lj=options.lj, coulomb=options.coulomb
                    ),
                )
                if assignment
                else None
            )
            return {
                "diagnostic_only": True,
                "typing": typing.payload,
                "charges": charges.payload,
                "assignment": assignment.payload if assignment else None,
                "model_complete": bool(
                    model and model.payload["model_definition_complete"]
                ),
                "model_diagnostics": model.payload["diagnostics"] if model else [],
                **FLAGS,
            }
        return inspect_pcff_operational_support(
            system,
            typing,
            resolution_policy=options.resolution_policy,
            special_pairs=special_pair_policy(lj=options.lj, coulomb=options.coulomb),
        )
    require(
        options.source_profile is not None,
        "Strict final-graph preparation needs an operational profile",
    )
    prepared = prepare_forcefield(system, ForceFieldRequest("pcff", options))
    return _finish_preparation(system, prepared)


def _finish_preparation(system, prepared):
    native = prepared.native_result
    data = unpack(native.assignment.json_text)
    require(
        data["parameter_coverage_complete"], "Missing/ambiguous required source row"
    )
    require(
        all(
            t["origin"] == "source_row" and t["source_rows"]
            for t in unpack(native.json_text)["terms"]
        ),
        "Strict mode cannot accept missing or policy-zero parameters",
    )
    result = PCFFFinalGraphPreparation(
        pack(_roles(system, native, prepared=prepared)), prepared
    )
    result.validate_integrity(system)
    return result


@boundary
def edit_polymer_topology(system, *, remove_hydrogens=(), add_crosslinks=()):
    """Return a copy after caller-specified bond edits; this is not curing/chemistry inference.

    Surviving IDs and coordinates are retained. Deleted H IDs are never reused.
    Chemical validity is enforced by subsequent native perception, not guessed.
    """
    system.validate()
    removed = tuple(remove_hydrogens)
    require(
        len(set(removed)) == len(removed) and all(type(i) is int for i in removed),
        "Invalid removed IDs",
    )
    links = tuple(tuple(p) for p in add_crosslinks)
    require(
        all(
            len(p) == 2 and all(type(i) is int for i in p) and p[0] != p[1]
            for p in links
        ),
        "Invalid crosslink endpoints",
    )
    require(len({tuple(sorted(p)) for p in links}) == len(links), "Duplicate crosslink")
    result = system.copy()
    for i in removed:
        require(
            i in result.topology.sites and result.topology.sites[i].element == "H",
            "Only explicit specified H removal supported",
        )
        result.topology.remove_site(i)
    result.coordinates = Coordinates(
        {i: system.coordinates.get(i) for i in result.topology.sites}
    )
    for i, j in links:
        require(
            i in result.topology.sites and j in result.topology.sites,
            "Unknown crosslink endpoint",
        )
        result.topology.add_bond(i, j, order=1, aromatic=False)
    old = result.metadata.setdefault("topology_edits", {})
    old.setdefault("crosslinks", []).extend([sorted(p) for p in links])
    old.setdefault("removed_hydrogens", []).extend(removed)
    old["policy"] = "explicit_stable_id_graph_edit_v1"
    result.validate()
    return result


@boundary
def reassign_pcff_polymer(previous_system, previous, system, options):
    """Validate old binding, then retype changed neighborhoods on the final graph.

    Global ring/component context and full integrity verification remain required.
    The receipt distinguishes local chemical decisions from the full publication
    check. No old interaction inventories or equilibrium dependencies are reused.
    """
    from .expanded import recognize

    previous.validate_integrity(previous_system)
    old, new = chemical_graph(previous_system), chemical_graph(system)
    old_atoms = {a["id"]: a for a in old["sites"]}
    new_atoms = {a["id"]: a for a in new["sites"]}
    old_b = {tuple(b["sites"]): b for b in old["bonds"]}
    new_b = {tuple(b["sites"]): b for b in new["bonds"]}
    changed = {
        i
        for i in old_atoms.keys() | new_atoms.keys()
        if old_atoms.get(i) != new_atoms.get(i)
    }
    for pair in old_b.keys() | new_b.keys():
        if old_b.get(pair) != new_b.get(pair):
            changed.update(pair)
    local = set(changed)
    for b in old["bonds"] + new["bonds"]:
        if set(b["sites"]) & changed:
            local.update(b["sites"])
    local &= new_atoms.keys()
    # Only profile v1 is authorized here; do not accidentally call another rule set.
    require(
        options.typing_profile == "island_pcff_source_graph_v1",
        "Local update profile unsupported",
    )
    entries, _environments, issues = recognize(
        new, defer_components=True, site_ids=local
    )
    require(not issues, "Invalid updated graph")
    # Reuse only validated source/type lookups; site inventories and every
    # equilibrium dependency are rebuilt from the new graph. Full native
    # validation independently re-resolves the candidate before publication.
    from island.forcefields import adopt_forcefield

    from .automatic import assign_automatic_pcff_charges
    from .class2 import PCFFClass2Result, derive
    from .model import define_pcff_model, special_pair_policy
    from .source import load_pcff_source

    source = load_pcff_source(options.source_path)
    previous_assignment = previous.prepared.native_result.assignment
    previous_assignment.validate_integrity(previous_system)
    require(
        source.identity == previous_assignment.source.identity,
        "Local parameter reuse requires the identical source",
    )
    prior = unpack(previous_assignment.json_text)
    require(
        prior["resolution_policy"] == options.resolution_policy
        if "resolution_policy" in prior
        else options.resolution_policy is None,
        "Local reuse resolution policy mismatch",
    )
    cache = {}
    for row in prior["assignments"]:
        if row["status"] == "assigned":
            cache[
                (
                    row.get("requested_family", row["family"]),
                    tuple(row["supplied_types"]),
                )
            ] = {
                k: deepcopy(v)
                for k, v in row.items()
                if k
                not in ("id", "family", "requested_family", "sites", "dependencies")
            }
    typing = type_pcff_atoms(system, source, profile=options.typing_profile)
    charges = assign_automatic_pcff_charges(
        system, typing, resolution_policy=options.resolution_policy
    )
    require(
        typing.complete and charges.complete, "Complete updated typing/charges required"
    )
    assignment = PCFFClass2Result(
        pack(
            derive(
                charges.payload,
                source,
                options.resolution_policy,
                resolution_cache=cache,
            )
        ),
        source,
    )
    assignment.validate_integrity(system)
    model = define_pcff_model(
        assignment,
        special_pairs=special_pair_policy(lj=options.lj, coulomb=options.coulomb),
    )
    profile = options.source_profile.profile()
    result = _finish_preparation(
        system, adopt_forcefield(system, "pcff", model, pcff_profile=profile)
    )
    actual = unpack(result.prepared.native_result.assignment.json_text)[
        "charge_record"
    ]["automatic_typing"]
    require(
        all(actual["assignments"].get(i) == label for i, label in entries.items()),
        "Local/full typing disagreement",
    )
    return result, {
        "schema": "island_pcff_local_update_receipt_v1",
        "previous_model": previous.prepared.native_result.identity,
        "new_model": result.prepared.native_result.identity,
        "changed_sites": sorted(changed),
        "locally_retyped_sites": sorted(local),
        "local_entries": entries,
        "reused_source_type_lookups": len(cache),
        "equilibrium_dependencies_rebuilt": True,
        "publication_check": "full native graph/charge/inventory/source validation",
        "no_cached_interaction_inventory": True,
        **FLAGS,
    }
