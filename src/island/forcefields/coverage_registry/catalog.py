"""Descriptive motifs on the final graph; these patterns never assign atom types."""

from dataclasses import dataclass

from island.forcefields.assignment_core.schema import (
    freeze_sequences,
    identity,
    require,
    unique,
)
from island.forcefields.assignment_core.topology import molecule, normalize


@dataclass(frozen=True)
class ChemicalMotif:
    id: str
    smarts: str
    description: str


@dataclass(frozen=True)
class MotifMatch:
    motif: str
    sites: tuple[int, ...]
    origin: str

    @property
    def identity(self):
        return identity(self)


DEFAULT_MOTIFS = (
    ChemicalMotif("aliphatic_ch", "[C;!a]", "Aliphatic carbon and attached H"),
    ChemicalMotif("aromatic", "[a]", "Aromatic atoms, including heteroaromatics"),
    ChemicalMotif("ether", "[#6;!$(C=O)]-[O;X2]-[#6;!$(C=O)]", "C-O-C ether"),
    ChemicalMotif("ester", "[#6]-[C](=[O])-[O]-[#6]", "Carboxylate ester"),
    ChemicalMotif("carbonyl", "[C]=[O]", "Carbonyl"),
    ChemicalMotif("amide", "[C](=[O])-[N]", "Amide"),
    ChemicalMotif("nitrile", "[C]#[N]", "Nitrile"),
    ChemicalMotif("halogenated_carbon", "[#6]-[F,Cl,Br,I]", "Carbon-halogen bond"),
    ChemicalMotif("sulfur", "[#16]", "Sulfur; oxidation state is not a typing claim"),
    ChemicalMotif("phosphorus", "[#15]", "Phosphorus; no generic parameter claim"),
    ChemicalMotif("silicon", "[#14]", "Silicon; no generic parameter claim"),
    ChemicalMotif("siloxane", "[#14]-[#8]-[#14]", "Si-O-Si environment"),
    ChemicalMotif("end_groups", "", "Builder-declared ends and explicit caps"),
    ChemicalMotif("crosslink_junctions", "", "Explicit crosslink or branch anchors"),
    ChemicalMotif("unclassified", "", "Atoms outside the descriptive catalog"),
)


@dataclass(frozen=True)
class ChemicalMotifCatalog:
    motifs: tuple[ChemicalMotif, ...] = DEFAULT_MOTIFS
    version: str = "island.j83.motifs.v1"

    def __post_init__(self):
        from rdkit import Chem

        freeze_sequences(self)
        unique([m.id for m in self.motifs], "duplicate_motif_identity")
        require(
            {"end_groups", "crosslink_junctions", "unclassified"}
            <= {m.id for m in self.motifs},
            "structural_motifs_required",
        )
        for m in self.motifs:
            require(bool(m.id and m.description), "incomplete_motif")
            require(
                not m.smarts or Chem.MolFromSmarts(m.smarts) is not None,
                "invalid_motif_smarts:" + m.id,
            )

    @property
    def identity(self):
        return identity(self)

    def detect(self, graph):
        from rdkit import Chem

        from island.builders.polymer import build_polymer_from_sequence

        # Bind topology metadata to the actual final graph, including all edits.
        require(normalize(graph.spec) == graph, "graph_spec_or_normalizer_mutation")
        mol = molecule(graph)
        matches = set()

        def add(name, indices, origin):
            indices = set(indices)
            indices |= {
                n.GetIdx()
                for i in tuple(indices)
                for n in mol.GetAtomWithIdx(i).GetNeighbors()
                if n.GetAtomicNum() == 1
            }
            if indices:
                matches.add(
                    MotifMatch(name, tuple(sorted(i + 1 for i in indices)), origin)
                )

        chemical_sites = set()
        for motif in self.motifs:
            if motif.smarts:
                query = Chem.MolFromSmarts(motif.smarts)
                # No match cap: a large polymer must not silently lose motifs.
                for hit in mol.GetSubstructMatches(query, uniquify=True, maxMatches=0):
                    add(motif.id, hit, "final_graph_smarts")
        for match in matches:
            chemical_sites.update(match.sites)

        spec = graph.spec
        baseline = build_polymer_from_sequence(
            dict(spec.repeat_units),
            spec.sequence or (spec.repeat_units[0][0],) * spec.degree_of_polymerization,
            generate_3d=False,
            add_hydrogens=False,
            tacticity=spec.tacticity,
            stereo_seed=spec.stereo_seed,
            head_end_group=spec.end_groups[0],
            tail_end_group=spec.end_groups[1],
        )
        # The unchanged J82 normalizer preserves these pre-edit site IDs.
        for role in ("head", "tail"):
            ends = {baseline.metadata["polymer"][role + "_site_id"] - 1}
            ends |= {
                i - 1
                for i, a in baseline.topology.sites.items()
                if a.metadata.get("end_group_role") == role
            }
            add("end_groups", ends, "builder_" + role + "_and_cap")
        for link in spec.crosslinks:
            add(
                "crosslink_junctions",
                (link.site1 - 1, link.site2 - 1),
                "explicit_crosslink",
            )
        for branch in spec.branches:
            add(
                "crosslink_junctions",
                (branch.anchor_site - 1,),
                "explicit_branch_anchor",
            )
        # Structural labels do not hide chemically unclassified atoms.
        unknown = set(graph.atom_order) - chemical_sites
        add("unclassified", (i - 1 for i in unknown), "outside_chemical_catalog")
        return tuple(sorted(matches, key=lambda m: (m.motif, m.sites, m.origin)))
