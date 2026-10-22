"""J34 explicit native component preparation for the disconnected linear fixture."""

import shutil

from island.forcefields import AmberToolsOptions, adopt_forcefield
from island.forcefields.ambertools.multimolecule import parameterize_multimolecule
from island.periodic import (
    PeriodicForceFieldConfig,
    periodic_chemical_view,
    prepare_periodic_forcefield,
)


def prepare_linear(fixture, root, runtime, runtime_identity):
    system = periodic_chemical_view(fixture.packed.system)
    native = parameterize_multimolecule(system, AmberToolsOptions('gaff2', 'am1bcc', amberhome=runtime['amberhome'], work_root=root, retain_success_artifacts=True, max_atoms=100, timeout_seconds=1800), root / 'native-components', runtime=runtime, runtime_identity=runtime_identity, box_lengths=fixture.packed.system.box.lengths)
    prepared = adopt_forcefield(system, 'gaff2', native)
    imported = native.imported_result
    config = PeriodicForceFieldConfig(family='GAFF2', source={'family': 'GAFF2', 'sha256': imported.source_sha256, 'source': imported.source}, box_lengths=fixture.packed.system.box.lengths, nonbonded_cutoff=5.0, lj_mixing_rule='lorentz_berthelot_12_6', one_four_scaling=(imported.nonbonded_policy.lj_scale_14, imported.nonbonded_policy.coulomb_scale_14), electrostatics_method='ewald', pme_tolerance=1e-5, provenance='J34 exact disjoint assembly of fresh complete GAFF2 molecule assignments', evidence=(fixture.packed.graph.identity, fixture.packed.plan.identity, runtime_identity))
    backend = prepare_periodic_forcefield(fixture.packed.system, fixture.packed.graph, prepared, config)
    shutil.copytree(root / 'native-components', root / 'source-assets/amber')
    return backend
