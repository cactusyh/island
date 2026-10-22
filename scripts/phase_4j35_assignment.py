"""GAFF2 parameters with explicit, independently sourced external charges."""

import shutil

from island.forcefields import adopt_forcefield
from island.forcefields.ambertools.provided_preparation import parameterize_provided
from island.periodic import (
    PeriodicForceFieldConfig,
    periodic_chemical_view,
    prepare_periodic_forcefield,
)


def prepare(fixture, root, runtime, runtime_identity, provided, reference_identity):
    native = parameterize_provided(fixture.packed.system, fixture.packed.graph, provided, root / 'native-components', runtime=runtime, runtime_identity=runtime_identity, expected_reference_identity=reference_identity)
    prepared = adopt_forcefield(periodic_chemical_view(fixture.packed.system), 'gaff2', native)
    imported = native.imported_result
    config = PeriodicForceFieldConfig(family='GAFF2', source={'family': 'GAFF2', 'sha256': imported.source_sha256, 'source': imported.source}, box_lengths=fixture.packed.system.box.lengths, nonbonded_cutoff=5.0, lj_mixing_rule='lorentz_berthelot_12_6', one_four_scaling=(imported.nonbonded_policy.lj_scale_14, imported.nonbonded_policy.coulomb_scale_14), electrostatics_method='ewald', pme_tolerance=1e-5, provenance='J35 explicit external OPLS-AA source charges plus GAFF2 parameters; no AM1-BCC claim', evidence=(fixture.packed.graph.identity, fixture.packed.plan.identity, runtime_identity, provided.identity, reference_identity))
    backend = prepare_periodic_forcefield(fixture.packed.system, fixture.packed.graph, prepared, config)
    shutil.copytree(root / 'native-components', root / 'source-assets/amber')
    return backend
