"""Fresh exact-benchmark native RESP audit; unchanged J36-C2 policy."""

import json
from pathlib import Path

from island.forcefields.ambertools.psi4_execution import reconstruct
from island.forcefields.ambertools.resp import file_hash
from island.forcefields.ambertools.resp_policy import create
from island.forcefields.ambertools.resp_precision import audit_in_memory
from scripts.phase_4j36_c2_reference import check
from scripts.phase_4j38_c3_input import build_input
from scripts.phase_4j38_fixture import build


def write(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')


def charges(root, kind, runtime):
    root = Path(root)
    original, native = root / ('quantum-' + kind), root / ('native-esp-' + kind)
    graph = build(kind)[0].packed.graph
    data = build_input(kind)
    worker = file_hash('scripts/phase_4j36_quantum_reference.py')
    execution = reconstruct(original, graph, expected_input=data,
                            runtime_identity=runtime['identity'], worker_sha256=worker)
    write(original / 'execution-record.json', execution)
    audit = audit_in_memory(native, graph, native_directory=original,
                           in_memory_worker_sha256=file_hash('scripts/phase_4j38_c3_native_esp.py'),
                           expected_identity=execution['identity'], expected_input=data,
                           runtime_identity=runtime['identity'], worker_sha256=worker)
    write(root / (kind + '-source-audit.json'), audit)
    independent = check(native, audit)
    write(root / (kind + '-independent.json'), independent)
    result = create(audit, graph, runtime, independent)
    result.save(root / (kind + '-charges.json'))
    return result
