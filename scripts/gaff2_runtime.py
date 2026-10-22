"""Capture or validate an explicitly selected, immutable GAFF2 runtime manifest."""

import argparse
import json
from pathlib import Path

from island.forcefields.ambertools.runtime import (
    capture_runtime,
    declared_environment,
    validate_runtime,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--capture', type=Path)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--identity')
    parser.add_argument('--activate', action='store_true', help='Explicitly select manifest AMBERHOME/bin for validation')
    args = parser.parse_args()
    if args.capture:
        if args.manifest.exists():
            raise ValueError('Refusing to overwrite a runtime manifest')
        record = capture_runtime(args.capture)
        args.manifest.parent.mkdir(parents=True, exist_ok=True)
        args.manifest.write_text(json.dumps(record, indent=2, sort_keys=True) + '\n')
    else:
        record = json.loads(args.manifest.read_text())
        validate_runtime(record, expected_identity=args.identity, environment=declared_environment(record) if args.activate else None)
    print(record['identity'])


if __name__ == '__main__':
    main()
