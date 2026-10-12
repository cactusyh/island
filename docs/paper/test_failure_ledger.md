# Test failure ledger

The full ordinary command `pytest -q` produced 1,462 passed, 180 skipped, 46
failed and 84 errors in the current environment. The suite is therefore not
green. The machine-readable ledger records the affected test groups, exception
classes, commands and evidence boundary in `test_failure_ledger.json`.

Most failures occur while importing the optional SciPy, OpenMM, ParmEd or
AmberTools stacks. The SciPy traceback shows an ABI mismatch between the Python
SciPy modules and `_rotation` extension. The GAFF/GAFF2 failures stop at the
prepared-bundle dependency boundary. These are environment/toolchain gates;
they are not treated as physical validation.

The PCFF near-miss and size-policy assertions remain listed as unresolved until
they can be reproduced in a clean supported environment. Tests were not
weakened and skips were not broadened.
