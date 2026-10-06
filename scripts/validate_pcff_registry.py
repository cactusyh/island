"""J9 polymer vertical slice, reusing the declared J8 independent harness.

The optional reference converter receives its own graph/type/charge inputs;
production coefficients are not exported as reference parameters.
"""

import sys

from validate_pcff_profile import main

if __name__ == "__main__":
    if "--declaration" not in sys.argv and "--child" not in sys.argv:
        sys.argv.extend(
            ["--declaration", "docs/evidence/phase_4j9_vertical_declaration.json"]
        )
    raise SystemExit(main())
