"""Analyze an existing saved run; creates reports, never runs simulation engines.

python examples/analyze_trajectory.py /path/to/run --output /path/to/new-report
Use --selection heavy or --weighting uniform to choose a different convention.
"""

from island.analysis.__main__ import main

if __name__ == "__main__":
    main()
