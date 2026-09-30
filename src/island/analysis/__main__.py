"""Thin CLI: python -m island.analysis RUN --output REPORT_DIRECTORY."""

import argparse

from island.analysis import (
    AnalysisOptions,
    GeometryOptions,
    analyze_workflow,
    export_analysis,
)
from island.exceptions import AnalysisError


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workflow")
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--selection", choices=("all", "heavy", "explicit"), default="all"
    )
    parser.add_argument("--site-ids", nargs="+", type=int)
    parser.add_argument("--endpoints", nargs=2, type=int)
    parser.add_argument("--weighting", choices=("mass", "uniform"), default="mass")
    for name in ("start-step", "end-step", "stride"):
        parser.add_argument(
            "--" + name, type=int, default=1 if name == "stride" else None
        )
    for name in ("start-time-ps", "end-time-ps"):
        parser.add_argument("--" + name, type=float)
    a = parser.parse_args()
    try:
        options = AnalysisOptions(
            GeometryOptions(a.selection, a.site_ids, a.weighting, a.endpoints),
            a.start_step,
            a.end_step,
            a.start_time_ps,
            a.end_time_ps,
            a.stride,
        )
        report = analyze_workflow(a.workflow, options)
        print(export_analysis(report, a.output))
        print(
            f"{report.payload['sample_count']} retained samples; not an equilibrium estimate"
        )
    except AnalysisError as error:
        parser.exit(2, f"analysis: {error}\n")


if __name__ == "__main__":
    main()
