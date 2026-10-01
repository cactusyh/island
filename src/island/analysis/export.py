"""Exclusive export directory, durable files, completion marker published last."""

import csv
import io
import shutil
from pathlib import Path

from island.exceptions import AnalysisError

from .workflow import AnalysisReport


def export_analysis(report, directory):
    from island.workflows import storage

    owned = False
    try:
        if type(report) is not AnalysisReport:
            raise AnalysisError("Expected AnalysisReport")
        # payload applies the same contract as validate_integrity(), before I/O.
        p = report.payload
        output = Path(directory).resolve()
        source = Path(p["source"]["directory"]).resolve()
        if output.is_relative_to(source) or source.is_relative_to(output):
            raise AnalysisError(
                "Analysis output must be separate from the source workflow"
            )
        output.mkdir(parents=True, exist_ok=False)
        owned = True
        table = io.StringIO(newline="")
        rows = []
        for frame in p["frames"]:
            g = frame["geometry"]
            row = {k: v for k, v in frame.items() if k != "geometry"}
            row.update(
                weighting=g["weighting"],
                selected_ids=" ".join(map(str, g["selected_ids"])),
                endpoint_ids=" ".join(map(str, g["endpoint_ids"] or [])),
                endpoint_convention=g["endpoint_convention"],
                Rg_angstrom=g["radius_of_gyration"],
                Re_angstrom=g["end_to_end_distance"],
                Re_over_Rg=g["re_over_rg"],
                kappa_squared=g["kappa_squared"],
                diagnostics="; ".join(g["diagnostics"]),
            )
            for i, axis in enumerate("xyz"):
                row[f"center_{axis}_angstrom"] = g["center"][i]
                row[f"lambda_{i + 1}_angstrom2"] = g["eigenvalues"][i]
                for j in range(i, 3):
                    row[f"G_{axis}{'xyz'[j]}_angstrom2"] = g["gyration_tensor"][i][j]
            rows.append(row)
        writer = csv.DictWriter(table, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
        files = {
            "report.json": storage.json_bytes(p),
            "frames.csv": table.getvalue().encode(),
        }
        for name, raw in files.items():
            storage.publish(output / name, raw)
        storage.publish(
            output / "COMPLETE.json",
            storage.json_bytes(
                {
                    "schema": "island.analysis-export.v1",
                    "files": {n: storage.checksum(raw) for n, raw in files.items()},
                }
            ),
        )
        return output
    except Exception as error:
        if owned:
            try:
                shutil.rmtree(output)
            except Exception as cleanup_error:  # noqa: BLE001 -- preserve original failure
                error.add_note(f"Export cleanup also failed: {cleanup_error}")
        raise AnalysisError(f"Analysis export failed: {error}") from error
