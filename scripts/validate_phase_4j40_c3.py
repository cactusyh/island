"""Validate J40-C3 preflight and pressure controls from retained C2 failures."""
import argparse
import json
from pathlib import Path

from island.analysis.periodic_preflight import (
    PeriodicPreflightPolicy,
    preflight_box,
    pressure_unavailable,
)


def audit(root: Path) -> dict:
    rows=[]
    for p in sorted(root.glob("**/uninterrupted/failure.json")):
        rec=json.loads(p.read_text()); frames=rec.get("frames",[])
        boxes=[f["box_vectors_nm"] for f in frames if f.get("box_vectors_nm")]
        lengths=[min(float(v[i][i]) for i in range(3)) for v in boxes] if boxes else []
        row={"path":str(p),"failure":rec,"initial_min_box_nm":lengths[0] if lengths else None,"minimum_box_nm":min(lengths) if lengths else None}
        if lengths:
            row["preflight_at_failure"] = preflight_box(initial_box_nm=[lengths[0]]*3, predicted_min_box_nm=[lengths[-1]]*3, cutoff_nm=.5, policy=PeriodicPreflightPolicy())
        rows.append(row)
    return {"schema":"island.j40-c3-failure-audit.v1","failures":rows,"pressure":pressure_unavailable(reason="OpenMM Reference retained histories expose no validated full virial pressure; barostat target is not measurement",target_pressure_bar=1.0),"scientific_acceptance":False}

def main():
    p=argparse.ArgumentParser();p.add_argument("--run",type=Path,required=True);p.add_argument("--output",type=Path,required=True);a=p.parse_args();a.output.write_text(json.dumps(audit(a.run),indent=2)+"\n");return 1
if __name__ == "__main__": raise SystemExit(main())
