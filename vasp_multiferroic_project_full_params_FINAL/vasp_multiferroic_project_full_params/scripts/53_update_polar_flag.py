#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""53_update_polar_flag.py

Update polar_candidate based on polar_kick relaxation:
- If polar_kick completed: analyze its final structure for inversion.
- Else: fall back to symmetry_summary polar_candidate.

Output:
  results/tables/polar_flag_summary.csv
    id, polar_candidate, source, kick_sg, kick_sgnum, kick_has_inv, collected_at
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, Any, List
import numpy as np

from _common import load_csv, write_csv, now, outcar_is_complete, get_verbose, log, vlog

RUNS_ROOT = Path("runs")
SYM_TABLE = Path("results/tables/symmetry_summary.csv")
OUT_TABLE = Path("results/tables/polar_flag_summary.csv")

def analyze(path: Path) -> Dict[str, Any]:
    out = {"kick_sg": None, "kick_sgnum": None, "kick_has_inv": None}
    try:
        from pymatgen.core import Structure
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
        s = Structure.from_file(str(path))
        sga = SpacegroupAnalyzer(s, symprec=1e-3, angle_tolerance=5)
        out["kick_sg"] = sga.get_space_group_symbol()
        out["kick_sgnum"] = int(sga.get_space_group_number())
        ops = sga.get_point_group_operations()
        out["kick_has_inv"] = any(np.allclose(op.rotation_matrix, -np.eye(3), atol=1e-6) for op in ops)
    except Exception:
        pass
    return out


VERBOSE = get_verbose()

def main():
    sym_rows = load_csv(SYM_TABLE)
    if not sym_rows:
        raise SystemExit(f"Missing {SYM_TABLE}. Run 50_symmetry_analyze.py first.")

    rows: List[Dict[str, Any]] = []
    for r0 in sym_rows:
        cid = r0["id"]
        kick_dir = RUNS_ROOT/cid/"polar_kick"
        outcar = kick_dir/"OUTCAR"

        if kick_dir.is_dir() and outcar_is_complete(outcar):
            cand = kick_dir/"CONTCAR" if (kick_dir/"CONTCAR").exists() else kick_dir/"POSCAR"
            info = analyze(cand)
            if info["kick_has_inv"] is True:
                polar = "N"
            elif info["kick_has_inv"] is False:
                polar = "Y"
            else:
                polar = "UNK"
            rows.append({"id": cid, "polar_candidate": polar, "source": "polar_kick", **info, "collected_at": now()})
        else:
            rows.append({
                "id": cid,
                "polar_candidate": r0.get("polar_candidate","UNK"),
                "source": "symmetry_only",
                "kick_sg": None,
                "kick_sgnum": None,
                "kick_has_inv": None,
                "collected_at": now(),
            })

    write_csv(OUT_TABLE, rows, ["id","polar_candidate","source","kick_sg","kick_sgnum","kick_has_inv","collected_at"])
    print(f"[OK] wrote {OUT_TABLE} ({len(rows)} rows)")

if __name__ == "__main__":
    main()
