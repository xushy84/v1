#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""50_symmetry_analyze.py

Analyze symmetry (space group + inversion) for each candidate.
Preference for structure source:
  1) runs/<id>/relax_tight/CONTCAR (or POSCAR)
  2) runs/<id>/static_gap/POSCAR (or CONTCAR)
  3) runs/<id>/mag_<GS>/CONTCAR (or POSCAR)
  4) runs/<id>/relax0/CONTCAR (or POSCAR)

Output:
  results/tables/symmetry_summary.csv
    id, sg_symbol, sg_number, has_inversion, polar_candidate, source, collected_at
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, Any, List
import numpy as np

from _common import ensure_dir, load_csv, write_csv, now, get_verbose

RUNS_ROOT = Path("runs")
MAG_GS_TABLE = Path("results/tables/mag_gs_summary.csv")
OUT_TABLE = Path("results/tables/symmetry_summary.csv")

def analyze_structure(struct_path: Path) -> Dict[str, Any]:
    out: Dict[str, Any] = {"sg_symbol": None, "sg_number": None, "has_inversion": None}
    try:
        from pymatgen.core import Structure
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
        s = Structure.from_file(str(struct_path))
        sga = SpacegroupAnalyzer(s, symprec=1e-3, angle_tolerance=5)
        out["sg_symbol"] = sga.get_space_group_symbol()
        out["sg_number"] = int(sga.get_space_group_number())
        ops = sga.get_point_group_operations()
        out["has_inversion"] = any(np.allclose(op.rotation_matrix, -np.eye(3), atol=1e-6) for op in ops)
    except Exception:
        pass
    return out

def pick_best_structure(cid: str, gs_state: str | None) -> tuple[Path | None, str]:
    for fn in ["CONTCAR", "POSCAR"]:
        p = RUNS_ROOT/cid/"relax_tight"/fn
        if p.exists():
            return p, f"relax_tight/{fn}"
    for fn in ["POSCAR", "CONTCAR"]:
        p = RUNS_ROOT/cid/"static_gap"/fn
        if p.exists():
            return p, f"static_gap/{fn}"
    if gs_state:
        for fn in ["CONTCAR", "POSCAR"]:
            p = RUNS_ROOT/cid/f"mag_{gs_state}"/fn
            if p.exists():
                return p, f"mag_{gs_state}/{fn}"
    for fn in ["CONTCAR", "POSCAR"]:
        p = RUNS_ROOT/cid/"relax0"/fn
        if p.exists():
            return p, f"relax0/{fn}"
    return None, ""


VERBOSE = get_verbose()

def main():
    ensure_dir(OUT_TABLE.parent)
    mag_rows = load_csv(MAG_GS_TABLE)
    mag_map = {r["id"]: r for r in mag_rows} if mag_rows else {}

    rows: List[Dict[str, Any]] = []
    for jobdir in sorted(RUNS_ROOT.glob("*")):
        cid = jobdir.name
        gs_state = mag_map.get(cid, {}).get("gs_state")
        p, src = pick_best_structure(cid, gs_state)
        if p is None:
            continue
        info = analyze_structure(p)
        if info["has_inversion"] is True:
            polar = "N"
        elif info["has_inversion"] is False:
            polar = "Y"
        else:
            polar = "UNK"
        rows.append({
            "id": cid,
            "sg_symbol": info["sg_symbol"],
            "sg_number": info["sg_number"],
            "has_inversion": info["has_inversion"],
            "polar_candidate": polar,
            "source": src,
            "collected_at": now(),
        })

    write_csv(OUT_TABLE, rows, ["id","sg_symbol","sg_number","has_inversion","polar_candidate","source","collected_at"])
    print(f"[OK] wrote {OUT_TABLE} ({len(rows)} rows)")

if __name__ == "__main__":
    main()
