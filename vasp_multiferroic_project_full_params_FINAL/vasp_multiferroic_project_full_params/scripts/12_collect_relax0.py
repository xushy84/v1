#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""12_collect_relax0.py

Collect relax0 status + basic metrics.
Enhanced:
  - n_atoms
  - reduced_formula (from CONTCAR if available else POSCAR)

Output:
  results/tables/relax0_summary.csv
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, Any, List
from _common import outcar_is_complete, parse_last_toten, parse_magmom_total, read_structure_volume_from_contcar, write_csv, now, get_verbose

RUNS_ROOT = Path("runs")
OUT_TABLE = Path("results/tables/relax0_summary.csv")

def read_formula_natoms(step: Path):
    try:
        from pymatgen.core import Structure
        p = step/"CONTCAR" if (step/"CONTCAR").exists() else step/"POSCAR"
        if not p.exists():
            return None, None
        s = Structure.from_file(str(p))
        return s.composition.reduced_formula.replace(" ",""), len(s)
    except Exception:
        return None, None


VERBOSE = get_verbose()

def main():
    rows: List[Dict[str, Any]] = []
    for jobdir in sorted(RUNS_ROOT.glob("*")):
        step = jobdir/"relax0"
        if not step.is_dir():
            continue
        outcar = step/"OUTCAR"
        contcar = step/"CONTCAR"
        complete = outcar_is_complete(outcar)
        formula, nat = read_formula_natoms(step)
        rows.append({
            "id": jobdir.name,
            "step": "relax0",
            "complete": "Y" if complete else "N",
            "formula": formula,
            "n_atoms": nat,
            "E_toten_eV": parse_last_toten(outcar) if complete else None,
            "M_tot": parse_magmom_total(outcar) if complete else None,
            "Volume_A3": read_structure_volume_from_contcar(contcar),
            "keep": "Y" if complete else "N",
            "collected_at": now(),
        })
    write_csv(OUT_TABLE, rows, ["id","step","complete","formula","n_atoms","E_toten_eV","M_tot","Volume_A3","keep","collected_at"])
    print(f"[OK] wrote {OUT_TABLE} ({len(rows)} rows)")

if __name__ == "__main__":
    main()
