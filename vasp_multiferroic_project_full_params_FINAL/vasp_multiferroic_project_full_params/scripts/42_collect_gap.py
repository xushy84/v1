#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""42_collect_gap.py

Collect static_gap results and write gap_summary.csv.
Primary source: vasprun.xml (pymatgen) for robust band gap and metal detection.
Fallback: attempt OUTCAR regex.

Enhanced outputs:
  - is_metal (bool from band structure when available)
  - n_kpoints (from vasprun when available)

Output:
  results/tables/gap_summary.csv
    id,complete,E_toten_eV,Eg_eV,metal_flag,is_metal,vbm_eV,cbm_eV,n_kpoints,keep,collected_at
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Any, List, Optional
import re

from _common import outcar_is_complete, parse_last_toten, write_csv, now, get_verbose, log, vlog

RUNS_ROOT = Path("runs")
OUT_TABLE = Path("results/tables/gap_summary.csv")

def parse_gap_from_vasprun(vasprun_xml: Path) -> Dict[str, Any]:
    out: Dict[str, Any] = {"Eg_eV": None, "is_metal": None, "vbm_eV": None, "cbm_eV": None, "n_kpoints": None}
    try:
        from pymatgen.io.vasp.outputs import Vasprun
        vr = Vasprun(str(vasprun_xml), parse_potcar_file=False)
        out["n_kpoints"] = int(getattr(vr, "nkpoints", None) or (len(vr.actual_kpoints) if vr.actual_kpoints else 0) or 0)
        bg = vr.get_band_structure()
        gap = bg.get_band_gap()
        out["Eg_eV"] = float(gap.get("energy")) if gap.get("energy") is not None else None
        out["is_metal"] = bool(bg.is_metal())
        if not out["is_metal"]:
            out["vbm_eV"] = float(bg.get_vbm()["energy"])
            out["cbm_eV"] = float(bg.get_cbm()["energy"])
        return out
    except Exception:
        return out

def parse_gap_from_outcar(outcar: Path) -> Dict[str, Any]:
    out: Dict[str, Any] = {"Eg_eV": None, "is_metal": None, "vbm_eV": None, "cbm_eV": None, "n_kpoints": None}
    txt = outcar.read_text(errors="ignore") if outcar.exists() else ""
    m = list(re.finditer(r"band\s+gap\s*[:=]\s*([-\d\.Ee+]+)", txt, re.IGNORECASE))
    if m:
        try:
            out["Eg_eV"] = float(m[-1].group(1))
        except Exception:
            pass
    return out

def metal_flag_from_gap(is_metal: Optional[bool], Eg_eV: Optional[float]) -> str:
    if is_metal is True:
        return "Y"
    if is_metal is False:
        return "N"
    if Eg_eV is not None:
        return "Y" if Eg_eV < 1e-3 else "N"
    return "UNK"


VERBOSE = get_verbose()

def main():
    rows: List[Dict[str, Any]] = []
    for jobdir in sorted(RUNS_ROOT.glob("*")):
        step = jobdir/"static_gap"
        if not step.is_dir():
            continue
        outcar = step/"OUTCAR"
        vasprun = step/"vasprun.xml"
        complete = outcar_is_complete(outcar)
        if not complete:
            rows.append({
                "id": jobdir.name,
                "complete": "N",
                "E_toten_eV": None,
                "Eg_eV": None,
                "metal_flag": "UNK",
                "is_metal": None,
                "vbm_eV": None,
                "cbm_eV": None,
                "n_kpoints": None,
                "keep": "N",
                "collected_at": now(),
            })
            continue

        E = parse_last_toten(outcar)
        if vasprun.exists():
            g = parse_gap_from_vasprun(vasprun)
        else:
            g = parse_gap_from_outcar(outcar)

        Eg = g.get("Eg_eV")
        is_metal = g.get("is_metal")
        vbm = g.get("vbm_eV")
        cbm = g.get("cbm_eV")
        nk = g.get("n_kpoints")

        rows.append({
            "id": jobdir.name,
            "complete": "Y",
            "E_toten_eV": E,
            "Eg_eV": Eg,
            "metal_flag": metal_flag_from_gap(is_metal, Eg),
            "is_metal": is_metal,
            "vbm_eV": vbm,
            "cbm_eV": cbm,
            "n_kpoints": nk,
            "keep": "Y",
            "collected_at": now(),
        })

    write_csv(
        OUT_TABLE,
        rows,
        fieldnames=[
            "id","complete","E_toten_eV","Eg_eV","metal_flag","is_metal","vbm_eV","cbm_eV","n_kpoints","keep","collected_at"
        ],
    )
    print(f"[OK] wrote {OUT_TABLE} ({len(rows)} rows)")

if __name__ == "__main__":
    main()
