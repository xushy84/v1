#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""82_collect_berryP.py

Parse Berry-phase polarization from OUTCAR and write polarization_summary.csv.

Output:
  results/tables/polarization_summary.csv
    id, complete, P_ux, P_uy, P_uz, P_uCcm2, note, collected_at

Parsing notes:
  VASP prints Berry-phase polarization in OUTCAR when LCALCPOL=.TRUE.
  Output formats vary. This script tries multiple patterns and falls back to UNK.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
import re
import math

from _common import outcar_is_complete, write_csv, now, get_verbose, log, vlog

RUNS_ROOT = Path("runs")
OUT_TABLE = Path("results/tables/polarization_summary.csv")

def parse_pol(outcar: Path) -> Tuple[Optional[float], Optional[float], Optional[float], str]:
    txt = outcar.read_text(errors="ignore") if outcar.exists() else ""

    # Pattern A: "Total polarization (eÅ)" or similar
    # Example variants include:
    #  "Total polarization (e Angstrom):    px  py  pz"
    pats = [
        r"Total\s+polarization[^\n]*\n\s*([-\d\.Ee+]+)\s+([-\d\.Ee+]+)\s+([-\d\.Ee+]+)",
        r"total\s+polarization[^\n]*\(e\s*\*\s*Angst\)[^\n]*\n\s*([-\d\.Ee+]+)\s+([-\d\.Ee+]+)\s+([-\d\.Ee+]+)",
    ]
    for pat in pats:
        m = list(re.finditer(pat, txt, re.IGNORECASE))
        if m:
            try:
                px, py, pz = map(float, m[-1].groups())
                return px, py, pz, "eA"
            except Exception:
                pass

    # Pattern B: "dipolmoment" (rare/older)
    m = list(re.finditer(r"dipolmoment\s*\(e\*Angst\)\s*:\s*([-\d\.Ee+]+)\s+([-\d\.Ee+]+)\s+([-\d\.Ee+]+)", txt, re.IGNORECASE))
    if m:
        try:
            px, py, pz = map(float, m[-1].groups())
            return px, py, pz, "eA"
        except Exception:
            pass

    return None, None, None, "UNK"


VERBOSE = get_verbose()

def main():
    rows: List[Dict[str, Any]] = []
    done=0; todo=0
    for jobdir in sorted(RUNS_ROOT.glob("*")):
        step = jobdir/"berryP"
        if not step.is_dir():
            continue
        outcar = step/"OUTCAR"
        complete = outcar_is_complete(outcar)
        if not complete:
            todo += 1
            rows.append({"id": jobdir.name, "complete":"N", "P_ux":None, "P_uy":None, "P_uz":None,
                         "P_uCcm2":None, "note":"incomplete", "collected_at": now()})
            continue

        px, py, pz, unit = parse_pol(outcar)

        # Convert eÅ to μC/cm^2 requires cell volume; we approximate by reading CONTCAR/POSCAR volume.
        P_uc = None
        note = f"unit={unit}"
        try:
            from pymatgen.core import Structure
            sp = step/"CONTCAR" if (step/"CONTCAR").exists() else step/"POSCAR"
            if sp.exists() and unit == "eA" and px is not None:
                s = Structure.from_file(str(sp))
                vol_A3 = s.volume
                # 1 eÅ / Å^3 = e / Å^2. Convert to C/m^2:
                # e = 1.602e-19 C, 1 Å^2 = 1e-20 m^2 => e/Å^2 = 1.602e1 C/m^2
                # Multiply by 100 to get μC/cm^2 (since 1 C/m^2 = 100 μC/cm^2)
                conv = 1.602176634e1 * 100.0  # μC/cm^2 per (eÅ/Å^3)
                pmod = math.sqrt(px*px + py*py + pz*pz)
                P_uc = conv * (pmod / vol_A3)
                note = f"unit=eA;volA3={vol_A3:.3f}"
        except Exception:
            pass

        done += 1
        rows.append({
            "id": jobdir.name,
            "complete":"Y",
            "P_ux": px,
            "P_uy": py,
            "P_uz": pz,
            "P_uCcm2": P_uc,
            "note": note,
            "collected_at": now(),
        })

    write_csv(OUT_TABLE, rows, ["id","complete","P_ux","P_uy","P_uz","P_uCcm2","note","collected_at"])
    log("INFO", f"Collected berryP: complete={done}, incomplete={todo}")
    log("INFO", f"Wrote table: {OUT_TABLE} ({len(rows)} rows)")

if __name__ == "__main__":
    main()
