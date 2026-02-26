#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""22_collect_mag_gs.py

Collect energies of mag_FM / mag_AFM1 / mag_LOW and pick the lowest as gs_state.
Enhanced outputs:
  - E(FM/AFM1/LOW)
  - dE_2nd in eV and meV per cell (for quick degeneracy check)
  - completeness flags per state

Output:
  results/tables/mag_gs_summary.csv
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, Any, List
from _common import outcar_is_complete, parse_last_toten, parse_magmom_total, write_csv, now, get_verbose, log, vlog

RUNS_ROOT = Path("runs")
OUT_TABLE = Path("results/tables/mag_gs_summary.csv")
STATES = ["FM","AFM1","LOW"]


VERBOSE = get_verbose()

def main():
    rows: List[Dict[str, Any]] = []
    for jobdir in sorted(RUNS_ROOT.glob("*")):
        cid = jobdir.name

        energies: Dict[str, float | None] = {}
        mags: Dict[str, float | None] = {}
        completes: Dict[str, bool] = {}

        for st in STATES:
            step = jobdir/f"mag_{st}"
            outcar = step/"OUTCAR"
            completes[st] = outcar_is_complete(outcar)
            energies[st] = parse_last_toten(outcar) if completes[st] else None
            mags[st] = parse_magmom_total(outcar) if completes[st] else None

        # pick GS
        gs_state, gs_E = None, None
        for st in STATES:
            e = energies.get(st)
            if e is None:
                continue
            if gs_E is None or e < gs_E:
                gs_E, gs_state = e, st

        # delta to 2nd
        sorted_es = sorted([e for e in energies.values() if e is not None])
        dE = (sorted_es[1] - sorted_es[0]) if len(sorted_es) >= 2 else None
        dE_meV = (dE * 1000.0) if dE is not None else None

        rows.append({
            "id": cid,
            "gs_state": gs_state,
            "E_gs_eV": gs_E,
            "dE_2nd_eV": dE,
            "dE_2nd_meV_cell": dE_meV,
            "M_gs": mags.get(gs_state) if gs_state else None,
            "E_FM_eV": energies.get("FM"),
            "E_AFM1_eV": energies.get("AFM1"),
            "E_LOW_eV": energies.get("LOW"),
            "FM_complete": "Y" if completes["FM"] else "N",
            "AFM1_complete": "Y" if completes["AFM1"] else "N",
            "LOW_complete": "Y" if completes["LOW"] else "N",
            "keep": "Y" if gs_state else "N",
            "collected_at": now(),
        })

    write_csv(
        OUT_TABLE,
        rows,
        ["id","gs_state","E_gs_eV","dE_2nd_eV","dE_2nd_meV_cell","M_gs",
         "E_FM_eV","E_AFM1_eV","E_LOW_eV",
         "FM_complete","AFM1_complete","LOW_complete","keep","collected_at"]
    )
    print(f"[OK] wrote {OUT_TABLE} ({len(rows)} rows)")

if __name__ == "__main__":
    main()
