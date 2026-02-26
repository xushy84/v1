#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""60_prepare_relax_tight.py

Prepare tight relax (ISIF=3) ONLY for shortlisted candidates.

Shortlist gates (default, can override via env):
  1) mag_gs.keep == 'Y' and gs_state exists
  2) gap_summary.complete == 'Y'
  3) gap_summary.metal_flag == 'N'
  4) Eg_eV >= EG_MIN  (env EG_MIN, default 0.5 eV)
  5) polar_flag_summary.polar_candidate == 'Y'  (after polar_kick if any)

Outputs:
  - runs/<id>/relax_tight/{POSCAR,POTCAR,INCAR,job_relax_tight.pbs}
  - results/tables/tight_shortlist.csv  (pass/fail + reasons + key metrics)

Notes:
  - This script never writes/copies KPOINTS (KSPACING used in INCAR).
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Dict, Any, List

from _common import ensure_dir, load_csv, write_csv, now, get_verbose, log, vlog

RUNS_ROOT = Path(os.environ.get("RUNS_ROOT", "runs"))

MAG_GS = Path("results/tables/mag_gs_summary.csv")
GAP = Path("results/tables/gap_summary.csv")
POLAR = Path("results/tables/polar_flag_summary.csv")
SYM = Path("results/tables/symmetry_summary.csv")

OUT_SHORTLIST = Path("results/tables/tight_shortlist.csv")

ENCUT = int(os.environ.get("ENCUT", "520"))
KSPACING = float(os.environ.get("KSPACING_TIGHT", "0.22"))
NCORE = int(os.environ.get("NCORE", "4"))
EG_MIN = float(os.environ.get("EG_MIN", "0.5"))

def _f(x):
    try:
        if x is None:
            return None
        s = str(x).strip()
        if s == "" or s.lower() == "nan":
            return None
        return float(s)
    except Exception:
        return None

def write_incar(path: Path):
    lines = [
        "SYSTEM = relax_tight",
        "ISTART = 0",
        "ICHARG = 2",
        f"ENCUT = {ENCUT}",
        "PREC = Accurate",
        "ALGO = Normal",
        "ISPIN = 2",
        "ISYM = 0",
        "EDIFF = 1E-6",
        "NELM = 200",
        "ISMEAR = 0",
        "SIGMA = 0.05",
        f"KSPACING = {KSPACING}",
        "KGAMMA = .FALSE.",
        "IBRION = 2",
        "NSW = 250",
        "ISIF = 3",
        "EDIFFG = -0.01",
        "LREAL = .FALSE.",
        "LASPH = .TRUE.",
        "ADDGRID = .TRUE.",
        "LWAVE = .FALSE.",
        "LCHARG = .FALSE.",
        f"NCORE = {NCORE}",
    ]
    path.write_text("\n".join(lines) + "\n")

def write_pbs(path: Path, jobname: str):
    vasp_bin = os.environ.get("VASP_BIN", "/public/apps/vasp/vasp_std")
    mpirun = os.environ.get("MPI_RUN", "mpirun")
    ppn = int(os.environ.get("PBS_PPN", "64"))
    wall = os.environ.get("PBS_WALLTIME_TIGHT", os.environ.get("PBS_WALLTIME", "72:00:00"))
    path.write_text(
        "#!/bin/bash\n"
        f"#PBS -N {jobname}\n"
        f"#PBS -l nodes=1:ppn={ppn}\n"
        f"#PBS -l walltime={wall}\n"
        "#PBS -j oe\n#PBS -V\n\n"
        "cd $PBS_O_WORKDIR\n"
        f"{mpirun} -np {ppn} {vasp_bin} > vasp.out\n"
    )


VERBOSE = get_verbose()

def main():
    mag_rows = load_csv(MAG_GS)
    gap_rows = load_csv(GAP)
    pol_rows = load_csv(POLAR)
    sym_rows = load_csv(SYM)

    mag = {r.get("id"): r for r in mag_rows}
    gap = {r.get("id"): r for r in gap_rows}
    pol = {r.get("id"): r for r in pol_rows}
    sym = {r.get("id"): r for r in sym_rows}

    rows_out: List[Dict[str, Any]] = []
    prepared = 0

    for jobdir in sorted(RUNS_ROOT.glob("*")):
        cid = jobdir.name

        mr = mag.get(cid, {})
        gr = gap.get(cid, {})
        pr = pol.get(cid, {})
        sr = sym.get(cid, {})

        reasons: List[str] = []

        vlog(VERBOSE, "INFO", f"Check {cid}")

        gs_state = mr.get("gs_state")
        if mr.get("keep") != "Y" or not gs_state:
            reasons.append("no_mag_gs")

        if gr.get("complete") != "Y":
            reasons.append("gap_incomplete")

        metal_flag = gr.get("metal_flag")
        if metal_flag != "N":
            reasons.append(f"metal_flag={metal_flag}")

        Eg = _f(gr.get("Eg_eV"))
        if Eg is None:
            reasons.append("Eg_missing")
        elif Eg < EG_MIN:
            reasons.append(f"Eg<{EG_MIN}")

        if pr.get("polar_candidate") != "Y":
            reasons.append(f"polar={pr.get('polar_candidate')}")

        passed = (len(reasons) == 0)
        if not passed:
            vlog(VERBOSE, "INFO", f"  FAIL {cid}: {reasons}")
        else:
            vlog(VERBOSE, "INFO", f"  PASS {cid}: preparing relax_tight")

        # Prepare if passed
        src_struct = ""
        if passed:
            base_dir = RUNS_ROOT/cid/f"mag_{gs_state}"
            base = base_dir/"CONTCAR" if (base_dir/"CONTCAR").exists() else base_dir/"POSCAR"
            if not base.exists():
                reasons.append("missing_gs_structure")
                passed = False
            elif not (base_dir/"POTCAR").exists():
                reasons.append("missing_POTCAR")
                passed = False
            else:
                out = RUNS_ROOT/cid/"relax_tight"
                ensure_dir(out)
                (out/"POSCAR").write_text(base.read_text())
                (out/"POTCAR").write_bytes((base_dir/"POTCAR").read_bytes())
                # ensure no KPOINTS
                kp = out/"KPOINTS"
                if kp.exists():
                    kp.unlink()
                write_incar(out/"INCAR")
                write_pbs(out/"job_relax_tight.pbs", f"{cid[:26]}_rt"[:36])
                prepared += 1
                src_struct = f"mag_{gs_state}/{base.name}"

        rows_out.append({
            "id": cid,
            "pass_tight": "Y" if passed else "N",
            "reasons": ";".join(reasons) if reasons else "",
            "gs_state": gs_state,
            "E_gs_eV": mr.get("E_gs_eV"),
            "dE_2nd_eV": mr.get("dE_2nd_eV"),
            "gap_complete": gr.get("complete"),
            "Eg_eV": gr.get("Eg_eV"),
            "metal_flag": metal_flag,
            "polar_candidate": pr.get("polar_candidate", sr.get("polar_candidate")),
            "sg_symbol": sr.get("sg_symbol"),
            "sg_number": sr.get("sg_number"),
            "has_inversion": sr.get("has_inversion"),
            "tight_source_structure": src_struct,
            "Eg_min_used": EG_MIN,
            "collected_at": now(),
        })

    write_csv(
        OUT_SHORTLIST,
        rows_out,
        fieldnames=[
            "id","pass_tight","reasons",
            "gs_state","E_gs_eV","dE_2nd_eV",
            "gap_complete","Eg_eV","metal_flag",
            "polar_candidate","sg_symbol","sg_number","has_inversion",
            "tight_source_structure","Eg_min_used","collected_at"
        ],
    )
    log("INFO", f"EG_MIN used = {EG_MIN} eV")
    log("INFO", f"Prepared relax_tight for {prepared} candidates")
    log("INFO", f"Wrote shortlist table: {OUT_SHORTLIST} ({len(rows_out)} rows)")

if __name__ == "__main__":
    main()
