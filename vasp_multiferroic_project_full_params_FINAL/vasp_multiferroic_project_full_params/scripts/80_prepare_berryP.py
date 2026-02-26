#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""80_prepare_berryP.py

Prepare Berry-phase polarization calculations for shortlisted candidates.

Inputs:
  results/tables/tight_shortlist.csv (pass_tight=Y) preferred
  fallback: results/tables/tight_summary.csv (complete=Y)

For each candidate:
  runs/<id>/berryP/
    POSCAR (from relax_tight/CONTCAR)
    POTCAR (from relax_tight/POTCAR)
    INCAR (LCALCPOL)
    job_berryP.pbs

Env:
  ENCUT (default 520)
  KSPACING_BERRY (default 0.22)
  NCORE (default 4)
  PBS_PPN (default 64)
  PBS_WALLTIME_BERRY (default 24:00:00)
  VASP_BIN, MPI_RUN

Notes:
  - Requires insulating systems. You should already have metal filtered before tight_shortlist.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Dict, Any, List

from _common import ensure_dir, load_csv, write_csv, now, get_verbose, log, vlog

RUNS_ROOT = Path(os.environ.get("RUNS_ROOT", "runs"))

SHORTLIST = Path("results/tables/tight_shortlist.csv")
TIGHT_SUMMARY = Path("results/tables/tight_summary.csv")
OUT_JOBS = Path("results/tables/berryP_jobs.csv")

ENCUT = int(os.environ.get("ENCUT", "520"))
KSPACING = float(os.environ.get("KSPACING_BERRY", "0.22"))
NCORE = int(os.environ.get("NCORE", "4"))

PBS_PPN = int(os.environ.get("PBS_PPN", "64"))
PBS_WALLTIME = os.environ.get("PBS_WALLTIME_BERRY", "24:00:00")

def pick_candidates() -> List[str]:
    if SHORTLIST.exists():
        rows = load_csv(SHORTLIST)
        return [r["id"] for r in rows if r.get("pass_tight") == "Y"]
    if TIGHT_SUMMARY.exists():
        rows = load_csv(TIGHT_SUMMARY)
        return [r["id"] for r in rows if r.get("complete") == "Y"]
    return []

def write_incar(path: Path):
    lines = [
        "SYSTEM = berryP_polarization",
        "ISTART = 0",
        "ICHARG = 2",
        f"ENCUT = {ENCUT}",
        "PREC = Accurate",
        "ALGO = Normal",
        "ISPIN = 2",
        "ISYM = 0",
        "EDIFF = 1E-7",
        "NELM = 250",
        "NSW = 0",
        "IBRION = -1",
        "ISMEAR = 0",
        "SIGMA = 0.05",
        "LREAL = .FALSE.",
        "LASPH = .TRUE.",
        "ADDGRID = .TRUE.",
        "LWAVE = .FALSE.",
        "LCHARG = .FALSE.",
        # Berry phase polarization
        "LCALCPOL = .TRUE.",
        # kmesh
        f"KSPACING = {KSPACING}",
        "KGAMMA = .FALSE.",
        f"NCORE = {NCORE}",
    ]
    path.write_text("\n".join(lines) + "\n")

def write_pbs(path: Path, jobname: str):
    vasp_bin = os.environ.get("VASP_BIN", "/public/apps/vasp/vasp_std")
    mpirun = os.environ.get("MPI_RUN", "mpirun")
    path.write_text(
        "#!/bin/bash\n"
        f"#PBS -N {jobname}\n"
        f"#PBS -l nodes=1:ppn={PBS_PPN}\n"
        f"#PBS -l walltime={PBS_WALLTIME}\n"
        "#PBS -j oe\n#PBS -V\n\n"
        "cd $PBS_O_WORKDIR\n"
        f"{mpirun} -np {PBS_PPN} {vasp_bin} > vasp.out\n"
    )


VERBOSE = get_verbose()

def main():
    cids = pick_candidates()
    log("INFO", f"Selected {len(cids)} candidates for berryP")
    log("INFO", f"KSPACING_BERRY={KSPACING}, NCORE={NCORE}")
    ensure_dir(OUT_JOBS.parent)
    rows_out: List[Dict[str, Any]] = []
    prepared = 0

    for cid in cids:
        base = RUNS_ROOT/cid/"relax_tight"
        pos = base/"CONTCAR" if (base/"CONTCAR").exists() else base/"POSCAR"
        pot = base/"POTCAR"
        if not pos.exists() or not pot.exists():
            rows_out.append({"id": cid, "status":"SKIP_MISSING_INPUTS", "note":"need relax_tight POSCAR/POTCAR", "collected_at": now()})
            continue

        out = RUNS_ROOT/cid/"berryP"
        ensure_dir(out)
        (out/"POSCAR").write_text(pos.read_text())
        (out/"POTCAR").write_bytes(pot.read_bytes())
        kp = out/"KPOINTS"
        if kp.exists():
            kp.unlink()
        write_incar(out/"INCAR")
        write_pbs(out/"job_berryP.pbs", f"{cid[:26]}_P"[:36])
        rows_out.append({"id": cid, "status":"PREPARED", "note":"LCALCPOL", "collected_at": now()})
        prepared += 1
        vlog(VERBOSE, "INFO", f"Prepared berryP for {cid}")

    write_csv(OUT_JOBS, rows_out, ["id","status","note","collected_at"])
    log("INFO", f"Prepared berryP for {prepared} candidates")
    log("INFO", f"Wrote jobs table: {OUT_JOBS} ({len(rows_out)} rows)")

if __name__ == "__main__":
    main()
