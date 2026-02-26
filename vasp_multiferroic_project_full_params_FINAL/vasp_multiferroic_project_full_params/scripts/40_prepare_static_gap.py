#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations
import os
from pathlib import Path
from _common import ensure_dir, load_csv, get_verbose

RUNS_ROOT = Path("runs")
MAG_GS_TABLE = Path("results/tables/mag_gs_summary.csv")

ENCUT = int(os.environ.get("ENCUT","520"))
KSPACING = float(os.environ.get("KSPACING_GAP","0.22"))
NCORE = int(os.environ.get("NCORE","4"))

def write_incar(path: Path):
    lines = [
        "SYSTEM = static_gap",
        "ISTART = 0",
        "ICHARG = 2",
        f"ENCUT = {ENCUT}",
        "PREC = Accurate",
        "ALGO = Normal",
        "ISPIN = 2",
        "ISYM = 0",
        "EDIFF = 1E-6",
        "NELM = 200",
        "NSW = 0",
        "IBRION = -1",
        "ISMEAR = 0",
        "SIGMA = 0.05",
        "LORBIT = 11",
        "LWAVE = .FALSE.",
        "LCHARG = .TRUE.",
        "LREAL = .FALSE.",
        "LASPH = .TRUE.",
        "ADDGRID = .TRUE.",
        f"KSPACING = {KSPACING}",
        "KGAMMA = .FALSE.",
        f"NCORE = {NCORE}",
    ]
    path.write_text("\n".join(lines)+"\n")

def write_pbs(path: Path, jobname: str):
    vasp_bin = os.environ.get("VASP_BIN","/public/apps/vasp/vasp_std")
    mpirun = os.environ.get("MPI_RUN","mpirun")
    ppn = int(os.environ.get("PBS_PPN","64"))
    wall = os.environ.get("PBS_WALLTIME","24:00:00")
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
    rows = load_csv(MAG_GS_TABLE)
    if not rows:
        raise SystemExit(f"Missing {MAG_GS_TABLE}. Run 22_collect_mag_gs.py first.")
    for r in rows:
        if r.get("keep") != "Y" or not r.get("gs_state"):
            continue
        cid = r["id"]
        gs = r["gs_state"]
        src = RUNS_ROOT/cid/f"mag_{gs}"
        pos = src/"CONTCAR" if (src/"CONTCAR").exists() else src/"POSCAR"
        pot = src/"POTCAR"
        if not pos.exists() or not pot.exists():
            continue
        out = RUNS_ROOT/cid/"static_gap"
        ensure_dir(out)
        (out/"POSCAR").write_text(pos.read_text())
        (out/"POTCAR").write_bytes(pot.read_bytes())
        kp = out/"KPOINTS"
        if kp.exists():
            kp.unlink()
        write_incar(out/"INCAR")
        write_pbs(out/"job_static_gap.pbs", f"{cid[:28]}_gap"[:36])
    print("[OK] prepared static_gap")

if __name__ == "__main__":
    main()
