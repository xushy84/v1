#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""51_prepare_polar_kick.py

For centrosymmetric candidates (has_inversion=True), apply a small deterministic kick
and run a short relax (ISIF=2) to see if inversion can be broken.

Inputs:
  results/tables/symmetry_summary.csv
  results/tables/mag_gs_summary.csv (for GS structure preference)

Outputs:
  runs/<id>/polar_kick/{POSCAR,POTCAR,INCAR,job_polar_kick.pbs}
"""

from __future__ import annotations
import os
from pathlib import Path
from typing import Optional
import numpy as np

from pymatgen.core import Structure

from _common import ensure_dir, load_csv, get_verbose

RUNS_ROOT = Path("runs")
SYM_TABLE = Path("results/tables/symmetry_summary.csv")
MAG_GS_TABLE = Path("results/tables/mag_gs_summary.csv")

ENCUT = int(os.environ.get("ENCUT","520"))
KSPACING = float(os.environ.get("KSPACING_KICK","0.25"))
NCORE = int(os.environ.get("NCORE","4"))

def pick_base_dir(cid: str, gs_state: Optional[str]) -> Path | None:
    if (RUNS_ROOT/cid/"relax_tight").is_dir():
        return RUNS_ROOT/cid/"relax_tight"
    if gs_state and (RUNS_ROOT/cid/f"mag_{gs_state}").is_dir():
        return RUNS_ROOT/cid/f"mag_{gs_state}"
    if (RUNS_ROOT/cid/"relax0").is_dir():
        return RUNS_ROOT/cid/"relax0"
    return None

def read_best_structure(base_dir: Path) -> Path | None:
    for fn in ["CONTCAR","POSCAR"]:
        p = base_dir/fn
        if p.exists():
            return p
    return None

def kick_structure(s: Structure, amp: float = 0.01) -> Structure:
    coords = s.frac_coords.copy()
    vec = np.array([0.37, 0.11, 0.23])
    for i in range(len(coords)):
        coords[i] += amp * vec * ((i % 5) - 2)
    return Structure(s.lattice, s.species, coords)

def write_incar(path: Path):
    lines = [
        "SYSTEM = polar_kick_relax",
        "ISTART = 0",
        "ICHARG = 2",
        f"ENCUT = {ENCUT}",
        "PREC = Accurate",
        "ALGO = Normal",
        "ISPIN = 2",
        "ISYM = 0",
        "EDIFF = 1E-5",
        "NELM = 160",
        "ISMEAR = 1",
        "SIGMA = 0.2",
        f"KSPACING = {KSPACING}",
        "KGAMMA = .FALSE.",
        "IBRION = 2",
        "NSW = 80",
        "ISIF = 2",
        "EDIFFG = -0.03",
        "LREAL = .FALSE.",
        "LASPH = .TRUE.",
        "ADDGRID = .TRUE.",
        "LWAVE = .FALSE.",
        "LCHARG = .FALSE.",
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
    sym_rows = load_csv(SYM_TABLE)
    if not sym_rows:
        raise SystemExit(f"Missing {SYM_TABLE}. Run 50_symmetry_analyze.py first.")
    mag_rows = load_csv(MAG_GS_TABLE)
    mag_map = {r["id"]: r for r in mag_rows} if mag_rows else {}

    nprep = 0
    for r in sym_rows:
        cid = r["id"]
        has_inv = str(r.get("has_inversion","")).lower()
        if has_inv not in ("true","1","yes"):
            continue

        gs_state = mag_map.get(cid, {}).get("gs_state")
        base_dir = pick_base_dir(cid, gs_state)
        if base_dir is None:
            continue
        spath = read_best_structure(base_dir)
        if spath is None:
            continue

        s = Structure.from_file(str(spath))
        s2 = kick_structure(s, amp=0.01)

        out = RUNS_ROOT/cid/"polar_kick"
        ensure_dir(out)
        s2.to(fmt="poscar", filename=str(out/"POSCAR"))

        if (base_dir/"POTCAR").exists():
            (out/"POTCAR").write_bytes((base_dir/"POTCAR").read_bytes())

        kp = out/"KPOINTS"
        if kp.exists():
            kp.unlink()

        write_incar(out/"INCAR")
        write_pbs(out/"job_polar_kick.pbs", f"{cid[:28]}_pk"[:36])
        nprep += 1

    print(f"[OK] prepared polar_kick for {nprep} candidates")

if __name__ == "__main__":
    main()
