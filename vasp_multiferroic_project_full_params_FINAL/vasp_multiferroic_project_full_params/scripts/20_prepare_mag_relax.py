#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations
import os
from pathlib import Path
from typing import List

from pymatgen.core import Structure
from _common import ensure_dir, pick_ids_from_table, get_verbose, log, vlog
from _vasp_inputs import sorted_structure_for_vasp

RUNS_ROOT = Path("runs")
RELAX0_TABLE = Path("results/tables/relax0_summary.csv")

MAG_STATES = ["FM","AFM1","LOW"]
MAG_SEED = {"Ti":1.0,"V":3.0,"Cr":5.0,"Mn":5.0,"Fe":5.0,"Co":3.0,"Ni":2.0,"Cu":1.0,"Mo":2.0,"W":2.0,"Ru":2.0,"Ir":1.0,"Ce":1.0,"Pr":3.0,"Nd":3.0,"Sm":5.0,"Eu":7.0,"Gd":7.0,"Tb":6.0,"Dy":5.0}

ENCUT = int(os.environ.get("ENCUT","520"))
KSPACING = float(os.environ.get("KSPACING_MAG","0.25"))
NCORE = int(os.environ.get("NCORE","4"))

def make_magmom(struct: Structure, mode: str) -> List[float]:
    base = [MAG_SEED.get(site.specie.symbol, 0.0) for site in struct.sites]
    if mode == "FM":
        return base
    if mode == "LOW":
        return [0.5*b for b in base]
    if mode == "AFM1":
        out, sign = [], 1.0
        for b in base:
            out.append(sign*b); sign *= -1.0
        return out
    return base

def write_incar(path: Path, magmom: List[float], mode: str):
    lines = [
        f"SYSTEM = mag_relax_{mode}",
        "ISTART = 0",
        "ICHARG = 2",
        f"ENCUT = {ENCUT}",
        "PREC = Accurate",
        "ALGO = Normal",
        "ISPIN = 2",
        "ISYM = 0",
        "EDIFF = 1E-5",
        "NELM = 180",
        "ISMEAR = 1",
        "SIGMA = 0.2",
        f"KSPACING = {KSPACING}",
        "KGAMMA = .FALSE.",
        "IBRION = 2",
        "NSW = 200",
        "ISIF = 2",
        "EDIFFG = -0.02",
        "LREAL = .FALSE.",
        "LASPH = .TRUE.",
        "ADDGRID = .TRUE.",
        "LWAVE = .FALSE.",
        "LCHARG = .FALSE.",
        f"NCORE = {NCORE}",
        "MAGMOM = " + " ".join(f"{m:.2f}" for m in magmom),
    ]
    path.write_text("\n".join(lines)+"\n")

def write_pbs(path: Path, jobname: str):
    vasp_bin = os.environ.get("VASP_BIN","/public/apps/vasp/vasp_std")
    mpirun = os.environ.get("MPI_RUN","mpirun")
    ppn = int(os.environ.get("PBS_PPN","64"))
    wall = os.environ.get("PBS_WALLTIME","48:00:00")
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
    ids = pick_ids_from_table(RELAX0_TABLE, "keep", "Y")
    if not ids:
        ids = [p.name for p in RUNS_ROOT.glob("*") if (p/"relax0").is_dir()]

    for cid in ids:
        relax0 = RUNS_ROOT/cid/"relax0"
        contcar = relax0/"CONTCAR"
        potcar = relax0/"POTCAR"
        if not contcar.exists() or not potcar.exists():
            continue
        s = Structure.from_file(str(contcar))
        s = sorted_structure_for_vasp(s)

        for st in MAG_STATES:
            out = RUNS_ROOT/cid/f"mag_{st}"
            ensure_dir(out)
            (out/"POSCAR").write_text(contcar.read_text())
            (out/"POTCAR").write_bytes(potcar.read_bytes())
            kp = out/"KPOINTS"
            if kp.exists(): kp.unlink()
            magmom = make_magmom(s, st)
            write_incar(out/"INCAR", magmom, st)
            write_pbs(out/"job_mag.pbs", f"{cid[:26]}_{st}"[:36])
    print(f"[OK] prepared mag relax for {len(ids)} ids")

if __name__ == "__main__":
    main()
