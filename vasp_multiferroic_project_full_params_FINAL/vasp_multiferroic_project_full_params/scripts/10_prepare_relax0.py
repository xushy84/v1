#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations

import os
from pathlib import Path
from typing import Dict, List, Tuple

from pymatgen.core import Structure
from pymatgen.io.vasp.inputs import Poscar

from _vasp_inputs import sorted_structure_for_vasp
from _common import ensure_dir, get_verbose

CIF_DIR = Path(os.environ.get("CIF_DIR", "data/cifs"))
RUNS_ROOT = Path(os.environ.get("RUNS_ROOT", "runs"))

POTCAR_ROOT = Path(os.environ.get("VASP_POTCAR_ROOT", "/public/home/xsy/VASP-POT/potpaw_pbe"))
VASP_BIN = os.environ.get("VASP_BIN", "/public/apps/vasp/vasp_std")
MPI_RUN = os.environ.get("MPI_RUN", "mpirun")
PBS_PPN = int(os.environ.get("PBS_PPN", "64"))
PBS_WALLTIME = os.environ.get("PBS_WALLTIME", "48:00:00")

ENCUT = int(os.environ.get("ENCUT", "520"))

KSPACING = float(os.environ.get("KSPACING_RELAX0", "0.30"))
NCORE = int(os.environ.get("NCORE", "4"))
EDIFF = os.environ.get("EDIFF_RELAX0", "1E-5")
EDIFFG = os.environ.get("EDIFFG_RELAX0", "-0.03")
NSW = int(os.environ.get("NSW_RELAX0", "180"))
ISIF = int(os.environ.get("ISIF_RELAX0", "2"))

MAG_SEED = {
    "Ti": 1.0, "V": 3.0, "Cr": 5.0, "Mn": 5.0, "Fe": 5.0,
    "Co": 3.0, "Ni": 2.0, "Cu": 1.0,
    "Mo": 2.0, "W": 2.0, "Ru": 2.0, "Ir": 1.0,
    "Ce": 1.0, "Pr": 3.0, "Nd": 3.0, "Sm": 5.0,
    "Eu": 7.0, "Gd": 7.0, "Tb": 6.0, "Dy": 5.0,
}
U_EFF = {
    "Cr": 3.5, "Mn": 4.0, "Fe": 4.0, "Co": 3.3,
    "Ni": 6.0, "V": 3.0, "Ti": 3.0,
    "Ce": 5.0, "Pr": 5.0, "Nd": 5.0, "Gd": 6.7,
}
D_ELEMENTS = {
    "Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
    "Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
    "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",
}
F_ELEMENTS = {
    "La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy",
    "Ho","Er","Tm","Yb","Lu",
    "Ac","Th","Pa","U","Np","Pu","Am",
}

POTCAR_PRIORITY_SUFFIX = ["_sv", "_pv", "_d", ""]
POTCAR_CACHE: Dict[Tuple[str, ...], bytes] = {}
AVAILABLE_POTCARS = {p.name for p in POTCAR_ROOT.iterdir() if p.is_dir()} if POTCAR_ROOT.exists() else set()

def find_best_potcar_symbol(el: str) -> str:
    priority = ["_sv","_pv",""] if (el in D_ELEMENTS or el in F_ELEMENTS) else POTCAR_PRIORITY_SUFFIX
    for suf in priority:
        name = el + suf
        if name in AVAILABLE_POTCARS:
            return name
    raise FileNotFoundError(f"No POTCAR for {el} under {POTCAR_ROOT}")

def concat_potcar(elements: List[str], out_path: Path):
    key = tuple(elements)
    if key in POTCAR_CACHE:
        out_path.write_bytes(POTCAR_CACHE[key])
        return
    blob = b""
    for el in elements:
        blob += (POTCAR_ROOT / find_best_potcar_symbol(el) / "POTCAR").read_bytes()
    POTCAR_CACHE[key] = blob
    out_path.write_bytes(blob)

def write_incar(path: Path, elements: List[str], magmom: List[float]):
    lines = [
        "SYSTEM = relax0_spin_U",
        "ISTART = 0",
        "ICHARG = 2",
        f"ENCUT = {ENCUT}",
        "PREC = Accurate",
        "ALGO = Normal",
        "ISPIN = 2",
        "ISYM = 0",
        f"EDIFF = {EDIFF}",
        "NELM = 160",
        "ISMEAR = 1",
        "SIGMA = 0.2",
        f"KSPACING = {KSPACING}",
        "KGAMMA = .FALSE.",
        "IBRION = 2",
        f"NSW = {NSW}",
        f"ISIF = {ISIF}",
        f"EDIFFG = {EDIFFG}",
        "LREAL = .FALSE.",
        "LASPH = .TRUE.",
        "ADDGRID = .TRUE.",
        "LWAVE = .FALSE.",
        "LCHARG = .FALSE.",
        f"NCORE = {NCORE}",
        "MAGMOM = " + " ".join(f"{m:.2f}" for m in magmom),
    ]
    any_u = any(el in U_EFF for el in elements)
    if any_u:
        lines += ["", "LDAU = .TRUE.", "LDAUTYPE = 2"]
        ldau_l, ldau_u = [], []
        for el in elements:
            u = U_EFF.get(el, 0.0)
            if el in F_ELEMENTS:
                ldau_l.append("3")
            elif el in D_ELEMENTS:
                ldau_l.append("2")
            else:
                ldau_l.append("-1")
            ldau_u.append(f"{u:.2f}")
        lines += [
            "LDAUL = " + " ".join(ldau_l),
            "LDAUU = " + " ".join(ldau_u),
            "LDAUJ = " + " ".join("0.00" for _ in elements),
            "LMAXMIX = 4",
        ]
    else:
        lines += ["", "LDAU = .FALSE."]
    path.write_text("\n".join(lines) + "\n")

def write_pbs(path: Path, jobname: str):
    path.write_text(
        "#!/bin/bash\n"
        f"#PBS -N {jobname}\n"
        f"#PBS -l nodes=1:ppn={PBS_PPN}\n"
        f"#PBS -l walltime={PBS_WALLTIME}\n"
        "#PBS -j oe\n#PBS -V\n\n"
        "cd $PBS_O_WORKDIR\n"
        f"{MPI_RUN} -np {PBS_PPN} {VASP_BIN} > vasp.out\n"
    )


VERBOSE = get_verbose()

def main():
    cifs = sorted(CIF_DIR.glob("*.cif"))
    if not cifs:
        raise SystemExit(f"No CIFs under {CIF_DIR}")
    ensure_dir(RUNS_ROOT)

    for idx, cif in enumerate(cifs):
        s = Structure.from_file(str(cif))
        s = sorted_structure_for_vasp(s)
        formula = s.composition.reduced_formula.replace(" ", "")
        job = f"{idx:04d}_{formula}"

        out = RUNS_ROOT / job / "relax0"
        ensure_dir(out)

        Poscar(s).write_file(out/"POSCAR")
        elem_order = list(Poscar(s).site_symbols)
        magmom = [MAG_SEED.get(site.specie.symbol, 0.0) for site in s.sites]
        write_incar(out/"INCAR", elem_order, magmom)

        kp = out/"KPOINTS"
        if kp.exists():
            kp.unlink()

        concat_potcar(elem_order, out/"POTCAR")
        write_pbs(out/"job_relax0.pbs", (job[:32] + "_r0"))
        print(f"[OK] {job} <- {cif.name}")

if __name__ == "__main__":
    main()
