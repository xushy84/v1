#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""70_prepare_phonon.py

Phonon force jobs preparation via phonopy finite-displacement.

Inputs:
  - results/tables/tight_shortlist.csv (preferred): use pass_tight=Y
  - fallback: results/tables/tight_summary.csv (complete=Y) if shortlist missing

For each selected candidate, generate a supercell and displaced structures using phonopy,
then create VASP single-point force calculations for each displaced supercell.

Directory layout:
  runs/<id>/phonon_force/
    phonopy.yaml / phonopy_disp.yaml (if generated)
    disp-001/ POSCAR INCAR POTCAR job_phonon_force.pbs
    disp-002/ ...

Outputs:
  results/tables/phonon_jobs.csv
    id, n_displacements, supercell, source_structure, status, note, collected_at

Requirements:
  - phonopy installed and available in Python (import phonopy)
    If phonopy is not available, the script will still write phonon_jobs.csv with status=NO_PHONOPY.

Tunable env vars:
  RUNS_ROOT (default runs)
  PH_SUPERCELL (default "2 2 2")  # phonopy supercell matrix
  PH_DISP (default 0.01)          # displacement in Angstrom
  ENCUT (default 520)
  KSPACING_PH (default 0.25)      # larger supercell -> can be coarser
  NCORE (default 4)
  PBS_PPN (default 64)
  PBS_WALLTIME_PH (default 48:00:00)
  VASP_BIN, MPI_RUN

Notes:
  - No KPOINTS file is created; KSPACING + KGAMMA=.FALSE. in INCAR.
  - Uses POTCAR copied from mag_<GS> or relax_tight when available.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Dict, Any, List, Tuple

from _common import ensure_dir, load_csv, write_csv, now, get_verbose, log

RUNS_ROOT = Path(os.environ.get("RUNS_ROOT", "runs"))

SHORTLIST = Path("results/tables/tight_shortlist.csv")
TIGHT_SUMMARY = Path("results/tables/tight_summary.csv")
MAG_GS = Path("results/tables/mag_gs_summary.csv")

OUT_TABLE = Path("results/tables/phonon_jobs.csv")

ENCUT = int(os.environ.get("ENCUT", "520"))
KSPACING = float(os.environ.get("KSPACING_PH", "0.25"))
NCORE = int(os.environ.get("NCORE", "4"))

PH_SUPERCELL = os.environ.get("PH_SUPERCELL", "2 2 2")  # string
PH_DISP = float(os.environ.get("PH_DISP", "0.01"))

PBS_PPN = int(os.environ.get("PBS_PPN", "64"))
PBS_WALLTIME = os.environ.get("PBS_WALLTIME_PH", "48:00:00")

def _parse_supercell(s: str) -> List[List[int]]:
    parts = [int(x) for x in s.split()]
    if len(parts) == 3:
        return [[parts[0],0,0],[0,parts[1],0],[0,0,parts[2]]]
    if len(parts) == 9:
        return [parts[0:3], parts[3:6], parts[6:9]]
    raise ValueError("PH_SUPERCELL must be 'a b c' or 9 integers")

def pick_candidates() -> List[str]:
    if SHORTLIST.exists():
        rows = load_csv(SHORTLIST)
        return [r["id"] for r in rows if r.get("pass_tight") == "Y"]
    if TIGHT_SUMMARY.exists():
        rows = load_csv(TIGHT_SUMMARY)
        return [r["id"] for r in rows if r.get("complete") == "Y"]
    return []

def pick_source_structure(cid: str) -> Tuple[Path | None, str]:
    # Prefer relax_tight CONTCAR; else mag_GS CONTCAR; else relax0 CONTCAR
    for fn in ["CONTCAR","POSCAR"]:
        p = RUNS_ROOT/cid/"relax_tight"/fn
        if p.exists():
            return p, f"relax_tight/{fn}"
    # mag GS
    gs_state = None
    if MAG_GS.exists():
        mag_rows = load_csv(MAG_GS)
        for r in mag_rows:
            if r.get("id") == cid:
                gs_state = r.get("gs_state")
                break
    if gs_state:
        for fn in ["CONTCAR","POSCAR"]:
            p = RUNS_ROOT/cid/f"mag_{gs_state}"/fn
            if p.exists():
                return p, f"mag_{gs_state}/{fn}"
    for fn in ["CONTCAR","POSCAR"]:
        p = RUNS_ROOT/cid/"relax0"/fn
        if p.exists():
            return p, f"relax0/{fn}"
    return None, ""

def pick_potcar(cid: str) -> Path | None:
    # Prefer relax_tight POTCAR; else mag_GS POTCAR; else relax0 POTCAR
    p = RUNS_ROOT/cid/"relax_tight"/"POTCAR"
    if p.exists():
        return p
    gs_state = None
    if MAG_GS.exists():
        for r in load_csv(MAG_GS):
            if r.get("id") == cid:
                gs_state = r.get("gs_state")
                break
    if gs_state:
        p = RUNS_ROOT/cid/f"mag_{gs_state}"/"POTCAR"
        if p.exists():
            return p
    p = RUNS_ROOT/cid/"relax0"/"POTCAR"
    return p if p.exists() else None

def write_incar(path: Path):
    lines = [
        "SYSTEM = phonon_force",
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
        "LWAVE = .FALSE.",
        "LCHARG = .FALSE.",
        "LREAL = .FALSE.",
        "LASPH = .TRUE.",
        "ADDGRID = .TRUE.",
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
    log("INFO", f"Selected {len(cids)} candidates for phonon_force")
    log("INFO", f"Supercell={PH_SUPERCELL}, disp={PH_DISP} Å, KSPACING_PH={KSPACING}, NCORE={NCORE}")
    ensure_dir(OUT_TABLE.parent)
    rows_out: List[Dict[str, Any]] = []

    # Try import phonopy
    try:
        from phonopy import Phonopy
        from phonopy.structure.atoms import PhonopyAtoms
    except Exception:
        for cid in cids:
            rows_out.append({
                "id": cid,
                "n_displacements": None,
                "supercell": PH_SUPERCELL,
                "source_structure": None,
                "status": "NO_PHONOPY",
                "note": "Install phonopy in python env to generate displacements.",
                "collected_at": now(),
            })
        write_csv(OUT_TABLE, rows_out, ["id","n_displacements","supercell","source_structure","status","note","collected_at"])
        print(f"[WARN] phonopy not available. wrote {OUT_TABLE} with NO_PHONOPY rows.")
        return

    from pymatgen.core import Structure

    scmat = _parse_supercell(PH_SUPERCELL)

    for cid in cids:
        spath, src = pick_source_structure(cid)
        potcar = pick_potcar(cid)
        out_root = RUNS_ROOT/cid/"phonon_force"
        ensure_dir(out_root)

        if spath is None or potcar is None:
            rows_out.append({
                "id": cid,
                "n_displacements": None,
                "supercell": PH_SUPERCELL,
                "source_structure": src or None,
                "status": "SKIP_MISSING_INPUTS",
                "note": "Missing structure or POTCAR.",
                "collected_at": now(),
            })
            continue

        s = Structure.from_file(str(spath))
        # Build PhonopyAtoms
        pa = PhonopyAtoms(
            symbols=[str(sp.specie) for sp in s.sites],
            cell=s.lattice.matrix,
            scaled_positions=s.frac_coords,
        )
        phonon = Phonopy(pa, scmat)
        phonon.generate_displacements(distance=PH_DISP)
        supercells = phonon.supercells_with_displacements

        # Write phonopy yaml (optional)
        try:
            phonon.save(out_root/"phonopy.yaml")
            phonon.save(out_root/"phonopy_disp.yaml")
        except Exception:
            pass

        # Prepare each displacement dir
        ndisp = 0
        for i, sc in enumerate(supercells, start=1):
            if sc is None:
                continue
            ndisp += 1
            dname = f"disp-{i:03d}"
            ddir = out_root/dname
            ensure_dir(ddir)

            # Write POSCAR in VASP format
            # Convert back to pymatgen Structure for easy POSCAR writing
            lattice = sc.cell
            frac = sc.scaled_positions
            symbols = sc.symbols
            ss = Structure(lattice, symbols, frac)
            ss.to(fmt="poscar", filename=str(ddir/"POSCAR"))

            # Copy POTCAR
            (ddir/"POTCAR").write_bytes(potcar.read_bytes())

            # Ensure no KPOINTS
            kp = ddir/"KPOINTS"
            if kp.exists():
                kp.unlink()

            write_incar(ddir/"INCAR")
            write_pbs(ddir/"job_phonon_force.pbs", f"{cid[:20]}_ph{i:03d}"[:36])

        rows_out.append({
            "id": cid,
            "n_displacements": ndisp,
            "supercell": PH_SUPERCELL,
            "source_structure": src,
            "status": "PREPARED",
            "note": f"disp={PH_DISP}A",
            "collected_at": now(),
        })

    write_csv(OUT_TABLE, rows_out, ["id","n_displacements","supercell","source_structure","status","note","collected_at"])
    log("INFO", f"Wrote phonon job table: {OUT_TABLE} ({len(rows_out)} rows)")

if __name__ == "__main__":
    main()
