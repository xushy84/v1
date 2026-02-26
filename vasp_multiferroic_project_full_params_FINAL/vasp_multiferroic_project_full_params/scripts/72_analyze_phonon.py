#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""72_analyze_phonon.py

Collect phonon forces from runs/<id>/phonon_force/disp-*/OUTCAR and run phonopy to compute phonons.

Outputs:
  results/tables/phonon_summary.csv
    id, status, imag_mode_flag, min_freq_THz, n_displacements, note, collected_at

Requirements:
  - phonopy installed.
  - 70_prepare_phonon.py must have generated displacements (phonon_force/disp-*/POSCAR).
  - OUTCAR must be complete and contain forces.

Env:
  PH_FORCE_SCALE (default 1.0)  # multiply forces if needed
"""

from __future__ import annotations

import os

import os
from pathlib import Path
from typing import Dict, Any, List, Tuple

from _common import ensure_dir, write_csv, now, outcar_is_complete, get_verbose, log, vlog

RUNS_ROOT = Path(os.environ.get("RUNS_ROOT", "runs"))
OUT_TABLE = Path("results/tables/phonon_summary.csv")

PH_FORCE_SCALE = float(os.environ.get("PH_FORCE_SCALE", "1.0"))

def parse_forces_from_outcar(outcar: Path) -> List[List[float]] | None:
    # Extract the last "TOTAL-FORCE (eV/Angst)" block.
    txt = outcar.read_text(errors="ignore")
    key = "TOTAL-FORCE (eV/Angst)"
    idx = txt.rfind(key)
    if idx < 0:
        return None
    block = txt[idx:].splitlines()
    forces = []
    # Skip header lines until dashed line then read coordinates+forces lines
    started = False
    for ln in block:
        if "----" in ln:
            started = True
            continue
        if not started:
            continue
        parts = ln.split()
        if len(parts) < 6:
            # end of block
            if forces:
                break
            continue
        try:
            fx, fy, fz = float(parts[3]), float(parts[4]), float(parts[5])
            forces.append([PH_FORCE_SCALE*fx, PH_FORCE_SCALE*fy, PH_FORCE_SCALE*fz])
        except Exception:
            if forces:
                break
            continue
    return forces if forces else None


VERBOSE = get_verbose()

def main():
    ensure_dir(OUT_TABLE.parent)
    rows_out: List[Dict[str, Any]] = []

    try:
        import phonopy
        from phonopy import Phonopy
        from phonopy.structure.atoms import PhonopyAtoms
        from phonopy.file_IO import parse_FORCE_SETS
    except Exception:
        rows_out.append({
            "id": "*",
            "status": "NO_PHONOPY",
            "imag_mode_flag": "UNK",
            "min_freq_THz": None,
            "n_displacements": None,
            "note": "Install phonopy in python env.",
            "collected_at": now(),
        })
        write_csv(OUT_TABLE, rows_out, ["id","status","imag_mode_flag","min_freq_THz","n_displacements","note","collected_at"])
        print(f"[WARN] phonopy not available. wrote {OUT_TABLE}.")
        return

    from pymatgen.core import Structure

    for jobdir in sorted(RUNS_ROOT.glob("*")):
        cid = jobdir.name
        root = jobdir/"phonon_force"
        if not root.is_dir():
            continue

        disp_dirs = sorted(root.glob("disp-*"))
        if not disp_dirs:
            rows_out.append({"id": cid, "status": "NO_DISP", "imag_mode_flag":"UNK", "min_freq_THz":None,
                             "n_displacements":0, "note":"No disp-* dirs", "collected_at": now()})
            continue

        # Try to build forces list aligned with disp order
        forces_all = []
        ok = True
        for d in disp_dirs:
            outcar = d/"OUTCAR"
            if not outcar_is_complete(outcar):
                ok = False
                break
            f = parse_forces_from_outcar(outcar)
            if f is None:
                ok = False
                break
            forces_all.append(f)

        if not ok:
            rows_out.append({"id": cid, "status": "INCOMPLETE", "imag_mode_flag":"UNK", "min_freq_THz":None,
                             "n_displacements":len(disp_dirs), "note":"missing OUTCAR/forces", "collected_at": now()})
            continue

        # Read the undeformed supercell POSCAR from first disp as reference
        # (phonopy itself should know the reference from yaml, but keep it robust)
        base_poscar = disp_dirs[0]/"POSCAR"
        try:
            s = Structure.from_file(str(base_poscar))
            pa = PhonopyAtoms(
                symbols=[str(sp.specie) for sp in s.sites],
                cell=s.lattice.matrix,
                scaled_positions=s.frac_coords,
            )
        except Exception:
            rows_out.append({"id": cid, "status":"BAD_POSCAR", "imag_mode_flag":"UNK", "min_freq_THz":None,
                             "n_displacements":len(disp_dirs), "note":"cannot read POSCAR", "collected_at": now()})
            continue

        # Load phonopy object from saved yaml if possible (preferred)
        phonon = None
        try:
            from phonopy import load
            y = root/"phonopy.yaml"
            yd = root/"phonopy_disp.yaml"
            if y.exists():
                phonon = load(str(y))
            elif yd.exists():
                phonon = load(str(yd))
        except Exception:
            phonon = None

        if phonon is None:
            rows_out.append({"id": cid, "status":"NO_YAML", "imag_mode_flag":"UNK", "min_freq_THz":None,
                             "n_displacements":len(disp_dirs), "note":"missing phonopy.yaml", "collected_at": now()})
            continue

        # Create FORCE_SETS in-memory
        # phonopy expects list of forces corresponding to supercells_with_displacements.
        try:
            phonon.forces = forces_all
            phonon.produce_force_constants()
            # Compute frequencies at Gamma only for quick stability flag
            mesh_str = os.environ.get("PH_MESH", "1 1 1")
            mesh = [int(x) for x in mesh_str.split()]
            phonon.run_mesh(mesh, with_eigenvectors=False, is_gamma_center=True)
            mesh = phonon.mesh
            freqs = mesh.frequencies  # shape (nq, nmodes)
            f0 = freqs[0]
            min_freq = float(min(f0))
            imag = "Y" if min_freq < -1e-3 else "N"
            rows_out.append({
                "id": cid,
                "status": "OK",
                "imag_mode_flag": imag,
                "min_freq_THz": min_freq,
                "n_displacements": len(disp_dirs),
                "note": f"mesh={os.environ.get('PH_MESH','1 1 1')} quick check",
                "collected_at": now(),
            })
        except Exception as e:
            rows_out.append({
                "id": cid,
                "status": "FAIL",
                "imag_mode_flag": "UNK",
                "min_freq_THz": None,
                "n_displacements": len(disp_dirs),
                "note": f"phonopy failed: {type(e).__name__}",
                "collected_at": now(),
            })

    write_csv(OUT_TABLE, rows_out, ["id","status","imag_mode_flag","min_freq_THz","n_displacements","note","collected_at"])
    print(f"[OK] wrote {OUT_TABLE} ({len(rows_out)} rows)")

if __name__ == "__main__":
    main()
