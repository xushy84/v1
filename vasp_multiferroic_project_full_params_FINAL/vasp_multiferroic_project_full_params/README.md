# VASP Multiferroic Screening Project (Params updated)

Put CIFs into `data/cifs/` (one structure per CIF).

Key fixes:
- NO KPOINTS files.
- Use KSPACING + KGAMMA=.FALSE. (avoids IBZKPT issues in VASP 5.4.1).
- NCORE=4 in all steps.

Steps with tuned params:
- relax0: ISIF=2, EDIFFG=-0.03, KSPACING=0.30
- mag_*: ISIF=2, EDIFFG=-0.02, KSPACING=0.25
- static_gap: NSW=0, ISMEAR=0, KSPACING=0.22
- polar_kick: short relax, KSPACING=0.25
- relax_tight: ISIF=3, EDIFF=1E-6, EDIFFG=-0.01, KSPACING=0.22

Submit scripts:
- 11_submit_relax0.sh
- 21_submit_mag_relax.sh (set STEP_DIR=mag_FM / mag_AFM1 / mag_LOW)
- 41_submit_static_gap.sh
- 52_submit_polar_kick.sh
- 61_submit_relax_tight.sh


## Added in v4
- 70/71/72: phonon finite displacement (phonopy) prepare/submit/analyze
- 80/81/82: Berry-phase polarization prepare/submit/collect


## Verbose output
Set `VERBOSE=1` to print detailed logs for most Python scripts, e.g.

    VERBOSE=1 python scripts/60_prepare_relax_tight.py

Bash submit scripts already print detailed stats and write logs.
