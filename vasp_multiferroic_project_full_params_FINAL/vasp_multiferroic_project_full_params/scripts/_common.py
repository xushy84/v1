#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations
from pathlib import Path
import os
from typing import Optional, Dict, Any, List
import re, csv
from datetime import datetime

def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def now() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def read_text_safe(p: Path) -> str:
    try:
        return p.read_text(errors="ignore")
    except FileNotFoundError:
        return ""

def outcar_is_complete(outcar: Path) -> bool:
    return "General timing and accounting" in read_text_safe(outcar)

def parse_last_toten(outcar: Path) -> Optional[float]:
    txt = read_text_safe(outcar)
    m = list(re.finditer(r"free\s+energy\s+TOTEN\s*=\s*([-\d\.Ee+]+)", txt))
    if not m:
        return None
    try:
        return float(m[-1].group(1))
    except Exception:
        return None

def parse_magmom_total(outcar: Path) -> Optional[float]:
    txt = read_text_safe(outcar)
    blocks = txt.split("magnetization (x)")
    if len(blocks) < 2:
        return None
    tail = blocks[-1]
    tot_lines = [ln for ln in tail.splitlines() if re.search(r"\btot\b", ln)]
    if not tot_lines:
        return None
    nums = re.findall(r"[-\d\.Ee+]+", tot_lines[-1])
    if not nums:
        return None
    try:
        return float(nums[-1])
    except Exception:
        return None

def read_structure_volume_from_contcar(contcar: Path) -> Optional[float]:
    try:
        from pymatgen.core import Structure
        s = Structure.from_file(str(contcar))
        return float(s.volume)
    except Exception:
        return None

def write_csv(path: Path, rows: List[Dict[str, Any]], fieldnames: List[str]) -> None:
    ensure_dir(path.parent)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k) for k in fieldnames})

def load_csv(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", newline="") as f:
        return list(csv.DictReader(f))

def pick_ids_from_table(table_path: Path, flag_col: str = "keep", flag_yes: str = "Y") -> List[str]:
    rows = load_csv(table_path)
    return [r["id"] for r in rows if r.get(flag_col,"").strip().upper()==flag_yes]

# =========================
# Logging helpers
# =========================
import sys
from datetime import datetime as _dt

def log(level: str, msg: str) -> None:
    ts = _dt.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{ts}] [{level.upper()}] {msg}", file=sys.stdout, flush=True)

def vlog(enabled: bool, level: str, msg: str) -> None:
    if enabled:
        log(level, msg)

def get_verbose() -> bool:
    v = str(os.environ.get("VERBOSE", "0")).strip().lower()
    return v in ("1","true","yes","y","on")
