#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations
from pymatgen.core import Structure

def sorted_structure_for_vasp(struct: Structure) -> Structure:
    sites = sorted(struct.sites, key=lambda s: (s.specie.symbol, *s.frac_coords))
    return Structure(struct.lattice, [s.specie for s in sites], [s.frac_coords for s in sites])
