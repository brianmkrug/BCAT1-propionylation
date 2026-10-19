#!/usr/bin/env python3
"""Superimpose the BCAT2 AlphaFold model onto BCAT1 using structurally
equivalent CA atoms, and report the fit."""
import numpy as np
import biotite.structure as struc
import biotite.structure.io.pdb as pdb

OUT = "/root/bc/data"


def load(acc):
    f = pdb.PDBFile.read(f"{OUT}/AF_{acc}.pdb")
    a = f.get_structure(model=1)
    return a[struc.filter_amino_acids(a)]


ref = load("P54687")   # BCAT1
mob = load("O15382")   # BCAT2

fitted, transform, ref_idx, mob_idx = struc.superimpose_homologs(
    ref, mob, max_iterations=10)

ca_r = ref[ref_idx][ref[ref_idx].atom_name == "CA"]
ca_m = fitted[mob_idx][fitted[mob_idx].atom_name == "CA"]
rmsd = struc.rmsd(ca_r.coord, ca_m.coord)
print(f"superimposed on {len(ca_r)} equivalent CA atoms, RMSD = {rmsd:.2f} A")
print("BCAT1 residues used:", ref[ref_idx].res_id.min(), "-", ref[ref_idx].res_id.max())
print("BCAT2 residues used:", mob[mob_idx].res_id.min(), "-", mob[mob_idx].res_id.max())

f = pdb.PDBFile()
f.set_structure(fitted)
f.write(f"{OUT}/AF_O15382_on_BCAT1.pdb")

# per-residue CA deviation for BCAT2 residues that have an equivalent
dev = np.linalg.norm(ca_r.coord - ca_m.coord, axis=1)
np.save(f"{OUT}/ca_dev.npy", np.vstack([ca_r.res_id, ca_m.res_id, dev]))
print("median per-residue CA deviation %.2f A; 90th pct %.2f A"
      % (np.median(dev), np.percentile(dev, 90)))
