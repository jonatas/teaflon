# /// script
# requires-python = ">=3.10, <3.13"
# dependencies = [
#     "pymol-open-source-whl",
# ]
# ///

import os
os.environ["PYOPENGL_PLATFORM"] = "osmesa"

import pymol
pymol.pymol_argv = ["pymol", "-cq"]
pymol.finish_launching()

from pymol import cmd

cmd.load("structures/AF-Q1JU72-F1-model_v6.cif", "deha")

# Get distances
n_term_dist = cmd.distance("dist_n", "resi 104+271+128", "resi 1 and name CA")
c_term_dist = cmd.distance("dist_c", "resi 104+271+128", "resi 304 and name CA")

print(f"[*] Distance from Active Site to N-terminus: {n_term_dist:.2f} Angstroms")
print(f"[*] Distance from Active Site to C-terminus: {c_term_dist:.2f} Angstroms")

cmd.quit()
