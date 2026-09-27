# /// script
# requires-python = ">=3.10, <3.13"
# dependencies = [
#     "pymol-open-source-whl",
# ]
# ///

import os
import sys

# Set environment variable for headless rendering
os.environ["PYOPENGL_PLATFORM"] = "osmesa"

import pymol # pytype: disable=import-error
pymol.pymol_argv = ["pymol", "-cq"]
pymol.finish_launching()

from pymol import cmd # pytype: disable=import-error

# Load the two domains
cmd.load("dehalogenase.pdb", "destroyer")
if cmd.count_atoms("destroyer") == 0:
    print("Error loading dehalogenase.pdb")
    cmd.quit()

cmd.load("hydrophobin.pdb", "hook")
if cmd.count_atoms("hook") == 0:
    print("Error loading hydrophobin.pdb")
    cmd.quit()

# Translate the hook 50 Angstroms to the right to simulate the flexible linker separation
cmd.translate([50, 0, 0], "hook")

# Styling
cmd.bg_color("white")
cmd.show_as("surface", "all")
cmd.set("transparency", 0.3)
cmd.show("cartoon", "all")

# Color Destroyer (Dehalogenase) red/orange
cmd.color("firebrick", "destroyer")
# Color Hook (Hydrophobin) blue
cmd.color("marine", "hook")

# Draw a dashed line to represent the GGGGSGGGGS linker
# Find C-term of destroyer and N-term of hook
cmd.select("c_term", "destroyer and name C and resi 298")
cmd.select("n_term", "hook and name N and resi 1")
cmd.distance("linker", "c_term", "n_term")
cmd.set("dash_color", "black")
cmd.set("dash_width", 3)
cmd.set("dash_gap", 0.5)

cmd.orient()
cmd.zoom("all", buffer=5)

cmd.set("ray_opaque_background", 1)
cmd.png("./output/teaflon_fusion_concept.png", width=1200, height=800, dpi=150)
cmd.save("teaflon_fusion.pse")
cmd.quit()
print("Successfully rendered fusion concept to ./output/teaflon_fusion_concept.png")
sys.exit(0)
