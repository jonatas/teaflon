# /// script
# requires-python = ">=3.10, <3.13"
# dependencies = [
#     "pymol-open-source-whl",
# ]
# ///

import os
import sys

os.environ["PYOPENGL_PLATFORM"] = "osmesa"

import pymol
pymol.pymol_argv = ["pymol", "-cq"]
pymol.finish_launching()

from pymol import cmd

# Load Hydrophobin (The Teflon Hook)
cmd.load("hydrophobin.pdb", "hook")
# Load PTFE Oligomer (Perfluorohexane)
cmd.load("ptfe_oligomer.sdf", "ptfe")

# Center and visualize
cmd.bg_color("white")

# Show Hydrophobin as surface showing electrostatics/hydrophobicity roughly
# We'll just color it by hydrophobicity (using a custom spectrum on resv, or just a solid color)
cmd.show_as("surface", "hook")
cmd.color("palegreen", "hook")
cmd.set("transparency", 0.2, "hook")
cmd.show("cartoon", "hook")

# Show PTFE as spheres
cmd.show_as("spheres", "ptfe")
# Color Carbon cyan, Fluorine light pink (or standard CPK)
cmd.color("cyan", "ptfe and elem C")
cmd.color("pink", "ptfe and elem F")

# Manually translate PTFE to sit on the surface of the hydrophobin to simulate docking
# (In a real pipeline, AutoDock Vina would output the coordinates)
cmd.translate([15, 5, 0], "ptfe")

cmd.orient()
cmd.zoom("all", buffer=3)

cmd.set("ray_opaque_background", 1)
cmd.png("./output/teaflon_docking_concept.png", width=1200, height=800, dpi=150)
cmd.save("teaflon_docking.pse")
cmd.quit()
print("Saved docking render to ./output/teaflon_docking_concept.png")
sys.exit(0)
