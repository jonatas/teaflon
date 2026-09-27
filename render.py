# /// script
# requires-python = ">=3.10, <3.13"
# dependencies = [
#     "pymol-open-source-whl",
# ]
# ///

import os
import sys

os.environ["PYOPENGL_PLATFORM"] = "osmesa"

import pymol # pytype: disable=import-error
pymol.pymol_argv = ["pymol", "-cq"]
pymol.finish_launching()

from pymol import cmd # pytype: disable=import-error

# Load the structures
cmd.load("structures/AF-Q1JU72-F1-model_v6.cif", "dehalogenase")
cmd.load("structures/AF-P79073-F1-model_v6.cif", "hydrophobin")

# Move them apart so they don't overlap
cmd.translate([50, 0, 0], "hydrophobin")

# Setup views
cmd.show_as("cartoon")

# Color by pLDDT confidence (b-factor column in AF2 CIFs)
# AlphaFold pLDDT colors:
# >90 Dark blue, >70 Light blue, >50 Yellow, <50 Orange/Red
cmd.color("blue", "b > 90")
cmd.color("cyan", "b < 90 and b > 70")
cmd.color("yellow", "b < 70 and b > 50")
cmd.color("orange", "b < 50")

# Center and setup rendering
cmd.zoom("all")
cmd.set("bg_rgb", "white")
cmd.set("ray_opaque_background", 1)

# Render and save
cmd.png("structures/visualization.png", width=1200, height=800, dpi=150)
cmd.save("structures/session.pse")

cmd.quit()
