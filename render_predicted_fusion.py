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

cmd.load("teaflon_fusion.pdb", "fusion")

# Color by pLDDT (B-factor column in ESMFold outputs)
# > 90 Blue, 70-90 Cyan, 50-70 Yellow, <50 Red
cmd.spectrum("b", "red_yellow_cyan_blue", "fusion", minimum=50, maximum=90)

cmd.bg_color("white")
cmd.show_as("cartoon", "fusion")
cmd.set("ray_opaque_background", 1)

cmd.orient()
cmd.zoom("all", buffer=5)
cmd.png("./output/teaflon_predicted_fusion.png", width=1200, height=800, dpi=150)
cmd.quit()
print("Saved predicted fusion render to ./output/teaflon_predicted_fusion.png")
