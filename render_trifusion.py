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

# 1. Load Dehalogenase (Destroyer)
cmd.load("dehalogenase.pdb", "destroyer")

# 2. Load Hydrophobin (Hook)
cmd.load("hydrophobin.pdb", "hook")
cmd.translate([50, 0, 0], "hook")

# 3. Build AMELX Biomineralization Tag (Tooth Builder)
cmd.fab("WPSTDKTKREEVD", "amelx_tag")
cmd.translate([80, 0, 0], "amelx_tag")

# Styling
cmd.bg_color("white")
cmd.show_as("surface", "destroyer")
cmd.show_as("surface", "hook")
cmd.show_as("sticks", "amelx_tag")
cmd.show("cartoon", "all")
cmd.set("transparency", 0.3)

# Coloring
cmd.color("firebrick", "destroyer")
cmd.color("marine", "hook")
cmd.color("purple", "amelx_tag")

# Connect them with visual dashed lines to simulate the linkers
cmd.select("c_term_destroyer", "destroyer and name C and resi 298")
cmd.select("n_term_hook", "hook and name N and resi 1")
cmd.distance("linker1", "c_term_destroyer", "n_term_hook")

cmd.select("c_term_hook", "hook and name C and resi 93")
cmd.select("n_term_amelx", "amelx_tag and name N and resi 1")
cmd.distance("linker2", "c_term_hook", "n_term_amelx")

cmd.set("dash_color", "black")
cmd.set("dash_width", 3)
cmd.set("dash_gap", 0.5)
cmd.hide("labels")

cmd.orient()
cmd.zoom("all", buffer=5)

cmd.set("ray_opaque_background", 1)
cmd.png("./output/teaflon_trifusion_concept.png", width=1200, height=800, dpi=150)
cmd.save("teaflon_trifusion.pse")
cmd.quit()
print("Saved TriFusion render to ./output/teaflon_trifusion_concept.png")
sys.exit(0)
