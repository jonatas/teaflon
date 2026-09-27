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

cmd.load("dehalogenase.pdb", "destroyer")
cmd.alter("destroyer", "chain='A'")

cmd.load("hydrophobin.pdb", "hook")
cmd.translate([50, 0, 0], "hook")
cmd.alter("hook", "chain='B'")

cmd.fab("WPSTDKTKREEVD", "amelx_tag")
cmd.translate([80, 0, 0], "amelx_tag")
cmd.alter("amelx_tag", "chain='C'")

cmd.create("trifusion_assembly", "destroyer or hook or amelx_tag")
cmd.save("teaflon_trifusion_assembly.pdb", "trifusion_assembly")
cmd.quit()
print("Saved teaflon_trifusion_assembly.pdb with distinct chains")
