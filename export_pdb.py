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

cmd.load('dehalogenase.pdb', 'destroyer')
cmd.load('hydrophobin.pdb', 'hook')
cmd.translate([50, 0, 0], 'hook')

os.makedirs('./output', exist_ok=True)
cmd.save('./output/teaflon_fusion.pdb', 'all')
cmd.quit()
sys.exit(0)
