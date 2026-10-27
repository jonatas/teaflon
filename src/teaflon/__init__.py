import os
import sys
import subprocess


def main() -> None:
    root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
    scripts = [
        "find_templates.py",
        "measure.py",
        "orca_dna_prep.py",
        "digital_twin_defluorination.py",
        "pg_bio_multiomics_pipeline/step1_genomics.py",
        "pg_bio_multiomics_pipeline/step2_folding.py",
        "pg_bio_multiomics_pipeline/step3_transcriptomics.py",
        "pg_bio_multiomics_pipeline/step4_visual_reporter.py",
        "pg_bio_multiomics_pipeline/step5_proteomics_ms.py",
        "pg_bio_multiomics_pipeline/step6_pan_organism_wetlab.py",
        "pg_bio_multiomics_pipeline/test_pipeline.py",
    ]
    for s in scripts:
        path = os.path.join(root, s)
        res = subprocess.run([sys.executable, path], cwd=root)
        if res.returncode != 0:
            sys.exit(res.returncode)

