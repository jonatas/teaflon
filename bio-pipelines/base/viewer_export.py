import os
import time
import psycopg
from rich.console import Console
from rich.table import Table
from rich import box

console = Console()
DB_URL = os.getenv("PG_BIO_URL", "postgresql://localhost:28818/bio_demo")
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "../output")


def export_merged_pdb_via_pg_bio():
    console.print("[bold cyan]🧊 pg_bio 3D Rigid-Body Domain Translator & Globular Gate (`translate_and_merge_pdb_chains`)[/]")
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    domain_a_path = os.path.join(REPO_ROOT, "dehalogenase.pdb")
    domain_b_path = os.path.join(REPO_ROOT, "hydrophobin.pdb")
    domain_a_pdb = open(domain_a_path).read()
    domain_b_pdb = open(domain_b_path).read()

    t0 = time.perf_counter()
    with psycopg.connect(DB_URL) as conn:
        with conn.cursor() as cur:
            cur.execute(
                """
                WITH merged AS (
                    SELECT translate_and_merge_pdb_chains(%s, %s, 28.0) AS pdb_str
                )
                SELECT m.pdb_str,
                       evaluate_globular_compactness(m.pdb_str)::text
                FROM merged m;
                """,
                (domain_a_pdb, domain_b_pdb),
            )
            pdb_str, glob_json = cur.fetchone()
            import json
            g = json.loads(glob_json)
            ca_count = g["ca_count"]
            mean_plddt = g["mean_plddt"]
            rg = g["radius_of_gyration"]
            compactness = g["compactness_ratio"]
            elongation = g["elongation"]

    elapsed_ms = (time.perf_counter() - t0) * 1000.0
    out_file = os.path.join(OUTPUT_DIR, "merged_domain_complex.pdb")
    with open(out_file, "w") as f:
        f.write(pdb_str)

    table = Table(
        title=f"In-Database 3D Rigid-Body PDB Assembly & Geometry Audit ({elapsed_ms:.2f} ms in pg_bio)",
        box=box.SIMPLE_HEAVY,
    )
    table.add_column("Metric", style="cyan")
    table.add_column("Value", style="bold green")
    table.add_row("Chain A + Chain B Cα Count", f"{ca_count} residues")
    table.add_row("Mean Confidence (pLDDT)", f"{mean_plddt:.2f}")
    table.add_row("Radius of Gyration (Rg)", f"{rg:.2f} Å")
    table.add_row("Flory Compactness Ratio", f"{compactness:.3f}")
    table.add_row("Principal Axis Elongation", f"{elongation:.3f}")
    table.add_row("Exported 3D Complex Path", out_file)
    console.print(table)
    return out_file


if __name__ == "__main__":
    export_merged_pdb_via_pg_bio()
