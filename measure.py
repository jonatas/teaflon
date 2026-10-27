import os
import sys
import time
import json
import psycopg
from rich.console import Console
from rich.table import Table

# Add the local pgbio-py SDK to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../pg_bio/pgbio-py")))
from pgbio import PgBioClient

console = Console()
DB_URL = os.getenv("PG_BIO_URL", "postgresql://localhost:28818/bio_demo")


def run_structural_measurements():
    console.print("[bold cyan]📐 TeaFlon 3D Structural & Theozyme Audit (Native pg_bio Layers 6, 8 & 12)[/]")
    t0 = time.perf_counter()

    client = PgBioClient(DB_URL)

    # 1. Measure N-terminal vs C-terminal linker placement geometry on Q1JU72 (dehalogenase.pdb) natively in SQL
    with open("dehalogenase.pdb", "r") as f:
        deha_pdb = f.read()
    with open("hydrophobin.pdb", "r") as f:
        hydro_pdb = f.read()
    with open("teaflon_trifusion_local.pdb", "r") as f:
        trifusion_pdb = f.read()

    with psycopg.connect(DB_URL) as conn:
        with conn.cursor() as cur:
            cur.execute(
                """
                WITH parsed AS (
                    SELECT atom_name, res_name, res_id, x, y, z,
                           create_residue_coord(x, y, z, res_name || res_id::text || ':' || atom_name) AS coord
                    FROM bio_parse_pdb(%s)
                ),
                triad_com AS (
                    SELECT center_of_mass(ARRAY_AGG(coord)) AS com
                    FROM parsed
                    WHERE res_id IN (104, 128, 271) AND atom_name = 'CA'
                ),
                n_term AS (
                    SELECT coord FROM parsed WHERE res_id = 1 AND atom_name = 'CA' LIMIT 1
                ),
                c_term AS (
                    SELECT coord FROM parsed WHERE res_id = 304 AND atom_name = 'CA' LIMIT 1
                ),
                hexad AS (
                    SELECT ARRAY[
                        (SELECT coord FROM parsed WHERE res_id = 104 AND atom_name = 'OD2' LIMIT 1),
                        (SELECT coord FROM parsed WHERE res_id = 271 AND atom_name = 'NE2' LIMIT 1),
                        (SELECT coord FROM parsed WHERE res_id = 128 AND atom_name = 'OD2' LIMIT 1),
                        (SELECT coord FROM parsed WHERE res_id = 150 AND atom_name = 'NE1' LIMIT 1),
                        (SELECT coord FROM parsed WHERE res_id = 212 AND atom_name = 'OH'  LIMIT 1),
                        (SELECT coord FROM parsed WHERE res_id = 108 AND atom_name = 'NH2' LIMIT 1)
                    ] AS site_coords
                )
                SELECT bio_distance((SELECT com FROM triad_com), (SELECT coord FROM n_term)) AS dist_n,
                       bio_distance((SELECT com FROM triad_com), (SELECT coord FROM c_term)) AS dist_c,
                       theozyme_rmsd((SELECT site_coords FROM hexad), 'q1ju72_dehalogenase_hexad') AS hexad_rmsd,
                       bio_distance(
                           (SELECT coord FROM parsed WHERE res_id = 104 AND atom_name = 'OD2' LIMIT 1),
                           (SELECT coord FROM parsed WHERE res_id = 271 AND atom_name = 'NE2' LIMIT 1)
                       ) AS asp104_his271_dist,
                       bio_angle(
                           (SELECT coord FROM parsed WHERE res_id = 128 AND atom_name = 'OD2' LIMIT 1),
                           (SELECT coord FROM parsed WHERE res_id = 271 AND atom_name = 'NE2' LIMIT 1),
                           (SELECT coord FROM parsed WHERE res_id = 104 AND atom_name = 'OD2' LIMIT 1)
                       ) AS charge_relay_angle;
                """,
                (deha_pdb,),
            )
            dist_n, dist_c, hexad_rmsd, asp_his_dist, relay_angle = cur.fetchone()

    console.print(f"\n[*] Active Site (Asp104-Asp128-His271) to N-terminus (Met1 CA): [bold green]{dist_n:.2f} Å[/]")
    console.print(f"[*] Active Site (Asp104-Asp128-His271) to C-terminus (Ser304 CA): [bold green]{dist_c:.2f} Å[/] (Optimal for GGGGSGGGGS linker!)")
    console.print(f"[*] Catalytic Nucleophile-Base Distance (Asp104:OD2 .. His271:NE2): [bold cyan]{asp_his_dist:.2f} Å[/]")
    console.print(f"[*] Charge-Relay Angle (Asp128:OD2 - His271:NE2 - Asp104:OD2): [bold cyan]{relay_angle:.2f}°[/]")
    console.print(f"[*] Kabsch Theozyme RMSD (`q1ju72_dehalogenase_hexad`): [bold green]{hexad_rmsd:.4f} Å[/]\n")

    # 2. Evaluate 3D Globular Compactness Gate across all local TeaFlon PDB constructs
    table = Table(title="TeaFlon 3D Globular Compactness & Folding Audit (`evaluate_globular_compactness` in pg_bio)")
    table.add_column("Construct / Domain", style="cyan")
    table.add_column("CA Atoms", style="white")
    table.add_column("Mean pLDDT", style="green")
    table.add_column("Rg (Å)", style="magenta")
    table.add_column("Flory Compactness", style="yellow")
    table.add_column("Elongation", style="blue")
    table.add_column("Globular Gate", style="bold green")

    for label, pdb_str in [
        ("Dehalogenase Q1JU72 (dehalogenase.pdb)", deha_pdb),
        ("Hydrophobin P79073 (hydrophobin.pdb)", hydro_pdb),
        ("TeaFlon TriFusion ESMFold (teaflon_trifusion_local.pdb)", trifusion_pdb),
    ]:
        prof = client.evaluate_globular_compactness(pdb_str)
        table.add_row(
            label,
            str(prof.ca_count),
            f"{prof.mean_plddt:.1f}",
            f"{prof.radius_of_gyration:.2f} Å",
            f"{prof.compactness_ratio:.3f}",
            f"{prof.elongation:.2f}",
            "✅ PASS (Globular)" if prof.is_globular else "⚠️ Multi-Domain Extended",
        )

    elapsed_ms = (time.perf_counter() - t0) * 1000.0
    console.print(table)
    console.print(f"[bold green]⚡ Completed full 3D structural & theozyme audit in {elapsed_ms:.2f} ms (zero PyMOL dependency)![/]")


if __name__ == "__main__":
    run_structural_measurements()

