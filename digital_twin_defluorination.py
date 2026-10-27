import os
import sys
import time
import numpy as np
import matplotlib.pyplot as plt
import psycopg
from rich.console import Console
from rich.table import Table

console = Console()
DB_URL = os.getenv("PG_BIO_URL", "postgresql://localhost:28818/bio_demo")


def simulate_ise_fluoride_release():
    console.print("[bold cyan]🧪 TEAFLON DIGITAL TWIN: Step 1 Defluorination & Fluorapatite Biomineralization Test[/]")
    console.print("=" * 78)
    console.print("Simulating Ion-Selective Electrode (ISE) [F⁻] & Fluorapatite Ca₅(PO₄)₃F yield natively in pg_bio (Layer 12).")

    t0 = time.perf_counter()
    with psycopg.connect(DB_URL) as conn:
        with conn.cursor() as cur:
            # Generate 49 time points (0 to 48 hours) across 4 enzymatic regimes natively in PostgreSQL
            cur.execute(
                """
                WITH grid AS (
                    SELECT g::float8 AS t_hr FROM generate_series(0, 48, 1) AS g
                )
                SELECT g.t_hr,
                       0.0::float8 AS f_ctrl,
                       wt.fluoride_ppm AS f_wt,
                       tf.fluoride_ppm AS f_teaflon,
                       tf.fluorapatite_mg_per_l AS fap_teaflon,
                       ex.fluoride_ppm AS f_extremophile,
                       ex.fluorapatite_mg_per_l AS fap_extremophile
                FROM grid g
                CROSS JOIN LATERAL simulate_defluorination_kinetics(1.5, 450.0, 200.0, 0.25, 10.0, g.t_hr) wt
                CROSS JOIN LATERAL simulate_defluorination_kinetics(12.0, 45.0, 1000.0, 5.0, 45.0, g.t_hr) tf
                CROSS JOIN LATERAL simulate_defluorination_kinetics(18.5, 32.0, 1000.0, 7.5, 50.0, g.t_hr) ex
                ORDER BY g.t_hr ASC;
                """
            )
            rows = cur.fetchall()

    elapsed_ms = (time.perf_counter() - t0) * 1000.0

    time_hours = np.array([r[0] for r in rows])
    fluoride_control = np.array([r[1] for r in rows])
    fluoride_natural = np.array([r[2] for r in rows])
    fluoride_teaflon = np.array([r[3] for r in rows])
    fap_teaflon = np.array([r[4] for r in rows])
    fluoride_extremo = np.array([r[5] for r in rows])
    fap_extremo = np.array([r[6] for r in rows])

    table = Table(title=f"In-Database ISE Fluoride & Fluorapatite Digital Twin ({elapsed_ms:.2f} ms in pg_bio)")
    table.add_column("Incubation Time", style="cyan")
    table.add_column("Wild-Type [F⁻] (ppm)", style="blue")
    table.add_column("TeaFlon v2 TriFusion [F⁻] (ppm)", style="red")
    table.add_column("TeaFlon Fluorapatite (mg/L)", style="magenta")
    table.add_column("Extremophile Swarm Chimera [F⁻] (ppm)", style="green")
    table.add_column("Extremophile Fluorapatite (mg/L)", style="bold green")

    for idx in (0, 6, 12, 24, 36, 48):
        r = rows[idx]
        table.add_row(
            f"{int(r[0]):2d} h",
            f"{r[2]:6.2f} ppm",
            f"{r[3]:6.2f} ppm",
            f"{r[4]:7.1f} mg/L",
            f"{r[5]:6.2f} ppm",
            f"{r[6]:7.1f} mg/L",
        )
    console.print(table)

    # Plotting the Digital Twin Sensor Data
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

    ax1.plot(time_hours, fluoride_control, label="Control (Buffer Only)", color="gray", linestyle="--")
    ax1.plot(time_hours, fluoride_natural, label="Wild-type Dehalogenase (No Hook)", color="blue")
    ax1.plot(time_hours, fluoride_teaflon, label="TeaFlon v2 TriFusion (Q1JU72+P79073+AMELX)", color="red", linewidth=2.5)
    ax1.plot(time_hours, fluoride_extremo, label="Extremophile Swarm Chimera (A1RTA7+Q7S3P5+AMELX)", color="darkgreen", linewidth=2.5, linestyle="-.")
    ax1.fill_between(time_hours, fluoride_teaflon, alpha=0.12, color="red")
    ax1.set_title("Fluoride ISE Sensor Digital Twin ([F⁻] Release)", pad=12)
    ax1.set_xlabel("Incubation Time (Hours)")
    ax1.set_ylabel("Released Fluoride Concentration [F⁻] (ppm)")
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=8.5)

    ax2.plot(time_hours, fap_teaflon, label="TeaFlon v2 Fluorapatite Ca₅(PO₄)₃F", color="crimson", linewidth=2.5)
    ax2.plot(time_hours, fap_extremo, label="Extremophile Chimera Fluorapatite Ca₅(PO₄)₃F", color="teal", linewidth=2.5, linestyle="-.")
    ax2.fill_between(time_hours, fap_extremo, alpha=0.12, color="teal")
    ax2.set_title("Stoichiometric Bioceramic Precipitation (AMELX Tag)", pad=12)
    ax2.set_xlabel("Incubation Time (Hours)")
    ax2.set_ylabel("Precipitated Fluorapatite Mass (mg/L)")
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=8.5)

    plt.tight_layout()
    output_path = "assets/step1_ise_simulation.png"
    os.makedirs("assets", exist_ok=True)
    os.makedirs("output", exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    console.print(f"\n[bold green]✅ Simulation Complete.[/] ISE & Fluorapatite chart saved to: [cyan]{output_path}[/]")


if __name__ == "__main__":
    simulate_ise_fluoride_release()

