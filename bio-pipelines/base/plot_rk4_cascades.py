import os
import time
import numpy as np
import matplotlib.pyplot as plt
import psycopg
from rich.console import Console
from rich.table import Table
from rich import box

console = Console()
DB_URL = os.getenv("PG_BIO_URL", "postgresql://localhost:28818/bio_demo")
OUTPUT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "../output"))

PIPELINES = [
    {
        "id": "ptfe-to-bioceramic",
        "title": "1. PTFE -> Fluorapatite Sink (37°C, pH 7.4)",
        "sub_label": "Surface C-F Sites (µM)",
        "int_label": "Free Fluoride [F⁻] (µM)",
        "prod_label": "Fluorapatite Ca₅(PO₄)₃F (µM eq)",
        "kcat1": 24.0,
        "km1": 380.0,
        "kcat2": 42.0,
        "km2": 30.0,
        "ki": 22.0,
        "boost": 8.5,
        "color": "#e63946",
    },
    {
        "id": "plastics-upcycling",
        "title": "2. PET Polyester -> PCA/Vanillin Sink (68°C, pH 8.0)",
        "sub_label": "PET Ester Bonds (µM)",
        "int_label": "MHET / Acidic Monomer (µM)",
        "prod_label": "Upcycled Sink Product (µM)",
        "kcat1": 25.5,
        "km1": 420.0,
        "kcat2": 38.0,
        "km2": 35.0,
        "ki": 18.0,
        "boost": 7.2,
        "color": "#2a9d8f",
    },
    {
        "id": "sugarcane-bagasse",
        "title": "3. Sugarcane Bagasse -> 2,5-FDCA Sink (55°C, pH 5.0)",
        "sub_label": "Cellulose/Xylan Chains (µM)",
        "int_label": "Cellobiose Intermediate (µM)",
        "prod_label": "G1P + 2,5-FDCA Sink (µM)",
        "kcat1": 23.0,
        "km1": 450.0,
        "kcat2": 40.0,
        "km2": 32.0,
        "ki": 15.0,
        "boost": 8.0,
        "color": "#d97706",
    },
    {
        "id": "corn-stover-bran",
        "title": "4. Corn Stover Diferulate -> Bio-Vanillin (50°C, pH 5.5)",
        "sub_label": "Diferulate Ester Bridges (µM)",
        "int_label": "Free Ferulic Acid (µM)",
        "prod_label": "Bio-Vanillin Sink (µM)",
        "kcat1": 24.5,
        "km1": 390.0,
        "kcat2": 41.0,
        "km2": 30.0,
        "ki": 20.0,
        "boost": 7.8,
        "color": "#457b9d",
    },
]


def run_rk4_multi_pipeline_benchmark():
    console.print("[bold cyan]📈 pg_bio Layer 14 RK4 Cascade Benchmark Across All 4 Feedstock Pipelines[/]")
    console.print("=" * 94)
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    t0 = time.perf_counter()
    results = {}

    with psycopg.connect(DB_URL) as conn:
        with conn.cursor() as cur:
            for p in PIPELINES:
                # 1. Full 4-Module Cascade (Module 1 Anchor Boost + Module 2 Scissor + Module 3 Equilibrium Sink)
                cur.execute(
                    """
                    SELECT time_hr, substrate_u_m, intermediate_u_m, sink_product_u_m, fractional_conversion
                    FROM simulate_cascade_kinetics(%s, %s, %s, %s, %s, %s, 1000.0, 48.0, 96)
                    ORDER BY time_hr ASC;
                    """,
                    (p["kcat1"], p["km1"], p["kcat2"], p["km2"], p["ki"], p["boost"]),
                )
                with_sink = np.array(cur.fetchall(), dtype=float)

                # 2. Control Without Module 3 Sink (kcat2 = 0.0 -> intermediate accumulates & inhibits E1 via Ki)
                cur.execute(
                    """
                    SELECT time_hr, substrate_u_m, intermediate_u_m, sink_product_u_m, fractional_conversion
                    FROM simulate_cascade_kinetics(%s, %s, 0.0, %s, %s, %s, 1000.0, 48.0, 96)
                    ORDER BY time_hr ASC;
                    """,
                    (p["kcat1"], p["km1"], p["km2"], p["ki"], p["boost"]),
                )
                no_sink = np.array(cur.fetchall(), dtype=float)

                # 3. Control Without Module 1 Anchor & Without Sink (boost = 1.0, kcat2 = 0.0)
                cur.execute(
                    """
                    SELECT time_hr, substrate_u_m, intermediate_u_m, sink_product_u_m, fractional_conversion
                    FROM simulate_cascade_kinetics(%s, %s, 0.0, %s, %s, 1.0, 1000.0, 48.0, 96)
                    ORDER BY time_hr ASC;
                    """,
                    (p["kcat1"], p["km1"], p["km2"], p["ki"]),
                )
                wild_type = np.array(cur.fetchall(), dtype=float)

                results[p["id"]] = {
                    "meta": p,
                    "with_sink": with_sink,
                    "no_sink": no_sink,
                    "wild_type": wild_type,
                }

    elapsed_ms = (time.perf_counter() - t0) * 1000.0

    summary_table = Table(
        title=f"48-Hour RK4 Cascade Kinetics Summary (12 Trajectories × 97 Steps Solved in {elapsed_ms:.2f} ms in pg_bio)",
        box=box.SIMPLE_HEAVY,
    )
    summary_table.add_column("Pipeline ID", style="bold cyan")
    summary_table.add_column("Peak Intermediate (With Sink)", justify="right", style="yellow")
    summary_table.add_column("48h Substrate Cleaved (Wild-Type)", justify="right", style="dim")
    summary_table.add_column("48h Substrate Cleaved (No Sink, Ki Stalled)", justify="right", style="magenta")
    summary_table.add_column("48h Sink Product Yield (Full 4-Module)", justify="right", style="bold green")
    summary_table.add_column("Sink Acceleration", justify="right", style="bold cyan")

    for p in PIPELINES:
        r = results[p["id"]]
        ws = r["with_sink"]
        ns = r["no_sink"]
        wt = r["wild_type"]

        peak_int_ws = np.max(ws[:, 2])
        wt_cleaved_pct = (1000.0 - wt[-1, 1]) / 10.0
        ns_cleaved_pct = (1000.0 - ns[-1, 1]) / 10.0
        ws_sink_pct = ws[-1, 4] * 100.0
        gain_vs_wt = ws_sink_pct / max(wt_cleaved_pct, 1e-6)

        summary_table.add_row(
            p["id"],
            f"{peak_int_ws:6.1f} µM",
            f"{wt_cleaved_pct:5.1f}%",
            f"{ns_cleaved_pct:5.1f}%",
            f"{ws_sink_pct:5.1f}% ({ws[-1, 3]:.0f} µM)",
            f"{gain_vs_wt:.2f}× vs WT",
        )

    console.print(summary_table)

    # Render 2x2 Multi-Feedstock RK4 Cascade Figure
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()

    for idx, p in enumerate(PIPELINES):
        ax = axes[idx]
        r = results[p["id"]]
        t = r["with_sink"][:, 0]
        ws = r["with_sink"]
        ns = r["no_sink"]
        wt = r["wild_type"]

        ax.plot(t, ws[:, 1], label=f"Remaining {p['sub_label']} (4-Module)", color="#264653", linestyle=":", linewidth=1.8)
        ax.plot(t, ws[:, 2], label=f"Transient {p['int_label']} (Drained by Sink)", color="#f4a261", linestyle="--", linewidth=2.0)
        ax.plot(t, ws[:, 3], label=f"{p['prod_label']} (Module 3 Sink)", color=p["color"], linewidth=2.8)
        ax.plot(t, 1000.0 - ns[:, 1], label="Cleaved w/o Sink (Ki Product-Stalled)", color="#6c757d", linestyle="-.", linewidth=1.8)
        ax.plot(t, 1000.0 - wt[:, 1], label="Wild-Type Scissor Alone (No Hook/Sink)", color="#adb5bd", linestyle="-", linewidth=1.5)

        ax.fill_between(t, ws[:, 3], alpha=0.12, color=p["color"])
        ax.set_title(p["title"], fontsize=11, fontweight="bold", pad=8)
        ax.set_xlabel("Incubation Time (Hours)")
        ax.set_ylabel("Concentration (µM eq)")
        ax.set_ylim(-20, 1050)
        ax.grid(True, alpha=0.28)
        ax.legend(fontsize=7.8, loc="center right")

    plt.suptitle(
        "In-Database 4th-Order Runge-Kutta (RK4) Enzymatic Cascade Digital Twin (pg_bio Layer 14)",
        fontsize=13,
        fontweight="bold",
        y=0.99,
    )
    plt.tight_layout()
    out_png = os.path.join(OUTPUT_DIR, "rk4_four_feedstock_cascades.png")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    console.print(f"[bold green]✅ Saved 4-Feedstock RK4 Cascade chart to:[/] [cyan]{out_png}[/]")


if __name__ == "__main__":
    run_rk4_multi_pipeline_benchmark()
