import json
import time
import urllib.request
from rich.console import Console
from rich.table import Table
from rich import box

console = Console()
BASE_URL = "http://127.0.0.1:8090"
PIPELINES = [
    "ptfe-to-bioceramic",
    "plastics-upcycling",
    "sugarcane-bagasse",
    "corn-stover-bran",
]


def min_interchain_distance(pdb_text: str) -> float:
    chain_a = []
    chain_b = []
    for line in pdb_text.splitlines():
        if (line.startswith("ATOM  ") or line.startswith("HETATM")) and len(line) >= 54:
            ch = line[21:22]
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                if ch == "A":
                    chain_a.append((x, y, z))
                elif ch == "B":
                    chain_b.append((x, y, z))
            except ValueError:
                pass
    if not chain_a or not chain_b:
        return 0.0
    max_ax = max(p[0] for p in chain_a)
    min_bx = min(p[0] for p in chain_b)
    return min_bx - max_ax


def run_validation_suite():
    console.print("[bold cyan]🔬 Validating pg_bio Layer 14 `swarm_web.py` Across All 4 Feedstock Pipelines[/]")
    console.print("=" * 98)

    t0 = time.perf_counter()
    table = Table(
        title="End-to-End Mathematical, Kinetic & Structural Validation of `http://127.0.0.1:8090`",
        box=box.SIMPLE_HEAVY,
    )
    table.add_column("Pipeline ID", style="bold cyan")
    table.add_column("Matrix Contrast (min..max)", justify="right", style="yellow")
    table.add_column("Top Pair (`Anchor × Scissor`)", style="white")
    table.add_column("Table vs KPI RK4 Yield", justify="right", style="bold green")
    table.add_column("RK4 Mass Conservation", justify="right", style="magenta")
    table.add_column("Simplex Σw", justify="right", style="blue")
    table.add_column("3D Clearance ΔX", justify="right", style="bold cyan")

    for pid in PIPELINES:
        # 1. Fetch /api/swarm
        swarm = json.loads(
            urllib.request.urlopen(f"{BASE_URL}/api/swarm?pipeline={pid}").read().decode()
        )
        assert len(swarm["anchors"]) == 12, f"Expected 12 anchors for {pid}"
        assert len(swarm["scissors"]) == 16, f"Expected 16 scissors for {pid}"
        assert len(swarm["matrix_cells"]) == 192, f"Expected 192 cells for {pid}"
        assert swarm["max_score"] > swarm["min_score"] + 0.08, "Expected dynamic Pareto contrast across buckets"

        top = swarm["top_pairs"][0]
        p = swarm["preset"]
        eff_boost = top["effective_boost"]

        # 2. Fetch /api/rk4 for the #1 pair and verify exact match with leaderboard yield + mass conservation
        qs = f"pipeline={pid}&kcat1={p['kcat1']}&km1={p['km1']}&kcat2={p['kcat2']}&ki={p['ki']}&boost={eff_boost}"
        rk4 = json.loads(urllib.request.urlopen(f"{BASE_URL}/api/rk4?{qs}").read().decode())

        assert abs(top["yield_pct"] - rk4["sink_48h_pct"]) < 1e-6, (
            f"Mismatch between leaderboard yield ({top['yield_pct']}) and RK4 KPI ({rk4['sink_48h_pct']})"
        )
        assert top["yield_um"] == rk4["sink_48h_um"], (
            f"Mismatch between leaderboard uM ({top['yield_um']}) and RK4 KPI ({rk4['sink_48h_um']})"
        )
        assert rk4["sink_48h_pct"] > rk4["nosink_48h_pct"] + 20.0, (
            "Full 4-module sink yield must significantly outperform Ki-stalled no-sink control"
        )

        # Verify S(t) + I(t) + P(t) = 1000.0 uM across all 49 time points (within rounding tolerance)
        max_mass_err = 0.0
        for s_v, i_v, p_v in zip(rk4["with_sink_sub"], rk4["with_sink_int"], rk4["with_sink_prod"]):
            err = abs((s_v + i_v + p_v) - 1000.0)
            if err > max_mass_err:
                max_mass_err = err
        assert max_mass_err <= 0.3, f"Mass conservation error {max_mass_err} uM exceeded 0.3 uM rounding tolerance"

        # 3. Verify Consortium Stoichiometry Simplex Sum (w1 + w2 + w3 == 1.00)
        st = swarm["stoichiometry"]
        w_sum = float(st["w1"]) + float(st["w2"]) + float(st["w3"])
        assert abs(w_sum - 1.0) < 1e-6, f"Simplex weights did not sum to 1.0: {w_sum}"

        # 4. Fetch /api/complex3d for the #1 pair and verify >= 28.0 A clearance & globular metrics
        c3d = json.loads(
            urllib.request.urlopen(
                f"{BASE_URL}/api/complex3d?anchor={top['anchor_id']}&scissor={top['scissor_id']}"
            ).read().decode()
        )
        dx_clearance = min_interchain_distance(c3d["merged_pdb"])
        assert abs(dx_clearance - 28.0) < 0.05, f"Expected 28.0 A X-clearance, got {dx_clearance}"
        assert c3d["globular"]["ca_count"] > 100, "Expected > 100 Ca atoms in merged complex"

        table.add_row(
            pid,
            f"{swarm['min_score']:.3f} .. {swarm['max_score']:.3f}",
            f"{top['anchor_id']} × {top['scissor_id']}",
            f"{top['yield_pct']:.1f}% == {rk4['sink_48h_pct']:.1f}% ✅",
            f"1000.0 ± {max_mass_err:.1f} µM ✅",
            f"{w_sum:.2f} ✅",
            f"{dx_clearance:.2f} Å ✅",
        )

    elapsed_ms = (time.perf_counter() - t0) * 1000.0
    console.print(table)
    console.print(f"[bold green]✅ All 4 Pipelines Passed 100% of Mathematical, Kinetic & 3D Validation Checks in {elapsed_ms:.1f} ms![/]")


if __name__ == "__main__":
    run_validation_suite()
