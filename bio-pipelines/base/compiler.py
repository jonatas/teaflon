import glob
import os
import time
import psycopg
import yaml
from rich.console import Console
from rich.table import Table
from rich import box

console = Console()
DB_URL = os.getenv("PG_BIO_URL", "postgresql://localhost:28818/bio_demo")
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
PIPELINES_ROOT = os.path.dirname(BASE_DIR)


def run_declarative_compiler(rebuild_atlas: bool = False):
    console.print("[bold cyan]🧬 pg_bio Declarative Bioprocess Compiler (Layer 7 + Layer 14 Rust Primitives)[/]")
    console.print("=" * 94)

    t0 = time.perf_counter()

    with psycopg.connect(DB_URL) as conn:
        with conn.cursor() as cur:
            cur.execute("SELECT to_regclass('public.bio_domain_atlas');")
            has_table = cur.fetchone()[0] is not None
            atlas_count = 0
            if has_table:
                cur.execute("SELECT COUNT(*) FROM bio_domain_atlas;")
                atlas_count = cur.fetchone()[0]

            if rebuild_atlas or atlas_count == 0:
                schema_sql = open(os.path.join(BASE_DIR, "schema.sql")).read()
                pareto_sql = open(os.path.join(BASE_DIR, "swarm_pareto.sql")).read()
                cur.execute(schema_sql)
                cur.execute(pareto_sql)
                conn.commit()

            # 1. Summarize Universal 4-Module Atlas across all 4 Feedstock Pipelines
            cur.execute(
                """
                SELECT pipeline_family,
                       module_role,
                       subtype,
                       COUNT(*) AS proteins,
                       COUNT(DISTINCT organism_name) AS organisms,
                       ROUND(AVG(mw_kda)::numeric, 1) AS avg_mw,
                       ROUND(AVG(gravy_hydropathy)::numeric, 3) AS avg_gravy
                FROM bio_domain_atlas
                GROUP BY pipeline_family, module_role, subtype
                ORDER BY pipeline_family, module_role;
                """
            )
            rows = cur.fetchall()

            cur.execute("SELECT COUNT(*), COUNT(DISTINCT uniprot_id), COUNT(DISTINCT organism_name) FROM bio_domain_atlas;")
            total_entries, unique_prots, unique_orgs = cur.fetchone()

            table = Table(
                title=f"Universal 4-Module Domain Atlas ({total_entries:,} Entries | {unique_prots:,} Unique Proteins | {unique_orgs:,} Organisms)",
                box=box.SIMPLE_HEAVY,
            )
            table.add_column("Pipeline Family", style="cyan")
            table.add_column("4-Module Role", style="yellow")
            table.add_column("Biocatalytic Subtype", style="white")
            table.add_column("Proteins", justify="right", style="bold green")
            table.add_column("Organisms", justify="right", style="magenta")
            table.add_column("Avg MW (kDa)", justify="right", style="blue")
            table.add_column("Avg GRAVY", justify="right", style="cyan")

            for r in rows:
                table.add_row(r[0], r[1], r[2], f"{r[3]:,}", f"{r[4]:,}", f"{r[5]}", f"{r[6]}")
            console.print(table)

            # 2. Parse declarative YAML specs & run native Rust Layer 14 functions:
            #    - pareto_score_domain_pair()
            #    - simulate_cascade_kinetics() (RK4 2-stage integrator with Ki inhibition & adsorption boost)
            #    - optimize_consortium_stoichiometry() (3-component simplex optimizer)
            yaml_files = sorted(glob.glob(os.path.join(PIPELINES_ROOT, "*/*.yaml")))
            cur.execute("TRUNCATE bio_chimeric_runs;")

            spec_table = Table(
                title="In-Database Pareto Domain Ranking & 48h RK4 Cascade Simulation (`pareto_score_domain_pair` + `simulate_cascade_kinetics`)",
                box=box.SIMPLE_HEAVY,
            )
            spec_table.add_column("Pipeline ID", style="bold cyan")
            spec_table.add_column("Reactor", style="white")
            spec_table.add_column("SMILES %", justify="right", style="green")
            spec_table.add_column("Top Pareto Pair (Anchor × Scissor)", style="yellow")
            spec_table.add_column("Pareto Score", justify="right", style="magenta")
            spec_table.add_column("48h RK4 Sink Yield", justify="right", style="bold green")

            for yf in yaml_files:
                spec = yaml.safe_load(open(yf))
                pid = spec["pipeline_id"]
                smiles = spec["feedstock"]["smiles"]
                ref_smiles = spec["feedstock"]["reference_analog_smiles"]
                temp_c = spec["reactor"]["temperature_c"]
                ph = spec["reactor"]["ph"]
                linker = spec["linkers"][0]
                host = spec["hosts"][0]
                sink_prod = spec["modules"]["module_3_sink"]["sink_product"]

                # Layer 7 + Layer 14 Pareto + Layer 14 RK4 Cascade in a single SQL query
                cur.execute(
                    """
                    WITH top_anchors AS (
                        SELECT uniprot_id, gravy_hydropathy, shannon_entropy
                        FROM bio_domain_atlas
                        WHERE pipeline_family = %s AND module_role = 'M1_Anchor'
                        ORDER BY gravy_hydropathy DESC, shannon_entropy DESC
                        LIMIT 15
                    ),
                    top_scissors AS (
                        SELECT uniprot_id, mw_kda, isoelectric_pi, shannon_entropy
                        FROM bio_domain_atlas
                        WHERE pipeline_family = %s AND module_role = 'M2_Scissor'
                          AND mw_kda BETWEEN 18.0 AND 65.0
                        ORDER BY shannon_entropy DESC
                        LIMIT 25
                    ),
                    ranked AS (
                        SELECT a.uniprot_id AS anchor_id,
                               s.uniprot_id AS scissor_id,
                               a.gravy_hydropathy,
                               s.mw_kda,
                               pareto_score_domain_pair(
                                   a.gravy_hydropathy,
                                   a.shannon_entropy,
                                   s.mw_kda,
                                   s.isoelectric_pi,
                                   s.shannon_entropy,
                                   %s
                               ) AS pareto_score
                        FROM top_anchors a
                        CROSS JOIN top_scissors s
                        ORDER BY pareto_score DESC
                        LIMIT 3
                    )
                    SELECT r.anchor_id,
                           r.scissor_id,
                           r.gravy_hydropathy,
                           r.mw_kda,
                           r.pareto_score,
                           smiles_to_fingerprint(%s) %% smiles_to_fingerprint(%s) AS tanimoto,
                           (
                               SELECT fractional_conversion
                               FROM simulate_cascade_kinetics(
                                   24.0 * r.pareto_score,
                                   75.0,
                                   36.0 * r.pareto_score,
                                   40.0,
                                   120.0,
                                   3.5 + GREATEST(r.gravy_hydropathy, 0.0) * 4.0,
                                   1000.0,
                                   48.0,
                                   48
                               )
                               ORDER BY time_hr DESC
                               LIMIT 1
                           ) AS rk4_48h_conversion
                    FROM ranked r;
                    """,
                    (pid, pid, ph, smiles, ref_smiles),
                )
                pareto_hits = cur.fetchall()

                for rank_idx, hit in enumerate(pareto_hits, start=1):
                    anchor_id, scissor_id, a_gravy, s_mw, p_score, tanimoto, rk4_conv = hit
                    cur.execute(
                        """
                        INSERT INTO bio_chimeric_runs (
                            pipeline_id, candidate_rank, anchor_uniprot, linker_sequence,
                            scissor_uniprot, sink_tag, host_code, anchor_gravy,
                            scissor_mw_kda, composite_pareto_score, predicted_48h_conversion
                        ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);
                        """,
                        (
                            pid,
                            rank_idx,
                            anchor_id,
                            linker,
                            scissor_id,
                            sink_prod[:40],
                            host,
                            float(a_gravy),
                            float(s_mw),
                            float(p_score),
                            float(rk4_conv),
                        ),
                    )

                top_hit = pareto_hits[0]
                spec_table.add_row(
                    pid,
                    f"{temp_c:.0f}°C / pH {ph}",
                    f"{top_hit[5]:.3f}",
                    f"{top_hit[0]} × {top_hit[1]}",
                    f"{top_hit[4]:.4f}",
                    f"{top_hit[6] * 100:.1f}%",
                )

            # 3. Run Consortium Stoichiometry Optimizer for Multi-Component Agro-Waste Feedstocks
            stoich_table = Table(
                title="Multi-Component Biomass Consortium Stoichiometry (`optimize_consortium_stoichiometry`)",
                box=box.SIMPLE_HEAVY,
            )
            stoich_table.add_column("Agro-Residue Feedstock", style="cyan")
            stoich_table.add_column("Feedstock Fractions (F1 : F2 : F3)", style="white")
            stoich_table.add_column("Optimal Enzyme Ratio (w1 : w2 : w3)", style="bold yellow")
            stoich_table.add_column("48h Conversion", justify="right", style="bold green")
            stoich_table.add_column("Balance Index", justify="right", style="magenta")

            for label, f1, f2, f3, k1, k2, k3 in [
                ("Sugarcane Bagasse (Cellulose:Xylan:Lignin)", 0.42, 0.25, 0.20, 14.0, 18.0, 9.5),
                ("Corn Stover & Bran (Cellulose:Arabinoxylan:Ferulate)", 0.36, 0.28, 0.18, 13.5, 16.5, 12.0),
            ]:
                cur.execute(
                    "SELECT * FROM optimize_consortium_stoichiometry(%s, %s, %s, %s, %s, %s, 48.0);",
                    (f1, f2, f3, k1, k2, k3),
                )
                w1, w2, w3, conv, bal = cur.fetchone()
                stoich_table.add_row(
                    label,
                    f"{int(f1*100)}% : {int(f2*100)}% : {int(f3*100)}%",
                    f"{w1:.2f} : {w2:.2f} : {w3:.2f}",
                    f"{conv * 100:.1f}%",
                    f"{bal:.4f}",
                )

            conn.commit()
            console.print(spec_table)
            console.print(stoich_table)

    elapsed_ms = (time.perf_counter() - t0) * 1000.0
    console.print(f"[bold green]⚡ Executed 4-Feedstock Declarative Compiler & Layer 14 Rust Engine in {elapsed_ms:.2f} ms![/]")


if __name__ == "__main__":
    run_declarative_compiler()
