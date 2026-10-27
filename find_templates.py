import sys
import os
import time
import psycopg
from rich.console import Console
from rich.table import Table

# Add the local pgbio-py SDK to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../pg_bio/pgbio-py")))

try:
    from pgbio import PgBioClient
except ImportError:
    print("Error: Could not import PgBioClient. Ensure pg_bio SDK is available.")
    sys.exit(1)

console = Console()


def main():
    db_url = os.getenv("PG_BIO_URL", "postgresql://localhost:28818/bio_demo")
    console.print(f"[bold cyan]Connecting to pg_bio at {db_url}...[/]")

    try:
        client = PgBioClient(db_url)
    except Exception as e:
        console.print(f"[bold red]Failed to connect to pg_bio: {e}[/]")
        return

    console.print("\n[bold yellow]--- Phase 1: Mining & Template Discovery using pg_bio (Full 2.34M Protein Space) ---[/]")

    # Canonical reference UniProt IDs for the 4 TeaFlon functional modules:
    # 1. Dehalogenase (Fluorine Scissor): Burkholderia sp. FA1 Fluoroacetate Dehalogenase (Q1JU72, 304 aa)
    # 2. Hydrophobin (Teflon Hook): Trichoderma reesei Class II Hydrophobin HFB2 (P79073, 86 aa)
    # 3. Amelogenin (Tooth Builder): Homo sapiens Amelogenin X isoform AMELX (Q99217, 191 aa)
    # 4. Visual Reporter (Red Bioceramic Watermark): Capsicum annuum Phytoene Desaturase PDS (P80093, 582 aa)
    targets = [
        ("Dehalogenase (Fluorine Scissor)", "Q1JU72", "%dehalogenase%"),
        ("Hydrophobin (Teflon Hook)", "P79073", "%hydrophobin%"),
        ("Amelogenin (Tooth Builder)", "Q99217", "%amelogen%"),
        ("Carotenoid Reporter (Red Watermark)", "P80093", "%phytoene%"),
    ]

    with psycopg.connect(db_url) as conn:
        with conn.cursor() as cur:
            for module_name, bait_id, family_ilike in targets:
                t0 = time.perf_counter()
                cur.execute(
                    "SELECT sequence, name, LENGTH(sequence) FROM proteins WHERE uniprot_id = %s;",
                    (bait_id,),
                )
                row = cur.fetchone()
                if not row:
                    console.print(f"[red]Bait {bait_id} not found in proteins table.[/]")
                    continue
                bait_seq, bait_name, bait_len = row

                # Count total homologs and distinct organisms across the 2.34M database for this module
                cur.execute(
                    """
                    SELECT count(*),
                           count(DISTINCT CASE
                               WHEN name LIKE '%%OS=%%' THEN split_part(split_part(name, 'OS=', 2), ' OX=', 1)
                               ELSE split_part(split_part(name, '_', 2), ' ', 1)
                           END)
                    FROM proteins
                    WHERE name ILIKE %s;
                    """,
                    (family_ilike,),
                )
                total_prots, total_orgs = cur.fetchone()

                # Query top structural + sequence hybrid homologs natively inside PostgreSQL using HNSW `<=>` + `<~>`
                cur.execute("SET LOCAL hnsw.ef_search = 100;")
                cur.execute(
                    """
                    WITH target AS (
                        SELECT embedding, sequence FROM proteins WHERE uniprot_id = %s
                    ),
                    candidates AS (
                        SELECT p.uniprot_id,
                               p.name,
                               p.sequence,
                               (p.embedding <=> (SELECT embedding FROM target)) AS cos_dist
                        FROM proteins p
                        WHERE p.uniprot_id != %s
                          AND p.embedding IS NOT NULL
                        ORDER BY p.embedding <=> (SELECT embedding FROM target) ASC
                        LIMIT 25
                    )
                    SELECT c.uniprot_id,
                           c.name,
                           c.cos_dist,
                           (ROW(p.embedding::vector, c.sequence)::bio_feature <~>
                            ROW((SELECT embedding::vector FROM target), (SELECT sequence FROM target))::bio_feature) AS hybrid_score,
                           round((molecular_weight(c.sequence) / 1000.0)::numeric, 2) AS mw_kda,
                           round(isoelectric_point(c.sequence)::numeric, 2) AS pi,
                           round(gravy_score(c.sequence)::numeric, 3) AS gravy,
                           round(sequence_entropy(c.sequence)::numeric, 2) AS entropy
                    FROM candidates c
                    JOIN proteins p ON p.uniprot_id = c.uniprot_id
                    ORDER BY hybrid_score ASC
                    LIMIT 5;
                    """,
                    (bait_id, bait_id),
                )
                matches = cur.fetchall()
                elapsed_ms = (time.perf_counter() - t0) * 1000.0

                table = Table(
                    title=f"{module_name} | Bait: {bait_id} ({bait_len} aa) | {total_prots:,} Annotated Proteins across {total_orgs:,} Organisms ({elapsed_ms:.1f} ms)"
                )
                table.add_column("Rank", style="dim")
                table.add_column("UniProt ID", style="cyan")
                table.add_column("Organism / Description", style="white")
                table.add_column("Cosine (<=>)", style="green")
                table.add_column("Hybrid (<~>)", style="bold green")
                table.add_column("MW (kDa)", style="magenta")
                table.add_column("pI", style="yellow")
                table.add_column("GRAVY", style="blue")

                for idx, m in enumerate(matches, 1):
                    uid, desc, c_dist, h_score, mw, pi, gravy, _ent = m
                    short_desc = desc[:62] + "..." if len(desc) > 62 else desc
                    table.add_row(
                        str(idx),
                        uid,
                        short_desc,
                        f"{float(c_dist):.4f}",
                        f"{float(h_score):.4f}",
                        f"{mw}",
                        f"{pi}",
                        f"{gravy}",
                    )
                console.print(table)

            # Also verify PgBioClient.find_homologues() SDK method directly on full Q1JU72 sequence
            q1_seq = client.find_homologues(query_sequence="MFEGFERRLVDVGDVTINCVVGGSGPALLLLHGFPQNLHMWARVAPLLANEYTVVCADLRGYGGSSKPVGAPDHANYSFRAMASDQRELMRTLGFERFHLVGHDRGGRTGHRMALDHPDSVLSLAVLDIIPTYVMFEEVDRFVARAYWHWYFLQQPAPYPEKVIGADPDTFYEGCLFGWGATGADGFDPEQLEEYRKQWRDPAAIHGSCCDYRAGGTIDFELDHGDLGRQVQCPALVFSGSAGLMHSLFEMQVVWAPRLANMRFASLPGGHFFVDRFPDDTARILREFLSDARSGIHQTERRES", limit=3)
            console.print(
                f"\n[bold green]✅ SDK Verification:[/] `client.find_homologues(Q1JU72_full_seq)` top hit = [cyan]{q1_seq[0].uniprot_id}[/] (distance = {q1_seq[0].embedding_distance:.4f})"
            )


if __name__ == "__main__":
    main()

