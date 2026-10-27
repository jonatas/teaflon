import os
import sys
import time
import psycopg
from rich.console import Console
from rich.table import Table

# Add the local pgbio-py SDK to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../pg_bio/pgbio-py")))
from pgbio import PgBioClient

console = Console()
DB_URL = os.getenv("PG_BIO_URL", "postgresql://localhost:28818/bio_demo")

WETLAB_HOSTS = [
    ("e_coli", "T7", "Escherichia coli BL21(DE3) (Mesophilic Bacterial Standard)"),
    ("b_subtilis", "T7", "Bacillus subtilis 168 (Gram+ Secretory Industrial Host)"),
    ("p_putida", "T7", "Pseudomonas putida KT2440 (Xenobiotic/Halocarbon Bioremediation Host)"),
    ("c_necator", "T7", "Cupriavidus necator H16 (Chemolithoautotrophic Bioplastic Host)"),
    ("s_cerevisiae", "GAL1", "Saccharomyces cerevisiae BY4741 (Eukaryotic Disulfide Folding Host)"),
    ("p_pastoris", "AOX1", "Komagataella phaffii (Pichia pastoris) (High-Yield Hydrophobin Secretion)"),
    ("h_salinarum", "PFDA", "Halobacterium salinarum NRC-1 (Extreme Halophile / 4M NaCl Host)"),
    ("s_acidocaldarius", "PFDA", "Sulfolobus acidocaldarius DSM 639 (Hyperthermophile 75°C / pH 2.5 Host)"),
    ("h_sapiens", "KRT14", "Homo sapiens Keratinocyte SafeLocus_B (Chr17 Bioceramic Nail Matrix)"),
]


def generate_orca_and_wetlab_cassettes(fasta_path: str, output_dir: str):
    console.print(f"[bold cyan]🧬 Reading TeaFlon TriFusion Protein FASTA:[/] {fasta_path}")
    with open(fasta_path, "r") as f:
        lines = f.readlines()
        protein_seq = "".join([l.strip() for l in lines if not l.startswith(">")])

    console.print(f"Loaded TriFusion Protein ([bold green]{len(protein_seq)}[/] Amino Acids)")
    client = PgBioClient(DB_URL)

    os.makedirs(output_dir, exist_ok=True)
    t0 = time.perf_counter()

    with psycopg.connect(DB_URL) as conn:
        with conn.cursor() as cur:
            cur.execute("""
                CREATE TABLE IF NOT EXISTS teaflon_expression_cassettes (
                    host_code VARCHAR(30) PRIMARY KEY,
                    host_description TEXT,
                    promoter_name VARCHAR(30),
                    coding_dna TEXT,
                    full_plasmid_dna TEXT,
                    length_bp INT,
                    gc_content FLOAT8,
                    cai_score FLOAT8,
                    restriction_sites_free BOOLEAN,
                    roundtrip_verified BOOLEAN
                );
                TRUNCATE teaflon_expression_cassettes;
            """)

            table = Table(title="TeaFlon Multi-Host Wet-Lab Expression Cassettes (Native pg_bio Layer 12)")
            table.add_column("Host Code", style="cyan")
            table.add_column("Expression Organism & Niche", style="white")
            table.add_column("Promoter", style="magenta")
            table.add_column("Length (bp)", style="green")
            table.add_column("GC %", style="yellow")
            table.add_column("CAI", style="bold green")
            table.add_column("BsaI Free", style="blue")
            table.add_column("Round-Trip", style="bold cyan")

            multi_fasta_path = os.path.join(output_dir, "teaflon_multihost_plasmids.fasta")
            orca_input_path = os.path.join(output_dir, "teaflon_orca_input.fasta")

            with open(multi_fasta_path, "w") as mf:
                for host_code, promoter, desc in WETLAB_HOSTS:
                    cassette = client.design_expression_cassette(
                        protein_seq=protein_seq,
                        host_organism=host_code,
                        promoter_type=promoter,
                    )
                    cur.execute("SELECT bio_translate_dna(%s);", (cassette.coding_dna,))
                    translated_aa = cur.fetchone()[0]
                    roundtrip_ok = (translated_aa == protein_seq)

                    cur.execute(
                        """
                        INSERT INTO teaflon_expression_cassettes
                        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s);
                        """,
                        (
                            host_code,
                            desc,
                            cassette.promoter_name,
                            cassette.coding_dna,
                            cassette.full_plasmid_dna,
                            cassette.length_bp,
                            cassette.gc_content,
                            cassette.cai_score,
                            cassette.restriction_sites_free,
                            roundtrip_ok,
                        ),
                    )

                    mf.write(
                        f">teaflon_trifusion_{host_code}|promoter={cassette.promoter_name}|bp={cassette.length_bp}|gc={cassette.gc_content:.4f}|cai={cassette.cai_score:.3f}\n"
                    )
                    mf.write(cassette.full_plasmid_dna + "\n")

                    if host_code == "e_coli":
                        with open(orca_input_path, "w") as of:
                            of.write(">chr_synthetic_teaflon_plasmid\n")
                            of.write(cassette.full_plasmid_dna + "\n")

                    table.add_row(
                        host_code,
                        desc,
                        cassette.promoter_name,
                        f"{cassette.length_bp:,}",
                        f"{cassette.gc_content * 100:.1f}%",
                        f"{cassette.cai_score:.3f}",
                        "✅ YES" if cassette.restriction_sites_free else "❌ NO",
                        "✅ EXACT" if roundtrip_ok else "❌ MISMATCH",
                    )

            conn.commit()

    elapsed_ms = (time.perf_counter() - t0) * 1000.0
    console.print(table)
    console.print(f"[bold green]⚡ Designed and verified 9 wet-lab expression plasmids natively in pg_bio in {elapsed_ms:.2f} ms![/]")
    console.print(f"Saved E. coli ORCA 3D chromatin input to: [cyan]{orca_input_path}[/]")
    console.print(f"Saved all 9 multi-host wet-lab plasmid FASTAs to: [cyan]{multi_fasta_path}[/]")


if __name__ == "__main__":
    fasta_file = "teaflon_trifusion.fasta"
    if not os.path.exists(fasta_file):
        print(f"Error: {fasta_file} not found.")
        sys.exit(1)

    generate_orca_and_wetlab_cassettes(fasta_file, "orca_validation")

