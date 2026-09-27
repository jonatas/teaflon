import sys
import os

# Add the local pgbio-py SDK to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../pg_bio/pgbio-py")))

try:
    from pgbio import PgBioClient
except ImportError:
    print("Error: Could not import PgBioClient. Ensure pg_bio SDK is available.")
    sys.exit(1)

def main():
    # Connect to the local pg_bio database
    # In a real setup, this would be an environment variable.
    # For now, we will use a safe default local url.
    db_url = os.getenv("PG_BIO_URL", "postgresql://localhost:28818/bio_demo")
    print(f"Connecting to pg_bio at {db_url}...")
    
    try:
        client = PgBioClient(db_url)
    except Exception as e:
        print(f"Failed to connect to pg_bio: {e}")
        print("Make sure pg_bio is running (docker-compose up -d in pg_bio directory).")
        return

    # Phase 1: Mining & Template Discovery
    print("\n--- Phase 1: Mining & Template Discovery using pg_bio ---")
    
    # We need:
    # 1. Dehalogenases (Fluorine Scissor)
    # 2. Hydrophobin or PETase-like domains (Teflon Hook)
    # 3. Amelogenin or Biomineralization scaffolds (Tooth Builder)
    
    # We will search for homologues using representative sequences.
    # For demonstration, we use some dummy or short sequences representing the domains
    # In practice, these would be the full FASTA sequences of known enzymes.
    
    targets = {
        "Dehalogenase (Fluorine Scissor)": "MFEGFERRLVD",  # Example short seq
        "Hydrophobin (Teflon Hook)": "VCPGLCCSQYG",      # Example short seq
        "Amelogenin (Tooth Builder)": "MPLPPHPGHPG"      # Example short seq
    }
    
    for name, seq in targets.items():
        print(f"\nSearching for {name} templates...")
        print(f"Query sequence: {seq}")
        try:
            homologues = client.find_homologues(query_sequence=seq, limit=3)
            if homologues:
                for i, match in enumerate(homologues):
                    print(f"  Match {i+1}: UniProt ID {match.uniprot_id} (Distance: {match.embedding_distance:.4f})")
            else:
                print("  No matches found in the current bio_demo database.")
        except Exception as e:
            print(f"  Error during search: {e}")
            print("  Note: pg_bio might need to be seeded with data first.")

if __name__ == "__main__":
    main()
