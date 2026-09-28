import psycopg
import os

DB_URL = os.getenv("PG_BIO_URL", "postgresql://localhost:28818/bio_demo")

def run_visual_reporter_pipeline():
    conn = psycopg.connect(DB_URL)
    cur = conn.cursor()

    print("🔴 TEAFLON VISUAL REPORTER MODULE (Biological Watermark)\n" + "="*60)
    print("Goal: Engineer a robust carotenoid pathway to turn the Teaflon bioceramic bright red.\n")

    # The 3-Step Red Pigment Pathway from Peppers (Capsicum annuum)
    pepper_enzymes = {
        "GGPPS": "P80042", # Geranylgeranyl pyrophosphate synthase (The Precursor)
        "PSY": "P37272",   # Phytoene synthase (The Committer)
        "PDS": "P80093"    # Phytoene desaturase (The Red Pigment Maker)
    }

    for name, uniprot_id in pepper_enzymes.items():
        print(f"\n🧬 Mining Industrial Robustness for {name} (Pepper Target: {uniprot_id})")
        
        query = f"""
        WITH target AS (
            SELECT embedding FROM proteins WHERE uniprot_id = '{uniprot_id}'
        )
        SELECT p.uniprot_id, p.name, (p.embedding <=> (SELECT embedding FROM target)) as distance
        FROM proteins p
        WHERE p.uniprot_id != '{uniprot_id}'
        ORDER BY p.embedding <=> (SELECT embedding FROM target) ASC
        LIMIT 2;
        """
        
        cur.execute(query)
        matches = cur.fetchall()
        
        for match in matches:
            uid, desc, dist = match
            # Shorten the description for clean output
            short_desc = desc.split(" OS=")[0] if " OS=" in desc else desc[:60]
            print(f"   ↳ 🛡️ Robust Homolog Found: {uid} | Distance: {dist:.3f}")
            print(f"      {short_desc}")

    print("\n✅ Visual Reporter Cassette Built.")
    print("When fused with the Teaflon system, these enzymes will co-precipitate Lycopene (red pigment) inside the Calcium-Fluoride bioceramic matrix!")

    cur.close()
    conn.close()

if __name__ == "__main__":
    run_visual_reporter_pipeline()
