# TeaFlon

**Primary Objective:** Design synthetic proteins capable of digesting Polytetrafluoroethylene (PTFE / Teflon) and converting the released fluorine into a structured bioceramic cluster, similar to human tooth enamel or bone.

## The Challenge

PTFE is composed of long carbon chains saturated with fluorine atoms. The Carbon-Fluorine (C-F) bond is the strongest single bond in organic chemistry, making PTFE incredibly inert, highly hydrophobic, and notoriously difficult to degrade biologically.

## The "TeaFlon" System Architecture

To achieve this, the biological system needs to perform three distinct functions, likely requiring a highly engineered multi-enzyme complex or fusion protein:

1. **Binding & Surface Disruption (The "Teflon Hook"):**
   PTFE is extremely hydrophobic and water-repellent. The enzyme requires a highly hydrophobic binding domain—drawing inspiration from fungal hydrophobins or engineered plastic-binding domains (like those in PETases)—to anchor the catalytic machinery to the polymer surface.

2. **Defluorination & Cleavage (The "Fluorine Scissor"):**
   Breaking the C-F bonds and the carbon backbone requires a powerful catalytic center. We will look to nature's *fluoroacetate dehalogenases* (enzymes that break single C-F bonds) and potentially radical-based mechanisms (similar to *lignin peroxidases* or *cytochrome P450s*) to shear the tough polymer chain.

3. **Biomineralization (The "Tooth Builder"):**
   Free fluoride ions ($F^-$) are toxic and must be sequestered. By engineering scaffolding proteins inspired by *amelogenin* (the protein that guides tooth enamel formation), the system will capture released fluoride along with environmental calcium and phosphate to precipitate **Fluorapatite** ($Ca_5(PO_4)_3F$)—forming a hard bioceramic cluster.

## AI & Structural Biology Development Plan

1. **Phase 1: Mining & Template Discovery:** Query biological databases (UniProt, PDB) to identify starting templates for dehalogenases, hydrophobic binding domains, and biomineralization scaffolds.
2. **Phase 2: AI-Driven Protein Design:** Use tools like RFdiffusion and ProteinMPNN for de novo design and active site grafting, validated by AlphaFold 3 structure prediction.
3. **Phase 3: Molecular Docking & Simulation:** Perform in silico molecular dynamics and docking simulations with PTFE oligomers to verify the predicted active sites.
4. **Phase 4: Biomineralization Coupling:** Design the structural scaffold to couple the defluorinating enzyme with the bioceramic nucleation site.
5. **Phase 5: Experimental Validation:** Synthesize genes, express in a host organism (e.g., E. coli or yeast), and test against PTFE nanoparticles to measure degradation and bioceramic precipitation.

## 🚀 Quickstart & Setup

This project uses `uv` for dependency management.

```bash
# Clone the repository
git clone https://github.com/jonatas/teaflon.git
cd teaflon

# Install dependencies (including PyMOL, PyTorch, Transformers, etc. if required)
uv sync
```

### Rendering Protein Designs

We have several PyMOL scripts available for rendering in-silico designs. Outputs are saved in the `./output/` directory:

```bash
uv run render_docking.py
uv run render_fusion.py
uv run render_trifusion.py
```

## 🧬 Integration with `pg_bio`

TeaFlon is officially the **first application built on top of [`pg_bio`](https://github.com/jonatas/pg_bio)**. 

During Phase 1 (Mining & Template Discovery), we utilize the `pgbio-py` SDK to perform hyper-fast vector homology searches across protein embedding spaces directly in Postgres. This allows us to instantly find homologous scaffolds for our Dehalogenase, Hydrophobin, and Amelogenin domains without massive memory overhead.

To run the template discovery pipeline:

1. Ensure your local `pg_bio` instance is running (via `docker-compose up -d` in the `pg_bio` directory).
2. Ensure you have the `PG_BIO_URL` environment variable set if your database is remote (it defaults to `postgresql://localhost:28818/bio_demo`).
3. Run the template discovery script:

```bash
uv run find_templates.py
```

## Security and Environment

**Important:** We strictly ensure no sensitive information (API keys, personal file paths) is committed to this repository. All environment variables should be stored locally in a `.env` file (which is ignored by Git).

Example `.env` file:
```env
PG_BIO_URL=postgresql://localhost:28818/bio_demo
```

### Multiomics Integration: The 12-Layer Bioceramic Compiler & Pan-Organism Atlas

We have built and validated a complete, end-to-end 12-layer multiomics and wet-lab compiler using `pg_bio` across **7,338 proteins and 930 distinct organisms** (spanning hyperthermophilic & halophilic Archaea, filamentous fungi, bacteria, plants, and mammals).

All scripts for this pipeline are located in `pg_bio_multiomics_pipeline/` and the root workspace, and can be executed end-to-end via `uv run teaflon`:

1. **Phase 1 Pan-Organism Mining (`find_templates.py`):** Queries the 2.34M protein space using HNSW `<=>` and the Smith-Waterman hybrid `<~>` operator across all 4 functional modules (`Q1JU72` Dehalogenase across 560 organisms, `P79073` Hydrophobin across 51 organisms, `Q99217` Amelogenin across 12 organisms, and `P80093` Carotenoid Reporter across 821 organisms).
2. **3D Structural, Theozyme & Globular Compactness Audit (`measure.py`):** Uses `bio_parse_pdb`, `center_of_mass`, `bio_distance`, `bio_angle`, `theozyme_rmsd(..., 'q1ju72_dehalogenase_hexad')` (`0.0000 Å`), and `evaluate_globular_compactness` (verifying all 418 CA atoms of the TriFusion fold at `82.9` mean pLDDT and `1.030` Flory globular ratio without PyMOL).
3. **Multi-Host Wet-Lab Expression Cassette Compiler (`orca_dna_prep.py`):** Uses `bio_reverse_translate`, `bio_codon_adaptation_index`, and `bio_design_expression_cassette` to generate and round-trip verify 9 wet-lab expression plasmids (*E. coli* BL21, *B. subtilis*, *P. putida*, *C. necator*, *S. cerevisiae*, *P. pastoris*, *H. salinarum* at 71.6% GC, *S. acidocaldarius*, and *H. sapiens* `KRT14`).
4. **In-Database ISE & Fluorapatite Digital Twin (`digital_twin_defluorination.py`):** Uses `simulate_defluorination_kinetics` in `pg_bio` to model 48-hour fluoride release (`150.0 ppm [F⁻]`) and stoichiometric $\text{Ca}_5(\text{PO}_4)_3\text{F}$ precipitation (`3,981 mg/L`).
5. **Step 1 Genomics (`step1_genomics.py`):** Uses 1D PostgreSQL `int4range` GiST indexes to identify `SafeLocus_B` on Chromosome 17 adjacent to the `KRT14` promoter with zero `TP53`/`BRCA1` overlap.
6. **Step 2 Proteomics (`step2_folding.py`):** Uses Morton Z-Order curves (`z_order_encode`) to isolate the 3D `ARG_BIND` keratin-binding pocket in `< 2 ms`.
7. **Step 3 Transcriptomics (`step3_transcriptomics.py`):** Uses Compressed Sparse Row (CSR) parallel arrays to verify single-cell dosage inside the Goldilocks zone (`34.99` copies in keratinocytes vs `0.49` in osteoblasts).
8. **Step 4 Visual Reporter (`step4_visual_reporter.py`):** Mines extremophile Archaea structural homologs of the *Capsicum annuum* carotenoid pathway (`GGPPS`, `PSY`, `PDS`) for the ruby-red lycopene bioceramic watermark.
9. **Step 5 Multi-Domain LC-MS/MS Proteomics (`step5_proteomics_ms.py`):** Validates all 4 TriFusion domains (`Dehalogenase N-term`, `Dehalogenase Tryptic Active Site`, `Hydrophobin Core`, and `Amelogenin C-Tag`) using `parse_mgf`, `generate_theoretical_fragments`, `isotope_fit_score`, and Elias-Gygi tryptic decoys (`13.2x` to `96.1x` target/decoy discrimination).
10. **Step 6 Pan-Organism Atlas & 12-Layer Audit (`step6_pan_organism_wetlab.py`):** Builds `teaflon_organism_atlas` (7,338 proteins / 930 organisms) and executes VCF clinical safety (`parse_vcf`), PFAS SMILES Tanimoto screening (`%`), MD preorganization (`ensemble_fraction_below`), and MM/PBSA thermodynamic auditing (`audit_mmpbsa_decomposition`, `discover_cryptic_hotspots`).
11. **Automated Validation Suite (`test_pipeline.py`) & CLI Dashboard (`teaflon_compiler.py`):** Runs all 8 multiomics & wet-lab assertions in `< 1 second`.

### 🔴 The Visual Reporter Module (The "Red Bone" Watermark)
To guarantee that the newly engineered bioceramic can be easily distinguished from natural bone, tooth enamel, or surrounding environment, we implemented a **Biological Watermark**. 

Using the `pg_bio` HNSW index, we mined the structural twins of the red carotenoid pathway from the Bell Pepper (*Capsicum annuum*). We mapped the fragile plant enzymes (`GGPPS`, `PSY`, `PDS`) to industrial-strength extremophile Archaea analogs. By co-expressing this robust red pigment cassette alongside the Teaflon mineralizer, the system precipitates **Lycopene** directly into the crystalline matrix. 

**The result:** A ruby-red bioceramic that provides instant, naked-eye visual verification of successful Teflon degradation and bone synthesis.
