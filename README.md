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
