-- =============================================================================
-- bio-pipelines/base/schema.sql
-- Universal 4-Module Synthetic Biology Compiler Schema for pg_bio
-- Supports:
--   1. ptfe-to-bioceramic        (PTFE/PFAS -> Fluorapatite Ca5(PO4)3F)
--   2. plastics-upcycling        (PET/PU -> Monomers / Vanillin / FDCA)
--   3. sugarcane-bagasse         (Cellulose/Xylan/Lignin -> C5/C6 & 2,5-FDCA)
--   4. corn-stover-bran          (Arabinoxylan Diferulate Bridges -> Bio-Vanillin)
-- =============================================================================

CREATE TABLE IF NOT EXISTS bio_substrates (
    pipeline_id VARCHAR(50) PRIMARY KEY,
    feedstock_name TEXT NOT NULL,
    target_bond VARCHAR(60) NOT NULL,
    bond_dissociation_kj_mol FLOAT8 NOT NULL,
    representative_smiles TEXT NOT NULL,
    reactor_temp_c FLOAT8 NOT NULL,
    reactor_ph FLOAT8 NOT NULL,
    equilibrium_sink_product TEXT NOT NULL,
    optical_beacon_readout TEXT NOT NULL
);

TRUNCATE bio_substrates CASCADE;

INSERT INTO bio_substrates (
    pipeline_id, feedstock_name, target_bond, bond_dissociation_kj_mol,
    representative_smiles, reactor_temp_c, reactor_ph,
    equilibrium_sink_product, optical_beacon_readout
) VALUES
(
    'ptfe-to-bioceramic',
    'PTFE / PFAS Fluoropolymers',
    'C-F (Aliphatic Fluorocarbon)',
    540.0,
    'FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F',
    37.0,
    7.4,
    'Fluorapatite Ca5(PO4)3F (via AMELX Tag)',
    'crcB Fluoride Riboswitch -> Ruby-Red Phytoene Desaturase (W9DPJ6)'
),
(
    'plastics-upcycling',
    'PET & Ester-Urethane Plastics',
    'Aromatic Ester / Urethane (C(=O)-O)',
    360.0,
    'OC(=O)c1ccc(C(=O)OCCO)cc1',
    68.0,
    8.0,
    'Ca-Terephthalate / Protocatechuate -> Vanillin or 2,5-FDCA',
    'PcaU Terephthalate Biosensor -> Violacein / Carotenoid Beacon'
),
(
    'sugarcane-bagasse',
    'Steam-Exploded Sugarcane Bagasse (42% Cellulose, 25% Xylan, 20% Lignin)',
    'beta-1,4-Glycosidic + Lignin beta-O-4 Ether',
    310.0,
    'OCC1OC(OC2C(O)C(O)C(OC2CO)O)C(O)C(O)C1O',
    55.0,
    5.0,
    'Cellobiose Phosphorylase + HMF/FDCA Oxidase -> 2,5-FDCA',
    'XylR / FerC Xylose-Ferulate Biosensor -> Lycopene Readout'
),
(
    'corn-stover-bran',
    'Corn Stover & Pericarp Bran (Arabinoxylan-Ferulate Crosslinks)',
    'Diferulate-Arabinoxylan Ester Bridge',
    335.0,
    'COc1cc(/C=C/C(=O)OCC2OC(O)C(O)C2O)ccc1O',
    50.0,
    5.5,
    'FDC1 + VpVAN Cascade -> Bio-Vanillin',
    'FerC Ferulate-Inducible Chromoprotein Beacon'
);

-- Universal 4-Module Domain Atlas populated across the 2.34M protein space
CREATE TABLE IF NOT EXISTS bio_domain_atlas (
    id SERIAL PRIMARY KEY,
    pipeline_family VARCHAR(40) NOT NULL,
    module_role VARCHAR(40) NOT NULL,     -- 'M1_Anchor', 'M2_Scissor', 'M3_Sink', 'M4_Beacon'
    subtype VARCHAR(50) NOT NULL,         -- e.g. 'Hydrophobin', 'CBM', 'PETase_Cutinase', 'Cellulase_LPMO', 'Feruloyl_Esterase'
    uniprot_id VARCHAR(20) NOT NULL,
    organism_name TEXT NOT NULL,
    protein_description TEXT,
    seq_length INT NOT NULL,
    mw_kda FLOAT8 NOT NULL,
    isoelectric_pi FLOAT8 NOT NULL,
    gravy_hydropathy FLOAT8 NOT NULL,
    shannon_entropy FLOAT8 NOT NULL,
    UNIQUE (pipeline_family, module_role, subtype, uniprot_id)
);

CREATE INDEX IF NOT EXISTS idx_bio_domain_family_role
    ON bio_domain_atlas (pipeline_family, module_role, subtype);
CREATE INDEX IF NOT EXISTS idx_bio_domain_gravy
    ON bio_domain_atlas (gravy_hydropathy DESC);

-- Chimeric Swarm Optimization & Pareto Front Runs
CREATE TABLE IF NOT EXISTS bio_chimeric_runs (
    run_id SERIAL PRIMARY KEY,
    pipeline_id VARCHAR(50) NOT NULL,
    candidate_rank INT NOT NULL,
    anchor_uniprot VARCHAR(20) NOT NULL,
    linker_sequence TEXT NOT NULL,
    scissor_uniprot VARCHAR(20) NOT NULL,
    sink_tag VARCHAR(40) NOT NULL,
    host_code VARCHAR(30) NOT NULL,
    anchor_gravy FLOAT8,
    scissor_mw_kda FLOAT8,
    composite_pareto_score FLOAT8,
    predicted_48h_conversion FLOAT8,
    created_at TIMESTAMPTZ DEFAULT NOW()
);
