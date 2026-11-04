-- =============================================================================
-- bio-pipelines/base/swarm_pareto.sql
-- Populates bio_domain_atlas across all 4 Feedstock Families from the 2.34M
-- protein table in bio_demo and exposes Pareto-optimal candidate views.
-- =============================================================================

TRUNCATE bio_domain_atlas;

INSERT INTO bio_domain_atlas (
    pipeline_family, module_role, subtype, uniprot_id, organism_name,
    protein_description, seq_length, mw_kda, isoelectric_pi,
    gravy_hydropathy, shannon_entropy
)
WITH candidate_domains AS (
    -- 1. Universal & Pipeline-Specific Surface Anchors (Module 1)
    SELECT 'ptfe-to-bioceramic'::varchar(40) AS fam,
           'M1_Anchor'::varchar(40) AS role,
           'Hydrophobin'::varchar(50) AS sub,
           uniprot_id, name, sequence
    FROM proteins
    WHERE name ILIKE '%hydrophobin%' AND sequence IS NOT NULL AND length(sequence) BETWEEN 40 AND 400

    UNION ALL
    SELECT 'plastics-upcycling', 'M1_Anchor', 'CBM_Hydrophobin', uniprot_id, name, sequence
    FROM proteins
    WHERE (name ILIKE '%carbohydrate-binding%' OR name ILIKE '%cellulose-binding%' OR name ILIKE '%hydrophobin%')
      AND sequence IS NOT NULL AND length(sequence) BETWEEN 40 AND 500

    UNION ALL
    SELECT 'sugarcane-bagasse', 'M1_Anchor', 'CBM_Cellulose_Anchor', uniprot_id, name, sequence
    FROM proteins
    WHERE (name ILIKE '%cellulose-binding%' OR name ILIKE '%carbohydrate-binding%' OR name ILIKE '%dockerin%')
      AND sequence IS NOT NULL AND length(sequence) BETWEEN 35 AND 600

    UNION ALL
    SELECT 'corn-stover-bran', 'M1_Anchor', 'CBM_Arabinoxylan_Anchor', uniprot_id, name, sequence
    FROM proteins
    WHERE (name ILIKE '%carbohydrate-binding%' OR name ILIKE '%xylan-binding%' OR name ILIKE '%cellulose-binding%')
      AND sequence IS NOT NULL AND length(sequence) BETWEEN 35 AND 600

    -- 2. Primary Bond Scissors (Module 2)
    UNION ALL
    SELECT 'ptfe-to-bioceramic', 'M2_Scissor', 'Dehalogenase', uniprot_id, name, sequence
    FROM proteins
    WHERE name ILIKE '%dehalogenase%' AND sequence IS NOT NULL AND length(sequence) BETWEEN 80 AND 600

    UNION ALL
    SELECT 'plastics-upcycling', 'M2_Scissor', 'PETase_Cutinase_Esterase', uniprot_id, name, sequence
    FROM proteins
    WHERE (name ILIKE '%cutinase%' OR name ILIKE '%petase%' OR name ILIKE '%terephthal%' OR name ILIKE '%arylesterase%')
      AND sequence IS NOT NULL AND length(sequence) BETWEEN 120 AND 650

    UNION ALL
    SELECT 'sugarcane-bagasse', 'M2_Scissor', 'Lignocellulose_3Enzyme_Swarm', uniprot_id, name, sequence
    FROM proteins
    WHERE (name ILIKE '%cellulase%' OR name ILIKE '%endoglucanase%' OR name ILIKE '%xylanase%' OR name ILIKE '%laccase%' OR name ILIKE '%dye-decolorizing peroxidase%' OR name ILIKE '%lytic polysaccharide monooxygenase%')
      AND sequence IS NOT NULL AND length(sequence) BETWEEN 100 AND 850

    UNION ALL
    SELECT 'corn-stover-bran', 'M2_Scissor', 'Feruloyl_Arabinofuranosidase', uniprot_id, name, sequence
    FROM proteins
    WHERE (name ILIKE '%feruloyl esterase%' OR name ILIKE '%arabinofuranosidase%')
      AND sequence IS NOT NULL AND length(sequence) BETWEEN 120 AND 850

    -- 3. Equilibrium Sinks (Module 3)
    UNION ALL
    SELECT 'ptfe-to-bioceramic', 'M3_Sink', 'AMELX_Biomineral_Scaffold', uniprot_id, name, sequence
    FROM proteins
    WHERE (name ILIKE '%amelogen%' OR name ILIKE '%enamelin%' OR name ILIKE '%biomineral%')
      AND sequence IS NOT NULL AND length(sequence) >= 30

    UNION ALL
    SELECT 'plastics-upcycling', 'M3_Sink', 'Protocatechuate_Dioxygenase', uniprot_id, name, sequence
    FROM proteins
    WHERE (name ILIKE '%protocatechuate%' OR name ILIKE '%catechol 1,2-dioxygenase%')
      AND sequence IS NOT NULL AND length(sequence) BETWEEN 120 AND 600

    UNION ALL
    SELECT 'sugarcane-bagasse', 'M3_Sink', 'Cellobiose_Phosphorylase_Oxidase', uniprot_id, name, sequence
    FROM proteins
    WHERE (name ILIKE '%cellobiose phosphorylase%' OR name ILIKE '%beta-glucosidase%' OR name ILIKE '%alcohol oxidase%')
      AND sequence IS NOT NULL AND length(sequence) BETWEEN 150 AND 900

    UNION ALL
    SELECT 'corn-stover-bran', 'M3_Sink', 'Ferulate_Decarboxylase_Vanillin', uniprot_id, name, sequence
    FROM proteins
    WHERE (name ILIKE '%ferulic acid decarboxylase%' OR name ILIKE '%phenolic acid decarboxylase%' OR name ILIKE '%vanillin%')
      AND sequence IS NOT NULL AND length(sequence) BETWEEN 100 AND 650

    -- 4. Optical Beacons (Module 4)
    UNION ALL
    SELECT 'all-pipelines', 'M4_Beacon', 'Carotenoid_Chromoprotein_Beacon', uniprot_id, name, sequence
    FROM proteins
    WHERE (name ILIKE '%phytoene desaturase%' OR name ILIKE '%phytoene synthase%' OR name ILIKE '%lycopene cyclase%')
      AND sequence IS NOT NULL AND length(sequence) BETWEEN 150 AND 650
)
SELECT DISTINCT ON (fam, role, sub, uniprot_id)
    fam,
    role,
    sub,
    uniprot_id,
    CASE
        WHEN name LIKE '%OS=%' THEN split_part(split_part(name, 'OS=', 2), ' OX=', 1)
        ELSE split_part(split_part(name, '_', 2), ' ', 1)
    END AS org_name,
    split_part(name, ' OS=', 1) AS prot_desc,
    length(sequence) AS seq_len,
    molecular_weight(sequence) / 1000.0 AS mw_kda,
    isoelectric_point(sequence) AS pi,
    gravy_score(sequence) AS gravy,
    sequence_entropy(sequence) AS entropy
FROM candidate_domains
WHERE sequence IS NOT NULL AND length(sequence) > 0;
