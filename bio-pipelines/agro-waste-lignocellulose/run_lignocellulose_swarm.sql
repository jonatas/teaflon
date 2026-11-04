-- =============================================================================
-- bio-pipelines/agro-waste-lignocellulose/run_lignocellulose_swarm.sql
-- Multi-Enzyme Consortium Stoichiometry & Chimera Search for Sugarcane Bagasse
-- (42% Cellulose, 25% Xylan, 20% Lignin) and Corn Stover Diferulate Bridges
-- =============================================================================

SELECT pipeline_id,
       candidate_rank,
       anchor_uniprot,
       linker_sequence,
       scissor_uniprot,
       sink_tag,
       host_code,
       anchor_gravy,
       scissor_mw_kda,
       composite_pareto_score,
       predicted_48h_conversion
FROM bio_chimeric_runs
WHERE pipeline_id IN ('sugarcane-bagasse', 'corn-stover-bran')
ORDER BY pipeline_id, candidate_rank ASC;
