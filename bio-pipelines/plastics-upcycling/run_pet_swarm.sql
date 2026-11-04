-- =============================================================================
-- bio-pipelines/plastics-upcycling/run_pet_swarm.sql
-- Swarm & Pareto Search for 68-72°C Tg-Tolerant PETase/Cutinase-CBM Chimeras
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
WHERE pipeline_id = 'plastics-upcycling'
ORDER BY candidate_rank ASC;
