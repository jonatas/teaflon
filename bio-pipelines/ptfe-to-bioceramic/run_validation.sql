-- =============================================================================
-- bio-pipelines/ptfe-to-bioceramic/run_validation.sql
-- Reproduces TeaFlon PTFE/PFAS -> Fluorapatite Ca5(PO4)3F Benchmarks in pg_bio
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
WHERE pipeline_id = 'ptfe-to-bioceramic'
ORDER BY candidate_rank ASC;
