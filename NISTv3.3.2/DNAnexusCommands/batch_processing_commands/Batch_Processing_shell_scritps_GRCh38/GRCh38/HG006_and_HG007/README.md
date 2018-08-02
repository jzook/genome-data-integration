## Scripts for processing GIAB Chinese Parents (HG006 and HG007)

- `giab_illumina_variant_calling_pipeline.sh` downloads bam files and splits by chromosome then runs callableLoci, Freebayes, and Sention variant caller on individual chromosomes.
- `run_illumina_chinese_paraents.sh` calls the Illumina pipeline script for the Mother and Father 100X and MatePair Illumina datasets.
- `get_depth.sh` extract coverage information from vcf files
- `coverage_stats_from_vcf.Rmd` calculate median coverage and max coverage threshold for use generating bed files with callable regions (callableLoci).
- `giab_10X_variant_calling_piepline.sh` dowloads and split bams by chromosome and haploid, calls variants on each haploid.
- `run_10X_chinese_parents.sh` calls the 10X pipeline for Mother and Father.
- `giab_complete_genomics_pipeline.sh` prepares Complete Genomics VCF for integration pipeline and performs liftover from GRCh37 to GRCh38. Liftover performed locally.
- `run_complete_genomics_chinese_parents.sh` run cg pipeline on Mother and Father CG vcfs.
- Callset table generation
  - `Make_callSet_tables.Rmd` - generate callset tables for HG006 and HG007 for both GRCh37 and GRCh38
    - `variant_callsets.tsv` list of variant callsets used in pipeline and annotation files
    - `filt_tbl_Asian_Trio_GRCh3{7,8}.tsv` filter table for HG006/HG007 GRCh37 and GRCh38 callsets.
- `run_integration_pipe_HG00{6,7}_GRCh3{7,8}.sh` script to run integration pipeline by chromosome.
