## Scripts for processing GIAB Chinese Parents (HG006 and HG007)

- `giab_illumina_variant_calling_pipeline.sh` downloads bam files and splits by chromosome then runs callableLoci, Freebayes, and Sention variant caller on individual chromosomes.
- `run_illumina_chinese_paraents.sh` calls the Illumina pipeline script for the Mother and Father 100X and MatePair Illumina datasets. 
