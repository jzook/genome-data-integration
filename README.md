genome-data-integration
=======================

Scripts used to integrate data from multiple genome sequencing datasets and form consensus variant calls.  The process and these scripts used to generate consensus calls for NA12878 are described in our manuscript "Integrating human sequence data sets provides a resource of benchmark SNP and indel genotype calls" at http://www.nature.com/nbt/journal/vaop/ncurrent/nbt.2835 and http://arxiv.org/abs/1307.4661.  The calls generated for NA12878 described in http://genomeinabottle.org/blog-entry/new-version-manuscript-and-highly-confident-genotypes-na12878 and available at ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/variant_calls/NIST. 

These scripts were used on Sun Grid Engine to parallelize processing.  They are provided in order to help with understanding the manuscript, but are not currently written to be easily adapted for use by others.  Rather, the resulting genotype calls for NA12878 are intended to be a resource for performance assessment.

Any questions can be directed to Justin Zook at NIST (jzook@nist.gov)