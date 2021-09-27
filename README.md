genome-data-integration
=======================

This repository contains scripts used to integrate data from multiple genome sequencing datasets and form high-confidence SNV, indel, and homozygous reference calls for Genome in a Bottle described in the preprint https://doi.org/10.1101/2020.07.24.212712.

For the NISTv4.2.1, we used the results for v3.3.2 plus new callsets for PacBio HiFi and 10x Genomics, in addition to better excluding structural variant and copy number variant regions. Most analyses were performed using apps or applets on DNAnexus The apps and applets are structured as:
dxapp.json specifies the input files and options, output files, and any dependencies that can be installed via apt.
src/code.sh contains the commands that are run
resources/ contains compiled binary files, scripts, and other files that are used in the applet

For the deprecated NISTv3.3.2, most analyses were performed using apps or applets on DNAnexus, except for mapping of all datasets and variant calling for Complete Genomics and Ion exome, since these steps were performed by others.  The apps and applets used in this work are included as directories under NISTv3.3.2.   They use an Ubuntu 12.04 machine on Amazon Web Services EC2.  The apps and applets are structured as:
dxapp.json specifies the input files and options, output files, and any dependencies that can be installed via apt.
src/code.sh contains the commands that are run
resources/ contains compiled binary files, scripts, and other files that are used in the applet

The commands were run per chromosome in parallel using the DNAnexus command line interface.  Note that some applets contain software that requires licenses for some or all uses, in particular GATK and Sentieon.

For deprecated NISTv2.19, the process and these scripts used to generate consensus calls for NA12878 were described in our manuscript "Integrating human sequence data sets provides a resource of benchmark SNP and indel genotype calls" at http://www.nature.com/nbt/journal/vaop/ncurrent/nbt.2835.  These scripts were used on Sun Grid Engine to parallelize processing.  They are provided in order to help with understanding the manuscript, but are not currently written to be easily adapted for use by others.  Rather, the resulting genotype calls for NA12878 are intended to be a resource for performance assessment.

Any questions can be submitted as issues to this repository or emailed to Justin Zook at NIST.