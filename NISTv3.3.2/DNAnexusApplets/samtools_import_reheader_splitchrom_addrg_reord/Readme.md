# SAMtools Mappings Split by chromosome

## What does this app do?

This app creates individual bam and index files for each chromosome if named as 1, 2, ... X, Y.

## What are typical use cases for this app?

This splits a bam file by chromosome so that subsequent processing steps can be faster, particularly for large bam files so that only one chromosome needs to be loaded

## What data are required for this app to run?

This app requires coordinate-sorted mappings in BAM format (`*.bam`) and index in BAI format (`*.bai`).

## What does this app output?

This app outputs the generated bam and index files with _1.bam and _1.bai for chromosome 1, and so on.

## How does this app work?

This app runs 'samtools view' and 'samtools index' from the SAMtools suite of tools. For more information, consult the SAMtools website at:

http://samtools.sourceforge.net/
