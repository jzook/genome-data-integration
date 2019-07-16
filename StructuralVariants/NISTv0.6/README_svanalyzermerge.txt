December 13, 2017

This directory contains an analysis of the similarities between structural variant 
calls in the "union" vcf file union_171212_refalt.sort.vcf.gz. Comparisons between
variants were made by the program "SVmerge.pl" version 0.2, resulting in the file of
distances union_171212_refalt.complete.svdist.gz, which is then given as input to
the clustering program SVcluster.pl. SVmerge.pl and SVcluster.pl are available in
the "SVanalyzer" package at https://github.com/nhansen/SVanalyzer.

SVmerge.pl aligns the widened haplotypes for the alternate alleles created by two
structural variants, and then aligns them to each other with a global aligner (in
this case, edlib's Needleman-Wunsch implementation), reporting edit distance
between the two sequences. SVmerge.pl reports three different measures of distance
between the variants:

RELSHIFT - This is the largest shift in coordinates between the two sequences within
the global alignment, divided by the size of the larger allele (REF or ALT), averaged
between the two variants, to give a "relative shift" within the alignment
between the two resulting haplotypes.

RELSIZEDIFF - This is the difference in size of the two structural variants (i.e.,
the difference between the number of bases deleted or inserted between the two),
divided by the size of the larger allele (REF or ALT), averaged between the two 
variants, giving a relative size difference between the two variants.

RELDIST - This is the edit distance between the two resultant haplotypes, divided by
size of the larger allele (REF or ALT), averaged between the two variants, giving a 
"relative edit distance" between the two resulting haplotypes.

Once these three distances were calculated by SVmerge.pl for all pairs of variants 
within a given distance of each other (in this case, 100,000 bases), clusters of 
variants were formed using the program SVcluster.pl version 0.2 by creating a graph
with SVs as the nodes, connected by edges if all three distance measures between
them were 0 (for "strict" clusters), or less than 1 (for "wide" clusters). Structural
variant clusters were then created by performing a depth-first-search of this graph
to determine connected components.

The cluster files named "union_171212_refalt.withsinglets.clusters.##.##.##.txt.gz
were created by merging calls that had distances less than 0.## for the three
distances described above, respectively. They also include lines for "singlet"
calls, i.e., calls which do not cluster with any other call according to the 
criteria used to create the file.

