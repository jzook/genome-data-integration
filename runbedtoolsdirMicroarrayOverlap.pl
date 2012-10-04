	#!/usr/bin/perl
    use strict;
    use warnings;

    my $directory = "/Volumes/SSD960/workspace/plots";
    my $directoryspace = "/Volumes/SSD960/workspace/plots/";
    my $status;
#java -jar -Xmx2g /Applications/bioinfo/GenomeAnalysisTKLite-2.0-35-gb3e7fbe/GenomeAnalysisTKLite.jar -R /Volumes/UHTS2/IonTorrentVallone/reference/PVallone_120406ENO1_ref.fasta -T UnifiedGenotyper -I /Volumes/SSD960/ValloneSOLiD/5500_12_07_10_FC1_2_03_30/5500_12_07_10_FC1_2_03-2-Idx_30-30.bam -stand_call_conf 20.0 -stand_emit_conf 2.0 -A AlleleBalance -A RMSMappingQuality -A QualByDepth -A FisherStrand -A DepthOfCoverage -o /Volumes/SSD960/ValloneSOLiD/5500_12_07_10_FC1_2_03_30/5500_12_07_10_FC1_2_03-2-Idx_30-30_GATK_conf202.vcf -nt 8
    opendir (DIR, $directory) or die $!;
#exit;
my $i=1;
    while (my $file = readdir(DIR)) {
    my $filestart="";
#		if ($file =~ /(.*vcf6.*(f|r|t))\.vcf$/) { $filestart = $directoryspace . $1; 
		if ($file =~ /(.*vcf11.*(f|r|t))\.vcf$/) { $filestart = $directoryspace . $1; 
#		if ($file =~ /(.*NoHomMapfilt_cons.*(f|r|t))\.vcf$/) { $filestart = $directoryspace . $1; 
#		if ($file =~ /(.*NoHomMapfilt_cons.*(f|r|t)_callable)\.vcf$/) { $filestart = $directoryspace . $1; 
		} else {next;}
		print $filestart . "\n";
		
		my $command;
#	        $command = "intersectBed -v -a ${filestart}.vcf -b /Volumes/SSD960/workspace/data/vcfannot/bed/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.20120117.defaultcallable_noGLMTXY.bed > ${filestart}_vcf6notcallable.vcf";
	        $command = "intersectBed -v -a ${filestart}.vcf -b /Volumes/SSD960/workspace/data/vcfannot/bed/union11callableonlymerged_nouncert_addcert.bed > ${filestart}_notcallable.vcf";
        print "$command\n";
        $status = system($command);
        if ($status != 0) {
            print "intersectBed failed!\n";
            exit(-1);
        }

#	        $command = "wc -l ${filestart}_vcf6notcallable.vcf >> ${directoryspace}vcf6notcallableconscallablecounts.txt";
	        $command = "wc -l ${filestart}_notcallable.vcf >> ${directoryspace}consnotcallablecounts.txt";
        print "$command\n";
        $status = system($command);
        if ($status != 0) {
            print "wc failed!\n";
        #    exit(-1);
        }

#	        $command = "intersectBed -a ${filestart}.vcf -b /Volumes/SSD960/workspace/data/vcfannot/bed/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.20120117.defaultcallable_noGLMTXY.bed > ${filestart}_vcf6callable.vcf";
	        $command = "intersectBed -a ${filestart}.vcf -b /Volumes/SSD960/workspace/data/vcfannot/bed/union11callableonlymerged_nouncert_addcert.bed > ${filestart}_callable.vcf";
        print "$command\n";
        $status = system($command);
        if ($status != 0) {
            print "intersectBed failed!\n";
            exit(-1);
        }
#	        $command = "wc -l ${filestart}_vcf6callable.vcf >> ${directoryspace}vcf6callableconscallablecounts.txt";
	        $command = "wc -l ${filestart}_callable.vcf >> ${directoryspace}conscallablecounts.txt";
        print "$command\n";
        $status = system($command);
        if ($status != 0) {
            print "wc failed!\n";
       #     exit(-1);
        }

$i++;
    }

closedir(DIR);
#/Users/jzook/tophat-1.4.1.OSX_x86_64/tophat -o /Volumes/Blank\ Disk/SEQC/MAY/SEQC_ILM_MAY_A_1_L01_ATCACG_AD0E8UACXX_tophat_T -T -G /Volumes/Blank\ Disk/TestRuns/human.ucscGenesERCCWithoutPolyA.gtf -p 8 /Volumes/Blank\ Disk/bowtieindex/hg19_ERCCWithoutPolyA_mito_rRNA /Volumes/Blank\ Disk/SEQC/MAY/SEQC_ILM_MAY_A_1_L01_ATCACG_AD0E8UACXX_R1.fastq /Volumes/Blank\ Disk/SEQC/MAY/SEQC_ILM_MAY_A_1_L01_ATCACG_AD0E8UACXX_R2.fastq
