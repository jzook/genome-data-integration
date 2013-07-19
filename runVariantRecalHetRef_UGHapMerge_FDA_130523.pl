#!/usr/bin/env perl
#run VariantRecalibrator for AB, SSEs, Alignment, and mapping biases
use POSIX;

#read in user parameters
my $filestart = shift @ARGV;    #filestart
my $dataset  = shift @ARGV; #reference fasta file
my $chrom = shift @ARGV;    #chromosome(s)
my $runNo = shift @ARGV;    #run number (combine allcall if it's 1)

my $bed;
if ($chrom eq "all") {
	$bed="";
} else { $bed="-L /data/results/justin/bedfiles/chr${chrom}.bed"; }

my $res = 1;

# open chrom list
if ($runNo==1 && $chrom eq "all") { #only run this part the first time
	unless ( open(CHROMS, "/home/justin.zook/GATK/runs/ChromNames.txt")) {
    	print "\nCannot open the file: ~/GATK/runs/ChromNames.txt! \n\n";
    	exit;
	}
	my $i=1; my $cats="";
	while (my $chromno=<CHROMS>) {
		if ($chromno =~ /(.*)\n/) {$chromno=$1;}
	     if ($i==1 && $res == 0) {
			my $samCmd = "gunzip -c ${filestart}_allcall_UGHapMerge_HetRef_${chromno}.vcf.gz > ${filestart}_Het_temp${i}.vcf";
			  print "$samCmd\n";
			 $res = system($samCmd);
			 $cats="cat ${filestart}_Het_temp${i}.vcf"
 	   	} elsif ($i>1 && $res == 0) {
			my $samCmd = "zgrep -v '^\\#' ${filestart}_allcall_UGHapMerge_HetRef_${chromno}.vcf.gz > ${filestart}_Het_temp${i}.vcf";
			  print "$samCmd\n";
			 $res = system($samCmd);
			 
			 $cats="$cats ${filestart}_Het_temp${i}.vcf"
 	   	}   
 	   	$i++;
	}
	close CHROMS ;
	     if ($res == 0) {
			 my $samCmd = "$cats | perl ~/scripts/vcfsorter.pl /scratch/justin.zook/references/human_g1k_v37.dict - > ${filestart}_allcall_UGHapMerge_HetRef_all.vcf";
			  print "$samCmd\n";
			 $res = system($samCmd);
			  $samCmd = "rm ${filestart}_Het_temp*.vcf";
			  print "$samCmd\n";
			 $res = system($samCmd);
 	   	} 

#$res=1;
     if ($res == 0) {
	 my $samCmd = " java -jar -Xmx4g  ~/GATK/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar -R /scratch/justin.zook/references/human_g1k_v37.fasta -T SelectVariants -V:variant ${filestart}_allcall_UGHapMerge_HetRef_${chrom}.vcf -selectType SNP -o ${filestart}_snp_allcall_UGHapMerge_HetRef_${chrom}.vcf -U LENIENT_VCF_PROCESSING";
	  print "$samCmd\n";
	 $res = system($samCmd);
    }   

     if ($res == 0) {
	 my $samCmd = " java -jar -Xmx4g  ~/GATK/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar -R /scratch/justin.zook/references/human_g1k_v37.fasta -T SelectVariants -V:variant ${filestart}_allcall_UGHapMerge_HetRef_${chrom}.vcf -selectType INDEL -o ${filestart}_indel_allcall_UGHapMerge_HetRef_${chrom}.vcf -U LENIENT_VCF_PROCESSING";
	  print "$samCmd\n";
	 $res = system($samCmd);
    }   
}

#$res=0;
	unless ( open(CHROMS, "/home/justin.zook/GATK/runs/ChromNames.txt")) {
    	print "\nCannot open the file: ~/GATK/runs/ChromNames.txt! \n\n";
    	exit;
	}
	my $i=1; my $cats="";
	while (my $chromno=<CHROMS>) {
		if ($chromno =~ /(.*)\n/) {$chromno=$1;}
	     if ($i==1 && $res == 0) {
			my $samCmd = "gunzip -c ${filestart}_${dataset}call_UGHapMerge_HetRefUncert_${chromno}.vcf.gz > ${filestart}_${dataset}call_Het_temp${i}.vcf";
			  print "$samCmd\n";
			 $res = system($samCmd);
			 $cats="cat ${filestart}_${dataset}call_Het_temp${i}.vcf"
 	   	} elsif ($i>1 && $res == 0) {
			my $samCmd = "zgrep -v '^\\#' ${filestart}_${dataset}call_UGHapMerge_HetRefUncert_${chromno}.vcf.gz > ${filestart}_${dataset}call_Het_temp${i}.vcf";
			  print "$samCmd\n";
			 $res = system($samCmd);
			 
			 $cats="$cats ${filestart}_${dataset}call_Het_temp${i}.vcf"
 	   	}   
 	   	$i++;
	}
	close CHROMS ;
	     if ($res == 0) {
			 my $samCmd = "$cats | perl ~/scripts/vcfsorter.pl /scratch/justin.zook/references/human_g1k_v37.dict - > ${filestart}_${dataset}call_UGHapMerge_HetRefUncert_all.vcf";
			  print "$samCmd\n";
			 $res = system($samCmd);
			  $samCmd = "rm ${filestart}_${dataset}call_Het_temp*.vcf";
			  print "$samCmd\n";
			 $res = system($samCmd);
 	   	} 

     if ($res == 0) {
	 my $samCmd = " java -jar -Xmx4g  ~/GATK/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar -R /scratch/justin.zook/references/human_g1k_v37.fasta -T SelectVariants -V:variant ${filestart}_${dataset}call_UGHapMerge_HetRefUncert_${chrom}.vcf -selectType SNP -o ${filestart}_snp_${dataset}call_UGHapMerge_HetRefUncert_${chrom}.vcf -U LENIENT_VCF_PROCESSING";
	  print "$samCmd\n";
	 $res = system($samCmd);
    }   
#$res=0;

#subtract SNPs from all records to get indels since otherwise some mixed sites are missed
     if ($res == 0) {
	 my $samCmd = " java -jar -Xmx4g  ~/GATK/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar -R /scratch/justin.zook/references/human_g1k_v37.fasta -T CombineVariants -V:all ${filestart}_${dataset}call_UGHapMerge_HetRefUncert_${chrom}.vcf -V:snp ${filestart}_snp_${dataset}call_UGHapMerge_HetRefUncert_${chrom}.vcf -o ${filestart}_${dataset}call_UGHapMerge_HetRefUncert_combsnp_${chrom}.vcf -U LENIENT_VCF_PROCESSING";
	  print "$samCmd\n";
	 $res = system($samCmd);
    }   
     if ($res == 0) {
	 my $samCmd = " java -jar -Xmx4g  ~/GATK/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar -R /scratch/justin.zook/references/human_g1k_v37.fasta -T SelectVariants -V:variant ${filestart}_${dataset}call_UGHapMerge_HetRefUncert_combsnp_${chrom}.vcf -select 'set == \"all\"' -o ${filestart}_indel_${dataset}call_UGHapMerge_HetRefUncert_${chrom}.vcf -U LENIENT_VCF_PROCESSING";
	  print "$samCmd\n";
	 $res = system($samCmd);
    }   


     if ($res == 0) {
	 my $samCmd = " export PATH=\$PATH:/usr/bin/";
	  print "$samCmd\n";
	 $res = system($samCmd);
    }   



$res=0;

#INDELS
 # SSEs
 if ($res == 0) {
 my $samCmd = "java -jar -Xmx4g  ~/GATK/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar -R /scratch/justin.zook/references/human_g1k_v37.fasta -T VariantRecalibrator -input ${filestart}_indel_${dataset}call_UGHapMerge_HetRefUncert_${chrom}.vcf -resource:consensusUGHap1,known=false,training=true,truth=true,prior=20.0 ${filestart}_indel_allcall_UGHapMerge_HetRef_${chrom}.vcf -resource:dbSNP129,known=true,training=false,truth=false,prior=5.0 /home/justin.zook/references/dbsnp_137.b37.excluding_sites_after_129.vcf -an FS -an BaseQRankSum -an NBQ -recalFile ${filestart}_indel_${dataset}call_UGHapMerge_Het_SSE.recal_1.6 -tranchesFile ${filestart}_indel_${dataset}call_UGHapMerge_Het_SSE.tranches_1.6 -rscriptFile ${filestart}_indel_${dataset}call_UGHapMerge_Het_SSE.plots_1.6.R -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 $bed -nt 1 --maxGaussians 1 -mode INDEL -U LENIENT_VCF_PROCESSING";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
#$res=0;
 # AB and QD
 if ($res == 0) {
 my $samCmd = "java -jar -Xmx4g  ~/GATK/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar -R /scratch/justin.zook/references/human_g1k_v37.fasta -T VariantRecalibrator -input ${filestart}_indel_${dataset}call_UGHapMerge_HetRefUncert_${chrom}.vcf -resource:consensusUGHap1,known=false,training=true,truth=true,prior=20.0 ${filestart}_indel_allcall_UGHapMerge_HetRef_${chrom}.vcf -resource:dbSNP129,known=true,training=false,truth=false,prior=5.0 /home/justin.zook/references/dbsnp_137.b37.excluding_sites_after_129.vcf -an QD -recalFile ${filestart}_indel_${dataset}call_UGHapMerge_Het_AB.recal_1.6 -tranchesFile ${filestart}_indel_${dataset}call_UGHapMerge_Het_AB.tranches_1.6 -rscriptFile ${filestart}_indel_${dataset}call_UGHapMerge_Het_AB.plots_1.6.R -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 $bed -nt 1 --percentBadVariants 0.08 --maxGaussians 1 -mode INDEL -U LENIENT_VCF_PROCESSING";
  print "$samCmd\n";
 $res = system($samCmd);
  
 if ($res != 0) {
 my $samCmd = "java -jar -Xmx4g  ~/GATK/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar -R /scratch/justin.zook/references/human_g1k_v37.fasta -T VariantRecalibrator -input ${filestart}_indel_${dataset}call_UGHapMerge_HetRefUncert_${chrom}.vcf -resource:consensusUGHap1,known=false,training=true,truth=true,prior=20.0 ${filestart}_indel_allcall_UGHapMerge_HetRef_${chrom}.vcf -resource:dbSNP129,known=true,training=false,truth=false,prior=5.0 /home/justin.zook/references/dbsnp_137.b37.excluding_sites_after_129.vcf -an QD -recalFile ${filestart}_indel_${dataset}call_UGHapMerge_Het_AB.recal_1.6 -tranchesFile ${filestart}_indel_${dataset}call_UGHapMerge_Het_AB.tranches_1.6 -rscriptFile ${filestart}_indel_${dataset}call_UGHapMerge_Het_AB.plots_1.6.R -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 $bed -nt 1 --percentBadVariants 0.3 --maxGaussians 1 -mode INDEL -U LENIENT_VCF_PROCESSING";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
 }

#$res=0;
 # Align
 if ($res == 0) {
 my $samCmd = "java -jar -Xmx4g  ~/GATK/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar -R /scratch/justin.zook/references/human_g1k_v37.fasta -T VariantRecalibrator -input ${filestart}_indel_${dataset}call_UGHapMerge_HetRefUncert_${chrom}.vcf -resource:consensusUGHap1,known=false,training=true,truth=true,prior=20.0 ${filestart}_indel_allcall_UGHapMerge_HetRef_${chrom}.vcf -resource:dbSNP129,known=true,training=false,truth=false,prior=5.0 /home/justin.zook/references/dbsnp_137.b37.excluding_sites_after_129.vcf -an HaplotypeScore -an ReadPosEndDist -an ReadPosRankSum -recalFile ${filestart}_indel_${dataset}call_UGHapMerge_Het_Align.recal_1.6 -tranchesFile ${filestart}_indel_${dataset}call_UGHapMerge_Het_Align.tranches_1.6 -rscriptFile ${filestart}_indel_${dataset}call_UGHapMerge_Het_Align.plots_1.6.R -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 $bed -nt 1 --maxGaussians 1 -mode INDEL -U LENIENT_VCF_PROCESSING";
  print "$samCmd\n";
 $res = system($samCmd);

 if ($res != 0) { #if it doesn't work, then increase percentBadVariants
 my $samCmd = "java -jar -Xmx4g  ~/GATK/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar -R /scratch/justin.zook/references/human_g1k_v37.fasta -T VariantRecalibrator -input ${filestart}_indel_${dataset}call_UGHapMerge_HetRefUncert_${chrom}.vcf -resource:consensusUGHap1,known=false,training=true,truth=true,prior=20.0 ${filestart}_indel_allcall_UGHapMerge_HetRef_${chrom}.vcf -resource:dbSNP129,known=true,training=false,truth=false,prior=5.0 /home/justin.zook/references/dbsnp_137.b37.excluding_sites_after_129.vcf -an HaplotypeScore -an ReadPosEndDist -an ReadPosRankSum -recalFile ${filestart}_indel_${dataset}call_UGHapMerge_Het_Align.recal_1.6 -tranchesFile ${filestart}_indel_${dataset}call_UGHapMerge_Het_Align.tranches_1.6 -rscriptFile ${filestart}_indel_${dataset}call_UGHapMerge_Het_Align.plots_1.6.R -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 $bed -nt 1  --percentBadVariants 0.08 --maxGaussians 1 -mode INDEL -U LENIENT_VCF_PROCESSING";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
}
#$res=1;
 # Mapping
 if ($res == 0) {
 my $samCmd = "java -jar -Xmx4g  ~/GATK/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar -R /scratch/justin.zook/references/human_g1k_v37.fasta -T VariantRecalibrator -input ${filestart}_indel_${dataset}call_UGHapMerge_HetRefUncert_${chrom}.vcf -resource:consensusUGHap1,known=false,training=true,truth=true,prior=20.0 ${filestart}_indel_allcall_UGHapMerge_HetRef_${chrom}.vcf -resource:dbSNP129,known=true,training=false,truth=false,prior=5.0 /home/justin.zook/references/dbsnp_137.b37.excluding_sites_after_129.vcf -an MQ -an MQRankSum -an MQ0Fraction -an DP -recalFile ${filestart}_indel_${dataset}call_UGHapMerge_Het_map.recal_1.6 -tranchesFile ${filestart}_indel_${dataset}call_UGHapMerge_Het_map.tranches_1.6 -rscriptFile ${filestart}_indel_${dataset}call_UGHapMerge_Het_map.plots_1.6.R -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 $bed -nt 1 --percentBadVariants 0.08 --maxGaussians 1 -mode INDEL -U LENIENT_VCF_PROCESSING";
  print "$samCmd\n";
 $res = system($samCmd);

 if ($res != 0) { #if it doesn't work, then remove MQ0Fraction
 my $samCmd = "java -jar -Xmx4g  ~/GATK/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar -R /scratch/justin.zook/references/human_g1k_v37.fasta -T VariantRecalibrator -input ${filestart}_indel_${dataset}call_UGHapMerge_HetRefUncert_${chrom}.vcf -resource:consensusUGHap1,known=false,training=true,truth=true,prior=20.0 ${filestart}_indel_allcall_UGHapMerge_HetRef_${chrom}.vcf -resource:dbSNP129,known=true,training=false,truth=false,prior=5.0 /home/justin.zook/references/dbsnp_137.b37.excluding_sites_after_129.vcf -an MQ -an MQRankSum -an DP -recalFile ${filestart}_indel_${dataset}call_UGHapMerge_Het_map.recal_1.6 -tranchesFile ${filestart}_indel_${dataset}call_UGHapMerge_Het_map.tranches_1.6 -rscriptFile ${filestart}_indel_${dataset}call_UGHapMerge_Het_map.plots_1.6.R -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 $bed -nt 1 --percentBadVariants 0.08 --maxGaussians 1 -mode INDEL -U LENIENT_VCF_PROCESSING";
  print "$samCmd\n";
 $res = system($samCmd);
 }   

 if ($res != 0) { #if it doesn't work, then remove MQ0Fraction
 my $samCmd = "java -jar -Xmx4g  ~/GATK/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar -R /scratch/justin.zook/references/human_g1k_v37.fasta -T VariantRecalibrator -input ${filestart}_indel_${dataset}call_UGHapMerge_HetRefUncert_${chrom}.vcf -resource:consensusUGHap1,known=false,training=true,truth=true,prior=20.0 ${filestart}_indel_allcall_UGHapMerge_HetRef_${chrom}.vcf -resource:dbSNP129,known=true,training=false,truth=false,prior=5.0 /home/justin.zook/references/dbsnp_137.b37.excluding_sites_after_129.vcf -an MQ -an DP -recalFile ${filestart}_indel_${dataset}call_UGHapMerge_Het_map.recal_1.6 -tranchesFile ${filestart}_indel_${dataset}call_UGHapMerge_Het_map.tranches_1.6 -rscriptFile ${filestart}_indel_${dataset}call_UGHapMerge_Het_map.plots_1.6.R -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 $bed -nt 1 --percentBadVariants 0.08 --maxGaussians 1 -mode INDEL -U LENIENT_VCF_PROCESSING";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
}

#$res=0;
 # SSE
 if ($res == 0) {
 my $samCmd = "java -jar -Xmx4g  ~/GATK/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar -R /scratch/justin.zook/references/human_g1k_v37.fasta -T ApplyRecalibration -input ${filestart}_indel_${dataset}call_UGHapMerge_HetRefUncert_${chrom}.vcf -recalFile ${filestart}_indel_${dataset}call_UGHapMerge_Het_SSE.recal_1.6 -tranchesFile ${filestart}_indel_${dataset}call_UGHapMerge_Het_SSE.tranches_1.6 -o ${filestart}_indel_${dataset}call_UGHapMerge_Het_SSE_recal90_${chrom}.vcf $bed -nt 1 --ts_filter_level 90.0 -mode INDEL -U LENIENT_VCF_PROCESSING";
  print "$samCmd\n";
 $res = system($samCmd);
 }   

 # AB and QD
 if ($res == 0) {
 my $samCmd = "java -jar -Xmx4g  ~/GATK/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar -R /scratch/justin.zook/references/human_g1k_v37.fasta -T ApplyRecalibration -input ${filestart}_indel_${dataset}call_UGHapMerge_HetRefUncert_${chrom}.vcf -recalFile ${filestart}_indel_${dataset}call_UGHapMerge_Het_AB.recal_1.6 -tranchesFile ${filestart}_indel_${dataset}call_UGHapMerge_Het_AB.tranches_1.6 -o ${filestart}_indel_${dataset}call_UGHapMerge_Het_AB_recal90_${chrom}.vcf $bed -nt 1 --ts_filter_level 90.0 -mode INDEL -U LENIENT_VCF_PROCESSING";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
#$res=0;
 # Align
 if ($res == 0) {
 my $samCmd = "java -jar -Xmx4g  ~/GATK/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar -R /scratch/justin.zook/references/human_g1k_v37.fasta -T ApplyRecalibration -input ${filestart}_indel_${dataset}call_UGHapMerge_HetRefUncert_${chrom}.vcf -recalFile ${filestart}_indel_${dataset}call_UGHapMerge_Het_Align.recal_1.6 -tranchesFile ${filestart}_indel_${dataset}call_UGHapMerge_Het_Align.tranches_1.6 -o ${filestart}_indel_${dataset}call_UGHapMerge_Het_Align_recal90_${chrom}.vcf $bed -nt 1 --ts_filter_level 90.0 -mode INDEL -U LENIENT_VCF_PROCESSING";
  print "$samCmd\n";
 $res = system($samCmd);
 }   

 # Mapping
 if ($res == 0) {
 my $samCmd = "java -jar -Xmx4g  ~/GATK/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar -R /scratch/justin.zook/references/human_g1k_v37.fasta -T ApplyRecalibration -input ${filestart}_indel_${dataset}call_UGHapMerge_HetRefUncert_${chrom}.vcf -recalFile ${filestart}_indel_${dataset}call_UGHapMerge_Het_map.recal_1.6 -tranchesFile ${filestart}_indel_${dataset}call_UGHapMerge_Het_map.tranches_1.6 -o ${filestart}_indel_${dataset}call_UGHapMerge_Het_map_recal90_${chrom}.vcf $bed -nt 1 --ts_filter_level 90.0 -mode INDEL -U LENIENT_VCF_PROCESSING";
  print "$samCmd\n";
 $res = system($samCmd);
 }   





#$res=1;


#SNPs
 # SSEs
 if ($res == 0) {
 my $samCmd = "java -jar -Xmx4g  ~/GATK/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar -R /scratch/justin.zook/references/human_g1k_v37.fasta -T VariantRecalibrator -input ${filestart}_snp_${dataset}call_UGHapMerge_HetRefUncert_${chrom}.vcf -resource:consensusUGHap1,known=false,training=true,truth=true,prior=20.0 ${filestart}_snp_allcall_UGHapMerge_HetRef_${chrom}.vcf -resource:dbSNP129,known=true,training=false,truth=false,prior=5.0 /home/justin.zook/references/dbsnp_137.b37.excluding_sites_after_129.vcf -an FS -an BaseQRankSum -an NBQ -recalFile ${filestart}_snp_${dataset}call_UGHapMerge_Het_SSE.recal_1.6 -tranchesFile ${filestart}_snp_${dataset}call_UGHapMerge_Het_SSE.tranches_1.6 -rscriptFile ${filestart}_snp_${dataset}call_UGHapMerge_Het_SSE.plots_1.6.R -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 $bed -nt 1 --maxGaussians 1 -U LENIENT_VCF_PROCESSING";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
#$res=0;
 # AB and QD
 if ($res == 0) {
 my $samCmd = "java -jar -Xmx4g  ~/GATK/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar -R /scratch/justin.zook/references/human_g1k_v37.fasta -T VariantRecalibrator -input ${filestart}_snp_${dataset}call_UGHapMerge_HetRefUncert_${chrom}.vcf -resource:consensusUGHap1,known=false,training=true,truth=true,prior=20.0 ${filestart}_snp_allcall_UGHapMerge_HetRef_${chrom}.vcf -resource:dbSNP129,known=true,training=false,truth=false,prior=5.0 /home/justin.zook/references/dbsnp_137.b37.excluding_sites_after_129.vcf -an QD -an ABHet -an OND -recalFile ${filestart}_snp_${dataset}call_UGHapMerge_Het_AB.recal_1.6 -tranchesFile ${filestart}_snp_${dataset}call_UGHapMerge_Het_AB.tranches_1.6 -rscriptFile ${filestart}_snp_${dataset}call_UGHapMerge_Het_AB.plots_1.6.R -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 $bed -nt 1 --percentBadVariants 0.08 --maxGaussians 1 -U LENIENT_VCF_PROCESSING";
  print "$samCmd\n";
 $res = system($samCmd);
  
 if ($res != 0) {
 my $samCmd = "java -jar -Xmx4g  ~/GATK/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar -R /scratch/justin.zook/references/human_g1k_v37.fasta -T VariantRecalibrator -input ${filestart}_snp_${dataset}call_UGHapMerge_HetRefUncert_${chrom}.vcf -resource:consensusUGHap1,known=false,training=true,truth=true,prior=20.0 ${filestart}_snp_allcall_UGHapMerge_HetRef_${chrom}.vcf -resource:dbSNP129,known=true,training=false,truth=false,prior=5.0 /home/justin.zook/references/dbsnp_137.b37.excluding_sites_after_129.vcf -an QD -an ABHet -an OND -recalFile ${filestart}_snp_${dataset}call_UGHapMerge_Het_AB.recal_1.6 -tranchesFile ${filestart}_snp_${dataset}call_UGHapMerge_Het_AB.tranches_1.6 -rscriptFile ${filestart}_snp_${dataset}call_UGHapMerge_Het_AB.plots_1.6.R -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 $bed -nt 1 --percentBadVariants 0.3 --maxGaussians 1 -U LENIENT_VCF_PROCESSING";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
 }

#$res=0;
 # Align
 if ($res == 0) {
 my $samCmd = "java -jar -Xmx4g  ~/GATK/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar -R /scratch/justin.zook/references/human_g1k_v37.fasta -T VariantRecalibrator -input ${filestart}_snp_${dataset}call_UGHapMerge_HetRefUncert_${chrom}.vcf -resource:consensusUGHap1,known=false,training=true,truth=true,prior=20.0 ${filestart}_snp_allcall_UGHapMerge_HetRef_${chrom}.vcf -resource:dbSNP129,known=true,training=false,truth=false,prior=5.0 /home/justin.zook/references/dbsnp_137.b37.excluding_sites_after_129.vcf -an HaplotypeScore -an ReadPosEndDist -an ReadPosRankSum -recalFile ${filestart}_snp_${dataset}call_UGHapMerge_Het_Align.recal_1.6 -tranchesFile ${filestart}_snp_${dataset}call_UGHapMerge_Het_Align.tranches_1.6 -rscriptFile ${filestart}_snp_${dataset}call_UGHapMerge_Het_Align.plots_1.6.R -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 $bed -nt 1 --maxGaussians 1 -U LENIENT_VCF_PROCESSING";
  print "$samCmd\n";
 $res = system($samCmd);

 if ($res != 0) { #if it doesn't work, then increase percentBadVariants
 my $samCmd = "java -jar -Xmx4g  ~/GATK/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar -R /scratch/justin.zook/references/human_g1k_v37.fasta -T VariantRecalibrator -input ${filestart}_snp_${dataset}call_UGHapMerge_HetRefUncert_${chrom}.vcf -resource:consensusUGHap1,known=false,training=true,truth=true,prior=20.0 ${filestart}_snp_allcall_UGHapMerge_HetRef_${chrom}.vcf -resource:dbSNP129,known=true,training=false,truth=false,prior=5.0 /home/justin.zook/references/dbsnp_137.b37.excluding_sites_after_129.vcf -an HaplotypeScore -an ReadPosEndDist -an ReadPosRankSum -recalFile ${filestart}_snp_${dataset}call_UGHapMerge_Het_Align.recal_1.6 -tranchesFile ${filestart}_snp_${dataset}call_UGHapMerge_Het_Align.tranches_1.6 -rscriptFile ${filestart}_snp_${dataset}call_UGHapMerge_Het_Align.plots_1.6.R -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 $bed -nt 1  --percentBadVariants 0.08 --maxGaussians 1 -U LENIENT_VCF_PROCESSING";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
}
#$res=1;
 # Mapping
 if ($res == 0) {
 my $samCmd = "java -jar -Xmx4g  ~/GATK/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar -R /scratch/justin.zook/references/human_g1k_v37.fasta -T VariantRecalibrator -input ${filestart}_snp_${dataset}call_UGHapMerge_HetRefUncert_${chrom}.vcf -resource:consensusUGHap1,known=false,training=true,truth=true,prior=20.0 ${filestart}_snp_allcall_UGHapMerge_HetRef_${chrom}.vcf -resource:dbSNP129,known=true,training=false,truth=false,prior=5.0 /home/justin.zook/references/dbsnp_137.b37.excluding_sites_after_129.vcf -an MQ -an MQRankSum -an MQ0Fraction -an DP -recalFile ${filestart}_snp_${dataset}call_UGHapMerge_Het_map.recal_1.6 -tranchesFile ${filestart}_snp_${dataset}call_UGHapMerge_Het_map.tranches_1.6 -rscriptFile ${filestart}_snp_${dataset}call_UGHapMerge_Het_map.plots_1.6.R -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 $bed -nt 1 --percentBadVariants 0.08 --maxGaussians 1 -U LENIENT_VCF_PROCESSING";
  print "$samCmd\n";
 $res = system($samCmd);

 if ($res != 0) { #if it doesn't work, then remove MQ0Fraction
 my $samCmd = "java -jar -Xmx4g  ~/GATK/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar -R /scratch/justin.zook/references/human_g1k_v37.fasta -T VariantRecalibrator -input ${filestart}_snp_${dataset}call_UGHapMerge_HetRefUncert_${chrom}.vcf -resource:consensusUGHap1,known=false,training=true,truth=true,prior=20.0 ${filestart}_snp_allcall_UGHapMerge_HetRef_${chrom}.vcf -resource:dbSNP129,known=true,training=false,truth=false,prior=5.0 /home/justin.zook/references/dbsnp_137.b37.excluding_sites_after_129.vcf -an MQ -an MQRankSum -an DP -recalFile ${filestart}_snp_${dataset}call_UGHapMerge_Het_map.recal_1.6 -tranchesFile ${filestart}_snp_${dataset}call_UGHapMerge_Het_map.tranches_1.6 -rscriptFile ${filestart}_snp_${dataset}call_UGHapMerge_Het_map.plots_1.6.R -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 $bed -nt 1 --percentBadVariants 0.08 --maxGaussians 1 -U LENIENT_VCF_PROCESSING";
  print "$samCmd\n";
 $res = system($samCmd);
 }   

 if ($res != 0) { #if it doesn't work, then remove MQ0Fraction
 my $samCmd = "java -jar -Xmx4g  ~/GATK/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar -R /scratch/justin.zook/references/human_g1k_v37.fasta -T VariantRecalibrator -input ${filestart}_snp_${dataset}call_UGHapMerge_HetRefUncert_${chrom}.vcf -resource:consensusUGHap1,known=false,training=true,truth=true,prior=20.0 ${filestart}_snp_allcall_UGHapMerge_HetRef_${chrom}.vcf -resource:dbSNP129,known=true,training=false,truth=false,prior=5.0 /home/justin.zook/references/dbsnp_137.b37.excluding_sites_after_129.vcf -an MQ -an DP -recalFile ${filestart}_snp_${dataset}call_UGHapMerge_Het_map.recal_1.6 -tranchesFile ${filestart}_snp_${dataset}call_UGHapMerge_Het_map.tranches_1.6 -rscriptFile ${filestart}_snp_${dataset}call_UGHapMerge_Het_map.plots_1.6.R -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 95.0 -tranche 90.0 $bed -nt 1 --percentBadVariants 0.08 --maxGaussians 1 -U LENIENT_VCF_PROCESSING";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
}

#$res=0;
 # SSE
 if ($res == 0) {
 my $samCmd = "java -jar -Xmx4g  ~/GATK/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar -R /scratch/justin.zook/references/human_g1k_v37.fasta -T ApplyRecalibration -input ${filestart}_snp_${dataset}call_UGHapMerge_HetRefUncert_${chrom}.vcf -recalFile ${filestart}_snp_${dataset}call_UGHapMerge_Het_SSE.recal_1.6 -tranchesFile ${filestart}_snp_${dataset}call_UGHapMerge_Het_SSE.tranches_1.6 -o ${filestart}_snp_${dataset}call_UGHapMerge_Het_SSE_recal90_${chrom}.vcf $bed -nt 1 --ts_filter_level 90.0 -U LENIENT_VCF_PROCESSING";
  print "$samCmd\n";
 $res = system($samCmd);
 }   

 # AB and QD
 if ($res == 0) {
 my $samCmd = "java -jar -Xmx4g  ~/GATK/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar -R /scratch/justin.zook/references/human_g1k_v37.fasta -T ApplyRecalibration -input ${filestart}_snp_${dataset}call_UGHapMerge_HetRefUncert_${chrom}.vcf -recalFile ${filestart}_snp_${dataset}call_UGHapMerge_Het_AB.recal_1.6 -tranchesFile ${filestart}_snp_${dataset}call_UGHapMerge_Het_AB.tranches_1.6 -o ${filestart}_snp_${dataset}call_UGHapMerge_Het_AB_recal90_${chrom}.vcf $bed -nt 1 --ts_filter_level 90.0 -U LENIENT_VCF_PROCESSING";
  print "$samCmd\n";
 $res = system($samCmd);
 }   
#$res=0;
 # Align
 if ($res == 0) {
 my $samCmd = "java -jar -Xmx4g  ~/GATK/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar -R /scratch/justin.zook/references/human_g1k_v37.fasta -T ApplyRecalibration -input ${filestart}_snp_${dataset}call_UGHapMerge_HetRefUncert_${chrom}.vcf -recalFile ${filestart}_snp_${dataset}call_UGHapMerge_Het_Align.recal_1.6 -tranchesFile ${filestart}_snp_${dataset}call_UGHapMerge_Het_Align.tranches_1.6 -o ${filestart}_snp_${dataset}call_UGHapMerge_Het_Align_recal90_${chrom}.vcf $bed -nt 1 --ts_filter_level 90.0 -U LENIENT_VCF_PROCESSING";
  print "$samCmd\n";
 $res = system($samCmd);
 }   

 # Mapping
 if ($res == 0) {
 my $samCmd = "java -jar -Xmx4g  ~/GATK/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar -R /scratch/justin.zook/references/human_g1k_v37.fasta -T ApplyRecalibration -input ${filestart}_snp_${dataset}call_UGHapMerge_HetRefUncert_${chrom}.vcf -recalFile ${filestart}_snp_${dataset}call_UGHapMerge_Het_map.recal_1.6 -tranchesFile ${filestart}_snp_${dataset}call_UGHapMerge_Het_map.tranches_1.6 -o ${filestart}_snp_${dataset}call_UGHapMerge_Het_map_recal90_${chrom}.vcf $bed -nt 1 --ts_filter_level 90.0 -U LENIENT_VCF_PROCESSING";
  print "$samCmd\n";
 $res = system($samCmd);
 }   

#$res=0;
 # Merge SNP and indel calls
 if ($res == 0) {
	my $samCmd = "grep -v '^\\#' ${filestart}_indel_${dataset}call_UGHapMerge_Het_SSE_recal90_${chrom}.vcf > ${filestart}_indel_${dataset}call_UGHapMerge_Het_SSE_recal90_${chrom}_nohead.vcf";
	  print "$samCmd\n";
	 $res = system($samCmd);
	 $samCmd = "cat ${filestart}_snp_${dataset}call_UGHapMerge_Het_SSE_recal90_${chrom}.vcf ${filestart}_indel_${dataset}call_UGHapMerge_Het_SSE_recal90_${chrom}_nohead.vcf | perl ~/scripts/vcfsorter.pl /scratch/justin.zook/references/human_g1k_v37.dict - > ${filestart}_${dataset}call_UGHapMerge_Het_SSE_recal90_${chrom}.vcf";
	 print "$samCmd\n";
	$res = system($samCmd);
	$samCmd = "rm ${filestart}_indel_${dataset}call_UGHapMerge_Het_SSE_recal90_${chrom}_nohead.vcf";
	print "$samCmd\n";
	 $res = system($samCmd);
} 

 if ($res == 0) {
	my $samCmd = "grep -v '^\\#' ${filestart}_indel_${dataset}call_UGHapMerge_Het_AB_recal90_${chrom}.vcf > ${filestart}_indel_${dataset}call_UGHapMerge_Het_AB_recal90_${chrom}_nohead.vcf";
	  print "$samCmd\n";
	 $res = system($samCmd);
	 $samCmd = "cat ${filestart}_snp_${dataset}call_UGHapMerge_Het_AB_recal90_${chrom}.vcf ${filestart}_indel_${dataset}call_UGHapMerge_Het_AB_recal90_${chrom}_nohead.vcf | perl ~/scripts/vcfsorter.pl /scratch/justin.zook/references/human_g1k_v37.dict - > ${filestart}_${dataset}call_UGHapMerge_Het_AB_recal90_${chrom}.vcf";
	 print "$samCmd\n";
	$res = system($samCmd);
	$samCmd = "rm ${filestart}_indel_${dataset}call_UGHapMerge_Het_AB_recal90_${chrom}_nohead.vcf";
	print "$samCmd\n";
	 $res = system($samCmd);
} 

 if ($res == 0) {
	my $samCmd = "grep -v '^\\#' ${filestart}_indel_${dataset}call_UGHapMerge_Het_Align_recal90_${chrom}.vcf > ${filestart}_indel_${dataset}call_UGHapMerge_Het_Align_recal90_${chrom}_nohead.vcf";
	  print "$samCmd\n";
	 $res = system($samCmd);
	 $samCmd = "cat ${filestart}_snp_${dataset}call_UGHapMerge_Het_Align_recal90_${chrom}.vcf ${filestart}_indel_${dataset}call_UGHapMerge_Het_Align_recal90_${chrom}_nohead.vcf | perl ~/scripts/vcfsorter.pl /scratch/justin.zook/references/human_g1k_v37.dict - > ${filestart}_${dataset}call_UGHapMerge_Het_Align_recal90_${chrom}.vcf";
	 print "$samCmd\n";
	$res = system($samCmd);
	$samCmd = "rm ${filestart}_indel_${dataset}call_UGHapMerge_Het_Align_recal90_${chrom}_nohead.vcf";
	print "$samCmd\n";
	 $res = system($samCmd);
} 

 if ($res == 0) {
	my $samCmd = "grep -v '^\\#' ${filestart}_indel_${dataset}call_UGHapMerge_Het_map_recal90_${chrom}.vcf > ${filestart}_indel_${dataset}call_UGHapMerge_Het_map_recal90_${chrom}_nohead.vcf";
	  print "$samCmd\n";
	 $res = system($samCmd);
	 $samCmd = "cat ${filestart}_snp_${dataset}call_UGHapMerge_Het_map_recal90_${chrom}.vcf ${filestart}_indel_${dataset}call_UGHapMerge_Het_map_recal90_${chrom}_nohead.vcf | perl ~/scripts/vcfsorter.pl /scratch/justin.zook/references/human_g1k_v37.dict - > ${filestart}_${dataset}call_UGHapMerge_Het_map_recal90_${chrom}.vcf";
	 print "$samCmd\n";
	$res = system($samCmd);
	$samCmd = "rm ${filestart}_indel_${dataset}call_UGHapMerge_Het_map_recal90_${chrom}_nohead.vcf";
	print "$samCmd\n";
	 $res = system($samCmd);
} 


