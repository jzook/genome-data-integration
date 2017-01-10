# Python script to prepare gatk-genotype-gvcfs-v3.5-anyref files

#Adjust file name for desired output name
out = open("160922_gatk-haplotypecaller-v35-anyref_10X_.sh", "w") 

#Set all of these values as appropriate 
HG = "HG001"		#genome
plat= "10XGenomics"    #platform
abv = "10X"        #abreviated platform name

chr = 1 
for i in range(22):  
	
	out.write("dx run -y GIAB:/Workflow/GATK_V3.5/gatk-haplotypecaller-v35-anyref " +
	 	
		"-isorted_bam=/"+HG+"/GRCh38/"+ plat +"/" + HG + "_GRCh38_" + abv + "_HP1_"  + str(chr) + ".bam " +
		"-isorted_bai=/"+HG+"/GRCh38/"+ plat +"/" + HG + "_GRCh38_" + abv + "_HP1_"  + str(chr) + ".bai " +
		"-ioutput_prefix=" + HG + "_" + str(chr) + "_" + abv + "_HP1_" + "GATKHC_gvcf " +
		"-iref=/assets/GRCh38hs38d1noalt.fasta-index.tar.gz " +
		"-iextra_options=\"-L chr" + str(chr) + " -stand_call_conf 2 -stand_emit_conf 2 -A BaseQualityRankSumTest -A ClippingRankSumTest -A Coverage -A FisherStrand -A LowMQ -A RMSMappingQuality -A ReadPosRankSumTest -A StrandOddsRatio -A HomopolymerRun -A TandemRepeatAnnotator\" " +
		"--destination=/" + HG+ "/" + "GRCh38/" + plat + "/GATKHC_output/" +	
		"\n") 
	
	chr = chr+1 
	

chr = 1 
for i in range(22):  
	
	out.write("dx run -y GIAB:/Workflow/GATK_V3.5/gatk-haplotypecaller-v35-anyref " +
	 	
		"-isorted_bam=/"+HG+"/GRCh38/"+ plat +"/" + HG + "_GRCh38_" + abv + "_HP2_"  + str(chr) + ".bam " +
		"-isorted_bai=/"+HG+"/GRCh38/"+ plat +"/" + HG + "_GRCh38_" + abv + "_HP2_"  + str(chr) + ".bai " +
		"-ioutput_prefix=" + HG + "_" + str(chr) + "_" + abv + "_HP2_" + "GATKHC_gvcf " +
		"-iref=/assets/GRCh38hs38d1noalt.fasta-index.tar.gz " +
		"-iextra_options=\"-L chr" + str(chr) + " -stand_call_conf 2 -stand_emit_conf 2 -A BaseQualityRankSumTest -A ClippingRankSumTest -A Coverage -A FisherStrand -A LowMQ -A RMSMappingQuality -A ReadPosRankSumTest -A StrandOddsRatio -A HomopolymerRun -A TandemRepeatAnnotator\" " +
		"--destination=/" + HG+ "/" + "GRCh38/" + plat + "/GATKHC_output/" +	
		"\n") 
	
	chr = chr+1 

out.close()	

