# Python script to prepare integration commands for chr 1-22
# This script was based on integration v3.3 calls and callset tables
# To run :  $ python Integration_batch_script.sh
# Integration script file output to directory script was run in

#Adjust file name for desired output name
out = open("161014_integration_hg3_v3.3.1.sh", "w") 

#Set all of these values as appropriate for current integration
int="v.3.3.1"                              #Integration version
HG="HG003" 								#genome number
intdate= "161014" 						#date integration is being run, yymmdd
CG="GS000037264" 						#identifier for Complete Genomics dataset
ion="20141120.NA24149" 					        #identifier for Ion dataset
callsets= "CG-IllFB-IllGATKHC-Ion-10X"  #callsets used for integration
ref= "/assets/hs37d5.fasta-index.tar.gz"
rtgsdf= "/assets/rtgsdf37.tar.gz"
refn= "/filtbeds/GRCh37/example_of_no_ref_regions_input_file_b37.bed"

#confirm these annotation files are being used for current integration, if not add or delete files and adjust in iannot section in BOTH loops
iannot1= "-iannotations=/Annotation_files/GATK_Annotations_160509.txt"
iannot2= "-iannotations=/Annotation_files/Freebayes_Annotations_160408.txt"
iannot3= "-iannotations=/Annotation_files/CG_Annotations_160408.txt"
iannot4= "-iannotations=/Annotation_files/Ion_Annotations_160408.txt"

#set path to desired callset table
cstable= "-icallsettable=/callset_tables/161014_HG003_Callset_tables_CG_Ilmn_Ion_10X_for_GRCh37_v3.3.1/HG003_RM8391_Datasets_CG_Illumina_Ion_10X_Files_GRCh37_"

#confirm bed files for filtering are being used for current integration, if not add or delete files and adjust in ifiltbed section in BOTH loops
ifiltbed1= "-ifiltbeds=/HG003/GRCh37/HG002_HG003_HG004_allsvs_merged.bed"
ifiltbed2= "-ifiltbeds=/filtbeds/GRCh37/AllRepeats_lt51bp_gt95identity_merged_slop5.bed.gz"
ifiltbed3= "-ifiltbeds=/filtbeds/GRCh37/AllRepeats_51to200bp_gt95identity_merged_slop5.bed.gz"
ifiltbed4= "-ifiltbeds=/filtbeds/GRCh37/AllRepeats_gt200bp_gt95identity_merged_sort.bed"
ifiltbed5= "-ifiltbeds=/filtbeds/GRCh37/hg19_self_chain_split_withalts_gt10k.bed.gz"
ifiltbed6= "-ifiltbeds=/filtbeds/GRCh37/SimpleRepeat_imperfecthomopolgt10_slop5.bed.gz"
ifiltbed7= "-ifiltbeds=/filtbeds/GRCh37/superdupsmerged_all_sort.bed"
ifiltbed8= "-ifiltbeds=/filtbeds/GRCh37/mm-2-merged.bed"

#Confirm the following callsets are being used for current integration: 
	 # CG
	 # Ilmn 150x150 300X GATKHC
	 # Ilmn 150x150 FB
	 # IonExome
	 # 10X	 		
     # Ilmn 250x250 GATKHC
	 # Ilmn 250x250 FB
	 # Ilmn MatePair GATKHC
	 # Ilmn MatePair FB
# If callsets vary delete and or add vcf/beds to BOTH loops
	 
# First for loop for chromosome 1-15 and has call for additional memory, adjust as neccessary
chr = 1 
for i in range(15):  
	
	out.write("dx run -y GIAB:/Workflow/nist-integration-v3.3.1-anyref "+
	 	
		"-ivcfs=/"+HG+"/GRCh37/Complete_Genomics/Integration_prepare_cg_output/vcfBeta-"+CG+"-ASM_"+str(chr)+".vcf.gz " + 
		"-ivcfs=/"+HG+"/GRCh37/Illumina/Illumina_hs37d5_150x150/GATKHC_output/"+HG+"_"+str(chr)+"_hs37d5_novoalign_Ilmn150bp300X_GATKHC.vcf.gz " +
		"-ivcfs=/"+HG+"/GRCh37/Illumina/Illumina_hs37d5_150x150/FreeBayes_output/"+HG+"_"+str(chr)+"_hs37d5_novoalign_Ilmn150bp300X_FB.vcf.gz " + 
		"-ivcfs=/"+HG+"/GRCh37/IonExome/Integration_prepare_ion_output/AmpliseqExome." + ion + "_" + str(chr) + ".vcf.gz " +
		"-ivcfs=/"+HG+"/GRCh37/10XGenomics/Integration_prepare_10X_output/"+HG+"_"+str(chr)+"_10X_GATKHCbyhaplo.vcf.gz "+
		"-ivcfs=/"+HG+"/GRCh37/Illumina/Illumina_hs37d5_250x250/GATKHC_output/"+HG+"_"+str(chr)+"_hs37d5_novoalign_Ilmn250x250_GATKHC.vcf.gz " +
		"-ivcfs=/"+HG+"/GRCh37/Illumina/Illumina_hs37d5_250x250/FreeBayes_output/"+HG+"_"+str(chr)+"_hs37d5_novoalign_Ilmn250x250_FB.vcf.gz " +
		"-ivcfs=/"+HG+"/GRCh37/Illumina/Illumina_hg19_MatePair/GATKHC_output/"+HG+"_"+str(chr)+"_hg19_bwa_mem_IlmnMatePair_GATKHC.vcf.gz " +
		"-ivcfs=/"+HG+"/GRCh37/Illumina/Illumina_hg19_MatePair/FreeBayes_output/"+HG+"_"+str(chr)+"_hg19_bwa_mem_IlmnMatePair_FB.vcf.gz " +
		#"-ivcfs=/"+HG+"/GRCh37/SOLID/GATKHC_output/"+ HG + "_" + str(chr) + "_hg19_solid5500_PE50x50bp_GATKHC.vcf.gz "+
		#"-ivcfs=/"+HG+"/GRCh37/SOLID/GATKHC_output/"+ HG + "_" + str(chr) + "_hg19_solid5500_SE75bp_GATKHC.vcf.gz "+


		"-ibeds=/"+HG+"/GRCh37/Complete_Genomics/Integration_prepare_cg_output/vcfBeta-"+CG+"-ASM_callable_"+str(chr)+".bed " +
		"-ibeds=/"+HG+"/GRCh37/Illumina/Illumina_hs37d5_150x150/Integration_prepare_GATKHC_v.3.3.1/" + HG + "_" + str(chr) +"_hs37d5_novoalign_Ilmn150bp300X_GATKHC_gvcf_callable.bed " +
		"-ibeds=/"+HG+"/GRCh37/Illumina/Illumina_hs37d5_250x250/Integration_prepare_GATKHC_v.3.3.1/" + HG + "_" + str(chr) +"_hs37d5_novoalign_Ilmn250x250_GATKHC_gvcf_callable.bed " +
		"-ibeds=/"+HG+"/GRCh37/Illumina/Illumina_hg19_MatePair/Integration_prepare_GATKHC_v.3.3.1/" + HG + "_" + str(chr) +"_hg19_bwa_mem_IlmnMatePair_GATKHC_gvcf_callable.bed " +
		"-ibeds=/"+HG+"/GRCh37/IonExome/Integration_prepare_ion_output/AmpliseqExome." + ion + "_callable_" + str(chr) + ".bed " +
		"-ibeds=/"+HG+"/GRCh37/10XGenomics/Integration_prepare_10X_output/"+HG+"_"+str(chr)+"_10X_GATKHCbyhaplo_callable.bed " +
		"-ibeds=/"+HG+"/GRCh37/Illumina/Illumina_hs37d5_150x150/CallableLoci_output/"+HG+"_"+str(chr)+"_hs37d5_novoalign_Ilmn150bp300X_callableloci.bed " +
		"-ibeds=/"+HG+"/GRCh37/Illumina/Illumina_hs37d5_250x250/CallableLoci_output/"+HG+"_"+str(chr)+"_hs37d5_novoalign_Ilmn250x250_callableloci.bed " +
		"-ibeds=/"+HG+"/GRCh37/Illumina/Illumina_hg19_MatePair/CallableLoci_output/"+HG+"_"+str(chr)+"_hg19_bwa_mem_IlmnMatePair_callableloci.bed "+
		 
		iannot1 + " " + iannot2 + " " + iannot3 + " " + iannot4 + " " +
		 
		cstable + str(chr) + ".txt" + " " +
		 
		ifiltbed1 + " " + ifiltbed2 + " " + ifiltbed3 + " " + ifiltbed4 + " " + ifiltbed5 + " " + ifiltbed6 + " " + ifiltbed7 + " " + ifiltbed8 + " " +
		
		"-iref=" + ref + " " + "-irtgsdf=" + rtgsdf + " " + "-irefn=" + refn + " " +
		
		"-ichrom=" + str(chr) + " " + "-iprefix="+ HG + "_GIAB_highconf_"+ callsets + "_" + str(chr) + "_" + int + " " + "--destination=/" + HG + "/GRCh37/Integration_"+ int + "_output/" + intdate + "_"+ callsets + "_" + int + "/ " +
		
		"--instance-type=mem3_hdd2_x2"+ "\n")
	
	chr = chr+1 

# Second for loop for smaller chromosomes 16-22, default memory used
chr = 16
for i in range(7):  
	
	out.write("dx run -y GIAB:/Workflow/nist-integration-v3.3.1-anyref "+
	 	
		"-ivcfs=/"+HG+"/GRCh37/Complete_Genomics/Integration_prepare_cg_output/vcfBeta-"+CG+"-ASM_"+str(chr)+".vcf.gz " + 
		"-ivcfs=/"+HG+"/GRCh37/Illumina/Illumina_hs37d5_150x150/GATKHC_output/"+HG+"_"+str(chr)+"_hs37d5_novoalign_Ilmn150bp300X_GATKHC.vcf.gz " +
		"-ivcfs=/"+HG+"/GRCh37/Illumina/Illumina_hs37d5_150x150/FreeBayes_output/"+HG+"_"+str(chr)+"_hs37d5_novoalign_Ilmn150bp300X_FB.vcf.gz " + 
		"-ivcfs=/"+HG+"/GRCh37/IonExome/Integration_prepare_ion_output/AmpliseqExome." + ion + "_" + str(chr) + ".vcf.gz " +
		"-ivcfs=/"+HG+"/GRCh37/10XGenomics/Integration_prepare_10X_output/"+HG+"_"+str(chr)+"_10X_GATKHCbyhaplo.vcf.gz "+
		"-ivcfs=/"+HG+"/GRCh37/Illumina/Illumina_hs37d5_250x250/GATKHC_output/"+HG+"_"+str(chr)+"_hs37d5_novoalign_Ilmn250x250_GATKHC.vcf.gz " +
		"-ivcfs=/"+HG+"/GRCh37/Illumina/Illumina_hs37d5_250x250/FreeBayes_output/"+HG+"_"+str(chr)+"_hs37d5_novoalign_Ilmn250x250_FB.vcf.gz " +
		"-ivcfs=/"+HG+"/GRCh37/Illumina/Illumina_hg19_MatePair/GATKHC_output/"+HG+"_"+str(chr)+"_hg19_bwa_mem_IlmnMatePair_GATKHC.vcf.gz " +
		"-ivcfs=/"+HG+"/GRCh37/Illumina/Illumina_hg19_MatePair/FreeBayes_output/"+HG+"_"+str(chr)+"_hg19_bwa_mem_IlmnMatePair_FB.vcf.gz " +
		#"-ivcfs=/"+HG+"/GRCh37/SOLID/GATKHC_output/"+ HG + "_" + str(chr) + "_hg19_solid5500_PE50x50bp_GATKHC.vcf.gz "+
		#"-ivcfs=/"+HG+"/GRCh37/SOLID/GATKHC_output/"+ HG + "_" + str(chr) + "_hg19_solid5500_SE75bp_GATKHC.vcf.gz "+


		"-ibeds=/"+HG+"/GRCh37/Complete_Genomics/Integration_prepare_cg_output/vcfBeta-"+CG+"-ASM_callable_"+str(chr)+".bed " +
		"-ibeds=/"+HG+"/GRCh37/Illumina/Illumina_hs37d5_150x150/Integration_prepare_GATKHC_v.3.3.1/" + HG + "_" + str(chr) +"_hs37d5_novoalign_Ilmn150bp300X_GATKHC_gvcf_callable.bed " +
		"-ibeds=/"+HG+"/GRCh37/Illumina/Illumina_hs37d5_250x250/Integration_prepare_GATKHC_v.3.3.1/" + HG + "_" + str(chr) +"_hs37d5_novoalign_Ilmn250x250_GATKHC_gvcf_callable.bed " +
		"-ibeds=/"+HG+"/GRCh37/Illumina/Illumina_hg19_MatePair/Integration_prepare_GATKHC_v.3.3.1/" + HG + "_" + str(chr) +"_hg19_bwa_mem_IlmnMatePair_GATKHC_gvcf_callable.bed " +
		"-ibeds=/"+HG+"/GRCh37/IonExome/Integration_prepare_ion_output/AmpliseqExome." + ion + "_callable_" + str(chr) + ".bed " +
		"-ibeds=/"+HG+"/GRCh37/10XGenomics/Integration_prepare_10X_output/"+HG+"_"+str(chr)+"_10X_GATKHCbyhaplo_callable.bed " +
		"-ibeds=/"+HG+"/GRCh37/Illumina/Illumina_hs37d5_150x150/CallableLoci_output/"+HG+"_"+str(chr)+"_hs37d5_novoalign_Ilmn150bp300X_callableloci.bed " +
		"-ibeds=/"+HG+"/GRCh37/Illumina/Illumina_hs37d5_250x250/CallableLoci_output/"+HG+"_"+str(chr)+"_hs37d5_novoalign_Ilmn250x250_callableloci.bed " +
		"-ibeds=/"+HG+"/GRCh37/Illumina/Illumina_hg19_MatePair/CallableLoci_output/"+HG+"_"+str(chr)+"_hg19_bwa_mem_IlmnMatePair_callableloci.bed "+
		 
		iannot1 + " " + iannot2 + " " + iannot3 + " " + iannot4 + " " +
		 
		cstable + str(chr) + ".txt" + " " +
		 
		ifiltbed1 + " " + ifiltbed2 + " " + ifiltbed3 + " " + ifiltbed4 + " " + ifiltbed5 + " " + ifiltbed6 + " " + ifiltbed7 + " " + ifiltbed8 + " " +
		
		"-iref=" + ref + " " + "-irtgsdf=" + rtgsdf + " " + "-irefn=" + refn + " " +
		
		"-ichrom=" + str(chr) + " " + "-iprefix="+ HG + "_GIAB_highconf_"+ callsets + "_" + str(chr) + "_" + int + " " + "--destination=/" + HG + "/GRCh37/Integration_"+ int + "_output/" + intdate + "_"+ callsets + "_" + int + "/ " +
		
		"\n") 
		
	chr = chr+1 	
out.close()	
