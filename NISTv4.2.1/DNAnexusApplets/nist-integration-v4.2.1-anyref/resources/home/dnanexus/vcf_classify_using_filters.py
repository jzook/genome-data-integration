#
# VcfClassifyUsingFilters_v3.pl - arbitrate uncertain genotype calls using filters and callable regions, determined differently for each technology
# v3.1 - if genotypes disagree between callsets and the arbitration decides to use 0/1, then make it uncertain, since CG and freebayes have this type of genotyping error
# v3.1 - remove callsettable column giving option to use only to confirm variants, since this can be done simply by using a notcallable bed containing the whole genome
#      - also filter implied homref records if another callset for the same dataset is filtered
# v3.1.1 - filter ion if it doesn't make a call at indels called by other methods, since it misses some
#        - don't report confident homref from any callset if it doesn't make a call at indels >10bp called by other methods, since none are perfectly sensitive to larger indels
# v3.3 - Changed not callable to callable bed file annotations
#      - added ADsum and alleleimbalance FILTER
#      - now take phasing and ID from GATK PGT and PID fields
# v3.3.1 - correct PS type to String in header
#        - enable GRCh38 or GRCh37
#        - filter sites where Ion misses indel even if all others are not callable
# v3.3.2 - filter sites that are within 50bp of another passing call but none of the callsets that support the 2 calls match, because some nearby conflicting calls from different callers were both considered high confidence if another callset from the same dataset was filtered
#        - transfer difficultregion annotation from input to output files
#		 - Fix interpretation of CG AD field for homozygous sites
import numpy
import glob
import re
import math

class ClassifyUsingFilter(object):

	#create header
	header = None

	genome_file = None
	genome_file_lines = None	
	fields = None
	len1 = None
	parser = None
	infile = None
	callsettable = None
	outfilestart = None
	
	vcfall_file = None
	callsettable_file = None
	output_allcalls = None
	output_arbitrated_file = None
	output_2_platforms_file = None
	output_test_file = None
	platform = None
	dataset = None
	callset = None
	fields = None
	callsettable_file_lines = None
	callsets = None

	x = None
	vcfall_file_lines = None
	vcfheader_line_idx = None
	agreecallsetsprev = None
	prevpasspos = None

	fields = None
	chrom = None
	pos = None
	ref = None
	alt = None
	qual = None
	info = None
	format_field = None
	formats = None
	callable_field = None
	filt = None
	callabletxt = None
	filttxt = None
	gt = None
	gq = None
	gl = None
	dp = None
	ad1 = None
	ad2 = None
	ad3 = None
	gtphased = None
	ps = None
	datasetfilttxt = None
	difficultregion = None
	infosplit = None

	gt1 = None
	agree = None
	disagree = None
	agreecallsets = None #list of callsets that agree with gt1 regardless of filter
	agreecallsetscnt = None #number of callsets that agree with gt1 regardless of filter
	agreedatasets = None #list of datasets that agree with gt1 regardless of filter
	agreedatasetscnt = None #number of datasets that agree with gt1 regardless of filter
	agreeplatforms = None #list of platforms that agree with gt1 regardless of filter
	agreeplatformscnt = None #number of platforms that agree with gt1 regardless of filter
	notcallablecnt = None #number of call sets that are not callable 
	gqsum = None #sum GQ's of the datasets with the final genotype
	DPSum = None #sum DP's of the datasets with the final genotype
	ADSum1 = None #sum AD's of the datasets with the final genotype and not filtered
	ADSum2 = None #sum AD's of the datasets with the final genotype and not filtered
	ADSum3 = None #sum AD's of the datasets with the final genotype and not filtered
	ADallSum1 = None #sum AD's of the datasets with the final genotype
	ADallSum2 = None #sum AD's of the datasets with the final genotype
	ADallSum3 = None #sum AD's of the datasets with the final genotype

	gt1nofilt = None #genotype of first callset without low cov or filter
	agreenofilt = None #number of callsets that agree with gt1nofilt and not filtered
	disagreenofilt = None #number of callsets that disagree with gt1nofilt and not filtered 
	notcallablefiltcnt = None #number of call sets with low cov or genotype with . or filter
	agreecallsetsnofilt = None #list of datasets that agree with gt1nofilt and not filtered
	agreecallsetsnofiltcnt = None #number of datasets that agree with gt1nofilt and not filtered
	agreedatasetsnofilt = None #list of datasets that agree with gt1nofilt and not filtered
	agreedatasetsnofiltcnt = None #number of datasets that agree with gt1nofilt and not filtered
	agreeplatformsnofilt = None #list of platforms that agree with gt1nofilt and not filtered
	agreeplatformsnofiltcnt = None #number of platforms that agree with gt1nofilt and not filtered
	filtout = None #output text to FILTER field
	arbitratetxt = None #output to INFO if arbitrated

	agreedatasetsnofilt = None #list of datasets that agree with gt1nofilt and not filtered; reset so that it can be used below to determine if AD's have been added from this dataset yet
	agreedatasetsnofiltcnt = None #number of datasets that agree with gt1nofilt and not filtered
	disagreecallsetsnofilt = None
	disagreedatasetsnofilt = None
	disagreeplatformsnofilt = None
	disagreecallsetsnofiltcnt = None
	disagreedatasetsnofiltcnt = None
	disagreeplatformsnofiltcnt = None

	currentcallsets = None
	matchcallset = None
	matchcallset = None
	prevpasspos = None
	agreecallsetsprev = None
			
	ADSumout = None
	ADallSumout = None
	disagreetxt = None
			
	qualout = None

	def __init__(self, unionvcf, CallsetTable, outputFileStart):
		self.header = ("##fileformat=VCFv4.2\n"
		"##FILTER=<ID=GQlessthan70,Description=\"Sum of GQ for datasets with this genotype less than 70\">\n"
		"##FILTER=<ID=allfilteredanddisagree,Description=\"All callsets have this call filtered or outside the callable regions and they have discordant genotypes or variant calls\">\n"
		"##FILTER=<ID=allfilteredbutagree,Description=\"All callsets have this call filtered or outside the callable regions but they have the same genotype\">\n"
		"##FILTER=<ID=discordantunfiltered,Description=\"Callsets with unfiltered calls have discordant genotypes or variant calls\">\n"
		"##FILTER=<ID=discordanthet,Description=\"Filtered calls where a passing call is het and a high GQ but filtered call is hom var, since often the het is wrong\">\n"
		"##FILTER=<ID=questionableindel,Description=\"Filtered calls where some callsets have a filtered indel larger than 10bp and another dataset has an implied homozygous reference call\">\n"
		"##FILTER=<ID=cgonly,Description=\"Filtered calls where only Complete Genomics had this call and it was completely missing from any other callset\">\n"
		"##FILTER=<ID=alleleimbalance,Description=\"Filtered calls where the net allele balance for unfiltered datasets is <0.2 or >0.8\">\n"
		"##FILTER=<ID=overlappingcall,Description=\"Filtered sites that are within 50bp of another passing call but none of the callsets that support the 2 calls match\">\n"
		"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth summed across all datasets, excluding MQ0 reads\">\n"
		"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Net Genotype quality across all datasets, calculated from GQ scores of callsets supporting the consensus GT, using only one callset from each dataset\">\n"
		"##FORMAT=<ID=ADALL,Number=R,Type=Integer,Description=\"Net allele depths across all datasets\">\n"
		"##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Net allele depths across all unfiltered datasets with called genotype\">\n"
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Consensus Genotype across all datasets with called genotype\">\n"
		"##FORMAT=<ID=PS,Number=1,Type=String,Description=\"Phase set in which this variant falls\">\n"
		"##INFO=<ID=DPSum,Number=1,Type=Integer,Description=\"Total read depth summed across all datasets, excluding MQ0 reads\">\n"
		"##INFO=<ID=platforms,Number=1,Type=Integer,Description=\"Number of different platforms for which at least one callset called this genotype, whether filtered or not\">\n"
		"##INFO=<ID=platformnames,Number=.,Type=String,Description=\"Names of platforms for which at least one callset called this genotype, whether filtered or not\">\n"
		"##INFO=<ID=platformbias,Number=.,Type=String,Description=\"Names of platforms that have reads containing a variant at this location, but the high-confidence call is homozygous reference, indicating that there is a potential bias.\">\n"
		"##INFO=<ID=datasets,Number=1,Type=Integer,Description=\"Number of different datasets for which at least one callset called this genotype, whether filtered or not\">\n"
		"##INFO=<ID=datasetnames,Number=.,Type=String,Description=\"Names of datasets for which at least one callset called this genotype, whether filtered or not\">\n"
		"##INFO=<ID=datasetsmissingcall,Number=.,Type=String,Description=\"Names of datasets that are missing a call or have an incorrect call at this location, and the high-confidence call is a variant\">\n"
		"##INFO=<ID=callsets,Number=1,Type=Integer,Description=\"Number of different callsets that called this genotype, whether filtered or not\">\n"
		"##INFO=<ID=callsetnames,Number=.,Type=String,Description=\"Names of callsets that called this genotype, whether filtered or not\">\n"
		"##INFO=<ID=varType,Number=1,Type=String,Description=\"Type of variant\">\n"
		"##INFO=<ID=filt,Number=.,Type=String,Description=\"List of callsets that had this call filtered.\">\n"
		"##INFO=<ID=callable,Number=.,Type=String,Description=\"List of callsets that had this call in a region with low coverage of high MQ reads.\">\n"
		"##INFO=<ID=difficultregion,Number=.,Type=String,Description=\"List of difficult region bed files containing this call.\">\n"
		"##INFO=<ID=arbitrated,Number=1,Type=String,Description=\"TRUE if callsets had discordant calls so that arbitration was needed.\">\n"
		"##INFO=<ID=callsetwiththisuniqgenopassing,Number=.,Type=String,Description=\"Callset that uniquely calls the PASSing genotype in GT when 2+ PASSing callsets support a different genotype.\">\n"
		"##INFO=<ID=callsetwithotheruniqgenopassing,Number=.,Type=String,Description=\"Callset that uniquely calls a PASSing genotype different from GT when 2+ PASSing callsets support the genotype in GT.\">\n")


		# self.genome_file = open("human.genome", "r")
		# self.genome_file_lines = genome_file.readlines()	
		# self.fields = line.split("\t")
		self.fields = None
		self.len1 = ""
		self.parser = None
		self.infile = unionvcf
		self.callsettable = CallsetTable
		self.outfilestart = outputFileStart
	
		# self.vcfall_file = open(vcfall, "r")
		# self.callsettable_file = open(callsettable, "r")
		self.output_allcalls = None
		self.output_arbitrated_file = None
		self.output_2_platforms_file = None
		self.output_test_file = None
		self.platform = None
		self.dataset = None
		self.callset = None
		self.callsettable_file_lines = None
		self.callsets = None

		self.x = None
		self.vcfall_file_lines = None
		self.vcfheader_line_idx = None
		self.agreecallsetsprev = None
		self.prevpasspos = None

		self.fields = None
		self.chrom = None
		self.pos = None
		self.id = None
		self.ref = None
		self.alt = None
		self.qual = None
		self.info = None
		self.format_field = None
		self.formats = None
		self.callable_field = None
		self.filt = None
		self.callabletxt = None
		self.filttxt = None
		self.gt = None
		self.gq = None
		self.gl = None
		self.dp = None
		self.ad1 = None
		self.ad2 = None
		self.ad3 = None
		self.gtphased = None
		self.ps = None
		self.datasetfilttxt = None
		self.difficultregion = None
		self.infosplit = None

		self.gt1 = None
		self.agree = None
		self.disagree = None
		self.agreecallsets = None #list of callsets that agree with gt1 regardless of filter
		self.agreecallsetscnt = None #number of callsets that agree with gt1 regardless of filter
		self.agreedatasets = None #list of datasets that agree with gt1 regardless of filter
		self.agreedatasetscnt = None #number of datasets that agree with gt1 regardless of filter
		self.agreeplatforms = None #list of platforms that agree with gt1 regardless of filter
		self.agreeplatformscnt = None #number of platforms that agree with gt1 regardless of filter
		self.notcallablecnt = None #number of call sets that are not callable 
		self.gqsum = None #sum GQ's of the datasets with the final genotype
		self.DPSum = None #sum DP's of the datasets with the final genotype
		self.ADSum1 = None #sum AD's of the datasets with the final genotype and not filtered
		self.ADSum2 = None #sum AD's of the datasets with the final genotype and not filtered
		self.ADSum3 = None #sum AD's of the datasets with the final genotype and not filtered
		self.ADallSum1 = None #sum AD's of the datasets with the final genotype
		self.ADallSum2 = None #sum AD's of the datasets with the final genotype
		self.ADallSum3 = None #sum AD's of the datasets with the final genotype

		self.gt1nofilt = None #genotype of first callset without low cov or filter
		self.agreenofilt = None #number of callsets that agree with gt1nofilt and not filtered
		self.disagreenofilt = None #number of callsets that disagree with gt1nofilt and not filtered 
		self.notcallablefiltcnt = None #number of call sets with low cov or genotype with . or filter
		self.agreecallsetsnofilt = None #list of datasets that agree with gt1nofilt and not filtered
		self.agreecallsetsnofiltcnt = None #number of datasets that agree with gt1nofilt and not filtered
		self.agreedatasetsnofilt = None #list of datasets that agree with gt1nofilt and not filtered
		self.agreedatasetsnofiltcnt = None #number of datasets that agree with gt1nofilt and not filtered
		self.agreeplatformsnofilt = None #list of platforms that agree with gt1nofilt and not filtered
		self.agreeplatformsnofiltcnt = None #number of platforms that agree with gt1nofilt and not filtered
		self.filtout = None #output text to FILTER field
		self.arbitratetxt = None #output to INFO if arbitrated

		self.agreedatasetsnofilt = None #list of datasets that agree with gt1nofilt and not filtered; reset so that it can be used below to determine if AD's have been added from this dataset yet
		self.agreedatasetsnofiltcnt = None #number of datasets that agree with gt1nofilt and not filtered
		self.disagreecallsetsnofilt = None
		self.disagreedatasetsnofilt = None
		self.disagreeplatformsnofilt = None
		self.disagreecallsetsnofiltcnt = None
		self.disagreedatasetsnofiltcnt = None
		self.disagreeplatformsnofiltcnt = None

		self.currentcallsets = None
		self.matchcallset = None
		self.matchcallset = None
		self.prevpasspos = None
		self.agreecallsetsprev = None
				
		self.ADSumout = None
		self.ADallSumout = None
		self.disagreetxt = None
				
		self.qualout = None

	def parse_input_files(self):
		try:
			self.genome_file = open("human.genome", "r")
		except:
			exit()
		self.genome_file_lines = self.genome_file.readlines()	
		for line in self.genome_file_lines:
			self.fields = line.split("\t")
			self.len1 = ""
			if re.search("(.*)\n", self.fields[1]): # if phased, convert from | to / for below analyses
				self.len1 = re.split("(.*)\n",self.fields[1])[1]
			else:
				self.len1 = self.fields[1]
			print(self.fields)
			self.header = self.header + "##contig=<ID=" + str(self.fields[0]) + ",length=" + str(self.fields[1].strip()) + ">\n"

		self.genome_file.close()
		self.header = self.header + "##fileDate=20160824\n#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	INTEGRATION\n"

		try:
			self.vcfall_file = open(self.infile, "r")
		except:
			exit()

		try:
			self.callsettable_file = open(self.callsettable, "r")
		except:
			exit()

		try:
			self.output_allcalls = open(self.outfilestart + "_ClassifyUsingFilters_allcalls.vcf", "w+")
		except:
			exit()

		try:
			self.output_arbitrated_file = open(self.outfilestart + "_ClassifyUsingFilters_arbitrated.vcf", "w+")
		except:
			exit()

		try:
			self.output_2_platforms_file = open(self.outfilestart + "_ClassifyUsingFilters_2platforms.vcf", "w+")
		except:
			exit()

		try:
			self.output_test_file = open(self.outfilestart + "_testout.txt", "w+")
		except:
			exit()


		#Read in information about each callset from callsettable
		self.platform = [""]*20
		self.dataset = [""]*20
		self.callset = [""]*20
		self.callsettable_file_lines = self.callsettable_file.readlines()
		self.callsets = 0

		cstable_idx = 0
		for line in self.callsettable_file_lines:
			if cstable_idx == 0:
				cstable_idx = cstable_idx + 1
				continue
			self.fields = line.split("\t")
			self.platform[self.callsets] = self.fields[0]
			self.dataset[self.callsets] = self.fields[1]
			self.callset[self.callsets] = self.fields[2]
			self.callsets = self.callsets + 1
			cstable_idx = cstable_idx + 1
		

		self.x = 0
		self.vcfall_file_lines = self.vcfall_file.readlines()
		self.vcfheader_line_idx = 0
		for line in self.vcfall_file_lines:
			self.vcfheader_line_idx = self.vcfheader_line_idx + 1
			if not(re.search("^\#\#", line)):
				self.x = 1 #skip header lines, but including line starting with #
			if self.x == 1:
				break
			

		self.output_allcalls.write(self.header)
		self.output_allcalls.flush()

		self.output_arbitrated_file.write(self.header)
		self.output_arbitrated_file.flush()

		self.output_2_platforms_file.write(self.header)
		self.output_2_platforms_file.flush()

		self.agreecallsetsprev = ""
		self.prevpasspos = -51

	def classify_and_filter(self):
		#Loop through union vcf and determine which calls are high confidence
		for vcf_line_idx in range(self.vcfheader_line_idx, len(self.vcfall_file_lines)):              
			# Split up the line into an array
			self.fields = self.vcfall_file_lines[vcf_line_idx].split("\t")
			#print("fields: " + str(self.fields))
			self.chrom = self.fields[0]
			self.pos = self.fields[1]
			self.id = self.fields[2]
			self.ref = self.fields[3]
			self.alt = self.fields[4]
			self.qual = self.fields[5]
			self.info = self.fields[7]
			self.format_field = self.fields[8]
			self.formats = self.format_field.split(":")
			self.callable_field = [0]*self.callsets
			self.filt = [0]*self.callsets
			self.callabletxt = ""
			self.filttxt = ""
			self.gt = [""]*self.callsets
			self.gq = ["-1.0"]*self.callsets
			self.gl = [""]*self.callsets
			self.dp = [0]*self.callsets
			self.ad1 = [0]*self.callsets
			self.ad2 = [0]*self.callsets
			self.ad3 = [0]*self.callsets
			self.gtphased = ""
			self.ps = "."
			self.minGQ = 70.0

			self.datasetfilttxt = ""
			self.difficultregion = ""

			self.infosplit = self.info.split(";")
			for k in range(0, len(self.infosplit)):
				if re.search("difficultregion", self.infosplit[k]):
					self.difficultregion=";" + str(self.infosplit[k])

			#Loop through callsets to determine if callable or filtered, and then find FORMAT annotations for each callset like GT, GQ, etc.
			for j in range(0, self.callsets):
				#is this call in the callset's callable regions?
				if re.search("CS_" + str(self.callset[j]) + "_callable", self.info): 
					self.callable_field[j] = 1
					if "DV" in self.callset[j]:
						self.minGQ = 25.0
					if self.callabletxt == "":
						self.callabletxt = ";callable=CS_" + str(self.callset[j]) + "_callable"
					else:
						self.callabletxt = self.callabletxt + "," + "CS_" + str(self.callset[j]) + "_callable"
				# 
				#is this call in the callset's filtered regions?
				if re.search("CS_" + str(self.callset[j]) + "_filt", self.info):
					self.filt[j] = 1
					if self.filttxt == "":
						self.filttxt = ";filt=CS_"+ str(self.callset[j]) + "_filt"
					else:
						self.filttxt = self.filttxt + ",CS_" + str(self.callset[j]) + "_filt"
					
					if self.datasetfilttxt == "":
						self.datasetfilttxt = ";datasetfilt=DS_" + str(self.dataset[j]) + "_filt"
					else:
						self.datasetfilttxt = self.datasetfilttxt + ",DS_" + str(self.dataset[j]) + "_filt"

				#find GT and GQ for this callset
				self.chars = self.fields[j+9].split(":")
				#print("chars: " + str(self.chars))
				if re.search("\.", self.chars[0]): #if no genotype, then don't need to find other fields
					self.gt[j] = "."
					self.gq[j] = "60.0"
					self.gl[j] = ""
				else:
					self.gt[j] = ""
					self.gq[j] = "60.0"
					self.gl[j] = ""
					self.revGT = 0 #are GT's not in increasing numerical order; if not, reverse AD's
					for k in range(0, len(self.chars)):
						#For GT field, format to (smaller number)/(larger number) to enable comparisons between callsets
						if self.formats[k] == "GT":
							if re.search("(.*)\n", self.chars[k]):
								self.gt[j] = re.split("(.*)\n", self.chars[k])[1]
							else: 
								self.gt[j] = self.chars[k]
							if re.search("(.*)\|(.*)", self.gt[j]): # if phased, convert from | to / for below analyses
								if j==0:
									self.gtphased = self.gt[j] #only keep phasing information from first callset
								self.gt[j] = re.split("(.*)\|(.*)", self.gt[j])[1] + "/" + re.split("(.*)\|(.*)", self.gt[j])[2]
							if re.search("(\d)\/(\d)", self.gt[j]): 
								if (re.split("(\d)\/(\d)", self.gt[j])[1] > re.split("(\d)\/(\d)", self.gt[j])[2]): 
									self.gt[j] = re.split("(\d)\/(\d)", self.gt[j])[2] + "/" + re.split("(\d)\/(\d)", self.gt[j])[1]
									self.revGT = 1  #order GT from smallest to largest number
							if re.search("\.", self.gt[j]): 
								self.gt[j] = "."
							if self.gt[j] == "1": 
								self.gt[j] = "1/1"

						#GATKHC outputs phased genotype in PGT
						if j == 0 and self.formats[k] == "PGT":
							if self.gt[j] == "1/1": 
								self.gtphased = "1|1" #account for bug in GATKHC that outputs 0|1 in PGT at hom ref sites
							else:
								if re.search("(.*)\n", self.chars[k]):  
									self.gtphased = re.split("(.*)\n", self.chars[k])[1]
								else:
									self.gtphased = self.chars[k]
							
						
						#get phase set from first callset (normally in PS but in PID for GATKHC)
						if j == 0 and (self.formats[k] == "PS" or self.formats[k] == "PID"): 
							# if "(.*)\n" in self.chars[k]:
							if re.search("(.*)\n", self.chars[k]):
								self.ps = re.split("(.*)\n", self.chars[k])[1]
							else:
								self.ps = self.chars[k]

						#get genotype quality
						if self.formats[k] == "GQ":
							if re.search("(.*)\n", self.chars[k]): 
								self.gq[j]=	re.split("(.*)\n", self.chars[k])[1]
							else:
								self.gq[j] = self.chars[k]
							if self.gq[j] == ".":
								self.gq[j] = "-1.0"

						
						#get coverage from DP
						if self.formats[k] == "DP":
							if re.search("(.*)\n", self.chars[k]): 
								self.dp[j] = re.split("(.*)\n",self.chars[k])[1]
							else:
								self.dp[j] = self.chars[k]
						
						#get allele depths, flip the order if we flipped the order of the GT field above
						if self.formats[k] == "AD":
							if re.search(".*,.*,.*", self.chars[k]): #if 2 ALT alleles, then 3 allele depths
								if re.search("(.*),(.*),(.*)\n", self.chars[k]):
									self.ad1[j] = re.split("(.*),(.*),(.*)\n", self.chars[k])[1]
									self.ad2[j] = re.split("(.*),(.*),(.*)\n", self.chars[k])[2]
									self.ad3[j] = re.split("(.*),(.*),(.*)\n", self.chars[k])[3]
								elif re.search("(.*),(.*),(.*)", self.chars[k]):
									self.ad1[j] = re.split("(.*),(.*),(.*)", self.chars[k])[1]
									self.ad2[j] = re.split("(.*),(.*),(.*)", self.chars[k])[2]
									self.ad3[j] = re.split("(.*),(.*),(.*)", self.chars[k])[3]
								if self.revGT==1:
									self.x = self.ad3[j]
									self.ad3[j] = self.ad2[j]
									self.ad2[j] = self.x
							else:
								if re.search("(.*),(.*)\n", self.chars[k]):
									self.ad1[j] = re.split("(.*),(.*)\n", self.chars[k])[1]
									self.ad2[j] = re.split("(.*),(.*)\n", self.chars[k])[2]
								elif re.search("(.*),(.*)", self.chars[k]):
									self.ad1[j] = re.split("(.*),(.*)", self.chars[k])[1]
									self.ad2[j] = re.split("(.*),(.*)", self.chars[k])[2]
								if self.revGT==1:
									self.x = self.ad1[j]
									self.ad1[j] = self.ad2[j]
									self.ad2[j] = self.x

			
			#first, check if all genotypes that are callable agree
			self.gt1 = "" #genotype of first callset without low cov
			self.agree = 0 #number of callsets that agree with gt1 regardless of filter
			self.disagree = 0 #number of callsets that disagree with gt1 regardless of filter 
			self.agreecallsets = "" #list of callsets that agree with gt1 regardless of filter
			self.agreecallsetscnt = 0 #number of callsets that agree with gt1 regardless of filter
			self.agreedatasets = "" #list of datasets that agree with gt1 regardless of filter
			self.agreedatasetscnt = 0 #number of datasets that agree with gt1 regardless of filter
			self.agreeplatforms = "" #list of platforms that agree with gt1 regardless of filter
			self.agreeplatformscnt = 0 #number of platforms that agree with gt1 regardless of filter
			self.notcallablecnt = 0 #number of call sets that are not callable 
			self.gqsum = 0 #sum GQ's of the datasets with the final genotype
			self.DPSum = 0 #sum DP's of the datasets with the final genotype
			self.ADSum1 = 0 #sum AD's of the datasets with the final genotype and not filtered
			self.ADSum2 = 0 #sum AD's of the datasets with the final genotype and not filtered
			self.ADSum3 = 0 #sum AD's of the datasets with the final genotype and not filtered
			self.ADallSum1 = 0 #sum AD's of the datasets with the final genotype
			self.ADallSum2 = 0 #sum AD's of the datasets with the final genotype
			self.ADallSum3 = 0 #sum AD's of the datasets with the final genotype

			#Loop through all callsets, check if they are callable, and check if callable callsets have the same genotype
			for j in range(0 , self.callsets):
				if self.gt1 == "": #Have we found any callable callsets yet?

					if ( self.callable_field[j] == 0 ) or ( not(self.gt[j] == ".") and float(self.gq[j]) < 20 ):  #is this not callable or (not homref and GQ<20)?
						self.notcallablecnt = self.notcallablecnt + 1
					else: #otherwise, it's callable
						self.gt1 = self.gt[j]
						self.agree = self.agree + 1
						self.agreecallsets = self.callset[j]
						self.agreecallsetscnt = self.agreecallsetscnt + 1
						self.agreedatasets= self.dataset[j]
						self.agreedatasetscnt = self.agreedatasetscnt + 1
						self.agreeplatforms = self.platform[j]
						self.agreeplatformscnt = self.agreeplatformscnt + 1
						self.gqsum = self.gqsum + float(self.gq[j])
				else: #if already found first callable callset, then compare to first callable callset's genotype
					if self.callable_field[j]==0 or (not(self.gt[j] == ".") and float(self.gq[j])<20): # #is this not callable or (not homref and GQ<20)?
						self.notcallablecnt = self.notcallablecnt + 1
					elif self.gt[j] == self.gt1: #otherwise, it's callable, so does it agree with the first callable callset's GT
						self.agree = self.agree + 1
						if (self.gqsum == 0 or not(re.search(str(self.dataset[j]), self.agreedatasets))): 
							self.gqsum = self.gqsum + float(self.gq[j]) #only add to gqsum once for each dataset
						if (not(re.search(str(self.dataset[j]), self.agreedatasets))):
							self.agreedatasets = self.agreedatasets + "," + self.dataset[j]
							self.agreedatasetscnt = self.agreedatasetscnt + 1
						if (not(re.search(str(self.callset[j]), self.agreecallsets))):
							self.agreecallsets = self.agreecallsets + "," + self.callset[j]
							self.agreecallsetscnt = self.agreecallsetscnt + 1
						if (not(re.search(str(self.platform[j]), self.agreeplatforms))):
							self.agreeplatforms = self.agreeplatforms + "," + self.platform[j]
							self.agreeplatformscnt = self.agreeplatformscnt + 1
					else: #otherwise, it disagrees with the first GT
						self.disagree = self.disagree + 1
			
			
			#Next, use filtering information in addition to callable information
			self.gt1nofilt = "" #genotype of first callset without low cov or filter
			self.agreenofilt = 0 #number of callsets that agree with gt1nofilt and not filtered
			self.disagreenofilt = 0 #number of callsets that disagree with gt1nofilt and not filtered 
			self.notcallablefiltcnt = 0 #number of call sets with low cov or genotype with . or filter
			self.agreecallsetsnofilt = "" #list of datasets that agree with gt1nofilt and not filtered
			self.agreecallsetsnofiltcnt = 0 #number of datasets that agree with gt1nofilt and not filtered
			self.agreedatasetsnofilt = "" #list of datasets that agree with gt1nofilt and not filtered
			self.agreedatasetsnofiltcnt = 0 #number of datasets that agree with gt1nofilt and not filtered
			self.agreeplatformsnofilt = "" #list of platforms that agree with gt1nofilt and not filtered
			self.agreeplatformsnofiltcnt = 0 #number of platforms that agree with gt1nofilt and not filtered
			self.filtout = "" #output text to FILTER field
			self.arbitratetxt = "" #output to INFO if arbitrated

			if (self.agree > 0 and self.disagree == 0):
				for j in range(0,self.callsets):
					if (self.gt[j] == self.gt1 and self.callable_field[j] > 0 and self.filt[j] == 0 and (self.gt[j] == "." or float(self.gq[j]) >= 20) and not((self.gt[j] == ".") and (re.search("DS_" + str(self.dataset[j]) + "_filt", self.datasetfilttxt)))):# { #is this callset's GT the same as the first callable callset's GT, is this one callable, is it not filtered, is GQ>=20, and (is it not (filtered in another callset for this dataset and homozygous reference))
						#TODO: Should probably reverse the logic here to match above and below, and also probably add the Ion indel filter from below
						self.agreenofilt = self.agreenofilt + 1
						if self.agreenofilt == 1:
							self.agreecallsetsnofilt = self.callset[j]
							self.agreecallsetsnofiltcnt = self.agreecallsetsnofiltcnt + 1
							self.agreedatasetsnofilt = self.dataset[j]
							self.agreedatasetsnofiltcnt = self.agreedatasetsnofiltcnt + 1
							self.agreeplatformsnofilt = self.platform[j]
							self.agreeplatformsnofiltcnt = self.agreeplatformsnofiltcnt + 1
						else:
							if (not(re.search(self.callset[j], self.agreecallsetsnofilt))):
								self.agreecallsetsnofilt = self.agreecallsetsnofilt + "," + self.callset[j]
								self.agreecallsetsnofiltcnt = self.agreecallsetsnofiltcnt + 1
							if (not(re.search(self.platform[j], self.agreeplatformsnofilt))):
								self.agreeplatformsnofilt = self.agreeplatformsnofilt + "," + self.platform[j]
								self.agreeplatformsnofiltcnt = self.agreeplatformsnofiltcnt + 1
								

				if self.gt1 == ".":
					self.gt1 = "0/0"
					self.filtout = "."
				
				if self.gt1 == "0/0" and self.agreenofilt>0 and (len(self.ref)>10 or len(self.alt)>10) or (self.agreeplatformsnofilt == "Ion" and (len(self.ref)>1 or len(self.alt)>1)):
					self.filtout = "questionableindel" #even if an unfiltered callset calls homozygous reference, make it uncertain if another callset has a filtered indel >10bp (or any indel if Ion misses it) because we found some callsets miss larger indels
				elif (self.gt1 == "0/0" and self.agreenofilt > 0):
					self.filtout = "."  #High-confidence homozygous reference if not a long indel
				elif self.agreenofilt > 0 and float(self.gqsum) < self.minGQ:
					self.filtout = "GQlessthan70" #if all callable, unfiltered datasets with this call sum to a GQ < 70, then make uncertain because there is insufficient support
				elif self.agreenofilt > 0:
					self.filtout = "PASS" #High-confidence variant if at least 1 unfiltered, callable callset
				else:
					self.filtout = "allfilteredbutagree" #otherwise, no unfiltered, callable callsets, so uncertain position
				

			elif self.agree > 0:
				#If all genotypes that were callable did not agree, then see if filtering information helps to arbitrate and determine which callset to trust
				self.gqsum=0
				for j in range(0,self.callsets):
					if (self.gt1nofilt == ""): #if haven't found first unfiltered call set yet
						if (self.callable_field[j]==0 or self.filt[j]>0 or (not(self.gt[j] == ".") and float(self.gq[j])<20) or ((self.gt[j] == ".") and (re.search("DS_" + str(self.dataset[j]) + "_filt", self.datasetfilttxt))) or ((self.gt[j] == ".") and (self.platform[j] == "Ion") and (len(self.ref)>1 or len(self.alt)>1))): #is this callset untrustworthy here
							self.notcallablefiltcnt = self.notcallablefiltcnt + 1
						else: #Otherwise, we trust it and get its genotype
							self.gt1nofilt = self.gt[j]
							self.agreenofilt = self.agreenofilt + 1
							self.agreecallsetsnofilt = self.callset[j]
							self.agreecallsetsnofiltcnt = self.agreecallsetsnofiltcnt + 1
							self.agreedatasetsnofilt = self.dataset[j]
							self.agreedatasetsnofiltcnt = self.agreedatasetsnofiltcnt + 1
							self.agreeplatformsnofilt = self.platform[j]
							self.agreeplatformsnofiltcnt = self.agreeplatformsnofiltcnt + 1
							self.gqsum = self.gqsum + float(self.gq[j])
				
					else:  #compare to first unfiltered genotype if it is callable and unfiltered
						if (self.callable_field[j]==0 or self.filt[j]>0 or (not(self.gt[j] == ".") and float(self.gq[j])<20) or ((self.gt[j] == ".") and (re.search("DS_" + str(self.dataset[j]) + "_filt", self.datasetfilttxt))) or ((self.gt[j] == ".") and (self.platform[j] == "Ion") and (len(self.ref)>1 or len(self.alt)>1))):  #is this callset untrustworthy here
							self.notcallablefiltcnt = self.notcallablefiltcnt + 1
						elif (self.gt[j] == self.gt1nofilt): #Otherwise, we trust it and compare to first trusted callset's genotype
							self.agreenofilt = self.agreenofilt + 1
							if (self.gqsum==0 or not(re.search(self.dataset[j], self.agreedatasetsnofilt))): 
								self.gqsum = self.gqsum + float(self.gq[j]) #only add to gqsum once for each dataset
							if (not(re.search(self.callset[j], self.agreecallsetsnofilt))):
								self.agreecallsetsnofilt = self.agreecallsetsnofilt + "," + self.callset[j]
								self.agreecallsetsnofiltcnt = self.agreecallsetsnofiltcnt + 1
							if (not(re.search(self.dataset[j], self.agreedatasetsnofilt))):
								self.agreedatasetsnofilt = self.agreedatasetsnofilt + "," + self.dataset[j]
								self.agreedatasetsnofiltcnt = self.agreedatasetsnofiltcnt + 1
							if (not(re.search(self.platform[j], self.agreeplatformsnofilt))):
								self.agreeplatformsnofilt = self.agreeplatformsnofilt + "," + self.platform[j]
								self.agreeplatformsnofiltcnt = self.agreeplatformsnofiltcnt + 1			
						else: #otherwise, this callset disagrees with the first callset
							self.disagreenofilt = self.disagreenofilt + 1
						
				if (self.gt1nofilt == "."):
					self.gt1 = "0/0"
				elif (not(self.gt1nofilt == "")):
					self.gt1 = self.gt1nofilt
				
				if (not(self.gt1nofilt == "")):
					self.agreecallsets = self.agreecallsetsnofilt
					self.agreecallsetscnt = self.agreecallsetsnofiltcnt
					self.agreedatasets = self.agreedatasetsnofilt
					self.agreedatasetscnt = self.agreedatasetsnofiltcnt
					self.agreeplatforms = self.agreeplatformsnofilt
					self.agreeplatformscnt = self.agreeplatformsnofiltcnt

				if (self.gt1 == "0/0" and self.agreenofilt>0 and ((len(self.ref)>10 or len(self.alt)>10) or (self.agreeplatformsnofilt == "Ion" and (len(self.ref)>1 or len(self.alt)>1)))):
					self.filtout= "questionableindel" #even if an unfiltered callset calls homozygous reference, make it uncertain if another callset has a filtered indel >10bp (or any indel if Ion misses it) because we found some callsets miss larger indels
				elif (self.gt1 == "0/0" and self.agreenofilt>0):
					self.filtout="." #High-confidence homozygous reference if not a long indel
				elif (self.agreenofilt>0 and float(self.gqsum)< self.minGQ):
					self.filtout = "GQlessthan70" #if all callable, unfiltered datasets with this call sum to a GQ < 70, then make uncertain because there is insufficient support
				elif (self.agreenofilt>0):
					self.filtout = "PASS" #High-confidence variant if at least 1 unfiltered, callable callset
				else:
					self.filtout = "allfilteredbutagree" #otherwise, no unfiltered, callable callsets, so uncertain position

				if (self.gt1 == "0/0" and self.agreenofilt>0 and self.disagreenofilt==0 and ((len(self.ref)>10 or len(self.alt)>10) or (self.agreeplatformsnofilt == "Ion" and (len(self.ref)>1 or len(self.alt)>1)))):
					self.filtout = "questionableindel" #even if an unfiltered callset calls homozygous reference, make it uncertain if another callset has a filtered indel >10bp (or any indel if Ion misses it) because we found some callsets miss larger indels
				elif (self.gt1 == "0/0" and self.agreenofilt>0 and self.disagreenofilt==0):
					self.filtout = "." #High-confidence homozygous reference if not a long indel
					self.arbitratetxt = ";arbitrated=TRUE"
				elif (self.agreenofilt>0 and self.disagreenofilt==0 and float(self.gqsum) < self.minGQ):
					self.filtout="GQlessthan70" #if all callable, unfiltered datasets with this call sum to a GQ < 70, then make uncertain because there is insufficient support
				elif (self.agreenofilt > 0 and self.disagreenofilt==0 and float(self.gqsum) >= self.minGQ):
					self.filtout = "PASS" #High-confidence variant if at least 1 unfiltered, callable callset
					self.arbitratetxt = ";arbitrated=TRUE"
				elif (self.agreenofilt==0):
					self.filtout = "allfilteredanddisagree" #no unfiltered, callable callsets, so uncertain position
				else:
					self.filtout = "discordantunfiltered" #otherwise, unfiltered callable callsets had discordant genotype calls
			
			else:
					#Otherwise, all callsets were not callable or had low GQ, so get the first genotype of a callset that is variant
					for j in range(0 , self.callsets):
						if (self.gt1 == ""):
							if (not(self.gt[j] == ".")):
								self.gt1 = self.gt[j]
								self.agree = self.agree + 1
								self.agreedatasets = self.dataset[j]
								self.agreedatasetscnt = self.agreedatasetscnt + 1
								self.agreeplatforms = self.platform[j]
								self.agreeplatformscnt = self.agreeplatformscnt + 1 
								self.gqsum = self.gqsum + float(self.gq[j])
						else:
							if (self.gt[j] == self.gt1):
								self.agree = self.agree + 1
								if (self.gqsum==0 or not(re.search(self.dataset[j], self.agreedatasets))):
									self.gqsum = self.gqsum + float(self.gq[j]) #only add to gqsum once for each dataset
								if not(re.search(self.callset[j], self.agreecallsets)):
									self.agreecallsets = self.agreecallsets + "," + self.callset[j]
									self.agreecallsetscnt = self.agreecallsetscnt + 1
								if not(re.search(self.dataset[j], self.agreedatasets)):
									self.agreedatasets = self.agreedatasets + "," + self.dataset[j]
									self.agreedatasetscnt = self.agreedatasetscnt + 1
								if not(re.search(self.platform[j], self.agreeplatforms)):
									self.agreeplatforms = self.agreeplatforms + "," + self.platform[j]
									self.agreeplatformscnt = self.agreeplatformscnt + 1
							else:
								self.disagree = self.disagree + 1

					if self.disagree > 0:
						self.filtout = "allfilteredanddisagree"
					else:
						self.filtout = "allfilteredbutagree"
			
			#We have now determined if the site is high-confidence, and if not, why not, so sum DP and AD for output annotations
			#sum DP for datasets that agree and count datasets and platforms that disagree and aren't filtered
			self.agreedatasetsnofilt = "" #list of datasets that agree with gt1nofilt and not filtered; reset so that it can be used below to determine if AD's have been added from this dataset yet
			self.agreedatasetsnofiltcnt = 0 #number of datasets that agree with gt1nofilt and not filtered
			self.disagreecallsetsnofilt = ""
			self.disagreedatasetsnofilt = ""
			self.disagreeplatformsnofilt = ""
			self.disagreecallsetsnofiltcnt = 0
			self.disagreedatasetsnofiltcnt = 0
			self.disagreeplatformsnofiltcnt = 0


			for j in range(0, self.callsets):
				if (self.gt[j] == self.gt1 or (self.gt[j] == "." and self.gt1 == "0/0")): #Does this callset agree with the output genotype
					#print("callsets[j]: " + str(j) + " dp[j]: " + str(self.dp[j]))
					if re.search("\d", str(self.dp[j])): 
						self.DPSum = self.DPSum + int(self.dp[j])  #add its DP to DPSum
					if (self.gqsum==0 or not(re.search(str(self.dataset[j]), self.agreedatasets))): 
						self.gqsum = self.gqsum + float(self.gq[j]) #only add to gqsum once for each dataset
					if ((int(self.ADallSum1) + int(self.ADallSum2) + int(self.ADallSum3))==0 or not(re.search(str(self.dataset[j]), self.agreedatasets))): 
						self.ADallSum1 = int(self.ADallSum1)+ (int(self.ad1[j]) if re.search("\d", str(self.ad1[j])) else 0) 
						self.ADallSum2 = int(self.ADallSum2)+ (int(self.ad2[j]) if re.search("\d", str(self.ad2[j])) else 0)
						self.ADallSum3 = int(self.ADallSum3)+ (int(self.ad3[j]) if re.search("\d", str(self.ad3[j])) else 0) #only add to AD sums once for each dataset
					if ((self.callable_field[j]>0 and self.filt[j]==0) and ((int(self.ADSum1) + int(self.ADSum2) + int(self.ADSum3))==0 or not(re.search(str(self.dataset[j]), self.agreedatasetsnofilt)))):
						self.ADSum1 = int(self.ADSum1) + (int(self.ad1[j]) if re.search("\d", str(self.ad1[j])) else 0)  
						self.ADSum2 = int(self.ADSum2) + (int(self.ad2[j]) if re.search("\d", str(self.ad2[j])) else 0)
						self.ADSum3 = int(self.ADSum3) + (int(self.ad3[j]) if re.search("\d", str(self.ad3[j])) else 0)
						if (not(re.search(str(self.dataset[j]), self.agreedatasetsnofilt))): 
								self.agreedatasetsnofilt = self.agreedatasetsnofilt + "," + self.dataset[j]
								self.agreedatasetsnofiltcnt = self.agreedatasetsnofiltcnt + 1
					if (not(re.search(str(self.callset[j]), self.agreecallsets))):
						self.agreecallsets = self.agreecallsets + "," + str(self.callset[j])
						self.agreecallsetscnt = self.agreecallsetscnt + 1
					if (not(re.search(str(self.dataset[j]), self.agreedatasets))):
						self.agreedatasets = self.agreedatasets + "," + str(self.dataset[j])
						self.agreedatasetscnt = self.agreedatasetscnt + 1
					if (not(re.search(str(self.platform[j]), self.agreeplatforms))):
						self.agreeplatforms = self.agreeplatforms + "," + str(self.platform[j])
						self.agreeplatformscnt = self.agreeplatformscnt + 1
				else:  #otherwise, this callset's GT is different from the output genotype
					if (self.agreeplatformsnofiltcnt == 1 and self.filtout == "PASS" and self.gt1 == "0/1" and self.gt[j] == "1/1" and float(self.gq[j])> self.minGQ): 
						self.filtout = "discordanthet" #fix problem with some CG and freebayes sites incorrectly calling hets when they should be homvar
					if (not(re.search(str(self.callset[j]), self.disagreecallsetsnofilt))):
						self.disagreecallsetsnofilt = self.disagreecallsetsnofilt + "," + str(self.callset[j])
						self.disagreecallsetsnofiltcnt = self.disagreecallsetsnofiltcnt + 1
					if (not(re.search(str(self.dataset[j]), self.disagreedatasetsnofilt))):
						self.disagreedatasetsnofilt = self.disagreedatasetsnofilt + "," + str(self.dataset[j])
						self.disagreedatasetsnofiltcnt = self.disagreedatasetsnofiltcnt + 1
					if (not(re.search(str(self.platform[j]), self.disagreeplatformsnofilt))):
						self.disagreeplatformsnofilt = self.disagreeplatformsnofilt + "," + str(self.platform[j])
						self.disagreeplatformsnofiltcnt = self.disagreeplatformsnofiltcnt + 1	
			
			if self.gt1 == "0/1" and self.filtout == "PASS" and (int(self.ADSum1) > (4 * int(self.ADSum2)) or int(self.ADSum2) > (4* int(self.ADSum1))):
				self.filtout = "alleleimbalance" #filter sites that have allele balance < 0.2 or > 0.8
			if self.agreeplatformscnt==1 and (re.search("CG", self.agreeplatforms)) and ("PASS" == self.filtout):
				self.filtout = "cgonly"
			self.gqsum = int(round(self.gqsum)) #round gqsum
			self.ADSumout = str(int(self.ADSum1)) + "," + str(int(self.ADSum2))
			self.ADallSumout = str(int(self.ADallSum1)) + "," + str(int(self.ADallSum2))
			if re.search(",",self.alt):
				self.ADSumout = str(int(self.ADSum1)) + "," + str(int(self.ADSum2)) + "," + str(int(self.ADSum3)) 
				self.ADallSumout = str(int(self.ADallSum1)) + "," + str(int(self.ADallSum2)) + "," + str(int(self.ADallSum3))

			
			#filter sites that are within 50bp of another passing call but none of the callsets that support the 2 calls match
			if (not(self.gt1 == "0/0") and self.filtout == "PASS" and (int(self.pos) - int(self.prevpasspos)) < 50 and (int(self.pos) - int(self.prevpasspos)) >= 0 ):
				self.currentcallsets = self.agreecallsets.split(",")
				self.matchcallset = 0
				for j in range(0,len(self.currentcallsets)):
					if re.search(self.currentcallsets[j], self.agreecallsetsprev):
						self.matchcallset = 1
				if self.matchcallset==0:
					self.filtout="overlappingcall" 
				
			if (not(self.gt1 == "0/0") and self.filtout == "PASS"):
				self.prevpasspos = self.pos
				self.agreecallsetsprev = self.agreecallsets

			
			#output if platform appears to have a bias causing reads with evidence for a false variant
			self.disagreetxt=""
			if self.gt1 == "0/0" and self.filtout == ".":
				self.DPSum = "."
				self.ADSumout = "."
				self.ADallSumout = "."
				self.gqsum = "."
				if re.search(",(..*)", self.disagreeplatformsnofilt):
					self.disagreetxt=";platformbias=" + re.split(",(..*)",self.disagreeplatformsnofilt)[1] #output platforms that have evidence for a variant at a high-confidence hom ref site, since these are potentially caused by bias
				
			elif not(self.gt1 == "0/0") and self.filtout == "PASS":
				if re.search(",(..*)", self.disagreedatasetsnofilt):
					self.disagreetxt = ";datasetsmissingcall=" + re.split(",(..*)", self.disagreedatasetsnofilt)[1] #output platforms that have an incorrect call at a high-confidence variant site
			elif self.filtout == "discordantunfiltered": #output if there is a callset that uniquely calls a PASSing genotype when 2+ PASSing callsets support a different genotype
				if self.agreecallsetsnofiltcnt > 1 and self.disagreecallsetsnofiltcnt == 1 and re.search(",(..*)", self.disagreecallsetsnofilt):
					self.disagreetxt = ";callsetwithotheruniqgenopassing=" + re.split(",(..*)", self.disagreecallsetsnofilt)[1]
				elif self.disagreecallsetsnofiltcnt > 1 and self.agreecallsetsnofiltcnt == 1 and re.search(",(..*)",self.agreecallsetsnofilt):
					self.disagreetxt = ";callsetwiththisuniqgenopassing=" + re.split(",(..*)", self.agreecallsetsnofilt)[1]

			
			self.qualout=5 #default QUAL to 5 if filtered
			if self.filtout == "PASS":
				self.qualout=50 #50 if passing variant
			elif self.filtout == ".":
				self.qualout=0 #0 if passing homref site
			if not(self.gtphased == "") and self.gt1 == self.gt[0]:
				self.gt1 = self.gtphased #keep phasing info from first callset

			self.output_allcalls.write(str(self.chrom) + "\t" + str(self.pos) + "\t" + str(self.id) + "\t" + str(self.ref) + "\t" + str(self.alt) + "\t" + str(self.qualout) + "\t" + str(self.filtout) + "\t" + "platforms=" + str(self.agreeplatformscnt) + ";platformnames=" + str(self.agreeplatforms) + ";datasets=" + str(self.agreedatasetscnt) + ";datasetnames=" + str(self.agreedatasets) + ";callsets=" + str(self.agreecallsetscnt) + ";callsetnames=" + str(self.agreecallsets) + str(self.disagreetxt) + str(self.callabletxt) + str(self.filttxt) + str(self.arbitratetxt) + str(self.difficultregion) + "\t" + "GT:PS:DP:ADALL:AD:GQ" + "\t" + str(self.gt1) + ":" + str(self.ps) + ":" + str(self.DPSum) + ":" + str(self.ADallSumout) + ":" + str(self.ADSumout) + ":" + str(self.gqsum) + "\n")
			if self.filtout == "PASS":
				self.output_arbitrated_file.write(str(self.chrom) + "\t" + str(self.pos) + "\t" + str(self.id) +"\t" + str(self.ref) + "\t" + str(self.alt) + "\t" + str(self.qualout) + "\t" + str(self.filtout) + "\t" + "platforms=" + str(self.agreeplatformscnt) + ";platformnames=" + str(self.agreeplatforms) + ";datasets=" + str(self.agreedatasetscnt) + ";datasetnames=" + str(self.agreedatasets) + ";callsets=" + str(self.agreecallsetscnt) + ";callsetnames=" + str(self.agreecallsets) + str(self.disagreetxt) + str(self.callabletxt) + str(self.filttxt) + str(self.arbitratetxt) + str(self.difficultregion) + "\t" + "GT:PS:DP:ADALL:AD:GQ" + "\t" + str(self.gt1) + ":" + str(self.ps) + ":" + str(self.DPSum) + ":" + str(self.ADallSumout) + ":" + str(self.ADSumout) + ":" + str(self.gqsum) + "\n")
			if self.filtout == "PASS" and self.agreeplatformsnofiltcnt > 1:
				self.output_2_platforms_file.write(str(self.chrom) + "\t" + str(self.pos) + "\t" + str(self.id) +"\t" + str(self.ref) + "\t" + str(self.alt) + "\t" + str(self.qualout) + "\t" + str(self.filtout) + "\t" + "platforms=" + str(self.agreeplatformscnt) + ";platformnames=" + str(self.agreeplatforms) + ";datasets=" + str(self.agreedatasetscnt) + ";datasetnames=" + str(self.agreedatasets) + ";callsets=" + str(self.agreecallsetscnt) + ";callsetnames=" + str(self.agreecallsets) + str(self.disagreetxt) + str(self.callabletxt) + str(self.filttxt) + str(self.arbitratetxt) + str(self.difficultregion) + "\t" + "GT:PS:DP:ADALL:AD:GQ" + "\t" + str(self.gt1) + ":" + str(self.ps) + ":" + str(self.DPSum) + ":" + str(self.ADallSumout) + ":" + str(self.ADSumout) + ":" + str(self.gqsum) + "\n")
			if self.agreeplatformsnofilt == "CG" and self.gt1 == "1/1":
				self.output_test_file.write(str(self.DPSum) + "\t" + "gqsum\n")				
		

	def run_filter_classification(self):
		self.parse_input_files()
		self.classify_and_filter()
		self.vcfall_file.close()
		self.callsettable_file.close()
		self.output_allcalls.close()
		self.output_arbitrated_file.close()
		self.output_2_platforms_file.close()
		self.output_test_file.close()