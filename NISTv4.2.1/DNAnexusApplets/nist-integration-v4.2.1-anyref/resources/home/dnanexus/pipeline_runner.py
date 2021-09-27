import os
import gzip
import shutil
import glob
import subprocess
import sys
import preprocess_combine_vcfs
import run_one_class_filter
import vcf_classify_using_filters
import vcf_one_class_filtering
import dxpy
import re


##Move all files from file arrays to the same data folder
class Integrator(object):
  file_list = None
  gz_file =  None
  unzipped_file = None
  genome_file = None
  genome_file_lines = None
  human_genome_file = None
  human_genome_bed_file = None
  nohomref_samplename_file_list = None
  vcf_list = None
  vars_file_list = None
  count = None
  tandemrepeatsbed = None
  techs_with_no_callable = None

  def __init__(self, callset_path_in, chrom_in, refn_in, prefix_in, tandemrepeatsbed_in, techs_with_no_callable_in):
    self.file_list = None
    self.gz_file =  None
    self.unzipped_file = None
    self.genome_file = None
    self.genome_file_lines = None
    self.human_genome_file = None
    self.human_genome_bed_file = None
    self.nohomref_samplename_file_list = None
    self.vcf_list = None
    self.vars_file_list = None
    self.count = None
    self.techs_with_no_callable = techs_with_no_callable_in
    print(callset_path_in + "\n")
    self.callset_table_path = "/home/dnanexus/in/callsettable/" + os.listdir("/home/dnanexus/in/callsettable")[0]
    print(chrom_in + "\n")
    self.chrom = chrom_in
    print(refn_in + "\n")
    self.refn = "/home/dnanexus/in/refn/" + os.listdir("/home/dnanexus/in/refn")[0]
    print(prefix_in + "\n")
    self.prefix = prefix_in
    print(str(tandemrepeatsbed_in) + "\n")
    self.tandemrepeatsbed = "/home/dnanexus/in/tandemrepeatsbed/" + os.listdir("/home/dnanexus/in/tandemrepeatsbed")[0]
    print(str(self.tandemrepeatsbed) + "\n")
    print(techs_with_no_callable_in + "\n")
    print(self.techs_with_no_callable + "\n")


  def import_files(self): 
    subprocess.check_call("mkdir data", shell=True)
    subprocess.check_call("mv ~/in/vcfs/*/* data", shell=True)
    subprocess.check_call("mv ~/in/beds/*/* data", shell=True)
    subprocess.check_call("mv ~/in/annotations/*/* data", shell=True)
    subprocess.check_call("mv ~/in/filtbeds/*/* .", shell=True)
    subprocess.check_call("mv ~/in/ref/* .", shell=True)
    subprocess.check_call("mv ~/in/rtgsdf/* .", shell=True)
    #unzip ref and rtgsdf files
    for v in glob.glob("*.tar.gz"):
        subprocess.check_call("tar -zxvf " + v, shell=True)


  #make genome and bed files for bedtools
  def preprocess(self): 
    self.genome_file = open("genome.dict", "r")
    self.genome_file_lines = self.genome_file.readlines()
    self.human_genome_file = open("human.genome", "w+")
    for i in range(1,26):
      split_line = self.genome_file_lines[i].split("\t")
      sn_split = split_line[1].split(":")[1]
      ln_split = split_line[2].split(":")[1]
      self.human_genome_file.write(sn_split + "\t" + ln_split + "\n")

    self.human_genome_file.flush()
    self.human_genome_file.close()
    self.genome_file.close()

    self.genome_file = open("genome.dict", "r")
    self.genome_file_lines = self.genome_file.readlines()
    self.human_genome_bed_file = open("human.genome.bed", "w+")
    for i in range(1,26):
      split_line = self.genome_file_lines[i].split("\t")
      sn_split = split_line[1].split(":")[1]
      ln_split = split_line[2].split(":")[1]
      self.human_genome_bed_file.write(sn_split + "\t" + "1" + "\t" + ln_split + "\n")

    self.human_genome_bed_file.flush()
    self.human_genome_bed_file.close()
    self.genome_file.close()

    #Process vcf and bed files (adding sample names, removing chr, removing homref sites)
    ppc_obj = preprocess_combine_vcfs.PreProcessCombine(cstable=self.callset_table_path, p = "data/")
    ppc_obj.preprocess()
			
    self.nohomref_samplename_file_list = glob.glob("data/*_nohomref_samplename.vcf.gz")
    self.nohomref_samplename_file_list.sort()
    self.vcf_list_str = " "
    for n in self.nohomref_samplename_file_list:
        self.vcf_list_str = self.vcf_list_str + n + " " 
    subprocess.check_call("vcfcombine " + str(self.vcf_list_str) + " > union.1.vcf", shell=True)

    i = 1
    callable_processed_file_list = glob.glob("data/*_callable_processed.bed")
    callable_processed_file_list.sort()
    for v in callable_processed_file_list:
      j = i + 1
      subprocess.check_call("vcfannotate -b " + v + " -k callable union." + str(i) + ".vcf > union." + str(j) + ".vcf", shell=True)
      subprocess.check_call("awk \'{ sum+=$3; sum-=$2 } END { print sum }\' " + v + "", shell=True)
      i = i + 1

    ##Annotate union vcf with difficult region bed files in the "difficultregion" INFO field (currently only used for information purposes)
    diff_bed_annot_file_list = glob.glob("*_diffbedannotcol.bed")
    diff_bed_annot_file_list.sort()
    for v in diff_bed_annot_file_list:
      j = i + 1
      subprocess.check_call("vcfannotate -b " + v + " -k difficultregion union." + str(i) + ".vcf > union." + str(j) + ".vcf", shell=True)
      subprocess.check_call("awk \'{ sum+=$3; sum-=$2 } END { print sum }\' " + v + "", shell=True)
      i = i + 1

    subprocess.check_call("cp union." + str(j) + ".vcf union_callableannotate.vcf", shell=True)

  def classify_and_intersect(self):
    ##Run first pass of integration to find training variants found in 2 technologies
    cuf_obj = vcf_classify_using_filters.ClassifyUsingFilter("union_callableannotate.vcf", self.callset_table_path, "union_callableannotate")
    cuf_obj.run_filter_classification()

    #bgzip and index integrated vcfs
    subprocess.check_call("bgzip union_callableannotate_ClassifyUsingFilters_2platforms.vcf", shell = True)
    subprocess.check_call("tabix -p vcf union_callableannotate_ClassifyUsingFilters_2platforms.vcf.gz", shell=True)

    subprocess.check_call("bgzip union_callableannotate_ClassifyUsingFilters_allcalls.vcf", shell=True)
    subprocess.check_call("tabix -p vcf union_callableannotate_ClassifyUsingFilters_allcalls.vcf.gz",shell=True)
    
    nohomref_samplename_file_list = glob.glob("data/*_nohomref_samplename.vcf.gz")
    nohomref_samplename_file_list.sort()
    for v in nohomref_samplename_file_list:
      subprocess.check_call("vcfintersect --intersect-vcf union_callableannotate_ClassifyUsingFilters_2platforms.vcf.gz -r genome.fa " + v + " > " + v + "_indivintersect2platforms.vcf", shell=True)
      subprocess.check_call("wc -l " + v + "_indivintersect2platforms.vcf", shell =True)
    

  def run_classifier(self): 
    ##Run perl script that runs another perl script that performs one-class filtering for each callset and generates filtering bed files for each callset (also generates filtered vcf files for vgraph method)
    roncf_obj = run_one_class_filter.RunOneClassFilter(self.callset_table_path, "data/")
    roncf_obj.run_classifier()

  def classify_final(self): 
    ##Annotate union vcf with each of the filtering bed files in the "oneclassfilt" INFO field
    subprocess.check_call("cp union_callableannotate.vcf union_callableannotate.1.vcf", shell=True)

    i=1
    samplename_oneclassfilter_file_list = glob.glob("data/*_samplename_oneclassfilter_filtered.bed")
    samplename_oneclassfilter_file_list.sort()
    for v in samplename_oneclassfilter_file_list:
      j = i + 1
      subprocess.check_call("vcfannotate -b " + v + " -k oneclassfilt union_callableannotate." + str(i) + ".vcf > union_callableannotate." + str(j) + ".vcf", shell=True)
      i = i + 1

    subprocess.check_call("cp union_callableannotate." + str(j) + ".vcf union_callableannotate_filterannotate.vcf", shell=True)

    ##Run final pass of integration to find high-confidence variants using not-callable bed and filtering bed for arbitartion
    cuf_obj2 = vcf_classify_using_filters.ClassifyUsingFilter("union_callableannotate_filterannotate.vcf", self.callset_table_path, "union_callableannotate_filterannotate")    
    cuf_obj2.run_filter_classification()

  def prep_bed_and_summarize(self):
    subprocess.check_call("bgzip union_callableannotate_filterannotate.vcf", shell=True)
    subprocess.check_call("tabix -p vcf union_callableannotate_filterannotate.vcf.gz", shell=True)

    subprocess.check_call("bgzip union_callableannotate_filterannotate_ClassifyUsingFilters_allcalls.vcf", shell=True)
    subprocess.check_call("tabix -p vcf union_callableannotate_filterannotate_ClassifyUsingFilters_allcalls.vcf.gz", shell=True)

    subprocess.check_call("bgzip union_callableannotate_filterannotate_ClassifyUsingFilters_arbitrated.vcf", shell=True)
    subprocess.check_call("tabix -p vcf union_callableannotate_filterannotate_ClassifyUsingFilters_arbitrated.vcf.gz", shell=True)

    subprocess.check_call("bgzip union_callableannotate_filterannotate_ClassifyUsingFilters_2platforms.vcf", shell=True)
    subprocess.check_call("tabix -p vcf union_callableannotate_filterannotate_ClassifyUsingFilters_2platforms.vcf.gz", shell=True)

    ##Calculate some stats about integrated variants
    subprocess.check_call("rtg vcfstats union_callableannotate_filterannotate_ClassifyUsingFilters_allcalls.vcf.gz", shell = True)

    subprocess.check_call("rtg vcfstats union_callableannotate_filterannotate_ClassifyUsingFilters_2platforms.vcf.gz", shell=True)

    subprocess.check_call("rtg vcfstats union_callableannotate_filterannotate_ClassifyUsingFilters_arbitrated.vcf.gz", shell=True)

    subprocess.check_call("zgrep GQlessthan70 union_callableannotate_filterannotate_ClassifyUsingFilters_allcalls.vcf.gz | wc -l", shell=True) 
    subprocess.check_call("zgrep allfilteredanddisagree union_callableannotate_filterannotate_ClassifyUsingFilters_allcalls.vcf.gz | wc -l", shell=True)
    subprocess.check_call("zgrep allfilteredbutagree union_callableannotate_filterannotate_ClassifyUsingFilters_allcalls.vcf.gz | wc -l", shell=True)
    subprocess.check_call("zgrep discordantunfiltered union_callableannotate_filterannotate_ClassifyUsingFilters_allcalls.vcf.gz | wc -l", shell=True)
    subprocess.check_call("zgrep discordanthet union_callableannotate_filterannotate_ClassifyUsingFilters_allcalls.vcf.gz | wc -l", shell=True)
    subprocess.check_call("zgrep questionableindel union_callableannotate_filterannotate_ClassifyUsingFilters_allcalls.vcf.gz | wc -l", shell=True)
    subprocess.check_call("zgrep cgonly union_callableannotate_filterannotate_ClassifyUsingFilters_allcalls.vcf.gz | wc -l", shell=True)
    subprocess.check_call("zgrep alleleimbalance union_callableannotate_filterannotate_ClassifyUsingFilters_allcalls.vcf.gz | wc -l", shell=True)

    ##Create high-confidence bed file
    ##first, find how many datasets cover each region
    # Selecting chrom
    self.g_file = open("human.genome", "r")
    self.g_file_lines = self.g_file.readlines()
    self.g_chrom_file = open("human.chrom.genome", "w+")
    for g in self.g_file_lines:
      split_line = g.split("\t")
      if split_line[0] == str(self.chrom):
        self.g_chrom_file.write(g)
        self.g_chrom_file.flush()
        break
    self.g_file.close()
    self.g_chrom_file.close()

    self.g_bed_file = open("human.genome.bed", "r")
    self.g_bed_file_lines = self.g_bed_file.readlines()
    self.g_bed_chrom_file = open("human.chrom.genome.bed", "w+")
    for gb in self.g_bed_file_lines:
      split_line = gb.split("\t")
      if split_line[0] == str(self.chrom):
        self.g_bed_chrom_file.write(gb)
        self.g_bed_chrom_file.flush()
        break
    self.g_bed_file.close()
    self.g_bed_chrom_file.close()


    self.vars = ""
    self.vars_file_list = glob.glob("data/*_callable_processed.bed")
    for v in self.vars_file_list:
      self.count = len(open(v, "r").readlines())
      if self.count > 1:
        self.vars = self.vars + " " + v

    subprocess.check_call("multiIntersectBed -i " + str(self.vars) + " | awk -v chrom=" + str(self.chrom) + " \'BEGIN {FS = OFS = \"\t\"} {if($1 == chrom || chrom == \"allbutY\" || chrom == \"all\" || $1 ~ /^#/) print }\' > callablemultinter.bed", shell=True)

    ##Create bed files with regions that are callable in at least 1 or at least 2 datasets
    subprocess.check_call("awk \'{if($4>1) print $0}\' callablemultinter.bed | mergeBed -i stdin > callablemultinter_gt1.bed", shell=True)
    subprocess.check_call("awk \'{if($4>0) print $0}\' callablemultinter.bed | mergeBed -i stdin > callablemultinter_gt0.bed", shell=True)
    subprocess.check_call("awk \'{ sum+=$3; sum-=$2 } END { print sum }\' callablemultinter_gt1.bed", shell=True)
    subprocess.check_call("awk \'{ sum+=$3; sum-=$2 } END { print sum }\' callablemultinter_gt0.bed", shell=True)

    #Create bed file with every potential variant that is not high-confidence +-50bp
    subprocess.check_call("gunzip -c union_callableannotate_filterannotate_ClassifyUsingFilters_allcalls.vcf.gz | awk \'{OFS=\"\t\"; if (!($7 ~ /^PASS/ || $7 ~ /^\./ || $1 ~ /^#/)) {print $1,$2-1,$2 + length($4),$4\"/\"$5}}\' | subtractBed -a human.chrom.genome.bed -b stdin | complementBed -i stdin -g human.chrom.genome | slopBed -i stdin -g human.chrom.genome -b 50 > union_callableannotate_filterannotate_ClassifyUsingFilters_filtered.bed", shell=True)

    subprocess.check_call("rtg vcfeval -b union_callableannotate_filterannotate_ClassifyUsingFilters_arbitrated.vcf.gz -c union_callableannotate_filterannotate_ClassifyUsingFilters_arbitrated.vcf.gz -o selfcompare -t rtgsdf", shell=True)

    subprocess.check_call("zgrep -v ^# selfcompare/fp.vcf.gz | awk \'BEGIN {FS = OFS = \"\t\"} {print $1,$2-50,$2+length($4)+50}\' | sed \'s/^X/23/;s/^Y/24/;s/^M/25/\' | sort -k1,1n -k2,2n | sed \'s/^23/X/;s/^24/Y/;s/^25/MT/\' | mergeBed -i stdin -d 50 > excludeoverlappingvars.bed", shell=True)

    subprocess.check_call("subtractBed -a callablemultinter_gt0.bed -b union_callableannotate_filterannotate_ClassifyUsingFilters_filtered.bed | subtractBed -a stdin -b excludeoverlappingvars.bed | subtractBed -a stdin -b " + str(self.refn) + "  > callablemultinter_gt0_nofilt_nooverlapvar.bed", shell=True)
    subprocess.check_call("subtractBed -a callablemultinter_gt1.bed -b union_callableannotate_filterannotate_ClassifyUsingFilters_filtered.bed | subtractBed -a stdin -b excludeoverlappingvars.bed | subtractBed -a stdin -b " + str(self.refn) + "  > callablemultinter_gt1_nofilt_nooverlapvar.bed", shell=True)
 
 
    subprocess.check_call("awk \'{ sum+=$3; sum-=$2 } END { print sum }\' callablemultinter_gt1_nofilt_nooverlapvar.bed", shell=True)
    subprocess.check_call("awk \'{ sum+=$3; sum-=$2 } END { print sum }\' callablemultinter_gt0_nofilt_nooverlapvar.bed", shell=True)

    ## Exclude fns if the call set is listed as making the call to fix vcfcombine issue
    subprocess.check_call("cp callablemultinter_gt0_nofilt_nooverlapvar.bed temp_highconf_0.bed", shell=True)
    vcf_file_name_without_prefix = ""
    callset_name = ""
    callable_bed_name_suffix = "_callable_processed.bed"
    filtered_bed_name_suffix = "_oneclassfilter_filtered.bed"
    previous_counter = 0
    current_counter = 1
    subprocess.check_call("ls data/", shell=True)
    vcfs_file_list = ""
    vcfs_file_list = glob.glob("data/*_nohomref_samplename.vcf.gz")
    for v in vcfs_file_list:
      if re.search(self.techs_with_no_callable, v) != None:
        continue
      delim1 = re.compile("VCF[0-9][0-9]_")
      vcf_file_name_without_prefix = delim1.split(v)[1]
      delim2 = "_nohomref_samplename.vcf.gz"
      callset_name = vcf_file_name_without_prefix.split(delim2)[0] 
      callable_bed_name = callset_name + callable_bed_name_suffix
      delim3 = ".vcf.gz"
      filter_bed_name_start = v.split(delim3)[0]
      filtered_bed_name = filter_bed_name_start + filtered_bed_name_suffix
      subprocess.check_call("rtg vcfeval -b selfcompare/tp.vcf.gz -c " + v + " -o " + callset_name + "_vs_highconf -t rtgsdf", shell=True)
      subprocess.check_call("bgzip -d " + callset_name + "_vs_highconf/fn.vcf.gz", shell=True)
      subprocess.check_call("bedtools intersect -header -a " + callset_name + "_vs_highconf/fn.vcf -b data/" + callable_bed_name + " > temp_FN_intersect_callable_" + str(current_counter) + ".vcf", shell=True)
      subprocess.check_call("bedtools subtract -a temp_FN_intersect_callable_" + str(current_counter) + ".vcf -b " + filtered_bed_name + " > temp_FN_intersect_callable_subtract_filtered_" + str(current_counter) + ".vcf", shell=True)
      subprocess.check_call("awk \'{FS=OFS=\"\t\"} {print $1,$2-50,$2+length($4)+50}\' temp_FN_intersect_callable_subtract_filtered_" + str(current_counter) + ".vcf | bedtools subtract -a temp_highconf_" + str(previous_counter) + ".bed -b stdin > temp_highconf_" + str(current_counter) + ".bed" , shell=True)
      current_counter = current_counter + 1
      previous_counter = previous_counter + 1
    
    # As a result of loop and increments, previous_counter will hold the value for temp_highconf
    subprocess.check_call("cp temp_highconf_" + str(previous_counter) + ".bed callablemultinter_gt0_nofilt_nooverlapvar.bed", shell=True)

    ## subtract partial tandem repeats
    subprocess.check_call("complementBed -i callablemultinter_gt0_nofilt_nooverlapvar.bed -g human.genome | intersectBed -wa -a " + str(self.tandemrepeatsbed) + " -b stdin | subtractBed -a callablemultinter_gt0_nofilt_nooverlapvar.bed -b stdin > callablemultinter_gt0_nofilt_nooverlapvar_removepartialrepeats.bed", shell=True)

    ##Calculate some stats about integrated variants

    subprocess.check_call("rtg vcffilter -i selfcompare/tp.vcf.gz --include-bed callablemultinter_gt0_nofilt_nooverlapvar_removepartialrepeats.bed -o selfcompare/tp_inhighconfbed.vcf.gz", shell=True)

    subprocess.check_call("rtg vcfstats selfcompare/tp_inhighconfbed.vcf.gz", shell=True)

    subprocess.check_call("rtg vcffilter -i selfcompare/tp.vcf.gz --exclude-bed callablemultinter_gt0.bed -o selfcompare/tp_notinuniongt0bed.vcf.gz", shell=True)

    subprocess.check_call("rtg vcffilter -i selfcompare/tp.vcf.gz --include-bed union_callableannotate_filterannotate_ClassifyUsingFilters_filtered.bed -o selfcompare/tp_infiltslop50bed.vcf.gz", shell=True)

    subprocess.check_call("rtg vcffilter -i selfcompare/tp.vcf.gz --include-bed excludeoverlappingvars.bed -o selfcompare/tp_inexcludeoverlappingvarsbed.vcf.gz", shell=True)

    subprocess.check_call("tar -zcf indivfilteredbeds.tar.gz data/*_callable_processed.bed", shell=True)

  #
  # Upload results
  #
  def upload_data(self):
    v = 1

    output = {}

    subprocess.check_call("cp union_callableannotate_filterannotate.vcf.gz " + self.prefix + "_annotated.vcf.gz", shell=True)
    annotated_vcf_output_file = dxpy.upload_local_file(self.prefix + "_annotated.vcf.gz")
    output["vcfanngz"] = dxpy.dxlink(annotated_vcf_output_file)

    subprocess.check_call("cp union_callableannotate_filterannotate.vcf.gz.tbi " + self.prefix + "_annotated.vcf.gz.tbi", shell=True)
    annotated_vcf_tbi_output_file = dxpy.upload_local_file(self.prefix + "_annotated.vcf.gz.tbi")
    output["vcfanntbi"] = dxpy.dxlink(annotated_vcf_tbi_output_file)

    subprocess.check_call("cp union_callableannotate_filterannotate_ClassifyUsingFilters_allcalls.vcf.gz " + self.prefix + "_all.vcf.gz", shell=True)
    all_vcf_output_file = dxpy.upload_local_file(self.prefix + "_all.vcf.gz")
    output["vcfallgz"] = dxpy.dxlink(all_vcf_output_file)

    subprocess.check_call("cp union_callableannotate_filterannotate_ClassifyUsingFilters_allcalls.vcf.gz.tbi " + self.prefix + "_all.vcf.gz.tbi", shell=True)
    all_vcf_tbi_output_file = dxpy.upload_local_file(self.prefix + "_all.vcf.gz.tbi")
    output["vcfalltbi"] = dxpy.dxlink(all_vcf_tbi_output_file)

    subprocess.check_call("cp selfcompare/tp.vcf.gz " + self.prefix + "_highconf.vcf.gz", shell=True)
    high_conf_vcf_output_file = dxpy.upload_local_file(self.prefix + "_highconf.vcf.gz")
    output["vcfhighconfgz"] = dxpy.dxlink(high_conf_vcf_output_file)    

    subprocess.check_call("cp selfcompare/tp.vcf.gz.tbi " + self.prefix + "_highconf.vcf.gz.tbi", shell=True)
    high_conf_vcf_tbi_output_file = dxpy.upload_local_file(self.prefix + "_highconf.vcf.gz.tbi")
    output["vcfhighconftbi"] = dxpy.dxlink(high_conf_vcf_tbi_output_file)    

    subprocess.check_call("cp union_callableannotate_filterannotate_ClassifyUsingFilters_2platforms.vcf.gz " + self.prefix + "_highconf_2tech.vcf.gz", shell=True)
    high_conf_2tech_vcf_output_file = dxpy.upload_local_file(self.prefix + "_highconf_2tech.vcf.gz")
    output["vcfhighconf2techgz"] = dxpy.dxlink(high_conf_2tech_vcf_output_file)    

    subprocess.check_call("cp union_callableannotate_filterannotate_ClassifyUsingFilters_2platforms.vcf.gz.tbi " + self.prefix + "_highconf_2tech.vcf.gz.tbi", shell=True)
    high_conf_2tech_vcf_tbi_output_file = dxpy.upload_local_file(self.prefix + "_highconf_2tech.vcf.gz.tbi")
    output["vcfhighconf2techtbi"] = dxpy.dxlink(high_conf_2tech_vcf_tbi_output_file)    

    subprocess.check_call("cp union_callableannotate_ClassifyUsingFilters_allcalls.vcf.gz " + self.prefix + "_nofilt_all.vcf.gz", shell=True)
    all_nofilt_vcf_output_file = dxpy.upload_local_file(self.prefix + "_nofilt_all.vcf.gz")
    output["vcfallnofiltgz"] = dxpy.dxlink(all_nofilt_vcf_output_file)   

    subprocess.check_call("cp union_callableannotate_ClassifyUsingFilters_allcalls.vcf.gz.tbi " + self.prefix + "_nofilt_all.vcf.gz.tbi", shell=True)
    all_nofilt_vcf_tbi_output_file = dxpy.upload_local_file(self.prefix + "_nofilt_all.vcf.gz.tbi")
    output["vcfallnofilttbi"] = dxpy.dxlink(all_nofilt_vcf_tbi_output_file)       

    subprocess.check_call("cp callablemultinter_gt0_nofilt_nooverlapvar_removepartialrepeats.bed " + self.prefix + "_highconf.bed", shell=True)
    bed1_output_file = dxpy.upload_local_file(self.prefix + "_highconf.bed")
    output["bed1"] = dxpy.dxlink(bed1_output_file)      

    subprocess.check_call("cp callablemultinter_gt1_nofilt_nooverlapvar.bed " + self.prefix + "_highconf_gt1.bed", shell=True)
    bed2_output_file = dxpy.upload_local_file(self.prefix + "_highconf_gt1.bed")
    output["bed2"] = dxpy.dxlink(bed2_output_file)          

    subprocess.check_call("cp callablemultinter_gt0.bed " + self.prefix + "_callablemultinter_gt0.bed", shell=True)
    bedcallablegt0_output_file = dxpy.upload_local_file(self.prefix + "_callablemultinter_gt0.bed")
    output["bedcallablegt0"] = dxpy.dxlink(bedcallablegt0_output_file)          

    subprocess.check_call("cp callablemultinter.bed " + self.prefix + "_callablemultinter_all.bed", shell=True)
    bedcallableall_output_file = dxpy.upload_local_file(self.prefix + "_callablemultinter_all.bed")
    output["bedcallableall"] = dxpy.dxlink(bedcallableall_output_file)

    subprocess.check_call("cp union_callableannotate_filterannotate_ClassifyUsingFilters_filtered.bed " + self.prefix + "_filteredsites.bed", shell=True)
    bedfilteredsites_output_file = dxpy.upload_local_file(self.prefix + "_filteredsites.bed")
    output["bedfilteredsites"] = dxpy.dxlink(bedfilteredsites_output_file)

    subprocess.check_call("cp indivfilteredbeds.tar.gz " + self.prefix + "_indivfilteredbeds.tar.gz", shell=True)
    bedindivfilteredsites_output_file = dxpy.upload_local_file(self.prefix + "_indivfilteredbeds.tar.gz")
    output["bedindivfilteredsites"] = dxpy.dxlink(bedindivfilteredsites_output_file)

    return output    

  def run_pipeline(self):
    import_files()
    preprocess()
    classify_and_intersect()
    run_classifier()
    classify_final()
    prep_bed_and_summarize()
    upload_data()
