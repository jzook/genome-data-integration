import argparse
import subprocess
import re
import sys
import vcf_one_class_filtering


class RunOneClassFilter(object):
  parser = None
  args = None

  path = None
  cstable_file = None

  cstable_file_lines = None
  cstable_file_header = None
  output_file_lines = None
  m = None
  vcfs = None
  cstable_line_head = []
  vcf_no_string = None
  line = None
  line_data = None
  vcf_prefix = None

  def __init__(self, cstable, p):
    # self.parser = argparse.ArgumentParser(description='Process and combine VCF files.', usage = "perl RunOneClassFilter.pl -cstable <txt> -p <FILEPATH>")

    # self.parser.add_argument('-cstable', metavar='txt', help='The input Tab-delimited text file containing columns Platform, Dataset, Callset, vcfAll.vcf.gz, callableBed.bed, and annotationsFile')
    # self.parser.add_argument('-p', metavar='FILEPATH', help='Path where vcf and bed files are located')

    # self.args = self.parser.parse_args()
    # if self.args.p is None:
    #   ##print("Please specify the path of the vcf and bed files.\n\n")
    #   exit()
    # self.path = self.args.p

    if p is None:
      ##print("Please specify the path of the vcf and bed files.\n\n")
      exit()
    self.path = p

    # try:
    #   self.cstable_file = open(self.args.cstable, "r")
    # except:
    #   ##print("Can't open the input table file. Please check if the file exists: " + args.cstable)
    #   exit()
    try:
      self.cstable_file = open(cstable, "r")
    except:
      ##print("Can't open the input table file. Please check if the file exists: " + cstable)
      exit()
      
    self.cstable_file_lines = self.cstable_file.readlines()
    self.cstable_file_header = ""
    self.output_file_lines = []
    self.m = 1
    self.vcfs = ""
    self.cstable_line_head = []
    self.vcf_no_string = None
    self.line = None
    self.line_data = None
    self.vcf_prefix = ""

  def process_callset_table(self): 
    ##print(self.cstable_file_lines)
    for li in self.cstable_file_lines:
            
      self.vcf_no_string = str(self.m)
      if self.m < 10:
        self.vcf_no_string = "0" + str(self.m)
          
      self.line = li.rstrip("\n")
      self.line_data = self.line.split("\t")
    
      if '#' in self.line:
        continue

      self.vcf_prefix = ""
      if re.search("(.*)\.vcf.gz", self.line_data[3]):
        self.vcf_prefix = re.split("(.*)\.vcf.gz", self.line_data[3])[1]
      else:
        ##print("vcfAll does not end in .vcf.gz\n")
        exit()
      if "none" == self.line_data[5]:
        subprocess.check_call("cp " + self.path + "VCF" + self.vcf_no_string + "_" + self.vcf_prefix + "_nohomref_samplename.vcf " + self.path + "VCF" + self.vcf_no_string + "_" + self.vcf_prefix + "_nohomref_samplename_oneclassfilter_filtered.vcf", shell=True)      
      else:
        #subprocess.check_call("perl VcfOneClassFiltering_v3.3.pl " + self.path + "VCF" + self.vcf_no_string + "_" + self.vcf_prefix + "_nohomref_samplename.vcf.gz_indivintersect2platforms.vcf " + self.path + "VCF" + self.vcf_no_string + "_" + self.vcf_prefix + "_nohomref_samplename.vcf " + self.path + self.line_data[5] + " " + self.path + "VCF" + self.vcf_no_string + "_" + self.vcf_prefix + "_nohomref_samplename_oneclassfilter " + self.line_data[2] + "\n", shell=True)
        print(self.path + "VCF" + self.vcf_no_string + "_" + self.vcf_prefix + "_nohomref_samplename.vcf.gz_indivintersect2platforms.vcf")
        print(self.path + "VCF" + self.vcf_no_string + "_" + self.vcf_prefix + "_nohomref_samplename.vcf")
        print(self.path + self.line_data[5])
        print(self.path + "VCF" + self.vcf_no_string + "_" + self.vcf_prefix + "_nohomref_samplename_oneclassfilter")
        print(self.line_data[2])
        ocf_obj = vcf_one_class_filtering.OneClassFiltering(self.path + "VCF" + self.vcf_no_string + "_" + self.vcf_prefix + "_nohomref_samplename.vcf.gz_indivintersect2platforms.vcf", self.path + "VCF" + self.vcf_no_string + "_" + self.vcf_prefix + "_nohomref_samplename.vcf", self.path + self.line_data[5], self.path + "VCF" + self.vcf_no_string + "_" + self.vcf_prefix + "_nohomref_samplename_oneclassfilter", self.line_data[2])
        ocf_obj.one_class_filtering()
        ##print("perl VcfOneClassFiltering_v3.3.pl " + self.path + "VCF" + self.vcf_no_string + "_" + self.vcf_prefix + "_nohomref_samplename.vcf.gz_indivintersect2platforms.vcf " + self.path + "VCF" + self.vcf_no_string + "_" + self.vcf_prefix + "_nohomref_samplename.vcf " + self.path + self.line_data[5] + " " + self.path + "VCF" + self.vcf_no_string + "_" + self.vcf_prefix + "_nohomref_samplename_oneclassfilter " + self.line_data[2] + "\n")
      self.m = self.m + 1    

    self.cstable_file.close()

  def run_classifier(self):
    self.process_callset_table()
