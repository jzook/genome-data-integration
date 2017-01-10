#!/usr/bin/env python

import os
import sys
import getopt
import bz2, gzip    


def main():

    
    # parse command line options and update variables
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:o:s:", [])
    except getopt.GetoptError:
        sys.exit(usage)

    for opt, arg in opts:
        if opt == "-h":
            sys.exit('Use me properly.')
        elif opt == '-i':
            INFILE = arg
        elif opt == '-o':
            OUTFILE = arg
            
    outf = open(OUTFILE, 'w')
    if "gz" in INFILE:
        inf = gzip.GzipFile(INFILE, 'r')
    elif "bz2" in INFILE:
        inf = bz2.BZ2File(INFILE, 'r')
    else:
        inf = open(INFILE, 'r')

    for line in inf:

        # Don't change header
        if line.startswith("#"):
            outf.write(line)
            continue

        # Remove SV, CNV, and MEI lines
        if "SVTYPE" in line or "<CGA_CNVWIN>" in line or "<INS:ME:" in line:
            continue

        # Keep remaining lines
        (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE) = line.rstrip().split("\t")

        # Change ALT for nocalls to "."
        if ALT == "<CGA_NOCALL>":
            ALT = "."

        # Copy FT value to FILTER for anything that's not a nocall
        else:
            sample = {}
            for (key,val) in zip(FORMAT.split(":"), SAMPLE.split(":")):
                sample[key] = val

            if sample['GT'] != "./." and sample['GT'] != ".|." and sample['GT']!=".":
                if "FT" in sample.keys():
                    FILTER = sample['FT']


        outf.write("\t".join( (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE) ) + "\n")

###################################################

if __name__ == "__main__":
    main()
                

