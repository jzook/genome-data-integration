#!/usr/bin/env python

import os
import dxpy
import subprocess
import gzip
import shutil
import glob
import random 
import sys
import re
import math
import argparse
import numpy
import sys
import pipeline_runner


@dxpy.entry_point('main')
def main(vcfs, beds,filtbeds, annotations, callsettable, ref, rtgsdf, refn, chrom, prefix, tandemrepeatsbed, techs_with_no_callable):

    all_inputs = dxpy.download_all_inputs(parallel = True)
    print(os.getcwd())
    print(all_inputs)
    print(chrom)
    print(prefix)
    integrator_obj = pipeline_runner.Integrator(all_inputs['callsettable_path'][0], chrom, all_inputs['refn_path'][0], prefix, tandemrepeatsbed, techs_with_no_callable)
    integrator_obj.import_files()
    integrator_obj.preprocess()
    integrator_obj.classify_and_intersect()
    integrator_obj.run_classifier()
    integrator_obj.classify_final()
    integrator_obj.prep_bed_and_summarize()
    output = integrator_obj.upload_data()

    return output

dxpy.run()