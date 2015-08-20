#!/usr/bin/env python

############################################################################
# Copyright (c) 2011-2013 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import sys
import getopt
import os
import logging
import shutil
import datetime
from time import gmtime, strftime

home_directory = os.path.abspath(os.path.dirname(os.path.realpath(__file__))) + '/'
ig_bin_directory = os.path.join(home_directory, "bin/")
python_src_directory = os.path.join(home_directory, "src/ig_tools/python_utils/")
ig_blast_directory = os.path.join(home_directory, "src/tools/igblast/")
inexact_repertoire_evaluator_src_directory = os.path.join(home_directory, "src/ig_tools/inexact_repertoire_evaluator/")
config_directory = os.path.join(home_directory, "configs/ig_tools/")
spades_py_scripts_directory = os.path.join(home_directory, "src/spades_pipeline/")
spades_release_bins = os.path.join(home_directory, "build/release/bin/")
somatic_mutation_search_src_dir = os.path.join(home_directory, "src/ig_tools/somatic_mutation_search") 

path_to_config_template = os.path.join(config_directory, "config.info.template")

sys.path.append(python_src_directory)
sys.path.append(inexact_repertoire_evaluator_src_directory)
sys.path.append(spades_py_scripts_directory)
sys.path.append(somatic_mutation_search_src_dir)

class PathToBins:
    repertoire_evaluator_tool = os.path.join(ig_bin_directory, "repertoire_evaluator")
    repertoire_comparer_tool = os.path.join(ig_bin_directory, "repertoire_comparer")
    ig_matcher_tool = os.path.join(ig_bin_directory, "ig_matcher")
    ig_kplus_vj_finder_tool = os.path.join(ig_bin_directory, "ig_matcher")

    run_repertoire_evaluator_tool = ig_bin_directory + "./repertoire_evaluator"
    run_repertoire_comparer_tool = ig_bin_directory + "./repertoire_comparer"
    run_ig_matcher_tool = ig_bin_directory + "./ig_matcher"
    run_ig_kplus_vj_finder_tool = ig_bin_directory + "./ig_kplus_vj_finder"

def PrintCommandLine(argv, log):
    command_line = " ".join([str(x) for x in argv] )
    log.info("\nCommand line: "+ command_line)

def ReadConfig():
    if not os.path.exists(path_to_config_template):
        print("ERROR: config file " + path_to_config_template + " was not found")
    f = open(path_to_config_template, "r")
    config_params = dict()
    for line in f.readlines():
        splits = line.split()
        config_params[splits[0]] = splits[1]
    return config_params

def RunIgblast():
    return os.path.join(ig_blast_directory, "/bin/igblastn")

def BlastDB():
    return os.path.join(home_directory, "src/tools/ncbi_blast_bin/filtered_imgt_blast_db/filtered_imgt_blastdb")

class Unbuffered:
   def __init__(self, stream):
           self.stream = stream
   def write(self, data):
           self.stream.write(data)
           self.stream.flush()
   def __getattr__(self, attr):
           return getattr(self.stream, attr)

def ErrorMsg(log):
    log.error("Something goes wrong. Please contact us and send .log file")
    sys.exit(1)

def AbnormalFinishMsg(log, program_name):
    log.info("Script " + program_name + " finished abnormally. Please contact us and send .log file")
