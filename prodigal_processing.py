# -*- coding: utf-8 -*-
#==============================================================================
# Copyright (C) 2017 Bryce L. Kille
# University of Illinois
# Department of Chemistry
#
# Copyright (C) 2017 Christopher J. Schwalen
# University of Illinois
# Department of Chemistry
#
# Copyright (C) 2017 Douglas A. Mitchell
# University of Illinois
# Department of Chemistry
#
# License: GNU Affero General Public License v3 or later
# Complete license availabel in the accompanying LICENSE.txt.
# or <http://www.gnu.org/licenses/>.
#
# This file is part of RODEO2.
#
# RODEO2 is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# RODEO2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#==============================================================================

import csv, subprocess, os, statistics
import logging
from rodeo_main import VERBOSITY

logger = logging.getLogger(__name__)
logger.setLevel(VERBOSITY)
# create console handler and set level
ch = logging.StreamHandler()
ch.setLevel(VERBOSITY)

# create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# add formatter to ch
ch.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)

def run_prodigal(record):
    try:
        prod_prefix = "tmp_files/" + record.query_short
        prod_os_file = open(prod_prefix+"os.txt", 'w')
        prod_use_file = open(prod_prefix+"prod.fasta", 'w')
        prod_use_file.write(">%s %s %s\n%s" % (record.query_accession_id, record.genus, record.species, 
        					record.cluster_sequence[record.window_start:record.window_end]))
        prod_use_file.close()
        if record.prod_window_end-record.prod_window_start < 100000:
            process = subprocess.Popen(["prodigal", "-i", prod_prefix+"prod.fasta", "-o", prod_prefix+"output.txt", 
            							"-p", "meta", "-s", prod_prefix+"orfs.tsv", "-q"], stdout=prod_os_file, stderr=prod_os_file)
            process.wait()
        else:
            prod_train_file = open(prod_prefix+"train.fasta", 'w')
            prod_train_file.write(">%s %s %s\n%s" % (record.query_accession_id, record.genus, record.species, 
            						record.cluster_sequence[record.prod_window_start:record.prod_window_end]))
            prod_train_file.close()
            process = subprocess.Popen(["prodigal", "-i", prod_prefix+"train.fasta", "-o", prod_prefix+"output.txt", 
            							"-t", prod_prefix+"train.txt", "-q"], stdout=prod_os_file, stderr=prod_os_file)
            process.wait()
            process = subprocess.Popen(["prodigal", "-i", prod_prefix+"prod.fasta", "-o", prod_prefix+"output.txt", 
            							"-t", prod_prefix+"train.txt", "-s", prod_prefix+"orfs.tsv", "-q"], stdout=prod_os_file, stderr=prod_os_file)
            process.wait()
             
        try:
            os.remove(prod_prefix+"train.fasta")
        except:
            pass
        try:
            os.remove(prod_prefix+"output.txt")
        except:
            pass
        try:
            os.remove(prod_prefix+"prod.fasta")
        except:
            pass
        try:
            os.remove(prod_prefix+"train.txt")
        except:
            pass
        try:
        	os.remove(prod_prefix+"os.txt")
        except:
        	pass
    except KeyboardInterrupt:
        try:
            os.remove(prod_prefix+"train.fasta")
        except:
            pass
        try:
            os.remove(prod_prefix+"output.txt")
        except:
            pass
        try:
            os.remove(prod_prefix+"prod.fasta")
        except:
            pass
        try:
            os.remove(prod_prefix+"train.txt")
        except:
            pass
        try:
            os.remove(prod_prefix+"orfs.tsv")
        except:
            pass
        try:
        	os.remove(prod_prefix+"os.txt")
        except:
        	pass
        logger.critical("SIGINT recieved during Prodigal")
        raise KeyboardInterrupt

def swarm_filter(output_dir, filtration_length=75, inclusion_length=150):
    input_file = open(output_dir + "/prodigal/prod_results.csv", 'r')
    output_fasta_file = open(output_dir + "/prodigal/filtered_sequences.txt", 'w')
    output_csv_file = open(output_dir + "/prodigal/filtered_prod_results.csv", 'w')

    sequence_data = input_file.readlines()
    output_csv_file.write(sequence_data[0])

    maxtotalscore = 0.0
    maxcodpotscore = 0.0
    maxstrtscore = 0.0
    maxrbsscore = 0.0
    maxupsscore = 0.0
    maxtypescore = 0.0
    adjtotalscore = 0.0
    adjcodpotscore = 0.0
    adjstrtscore = 0.0
    adjrbsscore = 0.0
    adjupsscore = 0.0
    adjtypescore = 0.0    

    for line in sequence_data[1:]:
        entry = line.split(",")
        if len(entry[16]) > filtration_length:
            continue
        maxtotalscore = max(maxtotalscore, float(entry[7]))
        maxcodpotscore = max(maxcodpotscore, float(entry[8]))
        maxstrtscore = max(maxstrtscore, float(entry[9]))
        maxrbsscore = max(maxrbsscore, float(entry[13]))
        maxupsscore = max(maxupsscore, float(entry[14]))
        maxtypescore = max(maxtypescore, float(entry[15]))    

    for line in sequence_data[1:]:
        entry = line.split(",")
        if len(entry[16]) > inclusion_length: 
            continue
        adjtotalscore = statistics.median([float(entry[7])/maxtotalscore, 0.0, 1.0])
        adjcodpotscore = statistics.median([float(entry[8])/maxcodpotscore, 0.0, 1.0])
        adjstrtscore = statistics.median([float(entry[9])/maxstrtscore, 0.0, 1.0])
        adjrbsscore = statistics.median([float(entry[13])/maxrbsscore, 0.0, 1.0])
        adjupsscore = statistics.median([float(entry[14])/maxupsscore, 0.0, 1.0])
        adjtypescore = statistics.median([float(entry[15])/maxtypescore, 0.0, 1.0])

        if adjtotalscore + adjcodpotscore + adjstrtscore + adjrbsscore + adjupsscore + adjupsscore + adjtypescore > 1.0:
            output_fasta_file.write(">%s\n%s" % (entry[0], entry[16]))
            output_csv_file.write(line)

def create_ssn(output_dir):
    process = subprocess.Popen(["diamond", "makedb", "--in", output_dir + "/prodigal/filtered_sequences.txt", "-d", output_dir + "/prodigal/align_db"])
    process.wait()
    process = subprocess.Popen(["diamond", "blastp", "-d", output_dir + "/prodigal/align_db", "-q", output_dir + "/prodigal/filtered_sequences.txt", "-o", output_dir + "/prodigal/blastp_results.tsv", "--min-score", "0"])
    process.wait()

    blastp_file = open(output_dir + "/prodigal/blastp_results.tsv", 'r')
    ssn_file = open(output_dir + "/prodigal/precursor.ssn", 'w')
    ssn_file.write("Node1\tNode2\tBitScore\n")

    blastp_results = blastp_file.readlines()
    for line in blastp_results:
    	result = line.split("\t")
    	if result[0] == result[1]:
    		continue
    	ssn_file.write("%s\t%s\t%s\n" % (result[0], result[1], result[11]))

    blastp_file.close()
    ssn_file.close()

    os.remove(output_dir + "/prodigal/align_db.dmnd")
    os.remove(output_dir + "/prodigal/blastp_results.tsv")
    	


    