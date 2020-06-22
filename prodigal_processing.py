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

def swarm_filter(output_dir, filtration_length=75, inclusion_length=120):
    input_file = open(output_dir + "/prodigal/prod_results.csv", 'r')
    output_fasta_file = open(output_dir + "/prodigal/filtered_sequences.txt", 'w')
    output_fasta_file2 = open(output_dir + "/prodigal/unfiltered_sequences.txt", 'w')
    output_csv_file = open(output_dir + "/prodigal/filtered_prod_results.csv", 'w')
    output_csv_file2 = open(output_dir + "/prodigal/unfiltered_prod_results.csv", "w")

    sequence_data = input_file.readlines()

    header = sequence_data[0].split(",")
    line_to_print = ",".join([header[0], header[1], header[2], header[3], header[4], header[5], header[6], header[16][:-1], "Length",
                            header[7], "Adj_" + header[7], header[8], "Adj_" + header[8], header[9], 
                            "Adj_" + header[9], header[10], header[11], header[12], header[13],
                            "Adj_" + header[13], header[14], "Adj_" + header[14], header[15], 
                            "Adj_" + header[15], "Sum_Adj_Scores\n"])

    output_csv_file.write(line_to_print)
    output_csv_file2.write(line_to_print)

    maxtotalscore = 0.1
    maxcodpotscore = 0.1
    maxstrtscore = 0.1
    maxrbsscore = 0.1
    maxupsscore = 0.1
    maxtypescore = 0.1
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

        logit = adjtotalscore + adjcodpotscore + adjstrtscore + adjrbsscore + adjupsscore + adjtypescore

        line_to_print = ",".join([entry[0], entry[1], entry[2], entry[3], entry[4], entry[5], entry[6], entry[16][:-1], str(len(entry[16])),
                                    entry[7], str(adjtotalscore), entry[8], str(adjcodpotscore), entry[9], 
                                    str(adjstrtscore), entry[10], entry[11], entry[12], entry[13],
                                    str(adjrbsscore), entry[14], str(adjupsscore), entry[15], 
                                    str(adjtypescore), str(logit)+"\n"])

        if logit > 1:
            output_fasta_file.write(">%s\n%s" % (entry[0], entry[16]))
            output_csv_file.write(line_to_print)

        output_csv_file2.write(line_to_print)
        output_fasta_file2.write(">%s\n%s" % (entry[0], entry[16]))

def create_ssn(output_dir):
    process = subprocess.Popen(["diamond", "makedb", "--in", output_dir + "/prodigal/filtered_sequences.txt", "-d", output_dir + "/prodigal/align_db"])
    process.wait()
    process = subprocess.Popen(["diamond", "blastp", "-d", output_dir + "/prodigal/align_db", "-q", output_dir + "/prodigal/filtered_sequences.txt", 
                                "-o", output_dir + "/prodigal/blastp_results.tsv", "--min-score", "0", "--more-sensitive"])
    process.wait()

    process = subprocess.Popen(["diamond", "makedb", "--in", output_dir + "/prodigal/unfiltered_sequences.txt", "-d", output_dir + "/prodigal/align_db2"])
    process.wait()
    process = subprocess.Popen(["diamond", "blastp", "-d", output_dir + "/prodigal/align_db2", "-q", output_dir + "/prodigal/unfiltered_sequences.txt", 
                                "-o", output_dir + "/prodigal/blastp_results2.tsv", "--min-score", "0", "--more-sensitive"])
    process.wait()

    blastp_file = open(output_dir + "/prodigal/blastp_results.tsv", 'r')
    ssn_file = open(output_dir + "/prodigal/filtered_precursor.ssn", 'w')
    ssn_file.write("Node1\tNode2\tBitScore\n")

    blastp_results = blastp_file.readlines()
    for line in blastp_results:
    	result = line.split("\t")
#    	if result[0] == result[1]:
#    		continue
    	ssn_file.write("%s\t%s\t%s\n" % (result[0], result[1], result[11]))

    blastp_file.close()
    ssn_file.close()

    blastp_file = open(output_dir + "/prodigal/blastp_results2.tsv", 'r')
    ssn_file = open(output_dir + "/prodigal/unfiltered_precursor.ssn", 'w')
    ssn_file.write("Node1\tNode2\tBitScore\n")

    blastp_results = blastp_file.readlines()
    for line in blastp_results:
        result = line.split("\t")
#        if result[0] == result[1]:
#            continue
        ssn_file.write("%s\t%s\t%s\n" % (result[0], result[1], result[11]))

    blastp_file.close()
    ssn_file.close()    

    os.remove(output_dir + "/prodigal/align_db.dmnd")
    os.remove(output_dir + "/prodigal/blastp_results.tsv")
    os.remove(output_dir + "/prodigal/align_db2.dmnd")
    os.remove(output_dir + "/prodigal/blastp_results2.tsv")


#def identify_precursor_clusters(output_dir):
#    edges = open(output_dir + "/prodigal/precursor.ssn", 'r').readlines()
#    nodes = []
#    for edge in edges:
#        data = edge.split("\t")
#        if data[0] not in nodes:
#            nodes.append(data[0])
#        if data[1] not in nodes:
#            nodes.append(data[1])
#
#    matrix = [[0 for i in range(0, len(nodes))] for j in range(0, len(nodes))]
#    for edge in edges:
#        data = edge.split("\t")
#        index1 = nodes.index(data[0])
#        index2 = nodes.index(data[1])
#        matrix[index1][index2] = data[3]
#
#    print(matrix)

    	


    