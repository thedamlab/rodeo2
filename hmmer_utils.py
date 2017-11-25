#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 21:38:20 2017

@author: bryce
"""

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

import os
import csv
import subprocess
import logging 
from rodeo_main import VERBOSITY

logger = logging.getLogger(__name__)
logger.setLevel(VERBOSITY)
# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(VERBOSITY)

# create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# add formatter to ch
ch.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)
pid_prefix = None    

#TODO try/except blocks for error checking processes
def _generate_fasta(accession):
    out_file = open("fasta_file.tmp.fasta", 'w+')
    esearch_process = subprocess.Popen(["/home/bryce/edirect/esearch", "-db", "protein", "-query", accession],
                                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    efetch_process = subprocess.Popen(["/home/bryce/edirect/efetch", "-format", "fasta"], 
                                      stdin=esearch_process.stdout, stdout=out_file)
    efetch_process.communicate()
    out_file.close()

def _generate_fasta_from_string(query_string):
    handle = open(pid_prefix+"fasta_file.tmp.fasta", "w+")
    handle.write(">temp_string\n" + query_string)

def _generate_hmmer(hmm_file):
    hmmer_process = subprocess.call(["hmmscan", "-o", pid_prefix+"hmm_out.tmp.tab", "--noali",
                                      "--domtblout", pid_prefix+"pFamInfo.tmp.tab", hmm_file,
                                      pid_prefix+"fasta_file.tmp.fasta"])    

#TODO try/except blocks for parsing to double check format
def _parse(n, e_cutoff):
    ret_list = []
    NUM_COLUMNS = 23
    pfam_accession = ''
    with open(pid_prefix+'pFamInfo.tmp.tab') as handle:
        for line in csv.reader(handle, delimiter=' '):
            if line[0][0] == '#':
                continue
            tmp_list = []
            for item in line:
                if item != '':
                    tmp_list.append(item)
            if tmp_list[1] != pfam_accession: #Check to make sure its not the same as the last one
                pfam_accession = tmp_list[1].split('.')[0]
                description = ' '.join(tmp_list[NUM_COLUMNS-1:])
                e_val = tmp_list[6]
                if float(e_val) < e_cutoff and len(ret_list) < n: #Will go through one extra line but oh well #TODO
                    ret_list.append((pfam_accession, description, float(e_val)))
                else:
                    break
    ret_list_final = []
    for acc, desc, e_val in ret_list:
        in_list = False
        for final_acc, _, _ in ret_list_final:
            if acc == final_acc:
                in_list = True
        if not in_list:
            ret_list_final.append((acc, desc, e_val))
    return ret_list_final

def get_hmmer_info(query, primary_hmm, cust_hmm, n=5, e_cutoff=.001, query_is_accession=False): #TODO handle lists of accessions
    """Returns top n hmmscan hits with e_values lower than e_cutoff"""
    try:
        global pid_prefix
        pid_prefix = "tmp_files/" + str(os.getpid())
        if query_is_accession:
            _generate_fasta(query)
        else:
            _generate_fasta_from_string(query)
        pfam_desc_list = []
        try:
            _generate_hmmer(primary_hmm)
            pfam_desc_list = _parse(n, e_cutoff) 
        except OSError:
            logger.error("Couldn't find %s" % (primary_hmm)) 
        
        for hmm in cust_hmm:
            try:
                _generate_hmmer(hmm)
                pfam_desc_list += _parse(n, e_cutoff) 
            except OSError:
                logger.error("Couldn't find %s" % (hmm))
        
        pfam_desc_list = sorted(pfam_desc_list, key=lambda entry: entry[2])
        os.remove(pid_prefix+'pFamInfo.tmp.tab')
        os.remove(pid_prefix+"fasta_file.tmp.fasta")
        os.remove(pid_prefix+"hmm_out.tmp.tab")
        return pfam_desc_list
    except KeyboardInterrupt:
        logger.critical("SIGINT recieved during HMMScan")
        raise KeyboardInterrupt
