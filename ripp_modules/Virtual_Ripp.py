#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 23 15:09:17 2017

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

import csv
import os
import re
import tempfile
import subprocess
import logging
from rodeo_main import VERBOSITY
VERBOSITY = logging.DEBUG

WEB_TOOL = False


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

def execute(commands, input=None):
    "Execute commands in a system-independent manner"

    if input is not None:
        stdin_redir = subprocess.PIPE
    else:
        stdin_redir = None

    try:
        proc = subprocess.Popen(commands, stdin=stdin_redir,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, 
                                shell=True)
        out, err = proc.communicate(input=input)
        retcode = proc.returncode
        return out, err, retcode
    except OSError, e:
        logger.error(e)
        raise e
        
        
class Virtual_Ripp(object):
    def __init__(self,
                 start,
                 end,
                 sequence,
                 upstream_sequence,
                 pfam_2_coords):
        self.start = start
        self.end = end
        self.sequence = sequence
        self.upstream_sequence = upstream_sequence
        #This line is here to ensure that the type is always set.
        #No peptide type should ever be left virtual. 
        self.peptide_type = 'virtual' 
        self.score  = 0
        self.pfam_2_coords = pfam_2_coords
        if start < end:
            self.direction = '+'
        else:
            self.direction = '-'
	self.valid_split = True #set to false if no valid split

#        self.set_leader_core()
#        self.set_monoisotopic_mass()
#        self.csv_columns = [self.leader, self.core, self.start, self.end]

    def run_svm(self, output_dir):
        if self.peptide_type == 'lasso':
            from ripp_modules.lasso.svm import svm_classify as svm
            from ripp_modules.lasso import CUTOFF as cutoff
        elif self.peptide_type == 'sacti':
            from ripp_modules.sacti.svm import svm_classify as svm
            from ripp_modules.sacti import CUTOFF as cutoff
        elif self.peptide_type == 'thio':
            from ripp_modules.thio.svm import svm_classify as svm
            from ripp_modules.thio import CUTOFF as cutoff
        elif self.peptide_type == 'lanthi':
            from ripp_modules.lanthi.svm import svm_classify as svm
            from ripp_modules.lanthi import CUTOFF as cutoff
            
        svm.run_svm()
        svm_output_reader = csv.reader(open("ripp_modules/" + self.peptide_type + "/svm/fitting_results.csv"))
        final_output_writer = csv.writer(open(output_dir + "/" + self.peptide_type + '/'\
                                              + self.peptide_type + "_features.csv", 'w'))
        features_reader = csv.reader(open(output_dir + "/" + self.peptide_type + '/'\
                                              + self.peptide_type + "/temp_features.csv"))
        header_row = features_reader.next() #skip header
        final_output_writer.writerow(header_row)
        
        #need to make ubiquitous 
        for row in features_reader:
            svm_output = svm_output_reader.next()[1]
            row.append(svm_output)
            if int(svm_output) == 1:
                row[-2] = int(row[-2]) + 10
            if int(row[-2]) > cutoff: #CUTOFF
                row[7] = 'Y'
            else:
                row[7] = 'N'
            final_output_writer.writerow(row)
        
    def run_fimo_simple(self, query_motif_file=None):
    #TODO change to temp file
        if not query_motif_file:
            query_motif_file = "ripp_modules/" + self.peptide_type + '/' + self.peptide_type + "_fimo.txt"
        pid = str(os.getpid())
        try:
            with open("tmp_files/" + pid + "FIMO.seq", 'w+') as tfile:
                tfile.write(">query\n%s" % (self.sequence))
            if WEB_TOOL:
                command = ["/home/ubuntu/meme/bin/fimo --text --verbosity 1 " + query_motif_file + ' ' + "tmp_files/" + pid + "FIMO.seq"]
            else:
                command = ["$HOME/meme/bin/fimo --text --thresh 0.01 --verbosity 1 " + query_motif_file + ' ' + "tmp_files/" + pid + "FIMO.seq"]
            try:
                out, err, retcode = execute(command)
            except OSError:
                logger.error("Could not run FIMO on %s" % (self.peptide_type))
                try:
                    os.remove(pid+"TEMP.seq")
                except OSError:
                    pass
                return ""
            if retcode != 0:
                logger.error('FIMO returned %d: %r while searching %r', retcode,
                                err, query_motif_file)
                return []
            try:
                os.remove("tmp_files/" + pid + "FIMO.seq")
            except OSError:
                    pass
        except KeyboardInterrupt:
            try:
                os.remove("tmp_files/" + pid + "FIMO.seq")
                return
            except OSError:
                pass
        return out
    
    def get_min_dist(self, coords_list):
        if coords_list == []:
            return None
        min_dist = abs(self.start-coords_list[0][0])
        for coord in coords_list:
            min_dist = min(abs(self.start-coord[0]), abs(self.end-coord[0]),
                           abs(self.start-coord[1]), abs(self.end-coord[1]),
                           min_dist)
        return min_dist
