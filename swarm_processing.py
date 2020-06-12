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

def swarm_filter(output_dir, filtration_length=75, inclusion_length=150):
	input_file = open(output_dir + "/prodigal/prod_results.csv", 'r')
	output_fasta_file = open(output_dir + "/prodigal/filtered_sequences.txt", 'w')
	output_csv_file = open(output_dir + "/prodigal/filtered_prod_results.csv", 'w')

	sequence_data = input_file.readlines()
	output_csv_file.write(sequence_data[0])
#Notes: total score float(entry[6]
#		cod pot score float(entry[7]
#		strt score float(entry [8]

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
		if len(entry[15]) > filtration_length:
			continue
		maxtotalscore = max(maxtotalscore, float(entry[6]))
		maxcodpotscore = max(maxcodpotscore, float(entry[7]))
		maxstrtscore = max(maxstrtscore, float(entry[8]))
		maxrbsscore = max(maxrbsscore, float(entry[12]))
		maxupsscore = max(maxupsscore, float(entry[13]))
		maxtypescore = max(maxtypescore, float(entry[14]))	

	for line in sequence_data[1:]:
		entry = line.split(",")
		if len(entry[15]) > inclusion_length: 
			continue
		adjtotalscore = statistics.median([float(entry[6])/maxtotalscore, 0.0, 1.0])
		adjcodpotscore = statistics.median([float(entry[7])/maxcodpotscore, 0.0, 1.0])
		adjstrtscore = statistics.median([float(entry[8])/maxstrtscore, 0.0, 1.0])
		adjrbsscore = statistics.median([float(entry[12])/maxrbsscore, 0.0, 1.0])
		adjupsscore = statistics.median([float(entry[13])/maxupsscore, 0.0, 1.0])
		adjtypescore = statistics.median([float(entry[14])/maxtypescore, 0.0, 1.0])

		if adjtotalscore + adjcodpotscore + adjstrtscore + adjrbsscore + adjupsscore + adjupsscore + adjtypescore > 1.0:
			output_fasta_file.write(">%s_Start:%s_End:%s\n%s" % (entry[0], entry[3], entry[4], entry[15]))
			output_csv_file.write(line)


	