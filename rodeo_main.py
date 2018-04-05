#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 20:34:38 2017

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

import logging
import multiprocessing
import shutil
import os
import argparse
import config_parser
import traceback
import sys
from shutil import copyfile

WEB_TOOL = False
if WEB_TOOL:
    RODEO_DIR = "/home/ubuntu/website/go/rodeo2/"
    os.chdir(RODEO_DIR)
VERSION = "2.1.2"
#VERBOSITY = logging.DEBUG
VERBOSITY = logging.INFO
QUEUE_CAP = "END_OF_QUEUE"
processes = []

def __main__():
    import nulltype_module
    import main_html_generator
    import ripp_html_generator
    from record_processing import fill_request_queue, Error_report
    import My_Record
    import record_processing
    
#==============================================================================
#     First we will handle the input, whether it be an accession, a list of acc
#   or a file which itself contains a list of accessions.
#==============================================================================
    parser = argparse.ArgumentParser("Main RODEO app.")
    parser.add_argument('query', type=str,
                        help='Accession number, genbank file or .txt file with an accession or .gbk query on each line') #accession # or gi
    parser.add_argument('-out', '--output_dir', type=str,
                        help='Name of output folder')
    parser.add_argument('-c', '--conf_file', nargs='*', default=[], 
                        help='Maximum size of potential ORF') 
    parser.add_argument('-hmm', '--custom_hmm', nargs='*', default=[], 
                        help='Maximum size of potential ORF') 
    parser.add_argument('-j', '--num_cores', type=int, default = 1,
                        help="Number of cores to use.")
#    parser.add_argument('-v', '--verbose', action='store_true', default=True,
#                        help="Include DEBUG output in stdout stream")
    parser.add_argument('-max', '--precursor_max', type=int,
                        help='Maximum size of potential ORF') #better word for potential?
    parser.add_argument('-min', '--precursor_min', type=int,
                        help='Minimum size of potential ORF') 
    parser.add_argument('-o', '--overlap', type=int, 
                        help='Maximum overlap of search with existing CDSs')
    parser.add_argument('-ft', '--fetch_type', type=str,
                        help='Type of window specification.\n' +
                        '\'cds\' will make the window +/- n CDSs from the query.\n' +
                        '\'nucs\' will make the window +/- n nucleotides from the query')
    parser.add_argument('-fn', '--fetch_n', type=int, 
                        help='The \'n\' variable for the -ft=orfs')
    parser.add_argument('-fd', '--fetch_distance', type=int, 
                        help='Number of nucleotides to fetch outside of window')
    parser.add_argument('-pt', '--peptide_types', nargs='*', default = [],
                        help='Type(s) of peptides to score.')
    parser.add_argument('-ea', '--evaluate_all', action='store_true', default=None,
                        help='Evaluate all duplicates if accession id corresponds to duplicate entries')
    parser.add_argument('-ex', '--exhaustive', action='store_true', default=None,
                        help="Score RiPPs even if they don't have a valid split site")
    parser.add_argument('-print', '--print_precursors', action='store_true', default=None,
                        help="Print precursors in HTML file")
    parser.add_argument('-w', '--web', action='store_true', default=False,
                        help="Only to use when running as a web tool")
    
    args, _ = parser.parse_known_args()
#==============================================================================
#     Set up logger
#==============================================================================
    
    logger = logging.getLogger("rodeo_main")
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
    
#==============================================================================
#   Handle configureations
#==============================================================================
    confs = []
    
    for conf in ['confs/default.conf'] + args.conf_file:
        try:
            confs.append(config_parser.parse_config_file(conf))
        except KeyboardInterrupt:
            raise KeyboardInterrupt
        except Exception as e:
            logger.error("Error with conf file %s" % (conf))
            print(e)
    master_conf = config_parser.merge_confs(confs)
    master_conf = config_parser.merge_conf_and_arg(master_conf, args)
    general_conf = master_conf['general']
    if WEB_TOOL:
        general_conf['variables']['pfam_dir'] = "/home/ubuntu/website/go/rodeo2/hmm_dir/Pfam-A.hmm"
        
    def pretty(d, indent=0):
        for key, value in d.items():
            print('\t' * indent + str(key))
            if isinstance(value, dict):
                pretty(value, indent+1)
            else:
                print('\t' * (indent+1) + str(value))
#    pretty(general_conf)
#==============================================================================
#   Set up output directory
#==============================================================================
    if args.output_dir == None:
        args.output_dir = args.query.split('.')[0] + "_rodeo_out"
    overwriting_folder = False
    try:
        os.mkdir(args.output_dir)
    except OSError:
        overwriting_folder = True
        shutil.rmtree(args.output_dir)
        os.mkdir(args.output_dir)
    
    try:
        os.mkdir(args.output_dir + '/confs')
        for conf_file in ['confs/default.conf'] + args.conf_file:
            try:
                name = conf_file.split('/')[-1]
                copyfile(conf_file, args.output_dir + '/confs/' + name)
            except:
                 logger.warning("Problem copying configuration file {}".format("conf_file"))
    except:
        logger.warning("Problem creating configuration copy directory")
    if overwriting_folder:
        logger.warning("Overwriting %s folder." % (args.output_dir))
    
    try: 
        os.mkdir("tmp_files")
    except OSError:
	pass
#==============================================================================
#   Check arguments        
#==============================================================================
    if not any(general_conf['variables']['fetch_type'].lower() == ft for ft in ['cds', 'nucs']):
        logger.critical("Invalid argument for -ft/-fetch_type")
        return None
    
    if 'sacti' in args.peptide_types or 'lanthi' in args.peptide_types:
        if not any ("tigr" in hmm_name.lower() for hmm_name in args.custom_hmm):
            logger.warn("Lanthi and/or sacti heuristics require TIGRFAM hmm. Make sure its location is specified with the -hmm or --custom_hmm flag.")
#==============================================================================
#   Set up queries/read query files   
#==============================================================================
    query = args.query
    queries = []
    
    if '.txt' == query[-4:]:
        try:
            input_handle = open(query)
            for line in input_handle:
                if line.strip() == "\n" or line.rstrip() == "": 
                    continue
                queries.append(line.rstrip())        
        except OSError:
            logger.critical("Could not find %s" % (query)) 
    else:
        queries.append(query)

#==============================================================================
#   Initializations and output file header writing        
#==============================================================================
    peptide_types = args.peptide_types
    output_dir = args.output_dir
    module = nulltype_module
    
    module.main_write_headers(output_dir)
    module.co_occur_write_headers(output_dir)
    main_html = open(output_dir + "/main_results.html", 'w')
    main_html_generator.write_header(main_html, master_conf)
    main_html_generator.write_table_of_contents(main_html, queries)
    ripp_modules = {}
    ripp_htmls = {}
    records = []
    
    for peptide_type in peptide_types:
        list_of_rows = []
        if peptide_type == "lasso":
            import ripp_modules.lasso.lasso_module as module
        elif peptide_type == "lanthi":
            import ripp_modules.lanthi.lanthi_module as module
        elif peptide_type == "sacti":
            import ripp_modules.sacti.sacti_module as module
        elif peptide_type == "thio":
            import ripp_modules.thio.thio_module as module
        else:
            logger.error("%s not in supported RiPP types" % (peptide_type))
            continue
        os.mkdir(output_dir + "/" + peptide_type)
        module.write_csv_headers(output_dir)
        ripp_modules[peptide_type] = module
        ripp_htmls[peptide_type] = open(output_dir + "/" + peptide_type + "/" + peptide_type + "_results.html", 'w')
        ripp_html_generator.write_header(ripp_htmls[peptide_type], master_conf, peptide_type)
        ripp_html_generator.write_table_of_contents(ripp_htmls[peptide_type], queries)

#==============================================================================
#   Set up paralellization (worker processes and fetch process)
#==============================================================================
    query_no = 0
    m = multiprocessing.Manager()
    processed_records_q = m.Queue(max(args.num_cores, 6))
    unprocessed_records_q = m.Queue(max(args.num_cores, 6))
    
    #Nest all in try in case of KeyboardInterrupt
    try:
        for i in range(args.num_cores):
            processes.append(multiprocessing.Process(target=record_processing.process_record_worker, args=(unprocessed_records_q, processed_records_q, args, master_conf, ripp_modules, i)))
            processes[-1].start()
        request_proc = multiprocessing.Process(target=fill_request_queue, args=(queries, processed_records_q, unprocessed_records_q, args, master_conf, ripp_modules,))
        request_proc.start()
        record = processed_records_q.get()
        queue_cap_count = 0
#==============================================================================
#       Main logic loop. Pull until you've pulled as many queue caps as there 
#       are worker processes.
#==============================================================================
        while queue_cap_count < args.num_cores:
            if record == QUEUE_CAP:
                queue_cap_count += 1
                if queue_cap_count == args.num_cores:
                    break
                record =  processed_records_q.get()
                continue
            query = record.query_accession_id
            query_no += 1
            logger.info("Writing output for query #%d.\t%s" % (query_no, query))
            if type(record) == Error_report:
                logger.error(("For %s:\t" + record.error_message) % (record.query))
                main_html_generator.write_failed_query(main_html, record.query, record.error_message)
                for peptide_type in peptide_types:
                    ripp_html_generator.write_failed_query(ripp_htmls[peptide_type], record.query, record.error_message)
                record = processed_records_q.get()
                continue
    
            module = nulltype_module
            for orf in record.intergenic_orfs:
                if orf.start < orf.end:
                    direction = "+"
                else:
                    direction = "-"
                row = [query, record.cluster_genus_species, record.cluster_accession, 
                       orf.start, orf.end, direction, orf.sequence]
                module.main_write_row(output_dir, row)
            for cds in record.CDSs:
                if cds.start < cds.end:
                    direction = "+"
                else:
                    direction = "-"
                row = [query, record.cluster_genus_species, record.cluster_accession,
                       cds.accession_id, cds.start, cds.end, direction]
                for pfam_acc, desc, e_val, name in cds.pfam_descr_list:
                    row += [pfam_acc, name, desc, e_val]
                module.co_occur_write_row(output_dir, row)
                
            main_html_generator.write_record(main_html, master_conf, record)
            
            for peptide_type in args.peptide_types:
                list_of_rows = []
                module = ripp_modules[peptide_type]
                for ripp in record.ripps[peptide_type]:
#                    if len(ripp.sequence) < master_conf[peptide_type]['variables']['precursor_min']:
#                        continue
#                    elif len(orf.sequence)  > master_conf[peptide_type]['variables']['precursor_max']:
#                        if not "M" in orf.sequence[2:]:
#                            continue
                    if not master_conf[peptide_type]['variables']['precursor_min'] <= len(ripp.sequence) <= master_conf[peptide_type]['variables']['precursor_max'] and \
                        not ("M" in ripp.sequence[-master_conf[peptide_type]['variables']['precursor_max']:]):
                        continue
                    list_of_rows.append(ripp.csv_columns)
                module.ripp_write_rows(args.output_dir, record.cluster_accession, #cluster acc or query acc?
                                       record.cluster_genus_species, list_of_rows)
            records.append(record)
            record = processed_records_q.get()    
            #END MAIN LOOP
        main_html.write("</html>")
        
        #Update score w SVM.
        
        try:
            for peptide_type in peptide_types:
                module = ripp_modules[peptide_type]
                module.run_svm(output_dir)
            My_Record.update_score_w_svm(output_dir, records)
        except KeyboardInterrupt:
            raise KeyboardInterrupt
        except IndexError:
            logger.error("No valid results. Input may be invalid or Genbank may not be responding")
        except ValueError as e:
            logger.critical("Value error when finalizing results. This most likely means that no results were obtained and that all queries failed.")
        except Exception as e:
            logger.error("Error running SVM")
            logger.error(e)
            traceback.print_exc(file=sys.stdout)
        
        for peptide_type in peptide_types: 
            for record in records:
                ripp_html_generator.write_record(ripp_htmls[peptide_type], master_conf, record, peptide_type)
            try:
                os.remove(output_dir + "/" + peptide_type + "/" + "temp_features.csv")
                pass
            except OSError:
                logger.debug("Temp feature file appears to be missing...")
         
        
    except KeyboardInterrupt as e:
        logger.critical("Keyboard interrupt recieved. Shutting down RODEO.")
        request_proc.join(5)
        for process in processes:
            process.join(5)
    
if __name__ == '__main__':
    __main__()

