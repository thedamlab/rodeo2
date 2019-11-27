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
#
#==============================================================================
# Special thanks goes to AntiSmash team, whose antiSMASH-rodeo repository at
# https://bitbucket.org/mmedema/antismash-rodeo/ provided the backbone code for 
# a great deal of the heuristic calculations.
#==============================================================================


import csv
import os
import re
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import numpy as np
from ripp_modules.VirtualRipp import VirtualRipp
import hmmer_utils

peptide_type = "lanthi_i"
CUTOFF = 20
index = 0

def write_csv_headers(output_dir):
    dir_prefix = output_dir + '/{}/'.format(peptide_type)
    if not os.path.exists(dir_prefix):
        os.makedirs(dir_prefix)
    svm_headers = 'PK,Classification, 1/5_frac_c, 1/5_frac_s_t_c, 2/5_frac_c, 2/5_frac_s_t_c, 3/5_frac_c, 3/5_frac_s_t, 4/5_frac_c, 4/5_frac_s_t, 4/5_frac_s_t_c, 5/5_frac_c, 5/5_frac_s_t, 5/5_frac_s_t_c, S+T, S+T+C, aliphatic, aromatic, charged, core_A, core_AC, core_AG, core_AH, core_AL, core_AM, core_AP, core_AR, core_AT, core_AV, core_AW, core_AxA, core_AxC, core_AxH, core_AxL, core_AxN, core_AxP, core_AxR, core_AxT, core_AxV, core_AxW, core_AxxA, core_AxxG, core_AxxH, core_AxxL, core_AxxP, core_AxxR, core_AxxT, core_AxxxA, core_AxxxH, core_AxxxL, core_AxxxN, core_AxxxP, core_AxxxR, core_AxxxxA, core_AxxxxC, core_AxxxxL, core_AxxxxM, core_AxxxxP, core_AxxxxR, core_AxxxxS, core_AxxxxW, core_AxxxxxA, core_AxxxxxL, core_AxxxxxN, core_AxxxxxP, core_AxxxxxR, core_AxxxxxT, core_AxxxxxW, core_C, core_CA, core_CC, core_CD, core_CE, core_CG, core_CI, core_CK, core_CN, core_CP, core_CR, core_CS, core_CT, core_CV, core_CW, core_CxC, core_CxD, core_CxE, core_CxG, core_CxK, core_CxN, core_CxP, core_CxR, core_CxS, core_CxT, core_CxV, core_CxxC, core_CxxD, core_CxxG, core_CxxK, core_CxxN, core_CxxQ, core_CxxR, core_CxxS, core_CxxT, core_CxxV, core_CxxY, core_CxxxC, core_CxxxD, core_CxxxG, core_CxxxH, core_CxxxI, core_CxxxK, core_CxxxN, core_CxxxR, core_CxxxS, core_CxxxT, core_CxxxxA, core_CxxxxC, core_CxxxxD, core_CxxxxG, core_CxxxxI, core_CxxxxN, core_CxxxxR, core_CxxxxS, core_CxxxxT, core_CxxxxxC, core_CxxxxxF, core_CxxxxxG, core_CxxxxxI, core_CxxxxxK, core_CxxxxxN, core_CxxxxxQ, core_CxxxxxR, core_CxxxxxS, core_CxxxxxT, core_D, core_DC, core_DD, core_DF, core_DG, core_DH, core_DN, core_DP, core_DR, core_DS, core_DT, core_DxC, core_DxE, core_DxG, core_DxN, core_DxR, core_DxT, core_DxxC, core_DxxG, core_DxxN, core_DxxS, core_DxxT, core_DxxxC, core_DxxxD, core_DxxxG, core_DxxxK, core_DxxxN, core_DxxxR, core_DxxxS, core_DxxxT, core_DxxxxC, core_DxxxxL, core_DxxxxN, core_DxxxxR, core_DxxxxS, core_DxxxxT, core_DxxxxY, core_DxxxxxC, core_DxxxxxD, core_DxxxxxE, core_DxxxxxN, core_DxxxxxR, core_DxxxxxT, core_ED, core_ER, core_ET, core_ExR, core_ExT, core_ExxC, core_ExxR, core_ExxT, core_ExxxD, core_ExxxN, core_ExxxxR, core_ExxxxT, core_ExxxxxC, core_ExxxxxH, core_ExxxxxN, core_ExxxxxR, core_ExxxxxT, core_F, core_FA, core_FF, core_FI, core_FL, core_FP, core_FQ, core_FR, core_FT, core_FxA, core_FxC, core_FxF, core_FxI, core_FxL, core_FxP, core_FxR, core_FxxC, core_FxxL, core_FxxN, core_FxxR, core_FxxS, core_FxxxC, core_FxxxF, core_FxxxI, core_FxxxL, core_FxxxR, core_FxxxS, core_FxxxV, core_FxxxW, core_FxxxxC, core_FxxxxG, core_FxxxxI, core_FxxxxL, core_FxxxxR, core_FxxxxxF, core_FxxxxxI, core_FxxxxxK, core_FxxxxxL, core_FxxxxxR, core_FxxxxxT, core_GC, core_GD, core_GE, core_GH, core_GL, core_GN, core_GP, core_GR, core_GS, core_GT, core_GW, core_GxA, core_GxC, core_GxL, core_GxN, core_GxP, core_GxR, core_GxS, core_GxT, core_GxV, core_GxxA, core_GxxC, core_GxxI, core_GxxL, core_GxxN, core_GxxP, core_GxxR, core_GxxT, core_GxxV, core_GxxxA, core_GxxxC, core_GxxxK, core_GxxxL, core_GxxxN, core_GxxxP, core_GxxxR, core_GxxxS, core_GxxxT, core_GxxxxA, core_GxxxxC, core_GxxxxH, core_GxxxxK, core_GxxxxL, core_GxxxxN, core_GxxxxP, core_GxxxxR, core_GxxxxS, core_GxxxxT, core_GxxxxxC, core_GxxxxxD, core_GxxxxxL, core_GxxxxxN, core_GxxxxxP, core_GxxxxxR, core_GxxxxxS, core_GxxxxxT, core_H, core_HA, core_HC, core_HL, core_HP, core_HR, core_HxA, core_HxG, core_HxR, core_HxxA, core_HxxC, core_HxxG, core_HxxL, core_HxxS, core_HxxxA, core_HxxxC, core_HxxxL, core_HxxxP, core_HxxxR, core_HxxxS, core_HxxxxA, core_HxxxxL, core_HxxxxR, core_HxxxxxK, core_HxxxxxL, core_HxxxxxP, core_HxxxxxR, core_IC, core_IF, core_IR, core_IT, core_IxC, core_IxL, core_IxR, core_IxS, core_IxT, core_IxxA, core_IxxC, core_IxxF, core_IxxT, core_IxxxR, core_IxxxS, core_IxxxT, core_IxxxxC, core_IxxxxF, core_IxxxxR, core_IxxxxT, core_IxxxxxC, core_IxxxxxF, core_IxxxxxL, core_IxxxxxR, core_IxxxxxT, core_KC, core_KT, core_KxC, core_KxG, core_KxL, core_KxR, core_KxS, core_KxT, core_KxxA, core_KxxC, core_KxxD, core_KxxR, core_KxxT, core_KxxxC, core_KxxxE, core_KxxxL, core_KxxxT, core_KxxxxC, core_KxxxxD, core_KxxxxE, core_KxxxxG, core_KxxxxR, core_KxxxxT, core_KxxxxxC, core_KxxxxxG, core_KxxxxxL, core_KxxxxxR, core_KxxxxxS, core_L, core_LA, core_LF, core_LH, core_LK, core_LL, core_LP, core_LR, core_LS, core_LV, core_LxF, core_LxG, core_LxH, core_LxL, core_LxM, core_LxN, core_LxP, core_LxQ, core_LxR, core_LxS, core_LxV, core_LxY, core_LxxF, core_LxxG, core_LxxH, core_LxxK, core_LxxL, core_LxxP, core_LxxR, core_LxxS, core_LxxV, core_LxxxA, core_LxxxF, core_LxxxH, core_LxxxL, core_LxxxP, core_LxxxR, core_LxxxS, core_LxxxxA, core_LxxxxC, core_LxxxxG, core_LxxxxH, core_LxxxxI, core_LxxxxL, core_LxxxxP, core_LxxxxR, core_LxxxxS, core_LxxxxV, core_LxxxxxA, core_LxxxxxF, core_LxxxxxH, core_LxxxxxI, core_LxxxxxL, core_LxxxxxM, core_LxxxxxQ, core_LxxxxxR, core_LxxxxxV, core_LxxxxxW, core_LxxxxxY, core_M, core_MC, core_MV, core_MW, core_MxR, core_MxxP, core_MxxR, core_MxxT, core_MxxxR, core_MxxxxR, core_N, core_NC, core_NG, core_NR, core_NS, core_NT, core_NxC, core_NxD, core_NxG, core_NxK, core_NxN, core_NxR, core_NxS, core_NxT, core_NxxC, core_NxxG, core_NxxQ, core_NxxR, core_NxxS, core_NxxT, core_NxxV, core_NxxxC, core_NxxxG, core_NxxxP, core_NxxxT, core_NxxxxA, core_NxxxxC, core_NxxxxE, core_NxxxxN, core_NxxxxP, core_NxxxxQ, core_NxxxxT, core_NxxxxxC, core_NxxxxxD, core_NxxxxxG, core_NxxxxxS, core_NxxxxxT, core_P, core_PA, core_PD, core_PF, core_PH, core_PL, core_PN, core_PP, core_PR, core_PS, core_PT, core_PV, core_PW, core_PxA, core_PxC, core_PxG, core_PxL, core_PxP, core_PxR, core_PxS, core_PxT, core_PxV, core_PxW, core_PxxA, core_PxxL, core_PxxP, core_PxxR, core_PxxS, core_PxxxA, core_PxxxG, core_PxxxL, core_PxxxP, core_PxxxR, core_PxxxS, core_PxxxT, core_PxxxxA, core_PxxxxG, core_PxxxxL, core_PxxxxP, core_PxxxxR, core_PxxxxV, core_PxxxxW, core_PxxxxxA, core_PxxxxxC, core_PxxxxxD, core_PxxxxxG, core_PxxxxxL, core_PxxxxxP, core_PxxxxxR, core_PxxxxxV, core_QC, core_QR, core_QT, core_QxR, core_QxxC, core_QxxF, core_QxxR, core_QxxxL, core_QxxxR, core_QxxxxC, core_QxxxxL, core_QxxxxR, core_QxxxxxC, core_QxxxxxH, core_QxxxxxL, core_QxxxxxR, core_R, core_RA, core_RC, core_RD, core_RE, core_RF, core_RG, core_RH, core_RI, core_RK, core_RL, core_RM, core_RN, core_RP, core_RQ, core_RR, core_RS, core_RT, core_RV, core_RW, core_RxA, core_RxC, core_RxD, core_RxE, core_RxF, core_RxG, core_RxH, core_RxK, core_RxL, core_RxN, core_RxP, core_RxQ, core_RxR, core_RxS, core_RxT, core_RxV, core_RxW, core_RxY, core_RxxA, core_RxxD, core_RxxF, core_RxxG, core_RxxH, core_RxxI, core_RxxL, core_RxxM, core_RxxP, core_RxxQ, core_RxxR, core_RxxS, core_RxxT, core_RxxV, core_RxxW, core_RxxY, core_RxxxA, core_RxxxC, core_RxxxD, core_RxxxE, core_RxxxF, core_RxxxG, core_RxxxH, core_RxxxI, core_RxxxL, core_RxxxM, core_RxxxP, core_RxxxQ, core_RxxxR, core_RxxxS, core_RxxxT, core_RxxxV, core_RxxxW, core_RxxxxA, core_RxxxxC, core_RxxxxD, core_RxxxxE, core_RxxxxF, core_RxxxxG, core_RxxxxH, core_RxxxxK, core_RxxxxL, core_RxxxxM, core_RxxxxP, core_RxxxxQ, core_RxxxxR, core_RxxxxS, core_RxxxxT, core_RxxxxV, core_RxxxxW, core_RxxxxxA, core_RxxxxxC, core_RxxxxxD, core_RxxxxxE, core_RxxxxxF, core_RxxxxxG, core_RxxxxxH, core_RxxxxxI, core_RxxxxxL, core_RxxxxxM, core_RxxxxxN, core_RxxxxxP, core_RxxxxxQ, core_RxxxxxR, core_RxxxxxS, core_RxxxxxT, core_RxxxxxV, core_RxxxxxW, core_RxxxxxY, core_SA, core_SC, core_SD, core_SG, core_SK, core_SN, core_SP, core_SR, core_SS, core_ST, core_SxA, core_SxC, core_SxD, core_SxG, core_SxK, core_SxL, core_SxM, core_SxN, core_SxP, core_SxR, core_SxS, core_SxW, core_SxxC, core_SxxD, core_SxxK, core_SxxL, core_SxxP, core_SxxR, core_SxxT, core_SxxxA, core_SxxxC, core_SxxxD, core_SxxxL, core_SxxxN, core_SxxxR, core_SxxxT, core_SxxxW, core_SxxxxA, core_SxxxxC, core_SxxxxK, core_SxxxxL, core_SxxxxN, core_SxxxxP, core_SxxxxR, core_SxxxxT, core_SxxxxW, core_SxxxxxA, core_SxxxxxC, core_SxxxxxL, core_SxxxxxN, core_SxxxxxR, core_SxxxxxT, core_SxxxxxW, core_T, core_TC, core_TD, core_TE, core_TG, core_TH, core_TI, core_TK, core_TN, core_TQ, core_TS, core_TT, core_TxA, core_TxC, core_TxD, core_TxE, core_TxG, core_TxI, core_TxK, core_TxN, core_TxR, core_TxS, core_TxT, core_TxV, core_TxY, core_TxxA, core_TxxC, core_TxxD, core_TxxF, core_TxxG, core_TxxI, core_TxxK, core_TxxN, core_TxxR, core_TxxS, core_TxxT, core_TxxV, core_TxxxA, core_TxxxC, core_TxxxD, core_TxxxE, core_TxxxF, core_TxxxG, core_TxxxI, core_TxxxK, core_TxxxN, core_TxxxQ, core_TxxxR, core_TxxxS, core_TxxxT, core_TxxxxC, core_TxxxxD, core_TxxxxE, core_TxxxxG, core_TxxxxI, core_TxxxxK, core_TxxxxN, core_TxxxxQ, core_TxxxxR, core_TxxxxS, core_TxxxxT, core_TxxxxxC, core_TxxxxxE, core_TxxxxxG, core_TxxxxxK, core_TxxxxxN, core_TxxxxxR, core_TxxxxxS, core_TxxxxxT, core_V, core_VC, core_VF, core_VG, core_VI, core_VL, core_VP, core_VR, core_VT, core_VW, core_VxA, core_VxC, core_VxG, core_VxK, core_VxL, core_VxP, core_VxR, core_VxV, core_VxxG, core_VxxH, core_VxxL, core_VxxR, core_VxxT, core_VxxV, core_VxxxA, core_VxxxG, core_VxxxL, core_VxxxP, core_VxxxR, core_VxxxxL, core_VxxxxP, core_VxxxxR, core_VxxxxV, core_VxxxxxA, core_VxxxxxC, core_VxxxxxL, core_VxxxxxP, core_VxxxxxR, core_VxxxxxS, core_VxxxxxV, core_W, core_WA, core_WP, core_WR, core_WxA, core_WxP, core_WxR, core_WxxP, core_WxxR, core_WxxS, core_WxxxA, core_WxxxP, core_WxxxR, core_WxxxxA, core_WxxxxR, core_WxxxxxA, core_WxxxxxL, core_WxxxxxR, core_YC, core_YL, core_YR, core_YT, core_YxC, core_YxR, core_YxxC, core_YxxL, core_YxxxL, core_YxxxxT, core_YxxxxxR, core_length, core_pI, frac_aliphatic, frac_aromatic, frac_charge, frac_neg, frac_polar, frac_pos, leader_length, neg, peptide_length, polar, pos'

    svm_headers = svm_headers.split(',')
    features_headers = ["Accession_id", "Genus/Species/Code", "Leader", "Core", "Start", "End", "Total Score", "Valid Precursor" ] + svm_headers
    features_csv_file = open(dir_prefix + "temp_features.csv", 'w')
    svm_csv_file = open("ripp_modules/{}/svm/fitting_set.csv".format(peptide_type), 'w')
    features_writer = csv.writer(features_csv_file)
    svm_writer = csv.writer(svm_csv_file)
    features_writer.writerow(features_headers)
    svm_writer.writerow(svm_headers)
    
        
class Ripp(VirtualRipp):
    def __init__(self, 
                 start, 
                 end, 
                 sequence,
                 upstream_sequence,
                 pfam_2_coords):
        super(Ripp, self).__init__(start, 
                                     end, 
                                     sequence,
                                     upstream_sequence,
                                     pfam_2_coords)
        self.peptide_type = peptide_type
        self.score = 0
        self.set_split()
#        self.set_monoisotopic_mass()
        self.csv_columns = [self.leader, self.core, self.start, self.end]
        self.CUTOFF = CUTOFF
    
    def get_fimo_score(self, fimo_file):
        fimo_output = self.run_fimo_simple(fimo_file)
        fimo_motifs = [int(line.partition("\t")[0]) for line in fimo_output.split("\n") if "\t" in line and line.partition("\t")[0].isdigit()]
        fimo_stops = []
        for line in fimo_output.split("\n"):
            line = line.split("\t")
            if len(line) > 3 and line[3].isdigit():
                fimo_stops.append(int(line[3]))
        fimo_scores = {int(line.split("\t")[0]): float(line.split("\t")[6]) for line in fimo_output.split("\n") if "\t" in line and line.partition("\t")[0].isdigit()}
        return fimo_motifs, fimo_scores, fimo_stops
    
    def set_split(self):
        """Try to identify cleavage site using regular expressions"""
        
        f_scores = self.get_fimo_score("ripp_modules/lanthi_i/lanthi_i_leader_fimo.txt")  
        l_stop = sorted([0] + f_scores[2])[-1]
        if len(f_scores[0]) > 0:
            self.score += 10
        regex1 = re.search('(.[S|T].{2,7}C)', self.sequence[l_stop:])
        regex2 = re.search('(G[G|A])', self.sequence[l_stop:])
    
        if regex1 and regex2:
            if regex2.end() < regex1.start() and regex2.end()+l_stop > 10:
                self.split_index = regex2.end()+l_stop
            elif regex1.start() + l_stop > 10:
                self.split_index = regex1.start()+l_stop
            else:
                self.split_index =  (int)((len(self.sequence)-l_stop)/2)+l_stop
        elif regex1 and not regex2:
            if regex1.start() + l_stop > 10:
                self.split_index = regex1.start()+l_stop
            else:
                self.split_index = (int)((len(self.sequence)-l_stop)/2)+l_stop
        elif regex2 and not regex1:
            if regex2.end() + l_stop > 10:
                self.split_index = regex2.end()+l_stop
            else:
                self.split_index = (int)((len(self.sequence)-l_stop)/2)+l_stop
        else:
            self.split_index = (int)((len(self.sequence)-l_stop)/2)+l_stop
        self.leader = self.sequence[:self.split_index]
        self.core = self.sequence[self.split_index:]
    
    def set_score(self, pfam_hmm, cust_hmm):
        scoring_csv_columns = []
        scoring_column_names = []
        fifth_size = len(self.sequence) // 5
        for i in range(5):
            scoring_csv_columns.append(self.sequence[i*fifth_size:(i+1)*fifth_size].count('C') / len(self.sequence))
            scoring_column_names.append("{}/5_frac_c".format(i+1))
        for i in [0, 1, 3, 4]:
            count = 0
            for c in "STC":
                count += self.sequence[i*fifth_size:(i+1)*fifth_size].count(c)
            scoring_csv_columns.append(count / len(self.sequence))
            scoring_column_names.append("{}/5_frac_s_t_c".format(i+1))
            
        for i in [2, 3, 4]:
            count = 0
            for c in "ST":
                count += self.sequence[i*fifth_size:(i+1)*fifth_size].count(c)
            scoring_csv_columns.append(count / len(self.sequence))
            scoring_column_names.append("{}/5_frac_s_t".format(i+1))
           
        scoring_csv_columns.append((self.core.count("S") + self.core.count("T")) / len(self.sequence))
        scoring_column_names.append("S+T")
        
        scoring_csv_columns.append((self.core.count("S") + self.core.count("T") + self.core.count("C"))/ len(self.sequence))
        scoring_column_names.append("S+T+C")
        
        aliphatic_core_count = self.core.count("A") + self.core.count("I") + self.core.count("M") + self.core.count("L") + self.core.count("V")
        scoring_csv_columns.append(aliphatic_core_count)
        scoring_column_names.append("aliphatic")
        scoring_csv_columns.append(aliphatic_core_count / len(self.core))
        scoring_column_names.append("frac_aliphatic")
        
        aromatic_core_count = self.core.count("F") + self.core.count("W") + self.core.count("Y")
        scoring_csv_columns.append(aromatic_core_count)
        scoring_column_names.append("aromatic")
        scoring_csv_columns.append(aromatic_core_count  / len(self.core))
        scoring_column_names.append("frac_aromatic")
        
        charged_count = self.core.count("D") + self.core.count("E") + self.core.count("K") + self.core.count("R")
        scoring_csv_columns.append(charged_count)
        scoring_column_names.append("charged")
        scoring_csv_columns.append(charged_count  / len(self.core))
        scoring_column_names.append("frac_charge")
        
        polar_count = self.core.count("C") + self.core.count("S") + self.core.count("T") + self.core.count("H") + self.core.count("N") + self.core.count("Q")
        scoring_csv_columns.append(polar_count)
        scoring_column_names.append("polar")
        scoring_csv_columns.append(polar_count  / len(self.core))
        scoring_column_names.append("frac_polar")
        
        pos = self.core.count("R") + self.core.count("K")
        scoring_csv_columns.append(pos)
        scoring_column_names.append("pos")
        scoring_csv_columns.append(pos / len(self.core))
        scoring_column_names.append("frac_pos")
        
        
        neg = self.core.count("D") + self.core.count("E")
        scoring_csv_columns.append(neg)
        scoring_column_names.append("neg")
        scoring_csv_columns.append(neg / len(self.core))
        scoring_column_names.append("frac_neg")
        
        scoring_csv_columns.append(len(self.leader))
        scoring_column_names.append("leader_length")
        
        scoring_csv_columns.append(len(self.sequence))
        scoring_column_names.append("peptide_length")
        
        scoring_csv_columns.append(len(self.core))
        scoring_column_names.append("core_length")
        
        core_seqs = ["core_CxxxC","core_C","core_TC","core_CxxxxC","core_TxxxxC","core_TxxxC","core_CxT","core_TxxC","core_CT","core_CxxT","core_SxC","core_TxxxxxC","core_TxxxxxT","core_DxC","core_SC","core_DxxxxxC","core_CxxxxxT","core_T","core_CxxxxxN","core_R","core_GxT","core_TxxxT","core_TxxT","core_SxxxxxC","core_TxD","core_CxxG","core_CxxxT","core_CK","core_CxxxxT","core_GxxC","core_CxS","core_DxxxxT","core_TxxxxT","core_CN","core_CxxxxN","core_TxG","core_NxxC","core_CxxxxxC","core_NxxxC","core_CxxC","core_CxC","core_CC","core_TxN","core_SxxxC","core_SxxC","core_CxxS","core_GxxN","core_DN","core_NT","core_CxxN","core_DxxC","core_TxxxN","core_CG","core_GC","core_CxxxxG","core_TxxG","core_L","core_GxxxT","core_NC","core_GxxxxC","core_TxT","core_DxxxxS","core_SxN","core_GxxxxxC","core_SxxxxT","core_SxxxxC","core_AR","core_CxN","core_RxR","core_KT","core_RxxR","core_RxxxR","core_RR","core_RS","core_TD","core_TxC","core_TxxN","core_TxxS","core_RxxA","core_TxV","core_RxxxxxR","core_AxxxxR","core_RxxxxR","core_IxxxxxC","core_AxR","core_DxxxC","core_TxxxxxG","core_GxxxC","core_RxxxA","core_CxxD","core_CxxxxxI","core_AxxR","core_SxxxR","core_AxxxxxR","core_RxxG","core_TT","core_RxG","core_SxxR","core_AxxxR","core_RxA","core_CxK","core_PxxxR","core_PxR","core_GxxR","core_PxxR","core_GxR","core_RP","core_P","core_KxxC","core_NxxxxxT","core_RxxxS","core_RxxxxxA","core_RA","core_RxxP","core_TN","core_RxP","core_N","core_RxxxxG","core_RxxS","core_GxC","core_CxxxxxS","core_RxxxxA","core_NxC","core_GR","core_RxxxxxP","core_CS","core_SxR","core_SxxxxR","core_SxxxT","core_TxxxxG","core_CxxxG","core_RxS","core_TxxxxN","core_GxxxR","core_PxxxxR","core_RxxxxS","core_SR","core_GxxxxR","core_RG","core_NxxxxC","core_RxxxP","core_NxxxxxG","core_RxxxG","core_RxxxxxS","core_CxxxD","core_KxC","core_PxxxxxR","core_RxxxxxG","core_SxxxxxR","core_GxxxxT","core_RxL","core_TI","core_RxxxxP","core_RL","core_LR","core_GxxxxxT","core_TG","core_GS","core_NxxT","core_KxxxxxC","core_TxxxxxN","core_LxxR","core_VR","core_PR","core_CxxxS","core_LxxxxxR","core_NxG","core_GxxxxxR","core_LxxxxR","core_TxxxxxS","core_GxxxxK","core_CxxxN","core_NxxxT","core_CxxxH","core_DxxxxxT","core_CxxxK","core_GxxxxxN","core_LxxxR","core_FxxxC","core_RxxV","core_NS","core_CR","core_MW","core_RxxxxL","core_SxxxxxT","core_GxN","core_TK","core_RxxL","core_W","core_AxP","core_LxR","core_DxxxxC","core_GxxxxN","core_KxxxxT","core_CxxxxS","core_KxT","core_DxxxS","core_TS","core_NxxxxxC","core_FxC","core_QxxxxC","core_RxV","core_RxxxxV","core_GxxT","core_RxxxL","core_VxxxxR","core_H","core_KxxxT","core_KxxT","core_CxD","core_VxR","core_GxxxS","core_LxxxxC","core_NxxxxT","core_PxxxxxC","core_VxxR","core_IxC","core_RxxxV","core_CxxxxxK","core_PA","core_TxK","core_TxE","core_GT","core_VxxxR","core_AC","core_A","core_DT","core_HR","core_CxG","core_FT","core_NxT","core_TxxxxxK","core_ExxxxT","core_DxT","core_LL","core_TxxR","core_NG","core_DxE","core_TxI","core_TxxxI","core_RxxT","core_RC","core_TxxxD","core_TxxxxK","core_IT","core_VT","core_YT","core_SP","core_RxxxxxL","core_PxxxxxP","core_RxxxxxV","core_SxxxxK","core_NxxxxP","core_RV","core_CxxK","core_SD","core_DC","core_PxxxP","core_LxxxL","core_IC","core_LP","core_GxxxN","core_AxL","core_VxxxxxC","core_CP","core_CxxR","core_PxxxxA","core_PxS","core_TxxxxxR","core_SxxxxN","core_AxxT","core_IxxxT","core_DxxT","core_TxA","core_PxL","core_CxxxR","core_HxR","core_VxxxxxR","core_AxxxxxN","core_DxxS","core_TxxxS","core_NxS","core_PxxA","core_RxxxxT","core_QR","core_PxxxxxA","core_QxxxxxC","core_LxxxxxL","core_AxxxxA","core_CI","core_PxxP","core_NxK","core_RxxxH","core_CxxxxR","core_DxxxxxD","core_D","core_AxxxA","core_TxxxR","core_PxC","core_TxxI","core_AxxxP","core_RxxH","core_SxxxD","core_KxxxxxG","core_PxxxxP","core_DD","core_AxxxxxP","core_LxxL","core_YxC","core_AxxA","core_AG","core_PS","core_LxxxxL","core_TxxxK","core_RxW","core_RxxxW","core_GxxxxxS","core_IxxT","core_LxxxxxA","core_ExxC","core_NxxxG","core_RT","core_TxxxxS","core_RxxxxxH","core_DF","core_HxxxxxK","core_WR","core_F","core_DxxxG","core_RxxxQ","core_SxxxL","core_RxxxT","core_RxxxxxQ","core_FxR","core_LxxxF","core_CxxxxxQ","core_VxxT","core_AxxxxxL","core_RxE","core_AxxxxP","core_RW","core_GP","core_WxxxxxR","core_QxxC","core_RxH","core_WxxxxR","core_YxxC","core_QxR","core_AP","core_CxxxxxR","core_CD","core_PxA","core_RxxxxH","core_PxxxxG","core_FxL","core_CA","core_HxxxxxR","core_PxxxxxV","core_FxxxR","core_V","core_VP","core_TxxxE","core_HxxxC","core_WxxR","core_AxxxxxA","core_RQ","core_RxxW","core_TxxxA","core_FxxL","core_VxC","core_AV","core_KxxxxxS","core_IxxxS","core_LxL","core_TxY","core_AxxxxL","core_DxxxT","core_VxK","core_DxxxK","core_CxxV","core_GxS","core_SxxxxxL","core_GxxxxP","core_CxxY","core_PxG","core_NxxxxxS","core_KxxD","core_CxxxxI","core_TxxxxR","core_QxxxxxR","core_ExxxxxR","core_ExxxxxC","core_IxxC","core_CxR","core_TxxV","core_WxR","core_SxxxN","core_M","core_LxxxxA","core_IxR","core_FxxxxxR","core_AL","core_GxxxxL","core_ExxT","core_DxxxxxR","core_TxxxxD","core_ED","core_LxG","core_CxE","core_IR","core_SxxxxL","core_RxxxxF","core_GxP","core_HxxxR","core_GxxxL","core_SxxP","core_DxxG","core_KxxxC","core_LxP","core_HxxxxR","core_ExxR","core_GxxI","core_LxS","core_FxxxxR","core_PxxxS","core_FL","core_AxxxxC","core_PxxL","core_ExR","core_VxA","core_QxxR","core_RxxxxxT","core_IxxxxC","core_SxP","core_GxA","core_PxxxxxL","core_PxxxxL","core_LxV","core_KxxxxD","core_DS","core_RxxxxxC","core_MxxT","core_GW","core_TQ","core_LS","core_AxxxL","core_LxxxxI","core_VL","core_HxxA","core_ER","core_DxxxR","core_LxxF","core_AxxxxxT","core_DP","core_VxxxxL","core_PxxxA","core_RxxxE","core_FP","core_KC","core_PH","core_TxxxG","core_LxxxxS","core_CE","core_FxxxxxL","core_VC","core_RxxxxD","core_VxG","core_HxxxxxP","core_HxxxP","core_IxxxxT","core_RxC","core_WxxxR","core_AxH","core_FF","core_SxxK","core_LxxxxP","core_RI","core_QxxxxR","core_LxH","core_LH","core_YC","core_FI","core_KxxxxG","core_GxV","core_RxxxxxE","core_QxxxxxL","core_RxxxxQ","core_RD","core_VxxxL","core_AxxxH","core_RxxxI","core_NxR","core_TE","core_LxxG","core_AxV","core_GxxxxS","core_FxxxxxT","core_PxP","core_NxxS","core_LxxxxxF","core_LV","core_RxxF","core_GD","core_GxxxP","core_SxxxxxN","core_ExxxxR","core_RxK","core_RH","core_WxA","core_ST","core_DR","core_TxxxxE","core_RxxxxxW","core_PxV","core_SxxxxP","core_WP","core_RxxxxxF","core_PxxxG","core_CxxxxxG","core_AxxL","core_RxxxxC","core_PT","core_NxxxP","core_TH","core_AH","core_QT","core_SxxxW","core_VxL","core_GxxxxxL","core_YR","core_FxxxxC","core_HxxxS","core_HxxL","core_VxxxxxV","core_GxxP","core_LxxxxxQ","core_CxP","core_DxR","core_KxxxE","core_HxxC","core_FxxxL","core_IF","core_QxxxxL","core_AxC","core_RxD","core_AxN","core_VxxxxxL","core_VxxG","core_FxxR","core_SxL","core_HxxxL","core_CxxxxxF","core_GH","core_WxxxxA","core_FxxxF","core_ExxxN","core_SxK","core_HxxxxA","core_RxxxxxD","core_GxxxxA","core_WA","core_TxR","core_PxxxL","core_KxxR","core_RxxxC","core_GE","core_DxxxxN","core_PxxxxW","core_WxxxxxA","core_KxxxxxL","core_RxQ","core_GL","core_LxxP","core_PD","core_TxxxxI","core_FxxC","core_TxxxQ","core_FxxS","core_KxxA","core_LxF","core_GxxL","core_RxT","core_AxxP","core_KxxxxC","core_AxxG","core_CxxxI","core_LxxV","core_RxF","core_GxxV","core_DxN","core_WxxxP","core_SxG","core_VxxxxP","core_CW","core_SxxxxA","core_MxxxxR","core_HxxS","core_GxxxxxP","core_ExxxD","core_HA","core_LxxxxH","core_FxxxxxF","core_NxxxxQ","core_SxxT","core_CxxxxA","core_FxxxxI","core_QxxxR","core_AxT","core_RxxxxW","core_KxxxxxR","core_SxW","core_SxxxA","core_YxxL","core_NxxV","core_KxS","core_IxxxxxL","core_RxxQ","core_ExT","core_PW","core_HxxxA","core_PL","core_KxR","core_GN","core_IxxxR","core_WxP","core_FxI","core_RxxxM","core_FR","core_PV","core_IxS","core_ExxxxxT","core_KxxxxR","core_RxxxxE","core_PF","core_HL","core_FxxxxL","core_RxxxF","core_SxS","core_WxxS","core_PxxxxxD","core_AxxxN","core_VxxxP","core_KxxxxE","core_CxV","core_MxR","core_HxxG","core_RK","core_RxY","core_LF","core_RM","core_VxxV","core_GxxxA","core_RN","core_FxxxI","core_RxxI","core_VxxL","core_FxF","core_AxxxxS","core_GxL","core_FxP","core_IxxxxR","core_IxT","core_DxxN","core_AxA","core_KxxxL","core_LxxxH","core_SxxxxxW","core_YL","core_RxxY","core_PP","core_KxG","core_RxxD","core_HxxxxL","core_PxxxT","core_TxxxF","core_AxxH","core_AW","core_IxxxxF","core_IxxxxxR","core_SxxL","core_LxxH","core_WxxxA","core_MxxR","core_YxR","core_QxxxL","core_RxxxxxN","core_RxxxD","core_SS","core_RxN","core_SxxD","core_LxxS","core_VxxxxxA","core_LxN","core_TxxD","core_LxxxP","core_DxxxxY","core_RE","core_ET","core_DxxxD","core_VxxxxV","core_VI","core_LxxxxxI","core_CxxQ","core_PN","core_MC","core_HC","core_SA","core_LxxxS","core_AxxxxxW","core_LxxxxV","core_NxxQ","core_VxxxG","core_VW","core_VxxxA","core_NxD","core_IxxxxxF","core_NxxxxA","core_LA","core_VxxxxxS","core_RxxM","core_YxxxxT","core_PxxxxV","core_VxxH","core_ExxxxxN","core_DxxxxxE","core_PxxxxxG","core_FxA","core_VxV","core_SxD","core_FxxxxxK","core_SxxxxW","core_PxT","core_SK","core_RF","core_RxxxxxM","core_AxxxxW","core_NxN","core_IxxF","core_TxS","core_PxW","core_TxxA","core_DH","core_TxxK","core_GxxxxxD","core_HxA","core_LxM","core_FxxxW","core_IxL","core_NxxG","core_KxL","core_SN","core_CxxxxD","core_AxxxxM","core_DG","core_MxxxR","core_FxxxxG","core_RxxxxxY","core_DxG","core_HxG","core_LxxxxG","core_LK","core_LxxxxxY","core_MV","core_LxY","core_TxxF","core_VxP","core_NxxxxN","core_GxxxK","core_IxxxxxT","core_PxxS","core_LxxxxxH","core_YxxxxxR","core_HP","core_AM","core_MxxP","core_HxxxxxL","core_NxxxxxD","core_VG","core_GxxxxH","core_SxA","core_LxxxA","core_NxxxxE","core_SxxxxxA","core_RxxxxM","core_NR","core_TxxxxxE","core_ExxxxxH","core_WxxP","core_GxxA","core_SxM","core_RxxxxxI","core_VF","core_SG","core_FQ","core_FxxN","core_DxxxN","core_CV","core_IxxA","core_QC","core_NxxR","core_TxxxxQ","core_FA","core_FxxxV","core_LxxxxxW","core_YxxxL","core_DxxxxL","core_RxxxxK","core_QxxF","core_DxxxxR","core_FxxxxxI","core_VxxxxxP","core_AxW","core_FxxxS","core_WxxxxxL","core_DxxxxxN","core_AT","core_LxQ","core_LxxxxxV","core_LxxxxxM","core_QxxxxxH","core_LxxK"]
#        
#        
        for seq in core_seqs:
            regex = re.compile(seq[5:].replace("x", "."))
            scoring_csv_columns.append(len(regex.findall(self.core)))
            scoring_column_names.append(seq)
        core_pi = ProteinAnalysis(self.core).isoelectric_point()
        scoring_csv_columns.append(core_pi)
        scoring_column_names.append("core_pI")
        
        if core_pi < 9:
            self.score += 4
        if self.core.count("C") >= 2:
            self.score += 4
            
        regex = re.compile("KL.L.K")
        if self.leader[-2:] == "GG" and len(regex.findall(self.leader)) > 0:
            self.score += 2
        
        scoring_csv_columns = [x for _, x in sorted(zip(scoring_column_names, scoring_csv_columns))]
#        print(sorted(scoring_column_names))
        self.csv_columns += [self.score] + scoring_csv_columns
        
                
  