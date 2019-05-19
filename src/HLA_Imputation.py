#-*- coding: utf-8 -*-

import os, sys, re
from os.path import join
import pandas as pd

from src.HLA_MultipleRefs import HLA_MultipleRefs
from src.redefineBPv1BH import redefineBP
from src.Panel_subset import Panel_Subset
from src.bgl2GC_trick_bgl import Bgl2GC


########## < Core Varialbes > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]

# __overlap__ = [3000, 4000, 5000]
__overlap__ = [3000]


class HLA_Imputation(object):

    def __init__(self, __MHC__, __REFERENCE__, _out, _hg, _p_LINKAGE2BEAGLE, _p_PLINK, *args, **kwargs):


        ###### < Loading information related to Multiple markers > ######

        df_reference_bim = pd.read_csv(__REFERENCE__+'.bim', sep='\s+', names=['Chr', 'Label', 'GD', 'POS', 'A1', 'A2'])
        # print(std_MAIN_PROCESS_NAME + 'Loaded reference bim file : \n{}'.format(df_reference_bim.head()))

        df_EXON_info = pd.read_csv('data/HLA_EACH_EXON_POSITIONS_hg{}.txt'.format(_hg), header=0, sep='\s+', usecols=[0, 1, 4], index_col=[1])
        # print(std_MAIN_PROCESS_NAME + 'Loaded Exon information : \n{}'.format(df_EXON_info.loc['exon2', :]))




        ###### < 'CONVERT_IN', 'IMPUTE', 'CONVERT_OUT' with multiple markers. > ######

        print(std_MAIN_PROCESS_NAME + "'CONVERT_IN', 'IMPUTE', 'CONVERT_OUT' with multiple markers.")

        ### Main iteration
        # for _exonN_ in ['exon2', 'exon3', 'exon4']:
        for _exonN_ in ['exon2']:

            __Multiple_Refs__ = HLA_MultipleRefs(_exonN_, df_EXON_info.loc[_exonN_, :].reset_index(drop=True).set_index('HLA'),
                                                 __REFERENCE__, df_reference_bim, _out, _hg, _p_PLINK)


            for _overlap_ in __overlap__:

                print(std_MAIN_PROCESS_NAME + "exonN: {} / Overlap: {}".format(_exonN_, _overlap_))


                # Starting with `MHC`+'.QC'
                