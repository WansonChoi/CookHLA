#-*- coding: utf-8 -*-

import os, sys, re
from os.path import join
from shutil import which
import argparse, textwrap


std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]



def CookQC():


    return 0



if __name__ == "__main__":

    ########## < Main parser > ##########

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    ###########################################################################################

        CookQC.py

        (Created by W. Choi.)



    ###########################################################################################
                                     '''),
                                     add_help=False
                                     )

    ### Common arguments to share over the modules.

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="Show this help message and exit\n\n", action='help')

    parser.add_argument("--input", "-i", help="\nCommon prefix of input files.\n\n", required=True)
    parser.add_argument("--reference", "-ref", help="\nPrefix of Reference files.\n\n", required=True)
    parser.add_argument("--out", "-o", help="\nOutput file name prefix\n\n", required=True)
    parser.add_argument("-hg", help="\nHuman Genome version(ex. 18, 19, 38)\n\n", choices=["18", "19", "38"],
                        metavar="HG", default='18')







    ##### < for Testing > #####

    # args = parser.parse_args(["--input", "data/Target/HM_CEU.FOUNDERS.filt",
    #                           "--out", "tests/_3_CookHLA/20190605_onlyAGM/_3_HM_CEU_T1DGC_REF",
    #                           "-ref", "data/HLA_PANEL/T1DGC/T1DGC_REF",
    #                           "-gm", "data/HLA_PANEL/Genetic_map/CEU_T1DGC.mach_step.avg.clpsB",
    #                           "-ae", "data/HLA_PANEL/Genetic_map/CEU_T1DGC.aver.erate",
    #                           "-an", "tests/HM_CEU_REF.bgl.phased.alleles.answer"])


    ##### < for Publish > #####
    args = parser.parse_args()
    print(args)

    CookQC()
