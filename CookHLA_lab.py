#-*- coding: utf-8 -*-

import os, sys, re
from os.path import join
import argparse, textwrap
import pandas as pd

from CookHLA import CookHLA



########## < Core Varialbes > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]



def CookHLA_lab(_args):

    INPUT = _args.input
    OUT = _args.out
    REFRENCE = _args.reference
    HG = _args.hg
    GeneticMap = _args.genetic_map
    AverageErate = _args.average_erate
    ANSWER = _args.answer
    ANSWER2 = _args.answer2
    HapMap_Map = _args.hapmap_map

    JAVA_MEM = _args.java_memory


    __accuracies__ = []


    # Output directory processing
    OUTPUT_dir = os.path.dirname(OUT)
    OUTPUT_prefix = os.path.basename(OUT)


    ###### < 2_Plain > ######
    print("<Imputation : _2_Plain>")

    OUT_2_Plain = join(OUTPUT_dir, '_2_Plain', OUTPUT_prefix+'.Plain')

    [t_HLA_Imptation_out, t_accuracy] = \
        CookHLA(INPUT, OUT_2_Plain, REFRENCE, _answer=ANSWER, _java_memory=JAVA_MEM)

    print("HLA_IMPUTATION_OUT : {}".format(t_HLA_Imptation_out))
    print("Accuracy : {}".format(t_accuracy))



    ###### _3_MM ######
    print("<Imputation : _3_MM>")

    OUT_3_MM = join(OUTPUT_dir, '_3_MM', OUTPUT_prefix+'.MM')




    ###### _4_AGM_HapMap_Map ######
    print("<Imputation : _4_AGM_HapMap_Map>")

    OUT_4_AGM_HapMap_Map = join(OUTPUT_dir, '_4_AGM_HapMap_Map', OUTPUT_prefix+'.AGM_HapMap_Map')



    ###### _5_AGM ######
    print("<Imputation : _5_AGM>")

    OUT_5_AGM = join(OUTPUT_dir, '_5_AGM', OUTPUT_prefix+'.AGM')



    ###### _6_AGM_MM ######
    print("<Imputation : _6_MM_AGM>")

    OUT_6_MM_AGM = join(OUTPUT_dir, '_6_MM_AGM', OUTPUT_prefix+'.MM.AGM')







    return 0



if __name__ == "__main__":

    ########## < Main parser > ##########

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    ###########################################################################################

        CookHLA_lab.py

        (Created by Wanson Choi.)
        
        To automate generating accuracy comparison table, this module has been introduced.
        Each cookhla imputation will be implemented as execution of functions and their accuracy 
        outputs will be summarized by pandas.
        
        Each imputations are
        
        _2_Plain
        _3_MM
        _4_AGM_HapMap_Map
        _5_AGM
        _6_MM_AGM
        
        To execute all above imputations, next arguments must be given to 'CookHLA_lab.py'
        
        1. --input/-i
        2. --out/-o
        3. --reference/-ref
        4. --genetic-map/-gm
        5. --average-erate/-ae
        6. --answer + --answer2(Fam column fixed)
        7. --hapmap-map/-hm



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


    # For Testing
    parser.add_argument("--genetic-map", "-gm", help="\nGenetic Map file.\n\n")
    parser.add_argument("--average-erate", "-ae", help="\nAverate error rate file.\n\n")
    parser.add_argument("--use-multiple-markers", "-ml", help="\nUsing multiple markers.\n\n", action='store_true')

    # parser.add_argument("--prephasing", "-pr", help="\nUtilizing prephasing strategy.\n\n", action='store_true')

    parser.add_argument("--answer", "-an", help="\nAnswer file to calculate imputation accuracy.\n\n")
    parser.add_argument("--answer2", "-an2", help="\nAnswer file to calculate imputation accuracy(Fam colum fixed).\n\n")

    # parser.add_argument("--multiprocess", "-mp", help="\nSetting parallel multiprocessing.\n\n", type=int, choices=[2,3,4,5,6,7,8,9], nargs='?', default=1, const=3)

    parser.add_argument("--java-memory", "-mem", help="\nMemory requried for beagle(ex. 12g).\n\n", default="2g")

    # parser.add_argument("--prephased", "-ph",
    #                     help="\n(For Testing Purpose) Passing prephased result manually to control the error rate of prephasing. "
    #                          "If given, Imputation will be done based on this phased file.\n\n")

    parser.add_argument("--hapmap-map", "-hm",
                        help="\n(For Testing Purpose) Hapmap Map(Adaptive Genetic Map).\n\n")





    ##### < for Testing > #####

    args = parser.parse_args(["--input", "data/Target/HM_CEU.FOUNDERS.filt",
                              "--out", "tests/CookHLA_lab/20190731/HM_CEU_T1DGC_REF",
                              "-ref", "data/HLA_PANEL/T1DGC/T1DGC_REF",
                              "-gm", "data/HLA_PANEL/Genetic_map/CEU_T1DGC.mach_step.avg.clpsB",
                              "-ae", "data/HLA_PANEL/Genetic_map/CEU_T1DGC.aver.erate",
                              "-an", "tests/HM_CEU_REF.bgl.phased.alleles.answer",
                              "-hm", "data/HapMap_Map.txt"])

    ## Only MM.
    # args = parser.parse_args(["--input", "data/Target/HM_CEU.FOUNDERS.filt",
    #                           "--out", "tests/_3_CookHLA/20190708_MM/_3_HM_CEU_T1DGC_REF",
    #                           "-ref", "data/HLA_PANEL/T1DGC/T1DGC_REF",
    #                           "-ml",
    #                           "-an", "tests/HM_CEU_REF.bgl.phased.FIDadj.alleles.answer",
    #                           "-mem", "4g"])




    ##### < for Publish > #####
    # args = parser.parse_args()
    print(args)

    CookHLA_lab(args)