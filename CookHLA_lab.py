#-*- coding: utf-8 -*-

import os, sys, re
from os.path import join
import argparse, textwrap
import pandas as pd
from time import time

from CookHLA import CookHLA



########## < Core Varialbes > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]



def CookHLA_lab(_args, _control_flags=(1,1,1,1,1)):

    # Arguments assignment.

    INPUT = _args.input
    OUT = _args.out
    REFRENCE = _args.reference
    HG = _args.hg
    GeneticMap = _args.genetic_map
    AverageErate = _args.average_erate
    ANSWER = _args.answer
    ANSWER2 = _args.answer2 if _args.answer2 else ANSWER
    HapMap_Map = _args.hapmap_map

    JAVA_MEM = _args.java_memory

    PREPHASING = _args.prephasing

    # Beagle5.1.
    __OVERLAP__ = _args.overlap




    # Accuracy outputs.
    __accuracies__ = {2: None,
                      3: None,
                      4: None,
                      5: None,
                      6: None}


    # Output directory processing
    OUTPUT_dir = os.path.dirname(OUT)
    OUTPUT_prefix = os.path.basename(OUT)


    # Control Flags
    if len(_control_flags) != 5:
        print(std_ERROR_MAIN_PROCESS_NAME + "Wrong '_control_flags' value. Please check it again.")
        sys.exit()

    if _control_flags == (1,1,1,1,1): # default value

        F_2_Plain = 1
        F_3_MM = 1
        F_4_AGM_HapMap_Map = 1
        F_5_AGM = 1
        F_6_MM_AGM = 1

    else:
        [F_2_Plain, F_3_MM, F_4_AGM_HapMap_Map, F_5_AGM, F_6_MM_AGM] = _control_flags



    ### Main Imputations.

    if F_2_Plain:

        ###### < 2_Plain > ######
        print(std_MAIN_PROCESS_NAME + "<Imputation : _2_Plain>")

        OUT_2_Plain = join(OUTPUT_dir, '_2_Plain', OUTPUT_prefix+'.Plain')

        time_start_2_Plain = time()

        [t_HLA_Imptation_out, t_accuracy] = \
            CookHLA(INPUT, OUT_2_Plain, REFRENCE, _answer=ANSWER, _java_memory=JAVA_MEM)

        time_end_2_Plain = time()
        print("Implementation time of _2_Plain : {}(min)".format((time_end_2_Plain - time_start_2_Plain)/60))

        # print("HLA_IMPUTATION_OUT : {}".format(t_HLA_Imptation_out))
        # print("Accuracy : {}".format(t_accuracy))

        __accuracies__[2] = t_accuracy



    if F_3_MM:

        ###### _3_MM ######
        print(std_MAIN_PROCESS_NAME + "<Imputation : _3_MM>")

        OUT_3_MM = join(OUTPUT_dir, '_3_MM', OUTPUT_prefix+'.MM')

        time_start_3_MM = time()

        [t_HLA_Imptation_out, t_accuracy] = \
            CookHLA(INPUT, OUT_3_MM, REFRENCE,
                    __use_Multiple_Markers=True, _MultP=_args.multiprocess,
                    _answer=ANSWER2, _java_memory=JAVA_MEM, f_prephasing=PREPHASING,
                    __overlap__=__OVERLAP__)

        time_end_3_MM = time()
        print("Implementation time of _3_MM : {}(min)".format((time_end_3_MM - time_start_3_MM)/60))


        __accuracies__[3] = t_accuracy



    if F_4_AGM_HapMap_Map:

        ###### _4_AGM_HapMap_Map ######
        print(std_MAIN_PROCESS_NAME + "<Imputation : _4_AGM_HapMap_Map>")

        OUT_4_AGM_HapMap_Map = join(OUTPUT_dir, '_4_AGM_HapMap_Map', OUTPUT_prefix+'.AGM_HapMap_Map')

        time_start_4_HapMap_Map = time()

        [t_HLA_Imptation_out, t_accuracy] = \
            CookHLA(INPUT, OUT_4_AGM_HapMap_Map, REFRENCE, _HapMap_Map=HapMap_Map,
                    _answer=ANSWER, _java_memory=JAVA_MEM)

        time_end_4_HapMap_Map = time()
        print("Implementation time of _4_HapMap_Map : {}(min)".format((time_end_4_HapMap_Map - time_start_4_HapMap_Map)/60))


        __accuracies__[4] = t_accuracy



    if F_5_AGM:

        ###### _5_AGM ######
        print(std_MAIN_PROCESS_NAME + "<Imputation : _5_AGM>")

        OUT_5_AGM = join(OUTPUT_dir, '_5_AGM', OUTPUT_prefix+'.AGM')

        time_start_5_AGM = time()

        [t_HLA_Imptation_out, t_accuracy] = \
            CookHLA(INPUT, OUT_5_AGM, REFRENCE, _AdaptiveGeneticMap=GeneticMap, _Average_Erate=AverageErate,
                    _answer=ANSWER, _java_memory=JAVA_MEM)

        time_end_5_AGM = time()
        print("Implementation time of _5_AGM : {}(min)".format((time_end_5_AGM - time_start_5_AGM)/60))


        __accuracies__[5] = t_accuracy



    if F_6_MM_AGM:

        ###### _6_MM_AGM ######
        print(std_MAIN_PROCESS_NAME + "<Imputation : _6_MM_AGM>")

        OUT_6_MM_AGM = join(OUTPUT_dir, '_6_MM_AGM', OUTPUT_prefix+'.MM.AGM')

        time_start_6_MM_AGM = time()

        [t_HLA_Imptation_out, t_accuracy] = \
            CookHLA(INPUT, OUT_6_MM_AGM, REFRENCE,
                    __use_Multiple_Markers=True, _MultP=_args.multiprocess,
                    _AdaptiveGeneticMap=GeneticMap, _Average_Erate=AverageErate,
                    _answer=ANSWER2, _java_memory=JAVA_MEM, f_prephasing=PREPHASING,
                    __overlap__=__OVERLAP__)

        time_end_6_MM_AGM = time()

        print("Implementation time of _6_MM_AGM : {}(min)".format((time_end_6_MM_AGM - time_start_6_MM_AGM)/60))


        __accuracies__[6] = t_accuracy


    # [Temporary Hard-coding]
    # __accuracies__[2] = 'tests/accuracy_ex_data/_2_HM_CEU_T1DGC_REF.Plain.overlap3000.MHC.HLA_IMPUTATION_OUT.alleles.accuracy'
    # __accuracies__[3] = 'tests/accuracy_ex_data/_3_HM_CEU_T1DGC_REF.MM.MHC.HLA_IMPUTATION_OUT.alleles.accuracy'
    # __accuracies__[4] = 'tests/accuracy_ex_data/_4_HM_CEU_T1DGC_REF.AGM_HapmapMap.MHC.HLA_IMPUTATION_OUT.alleles.accuracy'
    # __accuracies__[5] = 'tests/accuracy_ex_data/_5_HM_CEU_T1DGC_REF.AGM.MHC.HLA_IMPUTATION_OUT.alleles.accuracy'
    # __accuracies__[6] = 'tests/accuracy_ex_data/_6_HM_CEU_T1DGC_REF.MM.AGM.noprephasing.MHC.HLA_IMPUTATION_OUT.alleles.accuracy'
    # print(__accuracies__)


    __RETURN__ = CollectTable(__accuracies__)
    print(std_MAIN_PROCESS_NAME + "Accuracy Table:\n{}\n".format(__RETURN__))
    __RETURN__.to_csv(OUT+".ACCURACY_TABLE.txt", sep='\t', header=True, index=True)

    return __RETURN__



def CollectTable(__accuracies__):

    l_df = []
    l_label = []

    for i in range(2, 7):

        if __accuracies__[i] and os.path.exists(__accuracies__[i]):

            # measureAcc_v2
            df_temp = pd.read_csv(__accuracies__[i], sep='\s+', header=None, index_col=[0,1], names=['HLA', '4D', 'acc'])

            # measureAcc_v3.5
            # df_temp = pd.read_csv(__accuracies__[i], sep='\s+', header=None, index_col=0, names=['HLA', 'acc'])
            # print(df_temp)
            l_df.append(df_temp)
            l_label.append(Int2Label(i))


    df_acc = pd.concat(l_df, axis=1).applymap(lambda x : None if x == -1 else x).dropna()
    df_acc.columns = l_label
    # print(df_acc)

    ## Acquiring average accuracy.
    sr_mean = df_acc.mean(axis=0)
    # print(sr_mean)

    df_acc2 = df_acc.append(sr_mean, ignore_index=True)
    # print(df_acc2)

    ## New index
    idx_df_acc2 = df_acc.index.to_frame().loc[:, 'HLA'].tolist()
    idx_df_acc2.append('avg')
    df_acc2.index = idx_df_acc2
    df_acc2.index.name = 'HLA'
    # print(df_acc2)

    return df_acc2



def Int2Label(_x):

    if _x == 2:
        return '_2_Plain'
    elif _x == 3:
        return '_3_MM'
    elif _x == 4:
        return '_4_AGM_HapMap_Map'
    elif _x == 5:
        return '_5_AGM'
    elif _x == 6:
        return '_6_MM_AGM'
    else:
        return '-1'



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
    parser.add_argument("--genetic-map", "-gm", help="\nGenetic Map file.\n\n", required=True)
    parser.add_argument("--average-erate", "-ae", help="\nAverate error rate file.\n\n", required=True)
    # parser.add_argument("--use-multiple-markers", "-ml", help="\nUsing multiple markers.\n\n", action='store_true')

    # parser.add_argument("--prephasing", "-pr", help="\nUtilizing prephasing strategy.\n\n", action='store_true')

    parser.add_argument("--answer", "-an", help="\nAnswer file to calculate imputation accuracy.\n\n", required=True)
    parser.add_argument("--answer2", "-an2", help="\nAnswer file to calculate imputation accuracy(Fam colum fixed).\n\n")

    parser.add_argument("--multiprocess", "-mp", help="\nSetting parallel multiprocessing.\n\n", type=int, choices=[2,3,4,5,6,7,8,9], nargs='?', default=1, const=3)

    parser.add_argument("--java-memory", "-mem", help="\nMemory requried for beagle(ex. 12g).\n\n", default="2g")

    # parser.add_argument("--prephased", "-ph",
    #                     help="\n(For Testing Purpose) Passing prephased result manually to control the error rate of prephasing. "
    #                          "If given, Imputation will be done based on this phased file.\n\n")

    parser.add_argument("--hapmap-map", "-hm",
                        help="\n(For Testing Purpose) Hapmap Map(Adaptive Genetic Map).\n\n", required=True)

    parser.add_argument("--control-flags", "-cf",
                        help="\nBoolean sequence to nominate which imputations are to be done.\n\n", nargs=5, default=(1,1,1,1,1), type=int)

    parser.add_argument("--prephasing", "-pr", help="\nUtilizing prephasing strategy.\n\n", action='store_true')

    # Beagle5.1.
    parser.add_argument("--overlap", "-ol",
                        help="\n3 Overlap values(cM) for Beagle 5.1 implementation.\n\n", nargs=3, default=(4,8,12), type=int)




    ##### < for Testing > #####

    # args = parser.parse_args(["--input", "data/Target/HM_CEU.FOUNDERS.filt",
    #                           "--out", "tests/CookHLA_lab/20190731/HM_CEU_T1DGC_REF",
    #                           "-ref", "data/HLA_PANEL/T1DGC/T1DGC_REF",
    #                           "-gm", "data/HLA_PANEL/Genetic_map/CEU_T1DGC.mach_step.avg.clpsB",
    #                           "-ae", "data/HLA_PANEL/Genetic_map/CEU_T1DGC.aver.erate",
    #                           "-an", "data/answer/HM_CEU_REF.bgl.phased.alleles.answer",
    #                           "-an2", "data/answer/HM_CEU_REF.bgl.phased.FIDadj.alleles.answer",
    #                           "-hm", "data/HapMap_Map.txt",
    #                           "-mp", "9",
    #                           "-mem", "4g",
    #                           "-cf", "0", "0", "0", "0", "0"])

    ## Only MM.
    # args = parser.parse_args(["--input", "data/Target/HM_CEU.FOUNDERS.filt",
    #                           "--out", "tests/_3_CookHLA/20190708_MM/_3_HM_CEU_T1DGC_REF",
    #                           "-ref", "data/HLA_PANEL/T1DGC/T1DGC_REF",
    #                           "-ml",
    #                           "-an", "tests/HM_CEU_REF.bgl.phased.FIDadj.alleles.answer",
    #                           "-mem", "4g"])




    ##### < for Publish > #####
    args = parser.parse_args()
    print(args)


    ### Manual Argument Checking
    l_control_flags = args.control_flags

    checked = True

    for i in range(0, len(l_control_flags)):
        if not (l_control_flags[i] == 0 or l_control_flags[i] == 1):
            checked = False
            break

    if checked:
        l_control_flags = tuple(l_control_flags)
        # print(l_control_flags)
    else:
        print(std_ERROR_MAIN_PROCESS_NAME + "The value of '--control-flags/-cf' must be five 0 or 1s. Please check it again.")
        sys.exit()


    CookHLA_lab(args, _control_flags=l_control_flags)