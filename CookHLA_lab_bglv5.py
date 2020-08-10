#-*- coding: utf-8 -*-

import os, sys, re
from os.path import join, dirname, basename, exists
import argparse, textwrap
import pandas as pd
from time import time

from CookHLA import CookHLA



########## < Core Varialbes > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (basename(__file__))

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]



def CookHLA_lab_bglv5(_args, _control_flags=(1,1,1)):

    # Arguments assignment.

    TARGET = _args.target
    OUT = _args.out
    REFRENCE = _args.reference
    HG = _args.hg
    GeneticMap = _args.genetic_map
    AverageErate = _args.average_erate
    ANSWER = _args.answer
    ANSWER2 = _args.answer2 if _args.answer2 else ANSWER
    HapMap_Map = _args.hapmap_map

    JAVA_MEM = _args.java_memory
    N_threads = args.nthreads


    # Output directory processing
    OUTPUT_dir = os.path.dirname(OUT)
    OUTPUT_prefix = os.path.basename(OUT)



    ### Exception handlings

    # Check Control Flags
    if len(_control_flags) != 3:
        print(std_ERROR_MAIN_PROCESS_NAME + "Wrong '_control_flags' value. Please check it again.")
        sys.exit()

    if _control_flags == (1,1,1): # default value

        F_1_Vanilla = 1
        F_2_Vanilla_HapMapMap = 1
        F_3_MM_AGM = 1

    else:
        [F_1_Vanilla, F_2_Vanilla_HapMapMap, F_3_MM_AGM] = _control_flags


    # Check Input number
    if not (len(TARGET) == len(ANSWER)):
        print(std_ERROR_MAIN_PROCESS_NAME + "There must be same number of ANSWER entries. Please check '--answer' argument again.")
        sys.exit()

    if len(GeneticMap) != len(AverageErate):
        print(std_ERROR_MAIN_PROCESS_NAME + "There must be same number of Genetic map and Average Error Rate entries. "
                                            "Please check '--genetic-map' and '--average-erate' arguments again.")
        sys.exit()

    if len(GeneticMap) != len(REFRENCE):
        print(std_ERROR_MAIN_PROCESS_NAME + "There must be same number of Genetic map and Reference entries. "
                                            "Please check '--genetic-map' and '--reference' arguments again.")
        sys.exit()


    # Accuracy outputs.
    # __accuracies__ = {1: None,
    #                   2: None}

    __accuracies__ = {i: {j : [None, None, None] for j in range(len(REFRENCE))} for i in range(len(TARGET))}
    print(__accuracies__)




    print(TARGET)
    print(REFRENCE)
    print(GeneticMap)
    print(AverageErate)
    print(ANSWER)


    ##### < Main > #####

    for i in range(len(TARGET)):
        for j in range(len(REFRENCE)):

            Nth_label_TAR = '_{}_TAR'.format(i+1)
            Nth_label_REF = '_{}_REF'.format(j+1)

            print("\n\n\n\n")

            if existsTARGET(TARGET[i]) and existsREFERENCE(REFRENCE[j]):

                print(std_MAIN_PROCESS_NAME +
                      "CookHLA_lab for\n"
                      "{}: {}\n"
                      "{}: {}".format(Nth_label_TAR, TARGET[i],
                                      Nth_label_REF, REFRENCE[j]))

                if F_1_Vanilla:

                    ###### _1_Vanilla ######
                    print("\n\n<Imputation : _1_Vanilla>")

                    OUT_1_Vanilla = join(OUTPUT_dir, Nth_label_TAR+Nth_label_REF,
                                         '_1_Vanilla', OUTPUT_prefix + '.Vanilla')

                    time_start_1_Vanilla = time()

                    [t_HLA_Imptation_out, t_accuracy] = \
                        CookHLA(TARGET[i], OUT_1_Vanilla, REFRENCE[j],
                                __use_Multiple_Markers=False, _MultP=_args.multiprocess,
                                _answer=ANSWER[i], _java_memory=JAVA_MEM, f_BEAGLE5=True,
                                _nthreads=N_threads)

                    time_end_1_Vanilla = time()
                    print("Implementation time of 1_Vanilla : {}(min)".format(
                        (time_end_1_Vanilla - time_start_1_Vanilla) / 60))

                    __accuracies__[i][j][0] = t_accuracy


                if F_2_Vanilla_HapMapMap:

                    ###### _2_Vanilla_HapMap_Map ######
                    print("\n\n<Imputation : _2_Vanilla_HapMap_Map>")

                    OUT_2_Vanilla_HapMap_Map = join(OUTPUT_dir, Nth_label_TAR+Nth_label_REF,
                                                '_2_Vanilla_HapMap_Map', OUTPUT_prefix + '.Vanilla_HapMap_Map')

                    time_start_2_Vanilla_HapMapMap = time()

                    [t_HLA_Imptation_out, t_accuracy] = \
                        CookHLA(TARGET[i], OUT_2_Vanilla_HapMap_Map, REFRENCE[j],
                                __use_Multiple_Markers=False, _MultP=_args.multiprocess,
                                _HapMap_Map=HapMap_Map,
                                _answer=ANSWER[i], _java_memory=JAVA_MEM, f_BEAGLE5=True,
                                _nthreads=N_threads)

                    time_end_2_Vanilla_HapMapMap = time()
                    print("Implementation time of 2_Vanilla_HapMapMap : {}(min)".format(
                        (time_end_2_Vanilla_HapMapMap - time_start_2_Vanilla_HapMapMap) / 60))

                    __accuracies__[i][j][1] = t_accuracy


                if F_3_MM_AGM:

                    ###### _3_MM_AGM ######
                    print("\n\n<Imputation : _3_MM+AGM (full CookHLA)>")

                    OUT_3_MM_AGM = join(OUTPUT_dir, Nth_label_TAR+Nth_label_REF,
                                        '_3_MM_AGM', OUTPUT_prefix + '.MM+AGM')

                    time_start_3_MM_AGM = time()

                    [t_HLA_Imptation_out, t_accuracy] = \
                        CookHLA(TARGET[i], OUT_3_MM_AGM, REFRENCE[j],
                                __use_Multiple_Markers=True, _MultP=_args.multiprocess,
                                _AdaptiveGeneticMap=GeneticMap[j], _Average_Erate=AverageErate[j],
                                _answer=ANSWER[i], _java_memory=JAVA_MEM, f_BEAGLE5=True,
                                _nthreads=N_threads)

                    time_end_3_MM_AGM = time()

                    print("Implementation time of _3_MM_AGM : {}(min)".format((time_end_3_MM_AGM - time_start_3_MM_AGM) / 60))

                    __accuracies__[i][j][2] = t_accuracy

                CollectTable(join(OUTPUT_dir, Nth_label_TAR+Nth_label_REF+".ACCURACY_TABLE.txt"),
                             __accuracies__[i][j])

            else:
                print(std_WARNING_MAIN_PROCESS_NAME +
                      "CookHLA_lab for\n"
                      "Target: {}\n"
                      "Reference: {}\n"
                      "FAILED!\n".format(TARGET[i], REFRENCE[j]))
                os.makedirs(join(OUTPUT_dir, Nth_label_TAR+Nth_label_REF+'_failed'), exist_ok=True)


    # print(__accuracies__)

    # return __RETURN__



def CollectTable(_out, acc):

    # l_label = ['_1_Vanilla', '_2_Vanilla_HapMapMap', '_3_MM+AGM'] # Hard-coded.

    df_acc = pd.concat([pd.read_csv(item, sep='\s+', header=None, index_col=0, names=['HLA', 'acc'])
                        for item in acc], axis=1) \
                .applymap(lambda x: None if x <= 0 else x) \
                .dropna()

    # df_acc.columns = l_label

    sr_mean = df_acc.mean(axis=0)
    df_acc2 = df_acc.append(sr_mean, ignore_index=True)

    idx_df_acc2 = df_acc.index.tolist()
    idx_df_acc2.append('avg')
    df_acc2.index = idx_df_acc2
    df_acc2.index.name = 'HLA'


    df_acc2.to_csv(_out, sep='\t', header=True, index=True)

    return _out



def existsTARGET(_target):
    # print((exists(_target+'.bed') and exists(_target+'.bim') and exists(_target+'.fam')))
    return (exists(_target+'.bed') and exists(_target+'.bim') and exists(_target+'.fam'))

def existsREFERENCE(_reference):
    # print((exists(_reference+'.bed') and
    #         exists(_reference+'.bim') and
    #         exists(_reference+'.fam') and
    #         exists(_reference+'.FRQ.frq') and
    #         exists(_reference+'.bgl.phased') and
    #         exists(_reference+'.markers')))
    return (exists(_reference+'.bed') and
            exists(_reference+'.bim') and
            exists(_reference+'.fam') and
            exists(_reference+'.FRQ.frq') and
            exists(_reference+'.bgl.phased') and
            exists(_reference+'.markers'))





if __name__ == "__main__":

    ########## < Main parser > ##########

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    ###########################################################################################

        CookHLA_lab_bglv5.py

        (Created by Wanson Choi.)
        
        To automate generating accuracy comparison table, this module has been introduced.
        Each cookhla imputation will be implemented as execution of functions and their accuracy 
        outputs will be summarized by pandas.
        
        Adapted to work on the accuracy comparison job based on BEAGLE5. In BEAGLE5 framework, 
        Only (1) 'vanilla Beagle5 + hapmap genetic map' and (2) 'MM+AGM(full CookHLA)' will be
        investigated.
        
        (1) _1_Plain+HapMap_Map
        (2) _2_MM_AGM (full CookHLA)
        
        In addition, multiple targets and references can be used. Answer information for target
        data must be the same number of target data sets. Genetic maps must be the same number of
        reference data sets.
                
        To execute all above imputations, next arguments must be given to 'CookHLA_lab_bglv5.py'
        
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

    parser.add_argument("--target", "-t", help="\nTarget data set file prefix.\n\n", required=True, nargs='+')
    parser.add_argument("--reference", "-ref", help="\nReference data file prefix.\n\n", required=True, nargs='+')
    parser.add_argument("--out", "-o", help="\nOutput file name prefix\n\n", required=True)
    parser.add_argument("-hg", help="\nHuman Genome version(ex. 18, 19, 38)\n\n", choices=["18", "19", "38"],
                        metavar="HG", default='18')


    # For Testing
    parser.add_argument("--genetic-map", "-gm", help="\nGenetic Map file.\n\n", required=True, nargs='+')
    parser.add_argument("--average-erate", "-ae", help="\nAverate error rate file.\n\n", required=True, nargs='+')
    # parser.add_argument("--use-multiple-markers", "-ml", help="\nUsing multiple markers.\n\n", action='store_true')

    # parser.add_argument("--prephasing", "-pr", help="\nUtilizing prephasing strategy.\n\n", action='store_true')

    parser.add_argument("--answer", "-an", help="\nAnswer file to calculate imputation accuracy.\n\n", required=True, nargs='+')
    parser.add_argument("--answer2", "-an2", help="\nAnswer file to calculate imputation accuracy(Fam colum fixed).\n\n", nargs='+')

    parser.add_argument("--multiprocess", "-mp", help="\nSetting parallel multiprocessing.\n\n", type=int, choices=[2,3,4,5,6,7,8,9], nargs='?', default=1, const=3)

    parser.add_argument("--java-memory", "-mem", help="\nMemory requried for beagle(ex. 12g).\n\n", default="2g")
    parser.add_argument("--nthreads", "-nth",
                        help="\nThe number of threads for each BEAGLE implementation (default: 1).\n\n", default=1, type=int)

    parser.add_argument("--hapmap-map", "-hm",
                        help="\n(For Testing Purpose) Hapmap Map(Adaptive Genetic Map).\n\n", required=True)

    parser.add_argument("--control-flags", "-cf",
                        help="\nBoolean sequence to nominate which imputations are to be done.\n\n", nargs=2, default=(1,1), type=int)




    ##### < for Testing > #####

    # args = parser.parse_args(["-t", "example/1958BC", "example/1958BC",
    #                           "-o", "tests/CookHLA_lab_bglv5_test/JustPrefix",
    #                           "-ref", "example/HM_CEU_REF", "example/HM_CEU_REF",
    #                           "-gm", "example/AGM.1958BC+HM_CEU_REF.mach_step.avg.clpsB", "example/AGM.1958BC+HM_CEU_REF.mach_step.avg.clpsB",
    #                           "-ae", "example/AGM.1958BC+HM_CEU_REF.aver.erate", "example/AGM.1958BC+HM_CEU_REF.aver.erate",
    #                           "-an", "example/1958BC.answer.1454to1401.Marked.chped", "example/1958BC.answer.1454to1401.Marked.chped",
    #                           "-hm", "tests/HapMap_Map.txt",
    #                           "-mp", "3",
    #                           "-mem", "4g"
    #                           ])





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


    CookHLA_lab_bglv5(args, _control_flags=l_control_flags)