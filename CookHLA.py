#-*- coding: utf-8 -*-

import os, sys, re
import argparse, textwrap
from src.HLA_Imputation import HLA_Imputation

########## < Core Varialbes > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]

TOLERATED_DIFF = 0.15


def CookHLA(_input, _out, _reference, _geneticMap, _average_erate, _java_memory='2g',
            _p_src="./src", _p_dependency="./dependency", __save_intermediates=False,
            __use_Multiple_Markers=False):


    p_src = _p_src
    p_dependency = _p_dependency

    _p_plink = os.path.join(p_dependency, "plink")
    _p_beagle4 = os.path.join(p_dependency, "beagle4.jar")
    _p_linkage2beagle = os.path.join(p_dependency, "linkage2beagle.jar")
    _p_beagle2linkage = os.path.join(p_dependency, "beagle2linkage.jar")
    _p_beagle2vcf = os.path.join(p_dependency, "beagle2vcf.jar")
    _p_vcf2beagle = os.path.join(p_dependency, "vcf2beagle.jar")


    # Intermediate path.
    if not _out:
        print(std_ERROR_MAIN_PROCESS_NAME + 'The argument "{0}" has not been given. Please check it again.\n'.format("--out"))
        sys.exit()
    else:
        _out = _out if not _out.endswith('/') else _out.rstrip('/')
        if bool(os.path.dirname(_out)): os.makedirs(os.path.dirname(_out), exist_ok=True)


    ###### < Dependency Checking > ######

    ### External software

    ### Source files



    ###### < Bash Command Preparation > ######

    """
    1. plink
    2. beagle4
    3. linkage2beagle
    4. beagle2linkage
    5. beagle2vcf
    6. vcf2beagle
    7. excluding_target_snp_not_reference
    8. complete_header.R
    9. Doubling_vcf.R
    10. DP_min_selection.R
    """

    JAVATMP = _out+'.javatmpdir'
    os.makedirs(JAVATMP, exist_ok=True)


    # Memory representation check.

    p = re.compile(r'g|G$')

    if p.search(_java_memory):
        _java_memory = p.sub(repl="000m", string=_java_memory) # Gigabyte to Megabyte to use it in java.
    else:
        print(std_ERROR_MAIN_PROCESS_NAME + "Given memory for beagle4 is unappropriate.\n"
                                            "Please check '--java-memory/-mem' argument again.")
        sys.exit()


    PLINK = ' '.join([_p_plink, "--noweb", "--silent", '--allow-no-sex'])
    BEAGLE4 = ' '.join(["java", '-Djava.io.tmpdir={}'.format(JAVATMP), "-Xmx{}".format(_java_memory), "-jar", _p_beagle4])
    LINKAGE2BEAGLE = ' '.join(["java", '-Djava.io.tmpdir={}'.format(JAVATMP), "-Xmx{}".format(_java_memory), "-jar", _p_linkage2beagle])
    BEAGLE2LINKAGE = ' '.join(["java", '-Djava.io.tmpdir={}'.format(JAVATMP), "-Xmx{}".format(_java_memory), "-jar", _p_beagle2linkage])
    BEAGLE2VCF = ' '.join(["java", '-Djava.io.tmpdir={}'.format(JAVATMP), "-Xmx{}".format(_java_memory), "-jar", _p_beagle2vcf])
    VCF2BEAGLE = ' '.join(["java", '-Djava.io.tmpdir={}'.format(JAVATMP), "-Xmx{}".format(_java_memory), "-jar", _p_vcf2beagle])



    ###### < Control Flags > ######

    EXTRACT_MHC = 1
    FLIP = 1
    CONVERT_IN = 1
    IMPUTE = 1
    CONVERT_OUT = 1
    CLEAN_UP = 1



    print(std_MAIN_PROCESS_NAME + "CookHLA : Performing HLA imputation for '{}'\n"
                                  "- Java memory = {}(Mb)".format(_input, _java_memory))

    if __use_Multiple_Markers:
        print("- Using Multiple Markers.")
    if _geneticMap:
        print("- Using Genetic Map : {}.".format(_geneticMap))


    MHC = _out+'.MHC' # Prefix for MHC data.


    idx_process = 1


    if EXTRACT_MHC:

        print("[{}] Extracting SNPs from the MHC.".format(idx_process))

        idx_process += 1

    if FLIP:

        print("[{}] Performing SNP quality control.".format(idx_process))

        idx_process += 1


    ############################################################

    # if CONVERT_IN:
    #
    #     print("[{}] Converting data to beagle format.".format(idx_process))
    #     print("[{}] Converting data to reference_markers_Position.".format(idx_process))
    #     print("[{}] Converting data to target_markers_Position and extract not_including snp.".format(idx_process))
    #     print("[{}] Converting data to GC_change_beagle format.".format(idx_process))
    #     print("[{}] Converting data to vcf_format.".format(idx_process))
    #     print("[{}] Converting data to reference_phased.".format(idx_process))
    #
    #
    #     idx_process += 1
    #
    # if IMPUTE:
    #
    #     print("[{}] Performing HLA imputation (see {}.MHC.QC.imputation_out.log for progress).".format(idx_process, _out))
    #
    #     idx_process += 1
    #
    # if CONVERT_OUT:
    #
    #     print("[{}] Converting imputation vcf to beagle.".format(idx_process))
    #     print("[{}] Converting imputation GC_beagle to ori_beagle.".format(idx_process))
    #     print("[{}] Converting imputation genotypes to PLINK .ped format.".format(idx_process))
    #
    #
    #     idx_process += 1


    # This part will be taken by the instance of 'HLA_Imputation' class.

    ############################################################


    if CLEAN_UP:

        print("[{}] Clean Up.".format(idx_process))


        print("DONE!\n")

        idx_process += 1





    return 0




if __name__ == "__main__":

    ########## < Main parser > ##########

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    ###########################################################################################

        CookHLA.py

        (Created by Buhm Han.)



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


    # For publish
    # parser.add_argument("--genetic-map", "-gm", help="\nGenetic Map file.\n\n", required=True)
    # parser.add_argument("--average-erate", "-ae", help="\nAverate error rate file.\n\n", required=True)


    # For Testing
    parser.add_argument("--genetic-map", "-gm", help="\nGenetic Map file.\n\n")
    parser.add_argument("--average-erate", "-ae", help="\nAverate error rate file.\n\n")
    parser.add_argument("--use-multiple-markers", "-ml", help="\nUsing multiple markers.\n\n", action='store_true')



    parser.add_argument("--java-memory", "-mem", help="\nMemory requried for beagle(ex. 12g).\n\n", default="2g")





    ##### < for Testing > #####

    # args = parser.parse_args(["--imgt2sequence", "-imgt", "370", "-o", "TEST/TEST", "-hg", "18"])




    ##### < for Publish > #####
    args = parser.parse_args()
    print(args)

    CookHLA(args.input, args.out, args.reference, args.genetic_map, args.average_erate, _java_memory=args.java_memory,
            __use_Multiple_Markers=args.use_multiple_markers)