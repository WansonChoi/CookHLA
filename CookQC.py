#-*- coding: utf-8 -*-

import os, sys, re
from os.path import join
from shutil import which
import argparse, textwrap

import pandas as pd

from src.REFERENCE import REFERENCE
from MakeGeneticMap import MakeGeneticMap



std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]





def CookQC(_input, _reference, _out,
           _p_src='./src', _p_dependency='./dependency', _mem='2g'):


    if not (bool(re.match(r'\d[Mm]$', _mem)) or bool(re.match(r'\d[Gg]$', _mem))):
        print(std_ERROR_MAIN_PROCESS_NAME + "Wrong value for Memry('{}').\n"
                                            "Please check the '-mem' argument again.".format(_mem))
        sys.exit()


    # dependent software
    _p_plink = which('plink')
    _p_beagle4 = which('beagle')
    _p_linkage2beagle = os.path.join(_p_dependency, "linkage2beagle.jar")
    _p_beagle2linkage = os.path.join(_p_dependency, "beagle2linkage.jar")
    _p_beagle2vcf = os.path.join(_p_dependency, "beagle2vcf.jar")
    _p_vcf2beagle = os.path.join(_p_dependency, "vcf2beagle.jar")

    # Command
    PLINK = "{} --silent --allow-no-sex".format(_p_plink)
    LINKAGE2BEAGLE = 'java -jar {}'.format(_p_linkage2beagle)


    # Intermediate path.
    if not _out:
        print(std_ERROR_MAIN_PROCESS_NAME + 'The argument "{0}" has not been given. Please check it again.\n'.format("--out"))
        sys.exit()
    else:
        _out = _out if not _out.endswith('/') else _out.rstrip('/')
        if bool(os.path.dirname(_out)): os.makedirs(os.path.dirname(_out), exist_ok=True)

    OUTPUT_dir = os.path.dirname(_out)
    # OUTPUT_INPUT = os.path.join(OUTPUT_dir, os.path.basename(_input)) # Generated in output folder
    OUTPUT_REF = os.path.join(OUTPUT_dir, os.path.basename(_reference))




    """
    < Briefing >
    
    Given reference panel, CookQC performs additional QCs.
    
    - Perfect Phasing (suggested by B. Han.)
    
    """


    ########## < [0] Loading Reference panel data > ##########

    REF = REFERENCE(_reference, _which=(0,1,1,0,0,0))

    """
    if hasGD => Just collapse
    if not => MakeGenetic Map
    """

    # if REF.hasGD():
    #     # Just collapse
    #     pass
    #
    # else:
    #     pass


    toExclude_MKREF = REF.get_MKref_markers(_out=join(OUTPUT_dir, 'ToExclude.ExceptHLA.txt'))

    REF_only_Variants = REF.PLINK_subset(_toExclude=toExclude_MKREF, _out=OUTPUT_REF+'.ONLY_Variants')
    print(REF_only_Variants)


    GM = MakeGeneticMap(REF_only_Variants, REF.prefix, _out=join(OUTPUT_dir, 'AGM.{}+{}'.format(os.path.basename(REF_only_Variants), os.path.basename(REF.prefix))))

    print(GM)


    return 0




def PHASING():

    return 0






if __name__ == "__main__":

    ########## < Main parser > ##########

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    ###########################################################################################

        CookQC.py

        - Perform additional QC to reference panel data.



    ###########################################################################################
                                     '''),
                                     add_help=False
                                     )

    ### Common arguments to share over the modules.

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="Show this help message and exit\n\n", action='help')

    parser.add_argument("--input", "-i", help="\nCommon prefix of input files.\n\n")
    parser.add_argument("--reference", "-ref", help="\nPrefix of Reference files.\n\n", required=True)
    parser.add_argument("--out", "-o", help="\nOutput file name prefix\n\n", required=True)
    parser.add_argument("-hg", help="\nHuman Genome version(ex. 18, 19, 38)\n\n", choices=["18", "19", "38"],
                        metavar="HG", default='18')

    parser.add_argument("--java-memory", "-mem", help="\nMemory requried for beagle(ex. 12g).\n\n", default="2g")






    ##### < for Testing > #####

    args = parser.parse_args(["-ref", "tests/T1DGC/T1DGC_REF",
                              "-o", "tests/T1DGC_CookQC/T1DGC_REF.CookQC"])


    ##### < for Publish > #####
    # args = parser.parse_args()
    print(args)

    CookQC(args.input, args.reference, args.out)
