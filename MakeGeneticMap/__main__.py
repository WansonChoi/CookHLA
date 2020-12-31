#-*- coding: utf-8 -*-

import os, sys, re
import argparse, textwrap

from MakeGeneticMap.MakeGeneticMap import MakeGeneticMap



class CookHLA_MakeGeneticMap(object):

    """
    Main wrappwer for MakeGeneticMap.

    """

    def __init__(self, _input, _hg_input, _reference, _out):

        ## Exception Handling here

        # (1) _input

        # (2) _reference



        self.GeneticMap = MakeGeneticMap(_input, _hg_input, _reference, _out)



if __name__ == '__main__':

    ########## < Main parser > ##########

    parser = argparse.ArgumentParser(prog='MakeGeneticMap',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    ###########################################################################################

        MakeGeneticMap.py

        INPUTS: 
        1. Target PLINK genotype file (*.{bed,bim,fam})
        2. Reference panel (PLINK + Phased Beagle file (*.{bed,bim,fam,FRQ.frq,bgl.phased,markers})

        DEPENDENCIES: (download and place in the same folder as this script)
        1. PLINK (1.9)
        2. linkage2beagle and beagle2linkage (Beagle utilities for PED <-> Beagle format)
        3. STEP0_randomize_the_sample_about_fam_03_06_2017-COOK-V1.R
        4. bgl2GC_trick_bgl-v1.1COOK-02222017.R
        5. Panel-BGL2BED (Beagle <-> BED format)
        6. Panel-subset.py
        7. STEP4-buildMap.R
        8. STEP5-collapseHLA.R


    ###########################################################################################                                     
                                     '''),
                                     add_help=False
                                     )

    ### Common arguments to share over the modules.

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="Show this help message and exit\n\n", action='help')

    parser.add_argument("--input", "-i", help="\nCommon prefix of input files.\n\n", required=True)
    parser.add_argument("--human-genome", "-hg", help="\nHuman Genome version(ex. 18, 19, 38) of TARGET(INPUT) data, Not Reference.\n\n",
                        choices=["18", "19", "38"], metavar="HG", required=True)
    parser.add_argument("--reference", "-ref", help="\nPrefix of Reference files.\n\n", required=True)
    parser.add_argument("--out", "-o", help="\nOutput file name prefix\n\n", required=True)



    ##### < for Testing > #####

    # args = parser.parse_args(["--input", "/media/sf_VM_Shared/Projects/CookHLA/tests/MakeGeneticMap_seed/Korean",
    #                           "--out", "tests/MkGM_test/AGM.Korean+T1DGC_REF",
    #                           "-ref", "/media/sf_VM_Shared/Projects/CookHLA/tests/T1DGC/T1DGC_REF")



    ##### < for Publish > #####
    args = parser.parse_args()
    print(args)



    CookHLA_MakeGeneticMap(args.input, args.human_genome, args.reference, args.out)