#-*- coding: utf-8 -*-

import os, sys, re
from os.path import join
from shutil import which
import argparse, textwrap

from src.RUN_Bash import RUN_Bash
from src.MakeGeneticMap.Panel_subset import Panel_Subset

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]



def MakeGeneticMap(_input, _reference, _out,
                   _p_src="./src/MakeGeneticMap", _p_dependency="./dependency",
                   __save_intermediates=False):


    _p_plink = which("plink")
    _p_linkage2beagle = os.path.join(_p_dependency, "linkage2beagle.jar")
    _p_beagle2linkage = os.path.join(_p_dependency, "beagle2linkage.jar")
    _p_transpose = os.path.join(_p_dependency, "transpose.jar")
    _p_mach = os.path.join(_p_dependency, "mach1")


    PLINK = "{} --noweb --silent --allow-no-sex".format(_p_plink)
    LINKAGE2BEAGLE = 'java -jar {}'.format(_p_linkage2beagle)
    RANDOMIZE_FAM = 'Rscript {}/STEP0_randomize_the_sample_about_fam_03_06_2017-COOK-V1.R'.format(_p_src)
    BGL2GC_TRICK_BGL = 'Rscript {}/bgl2GC_trick_bgl-v1.1COOK-02222017.R'.format(_p_src)
    BGL2BED = os.path.join(_p_src, 'Panel-BGL2BED.sh')
    STEP4_buildMap = os.path.join(_p_src, 'STEP4-buildMap.R')
    STEP5_collapseHLA = os.path.join(_p_src, 'STEP5-collapseHLA.R')


    # Intermediate path.
    if not _out:
        print(std_ERROR_MAIN_PROCESS_NAME + 'The argument "{0}" has not been given. Please check it again.\n'.format("--out"))
        sys.exit()
    else:
        _out = _out if not _out.endswith('/') else _out.rstrip('/')
        if bool(os.path.dirname(_out)): os.makedirs(os.path.dirname(_out), exist_ok=True)


    OUTPUT_dir = os.path.dirname(_out)
    OUTPUT_INPUT = os.path.join(OUTPUT_dir, os.path.basename(_input)) # Generated in output folder
    OUTPUT_REF = os.path.join(OUTPUT_dir, os.path.basename(_reference))




    ###### < Control Flags > ######

    RANDOM = 1
    SUBSET_BGL = 1
    MAKING_MACH_INPUT = 1
    RUNMACH = 1
    BUILDING_MAP = 1
    Cleanup = 1



    if RANDOM:

        RUN_Bash('awk \'{print $1" "$2" ""0"" ""0"" "$5" "$6}\' %s > %s' % (_input+'.fam', OUTPUT_INPUT+'.trick.fam'))
        RUN_Bash('awk \'{print $1" "$2" ""0"" ""0"" "$5" "$6}\' %s > %s' % (_reference+'.fam', OUTPUT_REF+'.trick.fam'))

        RUN_Bash(RANDOMIZE_FAM + ' {} {}'.format(OUTPUT_INPUT+'.trick.fam', OUTPUT_INPUT+'.rearranged.fam'))
        RUN_Bash(RANDOMIZE_FAM + ' {} {}'.format(OUTPUT_REF+'.trick.fam', OUTPUT_REF+'.rearranged.fam'))

        RUN_Bash('rm {}'.format(OUTPUT_INPUT+'.trick.fam'))
        RUN_Bash('rm {}'.format(OUTPUT_REF+'.trick.fam'))


    if SUBSET_BGL:

        RUN_Bash('head -100 %s | awk \'{print $1" "$2}\' > %s' % (OUTPUT_INPUT+'.rearranged.fam', OUTPUT_INPUT+'.subset.samples'))
        RUN_Bash('head -100 %s | awk \'{print $1" "$2}\' > %s' % (OUTPUT_REF+'.rearranged.fam', OUTPUT_REF+'.subset.samples'))

        RUN_Bash(PLINK + ' --bfile {} --keep {} --recode --out {}'.format(_input, OUTPUT_INPUT+'.subset.samples', OUTPUT_INPUT+'.subset'))
        RUN_Bash(PLINK + ' --bfile {} --keep {} --make-bed --out {}'.format(_input, OUTPUT_INPUT+'.subset.samples', OUTPUT_INPUT+'.subset'))

        RUN_Bash("cut -d ' ' -f1-5,7- {} > {}".format(OUTPUT_INPUT+'.subset.ped', OUTPUT_INPUT+'.subset.nopheno.ped'))
        RUN_Bash('awk \'{print "M " $2}\' %s > %s' % (OUTPUT_INPUT+'.subset.map', OUTPUT_INPUT+'.subset.dat'))


        RUN_Bash(LINKAGE2BEAGLE + ' pedigree={} data={} beagle={} standard=true > {}'.format(
            OUTPUT_INPUT+'.subset.nopheno.ped', OUTPUT_INPUT+'.subset.dat',
            OUTPUT_INPUT + '.subset.bgl.phased', OUTPUT_INPUT+'.subset.bgl.phased.log'
        ))


        RUN_Bash('awk \'{print $2" "$4" "$5" "$6}\' %s > %s' % (OUTPUT_INPUT+'.subset.bim', OUTPUT_INPUT+'.subset.markers'))


        Panel_Subset(_reference, OUTPUT_REF+'.subset.samples', 'all', OUTPUT_REF+'.subset')

        RUN_Bash(BGL2GC_TRICK_BGL + ' {} {} {} {}'.format(
                 OUTPUT_INPUT + '.subset.bgl.phased', OUTPUT_INPUT+'.subset.markers',
                 OUTPUT_INPUT + '.subset.GCchange.bgl.phased', OUTPUT_INPUT + '.subset.GCchange.markers'))

        RUN_Bash(BGL2GC_TRICK_BGL + ' {} {} {} {}'.format(
            OUTPUT_REF + '.subset.bgl.phased', OUTPUT_REF+'.subset.markers',
            OUTPUT_REF + '.subset.GCchange.bgl.phased', OUTPUT_REF + '.subset.GCchange.markers'))



        RUN_Bash('rm {}'.format(OUTPUT_INPUT+'.rearranged.fam'))
        RUN_Bash('rm {}'.format(OUTPUT_INPUT+'.subset.samples'))
        RUN_Bash('rm {}'.format(OUTPUT_REF+'.rearranged.fam'))
        RUN_Bash('rm {}'.format(OUTPUT_REF+'.subset.samples'))
        RUN_Bash('rm {}'.format(OUTPUT_INPUT+'.subset.bgl.phased'))
        RUN_Bash('rm {}'.format(OUTPUT_INPUT+'.subset.markers'))
        RUN_Bash('rm {}'.format(OUTPUT_REF+'.subset.bgl.phased'))
        RUN_Bash('rm {}'.format(OUTPUT_REF+'.subset.markers'))




    if MAKING_MACH_INPUT:

        RUN_Bash(BGL2BED+' {} {} {} {}'.format(OUTPUT_INPUT+'.subset.GCchange', OUTPUT_INPUT+'.subset.GCchange', _p_beagle2linkage, _p_plink))
        RUN_Bash(PLINK+' --bfile {} --recode --out {}'.format(OUTPUT_INPUT+'.subset.GCchange', OUTPUT_INPUT+'.subset.GCchange'))

        RUN_Bash('awk \'{print "M", $2}\' %s > %s' % (OUTPUT_INPUT+'.subset.GCchange.map', OUTPUT_INPUT+'.subset.GCchange.dat'))    # -d
        RUN_Bash('cut -d " " -f1-5,7- %s > %s' % (OUTPUT_INPUT+'.subset.GCchange.ped', OUTPUT_INPUT+'.subset.GCchange.nophe.ped'))  # -p

        RUN_Bash('cat {} | java -jar {} > {}'.format(OUTPUT_REF+'.subset.GCchange.bgl.phased', _p_transpose, OUTPUT_REF+'.subset.GCchange.bgl.phased.tr'))
        RUN_Bash('cut -d " " -f1,2,6- %s | tail -n+3 > %s' % (OUTPUT_REF+'.subset.GCchange.bgl.phased.tr', OUTPUT_REF+'.subset.GCchange.haps')) # -h
        RUN_Bash('cut -d " " -f1 %s > %s' % (OUTPUT_REF+'.subset.GCchange.markers', OUTPUT_REF+'.subset.GCchange.haps.snps'))   # -s




    if RUNMACH:

        RUN_Bash(_p_mach+' -d {} -p {} -h {} -s {} --rounds 20 --greedy --prefix {}'.format(
            OUTPUT_INPUT + '.subset.GCchange.dat',
            OUTPUT_INPUT + '.subset.GCchange.nophe.ped',
            OUTPUT_REF + '.subset.GCchange.haps',
            OUTPUT_REF + '.subset.GCchange.haps.snps',
            _out+'.mach_step'
        ))




    if BUILDING_MAP:

        RUN_Bash(STEP4_buildMap+' {} {} {} {} {} > {}'.format(
            _out+'.mach_step.erate', _out+'.mach_step.rec',
            OUTPUT_REF+'.subset.GCchange.markers', _out+'.mach_step.gmap.avg', _out+'.mach_step.gmap.last',
            _out+'.aver.erate'
        ))

        RUN_Bash(STEP5_collapseHLA+' {} {} {}'.format(_out+'.mach_step.gmap.avg', _out+'.mach_step.avg.clpsA', _out+'.mach_step.avg.clpsB'))



    if Cleanup:

        RUN_Bash('rm {}'.format(OUTPUT_INPUT+'.subset.*'))
        RUN_Bash('rm {}'.format(OUTPUT_REF+'.subset.*'))





    return 0



if __name__ == "__main__":

    ########## < Main parser > ##########

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    ###########################################################################################

        MakeGeneticMap.py
        
        INPUTS: 
        1. Plink dataset
        2. Reference dataset (beagle format)

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
    parser.add_argument("--reference", "-ref", help="\nPrefix of Reference files.\n\n", required=True)
    parser.add_argument("--out", "-o", help="\nOutput file name prefix\n\n", required=True)
    parser.add_argument("-hg", help="\nHuman Genome version(ex. 18, 19, 38)\n\n", choices=["18", "19", "38"],
                        metavar="HG", default='18')







    ##### < for Testing > #####

    # args = parser.parse_args(["--input", "/media/sf_VM_Shared/Projects/CookHLA/tests/MakeGeneticMap_seed/Korean",
    #                           "--out", "tests/MkGM_test/AGM.Korean+T1DGC_REF",
    #                           "-ref", "/media/sf_VM_Shared/Projects/CookHLA/tests/T1DGC/T1DGC_REF",
    #                           "-hg", "18"])


    ##### < for Publish > #####
    args = parser.parse_args()
    print(args)

    MakeGeneticMap(args.input, args.reference, args.out)
