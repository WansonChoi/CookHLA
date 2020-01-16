#-*- coding: utf-8 -*-

import os, sys, re
import subprocess
from os.path import join
from shutil import which
import argparse, textwrap

import pandas as pd

from src.REFERENCE import REFERENCE
from MakeGeneticMap import MakeGeneticMap
from src.redefineBPv1BH import redefineBP
from src.bgl2GC_trick_bgl import Bgl2GC
from src.GC_tricked_bgl2ori_bgl import GCtricedBGL2OriginalBGL
from src.BGL2Alleles import BGL2Alleles


std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]

# Patterns
p_ONLY_4digit = re.compile(r'^\d{4},\d{4}$')
p_4to5digit = re.compile(r'^\d{4,5},\d{4,5}$')



def CookQC(_input, _reference, _out,
           _p_src='./src', _p_dependency='./dependency', _mem='2g',
           _given_AGM=None, _given_phased=None):


    if not (bool(re.match(r'\d+[Mm]$', _mem)) or bool(re.match(r'\d+[Gg]$', _mem))):
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
    LINKAGE2BEAGLE = 'java -Xmx{} -jar {}'.format(_mem, _p_linkage2beagle)
    BEAGLE2VCF = 'java -Xmx{} -jar {}'.format(_mem, _p_beagle2vcf)
    VCF2BEAGLE = 'java -Xmx{} -jar {}'.format(_mem, _p_vcf2beagle)


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

    myREF = REFERENCE(_reference, _which=(0,1,1,0,0,0))   # 'myREF' class instance.

    """
    It will be assumed that given reference panel doesn't have genetic distance.
    
    [1] Remove markers from 'MakeReference'.
    [2] Use 'MakeGeneticMap' to get genetic distance information. (ex. Target: T1DGC(No MakeReference markers) / Reference : T1DGC_REF)
        - Remove 'MakeReference' markers except HLA markers. (i.e. remove 'AA_', 'SNP_', 'INS_' not 'HLA_')
            - By '--extract' 'T1DGC_REF'.
    [3] Perform Phasing
        - PLINK to Beagle (linakge2beagle.jar)
        - GCtrick
    [4] Find wrongly phased samples
        - vcf(phased) to Beagle (vcf2beagle.jar)
        - BGL2Alleles.py
    
    [5] Split (1) wrongly phased samples PLINK dataset(K), (2) Rightly phased samples PLINK dataset(N-K)
        - Just remove samples in the output phased file. (Newly phased file). =>
        - CookHLA with prephasing strategy
        
    [6] Analyze Posterior probability of Haplotypes.
    
    """

    if bool(_given_phased) and os.path.exists(_given_phased):

        ### Using given Phased file.
        myBGL_Phased = _given_phased
        print(std_MAIN_PROCESS_NAME + "Using given Phased file('{}').".format(_given_phased))

    else:


        if bool(_given_AGM) and os.path.exists(_given_AGM+'.mach_step.avg.clpsB'):

            ### Using given AGM.
            GM = _given_AGM
            print("Using given Adaptive Genetic Map : '{}'".format(GM))

        else:
            # have to be newly created.

            ##### < Remove markers from 'MakeReference'. > #####
            toExclude_MKref = myREF.get_MKref_markers(_out=join(OUTPUT_dir, 'ToExclude.MKREF.txt'))

            ##### < [2] Use 'MakeGeneticMap' to get genetic distance information. > #####
            REF_only_Variants = myREF.PLINK_subset(_toExclude=toExclude_MKref,
                                                   _out=OUTPUT_REF + '.ONLY_Variants')  # Removing markers originated from 'MakeReference'.
            # print(REF_only_Variants)

            GM = MakeGeneticMap(REF_only_Variants, myREF.prefix, _out=join(OUTPUT_dir, 'AGM.{}+{}'.format(os.path.basename(REF_only_Variants), os.path.basename(myREF.prefix))))
            print(GM)


            # removal
            subprocess.call(['rm', REF_only_Variants+'.bed'])
            subprocess.call(['rm', REF_only_Variants+'.bim'])
            subprocess.call(['rm', REF_only_Variants+'.fam'])
            subprocess.call(['rm', REF_only_Variants+'.log'])
            subprocess.call(['rm', toExclude_MKref])




        df_clpsB = pd.read_csv(GM+'.mach_step.avg.clpsB', sep='\s+', header=None, names=['Chr', 'Label', 'GD', 'BP'])
        # print("df_clpsB :\n{}\n".format(df_clpsB))

        p_AA_SNP_INS = re.compile(r'AA_|SNP_|INS_')
        f_AA_SNP_INS = df_clpsB['Label'].str.match(p_AA_SNP_INS)

        df_clpsB = df_clpsB[~f_AA_SNP_INS]
        # print("df_clpsB :\n{}\n".format(df_clpsB))


        ToExtract_VariantsAndHLA = join(OUTPUT_dir, 'ToExtract.VariantsAndHLA.txt')
        df_clpsB['Label'].to_csv(ToExtract_VariantsAndHLA, header=False, index=False)

        REF_VariantsAndHLA = myREF.PLINK_subset(OUTPUT_REF+'.ONLY_Variants_HLA', _toExtract=ToExtract_VariantsAndHLA)

        df_bim_REF_VariantsAndHLA = pd.read_csv(REF_VariantsAndHLA+'.bim', sep='\s+', header=None, names=['Chr', 'Label', 'GD', 'BP', 'al1', 'al2'])

        # Just replacing
        df_bim_REF_VariantsAndHLA['GD'] = df_clpsB['GD'].reset_index(drop=True) # *** Potentially error-proun.
        df_bim_REF_VariantsAndHLA.to_csv(REF_VariantsAndHLA+'.bim', sep='\t', header=False, index=False)


        # # Using pd.merge()
        # df_merge0 = df_bim_REF_VariantsAndHLA.merge(df_clpsB, on='Label')
        # df_merge0 = df_merge0[['Chr_x', 'Label', 'GD_y', 'BP_x', 'al1', 'al2']]
        # df_merge0.to_csv(REF_VariantsAndHLA+'.bim', sep='\t', header=False, index=False)
        # print("df_merge0:\n{}\n".format(df_merge0))
        # print(df_clpsB.shape[0])
        # print(df_bim_REF_VariantsAndHLA.shape[0])


        # removal
        subprocess.call(['rm', ToExtract_VariantsAndHLA])
        del(df_bim_REF_VariantsAndHLA)

        # subprocess.call(['rm', GM+'.aver.erate'])
        # subprocess.call(['rm', GM+'.mach_step.avg.clpsA'])
        # subprocess.call(['rm', GM+'.mach_step.avg.clpsB'])
        # subprocess.call(['rm', GM+'.mach_step.erate'])
        # subprocess.call(['rm', GM+'.mach_step.gmap.avg'])
        # subprocess.call(['rm', GM+'.mach_step.gmap.last'])
        # subprocess.call(['rm', GM+'.mach_step.rec'])




        ##### < [3] Perform Phasing > #####

        subprocess.check_output([_p_plink, '--recode', '--keep-allele-order', '--bfile', REF_VariantsAndHLA, '--out', REF_VariantsAndHLA])

        # Making '*.markers'
        os.system(' '.join(['awk \'{print $2 " " $4 " " $5 " " $6}\'', REF_VariantsAndHLA+'.bim', '>', REF_VariantsAndHLA+'.markers']))
        # Making '*.nopheno.ped'
        os.system(' '.join(['awk \'{print "M " $2}\'', REF_VariantsAndHLA+'.map', '>', REF_VariantsAndHLA+'.dat']))
        # Making '*.dat'
        os.system(' '.join(["cut -d ' ' -f1-5,7-", REF_VariantsAndHLA+'.ped', '>', REF_VariantsAndHLA+'.nopheno.ped']))

        # Converting PLINK file to Beagle file with 'linkage2beagle.jar'
        os.system(' '.join([LINKAGE2BEAGLE,
                             'pedigree={}'.format(REF_VariantsAndHLA+'.nopheno.ped'),
                             'data={}'.format(REF_VariantsAndHLA+'.dat'),
                             'beagle={}'.format(REF_VariantsAndHLA+'.bgl'),
                             'standard=true', '>', REF_VariantsAndHLA+'.bgl.l2b.log']))

        myBGL = REF_VariantsAndHLA+'.bgl' # (***) the beagle file generated by the above linkage2beagle.jar
        myMarkers = REF_VariantsAndHLA+'.markers'

        # Removal
        subprocess.call(['rm', REF_VariantsAndHLA+'.nopheno.ped'])
        subprocess.call(['rm', REF_VariantsAndHLA+'.dat'])
        subprocess.call(['rm', REF_VariantsAndHLA+'.ped'])
        subprocess.call(['rm', REF_VariantsAndHLA+'.map'])


        ### GCtrick

        # 1. redefineBPv1BH.py
        redefineBP(myMarkers, REF_VariantsAndHLA+'.redefined.markers')
        # 2. GCchange trick
        [myBGL_GCtrick, myMarkers_redefined_GCtrcik] = \
            Bgl2GC(myBGL, REF_VariantsAndHLA+'.redefined.markers',
                   REF_VariantsAndHLA+'.GCtrick.bgl', REF_VariantsAndHLA+'.GCtrick.redefined.markers')

        # Removal
        subprocess.call(['rm', REF_VariantsAndHLA+'.redefined.markers'])


        ### beagle2vcf
        """
        usage: java -jar beagle2vcf.jar [chrom] [markers] [bgl] [missing] > [vcf]
        """
        command = ' '.join([BEAGLE2VCF, '6', myMarkers_redefined_GCtrcik, myBGL_GCtrick, '0', '>', REF_VariantsAndHLA+'.GCtrick.bgl.vcf'])
        os.system(command)


        ### Phasing
        """
        java -jar beagle.27Jan18.7e1.jar gl=test.27Jan18.7e1.vcf.gz out=out.gl
        """

        # print("Performing Phasing.")
        SP = subprocess.call([_p_beagle4, '-Xmx{}'.format(_mem),
                            'gt={}'.format(REF_VariantsAndHLA+'.GCtrick.bgl.vcf'),
                            'out={}'.format(REF_VariantsAndHLA+'.GCtrick.bgl.phased')],
                             stdout=open(REF_VariantsAndHLA+'.GCtrick.bgl.phased.vcf.log', 'w'))

        # Phasing output check.
        if not (SP == 0 and os.path.exists(REF_VariantsAndHLA+'.GCtrick.bgl.phased.vcf.gz')):
            print(std_ERROR_MAIN_PROCESS_NAME + "Phasing ('{}') failed.".format(REF_VariantsAndHLA+'.GCtrick.bgl.phased.vcf.gz'))
            return -1
        else:
            # removal
            subprocess.call(['rm', REF_VariantsAndHLA+'.GCtrick.bgl.vcf'])

        myBGL_Phased_vcfgz = REF_VariantsAndHLA+'.GCtrick.bgl.phased.vcf.gz'


        ### vcf2beagle
        """
        usage: cat [vcf file] | java -jar vcf2beagle.jar [missing] [prefix]
        """
        command = ' '.join(['gunzip -c', myBGL_Phased_vcfgz, '|', VCF2BEAGLE, '0', REF_VariantsAndHLA+'.GCtrick.phased'])
        os.system(command)

        if not os.path.exists(REF_VariantsAndHLA + '.GCtrick.phased.bgl.gz'):
            print(std_ERROR_MAIN_PROCESS_NAME + "Failed to generate Beagle file by vcf2beagle.jar.")
            sys.exit()
        else:
            command = ' '.join(['gunzip -c', REF_VariantsAndHLA + '.GCtrick.phased.bgl.gz', '>', REF_VariantsAndHLA + '.GCtrick.bgl.phased'])
            os.system(command)

            # removal
            subprocess.call(['rm', REF_VariantsAndHLA + '.GCtrick.phased.bgl.gz'])
            subprocess.call(['rm', REF_VariantsAndHLA + '.GCtrick.phased.int'])
            subprocess.call(['rm', REF_VariantsAndHLA + '.GCtrick.phased.markers'])

        myBGL_Phased = (REF_VariantsAndHLA + '.GCtrick.bgl.phased') # Generated by 'vcf2begale.jar' and gunziped.
        myBGL_Phased = GCtricedBGL2OriginalBGL(myBGL_Phased, myMarkers, REF_VariantsAndHLA+'.bgl.phased')




    ##### < [] Find wrongly phased samples. > #####

    myAlleles = BGL2Alleles(myBGL_Phased, OUTPUT_REF+'.alleles', 'all')
    # myAlleles = 'tests/T1DGC_CookQC/T1DGC_REF.alleles'

    wrong_samples = WronglyPhased(myAlleles, _out=join(OUTPUT_dir, 'WronglyPhased.samples.txt'), _fam=myREF.fam)
    # print(wrong_samples)

    # removal
    # subprocess.call(['rm', myAlleles])
    # subprocess.call(['rm', myBGL_Phased])




    ##### < [] Split (1) wrongly phased samples PLINK dataset(K), (2) Rightly phased samples PLINK dataset(N-K) > #####

    ### (1) K wrongly phased PLINK data. => (Target)
    toExclude_MKref = myREF.get_MKref_markers(_out=join(OUTPUT_dir, 'ToExclude.MKref.txt'))
    myREF_K = myREF.PLINK_subset(_toExclude=toExclude_MKref, _toKeep=wrong_samples, _out=OUTPUT_REF+'.WrongPhase')


    ### (2) N-K rightly phased PLINK data. => (Reference)
    myREF_N_K = myREF.PLINK_subset(_toRemove=wrong_samples, _out=OUTPUT_REF+'.RightPhase')
    
    # BGL subset





    ##### < Remove markers from 'MakeReference'. > #####


    return 0




def WronglyPhased(_alleles, _only4digit=True, _out=None, _fam=None):

    df_alleles = pd.read_csv(_alleles, sep='\s+', header=None, usecols=[1,2,4], names=['IID', 'HLA', '4digit']) \
                    .pivot('IID', 'HLA', '4digit')

    p_toUse = p_ONLY_4digit if _only4digit else p_4to5digit

    f = df_alleles.applymap(lambda x : bool(p_toUse.match(x))).apply(lambda x : x.all(), axis=1)

    df_wronglyPhased = df_alleles[~f]

    # To retrieve 'FID' info, *.fam file must be given.
    if isinstance(_fam, pd.DataFrame):

        f2 = _fam['IID'].isin(df_wronglyPhased.index.tolist())

        df_RETURN = _fam[f2].loc[:, ['FID', 'IID']]

    else:

        df_RETURN = pd.concat([df_wronglyPhased.index.to_series(), df_wronglyPhased.index.to_series()], axis=1)
        df_RETURN.columns = ['FID', 'IID']

    # print("df_RETURN:\n{}\n".format(df_RETURN))

    if bool(_out):
        df_RETURN.to_csv(_out, sep='\t', index=False, header=True)
        return _out
    else:
        return df_RETURN






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

    parser.add_argument("--given-phased", "-ph",
                        help="\n(For Testing Purpose) Passing prephased result manually(Not VCF, No GCtrick).\n"
                             "If given, the process will be done to this phased file.\n\n")

    parser.add_argument("--given-GM", "-gGM",
                        help="\n(For Testing Purpose) Passing the prefix of already created AGM.\n"
                             "If given, the process will be done to this AGM file.\n\n")





    ##### < for Testing > #####

    args = parser.parse_args(["-ref", "tests/T1DGC/T1DGC_REF",
                              "-o", "tests/T1DGC_CookQC/T1DGC_REF.CookQC",
                              '-mem', '12g',
                              '-ph', 'tests/T1DGC_CookQC/T1DGC_REF.ONLY_Variants_HLA.bgl.phased',
                              '-gGM', 'tests/T1DGC_CookQC/AGM.T1DGC_REF.ONLY_Variants+T1DGC_REF'])

    # args = parser.parse_args(["-ref", "tests/PAN-ASIAN/Pan-Asian_REF",
    #                           "-o", "tests/PAN-ASIAN_CookQC/Pan-Asian_REF.CookQC",
    #                           '-mem', '12g'])


    ##### < for Publish > #####
    # args = parser.parse_args()
    print(args)

    CookQC(args.input, args.reference, args.out, _mem=args.java_memory,
           _given_phased=args.given_phased, _given_AGM=args.given_GM)
