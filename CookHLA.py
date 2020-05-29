#-*- coding: utf-8 -*-

import os, sys, re
from os.path import join, exists
from shutil import which
import argparse, textwrap
from time import time

from src.HLA_Imputation import HLA_Imputation
from src.HLA_Imputation_GM import HLA_Imputation_GM


########## < Core Varialbes > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]

TOLERATED_DIFF = 0.15




def CookHLA(_input, _out, _reference, _hg='18', _AdaptiveGeneticMap=None, _Average_Erate=None, _java_memory='2g',
            _MultP=1, _answer=None, __save_intermediates=False, __use_Multiple_Markers=False, _p_src="./src",
            _p_dependency="./dependency", _given_prephased=None, f_prephasing=False, _HapMap_Map=None,
            __overlap__=(0.5, 1, 1.5), _window=5, _ne=1000000, _nthreads=1, f_measureAcc_v2=False):


    ### Argument exception

    if _given_prephased and not f_prephasing:
        print(std_ERROR_MAIN_PROCESS_NAME + "The arguments '--prephased/-ph' must be used with '--prephasing/nph'. Please check them again.")
        sys.exit()


    if (_AdaptiveGeneticMap and _Average_Erate) and _HapMap_Map:
        print(std_ERROR_MAIN_PROCESS_NAME + "The arguments '--hapmap-map/-hm' can't be used with '--genetic-map/-gm' and '--average-erate/-ae'. Please check them again.")
        sys.exit()


    p_src = _p_src
    p_dependency = _p_dependency

    # _p_plink = os.path.join(p_dependency, "plink")
    _p_plink = which('plink')
    # _p_beagle5 = os.path.join(p_dependency, "beagle4.jar")
    _p_beagle5 = which('beagle')    # replaced by the one of Anaconda(Bioconda).
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
    if not(bool(_p_plink) and exists(_p_plink)):
        sys.stderr.write(std_ERROR_MAIN_PROCESS_NAME + "PLINK(v1.9b) can't be found.\n")
        sys.exit()

    if not(bool(_p_beagle5) and exists(_p_beagle5)):
        sys.stderr.write(std_ERROR_MAIN_PROCESS_NAME + "Beagle4 can't be found.\n")
        sys.exit()

    if not(bool(_p_linkage2beagle) and exists(_p_linkage2beagle)):
        sys.stderr.write(std_ERROR_MAIN_PROCESS_NAME + "linkage2beagle can't be found in 'dependency/' folder.\n")
        sys.exit()

    if not(bool(_p_beagle2linkage) and exists(_p_beagle2linkage)):
        sys.stderr.write(std_ERROR_MAIN_PROCESS_NAME + "beagle2linkage can't be found in 'dependency/' folder.\n")
        sys.exit()

    if not(bool(_p_beagle2vcf) and exists(_p_beagle2vcf)):
        sys.stderr.write(std_ERROR_MAIN_PROCESS_NAME + "beagle2vcf can't be found in 'dependency/' folder.\n")
        sys.exit()

    if not(bool(_p_vcf2beagle) and exists(_p_vcf2beagle)):
        sys.stderr.write(std_ERROR_MAIN_PROCESS_NAME + "vcf2beagle can't be found in 'dependency/' folder.\n")
        sys.exit()



    ### Source files

    ### Given prephased result
    if _given_prephased and not exists(_given_prephased):
        print(std_ERROR_MAIN_PROCESS_NAME + "Given prephased result file('{}') doesn't exist. Please check '--prephase/-ph' argument again.\n".format(_given_prephased))
        sys.exit()


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

    OUTPUT_dir = os.path.dirname(_out)
    OUTPUT_dir_ref = join(OUTPUT_dir, os.path.basename(_reference))

    # Memory representation check.

    p = re.compile(r'g|G$')

    if p.search(_java_memory):
        _java_memory = p.sub(repl="000m", string=_java_memory) # Gigabyte to Megabyte to use it in java.
    else:
        print(std_ERROR_MAIN_PROCESS_NAME + "Given memory for beagle4 is unappropriate.\n"
                                            "Please check '--java-memory/-mem' argument again.")
        sys.exit()


    PLINK = ' '.join([_p_plink, "--noweb", "--silent", '--allow-no-sex'])
    BEAGLE5 = ' '.join([_p_beagle5, '-Djava.io.tmpdir={}'.format(JAVATMP), "-Xmx{}".format(_java_memory)])
    LINKAGE2BEAGLE = ' '.join(["java", '-Djava.io.tmpdir={}'.format(JAVATMP), "-Xmx{}".format(_java_memory), "-jar", _p_linkage2beagle])
    BEAGLE2LINKAGE = ' '.join(["java", '-Djava.io.tmpdir={}'.format(JAVATMP), "-Xmx{}".format(_java_memory), "-jar", _p_beagle2linkage])
    BEAGLE2VCF = ' '.join(["java", '-Djava.io.tmpdir={}'.format(JAVATMP), "-Xmx{}".format(_java_memory), "-jar", _p_beagle2vcf])
    VCF2BEAGLE = ' '.join(["java", '-Djava.io.tmpdir={}'.format(JAVATMP), "-Xmx{}".format(_java_memory), "-jar", _p_vcf2beagle])

    MERGE = os.path.join(p_src, 'merge_tables.pl')
    # PARSEDOSAGE = os.path.join(p_src, 'ParseDosage.csh')
    # BGL2BED = os.path.join(p_src, 'Panel-BGL2BED.sh')




    ###### < Adaptive Genetic Map checking > ######

    __use_GeneticMap = False # Check whether Adaptive Genetic Map is available.

    if _Average_Erate and _AdaptiveGeneticMap:

        if os.path.exists(_Average_Erate) and os.path.getsize(_Average_Erate) > 0 and \
                os.path.exists(_AdaptiveGeneticMap) and os.path.getsize(_AdaptiveGeneticMap) > 0:
            __use_GeneticMap = True     # Using Adaptive Genetic Map.

        else:
            if not os.path.exists(_Average_Erate):
                print(std_ERROR_MAIN_PROCESS_NAME + "The file ('{}') doesn't exist.\n"
                                                    "Please check '--average-erate/-ae' argument again.".format(_Average_Erate))
                sys.exit()
            elif os.path.getsize(_Average_Erate) == 0:
                print(std_ERROR_MAIN_PROCESS_NAME + "The file ('{}') doesn't contain anything.\n"
                                                    "Please check '--average-erate/-ae' argument again.".format(_Average_Erate))
                sys.exit()


            if not os.path.exists(_AdaptiveGeneticMap):
                print(std_ERROR_MAIN_PROCESS_NAME + "The file ('{}') doesn't exist.\n"
                                                    "Please check '--genetic-map/-gm' argument again.".format(_AdaptiveGeneticMap))
                sys.exit()
            elif os.path.getsize(_AdaptiveGeneticMap):
                print(std_ERROR_MAIN_PROCESS_NAME + "The file ('{}') doesn't contain anything.\n"
                                                    "Please check '--average-erate/-ae' argument again.".format(_Average_Erate))
                sys.exit()



    elif not (_Average_Erate or _AdaptiveGeneticMap):

        __use_GeneticMap = False     # No using Adaptive Genetic Map.

    else:
        print(std_ERROR_MAIN_PROCESS_NAME + "Either arguments '--genetic-map(-gm)' or '--average-erate(-ae)' wasn't given.\n"
                                            "Please check whether both of them are given or not.")
        sys.exit()



    ## Hapmap Map

    __use_HapMap_Map = False

    if _HapMap_Map:

        if os.path.exists(_HapMap_Map) and os.path.getsize(_HapMap_Map) > 0:
            __use_HapMap_Map = True
        else:
            print(std_ERROR_MAIN_PROCESS_NAME + "Given HapMap map file is wrong. Please check it again.")
            sys.exit()




    ###### < Control Flags > ######

    EXTRACT_MHC = 1
    FLIP = 1
    CLEAN_UP = 1




    ################## < MAIN > ##################

    print(std_MAIN_PROCESS_NAME + "CookHLA : Performing HLA imputation for '{}'\n"
                                  "- Java memory = {}(Mb)".format(_input, _java_memory))

    if __use_Multiple_Markers:
        print("- Using Multiple Markers.")

    if __use_GeneticMap:
        print("- Using Adaptive Genetic Map.")

    if __use_HapMap_Map:
        print("- (Test Purpose) Using HapMap Map")

    if _given_prephased:
        print("- (Test Purpose) Pre-phased result given.")

    if f_prephasing:
        print("- Prephasing.")




    MHC = _out+'.MHC' # Prefix for MHC data.


    idx_process = 1


    if EXTRACT_MHC:

        print("[{}] Extracting SNPs from the MHC.".format(idx_process))

        command = ' '.join([PLINK, '--bfile', _input, '--chr 6', '--from-mb 29 --to-mb 34', '--maf 0.025', '--make-bed', '--out', MHC])
        # print(command)
        os.system(command)

        """
        Input : `_input` (from argument)
        Output : `MHC` (SNPs in the MHC region)
        """

        idx_process += 1

    if FLIP:

        print("[{}] Performing SNP quality control.".format(idx_process))

        ### Identifying non-A/T non-C/G SNPs to flip
        command = ' '.join(['echo "SNP 	POS	A1	A2"', '>', _out+'.tmp1'])
        # print(command)
        os.system(command)

        command = ' '.join(['cut -f2,4-', MHC+'.bim', '>>', _out+'.tmp1']) # Cutting out columns of `_input`
        # print(command)
        os.system(command)

        command = ' '.join(['echo "SNP 	POSR	A1R	A2R"', '>', _out+'.tmp2'])
        # print(command)
        os.system(command)

        command = ' '.join(['cut -f2,4-', _reference+'.bim', '>>', _out+'.tmp2']) # Cutting out columns of `_reference`
        # print(command)
        os.system(command)

        command = ' '.join([MERGE, _out+'.tmp2', _out+'.tmp1', 'SNP', '|', 'grep -v -w NA', '>', _out+'.SNPS.alleles'])
        # print(command)
        os.system(command)

        if not __save_intermediates:
            os.system(' '.join(['rm', _out+'.tmp1']))


        command = ' '.join(["awk '{if ($3 != $6 && $3 != $7){print $1}}'", _out+'.SNPS.alleles', '>', _out+'.SNPS.toflip1']) # Acquiring suspected SNPs to flip.
        # print(command)
        os.system(command)

        command = ' '.join([PLINK, '--bfile', MHC, '--flip', _out+'.SNPS.toflip1', '--make-bed', '--out', MHC+'.FLP']) ### Flipping those suspected SNPs.
        # print(command)
        os.system(command)

        if not __save_intermediates:
            os.system(' '.join(['rm', MHC+'.bed']))
            os.system(' '.join(['rm', MHC+'.bim']))
            os.system(' '.join(['rm', MHC+'.log']))
            if __use_Multiple_Markers:
                os.system(' '.join(['rm', MHC+'.fam'])) # MHC+'.fam' will be used in 'CONVERT_OUT' when not using multiple markers.
            os.system(' '.join(['rm', _out+'.SNPS.alleles']))
            os.system(' '.join(['rm', _out+'.SNPS.toflip1']))

        # So far : `MHC+.FLP`



        ### Calculating allele frequencies
        command = ' '.join([PLINK, '--bfile', MHC+'.FLP', '--freq', '--out', MHC+'.FLP.FRQ'])
        # print(command)
        os.system(command)

        command = ' '.join(["sed 's/A1/A1I/g'", MHC+'.FLP.FRQ.frq', '|', "sed 's/A2/A2I/g'", '|', "sed 's/MAF/MAF_I/g'", '>', _out+'.tmp'])
        # print(command)
        os.system(command)


        command = ' '.join(['mv', _out+'.tmp', MHC+'.FLP.FRQ'])
        # print(command)
        os.system(command)

        command = ' '.join([MERGE, _reference+'.FRQ.frq', MHC+'.FLP.FRQ.frq', 'SNP', '|', 'grep -v -w NA', '>', _out+'.SNPS.frq'])
        # print(command)
        os.system(command)

        command = ' '.join(["sed 's/ /\t/g'", _out+'.SNPS.frq', '|', 'awk \'{if ($3 != $8){print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $9 "\t" $8 "\t" 1-$10 "\t*"}else{print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $8 "\t" $9 "\t" $10 "\t."}}\'',
                            '>', _out+'.SNPS.frq.parsed'])
        # print(command)
        os.system(command)

        if not __save_intermediates:
            os.system(' '.join(['rm', _out+'.SNPS.frq']))
            os.system(' '.join(['rm', MHC+'.FLP.FRQ.frq']))
            os.system(' '.join(['rm', MHC+'.FLP.FRQ', MHC+'.FLP.FRQ.log']))



        ### Finding A/T and C/G SNPs
        command = ' '.join(['awk \'{if (($2 == "A" && $3 == "T") || ($2 == "T" && $3 == "A") || ($2 == "C" && $3 == "G") || ($2 == "G" && $3 == "C")){if ($4 > $7){diff=$4 - $7; if ($4 > 1-$7){corrected=$4-(1-$7)}else{corrected=(1-$7)-$4}}else{diff=$7-$4;if($7 > (1-$4)){corrected=$7-(1-$4)}else{corrected=(1-$4)-$7}};print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" diff "\t" corrected}}\'',
                            _out+'.SNPS.frq.parsed', '>', _out+'.SNPS.ATCG.frq']) # 'ATCG' literally means markers with only A, T, C, G are extracted.
        # print(command)
        os.system(command)



        ### Identifying A/T and C/G SNPs to flip or remove
        command = ' '.join(["awk '{if ($10 < $9 && $10 < .15){print $1}}'", _out+'.SNPS.ATCG.frq', '>', _out+'.SNPS.toflip2'])
        # print(command)
        os.system(command)

        command = ' '.join(["awk '{if ($4 > 0.4){print $1}}'", _out+'.SNPS.ATCG.frq', '>', _out+'.SNPS.toremove'])
        # print(command)
        os.system(command)

        if not __save_intermediates:
            os.system(' '.join(['rm', _out+'.SNPS.ATCG.frq']))



        ### Identifying non A/T and non C/G SNPs to remove
        command = ' '.join(['awk \'{if (!(($2 == "A" && $3 == "T") || ($2 == "T" && $3 == "A") || ($2 == "C" && $3 == "G") || ($2 == "G" && $3 == "C"))){if ($4 > $7){diff=$4 - $7;}else{diff=$7-$4}; if (diff > \'%s\'){print $1}}}\'' % TOLERATED_DIFF,
                            _out+'.SNPS.frq.parsed', '>>', _out+'.SNPS.toremove'])
        # print(command)
        os.system(command)

        command = ' '.join(['awk \'{if (($2 != "A" && $2 != "C" && $2 != "G" && $2 != "T") || ($3 != "A" && $3 != "C" && $3 != "G" && $3 != "T")){print $1}}\'', _out+'.SNPS.frq.parsed', '>>', _out+'.SNPS.toremove'])
        # print(command)
        os.system(command)

        command = ' '.join(['awk \'{if (($2 == $5 && $3 != $6) || ($3 == $6 && $2 != $5)){print $1}}\'', _out+'.SNPS.frq.parsed', '>>', _out+'.SNPS.toremove'])
        # print(command)
        os.system(command)

        if not __save_intermediates:
            os.system(' '.join(['rm', _out+'.SNPS.frq.parsed']))


        ### Making QCd SNP file
        command = ' '.join([PLINK, '--bfile', MHC+'.FLP', '--geno 0.2', '--exclude', _out+'.SNPS.toremove', '--flip', _out+'.SNPS.toflip2', '--make-bed', '--out', MHC+'.QC'])
        # print(command)
        os.system(command)

        if not __save_intermediates:
            os.system(' '.join(['rm', _out+'.SNPS.toremove']))
            os.system(' '.join(['rm', _out+'.SNPS.toflip2']))
            os.system(' '.join(['rm', MHC+'.FLP.bed']))
            os.system(' '.join(['rm', MHC+'.FLP.bim']))
            os.system(' '.join(['rm', MHC+'.FLP.fam']))
            os.system(' '.join(['rm', MHC+'.FLP.log']))

        command = ' '.join([PLINK, '--bfile', MHC+'.QC', '--freq', '--out', MHC+'.QC.FRQ'])
        # print(command)
        os.system(command)


        command = ' '.join(["sed 's/A1/A1I/g'", MHC+'.QC.FRQ.frq', '|', "sed 's/A2/A2I/g'", '|', "sed 's/MAF/MAF_I/g'", '>', _out+'.tmp'])
        # print(command)
        os.system(command)

        command = ' '.join(['mv', _out+'.tmp', MHC+'.QC.FRQ.frq'])
        # print(command)
        os.system(command)

        command = ' '.join([MERGE, _reference+'.FRQ.frq', MHC+'.QC.FRQ.frq', 'SNP', '|', 'grep -v -w NA', '>', _out+'.SNPS.QC.frq'])
        # print(command)
        os.system(command)

        if not __save_intermediates:
            os.system(' '.join(['rm', MHC+'.QC.FRQ.frq']))
            os.system(' '.join(['rm', MHC+'.QC.FRQ.log']))



        command = ' '.join(['cut -f2', _out+'.SNPS.QC.frq', '|', "awk '{if (NR > 1){print $1}}'", '>', _out+'.SNPS.toinclude'])
        # print(command)
        os.system(command)

        if not __save_intermediates:
            os.system(' '.join(['rm', _out+'.SNPS.QC.frq']))

        command = ' '.join(['echo "SNP 	POS	A1	A2"', '>', _out+'.tmp1'])
        # print(command)
        os.system(command)

        command = ' '.join(['cut -f2,4-', MHC+'.QC.bim' ,'>>', _out+'.tmp1'])
        # print(command)
        os.system(command)


        command = ' '.join([MERGE, _out+'.tmp2', _out+'.tmp1', 'SNP', '|', 'awk \'{if (NR > 1){if ($5 != "NA"){pos=$5}else{pos=$2}; print "6\t" $1 "\t0\t" pos "\t" $3 "\t" $4}}\'',
                            '>', MHC+'.QC.bim'])
        # print(command)
        os.system(command)

        if not __save_intermediates:
            os.system(' '.join(['rm', _out+'.tmp1']))
            os.system(' '.join(['rm', _out+'.tmp2']))



        # Recoding QC'd file as ped
        command = ' '.join([PLINK, '--bfile', MHC+'.QC', '--extract', _out+'.SNPS.toinclude', '--make-bed', '--out', MHC+'.QC.reorder'])
        # print(command)
        os.system(command)

        command = ' '.join([PLINK, '--bfile', MHC+'.QC.reorder', '--recode', '--out', MHC+'.QC'])
        # print(command)
        os.system(command)

        if not __save_intermediates:
            # os.system(' '.join(['rm', MHC+'.QC.{bed,bim,fam,log}']))
            os.system(' '.join(['rm', MHC+'.QC.reorder.bed']))
            os.system(' '.join(['rm', MHC+'.QC.reorder.bim']))
            os.system(' '.join(['rm', MHC+'.QC.reorder.fam']))
            os.system(' '.join(['rm', MHC+'.QC.reorder.log']))
            os.system(' '.join(['rm', _out+'.SNPS.toinclude']))



        # Making SNP files (pre-beagle files)
        command = ' '.join(['awk \'{print "M " $2}\'', MHC+'.QC.map', '>', MHC+'.QC.dat'])
        # print(command)
        os.system(command)

        # command = ' '.join(['cut -f2', MHC+'.QC.map', '>', MHC+'.snps'])
        # print(command)
        # os.system(command)

        command = ' '.join(["cut -d ' ' -f1-5,7-", MHC+'.QC.ped', '>', MHC+'.QC.nopheno.ped'])
        # print(command)
        os.system(command)

        if not __save_intermediates:
            os.system(' '.join(['rm', MHC+'.QC.ped']))
            os.system(' '.join(['rm', MHC+'.QC.map']))


        idx_process += 1



    if __use_Multiple_Markers:

        # [3] Multiple Markers
        # [6] Multiple Markers + Adaptive Genetic Map
        __IMPUTE_OUT__ = HLA_Imputation(idx_process, MHC, _reference, _out, _hg, __overlap__, _window, _ne, _nthreads,
                                        _AdaptiveGeneticMap, _Average_Erate,
                                        LINKAGE2BEAGLE, BEAGLE2LINKAGE, BEAGLE2VCF, VCF2BEAGLE, PLINK, BEAGLE5,
                                        _answer=_answer, f_save_intermediates=__save_intermediates, _MultP=_MultP,
                                        _given_prephased=_given_prephased, f_prephasing=f_prephasing, f_measureAcc_v2=f_measureAcc_v2)


    elif not __use_Multiple_Markers:

        # [2] Plain
        # [4] Adaptive Genetic Map (HapMap)
        # [5] Adaptive Genetic Map
        __IMPUTE_OUT__ = HLA_Imputation_GM(idx_process, MHC, _reference, _out, _hg, _window, __overlap__[0], _ne, _nthreads,
                                           _AdaptiveGeneticMap, _Average_Erate, LINKAGE2BEAGLE, BEAGLE2LINKAGE, BEAGLE2VCF, VCF2BEAGLE,
                                           PLINK, BEAGLE5, _answer=_answer, f_save_intermediates=__save_intermediates,
                                           _HapMap_Map=_HapMap_Map, f_measureAcc_v2=f_measureAcc_v2)


    idx_process = __IMPUTE_OUT__.idx_process


    if CLEAN_UP:

        print("[{}] Clean Up.".format(idx_process))
        idx_process += 1


        if not __save_intermediates:
            os.system(' '.join(['rm', MHC + '.QC.nopheno.ped']))
            os.system(' '.join(['rm', MHC + '.QC.dat']))
            os.system(' '.join(['rm', MHC + '.QC.bed']))
            os.system(' '.join(['rm', MHC + '.QC.bim']))
            os.system(' '.join(['rm', MHC + '.QC.fam']))
            os.system(' '.join(['rm', MHC + '.QC.vcf']))
            os.system(' '.join(['rm', MHC + '.QC.GCchange.bgl']))
            os.system(' '.join(['rm', MHC + '.QC.GCchange.markers']))
            os.system(' '.join(['rm', MHC + '.QC.log']))
            os.system(' '.join(['rm', _out + '.bgl.log']))
            os.system(' '.join(['rm -rf', JAVATMP]))



        print("DONE!\n")





    return [__IMPUTE_OUT__.HLA_IMPUTATION_OUT, __IMPUTE_OUT__.accuracy]





if __name__ == "__main__":

    ########## < Main parser > ##########

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    ###########################################################################################

        CookHLA.py



    ###########################################################################################
                                     '''),
                                     add_help=False
                                     )

    ### Common arguments to share over the modules.

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="Show this help message and exit\n\n", action='help')

    parser.add_argument("--input", "-i", help="\nCommon prefix of Target Input files.\n\n", required=True)
    parser.add_argument("--reference", "-ref", help="\nPrefix of Reference files.\n\n", required=True)
    parser.add_argument("--out", "-o", help="\nOutput file name prefix\n\n", required=True)

    parser.add_argument("--genetic-map", "-gm", help="\nGenetic Map file.\n\n")
    parser.add_argument("--average-erate", "-ae", help="\nAverate error rate file.\n\n")
    # parser.add_argument("--use-multiple-markers", "-ml", help="\nUsing multiple markers.\n\n", action='store_true') => now default

    parser.add_argument("--prephasing", "-pr", help="\nUtilizing prephasing strategy.\n\n", action='store_true')

    parser.add_argument("--answer", "-an", help="\nAnswer file to calculate imputation accuracy.\n\n")

    parser.add_argument("--multiprocess", "-mp", help="\nSetting parallel multiprocessing.\n\n", type=int, choices=[2,3,4,5,6,7,8,9], nargs='?', default=1, const=3)

    parser.add_argument("--java-memory", "-mem", help="\nMemory requried for beagle(ex. 12g).\n\n", default="2g")

    parser.add_argument("--measureAcc_v2", "-macc_v2", help="\nCalculate accuracy with previous version module.\n\n", action='store_true')

    # Beagle5.1.
    parser.add_argument("--overlap", "-ol",
                        help="\n3 Overlap values(cM) for Beagle 5.1 implementation.\n\n", nargs=3, default=(0.5,1,1.5), type=float)
    parser.add_argument("--window", "-w",
                        help="\nWindow value(cM) for Beagle 5.1 implementation.\n\n", default=5, type=float)
    parser.add_argument("--effective-population-size", "-ne",
                        help="\nEffective population size value for Beagle 5.1 implementation.\n\n", default=1000000, type=int)
    parser.add_argument("--nthreads", "-nth",
                        help="\nThe number of theads to use in Beagle 5.1 implementation.\n\n", default=1, type=int)



    ##### < for Testing > #####

    # args = parser.parse_args(["--input", "data/Target/HM_CEU.FOUNDERS.filt",
    #                           "--out", "tests/_3_CookHLA/20190605_onlyAGM/_3_HM_CEU_T1DGC_REF",
    #                           "-ref", "data/HLA_PANEL/T1DGC/T1DGC_REF",
    #                           "-gm", "data/HLA_PANEL/Genetic_map/CEU_T1DGC.mach_step.avg.clpsB",
    #                           "-ae", "data/HLA_PANEL/Genetic_map/CEU_T1DGC.aver.erate",
    #                           "-an", "tests/HM_CEU_REF.bgl.phased.alleles.answer"])

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

    CookHLA_start = time()

    CookHLA(args.input, args.out, args.reference, "18", args.genetic_map, args.average_erate,
            _java_memory=args.java_memory, _MultP=args.multiprocess, _answer=args.answer, __use_Multiple_Markers=True,
            f_prephasing=args.prephasing, __overlap__=args.overlap, _ne=args.effective_population_size, f_measureAcc_v2=args.measureAcc_v2)

    CookHLA_end = time()

    CookHLA_time = (CookHLA_end - CookHLA_start)/60
    print(std_MAIN_PROCESS_NAME + "Total CookHLA time : {}(min)".format(CookHLA_time))
