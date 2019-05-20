#-*- coding: utf-8 -*-

import os, sys, re
from os.path import join
import argparse, textwrap
from src.HLA_Imputation import HLA_Imputation

########## < Core Varialbes > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]

TOLERATED_DIFF = 0.15


def CookHLA(_input, _out, _reference, _geneticMap, _average_erate, _hg, _java_memory='2g',
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
    BEAGLE4 = ' '.join(["java", '-Djava.io.tmpdir={}'.format(JAVATMP), "-Xmx{}".format(_java_memory), "-jar", _p_beagle4])
    LINKAGE2BEAGLE = ' '.join(["java", '-Djava.io.tmpdir={}'.format(JAVATMP), "-Xmx{}".format(_java_memory), "-jar", _p_linkage2beagle])
    BEAGLE2LINKAGE = ' '.join(["java", '-Djava.io.tmpdir={}'.format(JAVATMP), "-Xmx{}".format(_java_memory), "-jar", _p_beagle2linkage])
    BEAGLE2VCF = ' '.join(["java", '-Djava.io.tmpdir={}'.format(JAVATMP), "-Xmx{}".format(_java_memory), "-jar", _p_beagle2vcf])
    VCF2BEAGLE = ' '.join(["java", '-Djava.io.tmpdir={}'.format(JAVATMP), "-Xmx{}".format(_java_memory), "-jar", _p_vcf2beagle])

    MERGE = os.path.join(p_src, 'merge_tables.pl')
    PARSEDOSAGE = os.path.join(p_src, 'ParseDosage.csh')
    BGL2BED = os.path.join(p_src, 'Panel-BGL2BED.sh')



    ###### < Control Flags > ######

    EXTRACT_MHC = 1
    FLIP = 1
    CONVERT_IN = 1
    IMPUTE = 1
    CONVERT_OUT = 1
    CLEAN_UP = 0



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
            os.system(' '.join(['rm', MHC+'.{bed,bim,fam,log}']))
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
            os.system(' '.join(['rm', MHC+'.FLP.{bed,bim,fam,log}']))

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
            os.system(' '.join(['rm', MHC+'.QC.reorder.{bed,bim,fam,log}']))
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
            os.system(' '.join(['rm', MHC+'.QC.{ped,map}']))


        idx_process += 1


    ############################################################

    if CONVERT_IN:

        ### Testing Multiple Reference

        from src.HLA_Imputation import HLA_Imputation

        myImputation = HLA_Imputation(MHC, _reference, _out, _hg, LINKAGE2BEAGLE, BEAGLE2VCF, PLINK, BEAGLE4,
                                      __save_intermediates, idx_process)


        idx_process = myImputation.getIDX_PROCESS()

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
    parser.add_argument("-hg", help="\nHuman Genome version(ex. 18, 19, 38)\n\n", choices=["18", "19", "38"],
                        metavar="HG", default='18')


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

    CookHLA(args.input, args.out, args.reference, args.genetic_map, args.average_erate, args.hg,
            _java_memory=args.java_memory, __use_Multiple_Markers=args.use_multiple_markers)