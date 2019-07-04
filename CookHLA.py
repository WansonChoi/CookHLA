#-*- coding: utf-8 -*-

import os, sys, re
from os.path import join
import argparse, textwrap
from src.HLA_Imputation import HLA_Imputation
from src.HLA_MultipleRefs import HLA_MultipleRefs
from statistics import mean

########## < Core Varialbes > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]

TOLERATED_DIFF = 0.15


__ExonN__ = ['exon2', 'exon3', 'exon4']
# __ExonN__ = ['exon2']
# __ExonN__ = ['exon4']




def CookHLA(_input, _out, _reference, _hg='18', _geneticMap=None, _average_erate=None, _java_memory='2g', _MultP=1, _answer=None,
            __save_intermediates=False, __use_Multiple_Markers=False,
            _p_src="./src", _p_dependency="./dependency",):

    f_useGeneticMap = False

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
    CLEAN_UP = 0



    print(std_MAIN_PROCESS_NAME + "CookHLA : Performing HLA imputation for '{}'\n"
                                  "- Java memory = {}(Mb)".format(_input, _java_memory))




    ###### < Multiple Markers and Adaptive Genetic Map > ######

    # Multiple Markers?
    if __use_Multiple_Markers:
        print("- Using Multiple Markers.")

    # Adaptive Genetic Map?
    if _average_erate and _geneticMap:

        if os.path.exists(_average_erate) and os.path.exists(_geneticMap):
            f_useGeneticMap = True
            print("- Using Genetic Map : {}.".format(_geneticMap))
        else:
            if not os.path.exists(_average_erate):
                print(std_ERROR_MAIN_PROCESS_NAME + "The file ('{}') doesn't exist.\n"
                                                    "Please check '--average-erate/-ae' argument again.".format(_average_erate))
                sys.exit()

            if not os.path.exists(_geneticMap):
                print(std_ERROR_MAIN_PROCESS_NAME + "The file ('{}') doesn't exist.\n"
                                                    "Please check '--genetic-map/-gm' argument again.".format(_geneticMap))
                sys.exit()

    elif not (_average_erate or _geneticMap):

        f_useGeneticMap = False     # No using Adaptive Genetic Map.

    else:
        print(std_ERROR_MAIN_PROCESS_NAME + "Either arguments '--genetic-map(-gm)' or '--average-erate(-ae)' wasn't given.\n"
                                            "Please check whether both of them are given or not.")
        sys.exit()





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
            os.system(' '.join(['rm', MHC+'.{bed,bim,log}']))
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


    ### Testing Multiple Reference

    # Container for Imputation result.
    __IMPUTE_OUT__ = {}


    if __use_Multiple_Markers:

        if _MultP > 1:

            import multiprocessing as mp

            pool = mp.Pool(processes=_MultP)

            dict_Pool = {_exonN_: pool.apply_async(Main_Imputation, (idx_process, MHC, _reference, _out+'.{}'.format(_exonN_), _hg,
                                              LINKAGE2BEAGLE, BEAGLE2LINKAGE, BEAGLE2VCF, VCF2BEAGLE, PLINK, BEAGLE4,
                                              __save_intermediates,
                                              _average_erate, _geneticMap, f_useGeneticMap, __use_Multiple_Markers,
                                              _exonN_, _answer)) for _exonN_ in __ExonN__}

            pool.close()
            pool.join()


            __IMPUTE_OUT__ = {_exonN_: _OUT.get() for _exonN_, _OUT in dict_Pool.items()}

            # checking result
            # for k,v in __IMPUTE_OUT__.items():
            #     print("{} : {}".format(k, v))



        else:
            ### No Multiprocessing

            for _exonN_ in __ExonN__:

                ## Reference Panel for exon N.
                __ExonN_Refs__ = HLA_MultipleRefs(_exonN_, _reference, _out, _hg, BEAGLE2LINKAGE, PLINK, f_only_HLA=False)

                myImputation = HLA_Imputation(idx_process, MHC, __ExonN_Refs__.getOUTPUT(), _out+'.{}'.format(_exonN_), _hg,
                                              LINKAGE2BEAGLE, BEAGLE2LINKAGE, BEAGLE2VCF, VCF2BEAGLE,
                                              PLINK, BEAGLE4,
                                              __save_intermediates,
                                              _aver_erate=_average_erate, _Genetic_Map=_geneticMap,
                                              f_useGeneticMap = f_useGeneticMap, f_useMultipleMarkers=__use_Multiple_Markers,
                                              _exonN_=_exonN_, _answer=_answer)

                # Hardcoding for partial testing.
                # myImputation = HLA_Imputation(idx_process, MHC, 'tests/_3_CookHLA/20190530/T1DGC_REF.exon2',
                #                               _out + '.exon2', _hg,
                #                               LINKAGE2BEAGLE, BEAGLE2LINKAGE, BEAGLE2VCF, VCF2BEAGLE,
                #                               PLINK, BEAGLE4,
                #                               __save_intermediates,
                #                               _aver_erate=_average_erate, _Genetic_Map=_geneticMap,
                #                               f_useGeneticMap=f_useGeneticMap,
                #                               f_useMultipleMarkers=__use_Multiple_Markers,
                #                               _exonN_=_exonN_, _answer=_answer)

                __IMPUTE_OUT__ = {_exonN_: myImputation}

    else:

        __IMPUTE_OUT__ = HLA_Imputation(idx_process, MHC, _reference, _out, _hg,
                                      LINKAGE2BEAGLE, BEAGLE2LINKAGE, BEAGLE2VCF, VCF2BEAGLE, PLINK, BEAGLE4,
                                      __save_intermediates,
                                      _aver_erate=_average_erate, _Genetic_Map=_geneticMap,
                                      f_useGeneticMap = f_useGeneticMap, f_useMultipleMarkers=__use_Multiple_Markers,
                                      _answer=_answer)


    ############################################################



    ### Average of accuracies.

    if _answer:
        print("Raw Acc : \n{}".format(__IMPUTE_OUT__.accuracy))
        __MeanAccuracy__ = getMeanAccuracy(__IMPUTE_OUT__, __use_Multiple_Markers, f_useGeneticMap)



    if CLEAN_UP:

        print("[{}] Clean Up.".format(idx_process))


        if not __save_intermediates:
            os.system(' '.join(['rm', MHC + '.QC.nopheno.ped']))
            os.system(' '.join(['rm', MHC + '.QC.dat']))
            os.system(' '.join(['rm', MHC + '.QC.{bed,bim,fam,log}']))
            os.system(' '.join(['rm', _out + '.bgl.log']))



        print("DONE!\n")

        idx_process += 1





    return 0



def Main_Imputation(idx_process, MHC, _reference, _out, _hg,
                    LINKAGE2BEAGLE, BEAGLE2LINKAGE, BEAGLE2VCF, VCF2BEAGLE, PLINK, BEAGLE4,
                    __save_intermediates=False,
                    _average_erate=None, _geneticMap=None, f_useGeneticMap=False,
                    __use_Multiple_Markers=False, _exonN_=None, _answer=None):

    """
    Wrapper function for multiprocessing

    """

    ### Reference Panel for exon N.
    __ExonN_Refs__ = HLA_MultipleRefs(_exonN_, _reference, _out, _hg, BEAGLE2LINKAGE, PLINK)
    # __ExonN_Refs__ = 'tests/_3_CookHLA/20190523/T1DGC_REF.exon2'
    # print(__ExonN_Refs__.getOUTPUT())


    myImputation = HLA_Imputation(idx_process, MHC, __ExonN_Refs__.getOUTPUT(), _out, _hg,
                                  LINKAGE2BEAGLE, BEAGLE2LINKAGE, BEAGLE2VCF, VCF2BEAGLE, PLINK, BEAGLE4,
                                  __save_intermediates,
                                  _aver_erate=_average_erate, _Genetic_Map=_geneticMap,
                                  f_useGeneticMap = f_useGeneticMap, f_useMultipleMarkers=__use_Multiple_Markers,
                                  _exonN_=_exonN_, _answer=_answer)

    # Hardcoding for partial testing.
    # myImputation = HLA_Imputation(idx_process, MHC, 'tests/_3_CookHLA/20190523/T1DGC_REF.exon2', _out, _hg,
    #                               LINKAGE2BEAGLE, BEAGLE2LINKAGE, BEAGLE2VCF, VCF2BEAGLE,
    #                               PLINK, BEAGLE4,
    #                               __save_intermediates,
    #                               _aver_erate=_average_erate, _Genetic_Map=_geneticMap,
    #                               f_useGeneticMap=f_useGeneticMap, f_useMultipleMarkers=__use_Multiple_Markers,
    #                               _exonN_=_exonN_, _answer=_answer)



    return myImputation



def getMeanAccuracy(__IMPUTE_OUT__, __f_Multiple_Markers, __f_useGeneticMap):


    __RETURN__ = {_hla: None for _hla in HLA_names}



    if __f_useGeneticMap:

        __overlap__ = [3000, 4000, 5000]

        if __f_Multiple_Markers:
            for _hla in HLA_names:
                __RETURN__[_hla] = mean([__IMPUTE_OUT__[exonN].accuracy[_overlap_]['4D'][_hla] for exonN in __ExonN__ for _overlap_ in __overlap__])
        else:

            for _hla in HLA_names:
                __RETURN__[_hla] = mean([__IMPUTE_OUT__.accuracy[_overlap_]['4D'][_hla] for _overlap_ in __overlap__])


    else:

        if __f_Multiple_Markers:

            for _hla in HLA_names:

                __RETURN__[_hla] = mean([__IMPUTE_OUT__[exonN].accuracy['4D'][_hla] for exonN in __ExonN__])

        else:
            # Neither Adaptive Genetic map and Multiple Map
            __RETURN__ = __IMPUTE_OUT__.accuracy['4D']

    print(__RETURN__)

    return __RETURN__



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

    parser.add_argument("--answer", "-an", help="\nAnswer file to calculate imputation accuracy.\n\n")

    parser.add_argument("--multiprocess", "-mp", help="\nSetting parallel multiprocessing.\n\n", type=int, choices=[2,3], nargs='?', default=1, const=3)

    parser.add_argument("--java-memory", "-mem", help="\nMemory requried for beagle(ex. 12g).\n\n", default="2g")





    ##### < for Testing > #####

    # args = parser.parse_args(["--input", "data/Target/HM_CEU.FOUNDERS.filt",
    #                           "--out", "tests/_3_CookHLA/20190605_onlyAGM/_3_HM_CEU_T1DGC_REF",
    #                           "-ref", "data/HLA_PANEL/T1DGC/T1DGC_REF",
    #                           "-gm", "data/HLA_PANEL/Genetic_map/CEU_T1DGC.mach_step.avg.clpsB",
    #                           "-ae", "data/HLA_PANEL/Genetic_map/CEU_T1DGC.aver.erate",
    #                           "-an", "tests/HM_CEU_REF.bgl.phased.alleles.answer"])




    ##### < for Publish > #####
    args = parser.parse_args()
    print(args)

    CookHLA(args.input, args.out, args.reference, args.hg, args.genetic_map, args.average_erate,
            _java_memory=args.java_memory, __use_Multiple_Markers=args.use_multiple_markers, _MultP=args.multiprocess,
            _answer=args.answer)