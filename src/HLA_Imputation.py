# -*- coding: utf-8 -*-

import os, sys, re
import subprocess
from os.path import join, exists
import multiprocessing as mp
from time import time

# Both
# from src.GC_tricked_bgl2ori_bgl import GCtricedBGL2OriginalBGL
from src.RUN_Bash import RUN_Bash
from src.measureAccuracy import measureAccuracy

# Multiple Markers
# from src.redefineBPv1BH import redefineBP
from src.HLA_MultipleRefs import HLA_MultipleRefs

# HLA genotype calling
HLA_genotype_call_prephasing = 'src/9accuracy_pre.v2.csh' # Practically deprecated (2020.09.12.)
HLA_genotype_call_noprephasing = 'src/9accuracy_no_CI.v2.csh' # Updated to print confident score (2020.09.12.)

# Defined Error
from src.CookHLAError import CookHLAImputationError, CookHLAHLATypeCallError

# measureAcc_v3.5
from measureAcc.measureAccuracy import CookHLA_measureAcc
from measureAcc.src.ALLELES2HPED import ALLELES2HPED


########## < Core Varialbes > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
HLA_names_gen = ["A", "C", "B", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]
isClassI = {"A": True, "B": True, "C": True, "DPA1": False, "DPB1": False, "DQA1": False, "DQB1": False, "DRB1": False}


__EXON__ = ['exon2', 'exon3', 'exon4']
# __EXON__ = ['exon2']

__overlap__ = [3000, 4000, 5000]
# __overlap__ = [3000]



class HLA_Imputation(object):

    def __init__(self, idx_process, MHC, _reference, _out, _hg, _nthreads, _AdaptiveGeneticMap, _Average_Erate, _LINKAGE2BEAGLE,
                 _BEAGLE2LINKAGE, _BEAGLE2VCF, _VCF2BEAGLE, _PLINK, _BEAGLE4, _CSH, _answer=None, f_save_intermediates=False,
                 _MultP=1, _given_prephased=None, f_prephasing=False, f_remove_raw_IMP_results=False, f_measureAcc_v2=False):

        ### General
        self.idx_process = idx_process
        self.__save_intermediates = f_save_intermediates

        self.FLAG_AdaptiveGeneticMap = _AdaptiveGeneticMap and _Average_Erate # (***) Deciding whether to use Adaptive genetic map or not.


        # Prefixes
        self.OUTPUT_dir = os.path.dirname(_out)
        self.OUTPUT_dir_ref = join(self.OUTPUT_dir, os.path.basename(_reference))
        self.OUTPUT_dir_GM = join(self.OUTPUT_dir, os.path.basename(_AdaptiveGeneticMap)) if self.FLAG_AdaptiveGeneticMap else None


        # Result
        self.Exon234_Panel = None
        self.dict_ExonN_Panel = {_exonN: None for _exonN in __EXON__}
        self.dict_ExonN_AGM = {_exonN: None for _exonN in __EXON__}
        self.dict_IMP_Result = {_exonN: {_overlap: None for _overlap in __overlap__} for _exonN in __EXON__}
        self.accuracy = None
        self.HLA_IMPUTATION_OUT = None

        self.dict_DOUBLED_PHASED_RESULT = {_exonN: None for _exonN in __EXON__}
        self.dict_REF_PHASED_VCF = {_exonN: None for _exonN in __EXON__}

        # Dependencies
        self.LINKAGE2BEAGLE = _LINKAGE2BEAGLE
        self.BEAGLE2LINKAGE = _BEAGLE2LINKAGE
        self.BEAGLE2VCF = _BEAGLE2VCF
        self.VCF2BEAGLE = _VCF2BEAGLE
        self.PLINK = _PLINK
        self.BEAGLE4 = _BEAGLE4

        # created in 'CONVERT_IN'
        # self.refined_REF_markers = None # used in 'CONVERT_OUT'
        # self.refined_Genetic_Map = None # used in 'IMPUTE'
        # self.GCchangeBGL = None # used in 'CONVERT_OUT'

        # Adaptive Genetic Map
        self.__AGM__ = _AdaptiveGeneticMap if self.FLAG_AdaptiveGeneticMap else None
        self.__AVER__ = _Average_Erate if self.FLAG_AdaptiveGeneticMap else None





        ###### < Reference panel for Exon 2, 3, 4 > ######

        multiple_panels = HLA_MultipleRefs(_reference, self.OUTPUT_dir_ref, _hg,
                                           self.BEAGLE2LINKAGE, self.BEAGLE2VCF, self.PLINK,
                                           _MultP=_MultP,
                                           __AGM__=self.__AGM__, _out_AGM=self.OUTPUT_dir_GM)

        self.dict_ExonN_Panel = multiple_panels.ExonN_Panel
        self.Exon234_Panel = multiple_panels.EXON234_Panel

        self.dict_ExonN_AGM = multiple_panels.ExonN_AGM if self.FLAG_AdaptiveGeneticMap else {_exonN: None for _exonN in __EXON__}

        # [Temporary Hard-coding]
        # self.Exon234_Panel = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190716_BOTH/T1DGC_REF.exon234'
        #
        # self.dict_ExonN_Panel['exon2'] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190716_BOTH/T1DGC_REF.exon2'
        # self.dict_ExonN_Panel['exon3'] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190716_BOTH/T1DGC_REF.exon3'
        # self.dict_ExonN_Panel['exon4'] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190716_BOTH/T1DGC_REF.exon4'
        #
        # self.dict_ExonN_AGM['exon2'] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190716_BOTH/CEU_T1DGC.mach_step.avg.clpsB.exon2.txt'
        # self.dict_ExonN_AGM['exon3'] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190716_BOTH/CEU_T1DGC.mach_step.avg.clpsB.exon3.txt'
        # self.dict_ExonN_AGM['exon4'] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190716_BOTH/CEU_T1DGC.mach_step.avg.clpsB.exon4.txt'



        ###### < Main - 'CONVERT_IN', 'IMPUTE', 'CONVERT_OUT' > ######


        ### (1) CONVERT_IN

        IMPUTATION_INPUT = self.CONVERT_IN(MHC, self.Exon234_Panel, _out, _hg, _given_prephased=_given_prephased,
                                           f_prephasing=f_prephasing)
        # Only one time of pre-phasing with Exon234 reference panel.



        ### (2) Imputation

        if _MultP == 1:

            imputation_serial_start = time()

            ## Serial implementation of main.
            for _exonN in __EXON__:
                for _overlap in __overlap__:

                    self.dict_IMP_Result[_exonN][_overlap] = \
                        self.IMPUTE(MHC, _out, IMPUTATION_INPUT, self.dict_ExonN_Panel[_exonN] + '.phased.vcf',
                                    _overlap, _exonN, _nthreads, self.__AVER__, self.dict_ExonN_AGM[_exonN], f_prephasing=f_prephasing)

            imputation_serial_end = time()

            imputation_serial_time = (imputation_serial_end - imputation_serial_start)/60
            print(std_MAIN_PROCESS_NAME+"Total imputation time of Serial implementation: {}(min)\n".format(imputation_serial_time))

        else:

            ## Parallel implementation of main.

            imputation_parallel_start = time()

            pool = mp.Pool(processes=_MultP if _MultP <= 9 else 9)

            dict_Pool = {_exonN: {_overlap: pool.apply_async(self.IMPUTE, (MHC, _out, IMPUTATION_INPUT, self.dict_ExonN_Panel[_exonN] + '.phased.vcf', _overlap, _exonN, _nthreads, self.__AVER__, self.dict_ExonN_AGM[_exonN], f_prephasing))
                                  for _overlap in __overlap__}
                         for _exonN in __EXON__}

            pool.close()
            pool.join()


            for _exonN in __EXON__:
                for _overlap in __overlap__:
                    self.dict_IMP_Result[_exonN][_overlap] = dict_Pool[_exonN][_overlap].get()

            imputation_parallel_end = time()

            imputation_parallel_time = (imputation_parallel_end - imputation_parallel_start)/60
            print(std_MAIN_PROCESS_NAME + "Total imputation time of Parallel implementation (with {} core(s)): {}(min)\n".format(_MultP, imputation_parallel_time))

        self.idx_process += 1


        # [Temporary Hard-coding]
        # self.dict_IMP_Result['exon2'][3000] = '/home/wanson/Git_Projects/CookHLA/tests/_3_CookHLA/20190731_MM_AGM/HM_CEU_T1DGC_REF.MM.AGM.noprephasing.exon2.overlap3000.MHC.QC.double.imputation_out.vcf'
        # self.dict_IMP_Result['exon2'][4000] = '/home/wanson/Git_Projects/CookHLA/tests/_3_CookHLA/20190731_MM_AGM/HM_CEU_T1DGC_REF.MM.AGM.noprephasing.exon2.overlap4000.MHC.QC.double.imputation_out.vcf'
        # self.dict_IMP_Result['exon2'][5000] = '/home/wanson/Git_Projects/CookHLA/tests/_3_CookHLA/20190731_MM_AGM/HM_CEU_T1DGC_REF.MM.AGM.noprephasing.exon2.overlap5000.MHC.QC.double.imputation_out.vcf'
        # self.dict_IMP_Result['exon3'][3000] = '/home/wanson/Git_Projects/CookHLA/tests/_3_CookHLA/20190731_MM_AGM/HM_CEU_T1DGC_REF.MM.AGM.noprephasing.exon3.overlap3000.MHC.QC.double.imputation_out.vcf'
        # self.dict_IMP_Result['exon3'][4000] = '/home/wanson/Git_Projects/CookHLA/tests/_3_CookHLA/20190731_MM_AGM/HM_CEU_T1DGC_REF.MM.AGM.noprephasing.exon3.overlap4000.MHC.QC.double.imputation_out.vcf'
        # self.dict_IMP_Result['exon3'][5000] = '/home/wanson/Git_Projects/CookHLA/tests/_3_CookHLA/20190731_MM_AGM/HM_CEU_T1DGC_REF.MM.AGM.noprephasing.exon3.overlap5000.MHC.QC.double.imputation_out.vcf'
        # self.dict_IMP_Result['exon4'][3000] = '/home/wanson/Git_Projects/CookHLA/tests/_3_CookHLA/20190731_MM_AGM/HM_CEU_T1DGC_REF.MM.AGM.noprephasing.exon4.overlap3000.MHC.QC.double.imputation_out.vcf'
        # self.dict_IMP_Result['exon4'][4000] = '/home/wanson/Git_Projects/CookHLA/tests/_3_CookHLA/20190731_MM_AGM/HM_CEU_T1DGC_REF.MM.AGM.noprephasing.exon4.overlap4000.MHC.QC.double.imputation_out.vcf'
        # self.dict_IMP_Result['exon4'][5000] = '/home/wanson/Git_Projects/CookHLA/tests/_3_CookHLA/20190731_MM_AGM/HM_CEU_T1DGC_REF.MM.AGM.noprephasing.exon4.overlap5000.MHC.QC.double.imputation_out.vcf'



        ### (3) CONVERT_OUT

        self.HLA_IMPUTATION_OUT = self.CONVERT_OUT(self.dict_IMP_Result, MHC + '.HLA_IMPUTATION_OUT', _CSH, f_prephasing=f_prephasing)
        print(std_MAIN_PROCESS_NAME+'IMPUTATION_OUT:\n{}'.format(self.HLA_IMPUTATION_OUT))


        ## Acquiring accuracy

        if bool(_answer):

            print(std_MAIN_PROCESS_NAME + "Calculating accuracy of each HLA gene. (answer: '{}')".format(_answer))

            if not os.path.exists(_answer):
                print(std_WARNING_MAIN_PROCESS_NAME + "Given answer file doesn't exist. Please check '--answer/-an' argument again.\n"
                                                    "Skipping calculating imputation accuracy.")
            elif os.path.getsize(_answer) == 0:
                print(std_WARNING_MAIN_PROCESS_NAME + "Given answer file doesn't have any content. Please check '--answer/-an' argument again.\n"
                                                    "Skipping calculating imputation accuracy.")

            else:
                if f_measureAcc_v2:
                    # measureAcc_v2
                    self.accuracy = measureAccuracy(_answer, self.HLA_IMPUTATION_OUT, 'all',
                                                        outfile=self.HLA_IMPUTATION_OUT + '.accuracy', __only4digits=True)

                else:
                    # measureAcc_v3.5
                    measureAcc_start = time()

                    t = CookHLA_measureAcc(_answer, self.HLA_IMPUTATION_OUT, self.HLA_IMPUTATION_OUT)
                    self.accuracy = t.accuracy

                    measureAcc_end = time()

                    measureAcc_time = (measureAcc_end - measureAcc_start)/60
                    print("\nAccuracy : {}".format(self.accuracy))
                    print("measureAccuracy time: {}(min)\n".format(measureAcc_time))



        ### General Removal
        if not self.__save_intermediates:

            # 'Exon234 panel'
            RUN_Bash('rm {}'.format(self.Exon234_Panel+'.bed'))
            RUN_Bash('rm {}'.format(self.Exon234_Panel+'.bim'))
            RUN_Bash('rm {}'.format(self.Exon234_Panel+'.fam'))
            RUN_Bash('rm {}'.format(self.Exon234_Panel+'.FRQ.frq'))
            RUN_Bash('rm {}'.format(self.Exon234_Panel+'.markers'))
            RUN_Bash('rm {}'.format(self.Exon234_Panel+'.bgl.phased'))
            RUN_Bash('rm {}'.format(self.Exon234_Panel+'.GCchange.markers'))
            RUN_Bash('rm {}'.format(self.Exon234_Panel+'.GCchange.bgl.phased'))
            RUN_Bash('rm {}'.format(self.Exon234_Panel+'.phased.vcf'))
            RUN_Bash('rm {}'.format(self.Exon234_Panel+'.refined.markers')) # only in Exon234 panel


            # 'Exon 2,3,4 panel'
            for _exonN in __EXON__:
                RUN_Bash('rm {}'.format(self.dict_ExonN_Panel[_exonN] + '.bed'))
                RUN_Bash('rm {}'.format(self.dict_ExonN_Panel[_exonN] + '.bim'))
                RUN_Bash('rm {}'.format(self.dict_ExonN_Panel[_exonN] + '.fam'))
                RUN_Bash('rm {}'.format(self.dict_ExonN_Panel[_exonN] + '.FRQ.frq'))
                RUN_Bash('rm {}'.format(self.dict_ExonN_Panel[_exonN] + '.markers'))
                RUN_Bash('rm {}'.format(self.dict_ExonN_Panel[_exonN] + '.bgl.phased'))
                RUN_Bash('rm {}'.format(self.dict_ExonN_Panel[_exonN] + '.GCchange.markers'))
                RUN_Bash('rm {}'.format(self.dict_ExonN_Panel[_exonN] + '.GCchange.bgl.phased'))
                RUN_Bash('rm {}'.format(self.dict_ExonN_Panel[_exonN] + '.phased.vcf'))


            # 'Exon 2,3,4 AGM'
            RUN_Bash('rm {}'.format(multiple_panels.EXON234_AGM))
            for _exonN in __EXON__:
                RUN_Bash('rm {}'.format(self.dict_ExonN_AGM[_exonN]))


            # # 'CONVERT_IN'
            # RUN_Bash('rm {}'.format(MHC + '.QC.nopheno.ped'))
            # RUN_Bash('rm {}'.format(MHC + '.QC.dat'))


            # 'CONVERT_OUT'
            for _exonN in __EXON__:
                for _overlap in __overlap__:
                    for _hla in HLA_names:
                        RUN_Bash('rm {}'.format(self.dict_IMP_Result[_exonN][_overlap]+'.HLA_{}'.format(_hla)))
                        if f_remove_raw_IMP_results:
                            RUN_Bash('rm {}'.format(self.dict_IMP_Result[_exonN][_overlap]))
                            RUN_Bash('rm {}'.format(self.dict_IMP_Result[_exonN][_overlap].rstrip('.vcf') + '.log'))



    def CONVERT_IN(self, MHC, _reference, _out, _hg, _given_prephased=None, f_prephasing=False):


        if _given_prephased and f_prephasing:

            print("(Test Purpose) Given pre-phased result will be used. ('{}')".format(_given_prephased))

            ############### < Multiple Markers > ###############

            ### Phasing & Doubling (only on Target Sample.)

            # Phasing (If previously prephased result is given, then the process to make new phased result will be skipped.
            PHASED_RESULT = _given_prephased

            # Doubling
            DOUBLED_PHASED_RESULT = self.Doubling(MHC, PHASED_RESULT.rstrip('.vcf'))

            return DOUBLED_PHASED_RESULT




        OUTPUT_dir_Exon234_ref = join(self.OUTPUT_dir, os.path.basename(_reference))


        print("[{}] Converting data to beagle format.".format(self.idx_process))
        self.idx_process += 1

        RUN_Bash(self.LINKAGE2BEAGLE + ' pedigree={} data={} beagle={} standard=true > {}'.format(
            MHC + '.QC.nopheno.ped', MHC + '.QC.dat', MHC + '.QC.bgl', _out + '.bgl.log'))

        # if not self.__save_intermediates:
        #     os.system('rm {}'.format(MHC + '.QC.nopheno.ped'))
        #     os.system('rm {}'.format(MHC + '.QC.dat'))
        #     os.system('rm {}'.format(_out + '.bgl.log'))




        ### Converting data to reference_markers_Position (Dispersing same genomic position of some markers.)

        # RefinedMarkers = redefineBP(_reference + '.markers', OUTPUT_dir_Exon234_ref + '.refined.markers')
        RUN_Bash('cp {} {}'.format(_reference + '.markers', OUTPUT_dir_Exon234_ref + '.refined.markers'))
        RefinedMarkers = OUTPUT_dir_Exon234_ref + '.refined.markers'
        # self.refined_REF_markers = RefinedMarkers # => This will be used in 'CONVERT_OUT'.




        ### Converting data to target_markers_Position and extract not_including snp.

        RUN_Bash('awk \'{print $2" "$4" "$5" "$6}\' %s > %s' % (MHC + '.QC.bim', MHC + '.QC.markers'))

        RUN_Bash('Rscript src/excluding_snp_and_refine_target_position-v1COOK02222017.R {} {} {}'.format(
            MHC + '.QC.markers', RefinedMarkers, MHC + '.QC.pre.markers'
        ))
        if not self.__save_intermediates:
            os.system(' '.join(['rm', MHC + '.QC.markers']))

        RUN_Bash('mv {} {}'.format(MHC + '.QC.bgl', MHC + '.QC.pre.bgl.phased'))

        RUN_Bash("awk '{print $1}' %s > %s" % (MHC + '.QC.pre.markers', join(self.OUTPUT_dir, 'selected_snp.txt')))

        from src.Panel_subset import Panel_Subset
        qc_refined = Panel_Subset(MHC + '.QC.pre', 'all', join(self.OUTPUT_dir, 'selected_snp.txt'),
                                  MHC + '.QC.refined')

        if not self.__save_intermediates:
            RUN_Bash('rm {}'.format(MHC + '.QC.pre.bgl.phased'))
            RUN_Bash('rm {}'.format(MHC + '.QC.pre.markers'))
            RUN_Bash('rm {}'.format(join(self.OUTPUT_dir, 'selected_snp.txt')))




        ### Converting data to GC_change_beagle format.

        from src.bgl2GC_trick_bgl import Bgl2GC

        # target
        [GCchangeBGL, GCchangeMarkers] = Bgl2GC(MHC + '.QC.refined.bgl.phased',
                                                MHC + '.QC.refined.markers',
                                                MHC + '.QC.GCchange.bgl',
                                                MHC + '.QC.GCchange.markers')

        self.GCchangeBGL = GCchangeBGL  # it will be used in 'CONVERT_OUT' with Genetic Map

        # print("<Target GCchanged bgl and marker file>\n"
        #       "bgl : {}\n"
        #       "markers : {}".format(GCchangeBGL, GCchangeMarkers))

        # reference
        [GCchangeBGL_REF, GCchangeMarkers_REF] = Bgl2GC(_reference + '.bgl.phased', RefinedMarkers,
                                                        OUTPUT_dir_Exon234_ref + '.GCchange.bgl.phased',
                                                        OUTPUT_dir_Exon234_ref + '.GCchange.markers')
        # print("<Reference GCchanged bgl and marker file>\n"
        #       "bgl : {}\n"
        #       "markers : {}".format(GCchangeBGL_REF, GCchangeMarkers_REF))

        if not self.__save_intermediates:
            RUN_Bash('rm {}'.format(MHC + '.QC.refined.bgl.phased'))
            RUN_Bash('rm {}'.format(MHC + '.QC.refined.markers'))
            # RUN_Bash('rm {}'.format(RefinedMarkers))

            # os.system(' '.join(['rm', RefinedMarkers])) # => This will be used in 'CONVERT_OUT" when not using Multiple Markers.




        ### Converting data to vcf_format

        # target
        RUN_Bash(self.BEAGLE2VCF + ' 6 {} {} 0 > {}'.format(GCchangeMarkers, GCchangeBGL, MHC + '.QC.vcf'))

        MHC_QC_VCF_exonN = MHC + '.QC.vcf'


        # reference
        RUN_Bash(self.BEAGLE2VCF + ' 6 {} {} 0 > {}'.format(GCchangeMarkers_REF, GCchangeBGL_REF,
                                                            OUTPUT_dir_Exon234_ref + '.vcf'))

        reference_vcf = OUTPUT_dir_Exon234_ref + '.vcf'




        ### Converting data to reference_phased

        RUN_Bash('sed "s%/%|%g" {} > {}'.format(reference_vcf, OUTPUT_dir_Exon234_ref + '.phased.vcf'))

        REF_PHASED_VCF = OUTPUT_dir_Exon234_ref + '.phased.vcf'

        if not self.__save_intermediates:
            RUN_Bash('rm {}'.format(reference_vcf))

            # # if self.f_useMultipleMarkers:
            # if not self.f_useGeneticMap:
            #     os.system(' '.join(['rm {}'.format(GCchangeBGL)])) # 'GCchangeBGL' will be used in 'CONVERT_OUT'
            #     os.system(' '.join(['rm {}'.format(GCchangeMarkers_REF)]))  # 'GCchangeMarkers_REF' will be used in 'CONVERT_OUT'
            #     os.system(' '.join(['rm {}'.format(GCchangeMarkers)]))
            #     os.system(' '.join(['rm {}'.format(GCchangeBGL_REF)]))

        """
        (1) `MHC_QC_VCF_exonN` := MHC + '.QC.vcf',
        (2) `REF_PHASED_VCF` := OUTPUT_dir_Exon234_ref + '.phased.vcf'

        These two files are to be passed into Beagle phasing;
        """


        if f_prephasing:

            ############### < Multiple Markers > ###############

            ### Phasing & Doubling (only on Target Sample.)

            # Phasing
            PHASED_RESULT = self.Phasing(MHC, MHC_QC_VCF_exonN, REF_PHASED_VCF)

            # [Temporary Hardcoding for Phased Result]
            # PHASED_RESULT = "/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190716_BOTH/_3_HM_CEU_T1DGC_REF.MHC.QC.phasing_out_not_double"
            # print("[Temporary Hardcoding]Phased Result:\n{}".format(PHASED_RESULT))

            # Doubling
            DOUBLED_PHASED_RESULT = self.Doubling(MHC, PHASED_RESULT)



            return DOUBLED_PHASED_RESULT

        else:

            return MHC_QC_VCF_exonN



    def IMPUTE(self, MHC, _out, _IMPUTATION_INPUT, _REF_PHASED_VCF, _overlap, _exonN, _nthreads, _aver_erate, _Refined_Genetic_Map, f_prephasing=False):

        if os.path.getsize(_IMPUTATION_INPUT) == 0:
            print(std_ERROR_MAIN_PROCESS_NAME + "Input file for imputation('{}') contains nothing. Please check it again.".format(_IMPUTATION_INPUT))
            sys.exit()


        # print("[{}] Performing HLA imputation (see {}.MHC.QC.imputation_out.log for progress).".format(self.idx_process, _out))
        print("\n[{}] Performing HLA imputation({} / overlap:{}).".format(self.idx_process, _exonN, _overlap))
        # self.idx_process += 1



        raw_HLA_IMPUTATION_OUT = MHC + ('.QC.{}.{}.doubled.raw_imputation_out'.format(_exonN, _overlap) if f_prephasing else '.QC.{}.{}.raw_imputation_out'.format(_exonN, _overlap))




        if self.FLAG_AdaptiveGeneticMap: # With Adatpive Genetic Map

            """
            ### MM + AGM
            
            # prephasing
            java -jar beagle4.jar gt=$MHC.QC.phasing_out_double.vcf ref=$REFERENCE.phased.vcf out=$MHC.QC.double.imputation_out impute=true lowmem=true gprobs=true ne=10000 overlap=${OVERLAP} err=$aver_erate map=$geneticMap.refined.map
            
            # No-prephasing
            java -jar beagle4.jar gt=$MHC.QC.vcf                    ref=$REFERENCE.phased.vcf out=$MHC.QC.double.imputation_out impute=true lowmem=true gprobs=true ne=10000 overlap=${OVERLAP} err=$aver_erate map=$geneticMap.refined.map
            
            """


            # aver_erate
            with open(_aver_erate, 'r') as f:
                aver_erate = f.readline().rstrip('\n')



            command = '{} gt={} ref={} out={} impute=true lowmem=true gprobs=true ne=10000 overlap={} err={} map={} nthreads={}'.format(
                self.BEAGLE4, _IMPUTATION_INPUT, _REF_PHASED_VCF, raw_HLA_IMPUTATION_OUT, _overlap, aver_erate, _Refined_Genetic_Map, _nthreads)
            # print(command)

            try:
                f_log = open(raw_HLA_IMPUTATION_OUT+'.log', 'w')

                imputation_start = time()
                subprocess.run(re.split('\s+', command), check=True, stdout=f_log, stderr=f_log)
                imputation_end = time()

            except subprocess.CalledProcessError:
                raise CookHLAImputationError(std_ERROR_MAIN_PROCESS_NAME + "Imputation({} / overlap:{}) failed.\n".format(_exonN, _overlap))
                # sys.stderr.write(std_ERROR_MAIN_PROCESS_NAME + "Imputation({} / overlap:{}) failed.\n".format(_exonN, _overlap))
                # return -1
            else:
                # print(std_MAIN_PROCESS_NAME+"Imputation({} / overlap:{}) done.".format(_exonN, _overlap))
                # os.system("rm {}".format(raw_HLA_IMPUTATION_OUT+'.err.log'))
                f_log.close()

                imputation_time = (imputation_end - imputation_start)/60
                sys.stdout.write("Imputation({} / overlap:{}) time: {}(min)\n".format(_exonN, _overlap, imputation_time))



        else: # Without Adaptive Genetic Map

            """
            ### MM
            
            # prephasing
            java -jar beagle4.jar gt=$MHC.QC.phasing_out_double.vcf ref=$REFERENCE.phased.vcf out=$MHC.QC.double.imputation_out impute=true lowmem=true overlap=$OVERLAP gprobs=true
            
            # No-prephasing
            java -jar beagle4.jar gt=$MHC.QC.vcf                    ref=$REFERENCE.phased.vcf out=$MHC.QC.double.imputation_out impute=true lowmem=true overlap=$OVERLAP gprobs=true

            """


            command = '{} gt={} ref={} out={} impute=true lowmem=true overlap={} gprobs=true nthreads={}'.format(
                self.BEAGLE4, _IMPUTATION_INPUT, _REF_PHASED_VCF, raw_HLA_IMPUTATION_OUT, _overlap, _nthreads)
            # print(command)

            try:
                f_log = open(raw_HLA_IMPUTATION_OUT+'.log', 'w')

                imputation_start = time()
                subprocess.run(re.split('\s+', command), check=True, stdout=f_log, stderr=f_log)
                imputation_end = time()

            except subprocess.CalledProcessError:
                raise CookHLAImputationError(std_ERROR_MAIN_PROCESS_NAME + "Imputation({} / overlap:{}) failed.\n".format(_exonN, _overlap))

            else:
                # print(std_MAIN_PROCESS_NAME+"Imputation({} / overlap:{}) done.".format(_exonN, _overlap))
                f_log.close()

                imputation_time = (imputation_end - imputation_start)/60
                sys.stdout.write("Imputation({} / overlap:{}) time: {}(min)\n".format(_exonN, _overlap, imputation_time))



        RUN_Bash('gzip -d -f {}.vcf.gz'.format(raw_HLA_IMPUTATION_OUT))

        return raw_HLA_IMPUTATION_OUT + '.vcf'



    def CONVERT_OUT(self, _raw_IMP_Result, _out, _CSH, f_prephasing=False):


        print("[{}] Converting out imputation result(s).".format(self.idx_process, _out))
        self.idx_process += 1


        to_args = ' '.join([_raw_IMP_Result[_exon][_overlap] for _exon in __EXON__ for _overlap in __overlap__])

        if f_prephasing:
            # Practically deprecated (2020.09.12.)
            command = '{} {} {} {}'.format(_CSH, HLA_genotype_call_prephasing, to_args, _out)
        else:
            # Updated to print confident score (2020.09.12.)
            command = '{} {} {} {}'.format(_CSH, HLA_genotype_call_noprephasing, to_args, _out)
        # print(command)

        try:
            f_log = open(_out + '.HLATypeCall.log', 'w')
            subprocess.run(command.split(' '), check=True, stdout=f_log, stderr=f_log)

        except subprocess.CalledProcessError:
            raise CookHLAHLATypeCallError("HLA type Calling failed.\n")

        else:

            # Merge HLATypeCall results of each HLA gene
            with open(_out+'.alleles', 'w') as f_HLAtypeCall_merged:
                for hla in HLA_names:

                    HLAtypeCall_each = _raw_IMP_Result['exon2'][3000]+'.HLA_{}.alleles'.format(hla)

                    if exists(HLAtypeCall_each):
                        f_HLAtypeCall_each = open(HLAtypeCall_each, 'r')
                        f_HLAtypeCall_merged.writelines(f_HLAtypeCall_each.readlines())
                        f_HLAtypeCall_each.close()

                        # '*.alleles'
                        os.system("rm {}".format(HLAtypeCall_each))


            if exists(_out + '.alleles') and os.path.getsize(_out + '.alleles') > 0:

                # Generate HPED, too. (2020. 09. 12.)
                ALLELES2HPED(_out + '.alleles', _out)

                return _out + '.alleles'
            else:
                print(std_ERROR_MAIN_PROCESS_NAME + "Failed to perform final HLA genotype calling.")
                return '-1'



    def Phasing(self, MHC, MHC_QC_VCF, REF_PHASED_VCF):

        ### Phasing (only on Target Sample.)
        print("[{}] Performing pre-phasing".format(self.idx_process))
        self.idx_process += 1

        command = ' '.join([self.BEAGLE4, 'gt={} ref={} out={} impute=false > {}'.format(MHC_QC_VCF, REF_PHASED_VCF,
                                                                                         MHC + '.QC.phasing_out_not_double',
                                                                                         MHC + '.QC.phasing_out_not_double.vcf.log')])
        # print(command)

        if not os.system(command):
            if not self.__save_intermediates:
                os.system(' '.join(['rm', MHC_QC_VCF]))
                os.system(' '.join(['rm', MHC + '.QC.phasing_out_not_double.vcf.log']))
                os.system(' '.join(['rm', MHC + '.QC.phasing_out_not_double.log']))
        else:
            print(std_ERROR_MAIN_PROCESS_NAME + "Failed to Phasing.\n"
                                                "Please check log file('{}')".format(MHC + '.QC.phasing_out_not_double.vcf.log'))
            sys.exit()

        return MHC + '.QC.phasing_out_not_double'



    def Doubling(self, MHC, PHASED_RESULT):

        ### Target data doubling step.
        print("[{}] Performing Doubling".format(self.idx_process))
        self.idx_process += 1

        RUN_Bash('gzip -d -f {}'.format(PHASED_RESULT + '.vcf.gz'))

        RUN_Bash('grep ^## {} > {}'.format(PHASED_RESULT + '.vcf', PHASED_RESULT + '.vcf.header1'))  # Header part with '##'
        RUN_Bash('grep -v ^## {} | head -n 1 > {}'.format(PHASED_RESULT + '.vcf', PHASED_RESULT + '.vcf.header2'))  # Header part with '#'
        RUN_Bash('grep -v ^# {} > {}'.format(PHASED_RESULT + '.vcf', PHASED_RESULT + '.vcf.body'))  # Body part without '#' or '##'

        RUN_Bash('sed "s%#%%" {} > {}'.format(PHASED_RESULT + '.vcf.header2', PHASED_RESULT + '.vcf.noshop.header2'))
        RUN_Bash('cat {} {} > {}'.format(PHASED_RESULT + '.vcf.noshop.header2', PHASED_RESULT + '.vcf.body', PHASED_RESULT + '.tobeDoubled.vcf'))

        RUN_Bash('Rscript src/Doubling_vcf.R {} {}'.format(PHASED_RESULT + '.tobeDoubled.vcf', PHASED_RESULT + '.Doubled.pre.vcf'))

        if not os.path.exists(PHASED_RESULT + '.tobeDoubled.vcf'):
            print(std_ERROR_MAIN_PROCESS_NAME + "Doubled phased file('{}') can't be found(or wasn't generated at all.".format(PHASED_RESULT + '.tobeDoubled.pre.vcf'))
            sys.exit()

        RUN_Bash('sed "s%CHROM%#CHROM%" {} > {}'.format(PHASED_RESULT + '.Doubled.pre.vcf', PHASED_RESULT + '.Doubled.pre2.vcf'))

        RUN_Bash('cat {} {} > {}'.format(PHASED_RESULT + '.vcf.header1', PHASED_RESULT + '.Doubled.pre2.vcf', MHC + '.QC.phasing_out_Doubled.vcf'))

        if not self.__save_intermediates:
            # os.system(' '.join(['rm', PHASED_RESULT + '.vcf']))
            os.system(' '.join(['rm', PHASED_RESULT + '.vcf.header1']))
            os.system(' '.join(['rm', PHASED_RESULT + '.vcf.header2']))
            os.system(' '.join(['rm', PHASED_RESULT + '.vcf.noshop.header2']))
            os.system(' '.join(['rm', PHASED_RESULT + '.vcf.body']))
            os.system(' '.join(['rm', PHASED_RESULT + '.tobeDoubled.vcf']))
            os.system(' '.join(['rm', PHASED_RESULT + '.Doubled.pre.vcf']))
            os.system(' '.join(['rm', PHASED_RESULT + '.Doubled.pre2.vcf']))


        return MHC + '.QC.phasing_out_Doubled.vcf'


