#-*- coding: utf-8 -*-

import os, sys, re
from os.path import join
import multiprocessing as mp
from statistics import mean

from src.GC_tricked_bgl2ori_bgl import GCtricedBGL2OriginalBGL
from src.RUN_Bash import RUN_Bash
from src.measureAccuracy import measureAccuracy
from src.redefineBPv1BH import redefineBP
from src.HLA_MultipleRefs import HLA_MultipleRefs



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



class HLA_Imputation_MM(object):

    def __init__(self, idx_process, MHC, _reference, _out, _hg,
                 _LINKAGE2BEAGLE, _BEAGLE2LINKAGE, _BEAGLE2VCF, _VCF2BEAGLE, _PLINK, _BEAGLE4,
                 _answer=None, f_save_intermediates=False, _MultP=1):


        ### General
        self.idx_process = idx_process
        self.__save_intermediates = f_save_intermediates

        # Prefixes
        self.OUTPUT_dir = os.path.dirname(_out)
        self.OUTPUT_dir_ref = join(self.OUTPUT_dir, os.path.basename(_reference))
        self.IMP_Result_prefix = _out # when using only multiple markers.

        # Result
        # self.IMP_Result = None  # Finally consensed Imputation output ('*.imputed.alleles').
        self.dict_ExonN_Panel = {_exonN: None for _exonN in __EXON__}
        self.dict_IMP_Result = {_exonN: {_overlap: None} for _exonN in __EXON__ for _overlap in __overlap__}
        self.accuracy = {_exonN: {_overlap: None} for _exonN in __EXON__ for _overlap in __overlap__}

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




        ###### < Reference panel for Exon 2, 3, 4 > ######

        multiple_panels = HLA_MultipleRefs(_reference, self.OUTPUT_dir_ref, _hg, self.BEAGLE2LINKAGE, self.PLINK, _MultP=_MultP)
        self.dict_ExonN_Panel = multiple_panels.ExonN_Panel

        # [Temporary Hard-coding]
        # self.dict_ExonN_Panel['exon2'] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190709_MM/T1DGC_REF.exon2'
        # self.dict_ExonN_Panel['exon3'] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190709_MM/T1DGC_REF.exon3'
        # self.dict_ExonN_Panel['exon4'] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190709_MM/T1DGC_REF.exon4'



        ###### < Main - 'IMPUTE', 'CONVERT_OUT' > ######

        if _MultP == 1:

            ## Main iteration over 'overlap'
            for _exonN in __EXON__:

                ### (1) CONVERT_IN

                [DOUBLED_PHASED_RESULT, REF_PHASED_VCF] = self.CONVERT_IN(MHC, self.dict_ExonN_Panel[_exonN], _out, _hg, _exonN)


                for _overlap in __overlap__:
                    self.dict_IMP_Result[_exonN][_overlap] = self.IMPUTATION_MM(DOUBLED_PHASED_RESULT, REF_PHASED_VCF,
                                                                                _exonN, _overlap, MHC, self.dict_ExonN_Panel[_exonN], _out, _hg)


        else:

            ###### < 'CONVERT_IN' in parallel. > ######

            """
            - In Multiple Marker technique, Pre-phasing in 'CONVERT_IN' block takes tremendous time.
            - working process executed by multiprocessing can't call multiprocessing further.
            - 'CONVERT_IN' is independent to 'IMPUTE' block.

            """

            ### (1) CONVERT_IN (by exon 2,3,4)

            pool_prephasing = mp.Pool(processes=_MultP)
            dict_Pool_prephasing = {_exonN: pool_prephasing.apply_async(self.CONVERT_IN, (MHC, self.dict_ExonN_Panel[_exonN], _out, _hg, _exonN))
                                    for _exonN in __EXON__}

            pool_prephasing.close()
            pool_prephasing.join()


            for _exonN in __EXON__:

                [a, b] = dict_Pool_prephasing[_exonN].get()
                self.dict_DOUBLED_PHASED_RESULT[_exonN] = a
                self.dict_REF_PHASED_VCF[_exonN] = b

            # print("DOUBLED_PHASED_RESULT :\n{}".format(self.dict_DOUBLED_PHASED_RESULT))
            # print("REF_PHASED_VCF :\n{}".format(self.dict_REF_PHASED_VCF))

            # [Temporary Hardcoding - 'CONVERT_IN' in parallel.
            # self.dict_DOUBLED_PHASED_RESULT['exon2'] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190709_MM/set_aside/_3_HM_CEU_T1DGC_REF.MHC.exon2.QC.phasing_out_doubled.vcf'
            # self.dict_REF_PHASED_VCF['exon2'] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190709_MM/T1DGC_REF.exon2.phased.vcf'
            #
            # self.dict_DOUBLED_PHASED_RESULT['exon3'] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190709_MM/set_aside/_3_HM_CEU_T1DGC_REF.MHC.exon3.QC.phasing_out_doubled.vcf'
            # self.dict_REF_PHASED_VCF['exon3'] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190709_MM/T1DGC_REF.exon3.phased.vcf'
            #
            # self.dict_DOUBLED_PHASED_RESULT['exon4'] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190709_MM/set_aside/_3_HM_CEU_T1DGC_REF.MHC.exon4.QC.phasing_out_doubled.vcf'
            # self.dict_REF_PHASED_VCF['exon4'] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190709_MM/T1DGC_REF.exon4.phased.vcf'


            ### Main iteration

            pool = mp.Pool(processes=_MultP if _MultP <= 9 else 9)

            dict_Pool = {_exonN: {_overlap: pool.apply_async(self.IMPUTATION_MM, (self.dict_DOUBLED_PHASED_RESULT[_exonN], self.dict_REF_PHASED_VCF[_exonN], _exonN, _overlap, MHC, self.dict_ExonN_Panel[_exonN], _out, _hg))
                                  for _overlap in __overlap__}
                         for _exonN in __EXON__}

            pool.close()
            pool.join()


            for _exonN in __EXON__:
                for _overlap in __overlap__:
                    self.dict_IMP_Result[_exonN][_overlap] = dict_Pool[_exonN][_overlap].get()


        # [Temporary Hardcoding - measuring accuracy]
        # self.dict_IMP_Result['exon2'][3000] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190708_MM/_3_HM_CEU_T1DGC_REF.exon2.3000.imputed.alleles'
        # self.dict_IMP_Result['exon2'][4000] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190708_MM/_3_HM_CEU_T1DGC_REF.exon2.4000.imputed.alleles'
        # self.dict_IMP_Result['exon2'][5000] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190708_MM/_3_HM_CEU_T1DGC_REF.exon2.5000.imputed.alleles'
        # self.dict_IMP_Result['exon3'][3000] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190708_MM/_3_HM_CEU_T1DGC_REF.exon2.3000.imputed.alleles'
        # self.dict_IMP_Result['exon3'][4000] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190708_MM/_3_HM_CEU_T1DGC_REF.exon2.4000.imputed.alleles'
        # self.dict_IMP_Result['exon3'][5000] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190708_MM/_3_HM_CEU_T1DGC_REF.exon2.5000.imputed.alleles'
        # self.dict_IMP_Result['exon4'][3000] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190708_MM/_3_HM_CEU_T1DGC_REF.exon2.3000.imputed.alleles'
        # self.dict_IMP_Result['exon4'][4000] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190708_MM/_3_HM_CEU_T1DGC_REF.exon2.4000.imputed.alleles'
        # self.dict_IMP_Result['exon4'][5000] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190708_MM/_3_HM_CEU_T1DGC_REF.exon2.5000.imputed.alleles'


        # Getting Accuracy
        if _answer:

            if os.path.isfile(_answer):

                for _exonN in __EXON__:
                    for _overlap in __overlap__:
                        print("Imputation Result({}, {}) : {}".format(_exonN, _overlap, self.dict_IMP_Result[_exonN][_overlap]))
                        self.accuracy[_exonN][_overlap] = measureAccuracy(_answer, self.dict_IMP_Result[_exonN][_overlap], 'all', _out+'.{}.{}.imputed.alleles.accuracy'.format(_exonN, _overlap))

                self.getPosteriorAccuracy(self.accuracy, _out+'.imputed.alleles.accuracy.avg')

            else:
                print(std_WARNING_MAIN_PROCESS_NAME + "Answer file to get accuracy('{}') can't be found. Skipping accuracy calculation.\n"
                                                      "Please check '--answer/-an' argument again.".format(_answer))


        else:
            print(std_MAIN_PROCESS_NAME + "No answer file to calculate accuracy.")




        ### General Removal
        if not self.__save_intermediates:
            RUN_Bash('rm {}'.format(MHC + '.QC.nopheno.ped'))
            RUN_Bash('rm {}'.format(MHC + '.QC.dat'))
            # RUN_Bash('rm {}'.format(join(self.OUTPUT_dir, 'selected_snp.txt')))




    def CONVERT_IN(self, MHC, _reference, _out, _hg, _exonN):

        __MHC_exonN__ = MHC+'.{}'.format(_exonN)
        __out_exonN__ = _out+'.{}'.format(_exonN)
        __reference_exonN__ = os.path.basename(_reference)
        OUTPUT_dir_ref_exonN = join(self.OUTPUT_dir, __reference_exonN__)


        print("[{}] Converting data to beagle format.".format(self.idx_process))
        self.idx_process += 1


        RUN_Bash(self.LINKAGE2BEAGLE + ' pedigree={} data={} beagle={} standard=true > {}'.format(
            MHC + '.QC.nopheno.ped', MHC + '.QC.dat', __MHC_exonN__ + '.QC.bgl', __out_exonN__+'.bgl.log'))

        if not self.__save_intermediates:
            # os.system('rm {}'.format(MHC + '.QC.nopheno.ped'))
            # os.system('rm {}'.format(MHC + '.QC.dat'))
            os.system('rm {}'.format(__out_exonN__+'.bgl.log'))




        ### Converting data to reference_markers_Position (Dispersing same genomic position of some markers.)

        RefinedMarkers = redefineBP(_reference + '.markers', OUTPUT_dir_ref_exonN+'.refined.markers')
        # self.refined_REF_markers = RefinedMarkers # => This will be used in 'CONVERT_OUT'.




        ### Converting data to target_markers_Position and extract not_including snp.

        RUN_Bash('awk \'{print $2" "$4" "$5" "$6}\' %s > %s' % (MHC + '.QC.bim', __MHC_exonN__ + '.QC.markers'))

        RUN_Bash('Rscript src/excluding_snp_and_refine_target_position-v1COOK02222017.R {} {} {}'.format(
            __MHC_exonN__+'.QC.markers', RefinedMarkers, __MHC_exonN__+'.QC.pre.markers'
        ))
        if not self.__save_intermediates:
            os.system(' '.join(['rm', __MHC_exonN__ + '.QC.markers']))

        RUN_Bash('mv {} {}'.format(__MHC_exonN__+'.QC.bgl', __MHC_exonN__+'.QC.pre.bgl.phased'))

        RUN_Bash("awk '{print $1}' %s > %s" % (__MHC_exonN__+'.QC.pre.markers', join(self.OUTPUT_dir, 'selected_snp.{}.txt'.format(_exonN))))


        from src.Panel_subset import Panel_Subset
        qc_refined = Panel_Subset(__MHC_exonN__ + '.QC.pre', 'all', join(self.OUTPUT_dir, 'selected_snp.{}.txt'.format(_exonN)),
                                  __MHC_exonN__ + '.QC.refined')

        if not self.__save_intermediates:
            RUN_Bash('rm {}'.format(__MHC_exonN__ + '.QC.pre.bgl.phased'))
            RUN_Bash('rm {}'.format(__MHC_exonN__ + '.QC.pre.markers'))
            RUN_Bash('rm {}'.format(join(self.OUTPUT_dir, 'selected_snp.{}.txt'.format(_exonN))))




        ### Converting data to GC_change_beagle format.

        from src.bgl2GC_trick_bgl import Bgl2GC

        # target
        [GCchangeBGL, GCchangeMarkers] = Bgl2GC(__MHC_exonN__ + '.QC.refined.bgl.phased', __MHC_exonN__ + '.QC.refined.markers',
                                                __MHC_exonN__ + '.QC.GCchange.bgl', __MHC_exonN__ + '.QC.GCchange.markers')

        self.GCchangeBGL = GCchangeBGL # it will be used in 'CONVERT_OUT' with Genetic Map

        # print("<Target GCchanged bgl and marker file>\n"
        #       "bgl : {}\n"
        #       "markers : {}".format(GCchangeBGL, GCchangeMarkers))

        # reference
        [GCchangeBGL_REF, GCchangeMarkers_REF] = Bgl2GC(_reference + '.bgl.phased', RefinedMarkers,
                                                        OUTPUT_dir_ref_exonN + '.GCchange.bgl.phased',
                                                        OUTPUT_dir_ref_exonN + '.GCchange.markers')
        # print("<Reference GCchanged bgl and marker file>\n"
        #       "bgl : {}\n"
        #       "markers : {}".format(GCchangeBGL_REF, GCchangeMarkers_REF))

        if not self.__save_intermediates:

            RUN_Bash('rm {}'.format(__MHC_exonN__ + '.QC.refined.bgl.phased'))
            RUN_Bash('rm {}'.format(__MHC_exonN__ + '.QC.refined.markers'))
            # RUN_Bash('rm {}'.format(RefinedMarkers))

            # os.system(' '.join(['rm', RefinedMarkers])) # => This will be used in 'CONVERT_OUT" when not using Multiple Markers.




        ### Converting data to vcf_format

        # target
        RUN_Bash(self.BEAGLE2VCF + ' 6 {} {} 0 > {}'.format(GCchangeMarkers, GCchangeBGL, __MHC_exonN__+'.QC.vcf'))

        MHC_QC_VCF_exonN = __MHC_exonN__ + '.QC.vcf'


        # reference
        RUN_Bash(self.BEAGLE2VCF + ' 6 {} {} 0 > {}'.format(GCchangeMarkers_REF, GCchangeBGL_REF, OUTPUT_dir_ref_exonN + '.vcf'))

        reference_vcf = OUTPUT_dir_ref_exonN + '.vcf'




        ### Converting data to reference_phased

        RUN_Bash('sed "s%/%|%g" {} > {}'.format(reference_vcf, OUTPUT_dir_ref_exonN + '.phased.vcf'))

        REF_PHASED_VCF = OUTPUT_dir_ref_exonN + '.phased.vcf'

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
        (2) `REF_PHASED_VCF` := OUTPUT_dir_ref_exonN + '.phased.vcf'

        These two files are to be passed into Beagle phasing;
        """



        ############### < Multiple Markers > ###############


        ### Phasing & Doubling (only on Target Sample.)

        # Phasing
        PHASED_RESULT = self.Phasing(__MHC_exonN__, MHC_QC_VCF_exonN, REF_PHASED_VCF)

        # [Temporary Hardcoding for Phased Result]
        # PHASED_RESULT = "/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190708_MM/_3_HM_CEU_T1DGC_REF.MHC.exon2.QC.phasing_out_not_double"
        # print("[Temporary Hardcoding]Phased Result:\n{}".format(PHASED_RESULT))


        # Doubling
        DOUBLED_PHASED_RESULT = self.Doubling(__MHC_exonN__, PHASED_RESULT)



        return [DOUBLED_PHASED_RESULT, REF_PHASED_VCF]



    def IMPUTE(self, MHC, _out, _DOUBLED_PHASED_RESULT, _REF_PHASED_VCF, _overlap, _exonN):


        print("[{}] Performing HLA imputation (see {}.MHC.QC.imputation_out.log for progress).".format(self.idx_process, _out))
        self.idx_process += 1


        OUT = MHC + '.{}.{}.QC.doubled.imputation_out'.format(_exonN, _overlap)



        """
        # multiple_nomap
        java -jar beagle4.jar gt=$MHC.QC.phasing_out_double.vcf ref=$REFERENCE.phased.vcf out=$MHC.QC.double.imputation_out impute=true lowmem=true
        
        # multiple_map
        java -jar beagle4.jar gt=$MHC.QC.phasing_out_double.vcf ref=$REFERENCE.phased.vcf out=$MHC.QC.double.imputation_out impute=true lowmem=true gprobs=true ne=10000 overlap=5000 err=$aver_erate map=$geneticMap.refined.map  
        """

        command = '{} gt={} ref={} out={} impute=true lowmem=true gprobs=true ne=10000 overlap={}'.format(self.BEAGLE4, _DOUBLED_PHASED_RESULT, _REF_PHASED_VCF, OUT, _overlap)
        # print(command)
        if not os.system(command):
            if not self.__save_intermediates:
                # os.system(' '.join(['rm', OUT + '.log'])) # Imputation Log file will be saved.
                # os.system(' '.join(['rm', _DOUBLED_PHASED_RESULT]))
                # os.system(' '.join(['rm', _REF_PHASED_VCF]))
                pass # for temporarily
        else:
            print(std_ERROR_MAIN_PROCESS_NAME + "Imputation failed.")
            sys.exit()


        RUN_Bash('gzip -d -f {}.vcf.gz'.format(OUT))

        return OUT + '.vcf'



    def CONVERT_OUT(self, _raw_IMP_Result, _reference, _out, _overlap, _exonN):


        Prefix_raw_IMP_Result = _raw_IMP_Result.rstrip('.vcf')
        OUTPUT_dir_ref = join(self.OUTPUT_dir, os.path.basename(_reference))
        OUT = _out+'.{}.{}'.format(_exonN, _overlap)


        ### vcf2HLAVCF
        for _hla in HLA_names:
            command = 'grep HLA_%s %s > %s' % (_hla, _raw_IMP_Result, Prefix_raw_IMP_Result+'.EXON_VCF_HLA_{}.txt'.format(_hla))
            # print(command)
            os.system(command)

        # Is file empty? (Introduced by W. Choi)
        vcf2HLAVCF = {_hla: None for _hla in HLA_names}

        for _hla in HLA_names:
            if os.path.getsize(Prefix_raw_IMP_Result+'.EXON_VCF_HLA_{}.txt'.format(_hla)) > 0:
                # print(Prefix_raw_IMP_Result+'.EXON_VCF_HLA_{}.txt'.format(_hla))
                vcf2HLAVCF[_hla] = Prefix_raw_IMP_Result+'.EXON_VCF_HLA_{}.txt'.format(_hla)


        ### DP_min_selection.R
        for _hla in HLA_names:
            if vcf2HLAVCF[_hla]:
                command = 'Rscript src/DP_min_selection.R {} {}'.format(vcf2HLAVCF[_hla],
                                                                        Prefix_raw_IMP_Result+'.EXON_VCF_DP_MIN_VCF_HLA_{}.txt'.format(_hla))
                # print(command)
                if not os.system(command):
                    vcf2HLAVCF[_hla] = Prefix_raw_IMP_Result+'.EXON_VCF_DP_MIN_VCF_HLA_{}.txt'.format(_hla)
                    # Replacing('*.EXON_VCF_HLA{}.txt' -> '*.EXON_VCF_DP_MIN_VCF_HLA_{}.txt')

            else:
                print('No HLA_{} alleles. DP_min_seleciton for HLA_{} will be skipped.'.format(_hla, _hla))

            # Removing '*.EXON_VCF_HLA_{}.txt'
            os.system('rm {}'.format(Prefix_raw_IMP_Result + '.EXON_VCF_HLA_{}.txt'.format(_hla)))



        # keep vcf header
        command = 'grep "#" {} > {}'.format(_raw_IMP_Result, _raw_IMP_Result+'.header.txt')
        # print(command)
        os.system(command)

        # Concatenate body part in the genomic position order.
        to_cat = []
        for _hla in HLA_names_gen:
            if vcf2HLAVCF[_hla]:
                to_cat.append(vcf2HLAVCF[_hla])

        to_cat = ' '.join(to_cat)

        command = 'cat {} > {}'.format(to_cat, Prefix_raw_IMP_Result+'.DP_MIN_VCF_HLA_all.txt')
        # print(command)
        if not os.system(command):

            if not self.__save_intermediates:
                for _hla in HLA_names:
                    if vcf2HLAVCF[_hla]:
                        os.system('rm {}'.format(vcf2HLAVCF[_hla]))


        # complete_vcf_HLA
        command = 'cat {} {} > {}'.format(_raw_IMP_Result+'.header.txt', Prefix_raw_IMP_Result+'.DP_MIN_VCF_HLA_all.txt',
                                          Prefix_raw_IMP_Result+'.DP_MIN_VCF_HLA_all_with_header.txt')
        # print(command)
        if not os.system(command):
            # os.system('rm {}'.format(_raw_IMP_Result))
            os.system('rm {}'.format(_raw_IMP_Result+'.header.txt'))
            os.system('rm {}'.format(Prefix_raw_IMP_Result+'.DP_MIN_VCF_HLA_all.txt'))


        # Get Marker
        command = 'grep HLA {} > {}'.format(_reference+'.markers', OUTPUT_dir_ref+'.HLA.markers')
        # print(command)
        os.system(command)


        #####


        command = 'cat {} | {} 0 {}'.format(Prefix_raw_IMP_Result+'.DP_MIN_VCF_HLA_all_with_header.txt',
                                            self.VCF2BEAGLE, Prefix_raw_IMP_Result+'.DP_MIN_Beagle_HLA_all_with_header')
        # print(command)
        if not os.system(command):
            if not self.__save_intermediates:
                os.system('rm {}'.format(Prefix_raw_IMP_Result+'.DP_MIN_VCF_HLA_all_with_header.txt'))
                os.system('rm {}'.format(Prefix_raw_IMP_Result+'.DP_MIN_Beagle_HLA_all_with_header.int'))
                os.system('rm {}'.format(Prefix_raw_IMP_Result+'.DP_MIN_Beagle_HLA_all_with_header.markers'))


        command = 'gzip -d -f {}'.format(Prefix_raw_IMP_Result+'.DP_MIN_Beagle_HLA_all_with_header.bgl.gz')
        # print(command)
        os.system(command)


        command = 'head -1 {} > {}'.format(Prefix_raw_IMP_Result+'.DP_MIN_Beagle_HLA_all_with_header.bgl',
                                           Prefix_raw_IMP_Result+'.DP_MIN_Beagle_HLA_all_with_header.bgl.header.txt')
        # print(command)
        os.system(command)


        command = 'cat {} {} > {}'.format(Prefix_raw_IMP_Result+'.DP_MIN_Beagle_HLA_all_with_header.bgl.header.txt',
                                          Prefix_raw_IMP_Result+'.DP_MIN_Beagle_HLA_all_with_header.bgl',
                                          Prefix_raw_IMP_Result+'.DP_MIN_Beagle_HLA_all_with_header_with_fid.bgl') # (2019. 05. 28.) FID라고 적혀있는데 I행만 두 개 더 들어가는게 꺼림찍함.
        # print(command)
        if not os.system(command):
            if not self.__save_intermediates:
                os.system('rm {}'.format(Prefix_raw_IMP_Result+'.DP_MIN_Beagle_HLA_all_with_header.bgl'))
                os.system('rm {}'.format(Prefix_raw_IMP_Result+'.DP_MIN_Beagle_HLA_all_with_header.bgl.header.txt'))



        # gc_change_ori_bgl.

        GC_decodedBGL = GCtricedBGL2OriginalBGL(Prefix_raw_IMP_Result+'.DP_MIN_Beagle_HLA_all_with_header_with_fid.bgl',
                                                OUTPUT_dir_ref+'.HLA.markers', Prefix_raw_IMP_Result+'.DP_MIN_Beagle_HLA_all_with_header.notGC.bgl')
        # print(GC_decodedBGL)

        if not self.__save_intermediates:
            os.system('rm {}'.format(Prefix_raw_IMP_Result+'.DP_MIN_Beagle_HLA_all_with_header_with_fid.bgl'))
            # os.system('rm {}'.format(OUTPUT_dir_ref+'.HLA.markers')) # This file must be shared between overlap 3000, 4000, 5000 of exon_N


        from src.BGL2Alleles_for_merge import BGL2Alleles4Merge

        DOUBLE_ALLELES = BGL2Alleles4Merge(GC_decodedBGL, Prefix_raw_IMP_Result+'.DP_MIN_Beagle_HLA_all_double.alleles', 'all')
        # print(DOUBLE_ALLELES)

        if not self.__save_intermediates:
            os.system('rm {}'.format(GC_decodedBGL))




        # Double2Single

        command = 'Rscript src/Double_alleles_decoder.R {} {}'.format(DOUBLE_ALLELES, OUT+'.imputed.alleles') # single
        # print(command)
        if not os.system(command):
            if not self.__save_intermediates:
                os.system('rm {}'.format(DOUBLE_ALLELES))


        __RETURN__ = OUT+'.imputed.alleles'


        self.idx_process += 1
        return __RETURN__



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

        RUN_Bash('grep ^## {} > {}'.format(PHASED_RESULT + '.vcf', PHASED_RESULT + '.vcf.header1')) # Header part with '##'
        RUN_Bash('grep -v ^## {} | head -n 1 > {}'.format(PHASED_RESULT + '.vcf', PHASED_RESULT + '.vcf.header2')) # Header part with '#'
        RUN_Bash('grep -v ^# {} > {}'.format(PHASED_RESULT + '.vcf', PHASED_RESULT + '.vcf.body')) # Body part without '#' or '##'

        RUN_Bash('sed "s%#%%" {} > {}'.format(PHASED_RESULT + '.vcf.header2', PHASED_RESULT + '.vcf.noshop.header2'))
        RUN_Bash('cat {} {} > {}'.format(PHASED_RESULT + '.vcf.noshop.header2', PHASED_RESULT + '.vcf.body', PHASED_RESULT + '.tobeDoubled.vcf'))

        RUN_Bash('Rscript src/Doubling_vcf.R {} {}'.format(PHASED_RESULT + '.tobeDoubled.vcf', PHASED_RESULT + '.Doubled.pre.vcf'))

        if not os.path.exists(PHASED_RESULT + '.tobeDoubled.vcf'):
            print(std_ERROR_MAIN_PROCESS_NAME + "Doubled phased file('{}') can't be found(or wasn't generated at all.".format(PHASED_RESULT + '.tobeDoubled.pre.vcf'))
            sys.exit()

        RUN_Bash('sed "s%CHROM%#CHROM%" {} > {}'.format(PHASED_RESULT + '.Doubled.pre.vcf', PHASED_RESULT + '.Doubled.pre2.vcf'))


        RUN_Bash('cat {} {} > {}'.format(PHASED_RESULT + '.vcf.header1', PHASED_RESULT + '.Doubled.pre2.vcf', MHC+'.QC.phasing_out_Doubled.vcf'))

        if not self.__save_intermediates:
            # os.system(' '.join(['rm', PHASED_RESULT + '.vcf']))
            os.system(' '.join(['rm', PHASED_RESULT + '.vcf.header1']))
            os.system(' '.join(['rm', PHASED_RESULT + '.vcf.header2']))
            os.system(' '.join(['rm', PHASED_RESULT + '.vcf.noshop.header2']))
            os.system(' '.join(['rm', PHASED_RESULT + '.vcf.body']))
            os.system(' '.join(['rm', PHASED_RESULT + '.tobeDoubled.vcf']))
            os.system(' '.join(['rm', PHASED_RESULT + '.Doubled.pre.vcf']))
            os.system(' '.join(['rm', PHASED_RESULT + '.Doubled.pre2.vcf']))



        return MHC+'.QC.phasing_out_Doubled.vcf'



    def IMPUTATION_MM(self, _DOUBLED_PHASED_RESULT, _REF_PHASED_VCF, _exonN, _overlap, MHC, _reference, _out, _hg):


        ### (2) IMPUTE

        raw_IMP_Reuslt = self.IMPUTE(MHC, _out, _DOUBLED_PHASED_RESULT, _REF_PHASED_VCF, _overlap, _exonN)
        # print("raw Imputed Reuslt : \n{}".format(raw_IMP_Reuslt)) # *.vcf

        # [Temporary Hard-coding]
        # raw_IMP_Reuslt = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190708_MM/_3_HM_CEU_T1DGC_REF.MHC.exon2.3000.QC.doubled.imputation_out.vcf'


        ### (3) CONVERT_OUT

        self.IMP_Result = self.CONVERT_OUT(raw_IMP_Reuslt, _reference, _out, _overlap, _exonN)

        return self.IMP_Result



    def getPosteriorAccuracy(self, _accuracy, _out=None, __also2D=False):

        # print(_accuracy)

        dict_accuracy_4D = {_hla: [] for _hla in HLA_names}
        dict_accuracy_2D = {_hla: [] for _hla in HLA_names}

        for _hla in HLA_names:
            for _exonN in __EXON__:

                if _exonN == 'exon4' and not isClassI[_hla]:
                    continue
                else:
                    for _overlap in __overlap__:
                        dict_accuracy_4D[_hla].append(_accuracy[_exonN][_overlap]['4D'][_hla])

        # print(dict_accuracy_4D)
        dict_mean_accuracy_4D = {_hla: mean(dict_accuracy_4D[_hla]) for _hla in HLA_names}
        print("\nAverage of imputation accuracy :\n{}".format(dict_mean_accuracy_4D))


        if __also2D:

            for _hla in HLA_names:
                for _exonN in __EXON__:

                    if _exonN == 'exon4' and not isClassI[_hla]:
                        continue
                    else:
                        for _overlap in __overlap__:
                            dict_accuracy_2D[_hla].append(_accuracy[_exonN][_overlap]['2D'][_hla])

            # print(dict_accuracy_2D)
            dict_mean_accuracy_2D = {_hla: mean(dict_accuracy_2D[_hla]) for _hla in HLA_names}
            print(dict_mean_accuracy_2D)


        ### file write
        if _out:
            with open(_out, 'w') as f_out:
                for _hla in HLA_names:
                    f_out.write("%s\t4D\t%.5f\n"%(_hla, float(dict_mean_accuracy_4D[_hla])))
                    if __also2D:
                        f_out.write("%s\t2D\t%.5f\n" % (_hla, float(dict_mean_accuracy_2D[_hla])))

            return _out
        else:
            if not __also2D:
                return dict_mean_accuracy_4D
            else:
                return [dict_mean_accuracy_4D, dict_mean_accuracy_2D]
