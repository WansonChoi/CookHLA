#-*- coding: utf-8 -*-

import os, sys, re
from os.path import join
import pandas as pd

from src.GC_tricked_bgl2ori_bgl import GCtricedBGL2OriginalBGL


########## < Core Varialbes > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
HLA_names_gen = ["A", "C", "B", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]

# __overlap__ = [3000, 4000, 5000]
__overlap__ = [3000]


class HLA_Imputation_MM(object):

    def __init__(self, idx_process, MHC, _reference, _out, _hg,
                 _LINKAGE2BEAGLE, _BEAGLE2LINKAGE, _BEAGLE2VCF, _VCF2BEAGLE,
                 _PLINK, _BEAGLE4,
                 __save_intermediates=False,
                 _aver_erate=None, _Genetic_Map=None, f_useGeneticMap=False,
                 f_useMultipleMarkers=False, _exonN_ = None, _answer=None):


        ### Class variables
        self.idx_process = idx_process
        self.__save_intermediates = __save_intermediates

        self.f_useGeneticMap = f_useGeneticMap
        self.f_useMultipleMarkers = f_useMultipleMarkers

        self.exonN = _exonN_

        self.refined_Genetic_Map = None
        self.GCchangeBGL = None

        # 'CONVERT_OUT'
        self.refined_REF_markers = None

        # Result
        self.OUTPUT_dir = os.path.dirname(_out)
        self.OUTPUT_dir_ref = join(self.OUTPUT_dir, os.path.basename(_reference))
        self.OUTPUT_dir_GM = join(self.OUTPUT_dir, os.path.basename(_Genetic_Map)) if f_useGeneticMap else None

        self.raw_IMP_Reuslt = None
        self.IMP_Result_prefix = _out # when using only multiple markers.
        self.IMP_Result = None  # '*.imputed.alleles'

        self.accuracy = None







        ###### < Main - 'CONVERT_IN', 'IMPUTE', 'CONVERT_OUT' > ######


        """
        앞서 언급한 Multiple Marker와 Adaptive Genetic Map을 모두 활용하는 경우, 사실상 exonN에 대한 Genetic Map을 새로 떠야하는
        상황이기 때문에 CONVERT_IN() 함수 들어가기 전에 이쯤에서 Genetic Map을 만들어야 할 듯.
        """


        ### (1) CONVERT_IN

        [MHC_QC_VCF, REF_PHASED_VCF] = self.CONVERT_IN(MHC, _reference, _out, _hg,
                                                       _LINKAGE2BEAGLE, _BEAGLE2VCF, _PLINK, _BEAGLE4,
                                                       _aver_erate=_aver_erate, _Genetic_Map=_Genetic_Map)

        # self.CONVERT_IN(MHC, _reference, _out, _hg, _LINKAGE2BEAGLE, _BEAGLE2VCF, _PLINK, _BEAGLE4,
        #                 _aver_erate=_aver_erate, _Genetic_Map=_Genetic_Map)

        # print("Convert_IN :\n{}\n{}".format(MHC_QC_VCF, REF_PHASED_VCF))




        if self.f_useGeneticMap:


            self.IMP_Result = {_overlap_: None for _overlap_ in __overlap__}    # Imputation Results.

            if self.f_useMultipleMarkers:

                ### Imputation with both 'Adaptive Genetic Map' and 'Multiple Markers' (Main CookHLA).
                # No multiprocessing (Sequentially)

                for _overlap_ in __overlap__:

                    t = self.IMPUTE_GM(_overlap_, _out, MHC_QC_VCF, REF_PHASED_VCF, MHC, _reference, _aver_erate, _Genetic_Map,
                                       _BEAGLE4, _VCF2BEAGLE, _BEAGLE2LINKAGE, _PLINK)

                    self.IMP_Result[_overlap_] = t


            else:
                ### Imputation with only 'Adaptive Genetic Map'.
                # Multiprocessing

                import multiprocessing as mp

                pool = mp.Pool(processes=3)

                dict_Pool = {_overlap_: pool.apply_async(self.IMPUTE_GM, (_overlap_, _out, MHC_QC_VCF, REF_PHASED_VCF, MHC, _reference, _aver_erate, _Genetic_Map,
                                                                          _BEAGLE4, _VCF2BEAGLE, _BEAGLE2LINKAGE, _PLINK)) for _overlap_ in __overlap__}

                self.IMP_Result = {_overlap_: _OUT.get() for _overlap_, _OUT in dict_Pool.items()}



            print(self.IMP_Result)


        else:

            # Plain Single Implementation

            if f_useMultipleMarkers:
                print(std_MAIN_PROCESS_NAME + "exonN: {}".format(_exonN_))


            ### (2) IMPUTE

            # Temporary Hard coding
            # MHC_QC_VCF = 'tests/_3_CookHLA/20190523/temp/CONVERT_IN/_3_HM_CEU_T1DGC_REF.MHC.exon2.QC.phasing_out_not_double.doubled.vcf'
            # REF_PHASED_VCF = 'tests/_3_CookHLA/20190523/temp/CONVERT_IN/T1DGC_REF.exon2.phased.vcf'

            self.raw_IMP_Reuslt = self.IMPUTE(_out, MHC_QC_VCF, REF_PHASED_VCF, _BEAGLE4)
            print("raw Imputed Reuslt : {}".format(self.raw_IMP_Reuslt))


            ### (3) CONVERT_OUT

            # Temporary Hard-coding
            # print("Temporary Hard-coding for testing 'CONVERT_OUT'.")
            # self.raw_IMP_Reuslt = 'tests/_3_CookHLA/20190523/_3_HM_CEU_T1DGC_REF.exon2.QC.doubled.imputation_out.vcf'

            self.IMP_Result = self.CONVERT_OUT(MHC, _reference, _out, _VCF2BEAGLE, _BEAGLE2LINKAGE, _PLINK)

            print("\n\nImputation Result : {}".format(self.IMP_Result))



        ###### < Get Accuracy > ######

        if _answer:

            if not os.path.exists(_answer):
                print(std_WARNING_MAIN_PROCESS_NAME + "Given Answer file doesn't exist. Skipping calculating accuracy."
                                                      "Please check '--answer' argument again.")


            from src.measureAccuracy import measureAccuracy

            if isinstance(self.IMP_Result, str):
                self.accuracy = measureAccuracy(_answer, self.IMP_Result, 'all', outfile=_out+'.accuracy')
            elif isinstance(self.IMP_Result, dict):
                self.accuracy = {k: measureAccuracy(_answer, v, 'all') for k, v in self.IMP_Result.items()}




        ###### < Get Accuracy > ######

        # Removing
        # if not self.__save_intermediates:
        #     if self.f_useGeneticMap:
        #         os.system('rm {}'.format(self.refined_Genetic_Map))
        #         os.system('rm {}'.format(self.GCchangeBGL))
        #         os.system('rm {}'.format(self.refined_REF_markers))







    def CONVERT_IN(self, MHC, _reference, _out, _hg, _LINKAGE2BEAGLE, _BEAGLE2VCF, _PLINK, _BEAGLE4,
                   _aver_erate=None, _Genetic_Map=None):


        __MHC__ = MHC if not self.f_useMultipleMarkers else MHC+'.{}'.format(self.exonN)

        print("[{}] Converting data to beagle format.".format(self.idx_process))

        command = ' '.join(
            [_LINKAGE2BEAGLE, 'pedigree={}'.format(MHC + '.QC.nopheno.ped'), 'data={}'.format(MHC + '.QC.dat'),
             'beagle={}'.format(__MHC__ + '.QC.bgl'), 'standard=true', '>', __MHC__ + '.QC.bgl.log'])  # Making '*.bgl' file.
        # print(command)
        os.system(command)

        if not self.__save_intermediates:
            # os.system(' '.join(['rm', MHC + '.QC.nopheno.ped']))
            # os.system(' '.join(['rm', MHC + '.QC.dat']))
            os.system('rm {}'.format(__MHC__ + '.QC.bgl.log'))


        ### Converting data to reference_markers_Position (Dispersing same genomic position of some markers.)

        from src.redefineBPv1BH import redefineBP

        RefinedMarkers = redefineBP(_reference + '.markers', self.OUTPUT_dir_ref+'.refined.markers')
        self.refined_REF_markers = RefinedMarkers # => This will be used in 'CONVERT_OUT'.


        ### Converting data to target_markers_Position and extract not_including snp.

        command = ' '.join(['awk \'{print $2" "$4" "$5" "$6}\'', MHC + '.QC.bim', '>',
                            __MHC__ + '.QC.markers'])  # Making '*.markers' file.
        # print(command)
        os.system(command)

        # if not self.__save_intermediates:
        #     os.system(' '.join(['rm', MHC + '.QC.{bed,bim,fam,log}']))

        command = ' '.join(['Rscript src/excluding_snp_and_refine_target_position-v1COOK02222017.R',
                            __MHC__ + '.QC.markers', RefinedMarkers, __MHC__ + '.QC.pre.markers'])
        # print(command)
        os.system(command)

        if not self.__save_intermediates:
            os.system(' '.join(['rm', __MHC__ + '.QC.markers']))

        command = ' '.join(['mv', __MHC__ + '.QC.bgl', __MHC__ + '.QC.pre.bgl.phased'])
        # print(command)
        os.system(command)

        command = ' '.join(
            ["awk '{print $1}'", __MHC__ + '.QC.pre.markers', '>', os.path.join(self.OUTPUT_dir, 'selected_snp.txt')])
        # print(command)
        os.system(command)

        from src.Panel_subset import Panel_Subset
        qc_refined = Panel_Subset(__MHC__ + '.QC.pre', 'all', join(self.OUTPUT_dir, 'selected_snp.txt'), __MHC__ + '.QC.refined')
        # print(qc_refined) # Refined Beagle files are generated here.

        if not self.__save_intermediates:
            os.system(' '.join(['rm', __MHC__ + '.QC.pre.{bgl.phased,markers}']))
            os.system(' '.join(['rm', join(self.OUTPUT_dir, 'selected_snp.txt')]))


        ### Converting data to GC_change_beagle format.

        from src.bgl2GC_trick_bgl import Bgl2GC

        # target
        [GCchangeBGL, GCchangeMarkers] = Bgl2GC(__MHC__ + '.QC.refined.bgl.phased', __MHC__ + '.QC.refined.markers',
                                                __MHC__ + '.QC.GCchange.bgl', __MHC__ + '.QC.GCchange.markers')

        self.GCchangeBGL = GCchangeBGL # it will be used in 'CONVERT_OUT' with Genetic Map

        # print("<Target GCchanged bgl and marker file>\n"
        #       "bgl : {}\n"
        #       "markers : {}".format(GCchangeBGL, GCchangeMarkers))

        # reference
        [GCchangeBGL_REF, GCchangeMarkers_REF] = Bgl2GC(_reference + '.bgl.phased', RefinedMarkers,
                                                        self.OUTPUT_dir_ref + '.GCchange.bgl.phased',
                                                        self.OUTPUT_dir_ref + '.GCchange.markers')
        # print("<Reference GCchanged bgl and marker file>\n"
        #       "bgl : {}\n"
        #       "markers : {}".format(GCchangeBGL_REF, GCchangeMarkers_REF))

        if not self.__save_intermediates:
            os.system(' '.join(['rm', __MHC__ + '.QC.refined.{bgl.phased,markers}']))

            if self.f_useMultipleMarkers:
                os.system(' '.join(['rm', RefinedMarkers])) # => This will be used in 'CONVERT_OUT" when not using Multiple Markers.


        ### Converting data to vcf_format

        # target
        command = ' '.join([_BEAGLE2VCF, '6', GCchangeMarkers, GCchangeBGL, '0', '>', __MHC__ + '.QC.vcf'])
        # print(command)
        os.system(command)

        MHC_QC_VCF = __MHC__ + '.QC.vcf'


        # reference
        command = ' '.join([_BEAGLE2VCF, '6', GCchangeMarkers_REF, GCchangeBGL_REF, '0', '>', self.OUTPUT_dir_ref + '.vcf'])
        # print(command)
        os.system(command)

        reference_vcf = self.OUTPUT_dir_ref + '.vcf'


        ### Converting data to reference_phased

        command = ' '.join(['sed "s%/%|%g"', reference_vcf, '>', self.OUTPUT_dir_ref + '.phased.vcf'])
        # print(command)
        os.system(command)

        REF_PHASED_VCF = self.OUTPUT_dir_ref + '.phased.vcf'

        if not self.__save_intermediates:
            os.system(' '.join(['rm', reference_vcf]))
            # if self.f_useMultipleMarkers:
            if not self.f_useGeneticMap:
                os.system(' '.join(['rm {}'.format(GCchangeBGL)])) # 'GCchangeBGL' will be used in 'CONVERT_OUT'
                os.system(' '.join(['rm {}'.format(GCchangeMarkers_REF)]))  # 'GCchangeMarkers_REF' will be used in 'CONVERT_OUT'
                os.system(' '.join(['rm {}'.format(GCchangeMarkers)]))
                os.system(' '.join(['rm {}'.format(GCchangeBGL_REF)]))

        """
        (1) MHC + '.QC.vcf',
        (2) self.OUTPUT_dir_ref + '.phased.vcf'

        These two files are to be passed into Beagle phasing;
        """

        ### Adaptive Genetic Map

        if self.f_useGeneticMap:

            """
            awk '{print $1" "$2" "$3}' $geneticMap > $geneticMap.first
            awk '{print $2}' $REFERENCE.GCchange.markers > $geneticMap.second
            paste -d " " $geneticMap.first $geneticMap.second > $geneticMap.refined.map

            rm $geneticMap.first
            rm $geneticMap.second

            """

            REFINED_GENTIC_MAP = self.OUTPUT_dir_GM + ('.{}.refined.map'.format(self.exonN) if self.f_useMultipleMarkers else '.refined.map')

            # if self.f_useMultipleMarkers:
            #     # When using both 'Adaptive Genetic Map' and 'Multiple Markers'.
            #     # 'Adaptive Genetic Map' and 'ExonN reference markers file' have different number of rows.
            #     # This block leaves markers which 'Adaptive Genetic Map' and 'ExonN Reference' both have.
            #
            #     from src.ManualInnerJoin import ManualInnerJoin
            #
            #     ManualInnerJoin(self.exonN, _Genetic_Map, GCchangeMarkers_REF, REFINED_GENTIC_MAP)
            #
            # else:

            command = 'awk \'{print $1" "$2" "$3}\' %s > %s' % (_Genetic_Map, self.OUTPUT_dir_GM+'.first')
            # print(command)
            os.system(command)

            command = 'awk \'{print $2}\' %s > %s' % (GCchangeMarkers_REF, self.OUTPUT_dir_GM+'.second')
            # print(command)
            os.system(command)

            command = 'paste -d " " {} {} > {}'.format(self.OUTPUT_dir_GM+'.first', self.OUTPUT_dir_GM+'.second', REFINED_GENTIC_MAP)   # 이렇게 column bind시키는데는 당연히 *.first, *.second 파일의 row수가 같을 거라고 가정하는 상황.
            # print(command)
            os.system(command)


            """
            그냥 Adaptive Genetic Map만 활용할때는, _reference에 진짜 reference가 들어오고 Genetic Map또한 진짜 reference를 기반으로
            만들어졌기 때문에 위 3줄의 Refined Genetic Map을 만드는데 문제가 없음.
            
            그러나, Multiple Marker를 Adaptive Genetic Map과 같이 활용하는 경우, Multiple Marker를 활용하는 경우 _reference에 
            HLA만 남긴 전처리된 다른 reference를 활용하기 때문에(HLA_MultipleRefs.py), 진짜 reference를 바탕으로 만들어진 Genetic Map과
            row수부터가 달라서 저 위 3줄로 Refined Genetic Map을 만들기가 어려워짐.
            
            보아하니 GM컬럼의 값 또한 ascending order로 주어져야하고, 그 와중에 BP 또한 ascending order로 주어져야 하니 GCchangeMarkers는
            redefineMap을 거치고 온 상태이기 때문에 같이 쓰기는 힘듬(한쪽이 ascending 시키면 다른게 ascending형태로 준비되어질 수 없음.)
            
            그래서 MakeEXON234_Panel.py부터 다시 손봐야 할 것 같음.
            
            
            참고로 GM에서 HLA~, rs~하는 애들만 남겨도 ("GM", "BP")모두가 ascending되게 할 수는 없기 때문에 소용없음. 
            => 생각해보면 Make_EXON234_Panel.py 건든다고 이 부분이 해결 안되는거 아닌가? rs~, HLA~이런 마커들만 남긴 reference panel에 대해
                Genetic Map을 새로 만들어야 할 거 같은데.
            """



            if os.path.exists(REFINED_GENTIC_MAP):

                self.refined_Genetic_Map = REFINED_GENTIC_MAP

                # if not self.__save_intermediates:
                #     os.system('rm {}'.format(self.OUTPUT_dir_GM+'.first'))
                #     os.system('rm {}'.format(self.OUTPUT_dir_GM+'.second'))
                #     os.system('rm {}'.format(GCchangeMarkers_REF)) # (Genetic Map) *.GCchange.markers is removed here.

            else:
                print(std_ERROR_MAIN_PROCESS_NAME + "Failed to generate Refined Genetic Map.")
                sys.exit()




        ### Mutliple Markers

        if not self.f_useMultipleMarkers:

            __RETURN__ = [MHC_QC_VCF, REF_PHASED_VCF]

        else:

            ### Phasing & Doubling (only on Target Sample.)

            # Phasing
            PHASED_RESULT = self.Phasing(__MHC__, MHC_QC_VCF, REF_PHASED_VCF, _BEAGLE4)

            # Doubling
            DOUBLED_PHASED_RESULT = self.Doubling(PHASED_RESULT)


            __RETURN__ = [DOUBLED_PHASED_RESULT, REF_PHASED_VCF]



        self.idx_process += 1

        return __RETURN__



    def IMPUTE(self, _out, _MHC_QC_VCF, _REF_PHASED_VCF, _BEAGLE4, _overlap=None, _aver_erate=None, _Genetic_Map=None):


        print("[{}] Performing HLA imputation (see {}.MHC.QC.imputation_out.log for progress).".format(self.idx_process, _out))


        OUT = _out + ('.QC.doubled.imputation_out' if self.f_useMultipleMarkers else '.QC.imputation_out')


        if self.f_useGeneticMap:

            ### Using both 'Multiple Markers' and 'Adaptive Genetic Map'.

            """
            java -jar beagle4.jar gt=$MHC.QC.phasing_out_double.vcf ref=$REFERENCE.phased.vcf out=$MHC.QC.double.imputation_out impute=true lowmem=true 
                                    gprobs=true ne=10000 overlap=5000 err=$aver_erate map=$geneticMap.refined.map
            """

            with open(_aver_erate, 'r') as f:
                aver_erate = f.readline().rstrip('\n')

            command = '{} gt={} ref={} out={} impute=true lowmem=true gprobs=true ne=10000 overlap={} err={} map={}'.format(_BEAGLE4, _MHC_QC_VCF, _REF_PHASED_VCF, OUT, _overlap, aver_erate, self.refined_Genetic_Map)
            print(command)
            if not os.system(command):
                if not self.__save_intermediates:
                    os.system('rm {}'.format(OUT+'.log'))
                    # os.system('rm {}'.format(_MHC_QC_VCF))
                    # os.system('rm {}'.format(_REF_PHASED_VCF)) # These 2 files
            else:
                print(std_ERROR_MAIN_PROCESS_NAME + "Imputation with Geneticmap failed.")
                sys.exit()


        else:

            """
            java -jar beagle4.jar gt=$MHC.QC.phasing_out_double.vcf ref=$REFERENCE.phased.vcf out=$MHC.QC.double.imputation_out impute=true lowmem=true 
            """

            command = '{} gt={} ref={} out={} impute=true lowmem=true'.format(_BEAGLE4, _MHC_QC_VCF, _REF_PHASED_VCF, OUT)
            # print(command)
            if not os.system(command):
                if not self.__save_intermediates:
                    # os.system(' '.join(['rm', OUT + '.log'])) # Imputation Log file will be saved.
                    os.system(' '.join(['rm', _MHC_QC_VCF]))
                    os.system(' '.join(['rm', _REF_PHASED_VCF]))
            else:
                print(std_ERROR_MAIN_PROCESS_NAME + "Imputation failed.")
                sys.exit()



        command = 'gzip -d -f {}.vcf.gz'.format(OUT)
        # print(command)
        os.system(command)

        __RETURN__ = OUT + '.vcf'
        self.idx_process += 1

        return __RETURN__



    def CONVERT_OUT(self, MHC, _reference, _out, _VCF2BEAGLE, _BEAGLE2LINKAGE, _PLINK, _overlap=None, raw_IMP_Result=None):


        if not self.f_useMultipleMarkers:

            __MHC__ = MHC+'.overlap{}'.format(_overlap)

            ### Converting imputation result in vcf file to beagle format.

            command = 'cat {} | {} 0 {}'.format(raw_IMP_Result, _VCF2BEAGLE, __MHC__ + '.QC.imputation_GCchange')
            # print(command)
            if not os.system(command):
                if not self.__save_intermediates:
                    os.system('rm {}'.format(__MHC__ + '.QC.imputation_GCchange.int'))

            command = 'gzip -d -f {}'.format(__MHC__ + '.QC.imputation_GCchange.bgl.gz')
            # print(command)
            os.system(command)


            ### Decoding GC-encoded values in beagle file to original values.

            # Decoding GC-encoding
            GC_decodedBGL = GCtricedBGL2OriginalBGL(__MHC__ + '.QC.imputation_GCchange.bgl', self.refined_REF_markers, __MHC__ + '.QC.imputation_ori.bgl')
            # print('GC_decodedBGL : {}'.format(GC_decodedBGL))

            if not self.__save_intermediates:
                os.system('rm {}'.format(__MHC__ + '.QC.imputation_GCchange.bgl'))
                os.system('rm {}'.format(__MHC__ + '.QC.imputation_GCchange.markers'))

            #
            # HLA_IMPUTED_Result = self.IMP_Result_prefix + '.bgl.phased'
            HLA_IMPUTED_Result = _out + '.bgl.phased'

            command = 'Rscript src/complete_header.R {} {} {}'.format(self.GCchangeBGL, GC_decodedBGL, HLA_IMPUTED_Result)
            # print(command)
            if not os.system(command):
                if not self.__save_intermediates:
                    # os.system('rm {}'.format(__MHC__ + '.QC.GCchange.bgl'))
                    os.system('rm {}'.format(GC_decodedBGL))


            ### Converting decoded beagle file to PLINK file.

            # command = 'cat {} | {} {}'.format(HLA_IMPUTED_Result, _BEAGLE2LINKAGE, _out+'.tmp')
            # print(command)
            # os.system(command)
            #
            # command = "cut -d ' ' -f6- {} > {}".format(_out+'.tmp.ped', _out+'.tmp')
            # print(command)
            # os.system(command)
            #
            # command = "paste -d ' ' {} {} | tr -d '\015' > {}".format(__MHC__ + '.fam', _out + '.tmp', HLA_IMPUTED_Result + '.ped') # *.ped
            # print(command)
            # os.system(command)
            #
            # command = 'cut -f1-4 {} > {}'.format(_reference+'.bim', HLA_IMPUTED_Result+'.map') # *.map
            # print(command)
            # os.system(command)
            #
            # command = 'cp {} {}'.format(__MHC__ + '.fam', HLA_IMPUTED_Result + '.fam')
            # print(command)
            # os.system(command)
            #
            # # Create PLINK bed format
            # command = '{} --ped {} --map {} --make-bed --out {}'.format(_PLINK, HLA_IMPUTED_Result+'.ped', HLA_IMPUTED_Result+'.map', HLA_IMPUTED_Result)
            # print(command)
            # os.system(command)
            #
            # command = 'rm {}'.format(HLA_IMPUTED_Result + '.fam')
            # print(command)
            # os.system(command)
            #
            # command = 'cp {} {}'.format(_reference+'.bim', HLA_IMPUTED_Result+'.bim')
            # print(command)
            # os.system(command)


            # BGL2Allele.py

            from src.BGL2Alleles import BGL2Alleles

            __RETURN__ = BGL2Alleles(HLA_IMPUTED_Result, _out+'.imputed.alleles', 'all')




        else:

            Prefix_raw_IMP_Result = self.raw_IMP_Reuslt.rstrip('.vcf')

            ### vcf2HLAVCF
            for _hla in HLA_names:
                command = 'grep HLA_%s %s > %s' % (_hla, self.raw_IMP_Reuslt, Prefix_raw_IMP_Result+'.EXON_VCF_HLA_{}.txt'.format(_hla))
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



            # Set aside vcf header
            command = 'grep "#" {} > {}'.format(self.raw_IMP_Reuslt, self.raw_IMP_Reuslt+'.header.txt')
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
            command = 'cat {} {} > {}'.format(self.raw_IMP_Reuslt+'.header.txt', Prefix_raw_IMP_Result+'.DP_MIN_VCF_HLA_all.txt',
                                              Prefix_raw_IMP_Result+'.DP_MIN_VCF_HLA_all_with_header.txt')
            # print(command)
            if not os.system(command):
                # os.system('rm {}'.format(self.raw_IMP_Reuslt))
                os.system('rm {}'.format(self.raw_IMP_Reuslt+'.header.txt'))
                os.system('rm {}'.format(Prefix_raw_IMP_Result+'.DP_MIN_VCF_HLA_all.txt'))


            # Get Marker
            command = 'grep HLA {} > {}'.format(_reference+'.markers', self.OUTPUT_dir_ref+'.HLA.markers')
            # print(command)
            os.system(command)


            #####


            command = 'cat {} | {} 0 {}'.format(Prefix_raw_IMP_Result+'.DP_MIN_VCF_HLA_all_with_header.txt',
                                                _VCF2BEAGLE, Prefix_raw_IMP_Result+'.DP_MIN_Beagle_HLA_all_with_header')
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
                                                    self.OUTPUT_dir_ref+'.HLA.markers', Prefix_raw_IMP_Result+'.DP_MIN_Beagle_HLA_all_with_header.notGC.bgl')
            # print(GC_decodedBGL)

            if not self.__save_intermediates:
                os.system('rm {}'.format(Prefix_raw_IMP_Result+'.DP_MIN_Beagle_HLA_all_with_header_with_fid.bgl'))
                os.system('rm {}'.format(self.OUTPUT_dir_ref+'.HLA.markers'))


            from src.BGL2Alleles_for_merge import BGL2Alleles4Merge

            DOUBLE_ALLELES = BGL2Alleles4Merge(GC_decodedBGL, Prefix_raw_IMP_Result+'.DP_MIN_Beagle_HLA_all_double.alleles', 'all')
            # print(DOUBLE_ALLELES)

            if not self.__save_intermediates:
                os.system('rm {}'.format(GC_decodedBGL))




            # Double2Single

            command = 'Rscript src/Double_alleles_decoder.R {} {}'.format(DOUBLE_ALLELES, self.IMP_Result_prefix+'.imputed.alleles') # single
            # print(command)
            if not os.system(command):
                if not self.__save_intermediates:
                    os.system('rm {}'.format(DOUBLE_ALLELES))


            __RETURN__ = self.IMP_Result_prefix+'.imputed.alleles'


        self.idx_process += 1
        return __RETURN__



    def Phasing(self, MHC, MHC_QC_VCF, REF_PHASED_VCF, _BEAGLE4):


        ### Phasing & Doubling (only on Target Sample.)

        command = ' '.join([_BEAGLE4, 'gt={} ref={} out={} impute=false > {}'.format(MHC_QC_VCF, REF_PHASED_VCF,
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




    def Doubling(self, PHASED_RESULT):

        ### Target data doubling step.

        command = 'gzip -d -f {}'.format(PHASED_RESULT + '.vcf.gz')
        # print(command)
        os.system(command)

        command = 'grep ^## {} > {}'.format(PHASED_RESULT + '.vcf', PHASED_RESULT + '.vcf.header')
        # print(command)
        os.system(command)

        command = 'grep -v ^## {} > {}'.format(PHASED_RESULT + '.vcf', PHASED_RESULT + '.vcf.body')
        # print(command)
        os.system(command)

        from src.Doubling_vcf import Doubling_vcf

        DOUBLED_VCF_body = Doubling_vcf(PHASED_RESULT + '.vcf.body', PHASED_RESULT + '.doubled.vcf.body')
        # print(DOUBLED_VCF_body)

        command = 'cat {} {} > {}'.format(PHASED_RESULT + '.vcf.header', DOUBLED_VCF_body,
                                          PHASED_RESULT + '.doubled.vcf')
        # print(command)
        os.system(command)

        if not self.__save_intermediates:
            os.system(' '.join(['rm', PHASED_RESULT + '.vcf']))
            os.system(' '.join(['rm', PHASED_RESULT + '.vcf.header']))
            os.system(' '.join(['rm', PHASED_RESULT + '.vcf.body']))
            os.system(' '.join(['rm', PHASED_RESULT + '.doubled.vcf.body']))



        return PHASED_RESULT + '.doubled.vcf'




    def getIDX_PROCESS(self):
        return self.idx_process



    def getImputationResult(self):
        return self.IMP_Result


    def IMPUTE_GM(self, _overlap_, _out, MHC_QC_VCF, REF_PHASED_VCF, MHC, _reference, _aver_erate, _Genetic_Map,
                  _BEAGLE4,_VCF2BEAGLE, _BEAGLE2LINKAGE, _PLINK):


        OUT_prefix = _out

        if self.f_useMultipleMarkers:
            print(std_MAIN_PROCESS_NAME + "exonN: {} / Overlap: {}".format(self.exonN, _overlap_))
            OUT_prefix = OUT_prefix +'.{}.overlap{}'.format(self.exonN, _overlap_)
        else:
            print(std_MAIN_PROCESS_NAME + "Overlap: {}".format(_overlap_))
            OUT_prefix = OUT_prefix +'.overlap{}'.format(_overlap_)

        # self.IMP_Result_prefix = self.IMP_Result_prefix + '.overlap{}'.format(_overlap_)



        ### (2) IMPUTE

        t_raw_IMP_Reuslt = self.IMPUTE(OUT_prefix, MHC_QC_VCF, REF_PHASED_VCF, _BEAGLE4,
                                       _overlap=_overlap_, _aver_erate=_aver_erate, _Genetic_Map=_Genetic_Map)
        print("raw Imputed Reuslt : {}".format(t_raw_IMP_Reuslt))


        ### (3) CONVERT_OUT

        # temporary hard-coding
        # self.raw_IMP_Reuslt = 'tests/_3_CookHLA/20190529/_3_HM_CEU_T1DGC_REF.QC.imputation_out.vcf'

        t_IMP_Result = self.CONVERT_OUT(MHC, _reference, OUT_prefix, _VCF2BEAGLE, _BEAGLE2LINKAGE, _PLINK, _overlap=_overlap_,
                                        raw_IMP_Result=t_raw_IMP_Reuslt)
        print('t_IMP_Result : {}'.format(t_IMP_Result))


        return t_IMP_Result



    def ManualInnerJoin(self, _left, _right, _out):



        return 0