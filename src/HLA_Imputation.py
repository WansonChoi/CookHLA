#-*- coding: utf-8 -*-

import os, sys, re
from os.path import join

from src.GC_tricked_bgl2ori_bgl import GCtricedBGL2OriginalBGL


########## < Core Varialbes > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
HLA_names_gen = ["A", "C", "B", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]

# __overlap__ = [3000, 4000, 5000]
__overlap__ = [3000]


class HLA_Imputation(object):

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

        # 'CONVERT_OUT'
        self.refined_REF_markers = None

        # Result
        self.OUTPUT_dir = os.path.dirname(_out)
        self.OUTPUT_dir_ref = join(self.OUTPUT_dir, os.path.basename(_reference))

        self.raw_IMP_Reuslt = None
        self.IMP_Result_prefix = _out
        self.IMP_Result = None  # '*.imputed.alleles'

        self.accuracy = 0







        ###### < Main - 'CONVERT_IN', 'IMPUTE', 'CONVERT_OUT' > ######


        ### (1) CONVERT_IN

        [MHC_QC_VCF, REF_PHASED_VCF] = self.CONVERT_IN(MHC, _reference, _out, _hg,
                                                       _LINKAGE2BEAGLE, _BEAGLE2VCF, _PLINK, _BEAGLE4,
                                                       _aver_erate=_aver_erate, _Genetic_Map=_Genetic_Map)

        # self.CONVERT_IN(MHC, _reference, _out, _hg, _LINKAGE2BEAGLE, _BEAGLE2VCF, _PLINK, _BEAGLE4,
        #                 _aver_erate=_aver_erate, _Genetic_Map=_Genetic_Map)

        # print("Convert_IN :\n{}\n{}".format(MHC_QC_VCF, REF_PHASED_VCF))




        if f_useGeneticMap:

            for _overlap_ in __overlap__:

                if f_useMultipleMarkers:
                    print(std_MAIN_PROCESS_NAME + "exonN: {} / Overlap: {}".format(_exonN_, _overlap_))
                else:
                    print(std_MAIN_PROCESS_NAME + "Overlap: {}".format(_overlap_))

                self.IMP_Result_prefix = self.IMP_Result_prefix + '.overlap{}'.format(_overlap_)


                # Temporary Hard coding
                # MHC_QC_VCF = 'tests/_3_CookHLA/20190523/_3_HM_CEU_T1DGC_REF.MHC.QC.phasing_out_not_double.doubled.vcf'
                # REF_PHASED_VCF = 'tests/_3_CookHLA/20190523/T1DGC_REF.phased.vcf'

                ### (2) IMPUTE

                # self.raw_IMP_Reuslt = self.IMPUTE(_out, MHC_QC_VCF, REF_PHASED_VCF, _BEAGLE4, _aver_erate, _Genetic_Map)
                # print("raw Imputed Reuslt : {}".format(self.raw_IMP_Reuslt))


                ### (3) CONVERT_OUT

                # Temporary Hard-coding
                # self.raw_IMP_Reuslt = 'tests/_3_CookHLA/20190523/_3_HM_CEU_T1DGC_REF.QC.doubled.imputation_out.vcf'
                # self.refiend_REF_markers = 'tests/_3_CookHLA/20190523/T1DGC_REF.refined.markers'

                # self.IMP_Result = self.CONVERT_OUT(MHC, _reference, _out, _VCF2BEAGLE, _BEAGLE2LINKAGE, _PLINK)




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
            elif isinstance(self.IMP_Result, list):
                self.accuracy = [measureAccuracy(_answer, item, 'all', outfile=_out+'.accuracy') for item in self.IMP_Result]






    def CONVERT_IN(self, MHC, _reference, _out, _hg, _LINKAGE2BEAGLE, _BEAGLE2VCF, _PLINK, _BEAGLE4,
                   _aver_erate=None, _Genetic_Map=None):


        __MHC_exonN__ = MHC if not self.f_useMultipleMarkers else MHC+'.{}'.format(self.exonN)

        print("[{}] Converting data to beagle format.".format(self.idx_process))

        command = ' '.join(
            [_LINKAGE2BEAGLE, 'pedigree={}'.format(MHC + '.QC.nopheno.ped'), 'data={}'.format(MHC + '.QC.dat'),
             'beagle={}'.format(__MHC_exonN__ + '.QC.bgl'), 'standard=true', '>', __MHC_exonN__ + '.QC.bgl.log'])  # Making '*.bgl' file.
        # print(command)
        os.system(command)

        if not self.__save_intermediates:
            # os.system(' '.join(['rm', MHC + '.QC.nopheno.ped']))
            # os.system(' '.join(['rm', MHC + '.QC.dat']))
            os.system('rm {}'.format(__MHC_exonN__ + '.QC.bgl.log'))


        ### Converting data to reference_markers_Position (Dispersing same genomic position of some markers.)

        from src.redefineBPv1BH import redefineBP

        RefinedMarkers = redefineBP(_reference + '.markers', self.OUTPUT_dir_ref+'.refined.markers')
        self.refined_REF_markers = RefinedMarkers # => This will be used in 'CONVERT_OUT'.


        ### Converting data to target_markers_Position and extract not_including snp.

        command = ' '.join(['awk \'{print $2" "$4" "$5" "$6}\'', MHC + '.QC.bim', '>',
                            __MHC_exonN__ + '.QC.markers'])  # Making '*.markers' file.
        # print(command)
        os.system(command)

        # if not self.__save_intermediates:
        #     os.system(' '.join(['rm', MHC + '.QC.{bed,bim,fam,log}']))

        command = ' '.join(['Rscript src/excluding_snp_and_refine_target_position-v1COOK02222017.R',
                            __MHC_exonN__ + '.QC.markers', RefinedMarkers, __MHC_exonN__ + '.QC.pre.markers'])
        # print(command)
        os.system(command)

        if not self.__save_intermediates:
            os.system(' '.join(['rm', __MHC_exonN__ + '.QC.markers']))

        command = ' '.join(['mv', __MHC_exonN__ + '.QC.bgl', __MHC_exonN__ + '.QC.pre.bgl.phased'])
        # print(command)
        os.system(command)

        command = ' '.join(
            ["awk '{print $1}'", __MHC_exonN__ + '.QC.pre.markers', '>', os.path.join(self.OUTPUT_dir, 'selected_snp.txt')])
        # print(command)
        os.system(command)

        from src.Panel_subset import Panel_Subset
        qc_refined = Panel_Subset(__MHC_exonN__ + '.QC.pre', 'all', join(self.OUTPUT_dir, 'selected_snp.txt'), __MHC_exonN__ + '.QC.refined')
        # print(qc_refined) # Refined Beagle files are generated here.

        if not self.__save_intermediates:
            os.system(' '.join(['rm', __MHC_exonN__ + '.QC.pre.{bgl.phased,markers}']))
            os.system(' '.join(['rm', join(self.OUTPUT_dir, 'selected_snp.txt')]))


        ### Converting data to GC_change_beagle format.

        from src.bgl2GC_trick_bgl import Bgl2GC

        # target
        [GCchangeBGL, GCchangeMarkers] = Bgl2GC(__MHC_exonN__ + '.QC.refined.bgl.phased', __MHC_exonN__ + '.QC.refined.markers',
                                                __MHC_exonN__ + '.QC.GCchange.bgl', __MHC_exonN__ + '.QC.GCchange.markers')
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
            os.system(' '.join(['rm', __MHC_exonN__ + '.QC.refined.{bgl.phased,markers}']))

            if self.f_useMultipleMarkers:
                os.system(' '.join(['rm', RefinedMarkers])) # => This will be used in 'CONVERT_OUT" when not using Multiple Markers.


        ### Converting data to vcf_format

        # target
        command = ' '.join([_BEAGLE2VCF, '6', GCchangeMarkers, GCchangeBGL, '0', '>', __MHC_exonN__ + '.QC.vcf'])
        # print(command)
        os.system(command)

        MHC_QC_VCF = __MHC_exonN__ + '.QC.vcf'


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
            if self.f_useMultipleMarkers:
                os.system(' '.join(['rm {}'.format(GCchangeBGL)])) # 'GCchangeBGL' will be used in 'CONVERT_OUT'
            os.system(' '.join(['rm {}'.format(GCchangeMarkers)]))
            os.system(' '.join(['rm {}'.format(GCchangeBGL_REF)]))
            os.system(' '.join(['rm {}'.format(GCchangeMarkers_REF)]))

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

            command = 'awk \'{print $1" "$2" "$3}\' %s > %s' % (_Genetic_Map, _Genetic_Map+'.first')
            print(command)
            os.system(command)

            command = 'awk \'{print $2}\' %s > %s' % (self.OUTPUT_dir_ref + '.GCchange.markers', _Genetic_Map+'.second')
            print(command)
            os.system(command)

            command = 'paste -d " " {} {} > {}'.format(_Genetic_Map+'.first', _Genetic_Map+'.second', _Genetic_Map+'.refined.map')
            print(command)
            os.system(command)


            if os.path.exists(_Genetic_Map+'.refined.map'):

                self.refined_Genetic_Map = _Genetic_Map+'.refined.map'

                if not self.__save_intermediates:
                    os.system('rm {}'.format(_Genetic_Map+'.first'))
                    os.system('rm {}'.format(_Genetic_Map+'.second'))

            else:
                print(std_ERROR_MAIN_PROCESS_NAME + "Failed to generate Refined Genetic Map.")
                sys.exit()




        ### Mutliple Markers

        if not self.f_useMultipleMarkers:

            __RETURN__ = [MHC_QC_VCF, REF_PHASED_VCF]

        else:

            ### Phasing & Doubling (only on Target Sample.)

            # Phasing
            PHASED_RESULT = self.Phasing(__MHC_exonN__, MHC_QC_VCF, REF_PHASED_VCF, _BEAGLE4)

            # Doubling
            DOUBLED_PHASED_RESULT = self.Doubling(PHASED_RESULT)


            __RETURN__ = [DOUBLED_PHASED_RESULT, REF_PHASED_VCF]



        self.idx_process += 1

        return __RETURN__



    def IMPUTE(self, _out, _Doubled_VCF, _REF_PHASED_VCF, _BEAGLE4, _overlap=None, _aver_erate=None, _Genetic_Map=None):


        print("[{}] Performing HLA imputation (see {}.MHC.QC.imputation_out.log for progress).".format(self.idx_process, _out))


        OUT = _out + '.QC.doubled.imputation_out'


        if self.f_useGeneticMap:

            ### Using both 'Multiple Markers' and 'Adaptive Genetic Map'.

            """
            java -jar beagle4.jar gt=$MHC.QC.phasing_out_double.vcf ref=$REFERENCE.phased.vcf out=$MHC.QC.double.imputation_out impute=true lowmem=true 
                                    gprobs=true ne=10000 overlap=5000 err=$aver_erate map=$geneticMap.refined.map
            """

            with open(_aver_erate, 'r') as f:
                aver_erate = f.readline().rstrip('\n')

            command = '{} gt={} ref={} out={} impute=true lowmem=true gprobs=true ne=10000 overlap={} err={} map={}'.format(_BEAGLE4, _Doubled_VCF, _REF_PHASED_VCF, OUT, _overlap, aver_erate, self.refined_Genetic_Map)
            print(command)
            if not os.system(command):
                if not self.__save_intermediates:
                    os.system('rm {}'.format(OUT+'.log'))
                    os.system('rm {}'.format(_Doubled_VCF))
                    os.system('rm {}'.format(_REF_PHASED_VCF))
            else:
                print(std_ERROR_MAIN_PROCESS_NAME + "Imputation with Geneticmap failed.")
                sys.exit()


        else:

            """
            java -jar beagle4.jar gt=$MHC.QC.phasing_out_double.vcf ref=$REFERENCE.phased.vcf out=$MHC.QC.double.imputation_out impute=true lowmem=true 
            """

            command = '{} gt={} ref={} out={} impute=true lowmem=true'.format(_BEAGLE4, _Doubled_VCF, _REF_PHASED_VCF, OUT)
            # print(command)
            if not os.system(command):
                if not self.__save_intermediates:
                    # os.system(' '.join(['rm', OUT + '.log'])) # Imputation Log file will be saved.
                    os.system(' '.join(['rm', _Doubled_VCF]))
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



    def CONVERT_OUT(self, MHC, _reference, _out, _VCF2BEAGLE, _BEAGLE2LINKAGE, _PLINK):


        if not self.f_useMultipleMarkers:

            __MHC_exonN__ = MHC+'.{}'.format(self.exonN) if self.f_useMultipleMarkers else MHC

            ### Converting imputation result in vcf file to beagle format.

            command = 'cat {} | {} 0 {}'.format(self.raw_IMP_Reuslt, _VCF2BEAGLE, __MHC_exonN__ + '.QC.imputation_GCchange')
            # print(command)
            if not os.system(command):
                if not self.__save_intermediates:
                    os.system('rm {}'.format(__MHC_exonN__ + '.QC.imputation_GCchange.int'))

            command = 'gzip -d -f {}'.format(__MHC_exonN__ + '.QC.imputation_GCchange.bgl.gz')
            # print(command)
            os.system(command)


            ### Decoding GC-encoded values in beagle file to original values.

            # Decoding GC-encoding
            GC_decodedBGL = GCtricedBGL2OriginalBGL(__MHC_exonN__ + '.QC.imputation_GCchange.bgl', self.refined_REF_markers, __MHC_exonN__ + '.QC.imputation_ori.bgl')
            # print('GC_decodedBGL : {}'.format(GC_decodedBGL))

            if not self.__save_intermediates:
                os.system('rm {}'.format(self.refined_REF_markers))

            #
            HLA_IMPUTED_Result = os.path.join(self.OUTPUT_dir, 'HLA_IMPUTED_Result.{}'.format(os.path.basename(__MHC_exonN__)))

            command = 'Rscript src/complete_header.R {} {} {}'.format(__MHC_exonN__ + '.QC.GCchange.bgl', GC_decodedBGL, HLA_IMPUTED_Result + '.bgl.phased')
            print(command)
            os.system(command)


            ### Converting decoded beagle file to PLINK file.

            # command = 'cat {} | {} {}'.format(HLA_IMPUTED_Result, _BEAGLE2LINKAGE, _out+'.tmp')
            # print(command)
            # os.system(command)
            #
            # command = "cut -d ' ' -f6- {} > {}".format(_out+'.tmp.ped', _out+'.tmp')
            # print(command)
            # os.system(command)
            #
            # command = "paste -d ' ' {} {} | tr -d '\015' > {}".format(__MHC_exonN__ + '.fam', _out + '.tmp', HLA_IMPUTED_Result + '.ped') # *.ped
            # print(command)
            # os.system(command)
            #
            # command = 'cut -f1-4 {} > {}'.format(_reference+'.bim', HLA_IMPUTED_Result+'.map') # *.map
            # print(command)
            # os.system(command)
            #
            # command = 'cp {} {}'.format(__MHC_exonN__ + '.fam', HLA_IMPUTED_Result + '.fam')
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

            __RETURN__ = BGL2Alleles(HLA_IMPUTED_Result + '.bgl.phased', self.IMP_Result_prefix+'.imputed.alleles', 'all')




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