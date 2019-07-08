#-*- coding: utf-8 -*-

import os, sys, re
from os.path import join

from src.BGL2Alleles import BGL2Alleles
from src.GC_tricked_bgl2ori_bgl import GCtricedBGL2OriginalBGL
from src.RUN_Bash import RUN_Bash
from src.measureAccuracy import measureAccuracy



########## < Core Varialbes > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
HLA_names_gen = ["A", "C", "B", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]



class HLA_Imputation_GM(object):

    def __init__(self, idx_process, MHC, _reference, _out, _hg, _AdaptiveGeneticMap, _Average_Erate,
                 _LINKAGE2BEAGLE, _BEAGLE2LINKAGE, _BEAGLE2VCF, _VCF2BEAGLE, _PLINK, _BEAGLE4,
                 _answer=None, f_save_intermediates=False):


        ### Class variables

        # General
        self.idx_process = idx_process
        self.__save_intermediates = f_save_intermediates

        # Prefixes
        self.OUTPUT_dir = os.path.dirname(_out)
        self.OUTPUT_dir_ref = join(self.OUTPUT_dir, os.path.basename(_reference))
        self.OUTPUT_dir_GM = join(self.OUTPUT_dir, os.path.basename(_AdaptiveGeneticMap))
        self.HLA_IMPUTED_Result_MHC = join(os.path.dirname(MHC), "HLA_IMPUTED_Result.{}".format(os.path.basename(MHC)))

        # Result
        self.raw_IMP_Reuslt = None
        self.IMP_Result = None  # Final Imputation output ('*.imputed.alleles').
        self.accuracy = None

        # Dependencies
        self.LINKAGE2BEAGLE = _LINKAGE2BEAGLE
        self.BEAGLE2LINKAGE = _BEAGLE2LINKAGE
        self.BEAGLE2VCF = _BEAGLE2VCF
        self.VCF2BEAGLE = _VCF2BEAGLE
        self.PLINK = _PLINK
        self.BEAGLE4 = _BEAGLE4

        # Adaptive Genetic Map
        self.__AGM__ = _AdaptiveGeneticMap
        self.__AVER__ = _Average_Erate

        # created in 'CONVERT_IN'
        self.refined_REF_markers = None # used in 'CONVERT_OUT'
        self.refined_Genetic_Map = None # used in 'IMPUTE'
        self.GCchangeBGL = None # used in 'CONVERT_OUT'



        ###### < Main - 'CONVERT_IN', 'IMPUTE', 'CONVERT_OUT' > ######


        ### (1) CONVERT_IN

        # self.CONVERT_IN(MHC, _reference, _out, _hg, _aver_erate=self.__AVER__, _Genetic_Map=self.__AGM__)
        [MHC_QC_VCF, REF_PHASED_VCF] = self.CONVERT_IN(MHC, _reference, _out, _hg, _aver_erate=self.__AVER__, _Genetic_Map=self.__AGM__)

        # [Temporary Hard coding]
        # MHC_QC_VCF = "/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190605_onlyAGM/_3_HM_CEU_T1DGC_REF.MHC.QC.vcf"
        # REF_PHASED_VCF = "/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190605_onlyAGM/T1DGC_REF.phased.vcf"
        # self.refined_Genetic_Map = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190605_onlyAGM/CEU_T1DGC.mach_step.avg.clpsB.refined.map'
        # print("CONVERT_IN :\n{}\n{}".format(MHC_QC_VCF, REF_PHASED_VCF))



        ### (2) IMPUTE

        self.raw_IMP_Reuslt = self.IMPUTE(_out, MHC_QC_VCF, REF_PHASED_VCF, self.__AVER__, self.refined_Genetic_Map)

        # [Temporary Hard coding]
        # self.raw_IMP_Reuslt = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190605_onlyAGM/_3_HM_CEU_T1DGC_REF.QC.imputation_out.vcf'
        # self.refined_REF_markers = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190605_onlyAGM/T1DGC_REF.refined.markers'
        # self.GCchangeBGL = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190605_onlyAGM/_3_HM_CEU_T1DGC_REF.MHC.QC.GCchange.bgl'
        # print("raw Imputed Reuslt :\n{}".format(self.raw_IMP_Reuslt))


        ### (3) CONVERT_OUT

        self.IMP_Result = self.CONVERT_OUT(MHC, _reference, _out, self.raw_IMP_Reuslt)

        # [Temporary Hard coding]
        # self.IMP_Result = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190605_onlyAGM/HLA_IMPUTED_Result._3_HM_CEU_T1DGC_REF.MHC.imputed.alleles'
        # print("\n\nImputation Result : {}".format(self.IMP_Result))



        ###### < Get Accuracy > ######

        self.accuracy = self.getAccuracy(self.IMP_Result, _answer, self.HLA_IMPUTED_Result_MHC+'.imputed.alleles.accuracy')
        # print("Accuracy result :\n{}".format(self.accuracy))





    def CONVERT_IN(self, MHC, _reference, _out, _hg, _aver_erate, _Genetic_Map):


        print("[{}] Converting data to beagle format.".format(self.idx_process))

        RUN_Bash(self.LINKAGE2BEAGLE + ' pedigree={} data={} beagle={} standard=true > {}'.format(
            MHC + '.QC.nopheno.ped', MHC + '.QC.dat', MHC + '.QC.bgl', _out+'.bgl.log'))

        if not self.__save_intermediates:
            os.system('rm {}'.format(MHC + '.QC.nopheno.ped'))
            os.system('rm {}'.format(MHC + '.QC.dat'))
            os.system('rm {}'.format(_out+'.bgl.log'))




        ### Converting data to reference_markers_Position (Dispersing same genomic position of some markers.)

        from src.redefineBPv1BH import redefineBP

        RefinedMarkers = redefineBP(_reference + '.markers', self.OUTPUT_dir_ref+'.refined.markers')
        self.refined_REF_markers = RefinedMarkers # => This will be used in 'CONVERT_OUT'.




        ### Converting data to target_markers_Position and extract not_including snp.

        RUN_Bash('awk \'{print $2" "$4" "$5" "$6}\' %s > %s' % (MHC + '.QC.bim', MHC + '.QC.markers'))

        RUN_Bash('Rscript src/excluding_snp_and_refine_target_position-v1COOK02222017.R {} {} {}'.format(
            MHC+'.QC.markers', RefinedMarkers, MHC+'.QC.pre.markers'
        ))
        if not self.__save_intermediates:
            os.system(' '.join(['rm', MHC + '.QC.markers']))

        RUN_Bash('mv {} {}'.format(MHC+'.QC.bgl', MHC+'.QC.pre.bgl.phased'))

        RUN_Bash("awk '{print $1}' %s > %s" % (MHC+'.QC.pre.markers', join(self.OUTPUT_dir, 'selected_snp.txt')))


        from src.Panel_subset import Panel_Subset
        qc_refined = Panel_Subset(MHC + '.QC.pre', 'all', join(self.OUTPUT_dir, 'selected_snp.txt'), MHC + '.QC.refined')

        if not self.__save_intermediates:
            RUN_Bash('rm {}'.format(MHC + '.QC.pre.bgl.phased'))
            RUN_Bash('rm {}'.format(MHC + '.QC.pre.markers'))
            RUN_Bash('rm {}'.format(join(self.OUTPUT_dir, 'selected_snp.txt')))




        ### Converting data to GC_change_beagle format.

        from src.bgl2GC_trick_bgl import Bgl2GC

        # target
        [GCchangeBGL, GCchangeMarkers] = Bgl2GC(MHC + '.QC.refined.bgl.phased', MHC + '.QC.refined.markers',
                                                MHC + '.QC.GCchange.bgl', MHC + '.QC.GCchange.markers')

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

            RUN_Bash('rm {}'.format(MHC + '.QC.refined.bgl.phased'))
            RUN_Bash('rm {}'.format(MHC + '.QC.refined.markers'))
            # RUN_Bash('rm {}'.format(RefinedMarkers))

            # os.system(' '.join(['rm', RefinedMarkers])) # => This will be used in 'CONVERT_OUT" when not using Multiple Markers.




        ### Converting data to vcf_format

        # target
        RUN_Bash(self.BEAGLE2VCF + ' 6 {} {} 0 > {}'.format(GCchangeMarkers, GCchangeBGL, MHC+'.QC.vcf'))

        MHC_QC_VCF = MHC + '.QC.vcf'


        # reference
        RUN_Bash(self.BEAGLE2VCF + ' 6 {} {} 0 > {}'.format(GCchangeMarkers_REF, GCchangeBGL_REF, self.OUTPUT_dir_ref + '.vcf'))

        reference_vcf = self.OUTPUT_dir_ref + '.vcf'




        ### Converting data to reference_phased

        RUN_Bash('sed "s%/%|%g" {} > {}'.format(reference_vcf, self.OUTPUT_dir_ref + '.phased.vcf'))

        REF_PHASED_VCF = self.OUTPUT_dir_ref + '.phased.vcf'

        if not self.__save_intermediates:
            RUN_Bash('rm {}'.format(reference_vcf))

            # # if self.f_useMultipleMarkers:
            # if not self.f_useGeneticMap:
            #     os.system(' '.join(['rm {}'.format(GCchangeBGL)])) # 'GCchangeBGL' will be used in 'CONVERT_OUT'
            #     os.system(' '.join(['rm {}'.format(GCchangeMarkers_REF)]))  # 'GCchangeMarkers_REF' will be used in 'CONVERT_OUT'
            #     os.system(' '.join(['rm {}'.format(GCchangeMarkers)]))
            #     os.system(' '.join(['rm {}'.format(GCchangeBGL_REF)]))

        """
        (1) `MHC_QC_VCF` := MHC + '.QC.vcf',
        (2) `REF_PHASED_VCF` := self.OUTPUT_dir_ref + '.phased.vcf'

        These two files are to be passed into Beagle phasing;
        """




        ############### < Adaptive Genetic Map > ###############


        """
        awk '{print $1" "$2" "$3}' $geneticMap > $geneticMap.first
        awk '{print $2}' $REFERENCE.GCchange.markers > $geneticMap.second
        paste -d " " $geneticMap.first $geneticMap.second > $geneticMap.refined.map

        rm $geneticMap.first
        rm $geneticMap.second

        """

        REFINED_GENTIC_MAP = self.OUTPUT_dir_GM + '.refined.map'

        RUN_Bash('awk \'{print $1" "$2" "$3}\' %s > %s' % (_Genetic_Map, self.OUTPUT_dir_GM+'.first'))
        RUN_Bash('awk \'{print $2}\' %s > %s' % (GCchangeMarkers_REF, self.OUTPUT_dir_GM+'.second'))
        RUN_Bash('paste -d " " {} {} > {}'.format(self.OUTPUT_dir_GM+'.first', self.OUTPUT_dir_GM+'.second', REFINED_GENTIC_MAP))   # 이렇게 column bind시키는데는 당연히 *.first, *.second 파일의 row수가 같을 거라고 가정하는 상황.


        if os.path.exists(REFINED_GENTIC_MAP):

            self.refined_Genetic_Map = REFINED_GENTIC_MAP

            if not self.__save_intermediates:
                os.system('rm {}'.format(self.OUTPUT_dir_GM+'.first'))
                os.system('rm {}'.format(self.OUTPUT_dir_GM+'.second'))
                os.system('rm {}'.format(GCchangeMarkers_REF)) # (Genetic Map) *.GCchange.markers is removed here.

        else:
            print(std_ERROR_MAIN_PROCESS_NAME + "Failed to generate Refined Genetic Map.")
            sys.exit()



        __RETURN__ = [MHC_QC_VCF, REF_PHASED_VCF]



        self.idx_process += 1

        return __RETURN__



    def IMPUTE(self, _out, _MHC_QC_VCF, _REF_PHASED_VCF, _aver_erate, _Refined_Genetic_Map):


        print("[{}] Performing HLA imputation (see {}.MHC.QC.imputation_out.log for progress).".format(self.idx_process, _out))


        OUT = _out + '.QC.imputation_out'


        """
        beagle gt=$MHC.QC.vcf ref=$REFERENCE.phased.vcf out=$MHC.QC.imputation_out impute=true gprobs=true lowmem=true 
                    map=$geneticMap.refined.map ne=10000 overlap=5000 err=$aver_erate
        """

        with open(_aver_erate, 'r') as f:
            aver_erate = f.readline().rstrip('\n')

        command = '{} gt={} ref={} out={} impute=true gprobs=true lowmem=true map={} ne=10000 overlap=5000 err={} '.format(
            self.BEAGLE4, _MHC_QC_VCF, _REF_PHASED_VCF, OUT,
            _Refined_Genetic_Map, aver_erate)
        print(command)
        if not os.system(command):
            if not self.__save_intermediates:
                os.system('rm {}'.format(OUT+'.log'))
                # os.system('rm {}'.format(_MHC_QC_VCF))
                # os.system('rm {}'.format(_REF_PHASED_VCF)) # These 2 files
        else:
            print(std_ERROR_MAIN_PROCESS_NAME + "Imputation with Geneticmap failed.")
            sys.exit()


        RUN_Bash('gzip -d -f {}.vcf.gz'.format(OUT))

        __RETURN__ = OUT + '.vcf'
        self.idx_process += 1

        return __RETURN__



    def CONVERT_OUT(self, MHC, _reference, _out, raw_IMP_Result):


        ### Converting imputation result in vcf file to beagle format.

        RUN_Bash('cat {} | {} 0 {}'.format(raw_IMP_Result, self.VCF2BEAGLE, MHC + '.QC.imputation_GCchange'))

        if not self.__save_intermediates:
            os.system('rm {}'.format(MHC + '.QC.imputation_GCchange.int'))

        RUN_Bash('gzip -d -f {}'.format(MHC + '.QC.imputation_GCchange.bgl.gz'))



        ### Converting imputation GC_beagle to original beagle(Decoding GC-encoding).

        GC_decodedBGL = GCtricedBGL2OriginalBGL(MHC + '.QC.imputation_GCchange.bgl', self.refined_REF_markers, MHC+'.QC.imputation_ori.bgl')
        # print('GC_decodedBGL : {}'.format(GC_decodedBGL))

        if not self.__save_intermediates:
            os.system('rm {}'.format(MHC + '.QC.imputation_GCchange.bgl'))
            os.system('rm {}'.format(MHC + '.QC.imputation_GCchange.markers'))

        HLA_IMPUTED_Result_MHC = self.HLA_IMPUTED_Result_MHC

        RUN_Bash('Rscript src/complete_header.R {} {} {}'.format(self.GCchangeBGL, GC_decodedBGL, HLA_IMPUTED_Result_MHC+'.bgl.phased'))

        if not self.__save_intermediates:
            # os.system('rm {}'.format(MHC + '.QC.GCchange.bgl'))
            os.system('rm {}'.format(GC_decodedBGL))


        ### Converting imputation genotypes to PLINK .ped format.

        RUN_Bash('cat {} | {} {}'.format(HLA_IMPUTED_Result_MHC+'.bgl.phased', self.BEAGLE2LINKAGE, _out+'.tmp'))
        RUN_Bash("cut -d ' ' -f6- {} > {}".format(_out+'.tmp.ped', _out+'.tmp'))
        RUN_Bash("paste -d ' ' {} {} | tr -d '\015' > {}".format(MHC + '.fam', _out + '.tmp', HLA_IMPUTED_Result_MHC + '.ped')) # *.ped
        RUN_Bash('cut -f1-4 {} > {}'.format(_reference+'.bim', HLA_IMPUTED_Result_MHC+'.map')) # *.map
        RUN_Bash('cp {} {}'.format(MHC + '.fam', HLA_IMPUTED_Result_MHC + '.fam'))

        # Create PLINK bed format
        RUN_Bash('{} --ped {} --map {} --make-bed --out {}'.format(self.PLINK, HLA_IMPUTED_Result_MHC+'.ped', HLA_IMPUTED_Result_MHC+'.map', HLA_IMPUTED_Result_MHC))
        RUN_Bash('rm {}'.format(HLA_IMPUTED_Result_MHC + '.fam'))
        RUN_Bash('cp {} {}'.format(_reference+'.bim', HLA_IMPUTED_Result_MHC+'.bim'))


        if not self.__save_intermediates:
            RUN_Bash('rm {}'.format(_out+'.tmp'))
            RUN_Bash('rm {}'.format(_out+'.tmp.ped'))
            RUN_Bash('rm {}'.format(_out+'.tmp.dat'))

            RUN_Bash('rm {}'.format(self.raw_IMP_Reuslt))
            RUN_Bash('rm {}'.format(self.refined_REF_markers))
            RUN_Bash('rm {}'.format(self.GCchangeBGL))



        # BGL2Allele.py
        __RETURN__ = BGL2Alleles(HLA_IMPUTED_Result_MHC+'.bgl.phased', HLA_IMPUTED_Result_MHC+'.imputed.alleles', 'all')


        self.idx_process += 1
        return __RETURN__



    def getAccuracy(self, _IMT_Result, _answer, _out):

        if not _answer:
            return -1
        else:
            return measureAccuracy(_answer, _IMT_Result, 'all', _out)