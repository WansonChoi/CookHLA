#-*- coding: utf-8 -*-

import os, sys, re
from os.path import join

from src.GC_tricked_bgl2ori_bgl import GCtricedBGL2OriginalBGL
from src.RUN_Bash import RUN_Bash
from src.measureAccuracy import measureAccuracy
from src.redefineBPv1BH import redefineBP



########## < Core Varialbes > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
HLA_names_gen = ["A", "C", "B", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]

# __EXON__ = ['exon2', 'exon3', 'exon4']
__EXON__ = ['exon2']

# __overlap__ = [3000, 4000, 5000]
__overlap__ = [3000]



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
        self.raw_IMP_Reuslt = None
        self.IMP_Result = None  # Final Imputation output ('*.imputed.alleles').
        self.dict_IMP_Result = {_exonN: {_overlap: None} for _exonN in __EXON__ for _overlap in __overlap__}
        self.accuracy = None

        # Dependencies
        self.LINKAGE2BEAGLE = _LINKAGE2BEAGLE
        self.BEAGLE2LINKAGE = _BEAGLE2LINKAGE
        self.BEAGLE2VCF = _BEAGLE2VCF
        self.VCF2BEAGLE = _VCF2BEAGLE
        self.PLINK = _PLINK
        self.BEAGLE4 = _BEAGLE4

        # created in 'CONVERT_IN'
        self.refined_REF_markers = None # used in 'CONVERT_OUT'
        self.refined_Genetic_Map = None # used in 'IMPUTE'
        self.GCchangeBGL = None # used in 'CONVERT_OUT'







        ###### < Main - 'CONVERT_IN', 'IMPUTE', 'CONVERT_OUT' > ######

        if _MultP == 1:

            ## Single Core (No Multiprocessing)

            for _exonN in __EXON__:
                for _overlap in __overlap__:
                    self.dict_IMP_Result[_exonN][_overlap] = self.IMPUTATION_MM(_exonN, _overlap, MHC, _reference, _out, _hg)

        else:
            ## Multiprocessing

            # Under construction.

            pass







    def CONVERT_IN(self, MHC, _reference, _out, _hg):


        print("[{}] Converting data to beagle format.".format(self.idx_process))
        self.idx_process += 1


        RUN_Bash(self.LINKAGE2BEAGLE + ' pedigree={} data={} beagle={} standard=true > {}'.format(
            MHC + '.QC.nopheno.ped', MHC + '.QC.dat', MHC + '.QC.bgl', _out+'.bgl.log'))

        if not self.__save_intermediates:
            os.system('rm {}'.format(MHC + '.QC.nopheno.ped'))
            os.system('rm {}'.format(MHC + '.QC.dat'))
            os.system('rm {}'.format(_out+'.bgl.log'))




        ### Converting data to reference_markers_Position (Dispersing same genomic position of some markers.)

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



        ############### < Multiple Markers > ###############


        ### Phasing & Doubling (only on Target Sample.)

        # Phasing
        PHASED_RESULT = self.Phasing(MHC, MHC_QC_VCF, REF_PHASED_VCF)

        # [Temporary Hardcoding for Phased Result]
        # PHASED_RESULT = "/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190708_MM/_3_HM_CEU_T1DGC_REF.MHC.QC.phasing_out_not_double"
        # print("[Temporary Hardcoding]Phased Result:\n{}".format(PHASED_RESULT))


        # Doubling
        DOUBLED_PHASED_RESULT = self.Doubling(MHC, PHASED_RESULT)



        return [DOUBLED_PHASED_RESULT, REF_PHASED_VCF]



    def IMPUTE(self, _out, _MHC_QC_VCF, _REF_PHASED_VCF, _BEAGLE4, _overlap=None, _aver_erate=None, _Genetic_Map=None):


        print("[{}] Performing HLA imputation (see {}.MHC.QC.imputation_out.log for progress).".format(self.idx_process, _out))
        self.idx_process += 1


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


        RUN_Bash('cat {} {} > {}'.format(PHASED_RESULT + '.vcf.header1', PHASED_RESULT + '.Doubled.pre2.vcf', MHC+'.QC.phasing_out_doubled.vcf'))

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



    def IMPUTATION_MM(self, _exonN, _overlap, MHC, _reference, _out, _hg):

        ### (1) CONVERT_IN

        [DOUBLED_PHASED_RESULT, REF_PHASED_VCF] = self.CONVERT_IN(MHC, _reference, _out, _hg)


        # ### (2) IMPUTE
        #
        # self.raw_IMP_Reuslt = self.IMPUTE(_out, MHC_QC_VCF, REF_PHASED_VCF, _BEAGLE4)
        # print("raw Imputed Reuslt : {}".format(self.raw_IMP_Reuslt))
        #
        #
        # ### (3) CONVERT_OUT
        #
        # self.IMP_Result = self.CONVERT_OUT(MHC, _reference, _out, _VCF2BEAGLE, _BEAGLE2LINKAGE, _PLINK)
        # print("\n\nImputation Result : {}".format(self.IMP_Result))

        return 0



    def getAccuracy(self, _IMT_Result, _answer, _out):

        if not _answer:
            return -1
        else:
            return measureAccuracy(_answer, _IMT_Result, 'all', _out)