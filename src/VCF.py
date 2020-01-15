#-*- coding: utf-8 -*-

import os, sys, re
# import subprocess
# from shutil import which
# import random
import pandas as pd


std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]


class VCF():

    def __init__(self, _vcf):

        self.vcf_file = _vcf
        [self.l_headers, self.df_vcf, self.idx_FORMAT] = self.readVCF(_vcf)

        print(self.df_vcf.head(20))
        print(self.idx_FORMAT)




    def readVCF(self, _vcf):

        with open(_vcf) as f_vcf:

            l_header = []
            l_colnames = None
            # l_lines = []

            for line in f_vcf:
                if line.startswith('##'):
                    l_header.append(line)
                elif line.startswith('#'):
                    l_colnames = re.split('\s+', line.rstrip('\n'))
                    break

            df_vcf = pd.DataFrame([re.split('\s+', line.rstrip('\n')) for line in f_vcf], columns=l_colnames)

        return [l_header, df_vcf, l_colnames.index('FORMAT')]



if __name__ == "__main__":

    # print("Testing ")

    _vcf = '/media/sf_VirtualBox_Share/HM_CEU-T1DGC_REF/HM_CEU_T1DGC_REF.MM.AGM.MHC.QC.exon2.3000.raw_imputation_out.vcf'

    myVCF = VCF(_vcf)
