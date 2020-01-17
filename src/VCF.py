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

    def __init__(self, _vcf, _p_dependency='./dependency', _mem='2g'):

        self.filepath = _vcf
        [self.l_headers, self.df_vcf, self.idx_FORMAT] = self.readVCF(_vcf)

        # print(self.df_vcf.head(20))
        # print(self.idx_FORMAT)




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



    def DoubleVCF(self, _out=None):


        df_Right = self.df_vcf.iloc[:, (self.idx_FORMAT + 1):]
        # print(df_Right)

        l_temp = []

        for col in df_Right.iteritems():

            df_temp = col[1].str.split('|', expand=True)

            sr_left_haplotype = df_temp.iloc[:, 0].map(lambda x : '|'.join([x, x]))
            sr_right_haplotype = df_temp.iloc[:, 1].map(lambda x : '|'.join([x, x]))

            df_columns_Doubled = pd.concat([sr_left_haplotype, sr_right_haplotype], axis=1)
            df_columns_Doubled.columns = [col[0], col[0]+'_1']

            l_temp.append(df_columns_Doubled)


        df_Right_Doubled = pd.concat(l_temp, axis=1)
        df_RETURN = pd.concat([self.df_vcf.iloc[:, :(self.idx_FORMAT + 1)], df_Right_Doubled], axis=1)


        if bool(_out):
            return self.Export_VCF(df_RETURN, _out)
        else:
            return df_RETURN




    def Export_VCF(self, _vcf, _out):

        if isinstance(_vcf, pd.DataFrame):

            with open(_out, 'w') as f_vcf_out:
                f_vcf_out.writelines(self.l_headers)

            _vcf.to_csv(_out, sep='\t', header=True, index=False, mode='a')

            return _out

        else:
            print(std_ERROR_MAIN_PROCESS_NAME + "Given '{}' is not DataFrame.".format(_vcf))
            return -1



    def SubsetSamples(self, _toKeep, _toRemove):
        
        return 0




if __name__ == "__main__":

    # print("Testing ")

    _vcf = '/media/sf_VirtualBox_Share/CookHLA/tests/T1DGC_CookQC/T1DGC_REF.ONLY_Variants_HLA.GCtrick.bgl.phased.SUBSET_for_TEST.vcf'

    myVCF = VCF(_vcf)
    d = myVCF.DoubleVCF(_out='/media/sf_VirtualBox_Share/CookHLA/tests/T1DGC_CookQC/T1DGC_REF.ONLY_Variants_HLA.GCtrick.bgl.phased.SUBSET_for_TEST.DOUBLED.vcf')
    print(d)
