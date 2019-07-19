# -*- coding: utf-8 -*-

import os, sys, re
from os.path import join
import pandas as pd
import numpy as np

from src.RUN_Bash import RUN_Bash



########## < Core Varialbes > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
HLA_names2 = ["A", "B", "C", "DRB1", "DPA1", "DPB1", "DQA1", "DQB1"]
HLA_names_gen = ["A", "C", "B", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]
isClassI = {"A": True, "B": True, "C": True, "DPA1": False, "DPB1": False, "DQA1": False, "DQB1": False, "DRB1": False}
HLA_mid_position = {"A": 30019970, "B": 31431272, "C": 31346171, "DPA1": 33145064, "DPB1": 33157346, "DQA1": 32716284, "DQB1": 32739039, "DRB1": 32660042}


__EXON__ = ['exon2', 'exon3', 'exon4']
__overlap__ = [3000, 4000, 5000]


# Patterns
p_HLA = re.compile(r'^HLA_')
p_HLA2 = re.compile(r'^HLA_\w+_(\d{4,5})$')
p_suffix = re.compile(r'_exon\d$')



def HLA_Genotype_Call_v2(_dict_IMP_Result, _out=None):


    ###### < Variables > ######

    # Input
    dict_IMP_vcf = _dict_IMP_Result
    dict_IMP_vcf_LEFT = {_exonN: {_overlap: None for _overlap in __overlap__} for _exonN in __EXON__} # 0 ~ 8 columns (Marker information)
    dict_IMP_vcf_RIGHT = {_hla: [] for _hla in HLA_names} # 9 ~ columns (Patient genotype information.)

    # GP
    dict_IMP_vcf_RIGHT_mean = {_hla: None for _hla in HLA_names}


    # Called HLA (unphased) Genotype
    dict_HLA_Genotype = {_hla: None for _hla in HLA_names}

    # Output Prefix
    OUTPUT_dir = os.path.dirname(dict_IMP_vcf['exon2'][3000]) # Just using 1st one.



    ### Disintegrating each VCF file.
    for _exonN in __EXON__:
        for _overlap in __overlap__:

            # print("({}, {})".format(_exonN, _overlap))

            with open(dict_IMP_vcf[_exonN][_overlap], 'r') as f_vcf, \
                 open(dict_IMP_vcf[_exonN][_overlap]+'.body', 'w') as f_vcf_body:

                for line in f_vcf:

                    if not line.startswith('##'):
                        f_vcf_body.write(line)

            df_vcf_body = pd.read_csv(dict_IMP_vcf[_exonN][_overlap]+'.body', sep='\s+', header=0)
            # print(df_vcf_body.head())

            RUN_Bash('rm {}'.format(dict_IMP_vcf[_exonN][_overlap]+'.body'))



            ### Filtering out HLA allele binary markers.

            flag_HLA = df_vcf_body.iloc[:, 2].str.match(p_HLA)
            df_vcf_body_HLA = df_vcf_body.loc[flag_HLA, :]
            # print(df_vcf_body_HLA.head())

            dict_IMP_vcf_LEFT[_exonN][_overlap] = df_vcf_body_HLA.iloc[:, :9] # Left part of vcf file.

            # Columns to be used as Indexes
            ID = df_vcf_body_HLA.iloc[:, 2].apply(lambda x : p_suffix.sub(repl='', string=x))
            # print("Trimmed ID : \n{}".format(ID.head()))
            exon = pd.Series([_exonN for z in range(0, len(ID))], name='exon', index=ID.index)
            # print(exon.head())
            overlap = pd.Series([_overlap for z in range(0, len(ID))], name='overlap', index=ID.index)
            # print(overlap.head())

            df_idx = pd.concat([ID, exon, overlap], axis=1).set_index('ID')
            # print(df_idx)



            ### Extrating 'GP(Genotpe Probability)' or 'Dose' values.

            needed_columns = [i for i in range(9, df_vcf_body.shape[1])]

            df_vcf_body_HLA = pd.concat([ID, df_vcf_body_HLA.iloc[:, needed_columns]], axis=1).set_index('ID').astype(str).applymap(lambda x : x.split(':')[2].split(','))
            df_vcf_body_HLA = df_vcf_body_HLA.applymap(lambda x : list(map(float, x))).applymap(lambda x : x[0]+x[1]/2)

            df_vcf_body_HLA = pd.concat([df_idx, df_vcf_body_HLA], axis=1)
            # print("DataFrame of GP value :\n{}".format(df_vcf_body_HLA.head()))


            for _hla in HLA_names:

                if not isClassI[_hla] and _exonN == 'exon4':
                    continue

                df_vcf_body = df_vcf_body_HLA.filter(regex='^HLA_{}'.format(_hla), axis=0).reset_index('ID')

                if df_vcf_body.shape[0] > 0:
                    dict_IMP_vcf_RIGHT[_hla].append(df_vcf_body)





    ### Average posterior GP of 9 implementations.

    dict_IMP_vcf_RIGHT_byHLA = {_hla: pd.concat(dict_IMP_vcf_RIGHT[_hla], axis=0) for _hla in HLA_names}

    for _hla in HLA_names:

        l_HLA_mean = []

        for idx, df in dict_IMP_vcf_RIGHT_byHLA[_hla].groupby('ID'):

            # print("Group by 'ID':\n{}".format(df.head()))

            df_temp = df.iloc[:, 3:].mean()
            df_temp.name = idx
            # print("Avg of group by :\n{}".format(df_temp.head()))

            l_HLA_mean.append(df_temp)


        df_HLA_mean = pd.concat(l_HLA_mean, axis=1).transpose()
        df_HLA_mean.index.name = 'ID'
        # print("df_mean : \n{}".format(df_HLA_mean.head()))

        dict_IMP_vcf_RIGHT_mean[_hla] = df_HLA_mean
        df_HLA_mean.to_csv(join(OUTPUT_dir, 'Avg_Posterior.{}.txt'.format(_hla)), sep='\t', header=True, index=True)


        df_called = df_HLA_mean.apply(lambda x : Main_Call(x), axis=0)
        dict_HLA_Genotype[_hla] = df_called
        # print(df_called.head())




    ### HLA Genotype Call

    df_called = pd.DataFrame.from_dict(dict_HLA_Genotype)
    # print(df_called.head())


    _out = join(OUTPUT_dir, _out) if _out else join(OUTPUT_dir, 'HLA_IMPUTATION_OUT.alleles')

    with open(_out, 'w') as f_out:

        count = 0

        for row in df_called.iterrows():

            PID = row[0]
            genotypes = row[1].to_dict()

            for _hla in HLA_names2:
                f_out.write('\t'.join([PID, PID, _hla, ',', genotypes[_hla]]) + '\n')

            count += 1
            # if count > 5 : break;


    return _out






def Main_Call(_sr_EachPatient):

    # print(_sr_EachPatient.head())

    the_sorted = _sr_EachPatient.loc[_sr_EachPatient > 0].sort_values(ascending=False)
    index_4digits = the_sorted.index.to_series().str.extract(p_HLA2, expand=False)
    the_sorted.index = index_4digits
    # print("\n================================================================")
    # print("Substted and sorted df :\n{}".format(the_sorted.head()))
    # print(index_4digits)


    if the_sorted.shape[0] == 0:
        # Error: Everything is 0.
        return ','

    elif the_sorted.shape[0] == 1:
        # Homozygous
        return ','.join([the_sorted.index.tolist()[0], the_sorted.index.tolist()[0]])

    elif the_sorted.shape[0] == 2:
        # Definite Heterozygous
        return ','.join(the_sorted.index.tolist())


    elif the_sorted.shape[0] > 2:

        flag_best = the_sorted == the_sorted.max()

        df_Best = the_sorted.loc[flag_best]
        df_Rest = the_sorted.loc[~flag_best]

        # print("\nBest :\n{}".format(df_Best))
        # print("\nRest :\n{}".format(df_Rest))


        if df_Best.shape[0] == 0:
            # Actually this case is impossible.
            return ',' # Error


        ### Normal case
        if df_Best.shape[0] == 1:
            # Main call
            # (1) The one in `df_Best` has bigger probability than 2 times of the best one in the rest. => Homozygous
            # (2) The one in `df_Best` and the best one in the rest => Heterozygous.

            best_of_best = df_Best.iat[0]
            best_of_rest = df_Rest.iat[0]

            if best_of_best > 2*best_of_rest:
                # (1) Homozygous
                return ','.join([df_Best.index.tolist()[0], df_Best.index.tolist()[0]])
            else:
                # (2) Heterozygous
                return ','.join([df_Best.index.tolist()[0], df_Rest.index.tolist()[0]])

        ### Normal case
        if df_Best.shape[0] == 2:
            # Heterozygous
            # Two different entities with same probability
            return ','.join(df_Best.index.tolist())




        if df_Best.shape[0] > 2:
            # Actually Error case
            # More than 2 HLA alleles have same(even) best avg_posterior probability
            return ','.join(df_Best.index.tolist())


    else:
        # Every HLA allele binary marker has 0 avg_posterior probability
        return ','









if __name__ == '__main__':

    """
    HLA_Genotype_Call_v2.py
    
    made by Wanson Choi.
    
    
    Another HLA genotype calling method.
    When Prephasing and Doubling techniques aren't used, HLA genotype will be called based on this script.
     
    
    """

    ### Example

    dict_raw_IMP = {_exonN : {_overlap : None for _overlap in __overlap__} for _exonN in __EXON__}

    # in Ubuntu
    # dict_raw_IMP['exon2'][3000] = '/home/wanson/Git_Projects/CookHLA/tests/HLA_Genotype_Call_v2/MM_AGM/HM_CEU_T1DGC_REF.MM.AGM.noprephasing.MHC.exon2.3000.QC.doubled.imputation_out.vcf'
    # dict_raw_IMP['exon2'][4000] = '/home/wanson/Git_Projects/CookHLA/tests/HLA_Genotype_Call_v2/MM_AGM/HM_CEU_T1DGC_REF.MM.AGM.noprephasing.MHC.exon2.4000.QC.doubled.imputation_out.vcf'
    # dict_raw_IMP['exon2'][5000] = '/home/wanson/Git_Projects/CookHLA/tests/HLA_Genotype_Call_v2/MM_AGM/HM_CEU_T1DGC_REF.MM.AGM.noprephasing.MHC.exon2.5000.QC.doubled.imputation_out.vcf'
    # dict_raw_IMP['exon3'][3000] = '/home/wanson/Git_Projects/CookHLA/tests/HLA_Genotype_Call_v2/MM_AGM/HM_CEU_T1DGC_REF.MM.AGM.noprephasing.MHC.exon3.3000.QC.doubled.imputation_out.vcf'
    # dict_raw_IMP['exon3'][4000] = '/home/wanson/Git_Projects/CookHLA/tests/HLA_Genotype_Call_v2/MM_AGM/HM_CEU_T1DGC_REF.MM.AGM.noprephasing.MHC.exon3.4000.QC.doubled.imputation_out.vcf'
    # dict_raw_IMP['exon3'][5000] = '/home/wanson/Git_Projects/CookHLA/tests/HLA_Genotype_Call_v2/MM_AGM/HM_CEU_T1DGC_REF.MM.AGM.noprephasing.MHC.exon3.5000.QC.doubled.imputation_out.vcf'
    # dict_raw_IMP['exon4'][3000] = '/home/wanson/Git_Projects/CookHLA/tests/HLA_Genotype_Call_v2/MM_AGM/HM_CEU_T1DGC_REF.MM.AGM.noprephasing.MHC.exon4.3000.QC.doubled.imputation_out.vcf'
    # dict_raw_IMP['exon4'][4000] = '/home/wanson/Git_Projects/CookHLA/tests/HLA_Genotype_Call_v2/MM_AGM/HM_CEU_T1DGC_REF.MM.AGM.noprephasing.MHC.exon4.4000.QC.doubled.imputation_out.vcf'
    # dict_raw_IMP['exon4'][5000] = '/home/wanson/Git_Projects/CookHLA/tests/HLA_Genotype_Call_v2/MM_AGM/HM_CEU_T1DGC_REF.MM.AGM.noprephasing.MHC.exon4.5000.QC.doubled.imputation_out.vcf'


    # in OS X
    # dict_raw_IMP['exon2'][3000] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190718_MM_Noprephasing/HM_CEU_T1DGC_REF.MM.AGM.noprephasing.MHC.exon2.3000.QC.doubled.imputation_out.vcf'
    # dict_raw_IMP['exon2'][4000] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190718_MM_Noprephasing/HM_CEU_T1DGC_REF.MM.AGM.noprephasing.MHC.exon2.4000.QC.doubled.imputation_out.vcf'
    # dict_raw_IMP['exon2'][5000] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190718_MM_Noprephasing/HM_CEU_T1DGC_REF.MM.AGM.noprephasing.MHC.exon2.5000.QC.doubled.imputation_out.vcf'
    # dict_raw_IMP['exon3'][3000] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190718_MM_Noprephasing/HM_CEU_T1DGC_REF.MM.AGM.noprephasing.MHC.exon3.3000.QC.doubled.imputation_out.vcf'
    # dict_raw_IMP['exon3'][4000] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190718_MM_Noprephasing/HM_CEU_T1DGC_REF.MM.AGM.noprephasing.MHC.exon3.4000.QC.doubled.imputation_out.vcf'
    # dict_raw_IMP['exon3'][5000] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190718_MM_Noprephasing/HM_CEU_T1DGC_REF.MM.AGM.noprephasing.MHC.exon3.5000.QC.doubled.imputation_out.vcf'
    # dict_raw_IMP['exon4'][3000] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190718_MM_Noprephasing/HM_CEU_T1DGC_REF.MM.AGM.noprephasing.MHC.exon4.3000.QC.doubled.imputation_out.vcf'
    # dict_raw_IMP['exon4'][4000] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190718_MM_Noprephasing/HM_CEU_T1DGC_REF.MM.AGM.noprephasing.MHC.exon4.4000.QC.doubled.imputation_out.vcf'
    # dict_raw_IMP['exon4'][5000] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190718_MM_Noprephasing/HM_CEU_T1DGC_REF.MM.AGM.noprephasing.MHC.exon4.5000.QC.doubled.imputation_out.vcf'


    HLA_Genotype_Call_v2(dict_raw_IMP)


    ### Example `Main_Call()` function
    # df_temp = pd.read_csv('/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190718_BOTH_Noprephasing/Avg_Posterior.DRB1.txt', sep='\t', header=0, index_col=0)
    # print(df_temp.head(20))

    # the_Patient = df_temp.iloc[:, 0]
    # call = Main_Call(the_Patient)
    # print(call)

    # df_temp2 = df_temp.apply(lambda x : Main_Call(x), axis=0)
    # print(df_temp2)
    # df_temp2.to_csv('/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190718_BOTH_Noprephasing/Called.DRB1.txt', sep='\t', header=False, index=True)