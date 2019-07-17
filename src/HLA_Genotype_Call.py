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
# HLA_names = ["A"]
HLA_names_gen = ["A", "C", "B", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]
isClassI = {"A": True, "B": True, "C": True, "DPA1": False, "DPB1": False, "DQA1": False, "DQB1": False, "DRB1": False}
HLA_mid_position = {"A": 30019970, "B": 31431272, "C": 31346171, "DPA1": 33145064, "DPB1": 33157346, "DQA1": 32716284, "DQB1": 32739039, "DRB1": 32660042}


__EXON__ = ['exon2', 'exon3', 'exon4']
__overlap__ = [3000, 4000, 5000]


# Patterns
p_HLA = re.compile(r'^HLA_')
p_HLA2 = re.compile(r'^HLA_\w+_(\d{4,5})$')
p_suffix = re.compile(r'_exon\d$')



def HLA_Genotype_Call(_dict_IMP_Result, _feature='GP', __forCheck=True, _out=None):

    if _feature != 'DS' and _feature != 'GP' and _feature != 'BOTH':
        print(std_ERROR_MAIN_PROCESS_NAME + "Wrong measure to calling genotype.")
        sys.exit()


    if _feature == 'BOTH':
        _feature = 'DSxGP'


    ###### < Variables > ######

    # Input
    dict_IMP_vcf = _dict_IMP_Result
    dict_IMP_vcf_LEFT = {_exonN: {_overlap: None for _overlap in __overlap__} for _exonN in __EXON__} # 0 ~ 8 columns (Marker information)

    # GP
    dict_IMP_vcf_RIGHT_GP_raw = {_hla: [] for _hla in HLA_names} # 9 ~ columns (Patient genotype information.)
    dict_IMP_vcf_RIGHT_mean = {_hla: None for _hla in HLA_names}

    # DS
    dict_IMP_vcf_RIGHT_DS_raw = {_hla: [] for _hla in HLA_names} # 9 ~ columns (Patient genotype information.)
    dict_IMP_vcf_RIGHT_DS = {_hla: {_exonN: {_overlap: None for _overlap in __overlap__} for _exonN in __EXON__} for _hla in HLA_names}
    dict_IMP_vcf_RIGHT_DS_Score = {_hla: None for _hla in HLA_names}

    # DS x GP
    dict_DSxGP = {_hla: None for _hla in HLA_names}

    n_row = 0
    n_col = 0

    PID = None

    # Called HLA (unphased) Genotype
    dict_HLA_Single_allele = {_hla: None for _hla in HLA_names}

    dict_forCheck_byHLA = {_hla: None for _hla in HLA_names}

    # Prefix
    OUTPUT_dir = os.path.dirname(dict_IMP_vcf['exon2'][3000])
    # _out = _out if _out else join(OUTPUT_dir, 'HLA_IMPUTATION_OUT.{}.single.alleles'.format(_feature))
    # _out_rawcall = join(OUTPUT_dir, 'Receipt.{}.raw_call.txt'.format(_feature))



    ### Disintegrating each VCF file.
    for _exonN in __EXON__:
        for _overlap in __overlap__:

            # print("({}, {})".format(_exonN, _overlap))

            with open(dict_IMP_vcf[_exonN][_overlap], 'r') as f_vcf, \
                 open(dict_IMP_vcf[_exonN][_overlap]+'.header', 'w') as f_vcf_header, \
                 open(dict_IMP_vcf[_exonN][_overlap]+'.body', 'w') as f_vcf_body:

                for line in f_vcf:

                    if line.startswith('##'):
                        f_vcf_header.write(line)
                    else:
                        f_vcf_body.write(line)

            df_temp = pd.read_csv(dict_IMP_vcf[_exonN][_overlap]+'.body', sep='\s+', header=0)

            RUN_Bash('rm {}'.format(dict_IMP_vcf[_exonN][_overlap]+'.header'))
            RUN_Bash('rm {}'.format(dict_IMP_vcf[_exonN][_overlap]+'.body'))

            ### Filtering out HLA allele binary markers.

            flag_HLA = df_temp.iloc[:, 2].str.match(p_HLA)
            df_onlyHLA = df_temp.loc[flag_HLA, :]

            dict_IMP_vcf_LEFT[_exonN][_overlap] = df_onlyHLA.iloc[:, :9]

            # Columns to be used as Indexes
            ID = df_onlyHLA.iloc[:, 2].apply(lambda x : p_suffix.sub(repl='', string=x))
            # print("Trimmed ID : \n{}".format(ID.head()))
            exon = pd.Series([_exonN for z in range(0, len(ID))], name='exon', index=ID.index)
            # print(exon.head())
            overlap = pd.Series([_overlap for z in range(0, len(ID))], name='overlap', index=ID.index)
            # print(overlap.head())

            df_idx = pd.concat([ID, exon, overlap], axis=1).set_index('ID')
            # print(df_idx)



            ### Extrating 'GP(Genotype Probability)' or 'Dose' values.

            needed_columns = [i for i in range(9, df_temp.shape[1])]


            if _feature == 'GP' or _feature == 'DSxGP':

                df_onlyHLA_GP = pd.concat([ID, df_onlyHLA.iloc[:, needed_columns]], axis=1).set_index('ID').astype(str).applymap(lambda x : x.split(':')[2].split(','))
                df_onlyHLA_GP = df_onlyHLA_GP.applymap(lambda x : list(map(float, x))).applymap(lambda x : x[0]+x[1]/2)

                df_onlyHLA_GP = pd.concat([df_idx, df_onlyHLA_GP], axis=1)
                # print("DataFrame of GP value :\n{}".format(df_onlyHLA_GP.head()))


                for _hla in HLA_names:

                    if not isClassI[_hla] and _exonN == 'exon4':
                        continue

                    df_temp = df_onlyHLA_GP.filter(regex='^HLA_{}'.format(_hla), axis=0).reset_index('ID')

                    if df_temp.shape[0] > 0:
                        dict_IMP_vcf_RIGHT_GP_raw[_hla].append(df_temp)



            if _feature == 'DS' or _feature == 'DSxGP'\
                    :
                """
                Cf) VCF file works focused on ALT(Alternative allele). Meanwhile, The binary marker 'P(Present)' in CookHLA is located in REF(Reference).
                So, dose value of given vcf file is 'A(Absent)' marker's value and we will extract min value which corresponds to 'P(Present)'.
                """

                df_onlyHLA_DS = pd.concat([ID, df_onlyHLA.iloc[:, needed_columns]], axis=1).set_index('ID').astype(str).applymap(lambda x : x.split(':')[1]).astype(float)
                df_onlyHLA_DS = pd.concat([df_idx, df_onlyHLA_DS], axis=1)
                # print("DataFrame of DS value :\n{}".format(df_onlyHLA_DS.head()))


                for _hla in HLA_names:

                    if not isClassI[_hla] and _exonN == 'exon4':
                        continue

                    df_temp = df_onlyHLA_DS.filter(regex='^HLA_{}'.format(_hla), axis=0).reset_index('ID')

                    if df_temp.shape[0] > 0:
                        dict_IMP_vcf_RIGHT_DS_raw[_hla].append(df_temp)




    if _feature == 'DS' or _feature == 'DSxGP':

        dict_IMP_vcf_RIGHT_byHLA = {_hla: pd.concat(dict_IMP_vcf_RIGHT_DS_raw[_hla], axis=0) for _hla in HLA_names}

        for _hla in HLA_names:

            # print("Subsetted by HLA : \n{}".format(dict_IMP_vcf_RIGHT_byHLA[_hla].head()))

            ### Coverting DS value to score

            l_Dose_Score = []

            count = 0

            for idx, df in dict_IMP_vcf_RIGHT_byHLA[_hla].groupby(['exon', 'overlap']):

                # print("<{} : {}> :\n{}".format(count, idx, df.head()))

                (_exonN, _overlap) = idx

                df_Dose_Score = df.set_index(['ID', 'exon', 'overlap']).apply(lambda x : x.apply(lambda y : 1 if y == x.min() else 0), axis=0)
                # print(df_Dose_Score.head(20))

                dict_IMP_vcf_RIGHT_DS[_hla][_exonN][_overlap] = df_Dose_Score


            # Dimension Info (Using the data of the object in last iteration.)
            n_row = df_Dose_Score.shape[0]
            n_col = df_Dose_Score.shape[1]

            new_index = df.loc[:, 'ID']
            new_columns = df_Dose_Score.columns


            ### Accumulating Dose score over 9 implementation

            arr_acc = np.zeros((n_row, n_col), dtype=int)
            # print("arr_acc :\n{}".format(arr_acc))

            for _exonN in __EXON__:

                if not isClassI[_hla] and _exonN == 'exon4':
                    continue

                for _overlap in __overlap__:
                    # print(dict_IMP_vcf_RIGHT_DS[_hla][_exonN][_overlap].to_numpy())
                    arr_acc = arr_acc + dict_IMP_vcf_RIGHT_DS[_hla][_exonN][_overlap].to_numpy()


            df_acc = pd.DataFrame(arr_acc, index=new_index, columns=new_columns)
            dict_IMP_vcf_RIGHT_DS_Score[_hla] = df_acc
            # print("Vote score accumulated : \n{}".format(df_acc))
            df_acc.to_csv(join(OUTPUT_dir, 'Receipt.DS.Vote_Score.{}.txt'.format(_hla)), sep='\t', header=True, index=True)




            ### For checking (Sorted by ('exon', 'overlap') key pair.)

            l_forCheck = []

            for _exonN in __EXON__:

                if not isClassI[_hla] and _exonN == 'exon4':
                    continue

                for _overlap in __overlap__:
                    l_forCheck.append(dict_IMP_vcf_RIGHT_DS[_hla][_exonN][_overlap])


            df_forCheck = pd.concat(l_forCheck, axis=0)
            # print(df_forCheck.head())
            dict_forCheck_byHLA[_hla] = df_forCheck
            df_forCheck.to_csv(join(OUTPUT_dir, 'Receipt.DS.each_Imputation.{}.txt'.format(_hla)), sep='\t', header=True, index=True)




        # GenotypeCalling(dict_IMP_vcf_RIGHT_DS_Score, _out, _out_rawcall)




    if _feature == 'GP' or _feature == 'DSxGP':

        ### Acquring mean posterior GP

        dict_IMP_vcf_RIGHT_byHLA = {_hla: pd.concat(dict_IMP_vcf_RIGHT_GP_raw[_hla], axis=0) for _hla in HLA_names}

        for _hla in HLA_names:

            l_HLA_mean = []

            for idx, df in dict_IMP_vcf_RIGHT_byHLA[_hla].groupby('ID'):

                sr_temp = df.iloc[:, 3:].mean()
                sr_temp.name = idx

                l_HLA_mean.append(sr_temp)


            df_HLA_mean = pd.concat(l_HLA_mean, axis=1).transpose()
            df_HLA_mean.index.name = 'ID'
            # print("df_mean : \n{}".format(df_HLA_mean.head()))

            dict_IMP_vcf_RIGHT_mean[_hla] = df_HLA_mean
            df_HLA_mean.to_csv(join(OUTPUT_dir, 'Receipt.GP.Avg_Posterior.{}.txt'.format(_hla)), sep='\t', header=True, index=True)



        # GenotypeCalling(dict_IMP_vcf_RIGHT_mean, _out, _out_rawcall)



    ### HLA Genotype call

    if _feature == 'DS':
        __RETURN__ = GenotypeCalling(dict_IMP_vcf_RIGHT_DS_Score, join(OUTPUT_dir, 'HLA_IMPUTATION_OUT.DS.single.alleles'), join(OUTPUT_dir, 'Receipt.DS.raw_call.txt'))
    elif _feature == 'GP':
        __RETURN__ = GenotypeCalling(dict_IMP_vcf_RIGHT_mean, join(OUTPUT_dir, 'HLA_IMPUTATION_OUT.GP.single.alleles'), join(OUTPUT_dir, 'Receipt.GP.raw_call.txt'))
    elif _feature == 'DSxGP':

        ### Element-wise multiplication

        for _hla in HLA_names:

            df_DS = dict_IMP_vcf_RIGHT_DS_Score[_hla]
            df_GP = dict_IMP_vcf_RIGHT_mean[_hla]


            DSxGP = np.multiply(df_DS.to_numpy(), df_GP.to_numpy())

            dict_DSxGP[_hla] = pd.DataFrame(DSxGP, index=df_DS.index, columns=df_DS.columns) # Just using index and columns of `df_DS`.
            dict_DSxGP[_hla].to_csv(join(OUTPUT_dir, 'Receipt.DSxGP.{}.txt'.format(_hla)), sep='\t', header=True, index=True)



        GenotypeCalling(dict_IMP_vcf_RIGHT_DS_Score, join(OUTPUT_dir, 'HLA_IMPUTATION_OUT.DS.single.alleles'), join(OUTPUT_dir, 'Receipt.DS.raw_call.txt'))
        GenotypeCalling(dict_IMP_vcf_RIGHT_mean, join(OUTPUT_dir, 'HLA_IMPUTATION_OUT.GP.single.alleles'),join(OUTPUT_dir, 'Receipt.GP.raw_call.txt'))
        __RETURN__ = GenotypeCalling(dict_DSxGP, join(OUTPUT_dir, 'HLA_IMPUTATION_OUT.DSxGP.single.alleles'), join(OUTPUT_dir, 'Receipt.DSxGP.raw_call.txt'))




    return __RETURN__





def GenotypeCalling(_dict_aggregated, _out, _out_rawcall):


    _dict_HLA_Single_allele = {_hla: None for _hla in HLA_names}


    ### Call HLA genotype.

    with open(_out_rawcall, 'w') as f_out_rawcall:

        for _hla in HLA_names:

            # print("{} :".format(_hla))
            f_out_rawcall.write("{} :\n".format(_hla))

            curr_DOSE_Score = _dict_aggregated[_hla]
            # print(curr_DOSE_Score.head())

            PID = curr_DOSE_Score.columns.tolist()
            dict_Call = {}

            for j in range(0, curr_DOSE_Score.shape[1], 2):

                sr_chr_1st = curr_DOSE_Score.iloc[:, j]
                sr_chr_2nd = curr_DOSE_Score.iloc[:, j + 1]

                flag_max1 = sr_chr_1st == sr_chr_1st.max()
                flag_max2 = sr_chr_2nd == sr_chr_2nd.max()

                idx1 = sr_chr_1st.loc[flag_max1]
                if len(idx1) > 0:
                    idx1 = idx1.index.to_series().str.extract(p_HLA2, expand=False).tolist()
                else:
                    idx1 = []

                idx2 = sr_chr_2nd.loc[flag_max2]
                if len(idx2) > 0:
                    idx2 = idx2.index.to_series().str.extract(p_HLA2, expand=False).tolist()
                else:
                    idx2 = []

                ## Manual Correction for final genotyping

                if len(idx1) == len(idx2):

                    if len(idx1) == 1:
                        # Normal Calling
                        # ex) [['0701'], ['1101']]
                        idx1 = idx1[0]
                        idx2 = idx2[0]

                    elif len(idx2) > 1:

                        if set(idx1) == set(idx2):
                            # Mirror Symmetry
                            # ex) [['1301', '1302'], ['1301', '1302']]
                            t = idx1
                            idx1 = t[0]
                            idx2 = t[1]
                        else:
                            # ex) [['1301', '1302'], ['1301', '0104']]
                            # ex) DPB1 : NA11992, ['0401,0402', '0202,0401']
                            # Imputation error -> Just export them as it is, but it will be classified as error in 'measureAccuracy.py.'
                            idx1 = ','.join(idx1)
                            idx2 = ','.join(idx2)

                else:
                    # Assuming len(idx1) > len(idx2)
                    if len(idx1) < len(idx2):
                        t = idx1
                        idx1 = idx2
                        idx2 = t

                    # ex) [['4006', '4403', '5101'], ['4403']]
                    # ex) [['0201', '0301', '2601'], ['2601']]
                    # ex) [['0201', '1001'], ['0201', '0901', '1001']]

                    idx1 = list(set(idx1).difference(set(idx2)))

                    idx1 = ','.join(idx1)
                    idx2 = ','.join(idx2)

                # print("PID: {}, {}".format(PID[j], [idx1, idx2]))
                f_out_rawcall.write("{} : {}\n".format(PID[j], [idx1, idx2]))

                dict_Call[PID[j]] = ','.join([idx1, idx2])

            _dict_HLA_Single_allele[_hla] = dict_Call

        f_out_rawcall.write('\n')




    # Two chromosome => One patient
    PID = [PID[i] for i in range(0, len(PID), 2)]

    with open(_out, 'w') as f_out:

        for pid in PID:
            for _hla in HLA_names:
                line = '\t'.join([pid, pid, _hla, ',', _dict_HLA_Single_allele[_hla][pid]]) + '\n'
                # print(line)
                f_out.write(line)



    return _out






if __name__ == '__main__':

    """
    
    _feature := 'DS', 'GP' or 'BOTH'
    
    (1) Dose(DS) Score based Voting method
        - Receipt.each_Imputation
        - Receipt.Vote_Score
        - Receipt.raw_call
        - HLA_IMPUTATION_OUT.single.alleles
    
    (2) Average of posterior probability(GP) based method.
        - Receipt.Avg_Posterior
        - Receipt.raw_call
        - HLA_IMPUTATION_OUT.single.alleles
    
    """

    dict_raw_IMP = {_exonN : {_overlap : None for _overlap in __overlap__} for _exonN in __EXON__}

    # in Ubuntu
    # dict_raw_IMP['exon2'][3000] = '/home/wanson/Git_Projects/CookHLA/tests/_3_CookHLA/20190711_BOTH/_3_HM_CEU_T1DGC_REF.MHC.exon2.3000.QC.doubled.imputation_out.vcf'
    # dict_raw_IMP['exon2'][4000] = '/home/wanson/Git_Projects/CookHLA/tests/_3_CookHLA/20190711_BOTH/_3_HM_CEU_T1DGC_REF.MHC.exon2.4000.QC.doubled.imputation_out.vcf'
    # dict_raw_IMP['exon2'][5000] = '/home/wanson/Git_Projects/CookHLA/tests/_3_CookHLA/20190711_BOTH/_3_HM_CEU_T1DGC_REF.MHC.exon2.5000.QC.doubled.imputation_out.vcf'
    # dict_raw_IMP['exon3'][3000] = '/home/wanson/Git_Projects/CookHLA/tests/_3_CookHLA/20190711_BOTH/_3_HM_CEU_T1DGC_REF.MHC.exon3.3000.QC.doubled.imputation_out.vcf'
    # dict_raw_IMP['exon3'][4000] = '/home/wanson/Git_Projects/CookHLA/tests/_3_CookHLA/20190711_BOTH/_3_HM_CEU_T1DGC_REF.MHC.exon3.4000.QC.doubled.imputation_out.vcf'
    # dict_raw_IMP['exon3'][5000] = '/home/wanson/Git_Projects/CookHLA/tests/_3_CookHLA/20190711_BOTH/_3_HM_CEU_T1DGC_REF.MHC.exon3.5000.QC.doubled.imputation_out.vcf'
    # dict_raw_IMP['exon4'][3000] = '/home/wanson/Git_Projects/CookHLA/tests/_3_CookHLA/20190711_BOTH/_3_HM_CEU_T1DGC_REF.MHC.exon4.3000.QC.doubled.imputation_out.vcf'
    # dict_raw_IMP['exon4'][4000] = '/home/wanson/Git_Projects/CookHLA/tests/_3_CookHLA/20190711_BOTH/_3_HM_CEU_T1DGC_REF.MHC.exon4.4000.QC.doubled.imputation_out.vcf'
    # dict_raw_IMP['exon4'][5000] = '/home/wanson/Git_Projects/CookHLA/tests/_3_CookHLA/20190711_BOTH/_3_HM_CEU_T1DGC_REF.MHC.exon4.5000.QC.doubled.imputation_out.vcf'


    # in OS X
    dict_raw_IMP['exon2'][3000] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190717_BOTH_ph/HM_CEU_T1DGC_REF_ph.MHC.exon2.3000.QC.doubled.imputation_out.vcf'
    dict_raw_IMP['exon2'][4000] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190717_BOTH_ph/HM_CEU_T1DGC_REF_ph.MHC.exon2.4000.QC.doubled.imputation_out.vcf'
    dict_raw_IMP['exon2'][5000] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190717_BOTH_ph/HM_CEU_T1DGC_REF_ph.MHC.exon2.5000.QC.doubled.imputation_out.vcf'
    dict_raw_IMP['exon3'][3000] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190717_BOTH_ph/HM_CEU_T1DGC_REF_ph.MHC.exon3.3000.QC.doubled.imputation_out.vcf'
    dict_raw_IMP['exon3'][4000] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190717_BOTH_ph/HM_CEU_T1DGC_REF_ph.MHC.exon3.4000.QC.doubled.imputation_out.vcf'
    dict_raw_IMP['exon3'][5000] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190717_BOTH_ph/HM_CEU_T1DGC_REF_ph.MHC.exon3.5000.QC.doubled.imputation_out.vcf'
    dict_raw_IMP['exon4'][3000] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190717_BOTH_ph/HM_CEU_T1DGC_REF_ph.MHC.exon4.3000.QC.doubled.imputation_out.vcf'
    dict_raw_IMP['exon4'][4000] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190717_BOTH_ph/HM_CEU_T1DGC_REF_ph.MHC.exon4.4000.QC.doubled.imputation_out.vcf'
    dict_raw_IMP['exon4'][5000] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190717_BOTH_ph/HM_CEU_T1DGC_REF_ph.MHC.exon4.5000.QC.doubled.imputation_out.vcf'

    # HLA_Genotype_Call(dict_raw_IMP, _feature='GP')
    # HLA_Genotype_Call(dict_raw_IMP, _feature='DS')
    HLA_Genotype_Call(dict_raw_IMP, _feature='BOTH')
