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



def HLA_Genotype_Call(_dict_IMP_Result, _feature='GP', __as_VCF=False, __forCheck=True, _out=None):

    if _feature != 'DS' and _feature != 'GP':
        print(std_ERROR_MAIN_PROCESS_NAME + "Wrong measure to calling genotype.")
        sys.exit()


    ###### < Variables > ######

    # Input
    dict_IMP_vcf = _dict_IMP_Result
    dict_IMP_vcf_LEFT = {_exonN: {_overlap: None for _overlap in __overlap__} for _exonN in __EXON__} # 0 ~ 8 columns (Marker information)
    dict_IMP_vcf_RIGHT = {_hla: [] for _hla in HLA_names} # 9 ~ columns (Patient genotype information.)

    # GP
    dict_IMP_vcf_RIGHT_mean = {_hla: None for _hla in HLA_names}

    # DS
    dict_IMP_vcf_RIGHT_DS = {_hla: {_exonN: {_overlap: None for _overlap in __overlap__} for _exonN in __EXON__} for _hla in HLA_names}
    dict_IMP_vcf_RIGHT_DS_Score = {_hla: None for _hla in HLA_names}

    n_row = 0
    n_col = 0

    PID = None

    # Called HLA (unphased) Genotype
    dict_HLA_Single_allele = {_hla: None for _hla in HLA_names}

    dict_forCheck_byHLA = {_hla: None for _hla in HLA_names}

    # Prefix
    OUTPUT_dir = os.path.dirname(dict_IMP_vcf['exon2'][3000])
    _out = _out if _out else join(OUTPUT_dir, 'IMPUTATION_OUT.{}.single.alleles'.format(_feature))



    ### Disintegrating each VCF file.
    for _exonN in __EXON__:
        for _overlap in __overlap__:

            print("({}, {})".format(_exonN, _overlap))

            with open(dict_IMP_vcf[_exonN][_overlap], 'r') as f_vcf, \
                 open(dict_IMP_vcf[_exonN][_overlap]+'.header', 'w') as f_vcf_header, \
                 open(dict_IMP_vcf[_exonN][_overlap]+'.body', 'w') as f_vcf_body:

                for line in f_vcf:

                    if line.startswith('##'):
                        f_vcf_header.write(line)
                    else:
                        f_vcf_body.write(line)

            df_temp = pd.read_csv(dict_IMP_vcf[_exonN][_overlap]+'.body', sep='\s+', header=0)

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



            ### Extrating 'GP(Genotpe Probability)' or 'Dose' values.

            needed_columns = [i for i in range(9, df_temp.shape[1])]

            if _feature == 'GP':
                df_onlyHLA = pd.concat([ID, df_onlyHLA.iloc[:, needed_columns]], axis=1).set_index('ID').astype(str).applymap(lambda x : x.split(':')[2].split(','))
                df_onlyHLA = df_onlyHLA.applymap(lambda x : list(map(float, x))).applymap(lambda x : x[0]+x[1]/2)

                df_onlyHLA = pd.concat([df_idx, df_onlyHLA], axis=1)
                # print("DataFrame of GP value :\n{}".format(df_onlyHLA.head()))

            elif _feature == 'DS':
                """
                Cf) VCF file works focused on ALT(Alternative allele). Meanwhile, The binary marker 'P(Present)' in CookHLA is located in REF(Reference).
                So, dose value of given vcf file is 'A(Absent)' marker's value and we will extract min value which corresponds to 'P(Present)'.
                """
                df_onlyHLA = pd.concat([ID, df_onlyHLA.iloc[:, needed_columns]], axis=1).set_index('ID').astype(str).applymap(lambda x : x.split(':')[1]).astype(float)
                df_onlyHLA = pd.concat([df_idx, df_onlyHLA], axis=1)
                # print("DataFrame of DS value :\n{}".format(df_onlyHLA.head()))

            # print(df_onlyHLA.shape)

            for _hla in HLA_names:

                if not isClassI[_hla] and _exonN == 'exon4':
                    continue

                df_temp = df_onlyHLA.filter(regex='^HLA_{}'.format(_hla), axis=0).reset_index('ID')

                if df_temp.shape[0] > 0:
                    dict_IMP_vcf_RIGHT[_hla].append(df_temp)






    if _feature == 'DS':

        dict_IMP_vcf_RIGHT_byHLA = {_hla: pd.concat(dict_IMP_vcf_RIGHT[_hla], axis=0) for _hla in HLA_names}

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
            df_acc.to_csv(join(OUTPUT_dir, 'DOSE_Score.{}.txt'.format(_hla)), sep='\t', header=True, index=True)




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
            df_forCheck.to_csv(join(OUTPUT_dir, 'each_Imputation.{}.{}.txt'.format(_hla, _feature)), sep='\t', header=True, index=True)



        ### Call HLA genotype.

        for _hla in HLA_names:

            # print("{} :".format(_hla))

            curr_DOSE_Score = dict_IMP_vcf_RIGHT_DS_Score[_hla]
            # print(curr_DOSE_Score.head())

            PID = curr_DOSE_Score.columns.tolist()
            dict_Call = {}


            for j in range(0, curr_DOSE_Score.shape[1], 2):

                sr_chr_1st = curr_DOSE_Score.iloc[:, j]
                sr_chr_2nd = curr_DOSE_Score.iloc[:, j+1]

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

                dict_Call[PID[j]] = ','.join([idx1, idx2])

            dict_HLA_Single_allele[_hla] = dict_Call

        # Two chromosome => One patient
        PID = [PID[i] for i in range(0, len(PID), 2)]



        with open(_out, 'w') as f_out:

            for pid in PID:
                for _hla in HLA_names:

                    line = '\t'.join([pid, pid, _hla, ',', dict_HLA_Single_allele[_hla][pid]]) + '\n'
                    # print(line)
                    f_out.write(line)








        # if __as_VCF:
        #
        #     ### Output as VCF file. (Under construction)
        #
        #     for _hla in HLA_names:
        #
        #         ## mean DS
        #
        #         df_temp = dict_IMP_vcf_RIGHT_mean[_hla].transpose()
        #         # print("df_temp : \n{}".format(df_temp.head()))
        #         # print(df_temp.column)
        #
        #
        #         dict_temp = {
        #             '#CHROM': ['6' for z in range(0, df_temp.shape[0])],
        #             'POS': [HLA_mid_position[_hla] for z in range(0, df_temp.shape[0])],
        #             'ID': df_temp.index,
        #             'REF': ['G' for z in range(0, df_temp.shape[0])],
        #             'ALT': ['C' for z in range(0, df_temp.shape[0])],
        #             'QUAL': ['.' for z in range(0, df_temp.shape[0])],
        #             'FILTER': ['PASS' for z in range(0, df_temp.shape[0])],
        #             'INFO': ['.' for z in range(0, df_temp.shape[0])],
        #             'FORMAT': ['GT' for z in range(0, df_temp.shape[0])]
        #         }
        #
        #         df_temp_LEFT = pd.DataFrame.from_dict(dict_temp)
        #         pd.concat([df_temp_LEFT, df_temp.reset_index(drop=True)], axis=1).to_csv(join(OUTPUT_dir, 'EXON_VCF_HLA_{}.mean_DS.txt'.format(_hla)), sep='\t', header=True, index=False)
        #
        #
        #         ## GT of mean DS
        #         df_temp2 = dict_IMP_vcf_RIGHT_DS_Score[_hla]
        #         # print(df_temp2.head())
        #
        #         dict_temp = {
        #             '#CHROM': ['6' for z in range(0, df_temp.shape[0])],
        #             'POS': [HLA_mid_position[_hla] for z in range(0, df_temp.shape[0])],
        #             'ID': df_temp2.index,
        #             'REF': ['G' for z in range(0, df_temp.shape[0])],
        #             'ALT': ['C' for z in range(0, df_temp.shape[0])],
        #             'QUAL': ['.' for z in range(0, df_temp.shape[0])],
        #             'FILTER': ['PASS' for z in range(0, df_temp.shape[0])],
        #             'INFO': ['.' for z in range(0, df_temp.shape[0])],
        #             'FORMAT': ['GT' for z in range(0, df_temp.shape[0])]
        #         }
        #
        #         df_temp_LEFT = pd.DataFrame.from_dict(dict_temp)
        #
        #         pd.concat([df_temp_LEFT, df_temp2.reset_index(drop=True)], axis=1).to_csv(join(OUTPUT_dir, 'EXON_VCF_HLA_{}.mean_DS_GT.txt'.format(_hla)), sep='\t', header=True, index=False)







    elif _feature == 'GP':

        ### Acquring mean posterior GP

        dict_IMP_vcf_RIGHT_byHLA = {_hla: pd.concat(dict_IMP_vcf_RIGHT[_hla], axis=0) for _hla in HLA_names}

        for _hla in HLA_names:

            l_HLA_mean = []

            for idx, df in dict_IMP_vcf_RIGHT_byHLA[_hla].groupby('ID'):

                sr_temp = df.iloc[:, 3:].mean()
                sr_temp.name = idx

                l_HLA_mean.append(sr_temp)


            df_HLA_mean = pd.concat(l_HLA_mean, axis=1).transpose()
            df_HLA_mean.index.name = 'ID'
            print("df_mean : \n{}".format(df_HLA_mean.head()))

            dict_IMP_vcf_RIGHT_mean[_hla] = df_HLA_mean
            df_HLA_mean.to_csv(join(OUTPUT_dir, 'Mean_Posterior.{}.txt'.format(_hla)), sep='\t', header=True, index=True)


        if __as_VCF:

            # new Header

            representative_header = join(OUTPUT_dir, 'representative_header.txt')
            representative_body = join(OUTPUT_dir, 'representative_body.txt')

            with open(dict_IMP_vcf['exon2'][3000]+'.header', 'r') as f_header, open(representative_header, 'w') as f_header_out:

                l_to_remove = ['ID={}'.format(item) for item in ['AF', 'AR2', 'DR2', 'IMP', 'GT', 'DS']]

                for line in f_header:

                    flag_ommit = False

                    for to_remove in l_to_remove:
                        if re.search(to_remove, line):
                            flag_ommit = True
                            break

                    if not flag_ommit:
                        f_header_out.write(line)


            l_temp = []

            for _hla in HLA_names:

                df_temp = dict_IMP_vcf_RIGHT_mean[_hla]

                dict_temp = {
                    '#CHROM': ['6' for z in range(0, df_temp.shape[0])],
                    'POS': [HLA_mid_position[_hla] for z in range(0, df_temp.shape[0])],
                    'ID': df_temp.index,
                    'REF': ['G' for z in range(0, df_temp.shape[0])],
                    'ALT': ['C' for z in range(0, df_temp.shape[0])],
                    'QUAL': ['.' for z in range(0, df_temp.shape[0])],
                    'FILTER': ['PASS' for z in range(0, df_temp.shape[0])],
                    'INFO': ['.' for z in range(0, df_temp.shape[0])],
                    'FORMAT': ['GP' for z in range(0, df_temp.shape[0])]
                }

                df_temp_LEFT = pd.DataFrame.from_dict(dict_temp)
                # print(df_temp_LEFT.head())

                df_temp = pd.concat([df_temp_LEFT, df_temp.reset_index(drop=True)], axis=1)
                # print(df_temp.head())

                l_temp.append(df_temp)


            df_new_VCF_body = pd.concat(l_temp, axis=0)
            df_new_VCF_body.to_csv(representative_body, sep='\t', header=True, index=False)

            RUN_Bash('cat {} {} > {}'.format(representative_header, representative_body, join(OUTPUT_dir, 'posterior_probability.vcf')))
            RUN_Bash('rm {}'.format(representative_header))
            RUN_Bash('rm {}'.format(representative_body))


        ### Generating '*.alleles` file.

        # Patient_ID = dict_IMP_vcf_RIGHT_mean['A'].columns.tolist() # Using just first one.
        # # print(Patient_ID)
        #
        # for _hla in HLA_names:
        #
        #     # print(dict_IMP_vcf_RIGHT_mean[_hla].head())
        #
        #     HLA_alleles = dict_IMP_vcf_RIGHT_mean[_hla].index.tolist()
        #     max_probs = dict_IMP_vcf_RIGHT_mean[_hla].apply(max, axis=0).tolist()
        #
        #
        #
        #     dict_HLA_allele_called = {pid: None for pid in Patient_ID}
        #
        #     for j in range(0, dict_IMP_vcf_RIGHT_mean[_hla].shape[1]):
        #
        #         among_probs = dict_IMP_vcf_RIGHT_mean[_hla].iloc[:, j].tolist()
        #
        #         idx_maxes = [i for i in range(0, len(among_probs)) if among_probs[i] == max_probs[j]]
        #
        #         if len(idx_maxes) > 0:
        #
        #             l_temp = []
        #
        #             for idx in idx_maxes:
        #
        #                 m = p_HLA2.match(HLA_alleles[idx])
        #
        #                 if m:
        #                     l_temp.append(m.group(1))
        #
        #             dict_HLA_allele_called[Patient_ID[j]] = ','.join(l_temp)
        #         else:
        #             dict_HLA_allele_called[Patient_ID[j]] = ""
        #
        #     dict_HLA_Single_allele[_hla] = pd.DataFrame(dict_HLA_allele_called, index=[_hla])
        #
        #
        # # for _hla in HLA_names:
        # #     print("{}   :\n{}".format(_hla, dict_HLA_Single_allele[_hla]))
        #
        # ### Making '*.single.alleles' file.
        #
        # df_Single_alleles = pd.concat([dict_HLA_Single_allele[_hla] for _hla in HLA_names], axis=0)
        # print(df_Single_alleles)
        # df_Single_alleles.to_csv('single_allele.txt', sep='\t', header=True, index=False)
        #
        #
        #
        # n_col_single_allele = int(df_Single_alleles.shape[1]/2)
        #
        # l_df_one_patient = []
        #
        # for i in range(0, df_Single_alleles.shape[1], 2):
        #
        #     df_one_patient = df_Single_alleles.iloc[:, [i, i+1]].apply(lambda x : ','.join(x), axis=1)
        #
        #     print(df_one_patient)




        ### Call HLA genotype.

        for _hla in HLA_names:

            # print("{} :".format(_hla))

            curr_Posterior_avg = dict_IMP_vcf_RIGHT_mean[_hla]
            # print(curr_Posterior_avg.head())

            PID = curr_Posterior_avg.columns.tolist()
            dict_Call = {}


            for j in range(0, curr_Posterior_avg.shape[1], 2):

                sr_chr_1st = curr_Posterior_avg.iloc[:, j]
                sr_chr_2nd = curr_Posterior_avg.iloc[:, j+1]

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

                dict_Call[PID[j]] = ','.join([idx1, idx2])

            dict_HLA_Single_allele[_hla] = dict_Call

        # Two chromosome => One patient
        PID = [PID[i] for i in range(0, len(PID), 2)]



        with open(_out, 'w') as f_out:

            for pid in PID:
                for _hla in HLA_names:

                    line = '\t'.join([pid, pid, _hla, ',', dict_HLA_Single_allele[_hla][pid]]) + '\n'
                    # print(line)
                    f_out.write(line)









if __name__ == '__main__':

    """
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
    dict_raw_IMP['exon2'][3000] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_HM_CEU_T1DGC_REF.MHC.exon2.3000.QC.doubled.imputation_out.vcf'
    dict_raw_IMP['exon2'][4000] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_HM_CEU_T1DGC_REF.MHC.exon2.4000.QC.doubled.imputation_out.vcf'
    dict_raw_IMP['exon2'][5000] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_HM_CEU_T1DGC_REF.MHC.exon2.5000.QC.doubled.imputation_out.vcf'
    dict_raw_IMP['exon3'][3000] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_HM_CEU_T1DGC_REF.MHC.exon3.3000.QC.doubled.imputation_out.vcf'
    dict_raw_IMP['exon3'][4000] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_HM_CEU_T1DGC_REF.MHC.exon3.4000.QC.doubled.imputation_out.vcf'
    dict_raw_IMP['exon3'][5000] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_HM_CEU_T1DGC_REF.MHC.exon3.5000.QC.doubled.imputation_out.vcf'
    dict_raw_IMP['exon4'][3000] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_HM_CEU_T1DGC_REF.MHC.exon4.3000.QC.doubled.imputation_out.vcf'
    dict_raw_IMP['exon4'][4000] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_HM_CEU_T1DGC_REF.MHC.exon4.4000.QC.doubled.imputation_out.vcf'
    dict_raw_IMP['exon4'][5000] = '/Users/wansun/Git_Projects/CookHLA/tests/_3_HM_CEU_T1DGC_REF.MHC.exon4.5000.QC.doubled.imputation_out.vcf'

    HLA_Genotype_Call(dict_raw_IMP, _feature='GP')
    # AverageImpResults(dict_raw_IMP, _feature='DS', __as_VCF=True)








    # ### Acquiring mean
    #
    # for _hla in HLA_names:
    #
    #     if len(dict_IMP_vcf_RIGHT[_hla]) > 0:
    #
    #         df_acc = dict_IMP_vcf_RIGHT[_hla][0]
    #
    #         # accumulation of the values(GP or DS) of each result
    #         for i in range(1, len(dict_IMP_vcf_RIGHT[_hla])):
    #             df_acc = df_acc + dict_IMP_vcf_RIGHT[_hla][i]
    #
    #         dict_IMP_vcf_RIGHT_mean[_hla] = df_acc / (9 if isClassI[_hla] else 6) # dividing 9 if HLA is Class I else 6.
    #
    #         # print("\nHLA_{} mean:\n{}".format(_hla, dict_IMP_vcf_RIGHT_mean[_hla].head()))
    #
    #
    #
