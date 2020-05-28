#-*- coding: utf-8 -*-
import os, sys, re
import pandas as pd

from os.path import exists, join, dirname, basename

from measureAcc.NomenCleaner.NomenCleaner import HATK_NomenCleaner
from measureAcc.src.SieveCHPED import SieveCHPED
from measureAcc.src.ALLELES2HPED import ALLELES2HPED

from src.CookHLAError import CookHLAInputPreparationError

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
Meta_Info = ['FID', 'IID', 'PID', 'MID', 'Sex', 'Phe']

p_HLA_allele = re.compile(r'\w+\*\d{2,3}(:\d{2,3})*')




class CookHLA_measureAcc(object):

    def __init__(self, _answer, _imputed, _out, _allele_group='measureAcc/data/SameExon234.MERGED.nuc.txt'):

        """
        1. imputed *.alleles
        2. answer *.alleles

        (1) ALLELES2HPED.py

        (2) DRB1 1454 -> 1401

        (3) NomenCleaner

        (4) SieveCHEPD.py

        (5) measureAccuracy_v3.5

        """

        ## Exception Handling here.

        # (1) answer file
        if not exists(_answer):
            raise CookHLAInputPreparationError("Given answer file can't be found.('{}')".format(_answer))

        if not (_answer.endswith('.alleles') or _answer.endswith('.hped') or _answer.endswith('.Marked.chped')):
            # No *.chped alone
            raise CookHLAInputPreparationError("Given answer file must have file extension "
                                               "either '*.alleles', '*.hped' or '*.Marked.chped'.")

        # (2) imputed file
        if not exists(_imputed):
            raise CookHLAInputPreparationError("Given imputed file can't be found.('{}')".format(_imputed))

        if not (_imputed.endswith('.alleles') or _imputed.endswith('.hped') or _imputed.endswith('.Marked.chped')):
            # No *.chped alone
            raise CookHLAInputPreparationError("Given imputed file must have file extension "
                                               "either '*.alleles', '*.hped' or '*.Marked.chped'.")


        # (3) _allele_group file
        if not exists(_allele_group):
            raise CookHLAInputPreparationError("Given Allele Group file can't be found.('{}')".format(_allele_group))



        ## measureAccuracy_v3.5
        answer_Marked_chped = self.ConvertToMarkedCHPED(_answer, dirname(_out))
        imputed_Marked_chped = self.ConvertToMarkedCHPED(_imputed, dirname(_out))

        self.accuracy = measureAccuracy(answer_Marked_chped, imputed_Marked_chped, _out, _allele_group)




    def ConvertToMarkedCHPED(self, _f, _out_dir):

        if _f.endswith('.alleles'):
            # ALLELES2HPED
            t_out = join(_out_dir, re.sub(r'\.alleles$', '', basename(_f)))
            _f = ALLELES2HPED(_f, t_out, _f_HLA_DRB1_1454to1401=True)
            # print("alleles -> hped: {}".format(_f))

        if _f.endswith('.hped'):
            # (1) HLA_DRB1 1454 to 1401

            # (2) NomenCleaner
            t1_out = re.sub(r'hped$', 'imgt3320.4field', _f)
            t = HATK_NomenCleaner(_f, "measureAcc/NomenCleaner/HLA_ALLELE_TABLE.imgt3320.hat", '3320', t1_out,
                      __f_NoCaption=False, __leave_NotFound=False,
                      __oneF=False, __twoF=False, __threeF=False, __fourF=True, __Ggroup=False, __Pgroup=False)
            # print("hped -> chped: {}".format(t.chped))

            # (3) SieveCHPED
            _f = SieveCHPED(_f, t.chped, t1_out)
            # print("chped -> Marked.chped: {}".format(_f))

        if _f.endswith('.Marked.chped'):
            return _f
        else:
            return '-1'




def whichGroup(_allele, _allele_gruop):
    
    if (_allele == '0') or (_allele =='deprecated'):
        return _allele
    else:
    
        f_allele = _allele_gruop.map(lambda x : bool(re.search(pattern=re.escape(_allele), string=x)))

        if f_allele.any():
#             print(_allele_gruop[f_allele])
            return _allele_gruop[f_allele].iat[0]
        else:
            return '0'


        
def measureAccuracy(_answer, _imputed, _out=None, _allele_group=None):
    
    """
    version : v3.5


    By taking the answer and imputed files which went through NomenCleaner('*.chped'),
    calculate the accuracy of the imputed file.
    
    If allele group information is given, then some alleles that belong to the same group
    will be considered as the same ones and accuracy will be calculated.
    
    """
    
    df_answer = pd.read_csv(_answer, sep='\s+', header=None, dtype=str, names=Meta_Info+[hla+i for hla in HLA_names for i in ['_1', '_2']]) \
                    .drop(['PID', 'MID', 'Sex', 'Phe'], axis=1)
#     print("df_answer :\n{}\n".format(df_answer.head()))
    
    df_imputed = pd.read_csv(_imputed, sep='\s+', header=None, dtype=str, names=Meta_Info+[hla+i for hla in HLA_names for i in ['_1', '_2']]) \
                    .drop(['PID', 'MID', 'Sex', 'Phe'], axis=1)
#     print("df_imputed :\n{}\n".format(df_imputed.head()))


    if bool(_allele_group):
        
        df_allele_group = pd.read_csv(_allele_group, sep='\s+', header=0, dtype=str, usecols=[0])
#         print("df_allele_group :\n{}\n".format(df_allele_group))
        
        dict_allele_group = {HLA_names[i] : df_allele_group['Alleles'][df_allele_group['Alleles'].str.match(HLA_names[i])] \
                                 for i in range(len(HLA_names))}
        
#         for k, v in dict_allele_group.items():
#             print("{} :\n{}".format(k, v.head()))
              

    
    
    ### Splitting the columns of the given chped file based on HLA.
    
    dict_Table = {HLA_names[i] : None for i in range(len(HLA_names))}
    dict_acc = {HLA_names[i] : None for i in range(len(HLA_names))}
    
    
    
    ### Main iteration over each HLA.
    
    l_eachRow_log = []
    
    for i in range(len(HLA_names)):
#     for i in [0,3,4]:
        
        df_Table_byHLA = \
            df_imputed.filter(items=['FID', 'IID', '{}_1'.format(HLA_names[i]), '{}_2'.format(HLA_names[i])], axis=1).merge(
            df_answer.filter(items=['FID', 'IID', '{}_1'.format(HLA_names[i]), '{}_2'.format(HLA_names[i])], axis=1),
            left_on=['FID', 'IID'], right_on=['FID', 'IID'], how='left', suffixes=('_imputed', '_answer')
        ).fillna('0')
        
#         print("df_Table_byHLA(Before):\n{}".format(df_Table_byHLA))
#         print("dict_allele_group:\n{}".format(dict_allele_group[HLA_names[i]]))

    

        
        ### Filtering out rows where answer values are valid.
        f_answer1 = df_Table_byHLA['{}_1_answer'.format(HLA_names[i])].map(lambda x : bool(p_HLA_allele.match(x)))
        f_answer2 = df_Table_byHLA['{}_2_answer'.format(HLA_names[i])].map(lambda x : bool(p_HLA_allele.match(x)))
        
        f_BothValid = f_answer1 & f_answer2
        l_BothValid = f_BothValid.tolist()
#         print(f_BothValid)
        
        
        if f_BothValid.any():
            
            
            ### Applying allele group information

            if bool(_allele_group):

                df_temp = df_Table_byHLA[['{}_1_imputed'.format(HLA_names[i]), '{}_2_imputed'.format(HLA_names[i]), 
                                          '{}_1_answer'.format(HLA_names[i]), '{}_2_answer'.format(HLA_names[i])]] \
                            .applymap(lambda x : whichGroup(x, dict_allele_group[HLA_names[i]]))

                df_Table_byHLA = pd.concat([df_Table_byHLA[['FID', 'IID']], df_temp], axis=1)

#                 print("df_Table_byHLA(After):\n{}".format(df_Table_byHLA))
            
            
            
            count = 0
            
            correct = 0
            total = 0
            
            for idx in f_BothValid.index.tolist():
                
                t_correct = 0

                [FID, IID, imputed_1, imputed_2, answer_1, answer_2] \
                    = df_Table_byHLA.iloc[idx, :].tolist()

                set_imputed = set([imputed_1, imputed_2])
                set_answer = set([answer_1, answer_2])                
                
#                 print("set_imputed/answer : {} / {}".format(set_imputed, set_answer))


                
                ### Main classification
                
                if l_BothValid[idx]:
                    
                    
                    if 'deprecated' in set_imputed:
                        
                        set_imputed = set.difference(set_imputed, {'deprecated'})
                        
                        if len(set_imputed) == 0:
                            # Both imputed alleles are deprecated.
                            # No way to save this case.
                            # NO TOUCH.
                            
                            l_eachRow_log.append(df_Table_byHLA.iloc[idx, :].tolist() + ['NA'])
                            
                        elif len(set_imputed) == 1:

                            if not set.isdisjoint(set_imputed, set_answer):
                                t_correct += 1
                            else:
                                # t_correct += 0
                                pass

                            correct += t_correct
                            total += 1 # (***)
                            
                            l_eachRow_log.append([FID, IID, imputed_1, imputed_2, answer_1, answer_2, '{}*'.format(t_correct)])
                        
                    else:

                        if set.isdisjoint(set_imputed, set_answer):

                            # Nothing matched (Both alleles have been wrongly imputed.)
                            t_correct = 0

                        else:
                            # There is an element in intersection.

                            if set_imputed == set_answer:
                                # 2 alleles are both correctly imputed.
                                t_correct = 2

                            else:
                                t_correct = 1                                                        
                            
                            
                        correct += t_correct
                        total += 2

                        l_eachRow_log.append([FID, IID, imputed_1, imputed_2, answer_1, answer_2, '{}'.format(t_correct)])

                    count += 1
    #                 if count > 5 : break

    
                else:
                    # Answer sample doesn't have perfect two alleles.
                    # => Won't be considered in getting accuracy.
                    
                    if 'deprecated' in set_answer:
                        
                        set_answer = set.difference(set_answer, {'deprecated'})
                        
                        if len(set_answer) == 0:
                            # Both answer alleles are deprecated.
                            # No way to save this case.
                            # NO TOUCH.
                            
                            l_eachRow_log.append(df_Table_byHLA.iloc[idx, :].tolist() + ['NA'])
                            
                        elif len(set_answer) == 1:
                            
                            if not set.isdisjoint(set_answer, set_imputed):
                                # There is an element in common => one allele to be considered as answer.
                                t_correct = 1
                            else:
                                # t_correct += 0
                                pass
                            
                            
                            correct += t_correct
                            total += 1
                            
                            
                            l_eachRow_log.append([FID, IID, imputed_1, imputed_2, answer_1, answer_2, '{}*'.format(t_correct)])

                    else:
                        # NO TOUCH.
                        l_eachRow_log.append(df_Table_byHLA.iloc[idx, :].tolist() + ['NA'])
    
    
            acc = float(correct)/total
            dict_acc[HLA_names[i]] = acc

            print("correct : {} / total : {}".format(correct, total))
            print("HLA-{} acc : {}".format(HLA_names[i], acc))

            
           
        else:
            # Nothing to check.
            acc = -1
            dict_acc[HLA_names[i]] = acc
            
            
            
    sr_acc = pd.Series(dict_acc) # Final accuracy output
            
    if bool(_out):
            
        df_log_eachRow = pd.DataFrame(l_eachRow_log, columns=['FID', 'IID', 'imputed_1', 'imputed_2', 'answer_1', 'answer_2', 'correct'])
#         print("df_log:\n{}\n".format(df_log_eachRow))
        df_log_eachRow.to_csv(_out+'.accuracy.log', sep='\t', header=True, index=False)
        sr_acc.to_csv(_out+'.accuracy', sep='\t', header=False, index=True)
              
        return _out+'.accuracy'

    else:
        return sr_acc
    

    
if __name__ == '__main__':

    """
    measureAccuracy_v3.5.py
    """
    pass