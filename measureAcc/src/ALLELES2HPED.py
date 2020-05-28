#-*- coding: utf-8 -*-
import os, sys, re
import pandas as pd


HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
HLA_names2 = ["B", "A", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]

p_4digit = re.compile(r'\d{4,6}')
p_4digit_pair = re.compile(r'^(\d{4,5})?,(\d{4,5})?$')
p_2digit_pair = re.compile(r'\d{2}?,\d{2}?')



def ALLELES2HPED(_alleles, _out=None, _f_HLA_DRB1_1454to1401=False):
    
    df_alleles = pd.read_csv(_alleles, sep='\s+', header=None, dtype=str)
#     print("df_alleles :\n{}\n".format(df_alleles.head()))

    ### Finding out 4-digit allele column
    f_4digit = df_alleles.applymap(lambda x : (bool(p_4digit_pair.match(x)) and (not bool(re.match(r'^,$', x))))).apply(lambda x : x.any(), axis=0)
#     print(f_4digit)
    
    df_alleles = pd.concat([df_alleles.iloc[:, [0,1,2]], df_alleles.loc[:, f_4digit]], axis=1)
    df_alleles.columns = ['FID', 'IID', 'HLA', 'Allele']

    
    ### Bizarre alleles (especially for SNP2HLA)
#     print("df_alleles:\n{}\n".format(df_alleles))
    
    f_toExclude = ~df_alleles['Allele'].str.match(p_4digit_pair)
    
    if f_toExclude.any():
        df_toExclude = df_alleles[f_toExclude]
        print("[Heads-up] Next bizarre alleles will be set to '0,0' by force.\n{}\n".format(df_toExclude))
    
        if bool(_out):
            df_toExclude.to_csv(_out+'.Excluded.alleles', sep='\t', header=True, index=False)
            
        df_alleles.loc[f_toExclude, 'Allele'] = '0,0'
        print("Revised:\n{}\n".format(df_alleles[f_toExclude]))
    

    
    df_alleles = df_alleles.set_index(['FID', 'IID', 'HLA']).unstack('HLA')
    df_alleles.columns = df_alleles.columns.droplevel(level=0)
    # print("df_alleles:\n{}\n".format(df_alleles))
    
    ### processing not give HLA genes.
    given_HLAs = df_alleles.columns.tolist()
    given_index = df_alleles.index

#     print(given_HLAs)

    df_alleles = pd.concat([
        df_alleles[HLA_names[i]] if (HLA_names[i] in given_HLAs)  \
        else pd.Series(['0,0' for z in range(df_alleles.shape[0])], name=HLA_names[i], index=given_index) 
        for i in range(len(HLA_names))
    ], axis=1)

#     l_temp = []
    
#     for i in range(len(HLA_names)):
        
#         if HLA_names[i] in given_HLAs:
#             print(df_alleles[HLA_names[i]])
#             l_temp.append(df_alleles[HLA_names[i]])
#         else:
#             sr_temp = pd.Series(['0,0' for z in range(df_alleles.shape[0])], name=HLA_names[i], index=given_index)

#             l_temp.append(sr_temp)
            
#     df_alleles = pd.concat(l_temp, axis=1)
            
    
#     print("df_alleles :\n{}\n".format(df_alleles))
    
    
    l_temp = []

    for _c, _sr in df_alleles.iteritems():
        
        column_name = _c
        sr_column = _sr
        
        l_temp.append(sr_column.str.split(',', expand=True).applymap(lambda x : x if bool(p_4digit.match(x)) else '0'))
        
    df_alleles2 = pd.concat(l_temp, axis=1)
#     print("df_alleles2 :\n{}\n".format(df_alleles2.head()))


    if _f_HLA_DRB1_1454to1401:
        df_alleles2 = HLA_DRB1_1454to1401(df_alleles2, _alleles)
    
    
    df_Idx = df_alleles2.index.to_frame()
    df_Idx['PID'] = ['0' for z in range(df_Idx.shape[0])]
    df_Idx['MID'] = ['0' for z in range(df_Idx.shape[0])]
    df_Idx['Sex'] = ['0' for z in range(df_Idx.shape[0])]
    df_Idx['Phe'] = ['-9' for z in range(df_Idx.shape[0])]

#     print("df_Idx :\n{}\n".format(df_Idx.head()))
    
    
    df_alleles2 = pd.concat([df_Idx, df_alleles2], axis=1)
#     print("df_alleles2 :\n{}\n".format(df_alleles2.head()))

    
    
    if bool(_out):
        df_alleles2.to_csv(_out+'.hped', sep='\t', header=False, index=False)
        return _out + '.hped'
    else:
        return df_alleles2




def HLA_DRB1_1454to1401(_hped_right, _alleles):
    f_1451 = (_hped_right.iloc[:, [14, 15]] == '1454').apply(lambda x: x.any(), axis=1)

    if f_1451.any():
        print("HLA_DRB1*1454 will be considered as 1401. ('{}')".format(_alleles))

        df_HLA_DRB1 = _hped_right.iloc[:, [14, 15]].replace('1454', '1401')
        return pd.concat([_hped_right.iloc[:, :14], df_HLA_DRB1], axis=1)

    else:
        return _hped_right



    
if __name__ == '__main__':
    
    """
    ALLELES2HPED.py
    
    """

    pass
    
    
    # [_answer, _out] = sys.argv[1:]
    #
    # ALLELES2HPED(_answer, _out)
