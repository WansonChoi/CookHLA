import os, sys, re
import pandas as pd

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
Meta_Info = ['FID', 'IID', 'PID', 'MID', 'Sex', 'Phe']

p_4digit = re.compile(r'^\d{4,5}$')
p_Nfield = re.compile(r'^\w+\*\d{2,3}(:\d{2,3})+$')


def SieveCHPED(_hped, _chped, _out=None):
    
    ### HPED file
    __HPED__ = pd.read_csv(_hped, sep='\s+', header=None, dtype=str, index_col=[0,1,2,3,4,5])
#     print("__HPED__:\n{}\n".format(__HPED__.iloc[:5, :6]))
    
    
    
    ### CHPED file
    __CHPED__ = pd.read_csv(_chped, sep='\s+', header=None, dtype=str, index_col=[0,1,2,3,4,5])
#     print("__CHPED__:\n{}\n".format(__CHPED__.iloc[:5, :6]))
    
    
    
    if __HPED__.shape[0] != __CHPED__.shape[0]:
        print("# of rows of given HPED and CHPED file are diffenet.")
        return -1
    
    
    
    
    ##### < Main iteration > #####
    
    l_new_lines = []
    l_deprecated_ones = []
        
    gentr_hped = __HPED__.itertuples()
    gentr_chped = __CHPED__.itertuples()
    
    count = 0
    
    for i in range(__HPED__.shape[0]):
        
        line_hped = next(gentr_hped)
        line_chped = next(gentr_chped)
        
        l_temp = []
        
        for i in range(len(HLA_names)):

            idx_1 = (2*i + 1)
            idx_2 = idx_1 + 1

            f_al1 = False
            f_al2 = False
            
#             print("current HPED allele: {} / {}".format(line_hped[idx_1], line_hped[idx_2]))
#             print("current CHPED allele: {} / {}".format(line_chped[idx_1], line_chped[idx_2]))
            
            
            ##### Removed cases
            
            if line_hped[idx_1] != '0':

                f_al1 = bool(p_4digit.match(line_hped[idx_1])) and (line_chped[idx_1] == '0')
                
                if f_al1:
                    l_temp.append('deprecated')
                else:
                    l_temp.append(line_chped[idx_1])
                    
            else:
                l_temp.append('0')

            
            
            if line_hped[idx_2] != '0':
                
                f_al2 = bool(p_4digit.match(line_hped[idx_2])) and (line_chped[idx_2] == '0')
                
                if f_al2:
                    l_temp.append('deprecated')
                else:
                    l_temp.append(line_chped[idx_2])            
                    
            else:
                l_temp.append('0')
                        
                
                
            
            if f_al1 or f_al2:
                l_deprecated_ones.append(count)
                
            
        l_new_lines.append(l_temp)
            
#         print("\n============================\n")
        
        count += 1
#         if count > 5: break
            
            
    df_Marked_chped = pd.DataFrame(l_new_lines)
    df_Marked_chped.index = __CHPED__.index
#     print("df_Marked_chped:\n{}\n".format(df_Marked_chped))
    
    
    if len(l_deprecated_ones) > 0:
        
        if len(l_deprecated_ones) < 10:
            print("deprecated samples:\n{}\n".format(l_deprecated_ones))
        else:
            print("There are {} deprecated samples.".format(len(l_deprecated_ones)))
            
            
        if bool(_out):
            
#             pd.DataFrame(l_deprecated_ones).to_csv(_out+'.Marked.deprecated.txt', sep='\t', header=False, index=False)
            df_Marked_chped.iloc[l_deprecated_ones, :].to_csv(_out+'.Marked.deprecated.txt', sep='\t', header=False, index=True)
    else:
        print("There is no removed/deprecated alleles.")

                
    
    if bool(_out):
        df_Marked_chped.to_csv(_out+'.Marked.chped', sep='\t', header=False, index=True)
        return _out+'.Marked.chped'
    else:
        return df_Marked_chped
    
    

    
    
    
if __name__ == '__main__':
    
    """
    SieveCHPED.py
    
    """

    pass

    # [_hped, _chped, _out] = sys.argv[1:]
    #
    # SieveCHPED(_hped, _chped, _out)
