#-*- coding: utf-8 -*-
import sys, os, re
import pandas as pd
from statistics import mean


p_1stTwo = re.compile(r'^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)')

p_AA = re.compile(r'^AA_')
p_SNP = re.compile(r'^SNP_')
p_HLA = re.compile(r'^HLA_')
p_INS = re.compile(r'^INS_')



def Make_EXON234_AGM(_GM, _REF_exon234_markers, _out):

    df_GM = pd.read_csv(_GM, sep='\s+', header=None, names=['Chr', 'ID', 'GD', 'BP'], dtype=str)
    # print(df_GM.head())

    df_REF_markers = pd.read_csv(_REF_exon234_markers, sep='\s+', header=None, names=['ID', 'BP', 'A1', 'A2'], usecols=[0, 1], dtype=str)
    # print(df_REF_markers.head())


    ### Subsetting only SNP markers and HLA allele markers in Adaptive Genetic Map.

    ID = df_GM.iloc[:, 1]

    flag_AA = ID.str.match(p_AA)
    flag_SNP = ID.str.match(p_SNP)
    flag_HLA = ID.str.match(p_HLA)
    flag_INS = ID.str.match(p_INS)

    flag_toexclude = flag_AA | flag_SNP | flag_HLA | flag_INS

    df_GM_filtered = df_GM.loc[~flag_toexclude, :]
    # print(df_GM_filtered.head())
    # df_GM_filtered.to_csv(_out+'.onlyRS.txt', sep='\t', header=False, index=False)




    ### Merging Genetic Map and Ref_markers file.

    df_merged = df_GM_filtered.merge(df_REF_markers, left_on='ID', right_on='ID', how='outer').fillna('0')
    # print(df_merged.head(50))
    # df_merged.to_csv(_out+'.onlyRS.merged.txt', sep='\t', header=True, index=True)

    df_merged2 = pd.DataFrame(GEN_Trim_Columns(df_merged), columns=['Chr', 'ID', 'GD', 'BP']).sort_values('BP')
    # print(df_merged2.head())
    # df_merged2.to_csv(_out+'.onlyRS.merged2.txt', sep='\t', header=True, index=True)


    __RETURN__ = GEN_stitch_GD(df_merged2)
    # print(__RETURN__.head())
    __RETURN__.to_csv(_out, sep='\t', header=False, index=False)


    return _out



def GEN_Trim_Columns(_df_merged):

    for row in _df_merged.itertuples():

        _ID_ = row[2]
        _BP_ = row[5] if row[5] != '0' else row[4]

        if not p_HLA.match(_ID_):
            yield [row[1], _ID_, row[3], _BP_] # ex) '6	rs9277935	9.094223819100000	33268403'
        else:
            yield ['6', _ID_, row[3], _BP_]



def GEN_stitch_GD(_df):

    # print(_df.head())

    ## Version B: Collapse special variables into mid-point, while keeping other variables.

    epsilon = 1E-12
    n_rows = _df.shape[0]

    l_rows = [list(row) for row in _df.itertuples()]

    _GD_ = _df.iloc[:, 2].astype(float).tolist()
    # print(_GD_)


    l_new_GD = [0]

    i = 1

    while i < n_rows:

        # Start of the chunk(HLA marker region)
        if _GD_[i-1] != 0 and _GD_[i] == 0:

            start = i
            start_cap_SNP = start - 1

            end = i + 1

            # Finding end of HLA allele binary marker
            while _GD_[end] == 0:
                end += 1

            end -= 1

            end_cap_SNP = end + 1

            # print("End of HLA allele binary marker is : \n{}".format(_df.iloc[end, :]))
            # print("Next marker of this End point is : {}".format(_df.iloc[end+1, :]))


            mid = mean([float(l_rows[start_cap_SNP][3]), float(l_rows[end_cap_SNP][3])])
            acc = mid

            l_new_GD.append(mid) # First HLA allele binary marker's GD

            j = start + 1

            while j < end_cap_SNP:

                acc += epsilon
                l_new_GD.append(acc)

                j += 1

            i = end

        else:

            l_new_GD.append(_GD_[i])

        i += 1


    # print(len(l_new_GD))
    # print(n_rows)


    return pd.concat([_df.iloc[:, [0,1]].reset_index(drop=True), pd.Series(l_new_GD), _df.iloc[:, 3].reset_index(drop=True)], axis=1)








if __name__ == '__main__':

    """
    Input : (1) GeneticMap file, (2) (GCchanged) reference marker file, (3) Output file name.
    Output : Merged GeneticMap file
    """

    [_GeneticMap, _REF_Markers, _out] = sys.argv[1:]

    # ex)
    # [_GeneticMap, _REF_Markers, _out] = ['/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/CookHLA/data/HLA_PANEL/Genetic_map/CEU_T1DGC.mach_step.avg.clpsB',
    #                                      '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190710_BOTH/T1DGC_REF.exon234.markers',
    #                                      '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190710_BOTH/HLA_Collapsed.txt']

    Make_EXON234_AGM(_GeneticMap, _REF_Markers, _out)

