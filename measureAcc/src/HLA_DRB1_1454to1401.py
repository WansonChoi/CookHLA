#-*- coding: utf-8 -*-

import os, sys
import pandas as pd


def HPED_DRB1_1454to1401(_hped, _out=None):


    HPED = pd.read_csv(_hped, sep='\s+', header=None, dtype=str)
    # print("HPED:\n{}\n".format(HPED.head()))

    if HPED.shape[1] != 22:
        print("The number of columns of given HPED file('{}') has to be 22.".format(_hped))
        return -1



    df_DRB1_1454to1401 = HPED.iloc[:, [20, 21]].replace('1454', '1401')

    df_RETURN = pd.concat([HPED.iloc[:, :20], df_DRB1_1454to1401], axis=1)

    if bool(_out):

        df_RETURN.to_csv(_out+'.hped', sep='\t', header=False, index=False)
        return _out+'.hped'

    else:
        return df_RETURN







if __name__ == '__main__':


    [_hped, _out] = sys.argv[1:]

    HPED_DRB1_1454to1401(_hped, _out)
