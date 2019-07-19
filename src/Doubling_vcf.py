#-*- coding: utf-8 -*-

import os, sys, re
import pandas as pd


def Doubling_vcf(_input_vcf_body, _output_vcf_body):

    ### Loading input vcf file.

    __vcf_body__ = pd.read_csv(_input_vcf_body, sep='\s+', header=0, index_col=[0,1,2,3,4,5,6,7,8], dtype=str)
    # print(__vcf_body__.head())


    ### Doubling input vcf file.

    df_vcf_doubled = pd.concat([item[1].str.extract(r'^(\d)\|(\d)', expand=True).applymap(lambda x : '|'.join([x,x])) for item in __vcf_body__.iteritems()], axis=1)
    # print(df_vcf_doubled.head())


    ### Making Header part

    header_input = __vcf_body__.columns.tolist()
    header_output = ['_'.join([item, j]) for item in header_input for j in ['1', '2']]

    df_vcf_doubled.columns = header_output
    # print(df_vcf_doubled.head())

    df_vcf_doubled.to_csv(_output_vcf_body, sep='\t', header=True, index=True)



    return _output_vcf_body



if __name__ == '__main__':

    """
    Doubling each column(each sample) into two columns after 9th column.
    
    ex.
    884MF17 = [
        0|1, 
        0|1, 
        1|1, 
        0|1, 
        ...
    ]
    
    =>
    
    884MF17_1 = [
        0|0,
        0|0,
        1|1,
        0|0,
        ...
    ]
    
    and
    
    884MF17_2 = [
        1|1,
        1|1,
        1|1,
        1|1,
        ...
    ]
    
    
    """

    [input_vcf_body, output_vcf_body] = sys.argv[1:]

    Doubling_vcf(input_vcf_body, output_vcf_body)