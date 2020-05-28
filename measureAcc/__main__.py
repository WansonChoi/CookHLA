#-*- coding: utf-8 -*-

import os, sys, re


from measureAcc.measureAccuracy import CookHLA_measureAcc


if __name__ == '__main__':

    """
    measureAccuracy.py (v3.5)
    
    - takes (1) imputed result and (2) answer file in '*.Marked.chped' format
    - takes Allele group file.
        (Alleles classifed based on the nucleotide sequence of exon 2,3,(4) of each HLA gene)
    - calculates accuracy.
    
    """


    [_answer_Marked_chped, _imputed_Marked_chped, _out] = sys.argv[1:]

    CookHLA_measureAcc(_answer_Marked_chped, _imputed_Marked_chped, _out)