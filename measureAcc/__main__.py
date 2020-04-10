#-*- coding: utf-8 -*-

import os, sys

from NomenCleaner.__main__ import HATK_NomenCleaner
from measureAcc.SieveCHPED import SieveCHPED
from measureAcc.measureAccuracy import measureAccuracy


class CookHLA_measureAcc(object):

    def __init__(self):

        """
        1. imputed *.alleles
        2. answer *.alleles

        (1) ALLELES2HPED.py

        (2) DRB1 1454 -> 1401

        (3) NomenCleaner

        (4) SieveCHEPD.py

        (5) measureAccuracy_v3.5

        """

        pass



if __name__ == '__main__':

    """
    
    measureAccuracy.py (v3.5)
    
    
    """


    [_answer_Marked_chped, _imputed_Marked_chped, _out, _allele_group] = sys.argv[1:]

    # measureAccuracy(_answer_Marked_chped, _imputed_Marked_chped, _out, _allele_group)