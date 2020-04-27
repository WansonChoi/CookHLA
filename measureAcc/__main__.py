#-*- coding: utf-8 -*-

import os, sys, re
from os.path import exists, join, dirname, basename

from NomenCleaner.__main__ import HATK_NomenCleaner
from measureAcc.SieveCHPED import SieveCHPED
from measureAcc.ALLELES2HPED import ALLELES2HPED
from measureAcc.measureAccuracy import measureAccuracy

from src.CookHLAError import CookHLAInputPreparationError



class CookHLA_measureAcc(object):

    def __init__(self, _answer, _imputed, _out, _allele_group='measureAcc/SameExon234.MERGED.nuc.txt'):

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
            t = HATK_NomenCleaner(_f, "NomenCleaner/HLA_ALLELE_TABLE.imgt3320.hat", '3320', t1_out,
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