#-*- coding: utf-8 -*-

"""
2020. 12. 24.

"""

import os, sys, re


def getSampleNumbers(_fam):

    with open(_fam, 'r') as f_fam:
        return len(list(f_fam))



def fixLabel(_input, _reference, _out, _PLINK):

    """
    Assuming BP information is given properly, Replace the labels of target markers with the matched ones of reference markers.

    """

    command = ' '.join(
        [_PLINK, '--make-bed', '--bfile', _input, '--keep-allele-order', '--out', _out])
    # print(command)
    os.system(command)


    return _out