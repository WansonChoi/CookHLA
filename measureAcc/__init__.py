#-*- coding: utf-8 -*-

"""
measureAccuracy_v3.5

v1 : the original version introduced by Buhm Han.
v2 : calculates the only allele pair whose both answer alleles are perfectly available.

v3 : calculates accuracy based on (1) nucleotide or (2) amino acid sequence of exon 2,3,(4) in each HLA gene.
v3.5 : considers 'deprecated' HLA alleles due to NomenCleaner. In fact, it is common for HLA type data to contain
        HLA alleles in old nomenclature. Nomencleaner deprecates those alleles. So, With the module 'SieveCHPED.py',
        v3.5 considers two cases where 'deprecated' allele appears in (1) imputed allele pair and (2) answer allele pair.


Additionally, Based on the fact ('https://en.wikiversity.org/wiki/Genetics/Human_Leukocyte_Antigen/DRB1*14:54'),
HLA-DRB1*1401 will be considered to be same as HLA-DRB1*1454. However, HLA-DRB1*1451 will be modified to 1401 for
practical purpose. (applied from the v3).


"""

