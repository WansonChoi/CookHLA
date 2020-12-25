#-*- coding: utf-8 -*-

"""
2020. 12. 24.

"""

import os, sys, re
import pandas as pd


def getSampleNumbers(_fam):

    with open(_fam, 'r') as f_fam:
        return len(list(f_fam))


Ambiguous = ({'A', 'T'}, {'G', 'C'})
dict_Complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}


def UpdateInput(_bim_input, _bim_reference, _out):

    df_bim_input = pd.read_csv(_bim_input, sep='\s+', header=None, dtype=str, names=['Chr', 'Label', 'GD', 'BP', 'a1', 'a2'])
    #     print("bim_input:\n{}\n".format(df_bim_input))

    df_bim_reference = pd.read_csv(_bim_reference, sep='\s+', header=None, dtype=str, names=['Chr', 'Label', 'GD', 'BP', 'a1', 'a2'])
    #     print("bim_reference:\n{}\n".format(df_bim_reference))

    df_merge0 = df_bim_input[['Label', 'BP', 'a1', 'a2']].merge(df_bim_reference[['Label', 'BP', 'a1', 'a2']], on='BP') \
        .drop_duplicates('BP')
    # print("df_merge0:\n{}\n".format(df_merge0))
    #     df_merge0.to_csv(_out+'.merged', sep='\t', header=True, index=False)



    ### Main 0 - Extract
    sr_extract = df_merge0['Label_y']
    # print("sr_extract:\n{}\n".format(sr_extract))
    sr_extract.to_csv(_out + '.extract', header=False, index=False)



    ### Main 1 - Update Name
    df_update_name = df_merge0[['Label_x', 'Label_y']]
    # print("df_update_name:\n{}\n".format(df_update_name))
    df_update_name.to_csv(_out + '.update_name', sep='\t', header=False, index=False)



    ### Main iteration 2 - Update Alleles
    count = 0
    l_update_alleles = []
    for line in df_merge0.itertuples():

        #         print(line)
        [idx, Label_x, BP, a1_x, a2_x, Label_y, a1_y, a2_y] = line

        if a1_x == '0' and a2_x == '0':
            # Both alleles of Target are '0' => No way to save.
            continue

        if a1_y == '0' or a2_y == '0':
            # Either allele of Reference is '0' => No way to save.
            continue



        if a1_x == '0' or a2_x == '0':
            # Either one allele of target is '0' (ex. ('0', 'G') or ('A', '0'))

            l_temp = [Label_y, a1_x, a2_x]

            if a1_x == '0':  # (ex. ('0', 'G'))

                # Flip check
                if not (a2_x == a1_y or a2_x == a2_y):
                    a2_x = dict_Complement[a2_x]

                if a2_x == a1_y:  # (a1_x, a2_x) == ('0', 'G') and (a1_y, a2_y) == ('G', 'T')
                    l_temp.extend([a2_y, a1_y])
                elif a2_x == a2_y:  # (a1_x, a2_x) == ('0', 'G') and (a1_y, a2_y) == ('T', 'G')
                    l_temp.extend([a1_y, a2_y])
                else:
                    # Abnormal relationship
                    continue



            elif a2_x == '0':  # (ex. ('G', '0'))

                # Practically won't be here. Because plink usually set '0' allele as A1.

                # Flip check
                if not (a1_x == a1_y or a1_x == a2_y):
                    a1_x = dict_Complement[a1_x]

                if a1_x == a1_y:  # (a1_x, a2_x) == ('G', '0') and (a1_y, a2_y) == ('G', 'T')
                    l_temp.extend([a1_y, a2_y])
                elif a1_x == a2_y:  # (a1_x, a2_x) == ('G', '0') and (a1_y, a2_y) == ('T', 'G')
                    l_temp.extend([a2_y, a1_y])
                else:
                    # Abnormal relationship
                    continue

            l_update_alleles.append(l_temp)



        count += 1
    #         if count > 5: break

    df_update_alleles = pd.DataFrame(l_update_alleles)
    # print("df_update_alleles:\n{}\n".format(df_update_alleles))
    df_update_alleles.to_csv(_out + '.update_alleles', sep='\t', header=False, index=False)

    return (_out + '.extract', _out + '.update_name', _out + '.update_alleles')


def FixInput(_input, _reference, _out, _PLINK):

    """
    BP information is assumed to be given properly.

    """

    [t_extract, t_update_name, t_update_alleles] = UpdateInput(_input+'.bim', _reference+'.bim', _out+'.fix')

    command = ' '.join(
        [_PLINK, '--make-bed',
                '--bfile', _input,
                '--extract {}'.format(t_extract),
                '--update-name {}'.format(t_update_name),
                '--update-alleles {}'.format(t_update_alleles),
                '--out', _out])
    print(command)
    os.system(command)


    os.system('rm {}'.format(t_extract))
    os.system('rm {}'.format(t_update_name))
    os.system('rm {}'.format(t_update_alleles))


    return _out