#-*- coding: utf-8 -*-

"""
2020. 12. 24.

"""

import os, sys, re
import subprocess

import pandas as pd
from pyliftover import LiftOver

from src.CookHLAError import CookHLAInputPreparationError

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))


Ambiguous = ({'A', 'T'}, {'G', 'C'})
dict_Complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}


p_MKref = re.compile(r'^(AA_|HLA_|SNP_|INS_)') # To exclude 'MakeReference' markers.


def getSampleNumbers(_fam):

    with open(_fam, 'r') as f_fam:
        return len(list(f_fam))



def LiftDown_hg18(_bim, _hg, _out):

    HG_input = 'hg{}'.format(_hg)
    # print("HG: {}".format(HG_input))

    df_bim = pd.read_csv(_bim, sep='\s+', header=None, dtype=str, names=['Chr', 'Label', 'GD', 'BP', 'a1', 'a2'])
    # print("df_bim:\n{}\n".format(df_bim))


    ### Main Liftover ###

    if HG_input == 'hg38':

        """
        'hg38' -> 'hg19' -> 'hg18' is needed.
        The Liftover tool (by UCSC Genomics Institute) doesn't provide 'hg38' to 'hg18'.
        
        """

        lo_hg38_to_hg19 = LiftOver(HG_input, 'hg19')
        lo_hg19_to_hg18 = LiftOver('hg19', 'hg18')

        sr_hg19 = df_bim['BP'] \
            .astype(int) \
            .map(lambda x: lo_hg38_to_hg19.convert_coordinate('chr6', x)) \
            .map(lambda x: x[0][1] if len(x) > 0 else -1)
        # print("(hg19):\n{}\n".format(sr_hg19))

        sr_hg18 = sr_hg19 \
            .map(lambda x: lo_hg19_to_hg18.convert_coordinate('chr6', x)) \
            .map(lambda x: x[0][1] if len(x) > 0 else -1)
        # print("(hg18):\n{}\n".format(sr_hg18))


    else:

        lo = LiftOver(HG_input, 'hg18')  # Liftdown to hg18

        sr_hg18 = df_bim['BP'] \
            .astype(int) \
            .map(lambda x: lo.convert_coordinate('chr6', x)) \
            .map(lambda x: x[0][1] if len(x) > 0 else -1)


    df_bim['BP'] = sr_hg18 # Setting new BPs (Liftdown)



    ### Makrers that failed the Liftdown. ###

    f_failed = sr_hg18 == -1

    if f_failed.any():
        print(std_WARNING_MAIN_PROCESS_NAME + "Next markers of Target('{}') failed to Liftdown to hg18. These markers will be excluded."
              .format(_bim))
        print(df_bim[f_failed])

    # print("df_bim_hg18:\n{}\n".format(df_bim))
    df_bim.to_csv(_out, sep='\t', header=False, index=False)


    return _out



def UpdateInput(_bim_input, _bim_reference, _out):

    df_bim_input = pd.read_csv(_bim_input, sep='\s+', header=None, dtype=str, names=['Chr', 'Label', 'GD', 'BP', 'a1', 'a2'])
    #     print("bim_input:\n{}\n".format(df_bim_input))

    df_bim_reference = pd.read_csv(_bim_reference, sep='\s+', header=None, dtype=str, names=['Chr', 'Label', 'GD', 'BP', 'a1', 'a2'])
    f_MKref = df_bim_reference['Label'].str.match(p_MKref)
    df_bim_reference = df_bim_reference[~f_MKref]
    # print("bim_reference:\n{}\n".format(df_bim_reference))

    ## Merge based on BP (hg18)
    df_merge0 = df_bim_input[['Label', 'BP', 'a1', 'a2']].merge(df_bim_reference[['Label', 'BP', 'a1', 'a2']], on='BP') \
        .drop_duplicates('BP')
    # print("df_merge0:\n{}\n".format(df_merge0))
    #     df_merge0.to_csv(_out+'.merged', sep='\t', header=True, index=False)


    if df_merge0.shape[0] == 0:
        print(std_ERROR_MAIN_PROCESS_NAME + "Target('{}') and Reference('{}') don't have any markers matched by Base Position. Please check their Human Genome(hg) version again."
              .format(_bim_input, _bim_reference))
        sys.exit()


    ### Main 0 - Extract
    sr_extract = df_merge0['Label_x']
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


def FixInput(_input, _hg_input, _reference, _out, _PLINK):

    """
    BP information is assumed to be given properly.
    BP of the reference panel will be assumed to be hg 18.

    """

    if _hg_input != "18":

        # Liftdown target(input) HG to hg 18.
        INPUT_BIM_hg18 = LiftDown_hg18(_input+'.bim', _hg_input, _out+'.bim.LiftDown_hg18')

        [t_extract, t_update_name, t_update_alleles] = UpdateInput(INPUT_BIM_hg18, _reference+'.bim', _out+'.fix')

        _out = _out+'.LiftDown_hg18'



        # (1) Subset target markers to reference markers based on BP.

        command = ' '.join(
            [_PLINK, '--make-bed',
             '--bed', _input+'.bed',
             '--bim', INPUT_BIM_hg18,
             '--fam', _input+'.fam',
             '--extract {}'.format(t_extract),
             '--out', _out+'.subset'])
        # print(command)

        try:
            subprocess.run(command.split(), check=True, stdout=subprocess.DEVNULL)

        except subprocess.CalledProcessError:
            raise CookHLAInputPreparationError(
                std_ERROR_MAIN_PROCESS_NAME +
                "Subsetting target markers to reference markers based on BP failed. "
                "Please check the PLINK log file('{}').".format(_out+'.subset.log')
            )



        # (2) Replace a target marker label with the matched reference marker label.
        command = ' '.join(
            [_PLINK, '--make-bed',
             '--bfile', _out + '.subset',
             '--update-name {}'.format(t_update_name),
             '--update-alleles {}'.format(t_update_alleles),
             '--out', _out])
        # print(command)

        try:
            subprocess.run(command.split(), check=True, stdout=subprocess.DEVNULL)

        except subprocess.CalledProcessError:
            raise CookHLAInputPreparationError(
                std_ERROR_MAIN_PROCESS_NAME +
                "Replacing a target marker label with the matched reference marker label failed. "
                "Please check the PLINK log file('{}').".format(_out + '.log')
            )
        else:
            os.system('rm {}'.format(INPUT_BIM_hg18)) # Liftover


            os.system('rm {}'.format(t_extract))
            os.system('rm {}'.format(_out+'.subset.bed'))
            os.system('rm {}'.format(_out+'.subset.bim'))
            os.system('rm {}'.format(_out+'.subset.fam'))
            os.system('rm {}'.format(_out+'.subset.log'))

            os.system('rm {}'.format(t_update_name))
            os.system('rm {}'.format(t_update_alleles))




    else:

        # No liftover

        [t_extract, t_update_name, t_update_alleles] = UpdateInput(_input+'.bim', _reference+'.bim', _out+'.fix')



        # (1) Subset target markers to reference markers based on BP.

        command = ' '.join(
            [_PLINK, '--make-bed',
                    '--bfile', _input,
                    '--extract {}'.format(t_extract),
                    '--out', _out+'.subset'])
        # print(command)

        try:
            subprocess.run(command.split(), check=True, stdout=subprocess.DEVNULL)

        except subprocess.CalledProcessError:
            raise CookHLAInputPreparationError(
                std_ERROR_MAIN_PROCESS_NAME +
                "Subsetting target markers to reference markers based on BP failed. "
                "Please check the PLINK log file('{}').".format(_out+'.subset.log')
            )



        # (2) Replace a target marker label with the matched reference marker label.

        command = ' '.join(
            [_PLINK, '--make-bed',
                    '--bfile', _out+'.subset',
                    '--update-name {}'.format(t_update_name),
                    '--update-alleles {}'.format(t_update_alleles),
                    '--out', _out])
        # print(command)

        try:
            subprocess.run(command.split(), check=True, stdout=subprocess.DEVNULL)

        except subprocess.CalledProcessError:
            raise CookHLAInputPreparationError(
                std_ERROR_MAIN_PROCESS_NAME +
                "Replacing a target marker label with the matched reference marker label failed. "
                "Please check the PLINK log file('{}').".format(_out + '.log')
            )
        else:
            os.system('rm {}'.format(t_extract))
            os.system('rm {}'.format(_out+'.subset.bed'))
            os.system('rm {}'.format(_out+'.subset.bim'))
            os.system('rm {}'.format(_out+'.subset.fam'))
            os.system('rm {}'.format(_out+'.subset.log'))

            os.system('rm {}'.format(t_update_name))
            os.system('rm {}'.format(t_update_alleles))



    return _out