#-*- coding: utf-8 -*-

import os, sys, re
import pandas as pd

########## < Core Varialbes > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]



class HLA_MultipleRefs():

    def __init__(self, __exonN__, _df_EXON_info, __REFERENCE__, _df_reference_bim, _out, _hg,
                 _PLINK, _LINKAGE2BEAGLE, *args, **kwargs):

        """

        (1) __REFERENCE__.bed
        (2) __REFERENCE__.bim
        (3) __REFERENCE__.fam
        (4) __REFERENCE__.bgl.phased
        (5) __REFERENCE__.markers
        (6) __REFERENCE__.FRQ.frq

        """

        OUTPUT_dir = os.path.dirname(_out)

        # print(std_MAIN_PROCESS_NAME + "Locating HLA markers in {}".format(__exonN__))
        # print("Given Exon information : \n{}".format(_df_EXON_info))
        # print("Given Reference BIM file : \n{}".format(_df_reference_bim.head()))


        self.Di_HLA_markers = {HLA_names[i] : None for i in range(0, len(HLA_names))}

        for i in range(0, len(HLA_names)):
        # for i in range(0, 2):

            # print(HLA_names[i])
            p_HLA_markers = re.compile(r'^HLA_{}'.format(HLA_names[i]))
            f_HLA_markers = _df_reference_bim.iloc[:, 1].str.match(p_HLA_markers)

            df_temp = _df_reference_bim.loc[f_HLA_markers, 'Label'].reset_index(drop=True)
            # print(df_temp.head())
            sr_ExonPos = pd.Series([_df_EXON_info.loc[HLA_names[i], 'mid'] for z in range(0, df_temp.shape[0])], name='mid')
            # print(sr_ExonPos.head())

            self.Di_HLA_markers[HLA_names[i]] = pd.concat([df_temp, sr_ExonPos], axis=1)
            # print(df_temp.head())
            # print('================================')



        to_modify = _out+'.{}_mid.txt'.format(__exonN__)
        pd.concat([self.Di_HLA_markers[HLA_names[i]] for i in range(0, len(HLA_names))], axis=0).to_csv(to_modify, sep='\t', header=False, index=False)


        __RETURN__ = os.path.join(OUTPUT_dir, os.path.basename(__REFERENCE__)+'.{}'.format(__exonN__))
        command = ' '.join([_PLINK, '--make-bed --bfile {} --update-map {} --out {}'.format(__REFERENCE__, to_modify, __RETURN__)])
        # print(command)

        if not os.system(command):
            os.system(' '.join(['rm', to_modify]))
        else:
            print(std_ERROR_MAIN_PROCESS_NAME + "Failed to generate {} reference.".format(__exonN__))
            sys.exit()



        ### Converting ExonN reference to Beagle format. (*.bgl.phased, *.markers)

        # *.markers
        command = ' '.join(['awk \'{print $2"\t"$4"\t"$5"\t"$6}\'', '{} > {}'.format(__RETURN__ + '.bim', __RETURN__ + '.markers')])
        # print(command)
        os.system(command)


        # Plink bed to ped file.
        command = '{} --bfile {} --recode --out {}'.format(_PLINK, __RETURN__, __RETURN__)
        print(command)
        if not os.system(command):
            os.system('rm {}'.format(__RETURN__+'.{bed,bim,fam,log}'))


        command = ' '.join(['awk \'{print "M " $2}\'', __RETURN__+'.map', '>', __RETURN__+'.dat'])
        # print(command)
        if not os.system(command):
            os.system('rm {}'.format(__RETURN__ + '.map'))


        command = ' '.join(["cut -d ' ' -f1-5,7-", __RETURN__+'.ped', '>', __RETURN__+'.nopheno.ped'])
        # print(command)
        if not os.system(command):
            os.system('rm {}'.format(__RETURN__ + '.ped'))


        ### Linkage2Beagle
        command = '{} pedigree={} data={} beagle={} standard=true > {}'.format(_LINKAGE2BEAGLE, __RETURN__+'.nopheno.ped', __RETURN__+'.dat', __RETURN__+'.bgl', __RETURN__+'.bgl.log')
        print(command)
        if not os.system(command):
            os.system('rm {}'.format(__RETURN__ + '.{dat,nopheno.ped,markers}'))
            os.system('rm {}'.format(__RETURN__ + '.bgl.log'))





        self.MODIFIED_REF = __RETURN__

    def getOUTPUT(self):
        return self.MODIFIED_REF