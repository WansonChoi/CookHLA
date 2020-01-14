#-*- coding: utf-8 -*-

import os, sys, re
import subprocess
from shutil import which
import pandas as pd


std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]

# dependent software
_p_plink = which('plink')
_p_beagle4 = which('beagle')


class REFERENCE():

    def __init__(self, _reference, _which=(1,1,1,1,1,1)):

        """
        Class to take and manage 'MakeReference' reference panel.
        """

        ########## < Loading Reference panel data > ##########

        ### Main data
        self.prefix = _reference

        self.bed = None
        self.bim = None
        self.fam = None
        self.FRQ = None
        self.bgl_phased = None
        self.markers = None


        ### File existence check and loading.

        if _which[0]:

            # (1) *.bed
            if not os.path.exists(_reference+'.bed'):
                print(std_ERROR_MAIN_PROCESS_NAME + "Given Reference PLINK bed file('{}') can't be found.\n"
                                                    "Please check the '--reference/-ref' argument again.".format(_reference+'.bed'))
                sys.exit()


        if _which[1]:

            # (2) *.bim
            if not os.path.exists(_reference+'.bim'):
                print(std_ERROR_MAIN_PROCESS_NAME + "Given Reference PLINK bim file('{}') can't be found.\n"
                                                    "Please check the '--reference/-ref' argument again.".format(_reference+'.bim'))
                sys.exit()

            else:
                self.bim = pd.read_csv(_reference+'.bim', sep='\s+', header=None, names=['Chr', 'Label', 'GD', 'BP', 'al1', 'al2'])
                # print(self.bim.head())

        if _which[2]:

            # (3) *.fam
            if not os.path.exists(_reference+'.fam'):
                print(std_ERROR_MAIN_PROCESS_NAME + "Given Reference PLINK fam file('{}') can't be found.\n"
                                                    "Please check the '--reference/-ref' argument again.".format(_reference+'.fam'))
                sys.exit()

            else:
                self.fam = pd.read_csv(_reference+'.fam', sep='\s+', header=None, names=['FID', 'IID', 'PID', 'MID', 'Sex', 'Phe'])


        if _which[3]:

            # (4) *.bgl.phased
            if not os.path.exists(_reference+'.bgl.phased'):
                print(std_ERROR_MAIN_PROCESS_NAME + "Given Reference Beagle Phased file('{}') can't be found.\n"
                                                    "Please check the '--reference/-ref' argument again.".format(_reference+'.bgl.phased'))
                sys.exit()

            else:
                self.bgl_phased = pd.read_csv(_reference+'.bgl.phased', '\s+', header=None)


        if _which[4]:

            # (5) *.markers
            if not os.path.exists(_reference+'.markers'):
                print(std_ERROR_MAIN_PROCESS_NAME + "Given Reference Beagle Markers file('{}') can't be found.\n"
                                                    "Please check the '--reference/-ref' argument again.".format(_reference+'.markers'))
                sys.exit()

            else:
                self.markers = pd.read_csv(_reference+'.markers', sep='\s+', header=None)


        if _which[5]:

            # (6) *.FRQ.frq
            if not os.path.exists(_reference+'.FRQ.frq'):
                print(std_ERROR_MAIN_PROCESS_NAME + "Given Reference Allele frequency information file('{}') can't be found.\n"
                                                    "Please check the '--reference/-ref' argument again.".format(_reference+'.FRQ.frq'))
                sys.exit()

            else:
                self.FRQ = pd.read_csv(_reference+'.FRQ.frq', sep='\s+', header=0)



    ### class methods

    def get_MKref_markers(self, _out=None, _AA=True, _HLA=True, _SNP=True, _INS=True):

        l_toExtract = []

        if _AA:
            l_toExtract.append('AA_')
        if _HLA:
            l_toExtract.append('HLA_')
        if _SNP:
            l_toExtract.append('SNP_')
        if _INS:
            l_toExtract.append('INS_')

        if len(l_toExtract) == 0:
            return -1


        p_toExtract = re.compile(r'%s' % '|'.join(l_toExtract))
        f_toExtract = self.bim['Label'].str.match(p_toExtract)

        df_MKRef_markers = self.bim[f_toExtract]

        if df_MKRef_markers.shape[0] == 0:
            print(std_WARNING_MAIN_PROCESS_NAME + "There is no 'MakeReference' markers in given reference panel.")
            return -1
        else:
            if bool(_out):
                df_MKRef_markers['Label'].to_csv(_out, header=False, index=False)
                return _out
            else:
                return df_MKRef_markers



    def PLINK_subset(self, _out, _toKeep=None, _toRemove=None, _toExtract=None, _toExclude=None):

        if not (bool(_toKeep) or bool(_toRemove) or bool(_toExtract) or bool(_toExclude)):
            print("Nothing to subset. (Samples or Markers not given.)")
            return -1

        command = [
            _p_plink, '--make-bed',
            '--bfile', self.prefix,
            '--out', _out
        ]

        if bool(_toKeep):
            command.extend(['--keep', _toKeep])
        elif bool(_toRemove):
            command.extend(['--remove', _toRemove])

        if bool(_toExtract):
            command.extend(['--extract', _toExtract])
        elif bool(_toExclude):
            command.extend(['--exclude', _toExclude])

        # print(command)
        sub = subprocess.call(command, stdout=subprocess.DEVNULL)

        if sub == 0:
            return _out
        else:
            return -1



    def getArtificialGDBIM(self, _start_offset = 0, _interval=1E-5, __writeFile=False):

        t_bim = self.bim.copy()

        sr_GD = pd.Series([
            (_start_offset + _interval * i) for i in range(t_bim.shape[0])
        ], name='GD')

        t_bim['GD'] = sr_GD

        if __writeFile:
            # t_bim.to_csv(self.prefix+'.')
            pass
        else:
            return t_bim


