#-*- coding: utf-8 -*-

import os, sys, re

import pandas as pd


std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]

class REFERENCE():

    def __init__(self, _reference, _which=(1,1,1,1,1,1)):

        """
        Class to take and manage 'MakeReference' reference panel.
        """

        ########## < Loading Reference panel data > ##########

        ### Main data
        self.bed = None
        self.bim = None
        self.fam = None
        self.FRQ = None
        self.bgl_phased = None
        self.markers = None


        ### File existence check.

        # (1) *.bed
        if not os.path.exists(_reference+'.bed'):
            print(std_ERROR_MAIN_PROCESS_NAME + "Given Reference PLINK bed file('{}') can't be found.\n"
                                                "Please check the '--reference/-ref' argument again.".format(_reference+'.bed'))
            sys.exit()

        # (2) *.bim
        if not os.path.exists(_reference+'.bim'):
            print(std_ERROR_MAIN_PROCESS_NAME + "Given Reference PLINK bim file('{}') can't be found.\n"
                                                "Please check the '--reference/-ref' argument again.".format(_reference+'.bim'))
            sys.exit()

        # (3) *.fam
        if not os.path.exists(_reference+'.fam'):
            print(std_ERROR_MAIN_PROCESS_NAME + "Given Reference PLINK fam file('{}') can't be found.\n"
                                                "Please check the '--reference/-ref' argument again.".format(_reference+'.fam'))
            sys.exit()

        # (4) *.bgl.phased
        if not os.path.exists(_reference+'.bgl.phased'):
            print(std_ERROR_MAIN_PROCESS_NAME + "Given Reference Beagle Phased file('{}') can't be found.\n"
                                                "Please check the '--reference/-ref' argument again.".format(_reference+'.bgl.phased'))
            sys.exit()

        # (5) *.markers
        if not os.path.exists(_reference+'.markers'):
            print(std_ERROR_MAIN_PROCESS_NAME + "Given Reference Beagle Markers file('{}') can't be found.\n"
                                                "Please check the '--reference/-ref' argument again.".format(_reference+'.markers'))
            sys.exit()

        # (6) *.FRQ.frq
        if not os.path.exists(_reference+'.FRQ.frq'):
            print(std_ERROR_MAIN_PROCESS_NAME + "Given Reference Allele frequency information file('{}') can't be found.\n"
                                                "Please check the '--reference/-ref' argument again.".format(_reference+'.FRQ.frq'))
            sys.exit()
