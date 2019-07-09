#-*- coding: utf-8 -*-

import os, sys, re
import multiprocessing as mp

from src.Make_EXON234_Panel import Make_EXON234_Panel



########## < Core Varialbes > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
__EXON__ = ['exon2']
# __EXON__ = ['exon2', 'exon3', 'exon4']


# Patterns to use.
p_1st_two_columns = re.compile(r'^(\S+)\s+(\S+)') # First two columns (ex. 'P pedigree' or 'rs969931 29602876', 'M rs969931', ... )
p_ExonN = {'exon2' : re.compile(r'^HLA_\w+_\d+_exon2'),
           'exon3' : re.compile(r'^HLA_\w+_\d+_exon3'),
           'exon4' : re.compile(r'^HLA_\w+_\d+_exon4')}




class HLA_MultipleRefs():

    def __init__(self, __REFERENCE__, _out, _hg, _BEAGLE2LINKAGE, _PLINK, _MultP=1,
                 f_save_intermediates=False):

        """

        (1) __REFERENCE__.bed
        (2) __REFERENCE__.bim
        (3) __REFERENCE__.fam
        (4) __REFERENCE__.bgl.phased
        (5) __REFERENCE__.markers
        (6) __REFERENCE__.FRQ.frq

        """

        # general
        self.__save_intermediates = f_save_intermediates

        # dependency
        self.BEAGLE2LINKAGE = _BEAGLE2LINKAGE
        self.PLINK = _PLINK

        # Main panel data.
        self.EXON234_Panel = Make_EXON234_Panel(__REFERENCE__, _out+'.exon234', _BEAGLE2LINKAGE, _PLINK)
        # self.EXON234_Panel = "/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190708_MM/T1DGC_REF.exon234"    # [Temporary Hardcoding.]
        # print("[Temporary Hardcoding of EXON234_Panel] :\n{}".format(self.EXON234_Panel))
        self.ExonN_Panel = {_exonN : None for _exonN in __EXON__}



        ### Generating Exon 2, 3, 4 Panel

        if _MultP == 1:

            for _exonN in __EXON__:

                self.ExonN_Panel[_exonN] = self.Make_ExonN_Panel(_exonN, self.EXON234_Panel, _out+'.{}'.format(_exonN))

        else:

            ## Multiprocessing
            pool = mp.Pool(processes=_MultP)

            dict_Pool = {_exonN: pool.apply_async(
                self.Make_ExonN_Panel, (_exonN, self.EXON234_Panel, _out+'.{}'.format(_exonN))
            ) for _exonN in __EXON__}

            pool.close()
            pool.join()

            self.ExonN_Panel = {_exonN_: _OUT.get() for _exonN_, _OUT in dict_Pool.items()}



    def Make_ExonN_Panel(self, _exonN, _EXON234_Panel, _out):

        """
        (1) Make a regular expression to nominate HLA allele binary marker of exon 2, 3, 4
        (2) To make exon 2 panel, for example, subset out exon 3, 4 HLA allele binary markers. Maybe with Plink('--recode').
        (3) With subsetted *.ped and *.map files, make a *.nopheno and *.markers files.
        (4) Feed them to LINKAGE2BEAGLE => Beagle file generated.
        (5) Get a frequency file.
        """


        ### Subsetting `_exonN` HLA allele binary markers and SNP markers. ((1) + (2) above.)

        # Input files - (1) Exon234 beagle file, (2) Exon234 markers file.
        f_bgl_exon234 = open(_EXON234_Panel+'.bgl.phased', 'r')
        f_markers_exon234 = open(_EXON234_Panel+'.markers', 'r')

        # Output files
        f_out_bgl = open(_out+'.bgl.phased', 'w')
        f_out_markers = open(_out+'.markers', 'w')

        count = 0

        for line_bgl_exon234 in f_bgl_exon234:

            m = p_1st_two_columns.match(line_bgl_exon234)
            items_bgl_exon234 = [m.group(1), m.group(2)] # Skip exception where that pattern finds nothing.

            if items_bgl_exon234[0] != 'M':
                """
                - Header part of *.bgl file.
                
                'P pedigree'
                'I id'
                'fID father'
                'mID mother'
                'C gender'
                
                etc.
                
                """

                f_out_bgl.write(line_bgl_exon234) # Just forward the line to output file.

            else:
                """
                Main body
                
                'M   rs969931'
                'M  rs2745406'
                
                """

                # Each line of '*.markers' file.
                line_markers_exon234 = f_markers_exon234.readline()

                if re.match(pattern=r'^HLA_', string=items_bgl_exon234[1]):

                    if p_ExonN[_exonN].match(items_bgl_exon234[1]):
                        # Gather only `exonN` HLA allele binary markers.
                        f_out_bgl.write(line_bgl_exon234)
                        f_out_markers.write(line_markers_exon234)

                else:
                    # In case of normal SNP markers, just forward them.
                    f_out_bgl.write(line_bgl_exon234)
                    f_out_markers.write(line_markers_exon234)




            count += 1
            # if count > 10 : break


        f_bgl_exon234.close(); f_markers_exon234.close(); f_out_bgl.close(); f_out_markers.close();

        # Generated exonN Panel.
        bgl_exonN = _out+'.bgl.phased'
        markers_exonN = _out+'.markers'




        ### With those subsetted *.bgl and *.markers files, generate *.ped and *.map. ((3) + (4) above.)


        command = 'cat {} | {} {}'.format(bgl_exonN, self.BEAGLE2LINKAGE,
                                          _out + ".STEP4_tmp")  # *.ped, *.dat (cf. 'java -jar' is included in 'BEAGLE2LINKAGE'.)
        # print(command)
        if not os.system(command):
            # Remove
            if not self.__save_intermediates:
                os.system('rm {}'.format(_out + ".STEP4_tmp.dat"))  # *.dat file is unnecessary.

        command = 'cut -d \' \' -f-5 {} > {}'.format(_out + ".STEP4_tmp.ped",
                                                     _out + ".STEP4_tmp.ped.left")  # ['FID', 'IID', 'PID', 'MID', 'Sex']
        # print(command)
        os.system(command)

        command = 'cut -d \' \' -f6- {} > {}'.format(_out + ".STEP4_tmp.ped",
                                                     _out + ".STEP4_tmp.ped.right")  # genotype information part.
        # print(command)
        os.system(command)

        command = 'paste -d \' -9 \' {} /dev/null /dev/null /dev/null {} > {}'.format(_out + ".STEP4_tmp.ped.left",
                                                                                      _out + ".STEP4_tmp.ped.right",
                                                                                      _out + ".ped")
        # print(command)
        if not os.system(command):
            # Remove
            if not self.__save_intermediates:
                os.system('rm {}'.format(_out + ".STEP4_tmp.ped"))
                os.system('rm {}'.format(_out + ".STEP4_tmp.ped.left"))
                os.system('rm {}'.format(_out + ".STEP4_tmp.ped.right"))

        # (1) rsid, (2) bp, (3) allele1
        os.system(' '.join(["cut -d \' \' -f1", _out + ".markers", ">", _out + ".STEP4_map.rsid"]))

        os.system(' '.join(["cut -d \' \' -f2", _out + ".markers", ">", _out + ".STEP4_map.bp"]))

        os.system(' '.join(["cut -d \' \' -f3", _out + ".markers", ">", _out + ".STEP4_map.allele1"]))

        os.system(' '.join(
            ["paste -d \'6  0 \'", "/dev/null", "/dev/null", _out + ".STEP4_map.rsid", "/dev/null", "/dev/null",
             _out + ".STEP4_map.bp", ">", _out + ".map"]))

        # os.system(' '.join(
        #     ["paste -d \'   \'", _out + ".STEP4_map.rsid", _out + ".STEP4_map.bp", ">", _out + ".refallele"]))

        os.system(' '.join(
            ["paste -d \' \'", _out + ".STEP4_map.rsid", _out + ".STEP4_map.allele1", ">",
             _out + ".refallele"]))

        """
        (2019. 07. 09.)
        To make '*.refallele' file, I think right part is supposed to be '_out + ".STEP4_map.allele1"' not '_out + ".STEP4_map.bp"'
        """

        # bed, bim, fam files.
        command = ' '.join([self.PLINK, '--ped {} --map {} --make-bed --reference-allele {} --out {}'.format(
            _out + ".ped",
            _out + ".map",
            _out + ".refallele",
            _out
        )])
        # print(command)
        if not os.system(command):
            # Remove
            if not self.__save_intermediates:
                os.system('rm {}'.format(_out + ".STEP4_map.rsid"))
                os.system('rm {}'.format(_out + ".STEP4_map.bp"))
                os.system('rm {}'.format(_out + ".STEP4_map.allele1"))
                os.system('rm {}'.format(_out + ".ped"))
                os.system('rm {}'.format(_out + ".map"))
                os.system('rm {}'.format(_out + ".log"))
                os.system('rm {}'.format(_out + ".refallele"))

        # Allele Frequency file(*.frq)
        command = ' '.join([self.PLINK, '--bfile {} --keep-allele-order --freq --out {}'.format(_out, _out + ".FRQ")])
        # print(command)
        if not os.system(command):
            # Remove
            if not self.__save_intermediates:
                os.system('rm {}'.format(_out + ".FRQ.log"))


        return _out




    def removeOUTPUT(self, _prefix):

        # Removing Exon_N reference panel.
        os.system('rm {}'.format(_prefix + '.bed'))
        os.system('rm {}'.format(_prefix + '.bim'))
        os.system('rm {}'.format(_prefix + '.fam'))
        os.system('rm {}'.format(_prefix + '.bgl.phased'))
        os.system('rm {}'.format(_prefix + '.markers'))
        os.system('rm {}'.format(_prefix + '.FRQ.frq'))

        return 0
