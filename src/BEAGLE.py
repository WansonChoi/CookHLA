#-*- coding: utf-8 -*-

import os, sys, re
# import subprocess
# from shutil import which
# import random
import pandas as pd


std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]


p_header = re.compile(r'^FID\s+IID')
p_1stTwo = re.compile(r'^(\S+)\s+(\S+)\s+')



class BEAGLE():

    def __init__(self, _bgl, _markers=None, _idx_I=None, _idx_EndOfHeader=None, _Label='first'):

        print("New {} beagle object is created.".format(_Label))

        ### Class variables
        self.bgl = _bgl
        self.markers = _markers
        self.idx_I = _idx_I
        self.idx_EndOfHeader = _idx_EndOfHeader

        self.Label = _Label

        if isinstance(self.bgl, str) and os.path.exists(self.bgl):

            # ### Phased Beagle file.
            # if os.path.exists(_bgl + '.bgl.phased'):
            #     self.bgl = _bgl + '.bgl.phased'
            # else:
            #     print(std_ERROR_MAIN_PROCESS_NAME + "Given phased beagle file('{}') doesn't exist.".format(
            #         _bgl + '.bgl.phased'))
            #     sys.exit()

            # ### Beagle Markers file
            # if os.path.exists(_bgl + '.markers'):
            #     self.markers = _bgl + '.markers'
            # else:
            #     print(std_ERROR_MAIN_PROCESS_NAME + "Given beagle markers file('{}') doesn't exist.".format(_bgl + '.markers'))
            #     sys.exit()

            ### Find the index of 'I id ...' line (Individual ID line).
            with open(self.bgl, 'r') as f_BGL_PHASED:

                count = 0

                for line in f_BGL_PHASED:

                    m = p_1stTwo.match(line)

                    if m.group(1) == 'I':
                        self.idx_I = count
                    elif m.group(1) == 'M':
                        self.idx_EndOfHeader = count
                        break
                    else:
                        count += 1

            print('"I id" line is : {}'.format(self.idx_I))
            print('End of Header : {}'.format(self.idx_EndOfHeader))

            self.bgl = pd.read_csv(self.bgl, sep='\s+', header=None, dtype=str)
            print("df_BGL:\n{}\n".format(self.bgl.head(20)))



        elif isinstance(_bgl, pd.DataFrame):

            pass

    def __del__(self):
        print("Beagle {} obj died.".format(self.Label))



    def SubsetBGLPhased(self, _out=None, _toKeep=None, _toRemove=None, _toExtract=None, _toExclude=None):


        if not (bool(_toKeep) or bool(_toRemove) or bool(_toExtract) or bool(_toExclude)):
            print("Nothing to subset. (Samples or Markers not given.)")
            return -1


        if bool(_toKeep) and bool(_toRemove):
            print(std_ERROR_MAIN_PROCESS_NAME + "Keeping and Removing samples can't be done simultaneously.")
            return -1

        elif not (bool(_toKeep) or bool(_toRemove)):

            sr_Samples = None

        else:

            ## Sample list

            KeepOrRemove = _toKeep if bool(_toKeep) else _toRemove

            f_sample = open(KeepOrRemove, 'r')
            line1 = f_sample.readline()

            if bool(p_header.match(line1)):
                sr_Samples = pd.read_csv(KeepOrRemove, sep='\s+', header=0, dtype=str) \
                                .loc[:, 'IID']  # Intentionally use .loc[] to make it as 'Series'.
            else:
                sr_Samples = pd.read_csv(KeepOrRemove, sep='\s+', dtype=str, names=['FID', 'IID']) \
                                 .loc[:, 'IID']

            f_sample.close()
            # print(sr_Samples)


        if bool(_toExtract) and bool(_toExclude):
            print(std_ERROR_MAIN_PROCESS_NAME + "Extracting and Excluding markers can't be done simultaneously.")
            return -1

        elif not (bool(_toExtract) or bool(_toExclude)):

            sr_Markers = None

        else:

            ## Marker list
            ExtractOrExclude = _toExtract if bool(_toExtract) else _toExclude
            sr_Markers = pd.read_csv(ExtractOrExclude, header=None, dtype=str) \
                             .iloc[:, 0]  # Intentionally use .loc[] to make it 'Series'.



        ########## < Main subsetting > ##########

        ### Samples in columns

        if bool(_toKeep) != bool(_toRemove):

            if bool(_toKeep):

                ## Keeping given samples AND the 1st two columns.

                f_samples = self.bgl.iloc[self.idx_I, :].isin(sr_Samples)

                # Including first 2 columns(i.e. meta info.)
                f_samples.iat[0] = True
                f_samples.iat[1] = True

            else:

                ## Removing given samples

                f_samples = ~self.bgl.iloc[self.idx_I, :].isin(sr_Samples)  # The 1st two columns will be included.

        else:

            f_samples = pd.Series([True for z in range(self.bgl.shape[1])], index=self.bgl.columns)




        ### Markers in rows

        if bool(_toExtract) != bool(_toExclude):

            if bool(_toExtract):

                ## Extracting given markers AND the Header rows.

                f_markers = self.bgl.iloc[:, 1].isin(sr_Markers)

                for i in range(self.idx_EndOfHeader + 1):
                    f_markers.iat[i] = True

            else:

                ## Excluding given markers
                f_markers = ~self.bgl.iloc[:, 1].isin(sr_Markers)  # Header rows will be included.

        else:

            f_markers = pd.Series([True for z in range(self.bgl.shape[0])], index=self.bgl.index)


        df_RETURN = self.bgl.loc[f_markers, f_samples]



        if bool(_out):
            df_RETURN.to_csv(_out, sep=' ', header=False, index=False)
            return _out
        else:
            return df_RETURN.copy()



    def DoublingBGLPhased(self, _out=None, _Label=None):

        # IID = self.df_bgl_phased.iloc[self.idx_I, :]

        l_temp = []

        count = 0

        for col in self.bgl.iteritems():

            if count >= 2:

                sr_temp = col[1].copy()

                if not count % 2:
                    # chr1
                    sr_temp.iat[self.idx_I] = sr_temp.iat[self.idx_I] + '_1'
                else:
                    # chr2
                    sr_temp.iat[self.idx_I] = sr_temp.iat[self.idx_I] + '_2'

                # print(sr_temp.head())

                l_temp.append(sr_temp)
                l_temp.append(sr_temp)



            else:
                l_temp.append(col[1])

            count += 1
            # if count > 5: break;


        df_RETURN = pd.concat(l_temp, axis=1)
        # print("Doubled bgl phased file :\n{}\n".format(df_RETURN.iloc[:10, 2:10]))


        if bool(_out):
            df_RETURN.to_csv(_out, sep=' ', header=False, index=False)
            # return _out


        return BEAGLE(df_RETURN, _idx_I=self.idx_I, _idx_EndOfHeader=self.idx_EndOfHeader, _Label=_Label)




if __name__ == "__main__":

    _bgl = '/media/sf_VM_Shared/Projects/CookHLA/example/reference_panel/Pan-Asian_REF.100.example.bgl.phased'

    Pan_Asian_bgl_phased = BEAGLE(_bgl)
    print(Pan_Asian_bgl_phased.bgl)

    df_temp = Pan_Asian_bgl_phased.DoublingBGLPhased(_Label='Double').DoublingBGLPhased(_Label='DoubleDouble')
    print(df_temp.bgl.iloc[:, :5])