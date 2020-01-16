#-*- coding: utf-8 -*-

import os, sys, re
import argparse, textwrap

import pandas as pd


std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]


p_header = re.compile(r'^FID\s+IID')
p_1stTwo = re.compile(r'^(\S+)\s+(\S+)\s+')




def SubsetBGLPhased(_bgl, _out=None, _toKeep=None, _toRemove=None, _toExtract=None, _toExclude=None):

    if not (bool(_toKeep) or bool(_toRemove) or bool(_toExtract) or bool(_toExclude)):
        print("Nothing to subset. (Samples or Markers not given.)")
        return -1



    if bool(_toKeep) and bool(_toRemove):
        print(std_ERROR_MAIN_PROCESS_NAME + "Keeping and Removing samples can't be done simultaneously.")
        return -1

    elif not (bool(_toKeep) or bool(_toRemove)):

        sr_Samples = None
        pass

    else:

        ## Sample list

        KeepOrRemove = _toKeep if bool(_toKeep) else _toRemove

        f_sample = open(KeepOrRemove, 'r')
        line1 = f_sample.readline()

        if bool(p_header.match(line1)):
            sr_Samples = pd.read_csv(KeepOrRemove, sep='\s+', header=0, dtype=str).loc[:, 'IID'] # Intentionally use .loc[] to make it 'Series'.
        else:
            sr_Samples = pd.read_csv(KeepOrRemove, sep='\s+', dtype=str, names=['FID', 'IID']).loc[:, 'IID']

        f_sample.close()
        # print(sr_Samples)



    if bool(_toExtract) and bool(_toExclude):
        print(std_ERROR_MAIN_PROCESS_NAME + "Extracting and Excluding markers can't be done simultaneously.")
        return -1

    elif not (bool(_toExtract) or bool(_toExclude)):

        sr_Markers = None
        pass

    else:

        ## Marker list
        ExtractOrExclude = _toExtract if bool(_toExtract) else _toExclude
        sr_Markers = pd.read_csv(ExtractOrExclude, header=None, dtype=str).iloc[:, 0] # Intentionally use .loc[] to make it 'Series'.

        # print(sr_Markers)



    ### Find the index of 'I id ...' line (Individual ID line).
    with open(_bgl, 'r') as f_BGL:
        count = 0
        idx_I = 0
        idx_EndOfHeader = 0 # This will be used as `idx_EndOfHeader + 1` in the 'Extract' function.

        for line in f_BGL:

            m = p_1stTwo.match(line)

            if m.group(1) == 'I':
                idx_I = count
            elif m.group(1) == 'M':
                idx_EndOfHeader = count
                break
            else:
                count += 1

    # print('"I id" line is : {}'.format(idx_I))
    # print('End of Header : {}'.format(idx_EndOfHeader))

    df_BGL = pd.read_csv(_bgl, sep='\s+', header=None, dtype=str)
    # print("df_BGL:\n{}\n".format(df_BGL.head(20)))




    ########## < Main subsetting > ##########

    ### Samples in columns

    if bool(_toKeep) != bool(_toRemove):

        if bool(_toKeep):

            ## Keeping given samples AND the 1st two columns.

            f_samples = df_BGL.iloc[idx_I, :].isin(sr_Samples)

            # Including first 2 columns(i.e. meta info.)
            f_samples.iat[0] = True
            f_samples.iat[1] = True

        else:

            ## Removing given samples

            f_samples = ~df_BGL.iloc[idx_I, :].isin(sr_Samples)  # The 1st two columns will be included.

    else:

        f_samples = pd.Series([True for z in range(df_BGL.shape[1])], index=df_BGL.columns)




    ### Markers in rows

    if bool(_toExtract) != bool(_toExclude):

        if bool(_toExtract):

            ## Extracting given markers AND the Header rows.

            f_markers = df_BGL.iloc[:, 1].isin(sr_Markers)

            for i in range(idx_EndOfHeader + 1):
                f_markers.iat[i] = True

        else:

            ## Excluding given markers
            f_markers = ~df_BGL.iloc[:, 1].isin(sr_Markers)  # Header rows will be included.

    else:

        f_markers = pd.Series([True for z in range(df_BGL.shape[0])], index=df_BGL.index)


    df_RETURN = df_BGL.loc[f_markers, f_samples]



    if bool(_out):
        df_RETURN.to_csv(_out, sep=' ', header=False, index=False)
        return _out
    else:
        return df_RETURN



if __name__ == "__main__":

    ########## < Main parser > ##########

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    ###########################################################################################

        SubsetBGLPhased.py
        
        - Subset (1) samples in columns or (2) markers in rows
        - Argument '--keep/--remove' and '--extract/--exclude' is motivated by PLINK.


    ###########################################################################################
                                     '''),
                                     add_help=False
                                     )

    ### Common arguments to share over the modules.

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="Show this help message and exit\n\n", action='help')

    parser.add_argument("--bgl", help="\nBeagle file.\n\n", required=True)
    parser.add_argument("--out", "-o", help="\nOutput file name prefix\n\n", required=True)

    samples = parser.add_mutually_exclusive_group()
    samples.add_argument("--keep", help="\nSample list file to keep.\n\n")
    samples.add_argument("--remove", help="\nSample list file to remove.\n\n")

    markers = parser.add_mutually_exclusive_group()
    markers.add_argument("--extract", help="\nMarker list file to extract.\n\n")
    markers.add_argument("--exclude", help="\nMarker list file to exclude.\n\n")



    ##### < for Testing > #####

    # args = parser.parse_args(["--bgl", "/media/sf_VM_Shared/Projects/CookHLA/tests/T1DGC/T1DGC_REF.bgl.phased",
    #                           "--out", "/media/sf_VM_Shared/Projects/CookHLA/tests/T1DGC/T1DGC_REF.subsetsubset.bgl.phased",
    #                           "--remove", "/media/sf_VM_Shared/Projects/CookHLA/tests/T1DGC_CookQC/WronglyPhased.samples.txt"
    #                           ])

    # args = parser.parse_args(["--bgl", "/media/sf_VM_Shared/Projects/CookHLA/tests/T1DGC_CookQC/T1DGC_REF.ONLY_Variants_HLA.bgl.phased",
    #                           "--out", "/media/sf_VM_Shared/Projects/CookHLA/tests/T1DGC/T1DGC_REF.subsetsubset.bgl.phased",
    #                           "--extract", "/media/sf_VM_Shared/Projects/CookHLA/tests/T1DGC_CookQC/ToExclude.MKref.txt"
    #                           ])


    ##### < for Publish > #####
    args = parser.parse_args()
    print(args)


    SubsetBGLPhased(args.bgl, args.out, args.keep, args.remove, args.extract, args.exclude)