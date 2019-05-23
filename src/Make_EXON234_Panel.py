#-*- coding: utf-8 -*-

import os, sys, re
from src.redefineBPv1BH import redefineBP
from src.BGL2SortBGl import BGL2SortBGL_WS

########## < Core Varialbes > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))



def Make_EXON234_Panel(infile, outfile, BEAGLE2LINKAGE, PLINK, __save_intermediates=False):


    REF_base = os.path.basename(infile)
    OUTPUT_dir = os.path.dirname(outfile)
    OUTPUT_REF = os.path.join(OUTPUT_dir, REF_base)



    print("STEP1_Collect_SNP_HLA_DATA")

    # In STEP1, New *.markers file will be used just next step.
    command = "grep rs {} > {}".format(infile + ".markers", OUTPUT_REF+".STEP1_SNP.markers")
    # print(command)
    os.system(command)

    command = "grep \'HLA_[A-Z]_[0-9][0-9][0-9][0-9]\' {} > {}".format(infile + ".markers", OUTPUT_REF+".STEP1_class1_4dit.markers")
    # print(command)
    os.system(command)

    command = "grep \'HLA_[A-Z][A-Z][A-Z][0-9]_[0-9][0-9][0-9][0-9]\' {} > {}".format(infile + ".markers", OUTPUT_REF+".STEP1_class2_4dit.markers")
    # print(command)
    os.system(command)

    command = 'cat {} {} {} > {}'.format(OUTPUT_REF+".STEP1_SNP.markers", OUTPUT_REF+".STEP1_class1_4dit.markers",
                                         OUTPUT_REF+".STEP1_class2_4dit.markers", OUTPUT_REF+".STEP1_SNP_4dit.markers")
    # print(command)
    os.system(command)


    # Remove
    if not __save_intermediates:
        os.system('rm {}'.format(OUTPUT_REF+".STEP1_SNP.markers"))
        os.system('rm {}'.format(OUTPUT_REF+".STEP1_class1_4dit.markers"))
        os.system('rm {}'.format(OUTPUT_REF+".STEP1_class2_4dit.markers"))



    print("STEP2_EXON234_MARKERS")

    [outbgl, outmarker] = HLA2EXON234(OUTPUT_REF+".STEP1_SNP_4dit.markers",
                                      infile + ".bgl.phased", OUTPUT_REF+".STEP2_exon234.bgl.phased",
                                      infile + ".markers", OUTPUT_REF+".STEP2_exon234.markers")

    # Remove
    if not __save_intermediates:
        os.system('rm {}'.format(OUTPUT_REF+".STEP1_SNP_4dit.markers"))


    print("STEP3_SORT")

    # Dispersing genomic positions of given marker file (*.markers)
    refiend_outmarker = redefineBP(outmarker, OUTPUT_REF+".STEP3_refined.markers")
    # print(refiend_outmarker)


    # Sorting the dispersed marker file.
    command = 'sort -gk 2 {} > {}'.format(refiend_outmarker, outfile+'.markers')
    # print(command)
    if not os.system(command):
        # Remove
        if not __save_intermediates:
            os.system('rm {}'.format(outmarker))
            os.system('rm {}'.format(refiend_outmarker))


    # Sorting beagle file to the oreder of the above sorted markers file
    sorted_outbgl = BGL2SortBGL_WS(outfile+'.markers', outbgl, outfile + ".bgl.phased")
    # print(sorted_outbgl)
    if not os.path.exists(sorted_outbgl):
        print(std_ERROR_MAIN_PROCESS_NAME + "Failed to generate '{}'.".format(sorted_outbgl))
        sys.exit()
    else:
        # Remove
        if not __save_intermediates:
            os.system('rm {}'.format(outbgl))



    print("STEP_4_Make_plink_file")

    command = 'cat {} | {} {}'.format(sorted_outbgl, BEAGLE2LINKAGE, outfile + ".STEP4_tmp") # *.ped, *.dat (cf. 'java -jar' is included in 'BEAGLE2LINKAGE'.)
    # print(command)
    if not os.system(command):
        # Remove
        if not __save_intermediates:
            os.system('rm {}'.format(outfile + ".STEP4_tmp.dat")) # *.dat file is unnecessary.


    command = 'cut -d \' \' -f-5 {} > {}'.format(outfile + ".STEP4_tmp.ped", outfile + ".STEP4_tmp.ped.left") # ['FID', 'IID', 'PID', 'MID', 'Sex']
    # print(command)
    os.system(command)

    command = 'cut -d \' \' -f6- {} > {}'.format(outfile + ".STEP4_tmp.ped", outfile + ".STEP4_tmp.ped.right") # genotype information part.
    # print(command)
    os.system(command)


    command = 'paste -d \' -9 \' {} /dev/null /dev/null /dev/null {} > {}'.format(outfile + ".STEP4_tmp.ped.left", outfile + ".STEP4_tmp.ped.right", outfile + ".ped")
    # print(command)
    if not os.system(command):
        # Remove
        if not __save_intermediates:
            os.system('rm {}'.format(outfile + ".STEP4_tmp.ped"))
            os.system('rm {}'.format(outfile + ".STEP4_tmp.ped.left"))
            os.system('rm {}'.format(outfile + ".STEP4_tmp.ped.right"))


    # (1) rsid, (2) bp, (3) allele1
    os.system(' '.join(["cut -d \' \' -f1", outfile + ".markers", ">", outfile + ".STEP4_map.rsid"]))

    os.system(' '.join(["cut -d \' \' -f2", outfile + ".markers", ">", outfile + ".STEP4_map.bp"]))

    os.system(' '.join(["cut -d \' \' -f3", outfile + ".markers", ">", outfile + ".STEP4_map.allele1"]))


    os.system(' '.join(
        ["paste -d \'6  0 \'", "/dev/null", "/dev/null", outfile + ".STEP4_map.rsid", "/dev/null", "/dev/null",
         outfile + ".STEP4_map.bp", ">", outfile + ".map"]))

    os.system(' '.join(
        ["paste -d \' \'", outfile + ".STEP4_map.rsid", outfile + ".STEP4_map.bp", ">", outfile + ".refallele"]))


    # bed, bim, fam files.
    command = ' '.join([PLINK, '--ped {} --map {} --make-bed --reference-allele {} --out {}'.format(
        outfile + ".ped",
        outfile + ".map",
        outfile + ".refallele",
        outfile
    )])
    # print(command)
    if not os.system(command):
        # Remove
        if not __save_intermediates:
            os.system('rm {}'.format(outfile + ".STEP4_map.rsid"))
            os.system('rm {}'.format(outfile + ".STEP4_map.bp"))
            os.system('rm {}'.format(outfile + ".STEP4_map.allele1"))
            os.system('rm {}'.format(outfile + ".{ped,map}"))
            os.system('rm {}'.format(outfile + ".log"))
            os.system('rm {}'.format(outfile + ".refallele"))


    # Allele Frequency file(*.frq)
    command = ' '.join([PLINK, '--bfile {} --keep-allele-order --freq --out {}'.format(outfile, outfile + ".FRQ")])
    # print(command)
    if not os.system(command):
        # Remove
        if not __save_intermediates:
            os.system('rm {}'.format(outfile + ".FRQ.log"))


    return outfile



def HLA2EXON234(markerchoicefile, inbgl, outbgl, inmarker, outmarker):

    """
    Originally, this function was 'HLA2EXON234.py' source file.

    """


    selectMarkers = {}

    HLA_EXON2_POSITION = [['HLA_A', '30018647'], ['HLA_C', '31347489'], ['HLA_B', '31432578'], ['HLA_DRB1', '32659998'], ['HLA_DQA1', '32717189'], ['HLA_DQB1', '32740687'], ['HLA_DPA1', '33145518'], ['HLA_DPB1', '33156558']]
    HLA_EXON3_POSITION = [['HLA_A', '30019161'], ['HLA_C', '31346966'], ['HLA_B', '31432060'], ['HLA_DRB1', '32657452'], ['HLA_DQA1', '32717867'], ['HLA_DQB1', '32737862'], ['HLA_DPA1', '33144914'], ['HLA_DPB1', '33160845']]
    HLA_EXON4_POSITION = [['HLA_A', '30020015'], ['HLA_C', '31346103'], ['HLA_B', '31431210']]



    with open(markerchoicefile) as fin:
        for l in fin:
            c = l.split()
            selectMarkers[c[0]] = True

    with open(inbgl) as pf, open(outbgl, 'w') as of:
        for l in pf:
            c = l.split()
            header = c[:2]
            data = c[2:]
            if header[0] == "M" and header[1] not in selectMarkers:
                continue

            if "HLA" in l:
                header[1] = header[1] + '_exon2'
                newdata = []
                for j in range(len(data)):
                    newdata.append(data[j])
                of.write(' '.join(header + newdata) + '\n')
                header[1] = (header[1].replace("_exon2", ""))
                header[1] = header[1] + '_exon3'
                newdata = []
                for j in range(len(data)):
                    newdata.append(data[j])
                of.write(' '.join(header + newdata) + '\n')
                header[1] = (header[1].replace("_exon3", ""))

                if "DRB1" in l:
                    continue
                if "DQA1" in l:
                    continue
                if "DQB1" in l:
                    continue
                if "DPA1" in l:
                    continue
                if "DPB1" in l:
                    continue

                header[1] = header[1] + '_exon4'
                newdata = []
                for j in range(len(data)):
                    newdata.append(data[j])
                of.write(' '.join(header + newdata) + '\n')
                header[1] = (header[1].replace("_exon4", ""))

            if "HLA" in l:
                continue

            newdata = []
            for j in range(len(data)):
                newdata.append(data[j])
            of.write(' '.join(header + newdata) + '\n')



    with open(inmarker) as mf, open(outmarker, 'w') as of:
        for l in mf:
            c = l.split()
            rsid_bp = c[:2]
            allele = c[2:]

            if rsid_bp[0] not in selectMarkers:
                continue

            if "HLA" in l:
                rsid_bp[0] = rsid_bp[0] + '_exon2'
                for i in range(8):
                    if HLA_EXON2_POSITION[i][0] in rsid_bp[0]:
                        rsid_bp[1] = HLA_EXON2_POSITION[i][1]
                        of.write(' '.join(rsid_bp + allele) + '\n')
                        rsid_bp[0] = (rsid_bp[0].replace("_exon2", ""))
                    if HLA_EXON2_POSITION[i][0] in rsid_bp[0]:
                        continue

                rsid_bp[0] = rsid_bp[0] + '_exon3'
                for i in range(8):
                    if HLA_EXON3_POSITION[i][0] in rsid_bp[0]:
                        rsid_bp[1] = HLA_EXON3_POSITION[i][1]
                        of.write(' '.join(rsid_bp + allele) + '\n')
                        rsid_bp[0] = (rsid_bp[0].replace("_exon3", ""))
                    if HLA_EXON3_POSITION[i][0] in rsid_bp[0]:
                        continue

                rsid_bp[0] = rsid_bp[0] + '_exon4'

                for i in range(3):
                    if HLA_EXON4_POSITION[i][0] in rsid_bp[0]:
                        rsid_bp[1] = HLA_EXON4_POSITION[i][1]

                        of.write(' '.join(rsid_bp + allele) + '\n')
                        rsid_bp[0] = (rsid_bp[0].replace("_exon4", ""))
                    if HLA_EXON4_POSITION[i][0] in rsid_bp[0]:
                        continue

            if "HLA" in l:
                continue

            of.write(' '.join(rsid_bp + allele) + '\n')


    return [outbgl, outmarker]





if __name__ == '__main__':

    """
    < Make_EXON234_Panel.py >
    
    INPUT : (1) Prefix of Reference Panel, (2) Prefix of Output file, (3) Path to 'beagle2linkage.jar' file,
            (4) Path to Plink(v1.07) file.
            
    OUTPUT : Three copies of the reference panel of which the HLA markers have genomic position of middle point of Exon 2,3,4.
    
    """

    [infile, outfile, BEAGLE2LINKAGE, PLINK] = sys.argv[1:5]

    # ### Temporary Hardcoding
    # infile = '/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/CookHLA/data/HLA_PANEL/T1DGC/T1DGC_REF'
    # outfile = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190521_EXON234/Exon234_Panel_test2'
    #
    # p_beagle2linkage = '/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/dependency/beagle2linkage.jar'
    # BEAGLE2LINKAGE = 'java -jar '+p_beagle2linkage
    #
    # p_plink = '/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/dependency/plink107/osx/plink'
    # PLINK = ' '.join([p_plink, '--noweb --allow-no-sex'])

    Make_EXON234_Panel(infile, outfile, BEAGLE2LINKAGE, PLINK)