#-*- coding: utf-8 -*-

import os, sys, re
from src.redefineBPv1BH import redefineBP
from src.BGL2SortBGl import BGL2SortBGL_WS

########## < Core Varialbes > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
isClassI = {"A": True, "B": True, "C": True, "DPA1": False, "DPB1": False, "DQA1": False, "DQB1": False, "DRB1": False}

# Patterns to use.
p_HLA_2field = re.compile(r'^HLA_(\w+)_\d{4}')
# p = re.compile(r'^([A-Za-z0-9_-]+)\s+(\w+)') # Frist two columns (ex. 'P pedigree' or 'rs969931 29602876', ... )
p = re.compile(r'^(\S+)\s+(\S+)') # Frist two columns (ex. 'P pedigree' or 'rs969931 29602876', ... )


def Make_EXON234_Panel(__exonN__, infile, outfile, BEAGLE2LINKAGE, PLINK, __save_intermediates=False,
                       f_only_HLA=True):


    if __exonN__ not in ['exon2', 'exon3', 'exon4']:
        print(std_ERROR_MAIN_PROCESS_NAME + "Exon error.")
        sys.exit()


    REF_base = os.path.basename(outfile)
    OUTPUT_dir = os.path.dirname(outfile)
    OUTPUT_REF = os.path.join(OUTPUT_dir, REF_base)


    print(std_MAIN_PROCESS_NAME + "Generating Reference Panel for {}".format(__exonN__.capitalize()))

    if f_only_HLA:

        # print("STEP1_Collect_SNP_HLA_DATA")

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

    else:

        # In case there is no need to filter out only rs_id alleles and HLA alleles.

        command = 'cp {} {}'.format(infile + ".markers", OUTPUT_REF+".STEP1_SNP_4dit.markers")
        print(command)
        os.system(command)



    # print("STEP2_EXON234_MARKERS")

    # [outbgl, outmarker] = HLA2EXON234(OUTPUT_REF+".STEP1_SNP_4dit.markers",
    #                                   infile + ".bgl.phased", OUTPUT_REF+".STEP2_exon234.bgl.phased",
    #                                   infile + ".markers", OUTPUT_REF+".STEP2_exon234.markers")

    [outbgl, outmarker] = MakeExon234(__exonN__, OUTPUT_REF+".STEP1_SNP_4dit.markers",
                                      infile + ".bgl.phased", OUTPUT_REF+".STEP2_exon234.bgl.phased",
                                      infile + ".markers", OUTPUT_REF+".STEP2_exon234.markers")

    # Remove
    if not __save_intermediates:
        os.system('rm {}'.format(OUTPUT_REF+".STEP1_SNP_4dit.markers"))


    # print("STEP3_SORT")

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



    # print("STEP_4_Make_plink_file")

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
                # 이 블럭에서 2-digit HLA allele들이 걸려나가줌
                continue

            if "HLA" in l:

                # Exon 2
                header[1] = header[1] + '_exon2'
                newdata = []
                for j in range(len(data)):
                    newdata.append(data[j])
                of.write(' '.join(header + newdata) + '\n')

                # Exon 3
                header[1] = (header[1].replace("_exon2", ""))
                header[1] = header[1] + '_exon3'
                newdata = []
                for j in range(len(data)):
                    newdata.append(data[j])
                of.write(' '.join(header + newdata) + '\n')

                # Exon 4
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

            # "HLA_A_0101" 이런애들은 바로 위 if구문에서 마무리가 됨.
            # 여기는 사실상 "M rs1234" 이런 rs_id를 가지는 SNP들이 여기서 처리됨.
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
                continue # 뭐던동 marker의 label이 'HLA_A_0101' 이런식이면 여기서 마무리가 됨.

            of.write(' '.join(rsid_bp + allele) + '\n') # 여기도 얘가 rs_id가지는 marker들을 담당하는 부분.


    return [outbgl, outmarker]



def MakeExon234(__exonN__, markerchoicefile, inbgl, outbgl, inmarker, outmarker):



    HLA_POSITION = {
        'exon2': {'A': '30018647', 'C': '31347489', 'B': '31432578', 'DRB1': '32659998', 'DQA1': '32717189', 'DQB1': '32740687', 'DPA1': '33145518', 'DPB1': '33156558'},
        'exon3': {'A': '30019161', 'C': '31346966', 'B': '31432060', 'DRB1': '32657452', 'DQA1': '32717867', 'DQB1': '32737862', 'DPA1': '33144914', 'DPB1': '33160845'},
        'exon4': {'A': '30020015', 'C': '31346103', 'B': '31431210'}
    }



    ### Obtaining marker labels from sorted *.markers file(`sort_markers`)

    selectMarkers = []
    with open(markerchoicefile) as fin:
        selectMarkers = [p.match(string=l).group(1) for l in fin]

    # print(selectMarkers)


    ## Changing rows for exon N in beagle file.

    with open(inbgl) as pf, open(outbgl, 'w') as of:
        for l in pf:

            m = p.match(l)
            header = [m.group(1), m.group(2)]

            if header[0] == "M" and header[1] not in selectMarkers:
                # Excluding out 2-digit HLA alleles or Bizarre markers.
                continue


            m2 = p_HLA_2field.match(header[1])

            if m2:

                hla = m2.group(1)

                if __exonN__ == 'exon4' and not isClassI[hla]:
                    # new_line = l
                    continue
                else:
                    # HLA alleles (with 4-digit(2-field))
                    header2 = ' '.join([header[0], header[1]+'_{}'.format(__exonN__)])
                    new_line = p.sub(repl=header2, string=l)
                    # body = p.sub(repl='', string=l)
                    # new_line = header2 + body

                of.write(new_line)

            else:
                # Normal SNPs (ex. rs969931)
                of.write(l)


    ### Changin rows for exon N in markers file(*.markers)
    with open(inmarker) as mf, open(outmarker, 'w') as of:
        for l in mf:

            m = p.match(l)
            rsid_bp = [m.group(1), m.group(2)]
            # print(rsid_bp)

            if rsid_bp[0] not in selectMarkers:
                continue


            m2 = p_HLA_2field.match(rsid_bp[0])

            if m2:
                hla = m2.group(1)

                if __exonN__ == 'exon4' and not isClassI[hla]:
                    continue
                else:
                    new_marker_name = rsid_bp[0]+'_{}'.format(__exonN__)
                    new_bp = HLA_POSITION[__exonN__][hla]

                    new_line = p.sub(repl=' '.join([new_marker_name, new_bp]), string=l)
                    # print("Before : {}".format(l))
                    # print("After : {}".format(new_line))

                    of.write(new_line)
            else:
                of.write(l)


    return [outbgl, outmarker]


if __name__ == '__main__':

    """
    < Make_EXON234_Panel.py >
    
    INPUT : (0) 'exon2', 'exon3' or 'exon4', (1) Prefix of Reference Panel, (2) Prefix of Output file, 
            (3) Path to 'beagle2linkage.jar' file, (4) Path to Plink(v1.07) file.
            
    OUTPUT : Three copies of the reference panel of which the HLA markers have genomic position of middle point of Exon 2,3,4.
    
    """

    [__exonN__, infile, outfile, BEAGLE2LINKAGE, PLINK] = sys.argv[1:6]

    # ### Temporary Hardcoding
    # infile = '/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/CookHLA/data/HLA_PANEL/T1DGC/T1DGC_REF'
    # outfile = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190521_EXON234/T1DGC_REF_exon4'
    #
    # p_beagle2linkage = '/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/dependency/beagle2linkage.jar'
    # BEAGLE2LINKAGE = 'java -jar '+p_beagle2linkage
    #
    # p_plink = '/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/dependency/plink107/osx/plink'
    # PLINK = ' '.join([p_plink, '--noweb --allow-no-sex'])

    Make_EXON234_Panel(__exonN__, infile, outfile, BEAGLE2LINKAGE, PLINK)