#-*- coding: utf-8 -*-

import os, sys, re


########## < Core Varialbes > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))



def Make_EXON234_Panel(infile, outfile, BEAGLE2LINKAGE, plink):


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
    os.system('rm {}'.format(OUTPUT_REF+".STEP1_SNP.markers"))
    os.system('rm {}'.format(OUTPUT_REF+".STEP1_class1_4dit.markers"))
    os.system('rm {}'.format(OUTPUT_REF+".STEP1_class2_4dit.markers"))



    # print("STEP2_EXON234_MARKERS")
    #
    # os.system(' '.join(["python HLA2EXON234.py", infile + ".STEP1_SNP_4dit.markers", infile + ".bgl.phased",
    #                     infile + ".STEP2_exon234.bgl.phased", infile + ".markers", infile + ".STEP2_exon234.markers"]))
    #
    # command = 'python redefineBPv1BH.py ./data/HM_CEU_REF_temp.markers ./data/HM_CEU_REF_temp_refined.markers'
    #
    #
    #
    # print("STEP3_SORT")
    #
    # os.system(
    #     ' '.join(["python redefineBPv1BH.py", infile + ".STEP2_exon234.markers", infile + ".STEP3_refined.markers"]))
    #
    # os.system(' '.join(["sort -gk 2", infile + ".STEP3_refined.markers", '>', outfile + ".markers"]))
    #
    # os.system(' '.join(
    #     ["python BGL2SortBGl.py", outfile + ".markers", infile + ".STEP2_exon234.bgl.phased", outfile + ".bgl.phased"]))
    #
    #
    #
    # print("STEP_4_Make_plink_file")
    #
    # os.system(' '.join(["cat", outfile + ".bgl.phased", "|", "java -jar", BEAGLE2LINKAGE, outfile + ".STEP4_tmp"]))
    #
    # os.system(' '.join(["cut -d \' \' -f-5", outfile + ".STEP4_tmp.ped", ">", outfile + ".STEP4_tmp.ped.left"]))
    #
    # os.system(' '.join(["cut -d \' \' -f6-", outfile + ".STEP4_tmp.ped", ">", outfile + ".STEP4_tmp.ped.right"]))
    #
    # os.system(' '.join(["paste -d \' -9 \' ", outfile + ".STEP4_tmp.ped.left", "/dev/null", "/dev/null", "/dev/null",
    #                     outfile + ".STEP4_tmp.ped.right", ">", outfile + ".ped"]))
    #
    # os.system(' '.join(["cut -d \' \' -f1", outfile + ".markers", ">", outfile + ".STEP4_map.rsid"]))
    #
    # os.system(' '.join(["cut -d \' \' -f2", outfile + ".markers", ">", outfile + ".STEP4_map.bp"]))
    #
    # os.system(' '.join(["cut -d \' \' -f3", outfile + ".markers", ">", outfile + ".STEP4_map.allele1"]))
    #
    # os.system(' '.join(
    #     ["paste -d \'6  0 \'", "/dev/null", "/dev/null", outfile + ".STEP4_map.rsid", "/dev/null", "/dev/null",
    #      outfile + ".STEP4_map.bp", ">", outfile + ".map"]))
    #
    # os.system(' '.join(
    #     ["paste -d \' \'", outfile + ".STEP4_map.rsid", outfile + ".STEP4_map.bp", ">", outfile + ".refallele"]))
    #
    # os.system(' '.join(
    #     [plink, "--noweb --allow-no-sex", "--ped", outfile + ".ped", "--map", outfile + ".map", "--make-bed",
    #      "--reference-allele", outfile + ".refallele", "--out", outfile]))
    #
    # os.system(' '.join([plink, "--noweb --allow-no-sex", "--bfile", outfile, " --keep-allele-order", "--freq", "--out",
    #                     outfile + ".FRQ"]))

    return 0



if __name__ == '__main__':

    """
    < Make_EXON234_Panel.py >
    
    INPUT : (1) Prefix of Reference Panel, (2) Prefix of Output file, (3) Path to 'beagle2linkage.jar' file,
            (4) Path to Plink(v1.07) file.
            
    OUTPUT : Three copies of the reference panel of which the HLA markers have genomic position of middle point of Exon 2,3,4.
    
    """

    # [infile, outfile, BEAGLE2LINKAGE, plink] = sys.argv[1:5]

    ### Temporary Hardcoding
    infile = '/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/CookHLA/data/HLA_PANEL/T1DGC/T1DGC_REF'
    outfile = '/Users/wansun/Git_Projects/CookHLA/tests/_3_CookHLA/20190521_EXON234/Exon234_Panel_test'
    BEAGLE2LINKAGE = '/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/dependency/beagle2linkage.jar'
    plink = '/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/dependency/plink107/osx/plink'


    Make_EXON234_Panel(infile, outfile, BEAGLE2LINKAGE, plink)