# CookHLA

## (1) Introduction

CookHLA imputes HLA types of target patient genotype data based on a reference panel information.

[Improvement]

1. Local embedding
2. Adative Genetic Map

<br>
<br>

## (2) Requirements and Dependencies

Prepare OS X(Mac) or Linux operating system. CookHLA currently doesn't support Windows. 

> CookHLA were tested in the next specific operating systems. Linux: Ubuntu 19.04(Disco Dingo), Ubuntu 18.04.3 LTS(Bionic Beaver), CentOS_7, Linux Mint 19.2 Cinnamon(Tina). OS X : Catalina, Mojave

<!-- > CentOS 7 or Ubuntu is most recommended in case of using Linux. -->

> In case of using Catalina OS X, **Using default shell as 'Bash($)' not 'Zsh(%)'** is recommended. To change the default shell to Bash, Please reference this blog(https://www.howtogeek.com/444596/how-to-change-the-default-shell-to-bash-in-macos-catalina/).

<br>

Then, Download this project in somewhere directory of your OS X or Linux system. It will be assumed that 'git' command is already installed in your system.

```
$ git clone https://github.com/WansonChoi/CookHLA.git
$ cd CookHLA
```
<br>


The following requirements must be installed in your system.

- python>=3.6.x
- pandas>=0.25.3
- perl>=5.26.2
- R>=3.6.x
- Java
- csh(or tcsh)

Next external software must be prepared in 'dependency/' folder. Users have to do it themselves because of copyright issue.

- PLINK1.9b (https://www.cog-genomics.org/plink/)
- mach1 (http://csg.sph.umich.edu/abecasis/MACH/download/)
- beagle.27Jan18.7e1.jar (https://faculty.washington.edu/browning/beagle/b4_1.html - 'Download Beagle 4.1' section; **Rename it to 'Beagle4.jar'** after downloading it in 'dependency/' folder.)
- beagle.18May20.d20.jar (http://faculty.washington.edu/browning/beagle/beagle.html - 'Download Beagle 5.1' section; **Rename it to 'Beagle5.jar'** after downloading it in 'dependency/' folder.)
- beagle2linkage.jar (https://faculty.washington.edu/browning/beagle_utilities/utilities.html - 'File conversion utilities' section)
- beagle2vcf.jar (https://faculty.washington.edu/browning/beagle_utilities/utilities.html - 'File conversion utilities' section)
- linkage2beagle.jar (https://faculty.washington.edu/browning/beagle_utilities/utilities.html - 'File conversion utilities' section)
- vcf2beagle.jar (https://faculty.washington.edu/browning/beagle_utilities/utilities.html - 'File conversion utilities' section)
- transpose.jar (https://faculty.washington.edu/browning/beagle_utilities/utilities.html - 'File manipulation utilities' section)

<br>


Optionally, If Python 'Anaconda(https://www.anaconda.com/)' is installed in your system, you can create a new independent Python virtual environment for CookHLA with the provided YML file.


By using 'CookHLA_LINUX.yml' or 'CookHLA_OSX.yml' file in the project folder **depending on your operating system**, Create a new Python virtual environment.

```
$ conda env create -f CookHLA_OSX.yml          ## OS X(Mac)
$ conda env create -f CookHLA_LINUX.yml        ## Linux
```

The above command will generate a new Python virtual environment named 'CookHLA', which contains dependent Python packages, R and R libraries, independent to your original Python system. For more detailed explanation about Anaconda's managing Python virtual environment, Please check this reference(https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#create-env-file-manually).

If the new virtual environment has been succuessfully created, then activate it.

```
$ conda activate CookHLA
```


CookHLA can be implemented in this virtual environment. 


<br>


> (Tip) Type 'conda acitvate base' on your command line if you want to go back to your original Python system setting. (https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#deactivating-an-environment)


> (Tip) Type 'conda env remove -n CookHLA' in your command line if you want to remove this newly created virtual environment for CookHLA forever in your Anaconda. (https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#removing-an-environment)

<br>
<br>



## () Implementational Heads-up

- Currently, CookHLA **ONLY SUPPORTS hg18** for genotype data. Target and Reference panel data must be in hg18 coordinate.
<!-- - **Family ID and Individual ID, 'FID' and 'IID' in PLINK data, have to be same**. -->
- Again, After downloading Beagle4('beagle.27Jan18.7e1.jar') and Beagle5('beagle.18May20.d20.jar') software in the 'dependency/' folder, Don't forget to rename them as 'beagle4.jar' and 'beagle5.jar'.


<br>
<br>


## () Usage example

```
$ python CookHLA.py \
    -i example/1958BC \
    -o MyHLAImputation/1958BC+HM_CEU_REF \
    -ref example/HM_CEU_REF \
    -gm example/AGM.1958BC+HM_CEU_REF.mach_step.avg.clpsB \
    -ae example/AGM.1958BC+HM_CEU_REF.aver.erate \
    -mem 2g \
    # -mp 2   # The number of available cores for Multiprocessing.
```

Input can be divided into 3 major components. 

1. Target genotype data to impute ('-i'), 
2. Reference panel data to be used in imputation ('-ref'),
3. Adaptive genetic map created based on the input 1 and 2 ('-gm' and '-ae').

In the output, '\*.alleles' file contains the imputed HLA types of the target genotype.



<br>
<br>


## () How to generate a Reference panel

CookHLA uses the exactly same reference panel of SNP2HLA. SNP2HLA is distributed with the additional software **MakeReference**, which can generate a reference panel for both SNP2HLA and CookHLA. The official webpage of SNP2HLA is "http://software.broadinstitute.org/mpg/snp2hla/" and the manual of **MakeReference** can be found here, "http://software.broadinstitute.org/mpg/snp2hla/makereference_manual.html".

<br>
<br>


## () How to generate an Adaptive Genetic Map

To generate an Adaptive genetic map, You have to use **'MakeGeneticMap' module** in the CookHLA project folder.

```
$ python -m MakeGeneticMap \
    -i example/1958BC \
    -ref example/HM_CEU_REF \
    -o MyAGM/1958BC+HM_CEU_REF
```

This implementation will generate (1) '\*.aver.erate' and (2) '\*.mach_step.avg.clpsB' files, each of which can be passed into '-ae' and '-gm' arguments of CookHLA respectively.

It is recomended that **the number of markers of target genotype data should be less than equal to that of reference panel data**. Otherwise, it can cause an error in the MakeGeneticMap process.


<br>
<br>



## () Available Reference panels

Some researchers might not have a reference panel or HLA type data to generate one. In this case, they should use public reference panels listed below.

1. T1DGC panel (Conditionally Public)

    - Webpage Link: http://software.broadinstitute.org/mpg/snp2hla/ - 'Links' section.
    - Citation: Jia X, Han B, Onengut-Gumuscu S, et al. Imputing amino acid polymorphisms in human leukocyte antigens. PLoS One. 2013;8(6):e64683. Published 2013 Jun 6. doi:10.1371/journal.pone.0064683

2.  Korean panel

    - Webpage Link: https://sites.google.com/site/scbaehanyang/hla_panel
    - Citation: 
        1. Kim K, Bang SY, Lee HS, Bae SC. Construction and application of a Korean reference panel for imputing classical alleles and amino acids of human leukocyte antigen genes. PLoS One. 2014;9(11):e112546. Published 2014 Nov 14. doi:10.1371/journal.pone.0112546.
        2. Kim K, Bang SY, Yoo DH, et al. Imputing Variants in HLA-DR Beta Genes Reveals That HLA-DRB1 Is Solely Associated with Rheumatoid Arthritis and Systemic Lupus Erythematosus. PLoS One. 2016;11(2):e0150283. Published 2016 Feb 26. doi:10.1371/journal.pone.0150283.

3. Pan-Asian panel

    - Webpage Link: http://software.broadinstitute.org/mpg/snp2hla/
    - Citation: 
        1. Pillai NE, Okada Y, Saw WY, et al. Predicting HLA alleles from high-resolution SNP data in three Southeast Asian populations. Hum Mol Genet. 2014;23(16):4443-4451. doi:10.1093/hmg/ddu149,
        2.  Okada Y, Kim K, Han B, et al. Risk for ACPA-positive rheumatoid arthritis is driven by shared HLA amino acid polymorphisms in Asian and European populations. Hum Mol Genet. 2014;23(25):6916-6926. doi:10.1093/hmg/ddu387.

4. HAN China panel

    - Webpage Link: http://gigadb.org/dataset/100156
    - Citation: Zhou F, Cao H, Zuo X, et al. Deep sequencing of the MHC region in the Chinese population contributes to studies of complex disease. Nat Genet. 2016;48(7):740-746. doi:10.1038/ng.3576



5. 1000 Genome SNP + HLA type data.

    1. SNP
        - Webpage Link: http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/
        - Citation: 1000 Genomes Project Consortium, Auton A, Brooks LD, et al. A global reference for human genetic variation. Nature. 2015;526(7571):68-74. doi:10.1038/nature15393

    2. HLA
        - Webpage Link: https://www.internationalgenome.org/category/hla/
        - Citation: Abi-Rached L, Gouret P, Yeh JH, et al. Immune diversity sheds light on missing variation in worldwide genetic diversity panels. PLoS One. 2018;13(10):e0206512. Published 2018 Oct 26. doi:10.1371/journal.pone.0206512



The 5th 1000G related data are not a reference panel but SNP and HLA type data that can be used to generate a reference panel. For researchers who can't afford to perform their own pre-processing and text processing, We provide the example 1000G reference panels that were used in the CookHLA paper. You can find them in 'xxxx' folder.

In that folder, there are 6 1000G reference panels where each of them represents a super population specified in 1000G project (https://www.internationalgenome.org/category/population/), i.e. AFR, AMR, EAS, EUR and SAS including ALL.

1. 
2. 
3. 
4. 
5. 
6. 

<!-- These reference panels went through next QCs and pre-processing.

- hg18:29-34mb
- subsetted to have only common markers that are shared in the T1DGC reference panel.
- excluded the markers that can't be -->

Researchers can use these reference panels to implement CookHLA.

<br>
<br>


## () Citation

S. Cook, W. Choi, H. Lim, Y. Luo, K. Kim, X. Jia, S. Raychaudhuri and B. Han, CookHLA: Accurate Imputation of Human Leukocyte Antigens. Under review.

<br>
<br>

## () License
The CookHLA Software Code is freely available for non-commercial academic research use. If you would like to obtain a license to the Code for commercial use, please contact Wanson Choi (WC) at wansonchoi@snu.ac.kr and Buhm Han (BH) at buhm.han@snu.ac.kr. WE (Seungho Cook, WC, Hyunjoon Lim and BH) MAKE NO REPRESENTATIONS OR WARRANTIES WHATSOEVER, EITHER EXPRESS OR IMPLIED, WITH RESPECT TO THE CODE PROVIDED HERE UNDER. IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE WITH RESPECT TO CODE ARE EXPRESSLY DISCLAIMED. THE CODE IS FURNISHED "AS IS" AND "WITH ALL FAULTS" AND DOWNLOADING OR USING THE CODE IS UNDERTAKEN AT YOUR OWN RISK. TO THE FULLEST EXTENT ALLOWED BY APPLICABLE LAW, IN NO EVENT SHALL WE BE LIABLE, WHETHER IN CONTRACT, TORT, WARRANTY, OR UNDER ANY STATUTE OR ON ANY OTHER BASIS FOR SPECIAL, INCIDENTAL, INDIRECT, PUNITIVE, MULTIPLE OR CONSEQUENTIAL DAMAGES SUSTAINED BY YOU OR ANY OTHER PERSON OR ENTITY ON ACCOUNT OF USE OR POSSESSION OF THE CODE, WHETHER OR NOT FORESEEABLE AND WHETHER OR NOT WE HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES, INCLUDING WITHOUT LIMITATION DAMAGES ARISING FROM OR RELATED TO LOSS OF USE, LOSS OF DATA, DOWNTIME, OR FOR LOSS OF REVENUE, PROFITS, GOODWILL, BUSINESS OR OTHER FINANCIAL LOSS.


<br>
<br>
