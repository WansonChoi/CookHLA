# CookHLA

## (1) Introduction

CookHLA imputes HLA types of target patient genotype data based on a reference panel information. Currently, CookHLA **ONLY SUPPORTS hg18** for genotype data. Also, **Family ID and Individual ID, 'FID' and 'IID' in PLINK data, have to be same**.

<br>
<br>

## (2) Installation

First, Prepare OS X(Mac) or Linux operating system. CookHLA currently doesn't support Windows. It was checked that CookHLA can work in the next specific operating systems. (Linux environment, especially CentOS, is most recommended.)

- Linux : 
    <!-- - Ubuntu 19.04(Disco Dingo)
    - Ubuntu 18.04.3 LTS(Bionic Beaver) -->
    - CentOS_7
    - Linux Mint 19.2 Cinnamon(Tina)
- OS X : 
    - Catalina(**with Bash NOT Zsh**)
    - Mojave

    In case of using Catalina OS X, **Make sure your default shell is 'Bash($)' not 'Zsh(%)'**. To change the default shell to Bash, Please reference this blog(https://www.howtogeek.com/444596/how-to-change-the-default-shell-to-bash-in-macos-catalina/).

<br>

Then, Download this project in somewhere directory of your OS X or Linux system. It will be assumed that 'git' command is already installed in your system.

```
$ git clone https://github.com/WansonChoi/CookHLA.git
$ cd CookHLA
```
<br>

We strongly recommend using the latest version of 'Anaconda(or Miniconda)' to set up CookHLA.


1. install Anaconda or Miniconda.

    - Anaconda : (https://www.anaconda.com/)
    - Miniconda : (https://docs.conda.io/en/latest/miniconda.html)

<br>

2. Create a new independent Python virtual environment with the given YML file.

	By using 'CookHLA_LINUX.yml' or 'CookHLA_OSX.yml' file in the project folder **depending on your operating system**, Create a new Python virtual environment.
    
	```
	$ conda env create -f CookHLA_OSX.yml          ## OS X(Mac)
	$ conda env create -f CookHLA_LINUX.yml        ## Linux
	```
	
	The above command will generate a new Python virtual environment named 'CookHLA', which contains dependent Python packages, R and R libraries, independent to your original Python system. For more detailed explanation about Anaconda's managing Python virtual environment, Please check this reference(https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#create-env-file-manually).

	If the new virtual environment has been succuessfully installed, then activate it.

	```
	$ conda activate CookHLA
	```


    CookHLA will be implemented in this virtual environment. 


<br>


> (Tip) Type 'conda acitvate base' on your command line if you want to go back to your original Python system setting. (https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#deactivating-an-environment)


> (Tip) Type 'conda env remove -n CookHLA' in your command line if you want to remove this newly created virtual environment for CookHLA forever in your Anaconda. (https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#removing-an-environment)

<br>
<br>


## (3) Usage example

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