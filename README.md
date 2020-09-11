# CookHLA

## (1) Introduction

CookHLA imputes HLA types of target patient genotype data based on a reference panel information.

[Improvement]

1. Multiple Marker Strategy
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


<br>


Optionally, If Python 'Anaconda(https://www.anaconda.com/)' is installed in your system, you can create a new independent Python virtual environment with the provided YML file.


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

- Currently, CookHLA **ONLY SUPPORTS hg18** for genotype data.
<!-- - **Family ID and Individual ID, 'FID' and 'IID' in PLINK data, have to be same**. -->


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

<br>
<br>


## () How to generate a Reference panel



<br>
<br>

## () How to generate a Adaptive Genetic Map



<br>
<br>



## () Available Reference panels