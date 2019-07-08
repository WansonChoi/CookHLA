#-*- coding: utf-8 -*-
import subprocess

def RUN_Bash(_command, __print=False, __save_intermediates=False):


    if __print:
        print(_command)

    sb = subprocess.call(_command, shell=True)

    if not sb:
        return 0    # Success
    else:
        return -1   # Fail



if __name__ == '__main__':

    """
    Suite of modules related to bash execution.
    Primarily designed for the modules to be imported.
    """

    pass