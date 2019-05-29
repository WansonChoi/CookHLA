#-*- coding: utf-8 -*-
import sys, os

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]

def measureAccuracy(answerfile, predictfile, genes, outfile=None, __asSTDOUT = False):


    if isinstance(genes, list):
        if len(genes) == 0 or (len(genes) == 1 and genes[0] == "all"):
            genes = "A B C DRB1 DPA1 DPB1 DQA1 DQB1".split()
    elif isinstance(genes, str):
        if genes == 'all':
            genes = "A B C DPA1 DPB1 DQA1 DQB1 DRB1".split()


    # Accuracy dictionary
    __RETURN__ = {'2D' : {_hla: None for _hla in HLA_names},
                  '4D' : {_hla: None for _hla in HLA_names}}

    __asFileWrite = outfile and not __asSTDOUT

    if __asFileWrite:
        fo = open(outfile, 'w')


    for gene in genes:
        correct2d=0
        total2d=0
        correct4d=0
        total4d=0

        answers2d={}
        answers4d={}
        with open(answerfile) as fin:
            for l in fin:
                c=l.split()
                ID=c[0]+' '+c[1]
                if c[2] == gene and len(c) > 3:
                    answers2d[ID]=c[3].split(',')
                    answers4d[ID]=c[4].split(',')


        with open(predictfile) as fin:
            for l in fin:
                c=l.split()
                ID=c[0]+' '+c[1]
                if c[2] == gene:
                    predict2d=c[3].split(',')
                    predict4d=c[4].split(',')
                    if ID in answers2d:
                        (correct, total)=compare_and_score(predict2d, answers2d[ID])
                        correct2d+=correct
                        total2d+=total
                    if ID in answers4d:
                        (correct, total)=compare_and_score(predict4d, answers4d[ID])
                        correct4d+=correct
                        total4d+=total


        __RETURN__['2D'][gene] = float(correct2d)/total2d
        __RETURN__['4D'][gene] = float(correct4d)/total4d


        if __asSTDOUT:
            sys.stdout.write("%s\t2D\t%.5f\n"%(gene, float(correct2d)/total2d))
            sys.stdout.write("%s\t4D\t%.5f\n"%(gene, float(correct4d)/total4d))

        if __asFileWrite:
            fo.write("%s\t2D\t%.5f\n"%(gene, float(correct2d)/total2d))
            fo.write("%s\t4D\t%.5f\n"%(gene, float(correct4d)/total4d))

    if __asFileWrite:
        fo.close()


    return __RETURN__




def compare_and_score(predict, answer):
    ## return value is correct, total
    if answer[0] == '':
        answer[0]='answer1'
    if answer[1] == '':
        answer[1]='answer2'
    if predict[0] == '':
        predict[0]='predict1'
    if predict[1] == '':
        predict[1]='predict2'

    correct=max((answer[0]==predict[0])+(answer[1]==predict[1]), \
            (answer[0]==predict[1])+(answer[1]==predict[0]))
    return(correct,2)



if __name__ == "__main__":

    """
    < measureAccuracy.py >
    
    Module to get accuracy(%) of given '*.alleles' file.
    
    
    """

    answerfile = sys.argv[1]
    predictfile = sys.argv[2]
    genes = sys.argv[3:]

    measureAccuracy(answerfile, predictfile, genes, __asSTDOUT=True)

