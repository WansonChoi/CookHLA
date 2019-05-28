import sys, os, re

p_HLA_exon = re.compile(r'^HLA_\w+_(\d+)_\w+$')

def BGL2Alleles(bglfile, outfile, genes):

    if isinstance(genes, list):
        if len(genes) == 1 and genes[0] == "all":
            genes = "A B C DRB1 DPA1 DPB1 DQA1 DQB1".split()
    elif isinstance(genes, str):
        if genes == 'all':
            genes = "A B C DPA1 DPB1 DQA1 DQB1 DRB1".split()


    tmpfile=outfile+'.tmpfile'

    f = open(bglfile)

    # (Python 2.x.x)
    # FID = f.next().split()[2:]
    # IID = f.next().split()[2:]

    # (Python 3.x.x)
    FID = f.readline().split()[2:]
    IID = f.readline().split()[2:]
    f.close()

    # (Python 2.x.x)
    # N=len(IID)/2

    # (Python 3.x.x)
    N=int(len(IID)/2)

    alleles2d = {}
    alleles4d = {}

    for gene in genes:
      alleles2d[gene] = [[] for _ in range(N)] 
      reg = "HLA_%s_[0-9][0-9][0-9]?"%gene
#      print("Processing 2D of "+gene)
      os.system("egrep -w '%s' %s > %s"%(reg, bglfile, tmpfile))
      readAlleles(alleles2d[gene], tmpfile)
          
      alleles4d[gene] = [[] for _ in range(N)]
      reg = "HLA_%s_[0-9][0-9][0-9][0-9]"%gene
#      print("Processing 4D of "+gene)
      os.system("grep '%s' %s > %s"%(reg, bglfile, tmpfile))
      readAlleles(alleles4d[gene], tmpfile)

    ## Fill in empty string
    for i in range(N):
        for gene in genes:
            if len(alleles2d[gene][i]) == 1:
                alleles2d[gene][i].append('')
            elif len(alleles2d[gene][i]) == 0:
                alleles2d[gene][i].append('')
                alleles2d[gene][i].append('')
            if len(alleles4d[gene][i]) == 1:
                alleles4d[gene][i].append('')
            elif len(alleles4d[gene][i]) == 0:
                alleles4d[gene][i].append('')
                alleles4d[gene][i].append('')
      
    fo = open(outfile, "w")
    for i in range(N):
      for gene in genes:
        fo.write("%s\t%s\t%s\t"%(FID[2*i],IID[2*i],gene))
        fo.write(",".join(alleles2d[gene][i])+"\t")
        fo.write(",".join(alleles4d[gene][i])+"\n")
    fo.close()
                 
    os.system('rm %s'%tmpfile)


    return outfile


## Subroutin to read alleles
def readAlleles(alleles, tmpfile):
  for l in open(tmpfile):
    c = l.split()

    m = p_HLA_exon.match(string=c[1])

    if m:
        allele = c[1][c[1].rfind('_')+1:]
        presence = c[2:]
        for i in range(2*len(alleles)):
          if presence[i] == 'P':
              alleles[int(i/2)].append(allele)

    else:
        continue


          

if __name__ == "__main__":

    ## BGL file --> HLA allele converter.
    ## Buhm Han, 12/3/14
    ## Arguments:
    ## 1. Beagle File
    ## 2. Output File
    ## 3~ . Genes (DPB1, A, ...)

    [bglfile, outfile] = sys.argv[1:3]
    genes = sys.argv[3:]

    BGL2Alleles(bglfile, outfile, genes)
