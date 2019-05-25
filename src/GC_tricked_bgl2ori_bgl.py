import sys, os


def GCtricedBGL2OriginalBGL(bglfile, markerfile, bgloutfile):

    alleles={}
    with open(markerfile) as mrk:
        for l in mrk:
            [ID, pos, A1, A2]=l.split()
            alleles[ID]=(A1, A2)


    with open(bglfile) as bgl, open(bgloutfile, 'w') as bglout:
        for l in bgl:
            c=l.split()
            if c[0] != 'M':
                bglout.write(l)
            else: ## If 'M', we replace alleles.
                ID=c[1]
                assert (ID in alleles), "Marker ID %s is not in the marker file"%ID
                (A1, A2)=alleles[ID]
                data=c[2:]
                newdata=[A1 if x=="G" else (A2 if x=="C" else "0") for x in data]
                newl=' '.join(c[:2]+newdata)+'\n'
                bglout.write(newl)


    return bgloutfile



if __name__ == '__main__':

    """
    
    ## This code temporarilly replaces Allele1 and 2 with G and C, because
    ## beagle does not understand allele names such as P(presence) and A(absence)
    ## Usage: python [this file] [beagle file] [marker file] [beagle output file] [marker output file]
    ## v1.2: 2017-07-05. python porting. Buhm Han.

    """

    [bglfile, markerfile, bgloutfile] = sys.argv[1:]

    GCtricedBGL2OriginalBGL(bglfile, markerfile, bgloutfile)