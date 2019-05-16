#!/usr/bin/env python

## This code temporarilly replaces Allele1 and 2 with G and C, because
## beagle does not understand allele names such as P(presence) and A(absence)
## Usage: python [this file] [beagle file] [marker file] [beagle output file] [marker output file]
## v1.2: 2017-07-05. python porting. Buhm Han.

import sys, os
#from itertools import izip

[bglfile, markerfile, bgloutfile, markeroutfile]=sys.argv[1:]

alleles={}
with open(markerfile) as mrk, open(markeroutfile, 'w') as mrkout:
    for l in mrk:
        [ID, pos, A1, A2]=l.split()
        alleles[ID]=(A1, A2)
        newl=' '.join([ID, pos, "G", "C"])+'\n'
        mrkout.write(newl)

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
            newdata=["G" if x==A1 else ("C" if x==A2 else "0") for x in data]
            newl=' '.join(c[:2]+newdata)+'\n'
            bglout.write(newl)
