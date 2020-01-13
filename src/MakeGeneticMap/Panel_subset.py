#-*- coding: utf-8 -*-

import sys, os


def Panel_Subset(_panelprefix, _indvfile, _markerchoicefile, _outprefix):

    ## Read individuals and markers to select
    allIndividual=False
    selectIDs={}
    if _indvfile == "all":
        allIndividual=True
    else:
        with open(_indvfile) as fin:
            for l in fin:
                c=l.split()
                ID=c[0]+' '+c[1]
                selectIDs[ID]=True

    allMarker=False
    selectMarkers={}
    if _markerchoicefile == "all":
        allMarker=True
    else:
        with open(_markerchoicefile) as fin:
            for l in fin:
                c=l.split()
                selectMarkers[c[0]]=True

    ##-----------------------------
    ## Now perform filtering
    ##-----------------------------
    ## 1. "Tag" which columns to choose, first.
    with open(_panelprefix + '.bgl.phased') as pf:
        FID=pf.readline().split()
        IID=pf.readline().split()
        assert FID[0] == 'P' and IID[0] == 'I'
        IDs=[x+' '+y for (x,y) in zip(FID[2:], IID[2:])]
        tag=[False]*len(IDs)
        for i in range(len(IDs)):
            if allIndividual or IDs[i] in selectIDs:
                tag[i]=True

    ## 2. Filter phased file based on tag.
    with open(_panelprefix + '.bgl.phased') as pf, open(_outprefix + '.bgl.phased', 'w') as of:
        for l in pf:
            c=l.split()
            header=c[:2]
            data=c[2:]
            if header[0]=="M" and allMarker==False and header[1] not in selectMarkers:
                continue
            newdata=[]
            for j in range(len(data)):
                if tag[j]:
                    newdata.append(data[j])
            of.write(' '.join(header+newdata)+'\n')

    ## 3. Finally, filter markers.
    with open(_panelprefix + '.markers') as mf, open(_outprefix + '.markers', 'w') as of:
        for l in mf:
            [rsid, bp, A1, A2]=l.split()
            if allMarker==False and rsid not in selectMarkers:
                continue
            of.write(l)

    return _outprefix


if __name__ == '__main__':

    ####################################################################
    ## Select subset of individuals (and/or markers) from BGL format
    ## Launch: 12-2-2014
    ## Author: Buhm Han
    ## Arguments:
    ##   1: BGL prefix (X of X.bgl.phased and X.markers)
    ##   2: individual file (FID / IID pairs to select) <-- if "all", select all
    ##   3: marker choice file (markers to select) <-- if "all", select all
    ##   4: Output prefix
    ####################################################################

    [panelprefix, indvfile, markerchoicefile, outprefix] = sys.argv[1:]

    Panel_Subset(panelprefix, indvfile, markerchoicefile, outprefix)


