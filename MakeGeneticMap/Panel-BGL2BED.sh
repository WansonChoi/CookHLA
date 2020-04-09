#!/bin/bash

## BGL2PED
## Converts SNP2HLA beagle panel to plink format
## Launch: 11/28/14

## Arguments:
## $1: BGL prefix (X of X.bgl.phased and X.markers)
## $2: OUTPUT prefix
## $3: beagle2linkage.jar file
## $4: plink(v1.9) # added by Wanson Choi. (2020. 01. 14.)

BGL=$1
OUTPUT=$2
BEAGLE2LINKAGE=$3
PLINK=$4

cat $BGL.bgl.phased | java -jar $BEAGLE2LINKAGE $OUTPUT.tmp
cut -d ' ' -f-5 $OUTPUT.tmp.ped > $OUTPUT.ped.left
cut -d ' ' -f6- $OUTPUT.tmp.ped > $OUTPUT.ped.right
paste -d ' -9 ' $OUTPUT.ped.left /dev/null /dev/null /dev/null $OUTPUT.ped.right > $OUTPUT.ped

cut -d ' ' -f1 $BGL.markers > $OUTPUT.map.rsid
cut -d ' ' -f2 $BGL.markers > $OUTPUT.map.bp
cut -d ' ' -f3 $BGL.markers > $OUTPUT.map.allele1
paste -d '6  0 ' /dev/null /dev/null $OUTPUT.map.rsid /dev/null /dev/null $OUTPUT.map.bp > $OUTPUT.map
paste -d ' ' $OUTPUT.map.rsid $OUTPUT.map.allele1 > $OUTPUT.refallele

rm $OUTPUT.tmp.* $OUTPUT.map.* $OUTPUT.ped.*

${PLINK} --noweb --allow-no-sex --ped $OUTPUT.ped --map $OUTPUT.map --make-bed --reference-allele $OUTPUT.refallele --out $OUTPUT
${PLINK} --noweb --allow-no-sex --bfile $OUTPUT --keep-allele-order --freq --out $OUTPUT.FRQ

rm $OUTPUT.map $OUTPUT.ped $OUTPUT.log $OUTPUT.refallele $OUTPUT.FRQ.log



