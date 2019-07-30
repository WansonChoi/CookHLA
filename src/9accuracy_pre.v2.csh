set INPUT=$1
set INPUT2=$2
set INPUT3=$3
set INPUT4=$4
set INPUT5=$5
set INPUT6=$6
set INPUT7=$7
set INPUT8=$8
set INPUT9=$9

set out=$10

grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT.HLA_A
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT.HLA_B
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT.HLA_C
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT.HLA_DRB1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT.HLA_DQA1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT.HLA_DQB1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT.HLA_DPA1
grep "#CHROM" $INPUT | sed -e 's/^#//' > $INPUT.HLA_DPB1

grep "#CHROM" $INPUT2 | sed -e 's/^#//' > $INPUT2.HLA_A
grep "#CHROM" $INPUT2 | sed -e 's/^#//' > $INPUT2.HLA_B
grep "#CHROM" $INPUT2 | sed -e 's/^#//' > $INPUT2.HLA_C
grep "#CHROM" $INPUT2 | sed -e 's/^#//' > $INPUT2.HLA_DRB1
grep "#CHROM" $INPUT2 | sed -e 's/^#//' > $INPUT2.HLA_DQA1
grep "#CHROM" $INPUT2 | sed -e 's/^#//' > $INPUT2.HLA_DQB1
grep "#CHROM" $INPUT2 | sed -e 's/^#//' > $INPUT2.HLA_DPA1
grep "#CHROM" $INPUT2 | sed -e 's/^#//' > $INPUT2.HLA_DPB1

grep "#CHROM" $INPUT3 | sed -e 's/^#//' > $INPUT3.HLA_A
grep "#CHROM" $INPUT3 | sed -e 's/^#//' > $INPUT3.HLA_B
grep "#CHROM" $INPUT3 | sed -e 's/^#//' > $INPUT3.HLA_C
grep "#CHROM" $INPUT3 | sed -e 's/^#//' > $INPUT3.HLA_DRB1
grep "#CHROM" $INPUT3 | sed -e 's/^#//' > $INPUT3.HLA_DQA1
grep "#CHROM" $INPUT3 | sed -e 's/^#//' > $INPUT3.HLA_DQB1
grep "#CHROM" $INPUT3 | sed -e 's/^#//' > $INPUT3.HLA_DPA1
grep "#CHROM" $INPUT3 | sed -e 's/^#//' > $INPUT3.HLA_DPB1

grep "#CHROM" $INPUT4 | sed -e 's/^#//' > $INPUT4.HLA_A
grep "#CHROM" $INPUT4 | sed -e 's/^#//' > $INPUT4.HLA_B
grep "#CHROM" $INPUT4 | sed -e 's/^#//' > $INPUT4.HLA_C
grep "#CHROM" $INPUT4 | sed -e 's/^#//' > $INPUT4.HLA_DRB1
grep "#CHROM" $INPUT4 | sed -e 's/^#//' > $INPUT4.HLA_DQA1
grep "#CHROM" $INPUT4 | sed -e 's/^#//' > $INPUT4.HLA_DQB1
grep "#CHROM" $INPUT4 | sed -e 's/^#//' > $INPUT4.HLA_DPA1
grep "#CHROM" $INPUT4 | sed -e 's/^#//' > $INPUT4.HLA_DPB1

grep "#CHROM" $INPUT5 | sed -e 's/^#//' > $INPUT5.HLA_A
grep "#CHROM" $INPUT5 | sed -e 's/^#//' > $INPUT5.HLA_B
grep "#CHROM" $INPUT5 | sed -e 's/^#//' > $INPUT5.HLA_C
grep "#CHROM" $INPUT5 | sed -e 's/^#//' > $INPUT5.HLA_DRB1
grep "#CHROM" $INPUT5 | sed -e 's/^#//' > $INPUT5.HLA_DQA1
grep "#CHROM" $INPUT5 | sed -e 's/^#//' > $INPUT5.HLA_DQB1
grep "#CHROM" $INPUT5 | sed -e 's/^#//' > $INPUT5.HLA_DPA1
grep "#CHROM" $INPUT5 | sed -e 's/^#//' > $INPUT5.HLA_DPB1

grep "#CHROM" $INPUT6 | sed -e 's/^#//' > $INPUT6.HLA_A
grep "#CHROM" $INPUT6 | sed -e 's/^#//' > $INPUT6.HLA_B
grep "#CHROM" $INPUT6 | sed -e 's/^#//' > $INPUT6.HLA_C
grep "#CHROM" $INPUT6 | sed -e 's/^#//' > $INPUT6.HLA_DRB1
grep "#CHROM" $INPUT6 | sed -e 's/^#//' > $INPUT6.HLA_DQA1
grep "#CHROM" $INPUT6 | sed -e 's/^#//' > $INPUT6.HLA_DQB1
grep "#CHROM" $INPUT6 | sed -e 's/^#//' > $INPUT6.HLA_DPA1
grep "#CHROM" $INPUT6 | sed -e 's/^#//' > $INPUT6.HLA_DPB1

grep "#CHROM" $INPUT7 | sed -e 's/^#//' > $INPUT7.HLA_A
grep "#CHROM" $INPUT7 | sed -e 's/^#//' > $INPUT7.HLA_B
grep "#CHROM" $INPUT7 | sed -e 's/^#//' > $INPUT7.HLA_C
grep "#CHROM" $INPUT7 | sed -e 's/^#//' > $INPUT7.HLA_DRB1
grep "#CHROM" $INPUT7 | sed -e 's/^#//' > $INPUT7.HLA_DQA1
grep "#CHROM" $INPUT7 | sed -e 's/^#//' > $INPUT7.HLA_DQB1
grep "#CHROM" $INPUT7 | sed -e 's/^#//' > $INPUT7.HLA_DPA1
grep "#CHROM" $INPUT7 | sed -e 's/^#//' > $INPUT7.HLA_DPB1

grep "#CHROM" $INPUT8 | sed -e 's/^#//' > $INPUT8.HLA_A
grep "#CHROM" $INPUT8 | sed -e 's/^#//' > $INPUT8.HLA_B
grep "#CHROM" $INPUT8 | sed -e 's/^#//' > $INPUT8.HLA_C
grep "#CHROM" $INPUT8 | sed -e 's/^#//' > $INPUT8.HLA_DRB1
grep "#CHROM" $INPUT8 | sed -e 's/^#//' > $INPUT8.HLA_DQA1
grep "#CHROM" $INPUT8 | sed -e 's/^#//' > $INPUT8.HLA_DQB1
grep "#CHROM" $INPUT8 | sed -e 's/^#//' > $INPUT8.HLA_DPA1
grep "#CHROM" $INPUT8 | sed -e 's/^#//' > $INPUT8.HLA_DPB1

grep "#CHROM" $INPUT9 | sed -e 's/^#//' > $INPUT9.HLA_A
grep "#CHROM" $INPUT9 | sed -e 's/^#//' > $INPUT9.HLA_B
grep "#CHROM" $INPUT9 | sed -e 's/^#//' > $INPUT9.HLA_C
grep "#CHROM" $INPUT9 | sed -e 's/^#//' > $INPUT9.HLA_DRB1
grep "#CHROM" $INPUT9 | sed -e 's/^#//' > $INPUT9.HLA_DQA1
grep "#CHROM" $INPUT9 | sed -e 's/^#//' > $INPUT9.HLA_DQB1
grep "#CHROM" $INPUT9 | sed -e 's/^#//' > $INPUT9.HLA_DPA1
grep "#CHROM" $INPUT9 | sed -e 's/^#//' > $INPUT9.HLA_DPB1



grep HLA_A $INPUT >> $INPUT.HLA_A
grep HLA_B $INPUT >> $INPUT.HLA_B
grep HLA_C $INPUT >> $INPUT.HLA_C
grep HLA_DRB1 $INPUT >> $INPUT.HLA_DRB1
grep HLA_DQA1 $INPUT >> $INPUT.HLA_DQA1
grep HLA_DQB1 $INPUT >> $INPUT.HLA_DQB1
grep HLA_DPA1 $INPUT >> $INPUT.HLA_DPA1
grep HLA_DPB1 $INPUT >> $INPUT.HLA_DPB1

grep HLA_A $INPUT2 >> $INPUT2.HLA_A
grep HLA_B $INPUT2 >> $INPUT2.HLA_B
grep HLA_C $INPUT2 >> $INPUT2.HLA_C
grep HLA_DRB1 $INPUT2 >> $INPUT2.HLA_DRB1
grep HLA_DQA1 $INPUT2 >> $INPUT2.HLA_DQA1
grep HLA_DQB1 $INPUT2 >> $INPUT2.HLA_DQB1
grep HLA_DPA1 $INPUT2 >> $INPUT2.HLA_DPA1
grep HLA_DPB1 $INPUT2 >> $INPUT2.HLA_DPB1

grep HLA_A $INPUT3 >> $INPUT3.HLA_A
grep HLA_B $INPUT3 >> $INPUT3.HLA_B
grep HLA_C $INPUT3 >> $INPUT3.HLA_C
grep HLA_DRB1 $INPUT3 >> $INPUT3.HLA_DRB1
grep HLA_DQA1 $INPUT3 >> $INPUT3.HLA_DQA1
grep HLA_DQB1 $INPUT3 >> $INPUT3.HLA_DQB1
grep HLA_DPA1 $INPUT3 >> $INPUT3.HLA_DPA1
grep HLA_DPB1 $INPUT3 >> $INPUT3.HLA_DPB1


grep HLA_A $INPUT4 >> $INPUT4.HLA_A
grep HLA_B $INPUT4 >> $INPUT4.HLA_B
grep HLA_C $INPUT4 >> $INPUT4.HLA_C
grep HLA_DRB1 $INPUT4 >> $INPUT4.HLA_DRB1
grep HLA_DQA1 $INPUT4 >> $INPUT4.HLA_DQA1
grep HLA_DQB1 $INPUT4 >> $INPUT4.HLA_DQB1
grep HLA_DPA1 $INPUT4 >> $INPUT4.HLA_DPA1
grep HLA_DPB1 $INPUT4 >> $INPUT4.HLA_DPB1

grep HLA_A $INPUT5 >> $INPUT5.HLA_A
grep HLA_B $INPUT5 >> $INPUT5.HLA_B
grep HLA_C $INPUT5 >> $INPUT5.HLA_C
grep HLA_DRB1 $INPUT5 >> $INPUT5.HLA_DRB1
grep HLA_DQA1 $INPUT5 >> $INPUT5.HLA_DQA1
grep HLA_DQB1 $INPUT5 >> $INPUT5.HLA_DQB1
grep HLA_DPA1 $INPUT5 >> $INPUT5.HLA_DPA1
grep HLA_DPB1 $INPUT5 >> $INPUT5.HLA_DPB1

grep HLA_A $INPUT6 >> $INPUT6.HLA_A
grep HLA_B $INPUT6 >> $INPUT6.HLA_B
grep HLA_C $INPUT6 >> $INPUT6.HLA_C
grep HLA_DRB1 $INPUT6 >> $INPUT6.HLA_DRB1
grep HLA_DQA1 $INPUT6 >> $INPUT6.HLA_DQA1
grep HLA_DQB1 $INPUT6 >> $INPUT6.HLA_DQB1
grep HLA_DPA1 $INPUT6 >> $INPUT6.HLA_DPA1
grep HLA_DPB1 $INPUT6 >> $INPUT6.HLA_DPB1

grep HLA_A $INPUT7 >> $INPUT7.HLA_A
grep HLA_B $INPUT7 >> $INPUT7.HLA_B
grep HLA_C $INPUT7 >> $INPUT7.HLA_C
grep HLA_DRB1 $INPUT7 >> $INPUT7.HLA_DRB1
grep HLA_DQA1 $INPUT7 >> $INPUT7.HLA_DQA1
grep HLA_DQB1 $INPUT7 >> $INPUT7.HLA_DQB1
grep HLA_DPA1 $INPUT7 >> $INPUT7.HLA_DPA1
grep HLA_DPB1 $INPUT7 >> $INPUT7.HLA_DPB1

grep HLA_A $INPUT8 >> $INPUT8.HLA_A
grep HLA_B $INPUT8 >> $INPUT8.HLA_B
grep HLA_C $INPUT8 >> $INPUT8.HLA_C
grep HLA_DRB1 $INPUT8 >> $INPUT8.HLA_DRB1
grep HLA_DQA1 $INPUT8 >> $INPUT8.HLA_DQA1
grep HLA_DQB1 $INPUT8 >> $INPUT8.HLA_DQB1
grep HLA_DPA1 $INPUT8 >> $INPUT8.HLA_DPA1
grep HLA_DPB1 $INPUT8 >> $INPUT8.HLA_DPB1

grep HLA_A $INPUT9 >> $INPUT9.HLA_A
grep HLA_B $INPUT9 >> $INPUT9.HLA_B
grep HLA_C $INPUT9 >> $INPUT9.HLA_C
grep HLA_DRB1 $INPUT9 >> $INPUT9.HLA_DRB1
grep HLA_DQA1 $INPUT9 >> $INPUT9.HLA_DQA1
grep HLA_DQB1 $INPUT9 >> $INPUT9.HLA_DQB1
grep HLA_DPA1 $INPUT9 >> $INPUT9.HLA_DPA1
grep HLA_DPB1 $INPUT9 >> $INPUT9.HLA_DPB1



Rscript 9GP_merge.R $INPUT.HLA_A $INPUT2.HLA_A $INPUT3.HLA_A $INPUT4.HLA_A $INPUT5.HLA_A $INPUT6.HLA_A $INPUT7.HLA_A $INPUT8.HLA_A $INPUT9.HLA_A A
Rscript 9GP_merge.R $INPUT.HLA_B $INPUT2.HLA_B $INPUT3.HLA_B $INPUT4.HLA_B $INPUT5.HLA_B $INPUT6.HLA_B $INPUT7.HLA_B $INPUT8.HLA_B $INPUT9.HLA_B B
Rscript 9GP_merge.R $INPUT.HLA_C $INPUT2.HLA_C $INPUT3.HLA_C $INPUT4.HLA_C $INPUT5.HLA_C $INPUT6.HLA_C $INPUT7.HLA_C $INPUT8.HLA_C $INPUT9.HLA_C C
Rscript 9GP_merge.R $INPUT.HLA_DRB1 $INPUT2.HLA_DRB1 $INPUT3.HLA_DRB1 $INPUT4.HLA_DRB1 $INPUT5.HLA_DRB1 $INPUT6.HLA_DRB1 $INPUT7.HLA_DRB1 $INPUT8.HLA_DRB1 $INPUT9.HLA_DRB1 DRB1
Rscript 9GP_merge.R $INPUT.HLA_DQA1 $INPUT2.HLA_DQA1 $INPUT3.HLA_DQA1 $INPUT4.HLA_DQA1 $INPUT5.HLA_DQA1 $INPUT6.HLA_DQA1 $INPUT7.HLA_DQA1 $INPUT8.HLA_DQA1 $INPUT9.HLA_DQA1 DQA1
Rscript 9GP_merge.R $INPUT.HLA_DQB1 $INPUT2.HLA_DQB1 $INPUT3.HLA_DQB1 $INPUT4.HLA_DQB1 $INPUT5.HLA_DQB1 $INPUT6.HLA_DQB1 $INPUT7.HLA_DQB1 $INPUT8.HLA_DQB1 $INPUT9.HLA_DQB1 DQB1
Rscript 9GP_merge.R $INPUT.HLA_DPA1 $INPUT2.HLA_DPA1 $INPUT3.HLA_DPA1 $INPUT4.HLA_DPA1 $INPUT5.HLA_DPA1 $INPUT6.HLA_DPA1 $INPUT7.HLA_DPA1 $INPUT8.HLA_DPA1 $INPUT9.HLA_DPA1 DPA1
Rscript 9GP_merge.R $INPUT.HLA_DPB1 $INPUT2.HLA_DPB1 $INPUT3.HLA_DPB1 $INPUT4.HLA_DPB1 $INPUT5.HLA_DPB1 $INPUT6.HLA_DPB1 $INPUT7.HLA_DPB1 $INPUT8.HLA_DPB1 $INPUT9.HLA_DPB1 DPB1


### Merge the outputs from '9GP_merge.R'
# cat *.alleles > $INPUT.9.merged
cat ${INPUT}.HLA_A.alleles \
    ${INPUT}.HLA_B.alleles \
    ${INPUT}.HLA_C.alleles \
    ${INPUT}.HLA_DRB1.alleles \
    ${INPUT}.HLA_DQA1.alleles \
    ${INPUT}.HLA_DQB1.alleles \
    ${INPUT}.HLA_DPA1.alleles \
    ${INPUT}.HLA_DPB1.alleles > ${out}.alleles

# python measureAccuracy4D.py $ANSWER $INPUT.9.merged > $INPUT.9.accuracy
