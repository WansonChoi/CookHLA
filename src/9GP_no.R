arg <- commandArgs(trailingOnly = TRUE)
args <- arg[1]
args2 <- arg[2]
args3 <- arg[3]
args4 <- arg[4]
args5 <- arg[5]
args6 <- arg[6]
args7 <- arg[7]
args8 <- arg[8]
args9 <- arg[9]

gene <- arg[10]
HLA_EXON1 <- as.matrix(read.table(args,colClasses = 'character'))
HLA_EXON2 <- as.matrix(read.table(args2,colClasses = 'character'))
HLA_EXON3 <- as.matrix(read.table(args3,colClasses = 'character'))
HLA_EXON4 <- as.matrix(read.table(args4,colClasses = 'character'))
HLA_EXON5 <- as.matrix(read.table(args5,colClasses = 'character'))
HLA_EXON6 <- as.matrix(read.table(args6,colClasses = 'character'))
HLA_EXON7 <- as.matrix(read.table(args7,colClasses = 'character'))
HLA_EXON8 <- as.matrix(read.table(args8,colClasses = 'character'))
HLA_EXON9 <- as.matrix(read.table(args9,colClasses = 'character'))

a <- max(nrow(HLA_EXON1),nrow(HLA_EXON2),nrow(HLA_EXON3),nrow(HLA_EXON4),nrow(HLA_EXON5),nrow(HLA_EXON6),nrow(HLA_EXON7),nrow(HLA_EXON8),nrow(HLA_EXON9))
new_HLA_vcf<-matrix(0,a,ncol(HLA_EXON1))
new_HLA_vcf[,1:9]<-HLA_EXON1[1:a,1:9]

FID<-HLA_EXON1[1,10:ncol(HLA_EXON1)]
HLA<-matrix(0,ncol(HLA_EXON1)-9,7)
for(i in (2:a))
{
  for(j in (10:ncol(HLA_EXON1)))
  {
    EXON1 <- as.numeric(strsplit((strsplit(HLA_EXON1[i,j],":")[[1]][3]),",")[[1]][1])+as.numeric(strsplit((strsplit(HLA_EXON1[i,j],":")[[1]][3]),",")[[1]][2])/2
    EXON2 <- as.numeric(strsplit((strsplit(HLA_EXON2[i,j],":")[[1]][3]),",")[[1]][1])+as.numeric(strsplit((strsplit(HLA_EXON2[i,j],":")[[1]][3]),",")[[1]][2])/2
    EXON3 <- as.numeric(strsplit((strsplit(HLA_EXON3[i,j],":")[[1]][3]),",")[[1]][1])+as.numeric(strsplit((strsplit(HLA_EXON3[i,j],":")[[1]][3]),",")[[1]][2])/2
    EXON4 <- as.numeric(strsplit((strsplit(HLA_EXON4[i,j],":")[[1]][3]),",")[[1]][1])+as.numeric(strsplit((strsplit(HLA_EXON4[i,j],":")[[1]][3]),",")[[1]][2])/2
    EXON5 <- as.numeric(strsplit((strsplit(HLA_EXON5[i,j],":")[[1]][3]),",")[[1]][1])+as.numeric(strsplit((strsplit(HLA_EXON5[i,j],":")[[1]][3]),",")[[1]][2])/2
    EXON6 <- as.numeric(strsplit((strsplit(HLA_EXON6[i,j],":")[[1]][3]),",")[[1]][1])+as.numeric(strsplit((strsplit(HLA_EXON6[i,j],":")[[1]][3]),",")[[1]][2])/2
    if(i <=nrow(HLA_EXON7)){
      EXON7 <- as.numeric(strsplit((strsplit(HLA_EXON7[i,j],":")[[1]][3]),",")[[1]][1])+as.numeric(strsplit((strsplit(HLA_EXON7[i,j],":")[[1]][3]),",")[[1]][2])/2
      EXON8 <- as.numeric(strsplit((strsplit(HLA_EXON8[i,j],":")[[1]][3]),",")[[1]][1])+as.numeric(strsplit((strsplit(HLA_EXON8[i,j],":")[[1]][3]),",")[[1]][2])/2
      EXON9 <- as.numeric(strsplit((strsplit(HLA_EXON9[i,j],":")[[1]][3]),",")[[1]][1])+as.numeric(strsplit((strsplit(HLA_EXON9[i,j],":")[[1]][3]),",")[[1]][2])/2
      new_HLA_vcf[i,j] <- (max(EXON1,EXON2,EXON3)+max(EXON4,EXON5,EXON6)+max(EXON7,EXON8,EXON9))/3
    }
    else{
      new_HLA_vcf[i,j] <- (max(EXON1,EXON2,EXON3)+max(EXON4,EXON5,EXON6))/2
    }
  }
}

for(i in (10:ncol(new_HLA_vcf)))
{
  HLA[i-9,1] <- FID[i-9]
  HLA[i-9,2] <- FID[i-9]
  HLA[i-9,3] <- gene
  HLA[i-9,4] <- new_HLA_vcf[which(new_HLA_vcf[,i]==max(new_HLA_vcf[,i])),3][1]
  HLA[i-9,6] <- new_HLA_vcf[which(new_HLA_vcf[,i]==max(new_HLA_vcf[,i])),i][1]
  t <- 0
  for(j in (2:(nrow(new_HLA_vcf))))
  {
    if (new_HLA_vcf[j,3] != HLA[i-9,4] && new_HLA_vcf[j,i] > as.numeric(HLA[i-9,6])/2){
      HLA[i-9,5] <- new_HLA_vcf[j,3]
      HLA[i-9,7] <- new_HLA_vcf[j,i]
      t <- 1
      break
    }
  }
  if(t==0){
    HLA[i-9,5] <- HLA[i-9,4]
    HLA[i-9,7] <- HLA[i-9,6]
  }
}
for(i in (1:nrow(HLA)))
{
  i1 <- strsplit(HLA[i,4],"_")[[1]][3]
  i2 <- strsplit(HLA[i,5],"_")[[1]][3]
  i3 <- paste0(strsplit(i1,"")[[1]][1],strsplit(i1,"")[[1]][2])
  i4 <- paste0(strsplit(i2,"")[[1]][1],strsplit(i2,"")[[1]][2])
  HLA[i,4] <- paste0(i3,",",i4)
  HLA[i,5] <- paste0(i1,",",i2)
}

write.table(HLA,paste0(args,".alleles"),quote=F, col.names=F, row.names=F)
