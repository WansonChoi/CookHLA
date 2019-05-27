args<- commandArgs(trailingOnly = TRUE)


HLA_EXON<-as.matrix(read.table(args[1],colClasses = 'character'))


new_HLA_vcf<-matrix(0,nrow(HLA_EXON),ncol(HLA_EXON))
new_HLA_vcf[,1:9]<-HLA_EXON[1:(nrow(HLA_EXON)),1:9]
for(i in (1:(nrow(HLA_EXON))))
{
  for(j in (10:ncol(HLA_EXON)))
  {
    
    new_HLA_vcf[i,j]<-round(as.numeric(strsplit(HLA_EXON[i,j],":")[[1]][2]))   
    
  }
  
}



for(j in (10:ncol(new_HLA_vcf)))
{
  
  new_HLA_vcf[which(new_HLA_vcf[,j]==min(new_HLA_vcf[,j])),j]<-"0|0"
  
  
  
  
}


for(i in (1:(nrow(HLA_EXON))))
{
  for(j in (10:ncol(HLA_EXON)))
  {
    
    if(new_HLA_vcf[i,j]!="0|0")
    {
      new_HLA_vcf[i,j]<-"1|1"
      
    }
    
  }
  
}


new_HLA_vcf[,8]<-"."

new_HLA_vcf[,9]<-"GT"


final_file<-args[2]

write.table(new_HLA_vcf,final_file, sep = "\t",quote=F, col.names=F, row.names=F)

