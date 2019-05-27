##decode_double_vcf_alleles

args<- commandArgs(trailingOnly = TRUE)


input_double_alleles<-as.matrix(read.table(args[1],colClasses = 'character'))


alleles<-matrix(0,nrow(input_double_alleles)/2,ncol(input_double_alleles))

for(j in (1:(nrow(input_double_alleles)/16)-1))
{for(i in (1:8))
{ 
  
  alleles[8*j+i,1]<-input_double_alleles[8*2*j+i,1]
  alleles[8*j+i,2]<-input_double_alleles[8*2*j+i,1]
  alleles[8*j+i,3]<-input_double_alleles[8*2*j+i,3]
  alleles[8*j+i,4]<-paste(strsplit(input_double_alleles[8*2*j+i,4],",")[[1]][1],strsplit(input_double_alleles[8*2*j+i+8,4],",")[[1]][1],sep = ",")
  
  alleles[8*j+i,5]<-paste(strsplit(input_double_alleles[8*2*j+i,5],",")[[1]][1],strsplit(input_double_alleles[8*2*j+i+8,5],",")[[1]][1],sep = ",")
  
  
  
}
  
  
} 

final_file<-args[2]

write.table(alleles, final_file, sep = "\t",quote=F, col.names=F, row.names=F)

