args<- commandArgs(trailingOnly = TRUE)

input_phased_date<-as.matrix(read.table(args[1],colClasses = 'character'))


output_double_phased<-matrix(0,nrow(input_phased_date),9+2*(ncol(input_phased_date)-9))

#strsplit(input_phased_date[1,10],"[|]")[[1]][1]
#strsplit(input_phased_date[1,10],"[|]")[[1]][2]
#paste(strsplit(input_phased_date[1,10],"[|]")[[1]][1],"|",strsplit(input_phased_date[1,10],"[|]")[[1]][1],sep = "")
#
#


output_double_phased[,1:9]<-input_phased_date[,1:9]

for(k in (10:ncol(input_phased_date)))
{
  output_double_phased[1,2*k-10]<-input_phased_date[1,k]
  output_double_phased[1,2*k-9]<-paste(input_phased_date[1,k],"_1",sep = "")
  
}




for(j in (2:nrow(input_phased_date)))
{
  for(i in (10:ncol(input_phased_date)))
  {
    output_double_phased[j,2*i-10]<-paste(strsplit(input_phased_date[j,i],"[|]")[[1]][1],"|",strsplit(input_phased_date[j,i],"[|]")[[1]][1],sep = "")
    
    output_double_phased[j,2*i-9]<-paste(strsplit(input_phased_date[j,i],"[|]")[[1]][2],"|",strsplit(input_phased_date[j,i],"[|]")[[1]][2],sep = "")
    
    
    
  }
  
}



final_file<-args[2]

write.table(output_double_phased, final_file, sep = "\t",quote=F, col.names=F, row.names=F)
