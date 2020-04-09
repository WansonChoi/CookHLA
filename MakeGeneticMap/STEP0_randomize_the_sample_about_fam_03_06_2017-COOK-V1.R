args<- commandArgs(trailingOnly = TRUE)

FAM_DATA <- as.matrix(read.table(args[1],head=F))

now_FAM_DATA<-nrow(FAM_DATA)
col_FAM_DATA<-ncol(FAM_DATA)

randow_order<-sample(now_FAM_DATA)

line<-1

new_random_FAM_DATA<-matrix(0,now_FAM_DATA,col_FAM_DATA)

for(i in randow_order)
{
  
  new_random_FAM_DATA[line,]<-FAM_DATA[i,]
  
  line<-line+1
}  

result<-args[2]
write.table(new_random_FAM_DATA,result,sep = " " ,quote = FALSE,row.names = FALSE,col.names = FALSE)