args<- commandArgs(trailingOnly = TRUE)

target_MARKER_DATA <- as.matrix(read.table(args[1],head=F))



reference_MARKER_DATA <- as.matrix(read.table(args[2],head=F))

#finding M

row_target_MARKER_DATA<-nrow(target_MARKER_DATA)
col_target_MARKER_DATA<-ncol(target_MARKER_DATA)






row_reference_MARKER_DATA<-nrow(reference_MARKER_DATA)
col_reference_MARKER_DATA<-ncol(reference_MARKER_DATA)








target_row_in_ref<-matrix(0,row_target_MARKER_DATA,1)
b<-0


for(a in 1:row_target_MARKER_DATA)
{
  i<-which(reference_MARKER_DATA[,1]==target_MARKER_DATA[a,1])
  if(length(i)>0){
    target_row_in_ref[b]<-i
    b=b+1
  }
  
}

error<-which(target_row_in_ref==0)

target_row_in_ref<- target_row_in_ref[-error]

new_target_nrow<-length(target_row_in_ref)

new_data_marker<-matrix(0,new_target_nrow,4)

d<-1

for(c in (target_row_in_ref))
{
  
  new_data_marker[d,]<-reference_MARKER_DATA[c,]
  
  
  d<-d+1
  
}



new_data_marker_sorted<-new_data_marker[order(new_data_marker[,2]),]



save_new_markers<-args[3]

write.table(new_data_marker_sorted,save_new_markers,sep = " " ,quote = FALSE,row.names = FALSE,col.names = FALSE)

