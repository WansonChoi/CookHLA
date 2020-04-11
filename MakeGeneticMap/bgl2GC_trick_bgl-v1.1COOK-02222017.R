## TO CHANGE BEAGLE R PROJECT

args<- commandArgs(trailingOnly = TRUE)


beagle_DATA <- as.matrix(read.table(args[1]))
MARKER_DATA <- as.matrix(read.table(args[2]))


beagle_row <- (nrow(beagle_DATA))
beagle_col <- (ncol(beagle_DATA))
number_snp <- (nrow(MARKER_DATA))




first_order_chr <-c("G")
second_order_chr <-c("C")
thrid_order_chr <-c("G")





#finding M 

for(h in (1:beagle_row))
{
  if(beagle_DATA[h,1]=="M")
  {
    startSNP_row<-h
    break()
  }
  
  
}  

#checking not including snp

making_frist_line=beagle_row-number_snp

print("if this process is stop,please sort the marker and bgl by position(position is adjusted by redefineBPv1BH.py)")

for(t in (startSNP_row:beagle_row))
{
  stopifnot(beagle_DATA[t,2]==MARKER_DATA[t-making_frist_line,1])
  
  
}  

print("good snp order (marker and bgl_data) ")








for(g in (startSNP_row:beagle_row))
{
  snpline <- beagle_DATA[g,-(1:2)] 
  is.char1 = (snpline == MARKER_DATA[g-making_frist_line,3])
  is.char2 = (snpline == MARKER_DATA[g-making_frist_line,4])
  is.char3 = !is.char1 & !is.char2
  snpline[is.char1] <- first_order_chr
  snpline[is.char2] <- second_order_chr
  snpline[is.char3] <- thrid_order_chr
  beagle_DATA[g,-(1:2)] <- snpline 
  
  
  
}


#conversing GC overlap markers

for(k in (1:nrow(MARKER_DATA)))
{
  
  MARKER_DATA[,3]<-first_order_chr
  
  MARKER_DATA[,4]<-second_order_chr 
  
}

save_bgl<-args[3]
save_marker<-args[4]

write.table(beagle_DATA,save_bgl,sep = " " ,quote = FALSE,row.names = FALSE,col.names = FALSE)
write.table(MARKER_DATA,save_marker,sep = " " ,quote = FALSE,row.names = FALSE,col.names = FALSE)