args<- commandArgs(trailingOnly = TRUE)


ori<-as.matrix(read.table(args[1]),head=F)

target<-as.matrix(read.table(args[2]),head=F)

row_ori<-nrow(ori)
col_ori<-ncol(ori)


row_target<-nrow(target)
col_target<-ncol(target)


for(j in (1:row_ori))
{
  
  if(ori[j,1]=="M")
  {
    head_line<-j-1
    break()
  }
}  
  





new_bgl<-matrix(0,row_target+head_line-1,col_target)


for(b in(1:j))
{new_bgl[b,1]<-ori[b,1]
new_bgl[b,2]<-ori[b,2]
}

for(i in (3:col_target))
{
  for(k in (3:col_ori))
  {
    if(target[1,i]==ori[2,k])
    {
      for(l in (1:head_line))
      {
      new_bgl[l,i]<-ori[l,k]
      }
      
    }
  }
  
  
  
}

new_bgl_line<-1

for(t in (2:row_target))
{
  
  new_bgl[new_bgl_line+head_line,]<-target[t,]
  new_bgl_line<-new_bgl_line+1
  
}
final_file<-args[3]

write.table(new_bgl,final_file,sep = " " ,quote = FALSE,row.names = FALSE,col.names = FALSE)



