#!/usr/bin/env Rscript
args<- commandArgs(trailingOnly = TRUE)



erate=as.matrix(read.table(args[1],header=T))
snps=erate[,1]
n=length(snps)
avgrate=as.numeric(erate[,2])
lastrate=as.numeric(erate[,3])

write.table(as.matrix(c(mean(avgrate))), args[6], quote=F, row.names=F, col.names=F)
#cat(mean(avgrate), "\n")  # (2020. 01. 09.) Redirection has been deprecated.

rec=as.matrix(read.table(args[2],header=T))
avgrec=as.numeric(rec[,2])
lastrec=as.numeric(rec[,3])

bpos=as.matrix(read.table(args[3]))[,2]

myrec=avgrec
gmap=rep(0, n)
for (i in 2:n) {
    ## Assuming that MACH output is theta/|H|
    gmap[i]=gmap[i-1]+(-1/2)*log(1-myrec[i-1]) 
}

save_avg<-args[4]
save_last<-args[5]



write.table(cbind("6", snps, gmap, bpos), save_avg, quote=F, col.names=F, row.names=F)

myrec=lastrec
gmap=rep(0, n)
for (i in 2:n) {
    ## Assuming that MACH output is theta/|H|
    gmap[i]=gmap[i-1]+(-1/2)*log(1-myrec[i-1]) 
}
write.table(cbind("6", snps, gmap, bpos), save_last, quote=F, col.names=F, row.names=F)






