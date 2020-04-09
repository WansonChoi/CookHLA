#!/usr/bin/env Rscript
args<- commandArgs(trailingOnly = TRUE)


mapdata=as.matrix(read.table(args[1]))
snps=mapdata[,2]
gmap=as.numeric(mapdata[,3])
n=length(snps)
gdist=gmap[2:n]-gmap[1:(n-1)]

## Special variables are these.
is.special=grepl("HLA_", snps) | grepl("AA_", snps) | grepl("SNP_", snps)

epsilon=1E-12

## Version A: Collapse distances between special variables.
gmap2=rep(0,n)
for (i in 2:n) {
    if (is.special[i-1] && is.special[i]) {
        gmap2[i]=gmap2[i-1]+epsilon
    } else {
        gmap2[i]=gmap2[i-1]+gdist[i-1]
    }
}

save_avg.clpsA<-args[2]
save_avg.clpsB<-args[3]

write.table(cbind(mapdata[,1:2], format(gmap2, digits=15), mapdata[,4]), save_avg.clpsA, quote=F, col.names=F, row.names=F)

## Version B: Collapse special variables into mid-point, while keeping other variables. 
gmap3=gmap
k=2
while(k<=n) {
    if (!is.special[k-1] && is.special[k]) { ## if I entered into the region of special variables,
        begin=k 
        j=k
        while(is.special[j]) { 
            j=j+1
        }
        end=j-1 ## special variables ended here. 
        midpoint=(gmap3[begin-1]+gmap3[end+1])/2
        gmap3[begin]=midpoint
        for (m in (begin+1):end) {
            gmap3[m]=gmap3[m-1]+epsilon ## collapse those into mid-point 
        }
        k=end 
    } 
    k=k+1
}

write.table(cbind(mapdata[,1:2], format(gmap3, digits=15), mapdata[,4]), save_avg.clpsB, quote=F, col.names=F, row.names=F)



