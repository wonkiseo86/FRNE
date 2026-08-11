library(sde)
outiter<-1  ## Number of repetitions
initer<-100000 ## Number of calculated min-eigenvalues that used to calculate quantiles
numgrid<-1000  ## num.grid to evaluate BM
CV<-NULL
grid<-seq(0,1,0.001) ## quantiles

for (iiii in 1:5)
{

QTB<-NULL
for(jj in 1:outiter)
{

CR<-NULL
LEN<-numgrid
for(ii in 1:initer)
{
MM<-iiii
MMM<-NULL
IMMM=NULL
for (i in 1:MM)
{
aa=as.vector(BM(x=0, t0=0, T=1, N=LEN-1))
MMM=cbind(MMM,aa)
}
DMMM<-NULL
for (i in 1:ncol(MMM))
{
azz=MMM[,i]-mean(MMM[,i])
DMMM=cbind(DMMM, azz)
}

IMMM=NULL
for (i in 1:ncol(MMM))
{
IMMM=cbind(IMMM, cumsum(DMMM[,i])/LEN)
}

DMMMA=(t(DMMM)%*%(DMMM))/LEN
IMMMA=(t(IMMM)%*%(IMMM))/LEN  

AMMM=(DMMMA)%*%solve(IMMMA)
CR<-append(CR,sum(diag(AMMM)))
print(ii)
}
qqq<-quantile(CR,grid)


QTB=cbind(QTB,as.vector(qqq))

}

CV<-rbind(CV,rowSums(QTB)/ncol(QTB))

}
## CV contains critical vaues for quantiles in "grid" for R=1,...6





## Temperature anomalies
aa=c(11.73, 177.39, 1214.36, 3216.9, 7247.24)  #Test stats obtained from estimation
qqq[min(which(aa[5]<=CV[5,]))]
qqq[min(which(aa[4]<=CV[4,]))]
qqq[min(which(aa[3]<=CV[3,]))]
qqq[min(which(aa[2]<=CV[2,]))]
qqq[min(which(aa[1]<=CV[1,]))]




## GRP
aa=c(13.22,  102, 694.2, 1699.36, 3094.7)  #Test stats obtained from estimation
qqq[min(which(aa[5]<=CV[5,]))]
qqq[min(which(aa[4]<=CV[4,]))]
qqq[min(which(aa[3]<=CV[3,]))]
qqq[min(which(aa[2]<=CV[2,]))]
qqq[min(which(aa[1]<=CV[1,]))]