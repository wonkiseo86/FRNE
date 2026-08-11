
#####################################################
######## PREPARATION ################################
#####################################################

# Set working directory (change as needed)
setwd("path") 

# Load required library  
library(R.matlab)
library(sde)
library(robustbase)

# Define functions 
# numerical integration over grid
inner = function(f,g,grid){
  h = f*g
  return(sum((0.5*h[1:(length(grid)-1)] + 0.5*h[2:(length(grid))])*(grid[2] - grid[1])))
}

#####################################################
#####################################################
######## Some base parameters setup  ################
#####################################################
#####################################################	
LBF=LBFX

Scut=30
Non_dimX=3
Non_dimY=Non_dimX
decfac=0.8
decfac1=0.1
decfac2=0.8

bdw = 1/4
cuteigen=0.4

result_total=NULL
result_total_K=NULL
sparseset=c(1,0)
esignal=100       ## 10, 50, 100
ea=1; eb=Non_dimX+1; eside=Non_dimX+1

for (sparse in sparseset)
{
set.seed(12345)
	

TSET=c(100,200,400,800)
Kscheme=1
K_fix=4	# Not used.
addfac=0; ranst=1;	# For more general setups but not used in the paper. 

perturb=1;geotail=18;sparsecut=7

#####################################################
#####################################################
######## Simulation DGP and Monte Carlo repetition ##
#####################################################
#####################################################	
## Basic parameters setup
T_sim = max(TSET) ; burnin=150 ;  
KAPP = c(0,1)  # values of kappa 
maxiter=3000

RESULT1=matrix(0,nrow=maxiter,ncol=length(KAPP));RESULT1N=matrix(0,nrow=maxiter,ncol=length(KAPP));RESULT1S=matrix(0,nrow=maxiter,ncol=length(KAPP))
RESULT2=matrix(0,nrow=maxiter,ncol=length(KAPP));RESULT2N=matrix(0,nrow=maxiter,ncol=length(KAPP));RESULT2S=matrix(0,nrow=maxiter,ncol=length(KAPP))
RESULT3=matrix(0,nrow=maxiter,ncol=length(KAPP));RESULT3N=matrix(0,nrow=maxiter,ncol=length(KAPP));RESULT3S=matrix(0,nrow=maxiter,ncol=length(KAPP))
RESULT4=matrix(0,nrow=maxiter,ncol=length(KAPP));RESULT4N=matrix(0,nrow=maxiter,ncol=length(KAPP));RESULT4S=matrix(0,nrow=maxiter,ncol=length(KAPP))

RESULTK1=matrix(0,nrow=maxiter,ncol=length(KAPP))
RESULTK2=matrix(0,nrow=maxiter,ncol=length(KAPP))
RESULTK3=matrix(0,nrow=maxiter,ncol=length(KAPP))
RESULTK4=matrix(0,nrow=maxiter,ncol=length(KAPP))

lbnumber2=50
nt=200;  t = (0:(nt-1))/(nt-1)
LBF = matrix(NA,nrow = nt , ncol = lbnumber2)
for (i in 1:(lbnumber2/2)){
  LBF[,2*i-1] = sqrt(2)*sin(2*pi*i*t) /sqrt(inner(sqrt(2)*sin(2*pi*i*t),sqrt(2)*sin(2*pi*i*t),t))
  LBF[,2*i] = sqrt(2)*cos(2*pi*i*t)/sqrt(inner(sqrt(2)*cos(2*pi*i*t),sqrt(2)*cos(2*pi*i*t),t))
}

for(i in 2:lbnumber2){  
  for(j in 1:i)  { 
    if (j != i) {LBF[,i] = LBF[,i]-(inner(LBF[,i],LBF[,j],t)/inner(LBF[,j],LBF[,j],t))*LBF[,j]  }}}

LBF=LBF
LBF0=LBF


for (iter in 1:maxiter)
{

AR=append(c(runif(Non_dimX,-0.5,0.5)),c(runif(20-Non_dimX,0.5,0.9),decfac^(1:(lbnumber2-20))*runif(lbnumber2-20,-0.9,0.9)))

if(sparse==0){SD = append(c(rep(1,Non_dimX+sparsecut)),1*c(decfac^(1:(20-Non_dimX-sparsecut)),1*decfac^(20-Non_dimX-sparsecut+1) *(1:(lbnumber2-20))^(-2)))}
if(sparse==1){SD = append(c(rep(1,Non_dimX+sparsecut)),1*c(decfac1^(1:(20-Non_dimX-sparsecut)),1*decfac1^(20-Non_dimX-sparsecut+1) *(1:(lbnumber2-20))^(-2)))}
SD=sqrt(SD)
SIGU = sqrt(c(1*decfac^(0:(lbnumber2-1))))

bb=runif(lbnumber2,-1,1)*append(rep(1,Non_dimX+sparsecut),c(decfac2^(0:(lbnumber2-sparsecut-Non_dimX-1))))
LBFX=LBF; LBFY=LBF


## Data generation
score_sim=matrix(0,nrow=T_sim+burnin,ncol=lbnumber2)
for( j in 1:length(AR))
{
 aa=rnorm(1,0,SD[j])
 for (i in 2:(T_sim+burnin))
 {
  aa=append(aa,aa[i-1]*AR[j] + rnorm(1,0,SD[j]))
 }
 score_sim[,j] = aa
}
score_sim=score_sim[(burnin+1):(T_sim+burnin),]

Qhat_add <- t(score_sim[1:TSET[4],]-rowMeans(score_sim[1:TSET[4],]))%*%(score_sim[1:TSET[4],]-rowMeans(score_sim[1:TSET[4],]))/TSET[4]
eig_add3  <- sum(eigen(Qhat_add)$values)  
xeigensum3=sum(eig_add3)   
xeigensum=xeigensum3
magee= esignal  / 100 
strengdeg=(sqrt(magee*xeigensum*(1/eside)))
if(ranst==1){streng = strengdeg}
strengdeg2=(sqrt((1-magee)*xeigensum*(1/eside)))
streng2 = strengdeg2


for(j in 1:Non_dimX)
{
score_sim[,j]=cumsum(score_sim[,j])
}

### construct function time series ###
X_eig_fn=(LBFX)
X_sim_fn=matrix(0,nrow=nrow(LBFX),ncol=T_sim)
mux= LBFX[,1:10]%*%cbind(rnorm(10,0,1))
muy= LBFY[,1:10]%*%cbind(rnorm(10,0,1))
edex=c(sample(ea:eb,eside))

for(i in 1:T_sim)
{
ee=LBF[,edex]%*%cbind(rnorm(ncol(LBF[,edex]),0,1)) 
X_sim_fn[,i]=c(rbind(score_sim[i,])%*%t(X_eig_fn)) + mux +streng*ee 
}
X_sim_fn0=X_sim_fn


yscore_sim=matrix(0,nrow=T_sim,ncol=lbnumber2)


for (j in 1:lbnumber2)
{
yscore_sim[,j] = bb[j]*score_sim[,j] + rnorm(1,0,SIGU[j])
}

Y_eig_fn=(LBFY)
Y_sim_fn=matrix(0,nrow=nrow(LBFY),ncol=T_sim)
for(i in 1:T_sim)
{
ee=LBF[,edex]%*%cbind(rnorm(ncol(LBF[,edex]),0,1)) 
Y_sim_fn[,i]=c(rbind(yscore_sim[i,])%*%t(Y_eig_fn)) + muy +streng2*ee  
}

X_sim_fn0=X_sim_fn  
Y_sim_fn0=Y_sim_fn



for (TT_sim in TSET)
{
##temporal demeaning##
X_sim_fn = X_sim_fn0[,1:TT_sim]-rowMeans(X_sim_fn0[,1:TT_sim])
Y_sim_fn = Y_sim_fn0[,1:TT_sim]-rowMeans(Y_sim_fn0[,1:TT_sim])

######## estimation #################################
LBF=LBF0

  Xraw=t(X_sim_fn)
  T <- nrow(Xraw)
  X=t(Xraw)
  XX <- t(LBF[1:nt,])%*%X[1:nt,]*(t[2]-t[1])
  
  Yraw=t(Y_sim_fn)
  T <- nrow(Yraw)
  Y=t(Yraw)
  YY <- t(LBF[1:nt,])%*%Y[1:nt,]*(t[2]-t[1])
  
 
for (kap in c(1,0))
{ 
FC_X1=XX[,1:(T-kap)]
FC_X0=XX[,(kap+1):T]
FC_Y =YY[,(kap+1):T]

C_kap=FC_X0%*%t(FC_X1)/T
D_kap = t(C_kap)%*%C_kap

eval_D= eigen(D_kap)$values
evec_D= eigen(D_kap)$vectors

if (Kscheme==1){K_ind= max(Non_dimX + sum(eval_D[(Non_dimX+1):length(eval_D)]/sum(eval_D[(Non_dimX+1):length(eval_D)])>cuteigen*(TT_sim^(-bdw))),Non_dimX+1)}

Z1=t(t(FC_X1) %*% evec_D[,1:K_ind])
Z0=t(t(FC_X0) %*% evec_D[,1:K_ind])
Z1S=t(t(FC_X1) %*% evec_D[,(Non_dimX+1):K_ind])
Z0S=t(t(FC_X0) %*% evec_D[,(Non_dimX+1):K_ind])

C_kap_D=Z0%*%t(Z1)/T  
D_kap_inv=diag(eval_D[1:K_ind]^(-1))
CR_kap =  FC_Y %*% t(Z1)/T

fkap_N=CR_kap %*% C_kap_D%*%D_kap_inv %*% diag(c(rep(1,Non_dimX),rep(0,K_ind-Non_dimX)))

FC_resid=FC_Y
for (i in 1:ncol(FC_Y))
{
FC_resid[,i]=FC_Y[,i] - fkap_N %*% Z0[,i]
}


CR_kap2 = FC_resid %*% t(Z1)/T
{fkap_S=CR_kap2%*% C_kap_D%*%D_kap_inv %*% diag(c(rep(0,Non_dimX),rep(1,K_ind-Non_dimX)))}

## residuals from regression.
comp1=0
comp1N=0
atemindex=1:50
for ( j in 1:length(bb)) 
{
ateminput=rep(0,lbnumber2); ateminput[atemindex[j]]=1
atem=(LBF)%*%(fkap_N +fkap_S) %*% t(evec_D[,1:K_ind]) %*% cbind(ateminput)  
btem=Y_eig_fn[,j]*bb[j]
ctem=inner(atem-btem,atem-btem,t) 
comp1 = comp1+ ctem^2
if (j==Non_dimX){comp1N=comp1}
}
indexx=which(KAPP==kap)
if (TT_sim==TSET[1]){RESULT1[iter,indexx]=comp1; RESULT1N[iter,indexx]=comp1N; RESULT1S[iter,indexx]=comp1-comp1N; RESULTK1[iter,indexx]=K_ind}
if (TT_sim==TSET[2]){RESULT2[iter,indexx]=comp1; RESULT2N[iter,indexx]=comp1N; RESULT2S[iter,indexx]=comp1-comp1N; RESULTK2[iter,indexx]=K_ind}
if (TT_sim==TSET[3]){RESULT3[iter,indexx]=comp1; RESULT3N[iter,indexx]=comp1N; RESULT3S[iter,indexx]=comp1-comp1N; RESULTK3[iter,indexx]=K_ind}
if (TT_sim==TSET[4]){RESULT4[iter,indexx]=comp1; RESULT4N[iter,indexx]=comp1N; RESULT4S[iter,indexx]=comp1-comp1N; RESULTK4[iter,indexx]=K_ind}
}

}
incfac=1
if(iter%%10==0){print(round(c(sparse,geotail,decfac,decfac1,perturb,Non_dimX,magee,kap, bdw,incfac*colMeans(RESULT1[1:iter,]),incfac*colMeans(RESULT2[1:iter,]),incfac*colMeans(RESULT3[1:iter,]),incfac*colMeans(RESULT4[1:iter,])),digits=4))}
if(iter%%100==0){print(round(c(sparse,geotail,decfac,decfac1,perturb,Non_dimX,magee,kap,bdw, colMeans(RESULTK1[1:iter,]),colMeans(RESULTK2[1:iter,]),colMeans(RESULTK3[1:iter,]),incfac*colMeans(RESULTK4[1:iter,]),digits=4)))}
}


if (sparse==1){RESULT1_sp=RESULT1; RESULT1N_sp =RESULT1N; RESULT1S_sp=RESULT1S; RESULTK1_sp=RESULTK1; 
RESULT2_sp=RESULT2; RESULT2N_sp=RESULT2N; RESULT2S_sp=RESULT2S; RESULTK2_sp=RESULTK2; 
RESULT3_sp=RESULT3; RESULT3N_sp=RESULT3N; RESULT3S_sp=RESULT3S; RESULTK3_sp=RESULTK3; 
RESULT4_sp=RESULT4; RESULT4N_sp=RESULT4N; RESULT4S_sp=RESULT4S; RESULTK4_sp=RESULTK4}
}

print(round(c( incfac*colMeans(RESULT1[1:iter,]),incfac*colMeans(RESULT2[1:iter,]),incfac*colMeans(RESULT3[1:iter,]),incfac*colMeans(RESULT4[1:iter,])),digits=4))
print(round(c(incfac*colMeans(RESULT1_sp[1:iter,]),incfac*colMeans(RESULT2_sp[1:iter,]),incfac*colMeans(RESULT3_sp[1:iter,]),incfac*colMeans(RESULT4_sp[1:iter,])),digits=4))

###### Results ######
## kappa=0
print(round(c( mean(sqrt(RESULT1[1:iter,1])),mean(sqrt(RESULT2[1:iter,1])),mean(sqrt(RESULT3[1:iter,1])),mean(sqrt(RESULT4[1:iter,1])), mean(sqrt(RESULT1_sp[1:iter,1])),mean(sqrt(RESULT2_sp[1:iter,1])),mean(sqrt(RESULT3_sp[1:iter,1])),mean(sqrt(RESULT4_sp[1:iter,1])) ),digits=3))

## kappa=1
print(round(c( mean(sqrt(RESULT1[1:iter,2])),mean(sqrt(RESULT2[1:iter,2])),mean(sqrt(RESULT3[1:iter,2])),mean(sqrt(RESULT4[1:iter,2])), mean(sqrt(RESULT1_sp[1:iter,2])),mean(sqrt(RESULT2_sp[1:iter,2])),mean(sqrt(RESULT3_sp[1:iter,2])),mean(sqrt(RESULT4_sp[1:iter,2])) ),digits=3))
