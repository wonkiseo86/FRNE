# Required R packages (need to be installed first)
library(fda); library(geigen); library(sandwich); library("readxl"); library(GMCM); library(psych);

setwd("path")  

Case = 2

if (Case == 1) { # Loading GRP Gr FE2 data for 1951-2019 
  Grid <- read_excel("Data_FCRGT.xlsx", sheet = 1, range = "B1:GQ1", col_names = FALSE,col_types = "numeric"); nage <- unlist(Grid, use.names = FALSE);
  dmat0 <- read_excel("Data_FCRGT.xlsx", sheet = 1, range = "B1:GQ70", col_names = TRUE, col_types = "numeric"); dmat <- t(as.matrix(dmat0)); 
} else if (Case == 2) { # Loading Land Surface Temperature data for 1951-2019 
  Grid <- read_excel("Data_FCRGT.xlsx", sheet = 2, range = "B1:CRB1", col_names = FALSE,col_types = "numeric"); nage <- unlist(Grid, use.names = FALSE);
  dmat0 <- read_excel("Data_FCRGT.xlsx", sheet = 2, range = "B1:CRB70", col_names = TRUE, col_types = "numeric"); dmat <- t(as.matrix(dmat0)); # GTEMP f  
} 

# Cointegration Rank tests 
#1 parameter setting 
#Parameter restriction : adddim+smax <= numbasis, the below default setting is to replicate Table 12.
adddim = 2    # (ell-s0) : additional number of eigenfunctions used to construct projection.   
smax = 5      # staring value (s_max) stated in the top hypothesis
numbasis = 50 # Number of Bspline basis functions to represent functional observations. 
linear = 0    # 0 : with intercept only, 1 : with linear trend

#2 Data Preprocessing, Functional PCA
t=nage ; nobs=ncol(dmat); nrows =nrow(dmat); # x_mat = dmat;
if (Case == 1 || Case == 2 || Case == 3 || Case == 4 || Case == 5 || Case == 6) {
  x_mat = matrix(NA, nrow=nrows, ncol=nobs);
for (iter in 1:nobs){
  x_mat[,iter] = log(dmat[,iter]/exp(mean(log(dmat[,iter]))));  
}
} else if (Case == 7 || Case == 8 || Case == 9) {
  x_mat = dmat;
}    

basis_fn = create.bspline.basis(rangeval = c(min(nage),max(nage)),nbasis = numbasis)

if (linear == 0)
{
  xx_mat=x_mat-rowMeans(x_mat) # temporally demeaned functional data
   
  fd_xx =  (Data2fd(y=xx_mat,argvals = t,basisobj = basis_fn))
  hkmat = t(fd_xx$coefs)  
  hhtau = eigen(crossprod(hkmat),symmetric = TRUE) 
  
  fd_x =  (Data2fd(y=x_mat,argvals = t,basisobj = basis_fn))  
  
  h2kmat= t(fd_x$coefs -rowMeans(fd_x$coefs)) 
  fd_z = t(h2kmat%*%hhtau$vectors[,1:numbasis])
  fd_z = fd_z -rowMeans(fd_z) 
  
  #3 Test statistics 
  ind = rep(5,length(1:smax + adddim))
  test_vr = rep(NA,smax)
  test_vr_pca = rep(NA,smax)
  for (initer in 1:smax){
    index = ind[initer]
    fd_zz = as.matrix(fd_z[1:index,])
    kmat = t(fd_zz)
    cmat = t(kmat)%*%kmat
    smat = apply(kmat,2,cumsum)
    smat = t(smat) %*% smat
    tau = (nobs^2)*geigen(cmat,smat,symmetric = TRUE,only.values = TRUE)$values
    test_vr[initer] = sum(tau[1:initer])
  }

  
  Aresult=NULL
  for(i in 1:smax)
  {
    Aresult=rbind(Aresult,paste("VR test,", "s0 =", i , " : ", round(test_vr[i],digits=2)))
  }
  print(Aresult)
  
} else{
  if (linear == 1)
  {
    trend=1:nobs
    const=rep(1,length=nobs)
    VV=cbind(const,trend)
    
    u_mat=t(x_mat)-VV%*%solve(t(VV)%*%VV)%*%t(VV)%*%t(x_mat)
    xx_mat=t(u_mat)
    xx_mat= t(apply(t(xx_mat),2,cumsum))  
    
    fd_z = xx_mat
    
    ind = rep(5,length(1:smax + adddim))
    test_vr = rep(NA,smax)
    test_vr_pca = rep(NA,smax)
    for (initer in 1:smax){
      index = ind[initer]
      fd_zz = as.matrix(fd_z[1:index,])
      kmat = t(fd_zz)
      cmat = t(kmat)%*%kmat
      smat = apply(kmat,2,cumsum)
      smat = t(smat) %*% smat
      tau = (nobs^2)*geigen(cmat,smat,symmetric = TRUE,only.values = TRUE)$values
      test_vr[initer] = sum(tau[1:initer])
    }
    
    Aresult=NULL
    for(i in 1:smax)
    {
      Aresult=rbind(Aresult,paste("VR test,", "s0 =", i , " : ", round(test_vr[i],digits=2)))
    }
    print(Aresult)
  }
}

