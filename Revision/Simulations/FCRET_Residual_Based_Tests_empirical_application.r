# Required R packages (need to be installed first)
library(fda); library(geigen); library(sandwich); library("readxl"); library(GMCM); library(psych);

# Set working directory (change as needed)
setwd("Path, e.g., C:/Users/R_code")  

 
# Basic parameter setup 
Case = 7
if (Case==7){# Loading Residuals represented with Fourier basis
  Grid <- read_excel("Residual_Kappa0.xlsx", sheet = 1, range = "A1:A200", col_names = FALSE,col_types = "numeric"); nage <- unlist(Grid, use.names = FALSE);
  dmat0 <- read_excel("Residual_Kappa0.xlsx", sheet = 1, range = "A1:BQ200", col_names = FALSE, col_types = "numeric"); dmat <- t(as.matrix(dmat0));
} 

# Define functions 
# numerical integration over grid
inner = function(f,g,grid){
  h = f*g
  return(sum((0.5*h[1:(length(grid)-1)] + 0.5*h[2:(length(grid))])*(grid[2] - grid[1])))
}

# Data load and basic parameters
#1 parameter setting 
dmat <- t(as.matrix(dmat0)); nt = 101 ; t = (0:(nt-1))/(nt-1);lbnumber2=200; bnumbasis=102; lbnumber22=140;nobs=70
numbasis = lbnumber2 # Number of Bspline basis functions to represent functional observations. 
linear = 0    # 0 : with intercept only, 1 : with linear trend
 
# 3. Generating Fourier Basis 
LBF = matrix(NA, nrow = nt, ncol = lbnumber2)
for (i in 1:(lbnumber2 / 2)) {
  LBF[, 2*i - 1] = sqrt(2) * sin(2 * pi * i * t)
  LBF[, 2*i]     = sqrt(2) * cos(2 * pi * i * t)
}
# 4. Modified Gram-Schmidt
lb = LBF  
for (i in 1:lbnumber2) {  current_norm = sqrt(inner(lb[, i], lb[, i], t));   lb[, i] = lb[, i] / current_norm
  if (i < lbnumber2) {    for (j in (i + 1):lbnumber2) {      projection_scale = inner(lb[, j], lb[, i], t);       lb[, j] = lb[, j] - projection_scale * lb[, i] } }}
LBF = lb; dmat=LBF[,1:lbnumber22]%*%t(dmat[,1:lbnumber22]);

#2 Data Preprocessing, Functional PCA
x_mat=dmat
x_mat=x_mat
basis_fn = create.bspline.basis(rangeval = c(0,1),nbasis = bnumbasis,norder=3)

  xx_mat=x_mat  # temporally demeaned functional data
  
  fd_xx =  (Data2fd(y=xx_mat,argvals = t,basisobj = basis_fn))# 31 basis representation of partial sum of temporally demeaned functional data
  hkmat = t(fd_xx$coefs)  # 31 basis representation of partial sum of temporally demeaned functional data
  hhtau = eigen(crossprod(hkmat),symmetric = TRUE) 
  #Extract Eigenvalue of sample covariance operator of 31 basis representation of partial sum of demeaned functional data
  
  fd_x =  (Data2fd(y=x_mat,argvals = t,basisobj = basis_fn))  # 31 basis representation of raw functional data
  
  h2kmat= t(fd_x$coefs ) # 31 basis representation of raw functional data
 
 for(jjj in 3:1)
 {
  fd_z = t(h2kmat%*%hhtau$vectors[,1:bnumbasis]) 
  fd_z = fd_z 
  fd_z=fd_z
  domeigen=eigen(crossprod(t(fd_z)))$vectors[,1:jjj]
  fd_z=t(t(fd_z)%*%cbind(domeigen))
  fd_z=fd_z
  #3 Test statistics 
   kmat = t(fd_z);   cmat = t(kmat)%*%kmat;   smat = apply(kmat,2,cumsum)
   smat = t(smat) %*% smat
   tau = ((nobs^2)*geigen(cmat,smat,symmetric = TRUE,only.values = TRUE)$values)
   maxtau=max(tau)
   sumtau=sum(tau)
 print(sumtau)
   }
  
  