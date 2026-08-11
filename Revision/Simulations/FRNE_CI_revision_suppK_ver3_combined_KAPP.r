#####################################################
######## PREPARATION ################################
#####################################################

# Set working directory (change as needed)
setwd("Path, e.g., C:/Users/R_code")  

# Load required library  
library(R.matlab)
library(sde)
library(robustbase)

inner = function(f,g,grid){
  h = f*g
  return(sum((0.5*h[1:(length(grid)-1)] + 0.5*h[2:(length(grid))])*(grid[2] - grid[1])))
}


Kscheme=1
Kval_cut = 0.75

for (esignal in c(0,50,100))    ## 0, 50, 100
{
######
   

for (Non_dimX in c(2,3)) {
for (rho in c(0,0.2,0.4)){

Scut=30
Non_dimY=Non_dimX
decfac=0.8
decfac1=0.2
decfac2=0.8
zetacut=9

# bdw list
bdw_list = c(0.25, 0.3, 0.35)
cuteigen=1

sparseset=c(1,0)
ea=1; eb=Non_dimX+1; eside=Non_dimX+1

for (sparse in sparseset) {

set.seed(12345)
TSET=c(100,200,400,800,1600)
Kscheme=1
K_fix=4
addfac=0; ranst=1;
perturb=1;geotail=18;sparsecut=7

T_sim = max(TSET) ; burnin=150 ;  
KAPP = c(0,1,2)  
maxiter=3000

# basic setup
res_total <- list()
for (bd_val in bdw_list) {
  bd_str <- as.character(bd_val)
  res_total[[bd_str]] <- list(
    RESULT1=matrix(0, maxiter, length(KAPP)), RESULT2=matrix(0, maxiter, length(KAPP)),
    RESULT3=matrix(0, maxiter, length(KAPP)), RESULT4=matrix(0, maxiter, length(KAPP)),
    RESULT5=matrix(0, maxiter, length(KAPP)),
    RESULTK1=matrix(0, maxiter, length(KAPP)), RESULTK2=matrix(0, maxiter, length(KAPP)),
    RESULTK3=matrix(0, maxiter, length(KAPP)), RESULTK4=matrix(0, maxiter, length(KAPP)),
    RESULTK5=matrix(0, maxiter, length(KAPP))
  )
}

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

for (iter in 1:maxiter) {

  AR=append(c(runif(Non_dimX,-0.5,0.5)),c(runif(20-Non_dimX,0.5,0.9),decfac^(1:(lbnumber2-20))*runif(lbnumber2-20,-0.9,0.9)))
  if(sparse==0){SD = append(c(rep(1,Non_dimX+sparsecut)),1*c(decfac^(1:(20-Non_dimX-sparsecut)),1*decfac^(20-Non_dimX-sparsecut+1) *(1:(lbnumber2-20))^(-2)))}
  if(sparse==1){SD = append(c(rep(1,Non_dimX+sparsecut)),1*c(decfac1^(1:(20-Non_dimX-sparsecut)),1*decfac1^(20-Non_dimX-sparsecut+1) *(1:(lbnumber2-20))^(-2)))}
  SD=sqrt(SD)
  SIGU = sqrt(c(1*decfac^(0:(lbnumber2-1))))
  
  bb=runif(lbnumber2,-1,1)*append(rep(1,Non_dimX+sparsecut),c(decfac2^(0:(lbnumber2-sparsecut-Non_dimX-1))))
  LBFX=LBF; LBFY=LBF

  ## Data generation
  R_corr <- diag(lbnumber2); block_size <- Non_dimX+2; if (rho != 0) {R_corr[1:block_size, 1:block_size] <- toeplitz(rho^(0:(block_size - 1)))}
  Sigma <- outer(SD, SD) * R_corr
  L <- t(chol(Sigma))
  Z <- matrix(rnorm((T_sim + burnin) * lbnumber2), nrow = lbnumber2, ncol = T_sim + burnin)
  E <- t(L %*% Z) 
  
  score_sim <- matrix(0, nrow = T_sim + burnin, ncol = lbnumber2)
  for (j in 1:lbnumber2) { score_sim[, j] <- as.numeric(stats::filter(E[, j], filter = AR[j], method = "recursive")) }
  score_sim = score_sim[(burnin + 1):(T_sim + burnin), ]

 # tmp_score <- score_sim[1:TSET[4], ] - rowMeans(score_sim[1:TSET[4], ])
   tmp_score <- sweep(score_sim[1:TSET[4], ], 2, colMeans(score_sim[1:TSET[4], ]), "-")
  Qhat_add <- crossprod(tmp_score) / TSET[4]
  eig_add3  <- sum(eigen(Qhat_add)$values)  
  xeigensum3=sum(eig_add3)   
  xeigensum=xeigensum3
  magee= (esignal  / 100) * ((esignal  / 100)<=1)
  strengdeg=(sqrt(magee*xeigensum*(1/eside)))
  if(ranst==1){streng = strengdeg}
  strengdeg2=(sqrt((1-magee)*xeigensum*(1/eside)))
  streng2 = strengdeg2

  for(j in 1:Non_dimX) { score_sim[,j]=cumsum(score_sim[,j]) }

  X_eig_fn=(LBFX)
  X_sim_fn=matrix(0,nrow=nrow(LBFX),ncol=T_sim)
  mux= LBFX[,1:10]%*%cbind(rnorm(10,0,1))
  muy= LBFY[,1:10]%*%cbind(rnorm(10,0,1))
  edex=c(sample(ea:eb,eside))
  
  EE_X <- LBF[,edex] %*% matrix(rnorm(length(edex) * T_sim, 0, 1), nrow = length(edex), ncol = T_sim)
  X_sim_fn <- (X_eig_fn %*% t(score_sim)) + as.vector(mux) + (streng * EE_X)

  yscore_sim=matrix(0,nrow=T_sim,ncol=lbnumber2)
  for (j in 1:lbnumber2) { yscore_sim[,j] = bb[j]*score_sim[,j] + rnorm(T_sim,0,SIGU[j]) }

  Y_eig_fn=(LBFY)
  EE_Y <- LBF[,edex] %*% matrix(rnorm(length(edex) * T_sim, 0, 1), nrow = length(edex), ncol = T_sim)
  Y_sim_fn <- (Y_eig_fn %*% t(yscore_sim)) + as.vector(muy) + (streng2 * EE_Y)

  X_sim_fn0=X_sim_fn  
  Y_sim_fn0=Y_sim_fn

  for (TT_sim in TSET) {
    X_sim_fn = X_sim_fn0[,1:TT_sim]-rowMeans(X_sim_fn0[,1:TT_sim])
    Y_sim_fn = Y_sim_fn0[,1:TT_sim]-rowMeans(Y_sim_fn0[,1:TT_sim])

    LBF=LBF0
    Xraw=t(X_sim_fn); T_len <- nrow(Xraw); X=t(Xraw); XX <- t(LBF[1:nt,])%*%X[1:nt,]*(t[2]-t[1])
    Yraw=t(Y_sim_fn); Y=t(Yraw); YY <- t(LBF[1:nt,])%*%Y[1:nt,]*(t[2]-t[1])
    
    for (kap in KAPP) { 
      FC_X1=XX[,1:(T_len-kap)]
      FC_X0=XX[,(kap+1):T_len]
      FC_Y =YY[,(kap+1):T_len]
      
      C_kap=FC_X0%*%t(FC_X1)/T_len
      D_kap = t(C_kap)%*%C_kap
      eval_D= eigen(D_kap)$values
      evec_D= eigen(D_kap)$vectors

      for (bd_val in bdw_list) {
        bd_str <- as.character(bd_val)

   # 
        if (Kscheme==1){
          eval_DP=eval_D[(Non_dimX+1):length(eval_D)]*(eval_D[(Non_dimX+1):length(eval_D)]>0)
          c_val <- (1 - Kval_cut) / (TSET[1]^(-bd_val)) 
          Kval_cut_adj <- 1 - c_val * (TT_sim^(-bd_val))
          eval_temp=which(cumsum(eval_DP)/sum(eval_DP) >= Kval_cut_adj)
          if (length(eval_temp)==0){K_ind=1+Non_dimX}else{K_ind= max(Non_dimX+min(eval_temp),1+Non_dimX)}
        } else {
          eval_DP=eval_D[(Non_dimX+1):length(eval_D)]*(eval_D[(Non_dimX+1):length(eval_D)]>0)
          eval_temp=which(cumsum(eval_DP)/sum(eval_DP) >= 1-(TT_sim^(-bd_val))/(TSET[1]^(-bd_val)/(1-Kval_cut))) 
          if (length(eval_temp)==0){K_ind=1+Non_dimX}else{K_ind= Non_dimX+min(eval_temp)}
        }
		
		
        Z1=t(t(FC_X1) %*% evec_D[,1:K_ind])
        Z0=t(t(FC_X0) %*% evec_D[,1:K_ind])
        Z1S=t(t(FC_X1) %*% evec_D[,(Non_dimX+1):K_ind])
        Z0S=t(t(FC_X0) %*% evec_D[,(Non_dimX+1):K_ind])
        
        C_kap_D=Z0%*%t(Z1)/T_len  
        D_kap_inv=diag(eval_D[1:K_ind]^(-1))
        CR_kap =  FC_Y %*% t(Z1)/T_len
        fkap_N=CR_kap %*% C_kap_D%*%D_kap_inv %*% diag(c(rep(1,Non_dimX),rep(0,K_ind-Non_dimX)))
        
        FC_resid = FC_Y - (fkap_N %*% Z0)
        CR_kap2 = FC_resid %*% t(Z1)/T_len
        fkap_S=CR_kap2%*% C_kap_D%*%D_kap_inv %*% diag(c(rep(0,Non_dimX),rep(1,K_ind-Non_dimX)))

        ## CI inference
        RESID <- FC_Y - ((fkap_N + fkap_S) %*% t(evec_D[, 1:K_ind]) %*% FC_X0)

        if (K_ind!=(Non_dimX+1)){SD_kap_inv=diag(eval_D[(Non_dimX+1):K_ind]^(-1))}else{SD_kap_inv=(eval_D[(Non_dimX+1):K_ind]^(-1))}  
        SC_kap_D=Z0S%*%t(Z1S)/T_len    
        SC_D=Z0S%*%t(Z0S)/T_len  
        P_KS=t(evec_D[,(Non_dimX+1):K_ind])
        C_resid=(RESID)%*%t(RESID)/T_len

        ateminput=rep(0,lbnumber2); ateminput=c(1^(1:zetacut),1*(1:(50-(zetacut)))^(-2))  
        pert=LBFX[,1:length(ateminput)]%*%cbind(ateminput)
        atem=(LBF)%*%(fkap_N +fkap_S) %*% t(evec_D[,1:K_ind]) %*% (t(LBF)%*%pert)*(t[2]-t[1])
        btem=LBFY%*%cbind(ateminput*bb)

        psi=LBFY[,1]
        ctem=inner(atem-btem,psi,t) 

        theta_raw=(SD_kap_inv%*%t(SC_kap_D)%*%SC_D%*%SC_kap_D%*%SD_kap_inv%*%P_KS%*%ateminput) * (P_KS%*%ateminput) 
        theta=sum(theta_raw)
        
        IIF=t(LBF)%*%psi*(t[2]-t[1])
        comp1 = ctem - 1.96*sqrt(theta*t(IIF)%*%C_resid%*%IIF/T_len)
        comp2 = ctem + 1.96*sqrt(theta*t(IIF)%*%C_resid%*%IIF/T_len)

        result=1*(0>comp1 & 0<comp2)
        
        indexx=which(KAPP==kap)
        if (TT_sim==TSET[1]){res_total[[bd_str]]$RESULT1[iter,indexx]=result; res_total[[bd_str]]$RESULTK1[iter,indexx]=K_ind}
        if (TT_sim==TSET[2]){res_total[[bd_str]]$RESULT2[iter,indexx]=result; res_total[[bd_str]]$RESULTK2[iter,indexx]=K_ind}
        if (TT_sim==TSET[3]){res_total[[bd_str]]$RESULT3[iter,indexx]=result; res_total[[bd_str]]$RESULTK3[iter,indexx]=K_ind}
        if (TT_sim==TSET[4]){res_total[[bd_str]]$RESULT4[iter,indexx]=result; res_total[[bd_str]]$RESULTK4[iter,indexx]=K_ind}
        if (TT_sim==TSET[5]){res_total[[bd_str]]$RESULT5[iter,indexx]=result; res_total[[bd_str]]$RESULTK5[iter,indexx]=K_ind}
      } # End bdw loop
    } # End kap loop
  } # End TT_sim loop

 incfac=1
  if(iter%%100==0) {
    cat(sprintf("\n=== Iteration: %d (CI Coverage for all KAPPA) ===\n", iter))
    for (bd_val in bdw_list) {
      bd_str <- as.character(bd_val)
     
      res_means <- c(
        colMeans(res_total[[bd_str]]$RESULT1[1:iter, , drop=FALSE]),
        colMeans(res_total[[bd_str]]$RESULT2[1:iter, , drop=FALSE]),
        colMeans(res_total[[bd_str]]$RESULT3[1:iter, , drop=FALSE]),
        colMeans(res_total[[bd_str]]$RESULT4[1:iter, , drop=FALSE]),
        colMeans(res_total[[bd_str]]$RESULT5[1:iter, , drop=FALSE])
      )
      
      cov_vals <- round(c(sparse, rho, Non_dimX, esignal, bd_val, incfac * res_means), digits=4)
      
      if (iter %in% c(100, 1000, 2000)) {
        k_labels <- paste0("k", KAPP)
      
        t_labels <- as.vector(outer(k_labels, paste0("T", 1:5), function(x, y) paste(y, x, sep="_")))
        t_labels_fixed <- as.vector(t(outer(paste0("T", 1:5), k_labels, paste, sep="_")))
        
        names(cov_vals) <- c("Sp", "Rho", "dnX", "eSig", "bdw", t_labels_fixed)
      } else {
        names(cov_vals) <- NULL
      }
      print(cov_vals)
    }
  }
  
} # End iter loop

if (sparse==1){ res_total_sp <- res_total }
} # End sparse loop

for (bd_val in bdw_list) {
  bd_str <- as.character(bd_val)
  RESULT1 = res_total[[bd_str]]$RESULT1; RESULT2 = res_total[[bd_str]]$RESULT2; RESULT3 = res_total[[bd_str]]$RESULT3; RESULT4 = res_total[[bd_str]]$RESULT4; RESULT5 = res_total[[bd_str]]$RESULT5
  RESULTK1 = res_total[[bd_str]]$RESULTK1; RESULTK2 = res_total[[bd_str]]$RESULTK2; RESULTK3 = res_total[[bd_str]]$RESULTK3; RESULTK4 = res_total[[bd_str]]$RESULTK4; RESULTK5 = res_total[[bd_str]]$RESULTK5
  RESULT1_sp = res_total_sp[[bd_str]]$RESULT1; RESULT2_sp = res_total_sp[[bd_str]]$RESULT2; RESULT3_sp = res_total_sp[[bd_str]]$RESULT3; RESULT4_sp = res_total_sp[[bd_str]]$RESULT4; RESULT5_sp = res_total_sp[[bd_str]]$RESULT5
  RESULTK1_sp = res_total_sp[[bd_str]]$RESULTK1; RESULTK2_sp = res_total_sp[[bd_str]]$RESULTK2; RESULTK3_sp = res_total_sp[[bd_str]]$RESULTK3; RESULTK4_sp = res_total_sp[[bd_str]]$RESULTK4; RESULTK5_sp = res_total_sp[[bd_str]]$RESULTK5

  save_filename <- paste0("results/FRNENEW_KsuppadjBLOCKKAPP_esignal", esignal, "_rho", rho, "_Kscheme", Kscheme, "_Kcutmin", Kval_cut, "_dn", Non_dimX, "_CI_bdw", bd_val, "_defac1", decfac1, ".RData")
  save(file=save_filename, RESULT1,RESULT2,RESULT3,RESULT4,RESULT5,RESULT1_sp,RESULT2_sp,RESULT3_sp,RESULT4_sp,RESULT5_sp,RESULTK1,RESULTK2,RESULTK3,RESULTK4,RESULTK5,RESULTK1_sp,RESULTK2_sp,RESULTK3_sp,RESULTK4_sp,RESULTK5_sp)
}
}
}


} # End esignal loop


