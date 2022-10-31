Pred_ute <- function(X_unsamp,distM,post_samp1,post_samp2,r01,r02)
{
  K=nrow(post_samp1)
  R=nrow(post_samp2)
  Pz <- rep(0,K)
  for(m in 1:K){
    Ym <- Res_Clcpp(xdata=as.matrix(X_unsamp[,4:5]),distM=distM,para=post_samp1[j,],Z1=t(rnorm(nrow(X_unsamp))),Z2=t(rnorm(nrow(X_unsamp))),Z3=rnorm(nrow(X_unsamp)),v2=Unif_set[1:nrow(X_unsamp)],r01=r01,r02=r02)
    log_z <- rep(0,R)
    for(r in 1:R){
     log_z[r]<-sum(Log_likelihoodCpp(as.matrix(X_unsamp[,4:5]),as.matrix(Ym),distM,post_samp2[r,],Z=t(Z1[1:nrow(Ym),(1:B2)]),r01,r02))
    }
    Pz[m]<-mean(exp(log_z),na.omit=TRUE) 
  }
  pred_ute <- mean(log(Pz),na.omit=TRUE)
  pred_ute[is.na(pred_ute) | pred_ute < (-10000)]=(-10000)
  
  return(pred_ute)
}