exp.crit <- function(d,B)
{
  crit <- rep(0,B)
  crit_kld <- rep(0,B)
  crit_pred <- rep(0,B)
  ord <- order(d)
  locs1 <- which(Monitoring_sites$Site_No %in% d)
  locs <- locs1[ord]
  
  distM=dist_mat[locs,locs]
  X_data <- X.all[locs,]
  
  for (j in 1:B)
  {
    Y <- Res_Clcpp(xdata=X_data,distM=distM,para=theta_int[j,],Z1=t(Z_res[1:nrow(X_data),1]),Z2=t(Z_res[1:nrow(X_data),2]),Z3=Z_res[1:nrow(X_data),3],v2=Unif_set[1:nrow(X_data)],r01=r01,r02=r02)
    crit.out <- combine_ute(d=d,X=X_data,Y=Y,theta=theta_int[j,],distM=distM)
    crit[j] <- crit.out[[3]] 
    crit_kld[j] <- crit.out[[4]] 
    crit_pred[j] <- crit.out[[5]]
  } 
  Avg.crit <- mean(crit,na.rm=TRUE)
  Avg.crit.kld <- mean(crit_kld,na.rm=TRUE)
  Avg.crit.pred <- mean(crit_pred,na.rm=TRUE)
  
  return(c(Avg.crit,Avg.crit.kld,Avg.crit.pred))
}
