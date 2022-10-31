exp.crit.all <- function(d,B)
{
  crit=c()
  ord <- order(d)
  locs1 <- which(Monitoring_sites$Site_No %in% d)
  locs <- locs1[ord]
  
  distM=dist_mat[locs,locs]
  Loc_sel <- ST_uns[locs,]
  X_data <- X.all[locs,]

 	for(j in 1:B){
 	  Y <- Res_Clcpp(xdata=X_data,distM=distM,para=theta_int[j,],Z1=t(Z_res[1:nrow(X_data),1]),Z2=t(Z_res[1:nrow(X_data),2]),Z3=Z_res[1:nrow(X_data),3],v2=Unif_set[1:nrow(X_data)],r01=r01,r02=r02)
 	  crit.out <- combine_ute_all(d=d,X=X_data,distM=distM,Y=Y,theta=theta_int[j,])
 	  crit.dual <- crit.out[[3]]
 	  crit.kld <- crit.out[[4]] 
 	  crit.pred.App <- crit.out[[5]] 
    crit.pred <- crit.out[[6]]
    out.all<- c(crit.dual,crit.kld,crit.pred.App,crit.pred)
    crit <- c(crit,out.all)
  }
  return(crit)
}
