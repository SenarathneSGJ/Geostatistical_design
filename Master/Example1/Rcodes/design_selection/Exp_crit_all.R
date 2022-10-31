exp.crit.all <- function(d,B)
{
  distM=rdist(d)
  X1 <- ((d[,2]+.1)*5)-((d[,1]+.2)*5)
  X2 <- (d[,2]+.1)*(d[,1]+.2)
  X <- data.frame(X1,X2)

  crit=foreach(j = 1:B,.packages= c("mvtnorm","fields","corpcor","nlme","Matrix","matrixcalc","MaternEx1"),.errorhandling = 'stop',.export= ls(globalenv()),.combine = c) %dopar%  
	{
     Y = Res_Clcpp(xdata=as.matrix(X),distM=distM,para=theta_int[j,],Z1=t(Z_res[1:nrow(X),1]),Z2=t(Z_res[1:nrow(X),2]),Z3=Z_res[1:nrow(X),3],v2=Unif_set[1:nrow(X)],r01=r01,r02=r02)
     crit.out <- combine_ute(d,X,Y,theta=theta_int[j,],distM)
     crit.dual <- crit.out[[3]]
     crit.kld <- crit.out[[4]]
     crit.pred <- crit.out[[5]]
     c(crit.kld,crit.pred,crit.dual) 
  }
  return(crit)
}
