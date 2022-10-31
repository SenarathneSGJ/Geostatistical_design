exp.crit.var <- function(d,B)
{
  distM=rdist(d)
  X1 <- ((d[,2]+.1)*5)-((d[,1]+.2)*5)
  X2 <- (d[,2]+.1)*(d[,1]+.2)
  X <- data.frame(X1,X2)
  
  Y1.post.pred=matrix(ncol=nrow(Unsamp_X),nrow=B)
  Y2.post.pred=matrix(ncol=nrow(Unsamp_X),nrow=B)
  post.var=matrix(ncol=ncol(prior),nrow=B)
  
  for(j in 1:B)
	{
     Y = Res_Clcpp(xdata=as.matrix(X),distM=distM,para=theta_int[j,],Z1=t(Z_res[1:nrow(X),1]),Z2=t(Z_res[1:nrow(X),2]),Z3=Z_res[1:nrow(X),3],v2=Unif_set[1:nrow(X)],r01=r01,r02=r02)
     crit.out <- Post_pred(d,X,Y,theta=theta_int[j,],distM)
     post.var[j,] <- diag(crit.out[[2]])
     Y1.post.pred[j,] <- apply(crit.out[[3]],1,var,na.rm=T)
     Y2.post.pred[j,] <- log(apply(crit.out[[4]],1,var,na.rm=T))
    
  }
  
  Y1.post <- apply(Y1.post.pred,2,mean,na.rm=T)
  Y2.post <- apply(Y2.post.pred,2,mean,na.rm=T)
  return(list(post.var,Y1.post,Y2.post))
}
