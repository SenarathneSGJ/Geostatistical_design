log_post <- function(X,theta,Y,distM)
  {
    log.prior<- dmvnorm(x=theta,mean=mu_prior,sigma=Sigma_prior,log=TRUE)
    log.like <- sum(Log_likelihoodCpp(xdata=X,y=as.matrix(Y),distM=distM,para=theta,Z=t(Z1[1:nrow(X),]),r01=r01,r02=r02))
    
    Neg_log_post=-1*(log.prior + log.like)
    
    return(Neg_log_post)
  }