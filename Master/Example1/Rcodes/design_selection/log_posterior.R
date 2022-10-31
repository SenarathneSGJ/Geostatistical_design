log_post <- function(X,theta,Y,distM)
  {
    log.prior<- dmvnorm(x=theta,mean=mu_prior,sigma=Sigma_prior,log=TRUE)
    #log.like <- LogLike_Cl3(xdata=X,y=Y,para=theta,B=B2,distM=distM)
    log.like<- sum(Log_likelihoodCpp(xdata=as.matrix(X),y=as.matrix(Y),distM=distM,para=theta,Z=t(Z1[1:nrow(Y),(1:B2)]),r01,r02)) 
    
    Neg_log_post=-1*(log.prior + log.like)
    Neg_log_post[is.na(Neg_log_post) | Neg_log_post==Inf | Neg_log_post == -Inf]= 1e10
    #print(Neg_log_post)
    
    return(Neg_log_post)
  }