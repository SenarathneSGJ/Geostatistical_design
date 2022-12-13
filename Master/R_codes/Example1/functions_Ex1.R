# Combine_ute function returns the posterior mean, posterior covariance matrix and the three utility values
# for a given design=d, values for independent variables=X, data=Y, parameter vector=theta, and distance matrix=distM
combine_ute <- function(d,X,Y,theta,distM){
  
  Sol <- Laplace_approx(X,Y,theta,distM)
  mu_post <- Sol$par
  hessian1 <- optimHess(par=mu_post,distM=distM,fn=log_post,X=X,Y=Y)
  
  if(Sol$convergence==1)
  {
    print("error1")
  }
  
  if(!is.symmetric.matrix(hessian1)){
    print("error2")
    hess1 <- hessian1
    hess1[lower.tri(hess1,diag=F)] <- 0
    hessian1 <- hess1+t(hess1)
    if(!is.positive.semi.definite(hessian1,tol=1e-16) | det(hessian1)==0) {
      Sigma_post <- Sigma_prior
      log_det_post <- log(det(Sigma_post))
      kld.out <- 0.5*(trace_mat(iSigma_prior%*%Sigma_post) + t(mu_prior-mu_post)%*%iSigma_prior%*%(mu_prior-mu_post)-
					num_par + log_det_prior-log_det_post)
      post_samp1 <- rmvnorm(200,mu_post,Sigma_post)
      pred.out <- Pred_utecpp(X_unsamp=as.matrix(Unsamp_X[,4:5]),post1=as.matrix(post_samp1),distM=as.matrix(distM_uns),
					v2=Unif_set[1:nrow(Unsamp_X)],r01,r02)
      #pred.out <- Pred_ute.App(X_unsamp=Unsamp_X,distM=distM_uns,post_samp1=post_samp1,r01=r01,r02=r02)
      det.out <- kld.out+pred.out  
    }else{
      Sigma_post <- solve(hessian1)
      log_det_post <- log(det(Sigma_post))
      kld.out <- 0.5*(trace_mat(iSigma_prior%*%Sigma_post) + t(mu_prior-mu_post)%*%iSigma_prior%*%(mu_prior-mu_post)-
					num_par + log_det_prior-log_det_post)
      post_samp1 <- rmvnorm(200,mu_post,Sigma_post)
      pred.out <- Pred_utecpp(X_unsamp=as.matrix(Unsamp_X[,4:5]),post1=as.matrix(post_samp1),distM=as.matrix(distM_uns),
					v2=Unif_set[1:nrow(Unsamp_X)],r01,r02)
      #pred.out <- Pred_ute.App(X_unsamp=Unsamp_X,distM=distM_uns,post_samp1=post_samp1,r01=r01,r02=r02)
      det.out <- kld.out+pred.out  
    }
    
  }else if(!is.positive.semi.definite(hessian1,tol=1e-16) | det(hessian1)==0) {
    print("error3")
    Sigma_post <- Sigma_prior
    log_det_post <- log(det(Sigma_post))
    kld.out <- 0.5*(trace_mat(iSigma_prior%*%Sigma_post) + t(mu_prior-mu_post)%*%iSigma_prior%*%(mu_prior-mu_post)-
				num_par + log_det_prior-log_det_post)
    post_samp1 <- rmvnorm(200,mu_post,Sigma_post)
    pred.out <- Pred_utecpp(X_unsamp=as.matrix(Unsamp_X[,4:5]),post1=as.matrix(post_samp1),distM=as.matrix(distM_uns),
					v2=Unif_set[1:nrow(Unsamp_X)],r01,r02)
    #pred.out <- Pred_ute.App(X_unsamp=Unsamp_X,distM=distM_uns,post_samp1=post_samp1,r01=r01,r02=r02)
    det.out <- kld.out+pred.out   
  }else{
    Sigma_post <- solve(hessian1)
    log_det_post <- log(det(Sigma_post))
    kld.out <- 0.5*(trace_mat(iSigma_prior%*%Sigma_post) + t(mu_prior-mu_post)%*%iSigma_prior%*%(mu_prior-mu_post)-
				num_par + log_det_prior-log_det_post)
    post_samp1 <- (rmvnorm(200,mu_post,Sigma_post))
    pred.out <- Pred_utecpp(X_unsamp=as.matrix(Unsamp_X[,4:5]),post1=as.matrix(post_samp1),distM=as.matrix(distM_uns),
					v2=Unif_set[1:nrow(Unsamp_X)],r01,r02)
    #pred.out <- Pred_ute.App(X_unsamp=Unsamp_X,distM=distM_uns,post_samp1=post_samp1,r01=r01,r02=r02)
    det.out <- kld.out+pred.out
  }
  
  return(list(mu_post,Sigma_post,det.out,kld.out,pred.out)) 
}

###################
# combine_ute_both function returns the posterior mean, posterior covariance matrix, dual purpose utility value, estimation 
# utility value and the prediction utility values based on the two approximations and time taken to run each approxiation
# for a given design=d, values for independent variables=X, data=Y, parameter vector=theta, and the distance matrix=distM
combine_ute_both <- function(d,X,Y,theta,distM){
  
  Sol <- Laplace_approx(X,Y,theta,distM)
  mu_post <- Sol$par
  hessian1 <- optimHess(par=mu_post,distM=distM,fn=log_post,X=X,Y=Y)
  
  if(Sol$convergence==1)
  {
    print("error1")
  }
  
  if(!is.symmetric.matrix(hessian1)){
    print("error2")
    hess1 <- hessian1
    hess1[lower.tri(hess1,diag=F)] <- 0
    hessian1 <- hess1+t(hess1)
    if (!is.positive.semi.definite(hessian1,tol=1e-16) | det(hessian1)==0) {
      Sigma_post <- Sigma_prior
      kld.out <- NA
      pred.out <- NA
      pred.out2 <- NA
      det.out <- NA 
    }else{
      Sigma_post <- solve(hessian1)
      log_det_post <- log(det(Sigma_post))
      kld.out <- 0.5*(trace_mat(iSigma_prior%*%Sigma_post) + t(mu_prior-mu_post)%*%iSigma_prior%*%(mu_prior-mu_post)-
					num_par + log_det_prior-log_det_post)
      post_samp <- rmvnorm(100,mu_post,Sigma_post)
      t1 <- proc.time()
      pred.out <- Pred_ute.App(X_unsamp=Unsamp_X,distM=distM_uns,post_samp1=post_samp,r01=r01,r02=r02)
      t2 <- proc.time()
      time1 <- t2-t1
      
      post_samp2 <- rmvnorm(100,mu_post,Sigma_post)
      t3 <- proc.time()
      pred.out2 <- Pred_ute(X_unsamp=Unsamp_X,distM=distM_uns,post_samp1=post_samp,post_samp2=post_samp2,r01=r01,r02=r02)
      t4 <- proc.time()
      time2 <- t4-t3
      
      det.out <- kld.out+pred.out  
    }
    
  }else if (!is.positive.semi.definite(hessian1,tol=1e-16) | det(hessian1)==0) {
    print("error3")
    Sigma_post <- Sigma_prior
    kld.out <- NA
    pred.out <- NA
    pred.out2 <- NA
    det.out <- NA 
  }else{
    Sigma_post <- solve(hessian1)
    log_det_post <- log(det(Sigma_post))
    kld.out <- 0.5*(trace_mat(iSigma_prior%*%Sigma_post) + t(mu_prior-mu_post)%*%iSigma_prior%*%(mu_prior-mu_post)-
				num_par + log_det_prior-log_det_post)
    post_samp <- rmvnorm(100,mu_post,Sigma_post)
    t1 <- proc.time()
    pred.out <- Pred_ute.App(X_unsamp=Unsamp_X,distM=distM_uns,post_samp1=post_samp,r01=r01,r02=r02)
    t2=proc.time()
    time1 <- t2-t1
    
    post_samp2 <- rmvnorm(100,mu_post,Sigma_post)
    t3 <- proc.time()
    pred.out2 <- Pred_ute(X_unsamp=Unsamp_X,distM=distM_uns,post_samp1=post_samp,post_samp2=post_samp2,r01=r01,r02=r02)
    t4 <- proc.time()
    time2 <- t4-t3
    det.out <- kld.out+pred.out
  }
  
  return(list(mu_post,Sigma_post,det.out,kld.out,pred.out,pred.out2,time1,time2)) 
}

###################
# exp.crit function returns the utility values (u(d,y)) based on the selected utility function for a given design=d and 
# theta matrix=theta_int, where input B is the number of utility values cosidered to obtain the expected utilities.
exp.crit <- function(d,B)
{
  distM <- rdist(d)
  X1 <- ((d[,2]+.1)*5)-((d[,1]+.2)*5)
  X2 <- (d[,2]+.1)*(d[,1]+.2)
  X <- data.frame(X1,X2)
  
  crit <- foreach(j=1:B,.packages=c("mvtnorm","fields","corpcor","nlme","Matrix","matrixcalc","MaternEx1"),
			.errorhandling='stop',.export=ls(globalenv()),.combine=c) %dopar%  
    {
      Y <- Res_Clcpp(xdata=as.matrix(X),distM=distM,para=theta_int[j,],Z1=t(rnorm(1:nrow(X))),
			Z2=t(rnorm(1:nrow(X))),Z3=(rnorm(1:nrow(X))),v2=Unif_set[1:nrow(X)],r01=r01,r02=r02)
      crit.out <- combine_ute(d,X,Y,theta=theta_int[j,],distM)
      crit.out[[utility+2]] 
    }
  return(crit)
}

###################
# exp.crit.all function returns the three utility values (u(d,y)) for a given design=d and theta matrix=theta_int,
# where input B is the number of utility values cosidered to obtain the expected utilities.
exp.crit.all <- function(d,B)
{
  distM <- rdist(d)
  X1 <- ((d[,2]+.1)*5)-((d[,1]+.2)*5)
  X2 <- (d[,2]+.1)*(d[,1]+.2)
  X <- data.frame(X1,X2)
  
  crit <- foreach(j=1:B,.packages=c("mvtnorm","fields","corpcor","nlme","Matrix","matrixcalc","MaternEx1"),
				.errorhandling='stop',.export=ls(globalenv()),.combine=c) %dopar%  
    {
      Y <- Res_Clcpp(xdata=as.matrix(X),distM=distM,para=theta_int[j,],Z1=t(rnorm(1:nrow(X))),
				Z2=t(rnorm(1:nrow(X))),Z3=(rnorm(1:nrow(X))),v2=Unif_set[1:nrow(X)],r01=r01,r02=r02)
      crit.out <- combine_ute(d,X,Y,theta=theta_int[j,],distM)
      crit.dual <- crit.out[[3]]
      crit.kld <- crit.out[[4]]
      crit.pred <- crit.out[[5]]
      c(crit.kld,crit.pred,crit.dual) 
    }
  return(crit)
}

###################
# exp.crit.both.P function returns the dual purpose utility values, estimation utility values and the prediction utility 
# values based on the two approximations for a given design=d and theta matrix=theta_int, where input B is the number of 
# utility values cosidered to obtain the expected utilities.
exp.crit.both.P <- function(d,B)
{
  distM <- rdist(d)
  X1 <- ((d[,2]+.1)*5)-((d[,1]+.2)*5)
  X2 <- (d[,2]+.1)*(d[,1]+.2)
  X <- data.frame(X1,X2)
  
  crit <- foreach(j=1:B,.packages=c("mvtnorm","fields","corpcor","nlme","Matrix","matrixcalc","MaternEx1"),
			.errorhandling='stop',.export=ls(globalenv()),.combine=c) %dopar%  
    {
      Y <- Res_Clcpp(xdata=as.matrix(X),distM=distM,para=theta_int[j,],Z1=t(rnorm(1:nrow(X))),
			Z2=t(rnorm(1:nrow(X))),Z3=(rnorm(1:nrow(X))),v2=Unif_set[1:nrow(X)],r01=r01,r02=r02)
      crit.out <- combine_ute_both(d,X,Y,theta=theta_int[j,],distM)
      crit.dual <- crit.out[[3]]
      crit.kld <- crit.out[[4]]
      crit.pred.App <- crit.out[[5]]
      crit.pred <- crit.out[[6]]
      c(crit.dual,crit.kld,crit.pred.App,crit.pred) 
    }
  return(crit)
}

###################
# exp.crit.var function returns the posterior variances of the parameters, posterior prediction variances for response 1,
# and posterior prediction variances for response 2 for a given design=d and theta matrix=theta_int, where input B is  
# the number of datasets (Y sets) to evaluate the design.
exp.crit.var <- function(d,B)
{
  distM <- rdist(d)
  X1 <- ((d[,2]+.1)*5)-((d[,1]+.2)*5)
  X2 <- (d[,2]+.1)*(d[,1]+.2)
  X <- data.frame(X1,X2)
  
  Y1.post.pred <- matrix(ncol=nrow(Unsamp_X),nrow=B)
  Y2.post.pred <- matrix(ncol=nrow(Unsamp_X),nrow=B)
  post.var <- matrix(ncol=ncol(prior),nrow=B)
  
  for(j in 1:B)
  {
    Y <- Res_Clcpp(xdata=as.matrix(X),distM=distM,para=theta_int[j,],Z1=t(rnorm(1:nrow(X))),
			Z2=t(rnorm(1:nrow(X))),Z3=(rnorm(1:nrow(X))),v2=Unif_set[1:nrow(X)],r01=r01,r02=r02)
    crit.out <- Post_pred(d,X,Y,theta=theta_int[j,],distM)
    post.var[j,] <- diag(crit.out[[2]])
    Y1.post.pred[j,] <- apply(crit.out[[3]],1,var,na.rm=T)
    Y2.post.pred[j,] <- log(apply(crit.out[[4]],1,var,na.rm=T))
    
  }
  
  Y1.post <- apply(Y1.post.pred,2,mean,na.rm=T)
  Y2.post <- apply(Y2.post.pred,2,mean,na.rm=T)
  return(list(post.var,Y1.post,Y2.post))
}

###################
# Laplace_approx function returns the minimum negative log posterior value and its corresponding parameter vector and the
# Hessian matrix for given values for independent variables=X, data=Y, parameter vector=theta, and distance matrix=distM
Laplace_approx <- function(X,Y,theta,distM)
{
  #lp.approx <- optim(par=theta,X=X,Y=Y,distM=distM,fn=log_post,hessian=TRUE,method="L-BFGS-B",
  #					lower=c(apply(prior,2,min)),upper=c(apply(prior,2,max)))
  lp.approx <- nlminb(start=theta,X=X,distM=distM,Y=Y,objective=log_post,lower=c(apply(prior,2,min)),
				upper=c(apply(prior,2,max)),control=list(eval.max=10000,iter.max=10000,sing.tol=1e-18,rel.tol=1e-14)) 
  lp.approx
}

###################
# log_post function returns the negative log posterior value for given values for independent variables=X, data=Y,
# parameter vector=theta, and distance matrix=distM
log_post <- function(X,Y,theta,distM)
{
  log.prior <- dmvnorm(x=theta,mean=mu_prior,sigma=Sigma_prior,log=TRUE)
  log.like <- sum(Log_likelihoodCpp(xdata=as.matrix(X),y=as.matrix(Y),distM=distM,para=theta,
				Z=t(Z1[1:nrow(Y),(1:B2)]),r01,r02)) 
  
  Neg_log_post <- (-1)*(log.prior + log.like)
  Neg_log_post[is.na(Neg_log_post) | Neg_log_post==Inf | Neg_log_post==(-Inf)]=1e10
  return(Neg_log_post)
}

###################
# Post_pred function returns the posterior mean vector, posterior variances of the parameters, posterior prediction
# variances for response 1, and posterior prediction variances for response 2 for a given design=d, values for
# independent variables=X, data=Y, parameter vector=theta and distance matrix=distM
Post_pred <- function(d,X,Y,theta,distM){
  
  Sol <- Laplace_approx(X,Y,theta,distM)
  mu_post <- Sol$par
  hessian1 <- optimHess(par=mu_post,distM=distM,fn=log_post,X=X,Y=Y)
  
  if(Sol$convergence==1)
  {
    print("error1")
  }
  
  if(!is.symmetric.matrix(hessian1)){
    print("error2")
    hess1 <- hessian1
    hess1[lower.tri(hess1,diag=F)] <- 0
    hessian1 <- hess1+t(hess1)
    if(!is.positive.semi.definite(hessian1,tol=1e-16) | det(hessian1)==0) {
      Sigma_post <- Sigma_prior
      post_samp1 <- rmvnorm(100,mu_post,Sigma_post)
      Y1.mat=matrix(nrow=nrow(Unsamp_X),ncol=100)
      Y2.mat=matrix(nrow=nrow(Unsamp_X),ncol=100)
      for(i in 1:100){
        Y.post <- Res_Clcpp(xdata=as.matrix(Unsamp_X[,4:5]),distM=distM_uns,para=post_samp1[i,],
					Z1=t(rnorm(1:nrow(Unsamp_X))),Z2=t(rnorm(1:nrow(Unsamp_X))),Z3=(rnorm(1:nrow(Unsamp_X))),
					v2=Unif_set[1:nrow(Unsamp_X)],r01=r01,r02=r02)
        Y1.mat[,i] <- Y.post[,1]
        Y2.mat[,i] <- Y.post[,2]
      }
    }else{
      Sigma_post <- solve(hessian1)
      post_samp1 <- rmvnorm(100,mu_post,Sigma_post)
      Y1.mat <- matrix(nrow=nrow(Unsamp_X),ncol=100)
      Y2.mat <- matrix(nrow=nrow(Unsamp_X),ncol=100)
      for(i in 1:100){
        Y.post  <-  Res_Clcpp(xdata=as.matrix(Unsamp_X[,4:5]),distM=distM_uns,para=post_samp1[i,],
						Z1=t(rnorm(1:nrow(Unsamp_X))),Z2=t(rnorm(1:nrow(Unsamp_X))),Z3=(rnorm(1:nrow(Unsamp_X))),
						v2=Unif_set[1:nrow(Unsamp_X)],r01=r01,r02=r02)
        Y1.mat[,i] <- Y.post[,1]
        Y2.mat[,i] <- Y.post[,2]
      } 
    }
  }else if(!is.positive.semi.definite(hessian1,tol=1e-16) | det(hessian1)==0) {
    print("error3")
    Sigma_post <- Sigma_prior
    post_samp1 <- rmvnorm(100,mu_post,Sigma_post)
    Y1.mat <- matrix(nrow=nrow(Unsamp_X),ncol=100)
    Y2.mat <- matrix(nrow=nrow(Unsamp_X),ncol=100)
    for(i in 1:100){
      Y.post <- Res_Clcpp(xdata=as.matrix(Unsamp_X[,4:5]),distM=distM_uns,para=post_samp1[i,],
					Z1=t(rnorm(1:nrow(Unsamp_X))),Z2=t(rnorm(1:nrow(Unsamp_X))),Z3=(rnorm(1:nrow(Unsamp_X))),
					v2=Unif_set[1:nrow(Unsamp_X)],r01=r01,r02=r02)
      Y1.mat[,i] <- Y.post[,1]
      Y2.mat[,i] <- Y.post[,2]
    } 
  }else{
    Sigma_post <- solve(hessian1)
    post_samp1 <- rmvnorm(100,mu_post,Sigma_post)
    Y1.mat <- matrix(nrow=nrow(Unsamp_X),ncol=100)
    Y2.mat <- matrix(nrow=nrow(Unsamp_X),ncol=100)
    for(i in 1:100){
      Y.post <- Res_Clcpp(xdata=as.matrix(Unsamp_X[,4:5]),distM=distM_uns,para=post_samp1[i,],
					Z1=t(rnorm(1:nrow(Unsamp_X))),Z2=t(rnorm(1:nrow(Unsamp_X))),Z3=(rnorm(1:nrow(Unsamp_X))),
					v2=Unif_set[1:nrow(Unsamp_X)],r01=r01,r02=r02)
      Y1.mat[,i] <- Y.post[,1]
      Y2.mat[,i] <- Y.post[,2]
    }
  }
  
  return(list(mu_post,Sigma_post,Y1.mat,Y2.mat)) 
}

###################
# Pred_ute function returns the prediction utilty (using the first approximation) for a given design=d and data=Y,
# based on the inputs: values for the independent variables=X_unsamp,distance matrix=distM, 
# posterior sample 1=post_samp1, posterior sample 2=post_samp2, partial sill=r01, and range parameter=r02.
Pred_ute <- function(X_unsamp,distM,post_samp1,post_samp2,r01,r02)
{
  K <- nrow(post_samp1)
  R <- nrow(post_samp2)
  Pz <- rep(0,K)
  for(m in 1:K){
    Ym <- Res_Clcpp(xdata=as.matrix(X_unsamp[,4:5]),distM=distM,para=post_samp1[j,],Z1=t(rnorm(nrow(X_unsamp))),
			Z2=t(rnorm(nrow(X_unsamp))),Z3=rnorm(nrow(X_unsamp)),v2=Unif_set[1:nrow(X_unsamp)],r01=r01,r02=r02)
    log_z <- rep(0,R)
    for(r in 1:R){
      log_z[r] <- sum(Log_likelihoodCpp(as.matrix(X_unsamp[,4:5]),as.matrix(Ym),distM,post_samp2[r,],
					Z=t(Z1[1:nrow(Ym),(1:B2)]),r01,r02))
    }
    Pz[m] <- mean(exp(log_z),na.omit=TRUE) 
  }
  pred_ute <- mean(log(Pz),na.omit=TRUE)
  pred_ute[is.na(pred_ute) | pred_ute < (-10000)]=(-10000)
  
  return(pred_ute)
}

###################
# Pred_ute.App function returns the prediction utilty (using the second approximation) for a given design=d and data=Y,
# based on the input values: values for the independent variables=X_unsamp,distance matrix=distM, 
# posterior sample=post_samp1, partial sill=r01, and range parameter=r02.
Pred_ute.App <- function(X_unsamp,distM,post_samp1,r01,r02)
{
  
  K <- nrow(post_samp1)
  Avg_Hz <- matrix(nrow=K,ncol=nrow(X_unsamp))
  
  for(m in 1:K){
    Y1 <- matrix(ncol=nrow(X_unsamp),nrow=100)
    Y2 <- matrix(ncol=nrow(X_unsamp),nrow=100)
    for(j in 1:100){
      Ym <- Res_Clcpp(xdata=as.matrix(X_unsamp[,4:5]),distM=distM,para=post_samp1[j,],Z1=t(rnorm(nrow(X_unsamp))),
				Z2=t(rnorm(nrow(X_unsamp))),Z3=rnorm(nrow(X_unsamp)),v2=Unif_set[1:nrow(X_unsamp)],r01=r01,r02=r02)
      Y1[j,] <- Ym[,1]
      Y2[j,] <- Ym[,2]
    }
    Hz <- c()
    for(i in 1:ncol(Y1)){ 
      data.Y <- cbind(Y1[,i],Y2[,i])
      Sig.Y <- cov(data.Y)
      Hz[i] <- (ncol(data.Y)/2)*(1+log(2*pi))+log(det(Sig.Y))/2
    }
    Avg_Hz[m,] <- Hz
  }
  
  Avg_Hz[Avg_Hz==(-Inf) | Avg_Hz==Inf] <- 0
  pred_ute <- (-sum(colMeans(Avg_Hz,na.rm=T)))
  pred_ute[is.na(pred_ute) | pred_ute < (-1000)]=(-1000)
  pred_ute[pred_ute==Inf | pred_ute==(-Inf)]=NA
  
  return(pred_ute)
}

###################
# trace_mat function returns the trace value of an input matrix M.
trace_mat <- function(M)
{
  tr <- sum(diag(M))
  tr
}