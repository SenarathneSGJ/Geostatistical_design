# Combine_ute function returns the posterior mean, posterior covariance matrix and the three utility values
# for a given design=d, values for independent variables=X, data=Y, parameter vector=theta, and distance matrix=distM
combine_ute <- function(d,X,Y,theta,distM){
  
  Sol <- Laplace_approx(X,Y,theta,distM)
  mu_post <- Sol$par
  hessian1 <- optimHess(par=mu_post,fn=log_post,X=X,Y=Y,distM=distM)
  
  if(!is.symmetric.matrix(hessian1)){
    hess1 <- hessian1
    hess1[lower.tri(hess1,diag=F)] <- 0
    hessian1 <- hess1+t(hess1)
    if (!is.positive.semi.definite(hessian1,tol=1e-16) | det(hessian1)==0) {
      Sigma_post <- Sigma_prior
      kld.out <- NA
      pred.out <- NA
      det.out <- NA 
    }else{
      Sigma_post <- solve(hessian1)
      log_det_post <- log(det(Sigma_post))
      kld.out <- 0.5*(trace_mat(iSigma_prior%*%Sigma_post) + 
					        t(mu_prior-mu_post)%*%iSigma_prior%*%(mu_prior-mu_post)-
					        num_par + log_det_prior-log_det_post)
      post_samp <-  rmvnorm(500,mu_post,Sigma_post)
      pred.out <- Pred_utecpp(X_unsamp=X.all,distM=dist_mat,post1=post_samp,
                      v2=Unif_set[1:nrow(X.all)],r01=r01,r02=r02)
      det.out <- kld.out+pred.out  
    }
    
  }else if (!is.positive.semi.definite(hessian1,tol=1e-16) | det(hessian1)==0) {
    Sigma_post <- Sigma_prior
    kld.out <- NA
    pred.out <- NA
    det.out <- NA 
  }else{
    Sigma_post <- solve(hessian1)
    log_det_post=log(det(Sigma_post))
    kld.out <- 0.5*(trace_mat(iSigma_prior%*%Sigma_post) + 
                  t(mu_prior-mu_post)%*%iSigma_prior%*%(mu_prior-mu_post)-
				          num_par + log_det_prior-log_det_post)
    post_samp <-  rmvnorm(500,mu_post,Sigma_post)
    pred.out <- Pred_utecpp(X_unsamp=X.all,distM=dist_mat,post1=post_samp,
                            v2=Unif_set[1:nrow(X.all)],r01=r01,r02=r02)
    det.out <- kld.out+pred.out
  }
  
  return(list(mu_post,Sigma_post,det.out,kld.out,pred.out)) 
}

################
# combine_ute_all function returns the posterior mean, posterior covariance matrix, dual purpose utility value, estimation 
# utility value and the prediction utility values based on the two approximations and time taken to run each approxiation
# for a given design=d, values for independent variables=X, data=Y, parameter vector=theta, and the distance matrix=distM
combine_ute_all <- function(d,X,Y,theta,distM){
  
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
      det.out <- NA 
    }else{
      Sigma_post <- solve(hessian1)
      log_det_post <- log(det(Sigma_post))
      kld.out <- 0.5*(trace_mat(iSigma_prior%*%Sigma_post) + 
					        t(mu_prior-mu_post)%*%iSigma_prior%*%(mu_prior-mu_post)-
					        num_par + log_det_prior-log_det_post)
      post_samp <-  rmvnorm(100,mu_post,Sigma_post)
      t1 <- proc.time()
      pred.out <- Pred_ute.App(X_unsamp=X.all,distM=dist_mat,
                          post_samp1=post_samp,r01=r01,r02=r02)
      t2 <- proc.time()
      time1 <- t2-t1
      
      post_samp2 <-  rmvnorm(100,mu_post,Sigma_post)
      t3 <- proc.time()
      pred.out2 <- Pred_ute(X_unsamp=X.all,distM=dist_mat,
                        post_samp1=post_samp,post_samp2=post_samp2,r01=r01,r02=r02)
      t4 <- proc.time()
      time2 <- t4-t3
      
      det.out <- kld.out+pred.out  
    }
    
  }else if(!is.positive.semi.definite(hessian1,tol=1e-16) | det(hessian1)==0) {
    print("error3")
    Sigma_post <- Sigma_prior
    kld.out <- NA
    pred.out <- NA
    det.out <- NA 
  }else{
    Sigma_post <- solve(hessian1)
    log_det_post <- log(det(Sigma_post))
    kld.out <- 0.5*(trace_mat(iSigma_prior%*%Sigma_post) + 
                t(mu_prior-mu_post)%*%iSigma_prior%*%(mu_prior-mu_post)-
					      num_par + log_det_prior-log_det_post)
    post_samp <-  rmvnorm(100,mu_post,Sigma_post)
    t1 <- proc.time()
    pred.out <- Pred_ute.App(X_unsamp=X.all,distM=dist_mat,
                        post_samp1=post_samp,r01=r01,r02=r02)
    t2 <- proc.time()
    time1 <- t2-t1
    
    post_samp2 <-  rmvnorm(100,mu_post,Sigma_post)
    t3 <- proc.time()
    pred.out2 <- Pred_ute(X_unsamp=X.all,distM=dist_mat,
                    post_samp1=post_samp,post_samp2=post_samp2,r01=r01,r02=r02)
    t4 <- proc.time()
    time2 <- t4-t3
    det.out <- kld.out+pred.out
  }
  
  return(list(mu_post,Sigma_post,det.out,kld.out,pred.out,pred.out2,time1,time2)) 
}

################
# exp.crit function returns the expected utility values (u(d)) based on the selected utility function for a given design=d  
# and theta matrix=theta_int, where input B is the number of utility values cosidered to obtain the expected utilities.
exp.crit <- function(d,B)
{
  crit <- rep(0,B)
  crit_kld <- rep(0,B)
  crit_pred <- rep(0,B)
  ord <- order(d)
  locs1 <- which(Monitoring_sites$Site_No %in% d)
  locs <- locs1[ord]
  
  distM <- dist_mat[locs,locs]
  X_data <- X.all[locs,]
  
  for(j in 1:B)
  {
    Y <- Res_Clcpp(xdata=X_data,distM=distM,para=theta_int[j,],
          Z1=t(rnorm(1:nrow(X_data))),Z2=t(rnorm(1:nrow(X_data))),
          Z3=(rnorm(1:nrow(X_data))),v2=Unif_set[1:nrow(X_data)],r01=r01,r02=r02)
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

################
# exp.crit.all function returns the three utility values (u(d,y)) for a given design=d and theta matrix=theta_int,
# where input B is the number of utility values considered to obtain the expected utilities.
exp.crit.all <- function(d,B)
{
  crit <- list()
  ord <- order(d)
  locs1 <- which(Monitoring_sites$Site_No %in% d)
  locs <- locs1[ord]
  
  distM <- dist_mat[locs,locs]
  X_data <- X.all[locs,]
  
  for (j in 1:B)
  {
    Y <- Res_Clcpp(xdata=X_data,distM=distM,para=matrix(theta_int[j,]),
          Z1=t(rnorm(1:nrow(X_data))),Z2=t(rnorm(1:nrow(X_data))),
          Z3=(rnorm(1:nrow(X_data))),v2=Unif_set[1:nrow(X_data)],r01=r01,r02=r02)
    crit.out <- combine_ute(d=d,X=X_data,Y=Y,theta=theta_int[j,],distM=distM)
    crit[(j-1)*3+1] <- crit.out[[4]] 
    crit[(j-1)*3+2] <- crit.out[[5]] 
    crit[j*3] <- crit.out[[3]]
  } 
  return(crit)
}

################
# exp.crit.both.P function returns the dual purpose utility values, estimation utility values and the prediction utility 
# values based on the two approximations for a given design=d and theta matrix=theta_int, where input B is the number of 
# utility values cosidered to obtain the expected utilities.
exp.crit.both.P <- function(d,B)
{
  crit <- c()
  ord <- order(d)
  locs1 <- which(Monitoring_sites$Site_No %in% d)
  locs <- locs1[ord]
  
  distM <- dist_mat[locs,locs]
  Loc_sel <- ST_uns[locs,]
  X_data <- X.all[locs,]
  
  for(j in 1:B){
    Y <- Res_Clcpp(xdata=X_data,distM=distM,para=theta_int[j,],
          Z1=t(rnorm(1:nrow(X_data))),Z2=t(rnorm(1:nrow(X_data))),
          Z3=(rnorm(1:nrow(X_data))),v2=Unif_set[1:nrow(X_data)],r01=r01,r02=r02)
    crit.out <- combine_ute_all(d=d,X=X_data,distM=distM,Y=Y,theta=theta_int[j,])
    crit.dual <- crit.out[[3]]
    crit.kld <- crit.out[[4]] 
    crit.pred.App <- crit.out[[5]] 
    crit.pred <- crit.out[[6]]
    out.all <- c(crit.dual,crit.kld,crit.pred.App,crit.pred)
    crit <- c(crit,out.all)
  }
  return(crit)
}

################
# Laplace_approx function returns the minimum negative log posterior value and its corresponding parameter vector and the
# Hessian matrix for given values for independent variables=X, data=Y, parameter vector=theta, and distance matrix=distM
Laplace_approx <- function(X,Y,theta,distM)
{
  #lp.approx <- optim(par=theta,X=X,Y=Y,fn=log_post,hessian=TRUE)
  lp.approx <- nlminb(start=theta,X=X,distM=distM,Y=Y,objective=log_post,
                control=list(eval.max=5000,iter.max=5000,sing.tol=1e-18,rel.tol=1e-14),
                lower=c(apply(prior,2,min)),upper=c(apply(prior,2,max)))
  
  lp.approx
}

################
# log_post function returns the negative log posterior value for given values for independent variables=X, 
# parameter vector=theta, data=Y, and distance matrix=distM
log_post <- function(X,Y,theta,distM)
{
  log.prior <- dmvnorm(x=theta,mean=mu_prior,sigma=Sigma_prior,log=TRUE)
  log.like <- sum(Log_likelihoodCpp(xdata=X,y=as.matrix(Y),distM=distM,
                    para=theta,Z=t(Z1[1:nrow(X),]),r01=r01,r02=r02))
  
  Neg_log_post <- (-1)*(log.prior + log.like)
  
  return(Neg_log_post)
}

################
# Pred_ute function returns the prediction utilty (using the first approximation) for a given design=d and data=Y,
# based on the inputs: values for the independent variables=X_unsamp,distance matrix=distM, 
# posterior sample 1=post_samp1, posterior sample 2=post_samp2, partial sill=r01, and range parameter=r02.
Pred_ute <- function(X_unsamp,distM,post_samp1,post_samp2,r01,r02)
{
  
  K <- nrow(post_samp1)
  R <- nrow(post_samp2)
  
  Pz <- rep(0,K)
  for(m in 1:K){
    Ym <- Res_Clcpp(xdata=X_unsamp,distM=distM,para=post_samp1[m,],
            Z1=t(rnorm(nrow(X_unsamp))),Z2=t(rnorm(nrow(X_unsamp))),
            Z3=rnorm(nrow(X_unsamp)),v2=Unif_set[1:nrow(X_unsamp)],r01=r01,r02=r02) 
    log_z <- rep(0,R)
    for(r in 1:R){
      log_z[r] <- sum(Log_likelihoodCpp(xdata=X_unsamp,y=as.matrix(Ym),distM=distM,
                      para=post_samp2[r,],Z=t(Z1[1:nrow(X_unsamp),]),r01,r02))
    }
    Pz[m] <- mean(exp(log_z),na.omit=TRUE) 
  }
  pred_ute <- mean(log(Pz),na.omit=TRUE)
  pred_ute[is.na(pred_ute) | pred_ute < (-10000)] <- (-10000)
  
  return(pred_ute)
}

################
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
      Ym <- Res_Clcpp(xdata=X_unsamp,distM=distM,para=post_samp1[m,],
              Z1=t(rnorm(nrow(X_unsamp))),Z2=t(rnorm(nrow(X_unsamp))),
              Z3=rnorm(nrow(X_unsamp)),v2=Unif_set[1:nrow(X_unsamp)],r01=r01,r02=r02) 
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
  
  pred_ute[is.na(pred_ute) | pred_ute < (-1000)] <- (-1000)
  pred_ute[pred_ute==Inf | pred_ute==(-Inf)] <- NA
  
  return(pred_ute)
}

################
trace_mat <- function(M)
{
  tr <- sum(diag(M))
  tr
}