combine_ute <- function(d,X,Y,theta,distM){
  
  Sol <- Laplace_approx(X,Y,theta,distM)
  mu_post <- Sol$par
  hessian1=optimHess(par=mu_post,distM=distM,fn=log_post,X=X,Y=Y)
   
  if(Sol$convergence==1)
    {
  	  print("error1")
    }
  
  if(!is.symmetric.matrix(hessian1)){
      print("error2")
      hess1 <- hessian1
      hess1[lower.tri(hess1, diag = F)] <- 0
      hessian1 <- hess1+t(hess1)
      if(!is.positive.semi.definite(hessian1, tol=1e-16) | det(hessian1)==0) {
         Sigma_post <- Sigma_prior
         log_det_post=log(det(Sigma_post))
         kld.out <- 0.5*(trace_mat(iSigma_prior%*%Sigma_post) + t(mu_prior-mu_post)%*%iSigma_prior%*%(mu_prior-mu_post)-num_par + log_det_prior-log_det_post)
         set.seed(NAI*101)
         post_samp1 <-  rmvnorm(200,mu_post,Sigma_post)
         pred.out <- Pred_utecpp(X_unsamp=as.matrix(Unsamp_X[,4:5]),post1=as.matrix(post_samp1),distM=as.matrix(distM_uns),v2=Unif_set[1:nrow(Unsamp_X)],r01,r02)
         #pred.out <- Pred_ute.App(X_unsamp=Unsamp_X,distM=distM_uns,post_samp1=post_samp1,r01=r01,r02=r02)
         det.out <- kld.out+pred.out  
        }else{
		 Sigma_post <- solve(hessian1)
         log_det_post=log(det(Sigma_post))
         kld.out <- 0.5*(trace_mat(iSigma_prior%*%Sigma_post) + t(mu_prior-mu_post)%*%iSigma_prior%*%(mu_prior-mu_post)-num_par + log_det_prior-log_det_post)
         set.seed(NAI*101)
         post_samp1 <-  rmvnorm(200,mu_post,Sigma_post)
         pred.out <- Pred_utecpp(X_unsamp=as.matrix(Unsamp_X[,4:5]),post1=as.matrix(post_samp1),distM=as.matrix(distM_uns),v2=Unif_set[1:nrow(Unsamp_X)],r01,r02)
         #pred.out <- Pred_ute.App(X_unsamp=Unsamp_X,distM=distM_uns,post_samp1=post_samp1,r01=r01,r02=r02)
		 det.out <- kld.out+pred.out  
		}
      
  }else if(!is.positive.semi.definite(hessian1, tol=1e-16) | det(hessian1)==0) {
    print("error3")
    Sigma_post <- Sigma_prior
    log_det_post=log(det(Sigma_post))
    kld.out <- 0.5*(trace_mat(iSigma_prior%*%Sigma_post) + t(mu_prior-mu_post)%*%iSigma_prior%*%(mu_prior-mu_post)-num_par + log_det_prior-log_det_post)
    set.seed(NAI*101)
    post_samp1 <-  rmvnorm(200,mu_post,Sigma_post)
    pred.out <- Pred_utecpp(X_unsamp=as.matrix(Unsamp_X[,4:5]),post1=as.matrix(post_samp1),distM=as.matrix(distM_uns),v2=Unif_set[1:nrow(Unsamp_X)],r01,r02)
    #pred.out <- Pred_ute.App(X_unsamp=Unsamp_X,distM=distM_uns,post_samp1=post_samp1,r01=r01,r02=r02)
    det.out <- kld.out+pred.out   
  }else{
    Sigma_post <- solve(hessian1)
    log_det_post=log(det(Sigma_post))
    kld.out <- 0.5*(trace_mat(iSigma_prior%*%Sigma_post) + t(mu_prior-mu_post)%*%iSigma_prior%*%(mu_prior-mu_post)-num_par + log_det_prior-log_det_post)
    set.seed(NAI*101)
    post_samp1 <-  (rmvnorm(200,mu_post,Sigma_post))
    pred.out <- Pred_utecpp(X_unsamp=as.matrix(Unsamp_X[,4:5]),post1=as.matrix(post_samp1),distM=as.matrix(distM_uns),v2=Unif_set[1:nrow(Unsamp_X)],r01,r02)
    #pred.out <- Pred_ute.App(X_unsamp=Unsamp_X,distM=distM_uns,post_samp1=post_samp1,r01=r01,r02=r02)
    det.out <- kld.out+pred.out
  }
  
  return(list(mu_post,Sigma_post,det.out,kld.out,pred.out)) 
}