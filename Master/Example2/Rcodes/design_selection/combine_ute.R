combine_ute <- function(d,X,Y,theta,distM){
  
  Sol <- Laplace_approx(X,Y,theta,distM)
  mu_post <- Sol$par
  hessian1 <- optimHess(par=mu_post,fn=log_post,X=X,Y=Y,distM=distM)
   
  if(!is.symmetric.matrix(hessian1)){
      hess1 <- hessian1
      hess1[lower.tri(hess1, diag = F)] <- 0
      hessian1 <- hess1+t(hess1)
      if (!is.positive.semi.definite(hessian1, tol=1e-16) | det(hessian1)==0) {
        Sigma_post <- Sigma_prior
        kld.out <- NA
        pred.out <- NA
        det.out <- NA 
      }else{
        Sigma_post <- solve(hessian1)
        log_det_post <- log(det(Sigma_post))
        kld.out <- 0.5*(trace_mat(iSigma_prior%*%Sigma_post) + t(mu_prior-mu_post)%*%iSigma_prior%*%(mu_prior-mu_post)-num_par + log_det_prior-log_det_post)
		    post_samp <-  rmvnorm(500,mu_post,Sigma_post)
		    pred.out <- Pred_utecpp(X_unsamp=X.all,distM=dist_mat,post1=post_samp,v2=Unif_set[1:nrow(X.all)],r01=r01,r02=r02)
		    det.out <- kld.out+pred.out  
		}
      
  }else if (!is.positive.semi.definite(hessian1, tol=1e-16) | det(hessian1)==0) {
    Sigma_post <- Sigma_prior
    kld.out <- NA
	  pred.out <- NA
	  det.out <- NA 
  }else{
    Sigma_post <- solve(hessian1)
    log_det_post=log(det(Sigma_post))
    kld.out <- 0.5*(trace_mat(iSigma_prior%*%Sigma_post) + t(mu_prior-mu_post)%*%iSigma_prior%*%(mu_prior-mu_post)-num_par + log_det_prior-log_det_post)
    post_samp <-  rmvnorm(500,mu_post,Sigma_post)
    pred.out <- Pred_utecpp(X_unsamp=X.all,distM=dist_mat,post1=post_samp,v2=Unif_set[1:nrow(X.all)],r01=r01,r02=r02)
    det.out <- kld.out+pred.out
  }
  
  return(list(mu_post,Sigma_post,det.out,kld.out,pred.out)) 
}