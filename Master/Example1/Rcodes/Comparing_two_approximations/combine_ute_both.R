combine_ute_both <- function(d,X,Y,theta,distM){
  
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
      if (!is.positive.semi.definite(hessian1, tol=1e-16) | det(hessian1)==0) {
        Sigma_post <- Sigma_prior
        kld.out <- NA
        pred.out <- NA
        pred.out2 <- NA
        det.out <- NA 
      }else{
        Sigma_post <- solve(hessian1)
        log_det_post=log(det(Sigma_post))
        kld.out <- 0.5*(trace_mat(iSigma_prior%*%Sigma_post) + t(mu_prior-mu_post)%*%iSigma_prior%*%(mu_prior-mu_post)-num_par + log_det_prior-log_det_post)
        post_samp <-  rmvnorm(100,mu_post,Sigma_post)
        t1=proc.time()
        pred.out <- Pred_ute.App(X_unsamp=Unsamp_X,distM=distM_uns,post_samp1=post_samp,r01=r01,r02=r02)
		    t2=proc.time()
		    time1= t2-t1
		    
		    post_samp2 <-  rmvnorm(100,mu_post,Sigma_post)
		    t3=proc.time()
		    pred.out2 <- Pred_ute(X_unsamp=Unsamp_X,distM=distM_uns,post_samp1=post_samp,post_samp2=post_samp2,r01=r01,r02=r02)
		    t4=proc.time()
		    time2= t4-t3
		    
		    det.out <- kld.out+pred.out  
		}
      
  }else if (!is.positive.semi.definite(hessian1, tol=1e-16) | det(hessian1)==0) {
    print("error3")
    Sigma_post <- Sigma_prior
    kld.out <- NA
	  pred.out <- NA
	  pred.out2 <- NA
	  det.out <- NA 
  }else{
    Sigma_post <- solve(hessian1)
    log_det_post=log(det(Sigma_post))
    kld.out <- 0.5*(trace_mat(iSigma_prior%*%Sigma_post) + t(mu_prior-mu_post)%*%iSigma_prior%*%(mu_prior-mu_post)-num_par + log_det_prior-log_det_post)
    set.seed(NAI*101)
    post_samp <-  rmvnorm(100,mu_post,Sigma_post)
    t1=proc.time()
    pred.out <- Pred_ute.App(X_unsamp=Unsamp_X,distM=distM_uns,post_samp1=post_samp,r01=r01,r02=r02)
    t2=proc.time()
    time1= t2-t1
    
    post_samp2 <-  rmvnorm(100,mu_post,Sigma_post)
    t3=proc.time()
    pred.out2 <- Pred_ute(X_unsamp=Unsamp_X,distM=distM_uns,post_samp1=post_samp,post_samp2=post_samp2,r01=r01,r02=r02)
    t4=proc.time()
    time2= t4-t3
    det.out <- kld.out+pred.out
  }
  
  return(list(mu_post,Sigma_post,det.out,kld.out,pred.out,pred.out2,time1,time2)) 
}