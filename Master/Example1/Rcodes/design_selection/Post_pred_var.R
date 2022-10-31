Post_pred <- function(d,X,Y,theta,distM){
  
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
         post_samp1 <-  rmvnorm(100,mu_post,Sigma_post)
         Y1.mat=matrix(nrow=nrow(Unsamp_X),ncol=100)
         Y2.mat=matrix(nrow=nrow(Unsamp_X),ncol=100)
         for(i in 1:100){
           Y.post = Res_Clcpp(xdata=as.matrix(Unsamp_X[,4:5]),distM=distM_uns,para=post_samp1[i,],Z1=t(Z_res[1:nrow(Unsamp_X),1]),Z2=t(Z_res[1:nrow(Unsamp_X),2]),Z3=Z_res[1:nrow(Unsamp_X),3],v2=Unif_set[1:nrow(Unsamp_X)],r01=r01,r02=r02)
           Y1.mat[,i]=Y.post[,1]
           Y2.mat[,i]=Y.post[,2]
         }
        }else{
		    Sigma_post <- solve(hessian1)
		    post_samp1 <-  rmvnorm(100,mu_post,Sigma_post)
		    Y1.mat=matrix(nrow=nrow(Unsamp_X),ncol=100)
		    Y2.mat=matrix(nrow=nrow(Unsamp_X),ncol=100)
		    for(i in 1:100){
		      Y.post = Res_Clcpp(xdata=as.matrix(Unsamp_X[,4:5]),distM=distM_uns,para=post_samp1[i,],Z1=t(Z_res[1:nrow(Unsamp_X),1]),Z2=t(Z_res[1:nrow(Unsamp_X),2]),Z3=Z_res[1:nrow(Unsamp_X),3],v2=Unif_set[1:nrow(Unsamp_X)],r01=r01,r02=r02)
		      Y1.mat[,i]=Y.post[,1]
		      Y2.mat[,i]=Y.post[,2]
		    } 
		}
  }else if(!is.positive.semi.definite(hessian1, tol=1e-16) | det(hessian1)==0) {
    print("error3")
    Sigma_post <- Sigma_prior
    post_samp1 <-  rmvnorm(100,mu_post,Sigma_post)
    Y1.mat=matrix(nrow=nrow(Unsamp_X),ncol=100)
    Y2.mat=matrix(nrow=nrow(Unsamp_X),ncol=100)
    for(i in 1:100){
      Y.post = Res_Clcpp(xdata=as.matrix(Unsamp_X[,4:5]),distM=distM_uns,para=post_samp1[i,],Z1=t(Z_res[1:nrow(Unsamp_X),1]),Z2=t(Z_res[1:nrow(Unsamp_X),2]),Z3=Z_res[1:nrow(Unsamp_X),3],v2=Unif_set[1:nrow(Unsamp_X)],r01=r01,r02=r02)
      Y1.mat[,i]=Y.post[,1]
      Y2.mat[,i]=Y.post[,2]
    } 
  }else{
    Sigma_post <- solve(hessian1)
    post_samp1 <-  rmvnorm(100,mu_post,Sigma_post)
    Y1.mat=matrix(nrow=nrow(Unsamp_X),ncol=100)
    Y2.mat=matrix(nrow=nrow(Unsamp_X),ncol=100)
    for(i in 1:100){
      Y.post = Res_Clcpp(xdata=as.matrix(Unsamp_X[,4:5]),distM=distM_uns,para=post_samp1[i,],Z1=t(Z_res[1:nrow(Unsamp_X),1]),Z2=t(Z_res[1:nrow(Unsamp_X),2]),Z3=Z_res[1:nrow(Unsamp_X),3],v2=Unif_set[1:nrow(Unsamp_X)],r01=r01,r02=r02)
      Y1.mat[,i]=Y.post[,1]
      Y2.mat[,i]=Y.post[,2]
    }
  }
  
  return(list(mu_post,Sigma_post,Y1.mat,Y2.mat)) 
}