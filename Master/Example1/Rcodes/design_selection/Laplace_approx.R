Laplace_approx<- function(X,Y,theta,distM)
  {
   #lp.approx <- optim(par=theta,X=X, Y=Y,distM=distM, fn=log_post, hessian=TRUE,method="L-BFGS-B",lower=c(apply(prior,2,min)),upper=c(apply(prior,2,max)))
   lp.approx <-nlminb(start=theta,X=X,distM=distM, Y=Y, objective =log_post,lower=c(apply(prior,2,min)),upper=c(apply(prior,2,max)),control = list(eval.max=10000,iter.max=10000,sing.tol=1e-18,rel.tol=1e-14)) #control = list(eval.max=10000,iter.max=10000,sing.tol=1e-18,rel.tol=1e-14),
  
  lp.approx
}