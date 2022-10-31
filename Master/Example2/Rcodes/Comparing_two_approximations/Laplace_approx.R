Laplace_approx<- function(X,Y,theta,distM)
{
   lp.approx <-nlminb(start=theta,X=X,distM=distM, Y=Y, objective =log_post,control = list(eval.max=5000,iter.max=5000,sing.tol=1e-18,rel.tol=1e-14),lower=c(apply(prior,2,min)),upper=c(apply(prior,2,max)))
   lp.approx
}