Pred_ute.App <- function(X_unsamp,distM,post_samp1,v2,r01,r02)
  {
  
  K <- nrow(post_samp1)
  Avg_Hz <- matrix(nrow=K,ncol=nrow(X_unsamp))

  for(m in 1:K){
    Y1 <- matrix(ncol=nrow(X_unsamp),nrow=100)
    Y2 <- matrix(ncol=nrow(X_unsamp),nrow=100)
    for(j in 1:100){
      Ym <- response_Cl(xdata=X_unsamp,distM=distM,para=post_samp1[m,],Z1=cbind(rnorm(nrow(X_unsamp)),rnorm(nrow(X_unsamp))),Z2=c(rnorm(nrow(X_unsamp))),V2=v2)
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
  
  Avg_Hz[Avg_Hz== -Inf | Avg_Hz== Inf] <- 0
  pred_ute <- (-sum(colMeans(Avg_Hz,na.rm=T)))
  
  pred_ute[is.na(pred_ute) | pred_ute < (-1000)] <- (-1000)
  pred_ute[pred_ute== Inf | pred_ute == -Inf] <- NA
  
  return(pred_ute)
}