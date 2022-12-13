
# Parameters that can be changed based on the requirement
d_no <- 10   # number of design points
Dependence <- 0.2 #0.2=week, 0.5=moderate, 0.8=strong

#Required R-packages
library(MaternEx1) #Newly created package for the example
library(mvtnorm)
library(fields)
library(Matrix)
library(matrixcalc)
library(nlme)
library(corpcor)
library(doParallel)

source("R_codes/Example1/functions_Ex1.R")

#AI is a parameter(integer) that can be pass externally when running this code.
#AI <- Sys.getenv("casenumber") 
AI=1
NAI <- as.numeric(AI)

# Prior distribution
mu_prior1 <- c(5,-2.8,8,log(1.2),log(.7/Dependence),log(Dependence),log(1.5))  
Sigma_prior1 <- c(4,4,4,0.25,0.25,0.25,0.25)
mu_prior2 <- c(3.8,-0.5,-0.7,log(.6/Dependence),log(Dependence),log(0.25))   
Sigma_prior2 <- c(0.125,0.125,0.125,0.125,0.125,0.25)
mu_Ltau <- log(.7/.3)
Sigma_Ltau <- 0.25

mu_prior<- c(mu_prior1,mu_prior2,mu_Ltau)
Sigma_prior <- diag(c(Sigma_prior1,Sigma_prior2,Sigma_Ltau))
iSigma_prior <- solve(Sigma_prior)
log_det_prior=log(det(Sigma_prior))

set.seed(NAI*101) #set a seed value for reproduce the results
N <- 5000 # number of samples from the prior
prior <- rmvnorm(N,mu_prior,Sigma_prior) # Initial particle set

num_par<-length(mu_prior)
r01=0.001
r02=0.001

#Prediction locations
Unsamp_Lx<-rep(seq(0,1,.25),5)
Unsamp_Ly<-rep(seq(0,1,.25),each=5)

#Values of independent variables at prediction locations
X1 <- ((Unsamp_Ly+.1)*5)-((Unsamp_Lx+.2)*5)
X2 <- (Unsamp_Ly+.1)*(Unsamp_Lx+.2)

site_no=1:25
Unsamp_X=data.frame(site_no,Unsamp_Lx,Unsamp_Ly,X1,X2)
distM_uns=rdist(data.frame(Unsamp_Lx,Unsamp_Ly)) # distance matrix

# Parameters and particle sets for the simulations
set.seed(101)
tot_sites <- 1000
B1 <- 500
theta_int=rmvnorm(n=B1,mean=mu_prior, sigma = Sigma_prior)
B2 <- 1000
Z1 <- matrix(rnorm(tot_sites*B2), tot_sites)
Unif_set = runif(tot_sites,0,1)
Z_res <- matrix(rnorm(tot_sites*3),ncol=3)

# Equally spaced triangular design
x1=c(0,1/6,1/3,1/3,1/2,1/2,2/3,2/3,5/6,1)
y1=c(0,1/3,0,2/3,1/3,1,0,2/3,1/3,0)
EqSp.D=data.frame(X1=x1,X2=y1)

# A design with all the points on the boundary
x2=c(0,.5,1,1,0,1,0,0,.5,1)
y2=c(0,0,0.25,0,.25,.75,.75,1,1,1)
Bnd.D=data.frame(X1=x2,X2=y2)

# A design with all points close to the prediction locations 
x3=c(0,.5,.75,1,0,1,0.25,0,.5,1)
y3=c(0,0,0.5,0,.5,.5,.5,1,1,1)
PrC.D=data.frame(X1=x3,X2=y3)

# Estimation design
d.est<- paste("Results/Example1/Selected_designs/est_R",Dependence*10,"_D",d_no,".RData",sep="")
load(d.est)
Est.D <- out.data[[1]]
  
# Dual-purpose design
d.dual<- paste("Results/Example1/Selected_designs/dual_R",Dependence*10,"_D",d_no,".RData",sep="")
load(d.dual)
Dual.D <- out.data[[1]]

# Prediction design
d.pred<- paste("Results/Example1/Selected_designs/pred_R",Dependence*10,"_D",d_no,".RData",sep="")
load(d.pred)
Pred.D <- out.data[[1]]

D.list=list(EqSp.D,Bnd.D,PrC.D,Est.D,Dual.D,Pred.D)

Y1.pred=matrix(ncol=nrow(Unsamp_X),nrow=B1)
Y2.pred=matrix(ncol=nrow(Unsamp_X),nrow=B1)
for(k in 1:B1){
  prior_pred <- Res_Clcpp(xdata=as.matrix(Unsamp_X[,4:5]),distM=distM_uns,
                    para=theta_int[k,],Z1=t(Z_res[1:nrow(Unsamp_X),1]),
                    Z2=t(Z_res[1:nrow(Unsamp_X),2]),Z3=Z_res[1:nrow(Unsamp_X),3],
                    v2=Unif_set[1:nrow(Unsamp_X)],r01=r01,r02=r02)
  Y1.pred[k,] <- prior_pred[,1]
  Y2.pred[k,] <- prior_pred[,2]
}
Y1.var <- apply(Y1.pred,2,var,na.rm=T)
Y2.var <- log(apply(Y2.pred,2,var,na.rm=T))

data.new.Y1 <- data.frame(site_no,Y1.var)
data.new.Y2 <- data.frame(site_no,Y2.var)

for(i in 1:6){
  d.opt = D.list[[i]]
  Crit.out <- exp.crit.var(d=d.opt,B=B1)
  data.new.Y1[,i+2] <- Crit.out[[2]]
  data.new.Y2[,i+2] <- Crit.out[[3]]
  out[[i]] <- Crit.out[[1]]
}

colnames(data.new.Y1)[3:8] <-c("Y1.P.P1","Y1.E.P1","Y1.D.P1","Y1.Bnd.P1",
                                "Y1.Eqsp.P1","Y1.PrC.P1")
colnames(data.new.Y2)[3:8] <-c("log.Y2.P.P1","log.Y2.E.P1","log.Y2.D.P1",
                                "log.Y2.Bnd.P1","log.Y2.Eqsp.P1","log.Y2.PrC.P1")

save(data.new.Y1,
     file="Results/Suplementary_results/Compare_common_designs/Post_pred_var_Y1.RData")
save(data.new.Y2,
     file="Results/Suplementary_results/Compare_common_designs/Post_pred_var_Y2.RData")
save(out,
     file="Results/Suplementary_results/Compare_common_designs/Post_var.RData")
