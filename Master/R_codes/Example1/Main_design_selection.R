
# Parameters that can be changed based on the requirement
d_no <- 5   # number of design points
utility <- 1 # 1=dual_purpose_utility,2=Estimation_utility,3=prediction_utility 
Dependence <- 0.2 #0.2=week,0.5=moderate,0.8=strong
Des <- c("dual","est","pred")

#Required R-packages
library(MaternEx1) #Newly created package for the example
library(mvtnorm)
library(fields)
library(Matrix)
library(matrixcalc)
library(nlme)
library(corpcor)
library(acebayes)
library(doParallel)

source("R_codes/Example1/functions_Ex1.R")

#AI <- Sys.getenv("casenumber") 
#AI is a parameter that can be pass externally when running this code.
AI <- 1
NAI <- as.numeric(AI)
out_name <- paste("ClOut",NAI,".txt",sep="")

#Create a cluster in R
cl <- makeCluster(30,outfile=out_name)
registerDoParallel(cl)

#Change the library path
#.libPaths(c("/home/Rpack",.libPaths()))
#clusterCall(cl,function(x) .libPaths(x),.libPaths())

#Instead of clusterCall() function given above, can also use
#clusterEvalQ(cl, .libPaths("/home/Rpack"))

# Prior distribution
mu_prior1 <- c(5,-2.8,8,log(1.2),log(.7/Dependence),log(Dependence),log(1.5))  
Sigma_prior1 <- c(4,4,4,0.25,0.25,0.25,0.25)
mu_prior2 <- c(3.8,-0.5,-0.7,log(.6/Dependence),log(Dependence),log(0.25))   
Sigma_prior2 <- c(0.125,0.125,0.125,0.125,0.125,0.25)
mu_Ltau <- log(.7/.3)
Sigma_Ltau <- 0.25

mu_prior <- c(mu_prior1,mu_prior2,mu_Ltau)
Sigma_prior <- diag(c(Sigma_prior1,Sigma_prior2,Sigma_Ltau))
iSigma_prior <- solve(Sigma_prior)
log_det_prior <- log(det(Sigma_prior))

set.seed(NAI*101) #set a seed value for reproduce the results
N <- 5000 # number of samples from the prior
prior <- rmvnorm(N,mu_prior,Sigma_prior) # Initial particle set

num_par <- length(mu_prior)
r01 <- 0.001
r02 <- 0.001

#Prediction locations
Unsamp_Lx <- rep(seq(0,1,.25),5)
Unsamp_Ly <- rep(seq(0,1,.25),each=5)

#Values of independent variables at prediction locations
X1 <- ((Unsamp_Ly+.1)*5)-((Unsamp_Lx+.2)*5)
X2 <- (Unsamp_Ly+.1)*(Unsamp_Lx+.2)

site_no <- 1:25
Unsamp_X <- data.frame(site_no,Unsamp_Lx,Unsamp_Ly,X1,X2)
distM_uns <- rdist(data.frame(Unsamp_Lx,Unsamp_Ly)) # distance matrix

####### Coordinate exchange algorithm ##########
# Parameters and particle sets for the simulations
set.seed(101)
tot_sites <- 1000
B1 <- 2000
theta_int=rmvnorm(n=B1,mean=mu_prior,sigma=Sigma_prior)
B2 <- 1000
Z1 <- matrix(rnorm(tot_sites*B2),tot_sites)
Unif_set <- runif(tot_sites,0,1)
Z_res <- matrix(rnorm(tot_sites*3),ncol=3)

# Initial design
set.seed(100*NAI)
dint <- matrix(runif(d_no*2),ncol=2)

for(i in 1:5){
opt_data <- ace(exp.crit,start.d=dint,B=c(1000,300),Q=10,N1=1,N2=0,lower=0,upper=1)
opt_d <- opt_data$phase1.d
Trace.val <- opt_data$phase1.trace

out.data <- list(opt_d,Trace.val,opt_data) 
dint <- opt_data$phase1.d
}

f.name <- paste(Des[utility],"_R",Dependence*10,"_D",d_no,sep="")
v_name <- paste("Results/Example1/Selected_designs/",f.name,".RData",sep="")
save(out.data,file=v_name)

#######################################

