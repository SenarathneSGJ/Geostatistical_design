
# Parameters that can be changed based on the requirement
d_no <- 5   # number of design points

#Required R-packages
#Required R-packages
library(AirQualityRcpp) #Newly created package for the example
library(mvtnorm)
library(fields)
library(Matrix)
library(matrixcalc)
library(nlme)
library(corpcor)

source("Example2/Rcodes/design_selection/combine_ute.R")
source("Example2/Rcodes/design_selection/Exp_crit_all.R")
source("Example2/Rcodes/design_selection/Laplace_approx.R")
source("Example2/Rcodes/design_selection/log_posterior.R")
source("Example2/Rcodes/design_selection/trace_mat.R")

AI <- Sys.getenv("casenumber") #AI is a parameter(integer) that can be pass externally when running this code.(Eg: AI=1)
NAI <- as.numeric(AI)

load("Example2/Rcodes/design_selection/LP_post_New.RData")   # prior distribution
mu_prior <- LP_post[[1]]
Sigma_prior <- LP_post[[2]]   
iSigma_prior <- solve(Sigma_prior)
log_det_prior <- log(det(Sigma_prior))

set.seed(NAI*101) #set a seed value for reproduce the results
N <- 5000 # number of samples from the prior
prior <- rmvnorm(N,mu_prior,Sigma_prior) # Initial particle set

num_par <- ncol(prior)  # total number of parameters in the Copula model
r01 <- 0.001            # nugget effect parameter for model 1
r02 <- 0.01             # nugget effect parameter for model 2

###################################################
Monitoring_sites <- read.csv("Example2/Rcodes/design_selection/monitoring_sites_selected.csv") # prediction locations
Site_Locs <- data.frame(Monitoring_sites[,6:7]) 
dist1 <- rdist(Site_Locs)
dist_mat <- dist1/max(dist1)
tot_sites <- nrow(dist_mat)

St_loc=Site_Locs
St_loc[,1] <- St_loc[,1]-min(St_loc[,1])
St_loc[,1] <- St_loc[,1]/max(dist1)
St_loc[,2] <- St_loc[,2]-min(St_loc[,2])
St_loc[,2] <- St_loc[,2]/max(dist1)

Unsampling_loc <- Monitoring_sites

data2017 <- read.csv("Example2/Rcodes/design_selection/Meteorological_data2017.csv") # meteorological data at the prediction locations
X1 <- (-1*data2017$Loc_Y)/3000000
X2 <- data2017$Humidity/100
X3 <- data2017$Wind_speed/10
X.all <- cbind(X1,X2,X3)      # values of the covariates at the prediction locations

locs_uns <- which(Monitoring_sites$Site_No %in% Unsampling_loc$Site_No)
distM_uns=dist_mat[locs_uns,locs_uns]
ST_uns <- St_loc[locs_uns,]
X_unsamp <- X.all[locs_uns,]

set.seed(101)
B1 <- 100
theta_int <- rmvnorm(n=B1,mean=mu_prior, sigma = Sigma_prior)

B2 <- 1000
Z1 <- matrix(rnorm(tot_sites*B2), tot_sites)
Unif_set <- runif(tot_sites,0,1)  
Z_res <- matrix(rnorm(tot_sites*3),ncol=3)

Des=c("dual","est","pred")
d.spatial <- list()

for(i in 1:3){
  d.name=paste("Example2/Selected_design/",Des[i],"_d",d_no,".RData",sep="")
  load(d.name)
  d.opt = out[[2]]
  Crit.out <- exp.crit.all(d=d.opt,B=B1)
  d.spatial[[i]] <- Crit.out
}

v_name=paste("out",NAI,"D",d_no,".RData",sep="")
save(d.spatial, file = v_name)

