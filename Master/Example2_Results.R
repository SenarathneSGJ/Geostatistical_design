# Required R packages
library(ggplot2)
library(robustbase)
#Install the following package from the source file mapping_0.1.tar.gz or 
#download it from the link: https://www.flutterbys.com.au/stats/tut/tut5.4.html
library(mapping) 
library(sp)
library(ggplot2)
library(gridExtra)

## Figure 1 - Sampling/unsampling locations
Sites_all <- read.csv("Data/monitoring_sites_selected.csv")
sampling_sites <- c(6,9,24,25,26,29,34)

Unsamp_locs <- Sites_all[!(Sites_all$Site_No %in%sampling_sites),]
Site <- rep("Unsampling_sites",nrow(Unsamp_locs))
Unsamp_locs <- cbind(Unsamp_locs,Site)

Samp_locs <- Sites_all[(Sites_all$Site_No %in% sampling_sites),]
Site <- rep("Sampling_sites",nrow(Samp_locs))
Samp_locs <- cbind(Samp_locs,Site)

data_new <- rbind(Samp_locs,Unsamp_locs)

png("Results/Figures_and_Tables/Figure1.png")
P<-ggplot(fortify(qld.mainland),aes(y=lat,x=long)) + 
  geom_polygon(colour="#68694F",fill="#DFE58B")
P1 <- P+geom_point(data=data_new,mapping=aes(y=Latitude,x=Longitude,col=Site,shape=Site),size=2)+
  ylab("Latitude")+xlab("Longitude")+theme(legend.position="bottom")

P.new <- ggplot(data_new,aes(y=Latitude,x=Longitude,label=Site_No))  +
  geom_text(hjust=-0.8-runif(1,0,.3),vjust=0.7+runif(1,0,.3),size=1.8)+
  geom_point(mapping=aes(col=Site,shape=Site),size=2) +
  theme(axis.text.x=element_blank(),axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),axis.title.y=element_blank(),
        legend.position="none",panel.background=element_rect(fill="#DFE58B"))

C1 <- P.new+coord_cartesian(xlim=c(146.8,148.1),ylim=c(-20,-19.2))
C2 <- P.new+coord_cartesian(xlim=c(151.05,151.4),ylim=c(-23.95,-23.8))
C3 <- P.new+coord_cartesian(xlim=c(152.3,153.25),ylim=c(-28,-27.3))

grid.arrange(P1,C1,C2,C3,ncol=3,nrow=3,
             layout_matrix=cbind(c(1,1,1),c(1,1,1),c(2,3,4)))
dev.off()

### Compare two approximations for prediction loss function ###
# Figure 7
vt <- 1:100 #number of simulations
pred_approx <- c()
pred_act <- c()
for(i in 1:length(vt)){
  loc <- paste("Results/Example2/Simulation_results/Pred_Ute_approximations/out",vt[i],".RData",sep="")
  load(loc)
  ute_mat <- matrix(out.data,ncol=4,byrow=T)
  pred_approx <- c(pred_approx,ute_mat[,3])
  pred_act <- c(pred_act,ute_mat[,4])
}

png("Results/Figures_and_Tables/Figure7.png")
plot(pred_act,pred_approx,ylab=bquote("Loss values obtained from " ~~tilde(lambda)[P]),
	xlab=bquote("Loss values obtained from " ~hat(lambda)[P]),cex=0.5)
dev.off()

##### Selected Designs (Table 4) ######

loc.d <- c(5,7,10,15)
d1 <- matrix(nrow=4,ncol=3)

for(i in 1:4){
  path.dE <- paste("Results/Example2/Selected_design/est_d",loc.d[i],".RData",sep="")
  load(path.dE)
  d1[i,1] <- paste(sort(out[[2]]),collapse=',')
  
  path.dB <- paste("Results/Example2/Selected_design/dual_d",loc.d[i],".RData",sep="")
  load(path.dB)
  d1[i,2] <- paste(sort(out[[2]]),collapse=',')
  
  path.dP <- paste("Results/Example2/Selected_design/pred_d",loc.d[i],".RData",sep="")
  load(path.dP)
  d1[i,3] <- paste(sort(out[[2]]),collapse=',')
}  

colnames(d1) <- c("Estimation","Dual-purpose","Prediction")
write.csv(d1,"Results/Figures_and_Tables/Table4.csv")

##### Design plot in Example 2 ######
# Figure 8
Clust <- c("C1","C2","C3")
pl_name <- c("a","b","c","d","e","f","g","h","i","j","k","l")
loc.d <- c(5,7,10,15)
png("Results/Figures_and_Tables/Figure8.png")
par(oma=c(4,1,1,1),mar=c(6,4,1,1))
par(mfrow=c(4,3))

Site_Locs <- read.csv("Data/monitoring_sites_selected.csv")

for(i in 1:4){
  path.dE <- paste("Results/Example2/Selected_design/est_d",loc.d[i],".RData",sep="")
  load(path.dE)
  Est.Opt <- out[[2]]
  
  path.dB <- paste("Results/Example2/Selected_design/dual_d",loc.d[i],".RData",sep="")
  load(path.dB)
  EstnPred.Opt <- out[[2]]
  
  path.dP <- paste("Results/Example2/Selected_design/pred_d",loc.d[i],".RData",sep="")
  load(path.dP)
  Pred.Opt <- out[[2]]
  
  opt_loc.E <- Site_Locs[Site_Locs$Site_No %in% Est.Opt,]
  opt_loc.B <- Site_Locs[Site_Locs$Site_No %in% EstnPred.Opt,]
  opt_loc.P <- Site_Locs[Site_Locs$Site_No %in% Pred.Opt,]
  
  plot(opt_loc.B$Loc_X/10000,opt_loc.B$Loc_Y/10000,pch=18,cex.sub=0.8,cex.lab=0.8,cex.axis=0.8,col="green",cex=2.5,
		xlab="X-coordinate",ylab="Y-coordinate",sub=paste("(" ,pl_name[(i-1)*3+1],") C1,N=",loc.d[i],sep=""),
		ylim=range(-223,(max(Site_Locs$Loc_Y/10000)+1)),xlim=range(45,65))
  points(opt_loc.E$Loc_X/10000,opt_loc.E$Loc_Y/10000,pch=3,col="red",cex=2)
  points(opt_loc.P$Loc_X/10000,opt_loc.P$Loc_Y/10000,pch=1,col="blue",cex=2)
  
  plot(opt_loc.B$Loc_X/10000,opt_loc.B$Loc_Y/10000,pch=18,cex.sub=0.8,cex.lab=0.8,cex.axis=0.8,col="green",cex=2.5,
		xlab="X-coordinate",ylab="Y-coordinate",sub=paste("(" ,pl_name[(i-1)*3+2],") C2,N=",loc.d[i],sep=""),
		ylim=range(-266,-263),xlim=range(91,95))
  points(opt_loc.E$Loc_X/10000,opt_loc.E$Loc_Y/10000,pch=3,col="red",cex=2)
  points(opt_loc.P$Loc_X/10000,opt_loc.P$Loc_Y/10000,pch=1,col="blue",cex=2)
  
  plot(opt_loc.B$Loc_X/10000,opt_loc.B$Loc_Y/10000,pch=18,cex.sub=0.8,cex.lab=0.8,cex.axis=0.8,col="green",cex=2.5,
		xlab="X-coordinate",ylab="Y-coordinate",sub=paste("(" ,pl_name[(i-1)*3+3],") C3,N=",loc.d[i],sep=""),
		ylim=range(-315,-300),xlim=range(102,112))
  points(opt_loc.E$Loc_X/10000,opt_loc.E$Loc_Y/10000,pch=3,col="red",cex=2)
  points(opt_loc.P$Loc_X/10000,opt_loc.P$Loc_Y/10000,pch=1,col="blue",cex=2)
  
}

par(fig=c(0,1,0,1),oma=c(0.1,0,0,0),mar=c(0.2,2.5,0,2.5),new=TRUE)
plot(0,0,type="n",bty="n",xaxt="n",yaxt="n")
pch_types <- c(3,18,1)
text1 <- c("Estimation","Dual-purpose","Prediction")
plot.col <- c("red","green","blue")
legend("bottom",text1,x.intersp=1.8,ncol=3,text.width=strwidth("1,000,000,000,0000"),
       col=plot.col,pch=pch_types,pt.cex=c(2,2.5,2),cex=.8)
dev.off()

#### Design Comparison ####
loc.d <- c(5,7,10,15)
Des <- c("Estimation","Dual-purpose","Prediction")
Est_List <- list()  # -L(d,y) for estimation
Pred_List <- list() # -L(d,y) for prediction

# Prediction entropy (H(Z,theta|xi)) based on prior
Entro_pred_prior <- (-100.6421) 

# obtaining expected loss values 
N <- 100 # number of simulations
Est_Loss <- data.frame(matrix(ncol=3,nrow=N*12))  # L(d) for estimation
names(Est_Loss) <- c("Design","Design_points","Expected_loss")
Pred_Loss <- data.frame(matrix(ncol=3,nrow=N*12))  # L(d) for prediction
names(Pred_Loss) <- c("Design","Design_points","Expected_loss")

for(j in 1:4)
{
  E.loss <- c();P.loss <- c()
 for(i in 1:N){
  vec <- paste("Results/Example2/Simulation_results/Design_evaluation/D",loc.d[j],"/out",i,".RData",sep="")
  load(vec)
  
  ED <- unlist(d.spatial[[2]])
  E.loss[i] <- (-1)*mean(ED[seq(1,300,3)],na.rm=T)
  P.loss[i] <- (-1)*mean(ED[seq(2,300,3)]-Entro_pred_prior,na.rm=T)
  
  BD <- unlist(d.spatial[[1]])
  E.loss[i+N] <- (-1)*mean(BD[seq(1,300,3)],na.rm=T)
  P.loss[i+N] <- (-1)*mean(BD[seq(2,300,3)]-Entro_pred_prior,na.rm=T)
  
  PD <- unlist(d.spatial[[3]])
  E.loss[i+2*N] <- (-1)*mean(PD[seq(1,300,3)],na.rm=T)
  P.loss[i+2*N] <- (-1)*mean(PD[seq(2,300,3)]-Entro_pred_prior,na.rm=T)
  
  }
  Design <- rep(Des,each=N)
  Design_points <- rep(paste("n=",loc.d[j],sep=""),N*3)
  
  Est.summary <- data.frame(Design,Design_points,Expected_loss=E.loss)
  Est_Loss[((j-1)*3*N+1):(j*3*N),] <- Est.summary
  
  Pred.summary <- data.frame(Design,Design_points,Expected_loss=P.loss)
  Pred_Loss[((j-1)*3*N+1):(j*3*N),] <- Pred.summary
}

# Figure 9 - Expected loss for estimation
png("Results/Figures_and_Tables/Figure9.png")
p1 <- ggplot(Est_Loss,aes(x=factor(Design_points,level=c('n=5','n=7','n=10','n=15')),y=Expected_loss,
		fill=factor(Design,level=c('Estimation','Dual-purpose','Prediction'))))+
		geom_boxplot(outlier.shape=NA)
p1+labs(x="Number of design points",y=bquote("Expected loss for estimation" ~~(lambda[E])),fill='Design')+
		scale_fill_manual(values=c('#99CCFF','#3399FF','#CCFFFF'))+
		theme(legend.position="bottom",legend.margin=margin(0,0,0,0),legend.box.margin=margin(0,0,0,0))
dev.off()

# Figure 10 - Expected loss for prediction
png("Results/Figures_and_Tables/Figure10.png")
p2 <- ggplot(Pred_Loss,aes(x=factor(Design_points,level=c('n=5','n=7','n=10','n=15')),y=Expected_loss,
		fill=factor(Design,level=c('Estimation','Dual-purpose','Prediction'))))+
		geom_boxplot(outlier.shape=NA)
p2+labs(x="Number of design points",y=bquote("Expected loss for prediction" ~~(lambda[P])),fill='Design')+
		scale_fill_manual(values=c('#99CCFF','#3399FF','#CCFFFF'))+
		theme(legend.position="bottom",legend.margin=margin(0,0,0,0),legend.box.margin=margin(0,0,0,0))
dev.off()

###### Design Efficiency ######
Estimation_loss <- aggregate(Est_Loss$Expected_loss,list(Est_Loss$Design_points,Est_Loss$Design),FUN=mean) 
Estimation_loss <- Estimation_loss[order(Estimation_loss$Group.1),]

Prediction_loss <- aggregate(Pred_Loss$Expected_loss,list(Pred_Loss$Design_points,Pred_Loss$Design),FUN=mean) 
Prediction_loss <- Prediction_loss[order(Prediction_loss$Group.1),]

Est.vec <- rep(Estimation_loss$x[c(2,5,8,11)],each=3)
Pred.vec <- rep(Prediction_loss$x[c(3,6,9,12)],each=3)

Est.efficiency <- 1-((Est.vec-Estimation_loss$x)/Est.vec)
Pred.efficiency <- 1-((Pred.vec-Prediction_loss$x)/Pred.vec)

#Table 5
Data.loss <- data.frame(n.points=Estimation_loss[,1],Design=Estimation_loss[,2],Est.efficiency,Pred.efficiency)
write.csv(Data.loss,"Results/Figures_and_Tables/Table5.csv")
