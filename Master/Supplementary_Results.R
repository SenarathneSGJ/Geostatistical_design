#Required R packages
library(ggplot2)
library(ggpubr)
library(wesanderson)

##### Design plot in Example 2 ######
# Figure S1
Clust <- c("C3")
pl_name <- c("a","b")
loc.d <- c(5,7)
png("Results/Figures_and_Tables/Figure_S1.png")
par(oma=c(2,1,1,1),mar=c(6,4,1,1))
par(mfrow=c(1,2))

Monitoring_sites <- read.csv("Data/monitoring_sites_selected.csv")
C1C2_sites <- c(1:3,5,6,9,26,28,29,31)
sel.sites <- which(!(Monitoring_sites$Site_No %in% C1C2_sites))
Site_Locs <- Monitoring_sites[sel.sites,]

for(i in 1:2){
  path.dE <- paste("Results/Suplementary_results/C3_simulations/est_d",loc.d[i],".RData",sep="")
  load(path.dE)
  Est.Opt <- out[[2]]
  
  path.dB <- paste("Results/Suplementary_results/C3_simulations/dual_d",loc.d[i],".RData",sep="")
  load(path.dB)
  EstnPred.Opt <- out[[2]]
  
  path.dP <- paste("Results/Suplementary_results/C3_simulations/pred_d",loc.d[i],".RData",sep="")
  load(path.dP)
  Pred.Opt <- out[[2]]
  
  opt_loc.E <- Site_Locs[Site_Locs$Site_No %in% Est.Opt,]
  opt_loc.B <- Site_Locs[Site_Locs$Site_No %in% EstnPred.Opt,]
  opt_loc.P <- Site_Locs[Site_Locs$Site_No %in% Pred.Opt,]
  
  plot(opt_loc.B$Loc_X/10000,opt_loc.B$Loc_Y/10000,pch=18,cex.sub=0.8,cex.lab=0.8,cex.axis=0.8,col="green",cex=2.5,
		xlab="X-coordinate",ylab="Y-coordinate",sub=paste("(" ,pl_name[i],") C3,N=",loc.d[i],sep=""),
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

#### Design Evaluation ####
n1 <- 50
loc.d <- c(5,7)
Design <- c("Estimation","Prediction","Dual-purpose")

Est_List <- list()  # -L(d,y) for estimation
Pred_List <- list() # -L(d,y) for prediction

Entro_pred_prior <- (-58.64) # Prediction entropy (H(Z,theta|xi)) based on prior

for(j in 1:2)
{
  E.est <- c();P.est <- c();D.est <- c()
  E.pred <- c();P.pred <- c();D.pred <- c()
  for(i in 1:n1){
    vec <- paste("Results/Suplementary_results/C3_simulations/d",loc.d[j],"/out",i,".RData",sep="")
    load(vec)
    
    ED <- unlist(d.spatial[[1]])
    E.est[((i-1)*1000+1):(1000*i)] <- ED[seq(1,3000,3)]
    E.pred[((i-1)*1000+1):(1000*i)] <- ED[seq(2,3000,3)]-Entro_pred_prior
    
    PD <- unlist(d.spatial[[2]])
    P.est[((i-1)*1000+1):(1000*i)] <- PD[seq(1,3000,3)]
    P.pred[((i-1)*1000+1):(1000*i)] <- PD[seq(2,3000,3)]-Entro_pred_prior
    
    BD <- unlist(d.spatial[[3]])
    D.est[((i-1)*1000+1):(1000*i)] <- BD[seq(1,3000,3)]
    D.pred[((i-1)*1000+1):(1000*i)] <- BD[seq(2,3000,3)]-Entro_pred_prior
  }
  DE1 <- paste0("d",loc.d[j],"_est_E")
  DP1 <- paste0("d",loc.d[j],"_pred_E")
  DB1 <- paste0("d",loc.d[j],"_dual_E")
  
  Est_List[[((j-1)*3+1)]] <- E.est
  Est_List[[((j-1)*3+2)]] <- P.est
  Est_List[[3*j]] <- D.est
  names(Est_List) <- c(DE1,DP1,DB1)
  
  DE2 <- paste0("d",loc.d[j],"_est_P")
  DP2 <- paste0("d",loc.d[j],"_pred_P")
  DB2 <- paste0("d",loc.d[j],"_dual_P")
  
  Pred_List[[((j-1)*3+1)]] <- E.pred
  Pred_List[[((j-1)*3+2)]] <- P.pred
  Pred_List[[3*j]] <- D.pred
  names(Pred_List)=c(DE2,DP2,DB2)
}

# obtain expected loss values 
N=50 # number of simulations
Est_Loss <- data.frame(matrix(ncol=3,nrow=N*6))  # L(d) for estimation
names(Est_Loss)=c("Design","Design_points","Expected_loss")
Pred_Loss <- data.frame(matrix(ncol=3,nrow=N*6))  # L(d) for prediction
names(Pred_Loss)=c("Design","Design_points","Expected_loss")

d.points <- paste("n=",rep(c(5,7),each=3),sep="")
Des=rep(c("Estimation","Prediction","Dual-purpose"),2)

for(i in 1:6){
  Design=rep(Des[i],N)
  Design_points=rep(d.points[i],N)
  E_loss=P_loss=c()
  for(j in 1:N){
    set.seed(j*101)
    E_loss[j] <- (-1)*mean(Est_List[[i]][((j-1)*1000+1):(1000*j)],na.rm=T)
    P_loss[j] <- (-1)*mean(Pred_List[[i]][((j-1)*1000+1):(1000*j)],na.rm=T)
  }
  Est.summary=data.frame(Design,Design_points,Expected_loss=E_loss)
  Est_Loss[((i-1)*N+1):(i*N),]=Est.summary
  
  Pred.summary=data.frame(Design,Design_points,Expected_loss=P_loss)
  Pred_Loss[((i-1)*N+1):(i*N),]=Pred.summary
}

# Figure S2
png("Results/Figures_and_Tables/Figure_S2.png")

p1 <- ggplot(Est_Loss,aes(x=factor(Design_points,level=c('n=5','n=7')),y=Expected_loss,
		fill=factor(Design,level=c('Estimation','Dual-purpose','Prediction'))))+
  geom_boxplot(outlier.shape=NA)
p1+labs(x="Number of design points",y=bquote('Expected loss for estimation'~~(lambda[E])),fill='Design')+
  scale_fill_manual(values=c('#99CCFF','#3399FF','#CCFFFF'))+
  theme(legend.position="bottom",legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(0,0,0,0))
dev.off()

# Figure S3
png("Results/Figures_and_Tables/Figure_S3.png")
p2 <- ggplot(Pred_Loss,aes(x=factor(Design_points,level=c('n=5','n=7','n=10','n=15')),y=Expected_loss,
		fill=factor(Design,level=c('Estimation','Dual-purpose','Prediction'))))+
  geom_boxplot(outlier.shape=NA)
p2+labs(x="Number of design points",y=bquote('Expected loss for prediction'~~(lambda[P])),fill='Design')+
  scale_fill_manual(values=c('#99CCFF','#3399FF','#CCFFFF'))+
theme(legend.position="bottom",legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(0,0,0,0))
dev.off()

#### Design efficiency ####
#Table S1 
n.points <- paste("n=",rep(c(5,7),each=3),sep="")
Design=rep(c("Estimation","Prediction","Dual-purpose"),2)

Estimation_loss=Prediction_loss=c()
for(k in 1:6){
  Estimation_loss[k]=(-1)*mean(Est_List[[k]],na.rm=T)
  Prediction_loss[k]=(-1)*mean(Pred_List[[k]],na.rm=T)
}

Est.vec=rep(Estimation_loss[c(1,4)],each=3)
Pred.vec=rep(Prediction_loss[c(2,5)],each=3)

Est.efficiency=1-((Est.vec-Estimation_loss)/Est.vec)
Pred.efficiency=1-((Pred.vec-Prediction_loss)/Pred.vec)

Data.loss=data.frame(Design,n.points,Est.efficiency,Pred.efficiency)
write.csv(Data.loss,"Results/Figures_and_Tables/Table_S1.csv")

###### Compare two approximations for prediction loss ###

#Figure S4
Ute.Approx1=c()
Ute.Approx2=c()
for(i in 1:50){
  d.name=paste("Results/Suplementary_results/Efficiency_prediction_loss_approximations/Example1/out",i,".RData",sep="")
  load(d.name)
  Ute.Approx1[i]=as.numeric(crit.out[[7]][3])
  Ute.Approx2[i]=as.numeric(crit.out[[8]][3])
}

Ute=c(rep("ld.tilde",50),rep("ld.hat",50))
val.s=c(Ute.Approx1,Ute.Approx2)
data2=data.frame(Ute,val.s)

png("Results/Figures_and_Tables/Figure_S4.png")
p <- ggplot(data2,aes(x=Ute,y=val.s)) + 
  geom_boxplot()
p+ylab("Time (s)")+xlab("Approximated Loss Function")+
  scale_x_discrete(labels=c(bquote( ~hat(lambda)[P](d,y)),bquote( ~~tilde(lambda)[P](d,y))))
dev.off()

#Figure S5
Ute.Approx1=c()
Ute.Approx2=c()
for(i in 1:50){
  d.name=paste("Results/Suplementary_results/Efficiency_prediction_loss_approximations/Example2/out",i,".RData",sep="")
  load(d.name)
  Ute.Approx1[i]=as.numeric(crit.out[[7]][1])
  Ute.Approx2[i]=as.numeric(crit.out[[8]][1])
}

Ute=c(rep("ld.tilde",50),rep("ld.hat",50))
val.s=c(Ute.Approx1,Ute.Approx2)
data2=data.frame(Ute,val.s)
png("Results/Figures_and_Tables/Figure_S5.png")
p <- ggplot(data2,aes(x=Ute,y=val.s)) + 
  geom_boxplot()
p+ylab("Time (s)")+xlab("Approximated Loss Function")+
  scale_x_discrete(labels=c(bquote( ~hat(lambda)[P](d,y)),bquote( ~~tilde(lambda)[P](d,y))))
dev.off()

####### Compare theoretical design #######
Unsamp_Lx<-rep(seq(0,1,.25),5)
Unsamp_Ly<-rep(seq(0,1,.25),each=5)
data_UnS=data.frame(Unsamp_Lx,Unsamp_Ly)

# Equally spaced triangular design
x1=c(0,1/6,1/3,1/3,1/2,1/2,2/3,2/3,5/6,1)
y1=c(0,1/3,0,2/3,1/3,1,0,2/3,1/3,0)
d.type=c(rep("d",length(x1)),rep("p",length(Unsamp_Lx)))
X1=c(x1,Unsamp_Lx)
X2=c(y1,Unsamp_Ly)
EqSp.D=data.frame(d.type,X1,X2)

# A design with all the points on the boundary
x2=c(0,.5,1,1,0,1,0,0,.5,1)
y2=c(0,0,0.25,0,.25,.75,.75,1,1,1)
d.type=c(rep("d",length(x2)),rep("p",length(Unsamp_Lx)))
X1 <- c(x2,Unsamp_Lx)
X2 <- c(y2,Unsamp_Ly)
Bnd.D <-data.frame(d.type,X1,X2)

# A design with all points close to the prediction locations 
x3 <- c(0,.5,.75,1,0,1,0.25,0,.5,1)
y3 <- c(0,0,0.5,0,.5,.5,.5,1,1,1)
d.type <- c(rep("d",length(x3)),rep("p",length(Unsamp_Lx)))
X1 <- c(x3,Unsamp_Lx)
X2 <- c(y3,Unsamp_Ly)
PrC.D <- data.frame(d.type,X1,X2)

png("Results/Figures_and_Tables/Figure_S6.png")
p1 <- ggplot(data=EqSp.D,aes(X1,X2,col=d.type,shape=d.type))+xlab("X1")+ylab("X2")+
  labs(caption="(a) An equally spaced triangular design (Equal.SP)")+
  theme(plot.caption=element_text(hjust=0.5))+
  geom_point(size=3)+
  scale_shape_manual(labels=c("Design points","Prediction locations"),values=c(4,1))+
  scale_color_manual(labels=c("Design points","Prediction locations"),values=c("blue","black"))+
  theme(legend.title=element_blank())

p2 <- ggplot(data=Bnd.D,aes(X1,X2,col=d.type,shape=d.type))+xlab("X1")+ylab("X2")+
  labs(caption="(b)	 A design with all the points on the boundary (Boundary.P)")+
  theme(plot.caption=element_text(hjust=0.5))+
  geom_point(size=3)+
  scale_shape_manual(labels=c("Design points","Prediction locations"),values=c(4,1))+
  scale_color_manual(labels=c("Design points","Prediction locations"),values=c("blue","black"))+
  theme(legend.title=element_blank())

p3 <- ggplot(data=PrC.D,aes(X1,X2,col=d.type,shape=d.type))+xlab("X1")+ylab("X2")+
  labs(caption="(c) A design with all points close to the prediction locations (Close.Pred)")+
  theme(plot.caption=element_text(hjust=0.5))+
  geom_point(size=3)+
  scale_shape_manual(labels=c("Design points","Prediction locations"),values=c(4,1))+
  scale_color_manual(labels=c("Design points","Prediction locations"),values=c("blue","black"))+
  theme(legend.title=element_blank())

# Figure S6
ggarrange(p1,p2,p3,ncol=2,nrow=2,common.legend=TRUE,legend="bottom")
dev.off()

# Figure S7
load("Results/Suplementary_results/Compare_common_designs/Post_pred_var_Y1.RData")
png("Results/Figures_and_Tables/Figure_S7.png")
p1 <- ggplot(data.new.Y1,aes(x=Y1.var,y=Y1.Eqsp.P1)) +
  geom_point(col="blue") + ylim(0,175)+
  labs(subtitle="Equally.SP",x="prior pred. var.",y="mean(Posterior pred.var)")+
  geom_text(label=data.new.Y1$site_no,nudge_x=11,cex=3)

p2 <- ggplot(data.new.Y1,aes(x=Y1.var,y=Y1.Bnd.P1)) +
  geom_point(col="blue") + ylim(0,175)+
  labs(subtitle="Boundary.P",x="prior pred. var.",y="mean(Posterior pred.var)")+
  geom_text(label=data.new.Y1$site_no,nudge_x=11,cex=3)

p3 <- ggplot(data.new.Y1,aes(x=Y1.var,y=Y1.PrC.P1)) +
  geom_point(col="blue") + ylim(0,175)+
  labs(subtitle="Close.Pred",x="prior pred. var.",y="mean(Posterior pred.var)")+
  geom_text(label=data.new.Y1$site_no,nudge_x=11,cex=3)

p4 <- ggplot(data.new.Y1,aes(x=Y1.var,y=Y1.D.P1)) +
  geom_point(col="blue") + ylim(0,175)+ 
  labs(subtitle="Dual.Purpose",x="prior pred. var.",y="mean(Posterior pred.var)")+
  geom_text(label=data.new.Y1$site_no,nudge_x=11,cex=3)

p5 <- ggplot(data.new.Y1,aes(x=Y1.var,y=Y1.E.P1)) +
  geom_point(col="blue") + ylim(0,175)+ 
  labs(subtitle="Estimation",x="prior pred. var.",y="mean(Posterior pred.var)")+
  geom_text(label=data.new.Y1$site_no,cex=3,nudge_x=11)

p6 <- ggplot(data.new.Y1,aes(x=Y1.var,y=Y1.P.P1)) +
  geom_point(col="blue") + ylim(0,175)+ 
  labs(subtitle="Prediction",x="prior pred. var.",y="mean(Posterior pred.var)")+
  geom_text(label=data.new.Y1$site_no,cex=3,nudge_x=11)

ggarrange(p1,p2,p3,p4,p5,p6,ncol=3,nrow=2)
dev.off()

# Figure S8 
load("Results/Suplementary_results/Compare_common_designs/Post_pred_var_Y2.RData")
png("Results/Figures_and_Tables/Figure_S8.png")
p7 <- ggplot(data.new.Y2,aes(x=log.Y2.var.,y=log.Y2.Eqsp.P1.)) +
  geom_point(col="blue") + ylim(-10,7.5)+
  labs(subtitle="Equally.SP",x="prior pred. var.",y="mean(Posterior pred.var)")+
  geom_text(label=data.new.Y2$site_no,nudge_x=1,nudge_y=(-.5),cex=3)

p8 <- ggplot(data.new.Y2,aes(x=log.Y2.var.,y=log.Y2.Bnd.P1.)) +
  geom_point(col="blue") + ylim(-10,7.5)+
  labs(subtitle="Boundary.P",x="prior pred. var.",y="mean(Posterior pred.var)")+
  geom_text(label=data.new.Y2$site_no,nudge_x=1,nudge_y=(-.5),cex=3)

p9 <- ggplot(data.new.Y2,aes(x=log.Y2.var.,y=log.Y2.PrC.P1.)) +
  geom_point(col="blue") + ylim(-10,7.5)+
  labs(subtitle="Close.Pred",x="prior pred. var.",y="mean(Posterior pred.var)")+
  geom_text(label=data.new.Y2$site_no,nudge_x=1,nudge_y=(-.5),cex=3)

p10 <- ggplot(data.new.Y2,aes(x=log.Y2.var.,y=log.Y2.D.P1.)) +
  geom_point(col="blue") + ylim(-10,7.5)+
  labs(subtitle="Dual.Purpose",x="prior pred. var.",y="mean(Posterior pred.var)")+
  geom_text(label=data.new.Y2$site_no,nudge_x=1,nudge_y=(-.5),cex=3)

p11 <- ggplot(data.new.Y2,aes(x=log.Y2.var.,y=log.Y2.E.P1.)) +
  geom_point(col="blue") + ylim(-10,7.5)+
  labs(subtitle="Estimation",x="prior pred. var.",y="mean(Posterior pred.var)")+
  geom_text(label=data.new.Y2$site_no,nudge_x=1,nudge_y=(-.5),cex=3)

p12 <- ggplot(data.new.Y2,aes(x=log.Y2.var.,y=log.Y2.P.P1.)) +
  geom_point(col="blue") + ylim(-10,7.5)+
  labs(subtitle="Prediction",x="prior pred. var.",y="mean(Posterior pred.var)")+
  geom_text(label=data.new.Y2$site_no,nudge_x=1,nudge_y=(-.5),cex=3)

ggarrange(p7,p8,p9,p10,p11,p12,ncol=3,nrow=2)
dev.off()

# Figure S9
load("Results/Suplementary_results/Compare_common_designs/Post_var.RData")
png("Results/Figures_and_Tables/Figure_S9.png")
par(mfrow=c(2,3),mar=c(4,4,1,1))
parameter <- as.factor(c("B01","B11","B21","SdY1","R11","R21","K1","B02","B12","B22","R12","R22","K2","log_tau"))
des.C <- c('Equally.SP','Boundary.P','Close.Pred',
        'Estimation','Dual.Purpose',"Prediction")

# Prior distribution
Dependence <- 0.2
mu_prior1 <- c(5,-2.8,8,log(1.2),log(.7/Dependence),log(Dependence),log(1.5))  
Sigma_prior1 <- c(4,4,4,0.25,0.25,0.25,0.25)
mu_prior2 <- c(3.8,-0.5,-0.7,log(.6/Dependence),log(Dependence),log(0.25))   
Sigma_prior2 <- c(0.125,0.125,0.125,0.125,0.125,0.25)
mu_Ltau <- log(.7/.3)
Sigma_Ltau <- 0.25

mu_prior<- c(mu_prior1,mu_prior2,mu_Ltau)
Sigma_prior <- diag(c(Sigma_prior1,Sigma_prior2,Sigma_Ltau))

p <- list()
for(i in 1:6){
  des <- out[[i]]
  post <- apply(des,2,mean,na.rm=T)
  data2 <- data.frame(prior=diag(Sigma_prior),post)
  
  p[[i]] <- ggplot(data2,aes(x=prior,y=post)) +
    geom_point(size=2,aes(shape=parameter,col=parameter)) + ylim(0,1.75)+
    labs(subtitle=des.C[i],x="prior_variance",y="posterior_variance")+
    scale_color_hue(labels=c('B01'=expression(beta[10]),'B11'=expression(beta[11]),
                        'B21'=expression(beta[12]),'SdY1'=expression(sigma[1]),
                        'R11'=expression(gamma[11]),'R21'=expression(gamma[12]),
                        'K1'=expression(nu[1]),'B02'=expression(beta[20]),
                        'B12'=expression(beta[21]),'B22'=expression(beta[22]),
                        'R12'=expression(gamma[21]),'R22'=expression(gamma[22]),
                        'K2'=expression(nu[2]),'log_tau'=expression(tau)))+
    scale_shape_manual(labels=c('B01'=expression(beta[10]),'B11'=expression(beta[11]),
                        'B21'=expression(beta[12]),'SdY1'=expression(sigma[1]),
                        'R11'=expression(gamma[11]),'R21'=expression(gamma[12]),
                        'K1'=expression(nu[1]),'B02'=expression(beta[20]),
                        'B12'=expression(beta[21]),'B22'=expression(beta[22]),
                        'R12'=expression(gamma[21]),'R22'=expression(gamma[22]),
                        'K2'=expression(nu[2]),'log_tau'=expression(tau)),values=parameter)
}  

ggarrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],ncol=3,nrow=2,common.legend=TRUE,legend="bottom")
dev.off()
