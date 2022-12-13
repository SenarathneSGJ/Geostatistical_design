
### Compare two approximations for prediction loss function ###
#Figure 2

vt <- 1:100 #number of simulations
pred_approx <- c()
pred_act <- c()
for(i in 1:length(vt)){
  loc <- paste("Results/Example1/Simulation_results/Pred_Ute_approximations/out",vt[i],".RData",sep="")
  load(loc)
  ute_mat <- matrix(out.data,ncol=4,byrow=T)
  pred_approx <- c(pred_approx,ute_mat[,3])
  pred_act <- c(pred_act,ute_mat[,4])
}
png("Results/Figures_and_Tables/Figure2.png")
plot(pred_act,pred_approx,xlim=c(-200,-135),ylim=c(-190,-115),
		ylab=bquote("Loss values obtained from " ~~tilde(lambda)[P]),
		xlab=bquote("Loss values obtained from " ~hat(lambda)[P]),cex=0.5)
dev.off()

### Design plots in Example1 ###
library(ggplot2)
cor.X1X2 <- c("Week","Moderate","Strong")
d.crit <- c("est","dual","pred")
cor.val <- c("R2","R5","R8")
loss_func <- c("Estimation","Dual","Prediction")

# Figure 3
png("Results/Figures_and_Tables/Figure3.png")
par(mfrow=c(3,4),mar=c(4,3,1,1),mgp=c(2,.8,0))
for(i in 1:3){
plot(c(0,1),c(0,1),ann=F,bty='n',type='n',xaxt='n',yaxt='n')
text(x=0.5,y=0.5,cor.X1X2[i],cex=1.1,col="black")

  for(j in 1:3){
    vec <- paste("Results/Example1/selected_designs/",d.crit[j],"_",cor.val[i],"_D5.RData",sep="")
    load(vec)
    plot(out.data[[1]],xlab="X1",ylab="X2",sub=loss_func[j],pch=19,
         xlim=c(0,1),ylim=c(0,1),col='blue')
  }
}
dev.off()

## Figure 4
png("Results/Figures_and_Tables/Figure4.png")
par(mfrow=c(3,4),mar=c(4,3,1,1),mgp=c(2,.8,0))
for(i in 1:3){
  plot(c(0,1),c(0,1),ann=F,bty='n',type='n',xaxt='n',yaxt='n')
  text(x=0.5,y=0.5,cor.X1X2[i],cex=1.1,col="black")
  
  for(j in 1:3){
    vec <- paste("Results/Example1/selected_designs/",d.crit[j],"_",cor.val[i],"_D10.RData",sep="")
    load(vec)
    plot(out.data[[1]],xlab="X1",ylab="X2",sub=loss_func[j],pch=19,
         xlim=c(0,1),ylim=c(0,1),col='blue')
  }
}
dev.off()

#### Design Comparison ####
library(robustbase)
loc.d <- c(5,10)
Des <- c("Estimation","Dual-purpose","Prediction")
Est_List <- list()  # -L(d,y) for estimation
Pred_List <- list() # -L(d,y) for prediction

# Prediction entropy (H(Z,theta|xi)) based on prior under Week, moderate, Strong correlation
Entro_pred_prior <- c(-178.6816,-179.2311,-179.702)

# obtain expected loss values 
N <- 100 # number of simulations
Est_Loss <- data.frame(matrix(ncol=4,nrow=N*18))  # L(d) for estimation
names(Est_Loss) <- c("Design","Design_points","Dependence","Expected_loss")
Pred_Loss <- data.frame(matrix(ncol=4,nrow=N*18))  # L(d) for prediction
names(Pred_Loss) <- c("Design","Design_points","Dependence","Expected_loss")

for(j in 1:2)
{
  for(k in 1:3){
    E.loss <- c();P.loss <- c()
    for(i in 1:N){
      vec <- paste("Results/Example1/Simulation_results/Design_evaluation/",cor.val[k],"D",
				loc.d[j],"/out",i,".RData",sep="")
      load(vec)
  
      ED <- unlist(d.spatial[[2]])
      E.loss[i] <- (-1)*mean(ED[seq(1,600,3)],na.rm=T)
      P.loss[i] <- (-1)*mean(ED[seq(2,600,3)]-Entro_pred_prior[k],na.rm=T)
      
      BD <- unlist(d.spatial[[1]])
      E.loss[i+N] <- (-1)*mean(BD[seq(1,600,3)],na.rm=T)
      P.loss[i+N] <- (-1)*mean(BD[seq(2,600,3)]-Entro_pred_prior[k],na.rm=T)
      
      PD <- unlist(d.spatial[[3]])
      E.loss[i+2*N] <- (-1)*mean(PD[seq(1,600,3)],na.rm=T)
      P.loss[i+2*N] <- (-1)*mean(PD[seq(2,600,3)]-Entro_pred_prior[k],na.rm=T)
      
    }
    Design <- rep(Des,each=N)
    Design_points <- rep(paste("n <- ",loc.d[j],sep=""),3*N)
    Dependence <- rep(cor.X1X2[k],3*N)
    
    Est.summary <- data.frame(Design,Design_points,Dependence,Expected_loss=E.loss)
    Est_Loss[(((k-1)*3*N)+(900*(j-1))+1):((k*3*N)+(j-1)*900),] <- Est.summary
    
    Pred.summary <- data.frame(Design,Design_points,Dependence,Expected_loss=P.loss)
    Pred_Loss[(((k-1)*3*N)+(900*(j-1))+1):((k*3*N)+(j-1)*900),] <- Pred.summary
  }
}

# Figure 5- Expected loss for estimation
png("Results/Figures_and_Tables/Figure5.png")
p1 <- ggplot(Est_Loss,aes(x=factor(Design_points,level=c('n=5','n=10')),y=Expected_loss,
			fill=factor(Design,level=c('Estimation','Dual-purpose','Prediction'))))+
    geom_boxplot(outlier.shape=NA) + facet_wrap(.~factor(Dependence,level=c("Week","Moderate","Strong")),scales="free")
p1+labs(x="Number of design points",y=bquote("Expected loss for estimation" ~~(lambda[E])),fill='Design')+
  scale_fill_manual(values=c('#99CCFF','#3399FF','#CCFFFF'))+scale_y_continuous(limits=c(-23.5,-10))+ 
  theme(legend.position="bottom",legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,0))
dev.off()

# Figure 6 - Expected loss for prediction
png("Results/Figures_and_Tables/Figure6.png")
p2 <- ggplot(Pred_Loss,aes(x=factor(Design_points,level=c('n=5','n=10')),y=Expected_loss,
		fill=factor(Design,level=c('Estimation','Dual-purpose','Prediction'))))+geom_boxplot(outlier.shape=NA)+ 
		facet_wrap(.~factor(Dependence,level=c("Week","Moderate","Strong")),scales="free")
p2+labs(x="Number of design points",y=bquote("Expected loss for prediction" ~~(lambda[P])),fill='Design')+
	scale_fill_manual(values=c('#99CCFF','#3399FF','#CCFFFF'))+scale_y_continuous(limits=c(-27.0,0))+ 
	theme(legend.position="bottom",legend.margin=margin(0,0,0,0),legend.box.margin=margin(0,0,0,0))
dev.off()

###### Design Efficiency ######
Estimation_loss <- aggregate(Est_Loss$Expected_loss,list(Est_Loss$Design_points,Est_Loss$Dependence,
					Est_Loss$Design),FUN=mean) 
Estimation_loss <- Estimation_loss[order(Estimation_loss[,1],Estimation_loss[,2],Estimation_loss[,3],decreasing=T),]

Prediction_loss <- aggregate(Pred_Loss$Expected_loss,list(Pred_Loss$Design_points,Pred_Loss$Dependence,
					Pred_Loss$Design),FUN=mean) 
Prediction_loss <- Prediction_loss[order(Prediction_loss[,1],Prediction_loss[,2],Prediction_loss[,3],decreasing=T),]

Est.vec <- rep(Estimation_loss$x[c(2,5,8,11,14,17)],each=3)
Pred.vec <- rep(Prediction_loss$x[c(1,4,7,10,13,16)],each=3)

Est.efficiency <- 1-((Est.vec-Estimation_loss$x)/Est.vec)
Pred.efficiency <- 1-((Pred.vec-Prediction_loss$x)/Pred.vec)

#Table3
Data.loss <- data.frame(n.points=Estimation_loss[,1],Dependece=Estimation_loss[,2],Design=Estimation_loss[,3],
				Est.efficiency,Pred.efficiency)
write.csv(Data.loss,"Results/Figures_and_Tables/Table3.csv")

