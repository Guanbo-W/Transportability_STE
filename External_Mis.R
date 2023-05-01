##### Simulation code for subgroup transportability
# codes for obtaining the balancing intercepts are borrowed from https://github.com/pzivich/publications-code/blob/master/ReBalancingIntercept/balance_intercept.R 
EST=function(n, n_m, x_tilde){
  Data=Gen(n, n_m)
  Y_R=Data$Y[Data$R==1] 
  R=Data$R
  n_internal=dim(Data[R==1,])[1]
  ##### nuisance parameters for psi_{1}(x_tilde):E(Y1|X1=x_tilde,R=0)
  inverse_weight=(length(which(Data$X1==x_tilde& Data$R==0))/n)^(-1)
  # mu=E(Y|A=1,X, R=1)
  model_mu=gam(Y~A+as.factor(X1)+s(X2)+s(X3)+s(X4)+s(X5)+s(X6)+s(X7)+lo(X8)+lo(X9)+lo(X10)+A*as.factor(X1)+A*X2+A*X3+A*X4+A*X5,data=Data[R==1,])  
  newdata=data.frame(rep(1,n),as.factor(Data$X1),Data$X2,Data$X3,Data$X4,Data$X5,Data$X6,Data$X7,Data$X8,Data$X9,Data$X10,as.factor(Data$X1),Data$X2,Data$X3,Data$X4,Data$X5)
  names(newdata)=c("A","X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","A1X1","A1X2","A1X3","A1X4","A1X5") 
  g=predict.Gam(model_mu,newdata = newdata)
  
  mis_model_mu=lm(Y~A+as.factor(X1)+X2+X3+X4+X5,data=Data[R==1,])  
  newdata=data.frame(rep(1,n),as.factor(Data$X1),Data$X2,Data$X3,Data$X4,Data$X5)
  names(newdata)=c("A","X1","X2","X3","X4","X5") 
  mis_g=predict(mis_model_mu,newdata = newdata)
  
  # eta=Pr(A=1|X, R=1)=sum_{s=1}^{3}Pr(A=a|X, S=s, R=1)Pr(S=s|X, R=1)=model_a*model_s
  model_s=multinom(as.factor(S)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,data=Data[R==1,],trace=F)
  s_prob=predict(model_s,newdata=Data[R==1,],"probs")
  q=s_prob[,1]
  model_a1=glm(A~as.factor(X1)+X2+X3+X4+X5+X6+X7+X8+X9+X10,family=binomial,data=Data[Data$S==1 & R==1,]) # not as.factor(S) nor as.factor(X1), a bit cheating here, but align with the data generating mechnism
  model_a2=glm(A~as.factor(X1)+X2+X3+X4+X5+X6+X7+X8+X9+X10,family=binomial,data=Data[Data$S==2 & R==1,])
  model_a3=glm(A~as.factor(X1)+X2+X3+X4+X5+X6+X7+X8+X9+X10,family=binomial,data=Data[Data$S==3 & R==1,])
  A.pre_s1=predict.glm(model_a1,newdata=Data[R==1,],type="response")
  A.pre_s2=predict.glm(model_a2,newdata=Data[R==1,],type="response")
  A.pre_s3=predict.glm(model_a3,newdata=Data[R==1,],type="response")
  eta=A.pre_s1*s_prob[,1]+A.pre_s2*s_prob[,2]+A.pre_s3*s_prob[,3]
  e=eta
  
  mis_model_s=multinom(as.factor(S)~X1+X2+X3+X4+X5,data=Data[R==1,],trace=F)
  mis_s_prob=predict(mis_model_s,newdata=Data[R==1,],"probs")
  mis_q=mis_s_prob[,1]
  mis_model_a1=glm(A~as.factor(X1)+X2+X3+X4+X5,family=binomial,data=Data[Data$S==1 & R==1,]) # not as.factor(S) nor as.factor(X1), a bit cheating here, but align with the data generating mechnism
  mis_model_a2=glm(A~as.factor(X1)+X2+X3+X4+X5,family=binomial,data=Data[Data$S==2 & R==1,])
  mis_model_a3=glm(A~as.factor(X1)+X2+X3+X4+X5,family=binomial,data=Data[Data$S==3 & R==1,])
  mis_A.pre_s1=predict.glm(mis_model_a1,newdata=Data[R==1,],type="response")
  mis_A.pre_s2=predict.glm(mis_model_a2,newdata=Data[R==1,],type="response")
  mis_A.pre_s3=predict.glm(mis_model_a3,newdata=Data[R==1,],type="response")
  mis_eta=mis_A.pre_s1*mis_s_prob[,1]+mis_A.pre_s2*mis_s_prob[,2]+mis_A.pre_s3*mis_s_prob[,3]
  mis_e=mis_eta
  
  # p=Pr(R=1|X)
  model_p=glm(R~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,data=Data,family = "binomial")
  p=predict.glm(model_p,newdata=Data[R==1,],type="response")
  
  mis_model_p=glm(R~X1+X2+X3+X4+X5,data=Data,family = "binomial")
  mis_p=predict.glm(mis_model_p,newdata=Data[R==1,],type="response")
  
  ##### estimation
  
  est_phi=inverse_weight/n*(sum(g[which(Data$X1==x_tilde& Data$R==0)])+sum(second_piece[which(Data[R==1,]$X1==x_tilde& Data[R==1,]$A==1)]))
  ##### estimation
  G_comp_v=mean(g[which(Data$X1==x_tilde& Data$R==0)])
  G_comp_x=mean(mis_g[which(Data$X1==x_tilde& Data$R==0)])
  IPTW_core_v=((1-p)/(e*p))*Y_R
  IPTW_core_x=((1-mis_p)/(mis_e*mis_p))*Y_R
  IPTW_v=inverse_weight/n*(sum(IPTW_core_v[which(Data[R==1,]$X1==x_tilde& Data[R==1,]$A==1)]))
  IPTW_x=inverse_weight/n*(sum(IPTW_core_x[which(Data[R==1,]$X1==x_tilde& Data[R==1,]$A==1)]))
  second_piece_v_v=((1-p)/(e*p))*(Y_R-g[which(Data$R==1)])
  second_piece_v_x=((1-mis_p)/(mis_e*mis_p))*(Y_R-g[which(Data$R==1)])
  second_piece_x_v=((1-p)/(e*p))*(Y_R-mis_g[which(Data$R==1)])
  second_piece_x_x=((1-mis_p)/(mis_e*mis_p))*(Y_R-mis_g[which(Data$R==1)])
  est_phi_muv_qv=inverse_weight/n*(sum(g[which(Data$X1==x_tilde& Data$R==0)])+sum(second_piece_v_v[which(Data[R==1,]$X1==x_tilde& Data[R==1,]$A==1)]))
  est_phi_muv_qx=inverse_weight/n*(sum(g[which(Data$X1==x_tilde& Data$R==0)])+sum(second_piece_v_x[which(Data[R==1,]$X1==x_tilde& Data[R==1,]$A==1)]))
  est_phi_mux_qv=inverse_weight/n*(sum(mis_g[which(Data$X1==x_tilde& Data$R==0)])+sum(second_piece_x_v[which(Data[R==1,]$X1==x_tilde& Data[R==1,]$A==1)]))
  est_phi_mux_qx=inverse_weight/n*(sum(mis_g[which(Data$X1==x_tilde& Data$R==0)])+sum(second_piece_x_x[which(Data[R==1,]$X1==x_tilde& Data[R==1,]$A==1)]))
  ##### Results
  results=c(est_phi_muv_qv,G_comp_v,IPTW_v,
            est_phi_muv_qx,G_comp_v,IPTW_x,
            est_phi_mux_qv,G_comp_x,IPTW_v,
            est_phi_mux_qx,G_comp_x,IPTW_x)
  return(results)
}

###########------- get the plot for one sample size
n=10000
#n=100000
n_m=1000
# n_m=2000 
# n_m=5000

true3=GetTrue(n,n_m,3)

est=NULL
for (i in 1:nsim) {
  temp=EST(n, n_m, 3)
  est=rbind(est,temp)
}

bias_mu_q_eta=colMeans(est[,c(1:3)])-true3[2] # bias of DR, G_comp, IPTW
sd_mu_q_eta=apply(est[,c(1:3)],2,sd) # Monte-Carlo sd of DR, G_comp, IPTW

bias_mu=colMeans(est[,c(4:6)])-true3[2] # bias of DR, G_comp, IPTW
sd_mu=apply(est[,c(4:6)],2,sd) # Monte-Carlo sd of DR, G_comp, IPTW

bias_q_eta=colMeans(est[,c(7:9)])-true3[2] # bias of DR, G_comp, IPTW
sd_q_eta=apply(est[,c(7:9)],2,sd) # Monte-Carlo sd of DR, G_comp, IPTW

bias_wrong=colMeans(est[,c(10:12)])-true3[2] # bias of DR, G_comp, IPTW
sd_wrong=apply(est[,c(10:12)],2,sd) # Monte-Carlo sd of DR, G_comp, IPTW

#### plots
plot( c(1,4), las = 2, xlim =c(0.6, 4.4), ylim = c(-2, 2),  type="n",
      xlab=TeX(r'(Correct Model(s))',bold=TRUE), ylab=TeX(r'(Bias$\pm$SD)', bold=TRUE), 
      main="Comparison of different model performance, subgroup=3, n_m=1000", xaxt="n")
axis(1, at=1:4, labels=c(TeX(r'($(\mu, \eta, q)$)', bold=FALSE),TeX(r'($(\mu)$)', bold=FALSE),TeX(r'($(\eta, q)$)', bold=FALSE),TeX(r'(None)', bold=FALSE)) )
for (i in 1:3){
  points(0.8+0.1*i,bias_mu_q_eta[i],pch=14+i, cex=2)
  segments(0.8+0.1*i, bias_mu_q_eta[i]-sd_mu_q_eta[i], 0.8+0.1*i, bias_mu_q_eta[i]+sd_mu_q_eta[i], lwd = 2 )
}
for (i in 1:3){
  points(1.8+0.1*i,bias_mu[i],pch=14+i, cex=2)
  segments(1.8+0.1*i, bias_mu[i]-sd_mu[i], 1.8+0.1*i, bias_mu[i]+sd_mu[i], lwd = 2 )
}
for (i in 1:3){
  points(2.8+0.1*i,bias_q_eta[i],pch=14+i, cex=2)
  segments(2.8+0.1*i, bias_q_eta[i]-sd_q_eta[i], 2.8+0.1*i, bias_q_eta[i]+sd_q_eta[i], lwd = 2 )
}
for (i in 1:3){
  points(3.8+0.1*i,bias_wrong[i],pch=14+i, cex=2)
  segments(3.8+0.1*i, bias_wrong[i]-sd_wrong[i], 3.8+0.1*i, bias_wrong[i]+sd_wrong[i], lwd = 2 )
}
abline(h=0)

legend(0.5, -1.3, c('Doubly robust', 'Regression', 'IPTW'),  horiz = F, bty = 'n', cex=1.5,pch=c(15,16,17))

################ change sample size, additional simulation in the appendix
SIM=function(n,n_m){
  true3=GetTrue(n,n_m,3)
  est=NULL
  for (i in 1:nsim) {
    temp=EST(n, n_m, 3)
    est=rbind(est,temp)
  }
  bias_mu_q_eta=colMeans(est[,c(1:3)])-true3[1] # bias of DR, G_comp, IPTW
  sd_mu_q_eta=apply(est[,c(1:3)],2,sd) # Monte-Carlo sd of DR, G_comp, IPTW
  
  bias_mu=colMeans(est[,c(4:6)])-true3[1] # bias of DR, G_comp, IPTW
  sd_mu=apply(est[,c(4:6)],2,sd) # Monte-Carlo sd of DR, G_comp, IPTW
  
  bias_q_eta=colMeans(est[,c(7:9)])-true3[1] # bias of DR, G_comp, IPTW
  sd_q_eta=apply(est[,c(7:9)],2,sd) # Monte-Carlo sd of DR, G_comp, IPTW
  
  bias_wrong=colMeans(est[,c(10:12)])-true3[1] # bias of DR, G_comp, IPTW
  sd_wrong=apply(est[,c(10:12)],2,sd) # Monte-Carlo sd of DR, G_comp, IPTW
  
  return(round(rbind(cbind(bias_mu_q_eta,sd_mu_q_eta),cbind(bias_mu,sd_mu),cbind(bias_q_eta,sd_q_eta),cbind(bias_wrong,sd_wrong)),2))
}

results=NULL
for (n in c(10000,100000)) {
  for (n_m in c(1000,2000,5000)) {
    result=SIM(n,n_m)
    results=cbind(result,results)
  }
}
results
