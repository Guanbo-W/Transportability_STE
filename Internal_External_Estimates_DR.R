###### To show Biases, coverage of simultaneous confidence bands, coverage of pointwise confidence intervals,theoretical standard deviation and Monte-Carlo standard deviation of for each of the subgroup effects using the propsed doubly robust estimator. 
qtmax <- function(p, B, alpha){
  tmaxs <- apply(abs(matrix(rnorm(p*B), nrow = p, ncol = B)), 2, max)
  return(quantile(tmaxs, 1-alpha))
}
EST=function(n, n_m, x_tilde){
  Data=Gen(n, n_m)
  #Data$X1_1=ifelse(I(Data$X1==1),1,0);Data$X1_2=ifelse(I(Data$X1==2),1,0);Data$X1_3=ifelse(I(Data$X1==3),1,0);Data$X1_4=ifelse(I(Data$X1==4),1,0);Data$X1_5=ifelse(I(Data$X1==5),1,0)
  #X=Data[,1:10];R=Data$R;S=Data$S;A=Data$A;Y=Data$Y;S_R=Data$S[R==1];A_R=Data$A[R==1];
  Y_R=Data$Y[Data$R==1] 
  
  #Data[R==1,]=Data[R==1,]
  R=Data$R
  n_internal=dim(Data[R==1,])[1]
  #length(Y_R)
  #n_internal
  ################################# for psi
  ##### nuisance parameters for psi_{1,1}(x_tilde):E(Y1|X1=x_tilde,S=1)
  # inverse_weight={1/n_internal sum_{i=1}^{n_internal}I(X1=x_tilde, S=1)}^(-1)
  inverse_weight=(length(which(Data[R==1,]$X1==x_tilde& Data[R==1,]$S==1))/n_internal)^(-1)
  # mu=E(Y|A=1,X, R=1)
  model_mu=gam(Y~A+as.factor(X1)+s(X2)+s(X3)+s(X4)+s(X5)+s(X6)+s(X7)+lo(X8)+lo(X9)+lo(X10)+A*as.factor(X1)+A*X2+A*X3+A*X4+A*X5,data=Data[R==1,])  
  newdata=data.frame(rep(1,n_internal),as.factor(Data[R==1,]$X1),Data[R==1,]$X2,Data[R==1,]$X3,Data[R==1,]$X4,Data[R==1,]$X5,Data[R==1,]$X6,Data[R==1,]$X7,Data[R==1,]$X8,Data[R==1,]$X9,Data[R==1,]$X10,as.factor(Data[R==1,]$X1),Data[R==1,]$X2,Data[R==1,]$X3,Data[R==1,]$X4,Data[R==1,]$X5)
  names(newdata)=c("A","X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","A1X1","A1X2","A1X3","A1X4","A1X5") 
  mu=predict.Gam(model_mu,newdata = newdata)
  
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
  ##### estimation
  second_piece=q/eta*(Y_R-mu)
  est_psi=inverse_weight/n_internal*(sum(mu[which(Data[R==1,]$X1==x_tilde& Data[R==1,]$S==1)])+sum(second_piece[which(Data[R==1,]$X1==x_tilde& Data[R==1,]$A==1)]))
  ##### pointwise SD and 95% confidence interval
  IF=inverse_weight*((mu-est_psi)*ifelse(I(Data[R==1,]$X1==x_tilde& Data[R==1,]$S==1),1,0)+second_piece*ifelse(I(Data[R==1,]$X1==x_tilde& Data[R==1,]$A==1),1,0))
  #ATE_uncentered_IF=(n_internal/length(which(Data[R==1,]$S==1)))*((mu)*ifelse(I( Data[R==1,]$S==1),1,0)+second_piece*ifelse(I(Data[R==1,]$A==1),1,0))
  sd_psi=sqrt(var(IF)/(n_internal))
  ci_psi_point=c(est_psi+qnorm(0.05/2)*sd_psi,est_psi+qnorm(1-0.05/2)*sd_psi)
  ##### simultaneous confidence bands
  quan_cb=qtmax(5,100000,0.05)
  cb_psi=c(est_psi-quan_cb*sd_psi,est_psi+quan_cb*sd_psi)
  ##### Results
  results1=c(est_psi,sd_psi,ci_psi_point, cb_psi)
  ################################# for phi
  ##### nuisance parameters for phi_{1,1}(x_tilde):E(Y1|X1=x_tilde,S=1)
  # inverse_weight={1/n_internal sum_{i=1}^{n_internal}I(X1=x_tilde, S=1)}^(-1)
  inverse_weight=(length(which(Data$X1==x_tilde& Data$R==0))/n)^(-1)
  # g=mu
  newdata=data.frame(rep(1,n),as.factor(Data$X1),Data$X2,Data$X3,Data$X4,Data$X5,Data$X6,Data$X7,Data$X8,Data$X9,Data$X10,as.factor(Data$X1),Data$X2,Data$X3,Data$X4,Data$X5)
  names(newdata)=c("A","X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","A1X1","A1X2","A1X3","A1X4","A1X5") 
  g=predict.Gam(model_mu,newdata = newdata)
  # p=Pr(R=1|X)
  model_p=glm(R~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10,data=Data,family = "binomial")
  p=predict.glm(model_p,newdata=Data[R==1,],type="response")
  # e=Pr(A=a|X, R=1)
  e=eta
  ##### estimation
  second_piece=((1-p)/(e*p))*(Y_R-mu)
  est_phi=inverse_weight/n*(sum(g[which(Data$X1==x_tilde& Data$R==0)])+sum(second_piece[which(Data[R==1,]$X1==x_tilde& Data[R==1,]$A==1)]))
  ##### pointwise SD and 95% confidence interval
  u=rep(0,n)
  u[which(Data$R==1)]=second_piece*ifelse(I(Data[R==1,]$X1==x_tilde& Data[R==1,]$A==1),1,0)
  IF=inverse_weight*((g-est_phi)*ifelse(I(Data$X1==x_tilde& Data$R==0),1,0)+u)
  sd_phi=sqrt(var(IF)/(n))
  ci_phi_point=c(est_phi+qnorm(0.05/2)*sd_phi,est_phi+qnorm(1-0.05/2)*sd_phi)
  ##### simultaneous confidence bands
  quan_cb=qtmax(5,100000,0.05)
  cb_phi=c(est_phi-quan_cb*sd_phi,est_phi+quan_cb*sd_phi)
  ##### Results
  results2=c(est_phi,sd_phi,ci_phi_point, cb_phi)
  results=c(results1,results2)
  return(results)
}

source("Data_Generation.R")
true1=GetTrue(n,n_m,1)
true2=GetTrue(n,n_m,2)
true3=GetTrue(n,n_m,3)
true4=GetTrue(n,n_m,4)
true5=GetTrue(n,n_m,5)

est=NULL
for (i in 1:nsim) {
  temp=EST(n, n_m, 1)
  est=rbind(est,temp)
}
coverage_CI_1_psi=mean(ifelse(est[,3]<true1[1] & true1[1]<est[,4],1,0))
coverage_CB_1_psi=mean(ifelse(est[,5]<true1[1] & true1[1]<est[,6],1,0))
coverage_CI_1_phi=mean(ifelse(est[,9]<true1[2] & true1[2]<est[,10],1,0))
coverage_CB_1_phi=mean(ifelse(est[,11]<true1[2] & true1[2]<est[,12],1,0))
est1=colMeans(est) # est, sd, CI
sd1=sd(est[,1]) # Monte-Carlo sd for psi
sd1_phi=sd(est[,7])

est=NULL
for (i in 1:nsim) {
  temp=EST(n, n_m, 2)
  est=rbind(est,temp)
}
coverage_CI_2_psi=mean(ifelse(est[,3]<true2[1] & true2[1]<est[,4],1,0))
coverage_CB_2_psi=mean(ifelse(est[,5]<true2[1] & true2[1]<est[,6],1,0))
coverage_CI_2_phi=mean(ifelse(est[,9]<true2[2] & true2[2]<est[,10],1,0))
coverage_CB_2_phi=mean(ifelse(est[,11]<true2[2] & true2[2]<est[,12],1,0))
est2=colMeans(est) # est, sd, CI
sd2=sd(est[,1])
sd2_phi=sd(est[,7])

est=NULL
for (i in 1:nsim) {
  temp=EST(n, n_m, 3)
  est=rbind(est,temp)
}
coverage_CI_3_psi=mean(ifelse(est[,3]<true3[1] & true3[1]<est[,4],1,0))
coverage_CB_3_psi=mean(ifelse(est[,5]<true3[1] & true3[1]<est[,6],1,0))
coverage_CI_3_phi=mean(ifelse(est[,9]<true3[2] & true3[2]<est[,10],1,0))
coverage_CB_3_phi=mean(ifelse(est[,11]<true3[2] & true3[2]<est[,12],1,0))
est3=colMeans(est) # est, sd, CI
sd3=sd(est[,1])
sd3_phi=sd(est[,7])

est=NULL
for (i in 1:nsim) {
  temp=EST(n, n_m, 4)
  est=rbind(est,temp)
}
coverage_CI_4_psi=mean(ifelse(est[,3]<true4[1] & true4[1]<est[,4],1,0))
coverage_CB_4_psi=mean(ifelse(est[,5]<true4[1] & true4[1]<est[,6],1,0))
coverage_CI_4_phi=mean(ifelse(est[,9]<true4[2] & true4[2]<est[,10],1,0))
coverage_CB_4_phi=mean(ifelse(est[,11]<true4[2] & true4[2]<est[,12],1,0))
est4=colMeans(est) # est, sd, CI
sd4=sd(est[,1])
sd4_phi=sd(est[,7])

est=NULL
for (i in 1:nsim) {
  temp=EST(n, n_m, 5)
  est=rbind(est,temp)
}
coverage_CI_5_psi=mean(ifelse(est[,3]<true5[1] & true5[1]<est[,4],1,0))
coverage_CB_5_psi=mean(ifelse(est[,5]<true5[1] & true5[1]<est[,6],1,0))
coverage_CI_5_phi=mean(ifelse(est[,9]<true5[2] & true5[2]<est[,10],1,0))
coverage_CB_5_phi=mean(ifelse(est[,11]<true5[2] & true5[2]<est[,12],1,0))
est5=colMeans(est) # est, sd, CI
sd5=sd(est[,1])
sd5_phi=sd(est[,7])

bias1=est1[1]-true1[1];bias2=est2[1]-true2[1];bias3=est3[1]-true3[1];bias4=est4[1]-true4[1];bias5=est5[1]-true5[1]
bias=c(bias1,bias2,bias3,bias4,bias5)
tes.lci=c(est1[3]-est1[1],est2[3]-est2[1],est3[3]-est3[1],est4[3]-est4[1],est5[3]-est5[1])
tes.uci=c(est1[4]-est1[1],est2[4]-est2[1],est3[4]-est3[1],est4[4]-est4[1],est5[4]-est5[1])
tes.lcb=c(est1[5]-est1[1],est2[5]-est2[1],est3[5]-est3[1],est4[5]-est4[1],est5[5]-est5[1])
tes.ucb=c(est1[6]-est1[1],est2[6]-est2[1],est3[6]-est3[1],est4[6]-est4[1],est5[6]-est5[1])
sd=c(est1[2],est2[2],est3[2],est4[2],est5[2])
MCMC_sd=c(sd1,sd2,sd3,sd4,sd5)
coverage_CI_psi=c(coverage_CI_1_psi,coverage_CI_2_psi,coverage_CI_3_psi,coverage_CI_4_psi,coverage_CI_5_psi)
coverage_CB_psi=c(coverage_CB_1_psi,coverage_CB_2_psi,coverage_CB_3_psi,coverage_CB_4_psi,coverage_CB_5_psi)
#### plots
par(mar = c(5, 4, 4, 4) + 0.3)
plot( c(1,5), las = 2, xlim =c(0.6, 5.4), ylim = c(-2, 2),  type="n",xlab=TeX(r'(Subgroup)', bold=TRUE), 
      ylab=TeX(r'(Bias\pm$SD)', bold=TRUE), xaxt="n")
axis(1, at=1:5, labels=c("1","2","3","4","5"))
for (i in 1:5){
  # rect(i-0.2, tes.lci[i], i+0.2,  tes.uci[i], col = NA,  border = "red", lwd = 3)  
  # rect(i-0.2, tes.lcb[i], i+0.2, tes.ucb[i], col = NA,  border = 4, lwd = 3 )
  segments(i-0.2, bias[i], i+0.2, bias[i], lwd = 5 )
  segments(i-0.05, bias[i]-sd[i], i-0.05, bias[i]+sd[i], lwd = 3 ,col="yellow")
  segments(i+0.05, bias[i]-MCMC_sd[i], i+0.05, bias[i]+MCMC_sd[i], lwd = 3 ,col="green")
}

abline(h=0)
par(new = TRUE)
plot(c(1:5), coverage_CI_psi, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "",ylim=c(0.7,1),col=2,lwd = 3)
lines(c(1:5),coverage_CB_psi,type = "l", lty =1,col=4,lwd = 3)
axis(side=4, ylim=c(0.7,1),at = c(0.9,0.95,1))
mtext(TeX(r'(Coverage)', bold=TRUE), side=4, line=3)

legend(2.5, 0.75, c('Bias', 'Coverage of simultaneous confidence band', 'Coverage of pointwise confidence interval',"Theoretical standard deviation", "Monte-Carlo standard deviation"), col = c(1,4,2,"yellow","green"), lwd = c(5,3,3,3,3), horiz = F, bty = 'n', cex=1)

bias1=est1[7]-true1[2];bias2=est2[7]-true2[2];bias3=est3[7]-true3[2];bias4=est4[7]-true4[2];bias5=est5[7]-true5[2]
bias=c(bias1,bias2,bias3,bias4,bias5)
tes.lci=c(est1[9]-est1[7],est2[9]-est2[7],est3[9]-est3[7],est4[9]-est4[7],est5[9]-est5[7])
tes.uci=c(est1[10]-est1[7],est2[10]-est2[7],est3[10]-est3[7],est4[10]-est4[7],est5[10]-est5[7])
tes.lcb=c(est1[11]-est1[7],est2[11]-est2[7],est3[11]-est3[7],est4[11]-est4[7],est5[11]-est5[7])
tes.ucb=c(est1[12]-est1[7],est2[12]-est2[7],est3[12]-est3[7],est4[12]-est4[7],est5[12]-est5[7])
sd=c(est1[8],est2[8],est3[8],est4[8],est5[8])
MCMC_sd=c(sd1_phi,sd2_phi,sd3_phi,sd4_phi,sd5_phi)
coverage_CI_phi=c(coverage_CI_1_phi,coverage_CI_2_phi,coverage_CI_3_phi,coverage_CI_4_phi,coverage_CI_5_phi)
coverage_CB_phi=c(coverage_CB_1_phi,coverage_CB_2_phi,coverage_CB_3_phi,coverage_CB_4_phi,coverage_CB_5_phi)
#### plots
par(mar = c(5, 4, 4, 4) + 0.3)
plot( c(1,5), las = 2, xlim =c(0.6, 5.4), ylim = c(-2, 2),  type="n",xlab=TeX(r'(Subgroup)', bold=TRUE), 
      ylab=TeX(r'(Bias$\pm$SD)', bold=TRUE), xaxt="n")
axis(1, at=1:5, labels=c("1","2","3","4","5"))
for (i in 1:5){
  # rect(i-0.2, tes.lci[i], i+0.2,  tes.uci[i], col = NA,  border = "red", lwd = 3)  
  # rect(i-0.2, tes.lcb[i], i+0.2, tes.ucb[i], col = NA,  border = 4, lwd = 3 )
  segments(i-0.2, bias[i], i+0.2, bias[i], lwd = 5 )
  segments(i-0.05, bias[i]-sd[i], i-0.05, bias[i]+sd[i], lwd = 3 ,col="yellow")
  segments(i+0.05, bias[i]-MCMC_sd[i], i+0.05, bias[i]+MCMC_sd[i], lwd = 3 ,col="green")
}
abline(h=0)
par(new = TRUE)
plot(c(1:5), coverage_CI_phi, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "",ylim=c(0.7,1),col=2,lwd = 3)
lines(c(1:5),coverage_CB_phi,type = "l", lty =1,col=4,lwd = 3)
axis(side=4, ylim=c(0.7,1),at = c(0.9,0.95,1))
mtext(TeX(r'(Coverage)', bold=TRUE), side=4, line=3)

legend(2.5, 0.75, c('Bias', 'Coverage of simultaneous confidence band', 'Coverage of pointwise confidence interval',"Theoretical standard deviation", "Monte-Carlo standard deviation"), col = c(1,4,2,"yellow","green"), lwd = c(5,3,3,3,3), horiz = F, bty = 'n', cex=1)

