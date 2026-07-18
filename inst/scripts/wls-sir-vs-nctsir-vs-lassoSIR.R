###boucle -------
rm(list = ls())
setwd("~/Desktop/phd /R codes /VC-pcr codes ")
source("~/Desktop/phd /R codes /myfunctions.R")

#source("~/Desktop/phd /R codes /VC-pcr codes /all VCPCR function.R")
#q=20;p=200; n=100 #;n=50
#Beta= c(rep(1,5),rep(0,5),rep(-1,5),rep(0,5),rep(1,5),rep(0,5),rep(-1,5),rep(0,5),rep(0,160))
#S<-which(Beta!=0)    

#q=floor(sqrt(p)) 
#h = max(2, round(nrow(X)/5)
n=500; p=700;K=50;H=10;q=20; alpha=c(0.2, 0.3, 0.7) ;rho = 0.5 #(0.3, 0.4, 0.8)----
#n=300; p=500; q= 10;H=5 ;alpha=c (0.3,0.5,0.7)  ; K= 50;rho = 0.5 #beta= 0.1 
#n=1200; p=1500; q=70;H=5  ; k= 50;alpha=c(0.2, 0.4,0.9);   
#n=1200; p=1500; q=70;H=5  ; k= 50;alpha=c(0.5, 0.6, 0.9); beta= 0.1 
#n=1200; p=1500; q=40;H=5 ; alpha=c(0.3,0.4,0.8) ; k= 50;rho = 0.5  #beta= 0.1 
#Beta= c(rep(0.1,q),rep(0,p-q));tho=0.5;-----
Beta= c(rep(1,q),rep(0,p-q)); tau=0.1
#values <- c(-2, -1, 1, 2)
#random_values <- sample(values, q, replace = TRUE)
#Beta<-c(random_values ,rep(0,p-q)); tau=0.1  
#r = c("lassoSIR"=0,"plsir"=0,"palsir"=0, "pralsir"=0,"psir"=0, "pcsir"=0,"nctsir"=0)

r = c("lassoSIR"=0,"nctsir"=0,"wls.sir"=0, "LASSO"=0, "SirCHOMP"=0)
result = list("abscos"=r,"Cor"=r ,"Time"=r, "TPR"=r, "FPR"=r ,"FNR"=r)
repetition = 1
abscos=matrix(nrow=5,ncol=repetition)
TPR=matrix(nrow=5,ncol=repetition)
FPR=matrix(nrow=5,ncol=repetition)
FNR=matrix(nrow=5,ncol=repetition)
Cor=matrix(nrow=5,ncol=repetition)
Time=matrix(nrow=5,ncol=repetition)

#### estimation with normalisation and without round------
iteration = 1
set.seed(1)

#set.seed(1234)

while (iteration <= repetition) {
  # X=X_rand(n=n,p=p,q=q,type = "homog",rho)
# X=X_rand(n=n,p=p,q=q,type = "autoR",rho)
X=X_rand(n=n,p=p,q=q,type = "block",alpha=alpha)
  #X<- mvtnorm::rmvnorm(n, rep(0,p), diag(p))
  t=X%*%Beta
  epsilon=rnorm(n,mean=0,sd=1)
  # y=as.vector(1+exp(t/sqrt(3))+epsilon )  ## expo
  sig=1
 # y= as.vector((X[,s[1]]+X[,s[2]]+1,5*X[,s[3]]+1.2*X[,s[4]]) /(0.5+(X[,s[5]]+1.2*X[,s[6]]+1)^2) + sig*epsilon)
 # X=mvrnorm(n, rep(0,p), diag(1,p))
#  y=X[,1] + X[,2] + X[,3] + X[,4] + X[,5] + X[,6] + rnorm(n)
  
#   y=as.vector(exp(t)+epsilon )  ## expo
# y=as.vector(t+epsilon) # linear cas 
  y=as.vector(t+ sig*epsilon  )            ### model 1 
 # y= as.vector(((t)^3)/2 + epsilon )    ### model 2
  #y=as.vector(exp(t)+epsilon)         ### model 3 
  #
 # y=as.vector(sin(pi*t/2)+epsilon )### model 4
  
      #y=as.vector(sin(t)*exp(t)+epsilon )### model 4
  # y=t/sqrt(3)+ epsilon
  
  ####### LASSOSIR-
  T1<-Sys.time() 
  sir.lasso <- LassoSIR( X, y, H, no.dim=1,solution.path=FALSE, categorical=FALSE, nfolds=10,screening=FALSE)
  beta.lassosir <- as.vector(sir.lasso$beta)
 #beta.lassosir<-normalize_vector(beta.lassosir)
  T2<-Sys.time()
  Time[1,iteration]<-round(difftime(T2, T1, units = "secs"),3)
  abscos[1,iteration]=round(abs(costeta(Beta,beta.lassosir)),3)
  #Dist[1,iteration]=round(Dproj(Beta,beta.lassosir),3)
  if(var(beta.lassosir)!=0) {Cor[1,iteration]=round(abs(cor(X%*%beta.lassosir,t)),3)} else {Cor[1,iteration]=0 };
  TPR[1,iteration]=round(TPR.funct(Beta,beta.lassosir),3)
  FPR[1,iteration]=round(FPR.funct(Beta,beta.lassosir),3)
  FNR[1,iteration]=round(FNR.funct(Beta,beta.lassosir),3)
  rm(T1,T2)
  
  
  ###SIR CHOMP --
  T1<-Sys.time() 
  beta.chomp<-as.numeric(AdaptCHOMPwithPIC(X, y, d = 1, gamma.pow = 2))

    T2<-Sys.time()
  Time[2,iteration]<-round(difftime(T2, T1, units = "secs"),3)
  abscos[2,iteration]=round(abs(costeta(Beta,beta.chomp)),3)
  #Dist[2,iteration]=round(Dproj(Beta,beta.svd),3)
  if(var(beta.chomp)!=0) {Cor[2,iteration]=round(abs(cor(X%*%beta.chomp,t)),3)} else {Cor[2,iteration]=0 };
  FPR[2,iteration]=round(FPR.funct(Beta,beta.chomp),3)
  FNR[2,iteration]=round(FNR.funct(Beta,beta.chomp),3)
  TPR[2,iteration]=round(TPR.funct(Beta,beta.chomp),3)
  rm(T1,T2)
  
  ###SIR+NCT natural cannonical tresholding direct --
  T1<-Sys.time() 
  # tau_min<-nct_cv_opt(x,y,k,tau_num =10, threshold="soft")
  #beta.nct<-as.vector(sir.nct(X,y,K,tau, threshold="soft"))
  beta.nct<-normalize_vector(as.vector(sir.nct(X,y,K,tau, threshold="soft")))
  #beta.nct <- beta.nct/sqrt( sum( beta.nct^2 ))#normaliser
  T2<-Sys.time()
  Time[3,iteration]<-round(difftime(T2, T1, units = "secs"),3)
  abscos[3,iteration]=round(abs(costeta(Beta,beta.nct)),3)
  #Dist[3,iteration]=round(Dproj(Beta,beta.nct),3)
  if(var(beta.nct)!=0) {Cor[3,iteration]=round(abs(cor(X%*%beta.nct,t)),3)} else {Cor[3,iteration]=0 };
  TPR[3,iteration]=round(TPR.funct(Beta,beta.nct),3)
  FPR[3,iteration]=round(FPR.funct(Beta, beta.nct),3)
  FNR[3,iteration]=round(FNR.funct(Beta, beta.nct),3)
   rm(T1,T2)
  
  ###SIR+WLS  -
  T1<-Sys.time() 
  wlsir<-wls.sir(X,as.vector(y))
  T2<-Sys.time()
  Time[4,iteration]<-round(difftime(T2, T1, units = "secs"),3)
  
  beta.wls<-as.vector(Re(wlsir$betahat))
  Time[4,iteration]<-round(difftime(T2, T1, units = "secs"),3)
  
  abscos[4,iteration]=round(abs(costeta(Beta, beta.wls)),3)
 # Dist[4,iteration]=round(Dproj(Beta, beta.wls),4)
  if(var(beta.nct)!=0) {Cor[4,iteration]=round(abs(cor(X%*%beta.wls,t)),3)} else {Cor[4,iteration]=0 };
  TPR[4,iteration]=round(TPR.funct(Beta,beta.wls),3)
  FPR[4,iteration]=round(FPR.funct(Beta,beta.wls),3)
  FNR[4,iteration]=round(FNR.funct(Beta,beta.wls),3)

  ###lasso  -
  T1<-Sys.time() 
  #beta.lasso<-as.vector(Predict(X,y,method="LASSO"))
 beta.lasso<-normalize_vector(as.vector(Predict(X,y,method="LASSO")))
  T2<-Sys.time()
  Time[5,iteration]<-round(difftime(T2, T1, units = "secs"),3)
  abscos[5,iteration]=round(abs(costeta(Beta, beta.lasso)),3)
  #Dist[5,iteration]=round(Dproj(Beta, beta.lasso),5)
  if(var(beta.lasso)!=0) {Cor[5,iteration]=round(abs(cor(X%*% beta.lasso,t)),3)} else {Cor[5,iteration]=0 };
  TPR[5,iteration]=round(TPR.funct(Beta,beta.lasso),3)
  FPR[5,iteration]=round(FPR.funct(Beta, beta.lasso),3)
  FNR[5,iteration]=round(FNR.funct(Beta, beta.lasso),3)
  
  iteration = iteration + 1
}
result$abscos["lassoSIR"]=sum(abscos[1,]) ;result$Cor["lassoSIR"]=sum(Cor[1,]) ;#result$Dist["lassoSIR"]=sum(Dist[1,])
result$Time["lassoSIR"]=sum(Time[1,]);result$FNR["lassoSIR"]=sum(FNR[1,]) 
result$TPR["lassoSIR"]=sum(TPR[1,]);result$FPR["lassoSIR"]=sum(FPR[1,]) 

result$abscos["SirCHOMP"]=sum(abscos[2,]) ;result$Cor["SirCHOMP"]=sum(Cor[2,]) #result$Dist["SirCHOMP"]=sum(Dist[2,]);
result$Time["SirCHOMP"]=sum(Time[2,]); result$FNR["SirCHOMP"]=sum(FNR[2,]) 
result$TPR["SirCHOMP"]=sum(TPR[2,]);result$FPR["SirCHOMP"]=sum(FPR[2,]) 


result$abscos["nctsir"]=sum(abscos[3,]) ;result$Cor["nctsir"]=sum(Cor[3,]) #result$Dist["nctsir"]=sum(Dist[3,]);
result$Time["nctsir"]=sum(Time[3,]);result$FNR["nctsir"]=sum(FNR[3,]) 
result$TPR["nctsir"]=sum(TPR[3,]);result$FPR["nctsir"]=sum(FPR[3,]) 


result$abscos["wls.sir"]=sum(abscos[4,]) ;result$Cor["wls.sir"]=sum(Cor[4,]) #result$Dist["wls.sir"]=sum(Dist[4,]);
result$Time["wls.sir"]=sum(Time[4,]) ;result$FNR["wls.sir"]=sum(FNR[4,]) 
result$TPR["wls.sir"]=sum(TPR[4,]);result$FPR["wls.sir"]=sum(FPR[4,]) 

result$abscos["LASSO"]=sum(abscos[5,]) ;result$Cor["LASSO"]=sum(Cor[5,]) #result$Dist["LASSO"]=sum(Dist[5,]);
result$Time["LASSO"]=sum(Time[5,]);result$FNR["LASSO"]=sum(FNR[5,]) 
result$TPR["LASSO"]=sum(TPR[5,]);result$FPR["LASSO"]=sum(FPR[5,]) 


res=lapply(result,  function(x) x/repetition)
res

#ratio[1,2]=round((res$FPR["lassoSIR"]+res$FNR["lassoSIR"])/(res$FPR["DDPCA.spac.sir"]+res$FNR["DDPCA.spac.sir"]),3)

####------#les plot ------
t1=X%*%beta.lassosir
t4=X%*%beta.wls
t5=X%*%beta.lasso
#t=x%*%b1+(x%*%b2)^2
#t1=x%*%beta[,1]+(x%*%beta[,2])^2
#t2=x%*%dr$evectors[,1]+(x%*%dr$evectors[,2])^2
y1= as.vector(((t1)^3)/2)    ### model 2
y4= as.vector(((t4)^3)/2)    ### model 2
y5= as.vector(((t5)^3)/2)    ### model 2
#par(mgp=c(2.5, 1, 0))
par(mfrow = c(1, 3))

# Set up the plot with the first series
plot(t, y, col="red") #main="Comparison of Series", pch=16, cex=0.5, xlab="X-Axis Label", ylab="Y-Axis Label")
# Add other series to the plot using points() function
points(t1, y1, col="blue")##. ### t1 lasso sir 

plot(t, y, col="red") #main="Comparison of Series", pch=16, cex=0.5, xlab="X-Axis Label", ylab="Y-Axis Label")
points(t4, y4, col="green")#. #### t4 wls SIR 

plot(t, y, col="red") #main="Comparison of Series", pch=16, cex=0.5, xlab="X-Axis Label", ylab="Y-Axis Label")
points(t5, y5, col="orange")#.  #### t5 lasso 

# Add legend
legend("topright", legend=c("t1", "t4", "t5"), col=c("blue", "green", "orange"), pch=16)



# Set up the plot with the first series
plot(sort(t), sort(y), col="red",type="l") #main="Comparison of Series", pch=16, cex=0.5, xlab="X-Axis Label", ylab="Y-Axis Label")
#plot(t,y,pch=".")
# Add other series to the plot using points() function
lines(sort(t1), sort(y1), col="blue")
lines(sort(t4),sort( y4), col="green")
lines(sort(t5),sort(y5), col="orange")
#lines(sort(t), sort(y))
####
par(mgp=c(2.5, 1, 0))
plot(sort(t),sort(t), type='l', xlab=expression(X~beta), ylab=expression(X~hat(beta)),las=2)
#lines(sort(t), sort(t))
lines(sort(t), sort(t1), col="blue")
lines(sort(t),sort( t4), col="green")
lines(sort(t),sort(t5), col="orange")# Add legend

legend("topright", legend=c("real","lassoSIR", "wls.SIR", "LASSO"), col=c("black","blue", "green", "orange"), pch=16, x.intersp = 0.5, y.intersp = 0.5, inset = c(0.5, 0.1))
######
par(mgp=c(2.5, 1, 0))
par(mfrow = c(2, 2))
#barplot(val,ylim=c(0,0.7),col="blue",xlab="",cex.axis = 0.7,main = "valeurs propres"); box()
#plot(t,t1,main = "t1 en fonction de t")
plot(t,y,col="red",main="real y ", pch=16, cex =0.5, xlab=expression(X~beta),las=2)
plot(t1,y,col="blue",main="lassosir   ",pch=16, cex =0.5, xlab=expression(X~hat(beta)),las=2)
#plot(t2,y,col="brown",main="nctSIR",pch=16, cex =0.5
#plot(t3,y,main = " sirSVD",pch=16, cex =0.5)
plot(t4,y,col="green",main="wls.SIR  ",pch=16, cex =0.5,xlab=expression(X~hat(beta)),las=2)
plot(t5,y,col="orange",main="lasso  ",pch=16, cex =0.5,xlab=expression(X~hat(beta)),las=2)


#### with ggplots -----
par(mfrow = c(2, 2))


ggplot() +
  geom_point(aes(x = t, y = y), color = "red", shape = 16, size = 2) +
  labs(title = "Real y", x = expression(X~beta), y = "y") +
  theme_minimal()

# Create a ggplot for the second plot
ggplot() +
  geom_point(aes(x = t1, y = y), color = "blue", shape = 16, size = 2) +
  labs(title = "lassosir", x = expression(X~hat(beta)), y = "y") +
  theme_minimal()

# Create a ggplot for the third plot
ggplot() +
  geom_point(aes(x = -t4, y = y), color = "green", shape = 16, size = 2) +
  labs(title = "wls.SIR", x = expression(X~hat(beta)), y = "y") +
  theme_minimal()

# Create a ggplot for the fourth plot
ggplot() +
  geom_point(aes(x = t5, y = y), color = "orange", shape = 16, size = 2) +
  labs(title = "lasso", x = expression(X~hat(beta)), y =  "y") +
  theme_minimal()

# Add legend to the last plot
last_plot() +
  theme(legend.position = "topright") +
  scale_color_manual(
    values = c("red", "blue", "green", "orange"),
    name = "Legend",
    labels = c("real", "lassoSIR", "wls.SIR", "LASSO")
  )




par(mgp=c(2.5, 1, 0))
par(mfrow = c(1, 1))
plot(sort(t),sort(t), type='l',asp=1,  xlab=expression(X~beta), ylab=expression(X~hat(beta)),las=2)
points(t,t1,col="blue",main="lassosir   ",pch=16, cex =0.5, xlab=expression(X~beta), ylab=expression(X~hat(beta)),las=2)
#plot(t2,y,col="brown",main="nctSIR",pch=16, cex =0.5
#plot(t3,y,main = " sirSVD",pch=16, cex =0.5)
plot(t,t, type='l', xlab=expression(X~beta), ylab=expression(X~hat(beta)),las=2)
points(t,t4,col="green",main="wls.SIR  ",pch=16, cex =0.5,xlab=expression(X~beta), ylab=expression(X~hat(beta)),las=2)

plot(sort(t),sort(t), type='l', xlab=expression(X~beta), ylab=expression(X~hat(beta)),las=2)
points(t,t5,col="orange",main="lasso  ",pch=16, cex =0.5,xlab=expression(X~beta), ylab=expression(X~hat(beta)),las=2)



######### in same ggplot ----
plot_real <-ggplot() +
  geom_point(aes(x = t, y = y), color = "red", shape = 16, size = 2) +
  labs(title = " real direction", x = expression(X~beta), y = "y") +
  theme_minimal()

plot_lassosir <-ggplot() +
  geom_point(aes(x = t1, y = y), color = "blue", shape = 16, size = 2) +
  labs(title = "lassosir", y= "y", x= expression(X~hat(beta))) +
  theme_minimal()


# Create a ggplot for the second plot
plot_wlsSIR <-ggplot() +
  geom_point(aes(x = t4, y = y), color = "green", shape = 16, size = 2) +
  labs(title = "wls.SIR", x = expression(X~hat(beta)), y = "y") +
  theme_minimal()

# Create a ggplot for the third plot
plot_lasso <-ggplot() +
  geom_point(aes(x = t5, y = y), color = "orange", shape = 16, size = 2) +
  labs(title = "lasso", x = expression(X~hat(beta)) , y = "y") +
  theme_minimal()

grid.arrange(plot_real, plot_lassosir, plot_wlsSIR, plot_lasso,
             ncol = 2)


#T1<-Sys.time() ------ normal estimation -----
iteration = 1
set.seed(1)
while (iteration <= repetition) {
  # X=X_rand(n=n,p=p,q=q,type = "homog",rho)
  X=X_rand(n=n,p=p,q=q,type = "autoR",rho)
  #X=X_rand(n=n,p=p,q=q,type = "block",alpha=alpha)
  # x = matrix(rnorm(n*p), nrow = n, ncol = p)
  t=X%*%Beta
  epsilon=rnorm(n,mean=0,sd=1)
  # y=as.vector(1+exp(t/sqrt(3))+epsilon )  ## expo
  
  #y=as.vector(exp(t)+epsilon )  ## expo
  # y=as.vector(t+epsilon) # linear cas 
  y=as.vector(t+ epsilon  )            ### model 1 
  # y= as.vector(((t)^3)/2 + epsilon )    ### model 2
  #  y=as.vector(exp(t)+epsilon)         ### model 3 
  #
  #    y=as.vector(sin(t)*exp(t)+epsilon )### model 4
  # y=t/sqrt(3)+ epsilon
  
  ####### LASSOSIR-
  T1<-Sys.time() 
  sir.lasso <- LassoSIR( X, y, H, no.dim=1,solution.path=FALSE, categorical=FALSE, nfolds=10,screening=FALSE)
  beta.lassosir <- as.vector(sir.lasso$beta)
  T2<-Sys.time()
  Time[1,iteration]<-round(difftime(T2, T1, units = "secs"),3)
  abscos[1,iteration]=round(abs(costeta(Beta,beta.lassosir)),3)
  # Dist[1,iteration]=round(Dproj(Beta,beta.lassosir),3)
  if(var(beta.lassosir)!=0) {Cor[1,iteration]=round(abs(cor(X%*%beta.lassosir,t)),3)} else {Cor[1,iteration]=0 };
  TPR[1,iteration]=round(TPR.funct(Beta,round(beta.lassosir,4)),3)
  FPR[1,iteration]=round(FPR.funct(Beta,round(beta.lassosir,4)),3)
  FNR[1,iteration]=round(FNR.funct(Beta,round(beta.lassosir,4)),3)
  rm(T1,T2)
  
  
  ###SIR+PCA direct --
  T1<-Sys.time() 
  beta.svd<-as.vector(sir.svd(X,y,K))
  #beta.svd <- beta.svd/sqrt( sum( beta.svd^2 ))#normaliser
  T2<-Sys.time()
  Time[2,iteration]<-round(difftime(T2, T1, units = "secs"),3)
  abscos[2,iteration]=round(abs(costeta(Beta,beta.svd)),3)
  #Dist[2,iteration]=round(Dproj(Beta,beta.svd),3)
  if(var(beta.svd)!=0) {Cor[2,iteration]=round(abs(cor(X%*%beta.svd,t)),3)} else {Cor[2,iteration]=0 };
  FPR[2,iteration]=round(FPR.funct(Beta,round(beta.svd,4)),3)
  TPR[2,iteration]=round(TPR.funct(Beta,round(beta.svd,4)),3)
  FNR[2,iteration]=round(FNR.funct(Beta,round(beta.svd,4)),3)
  
  rm(T1,T2)
  
  ###SIR+NCT natural cannonical tresholding direct 
  T1<-Sys.time() 
  # tau_min<-nct_cv_opt(x,y,k,tau_num =10, threshold="soft")
  beta.nct<-as.vector(sir.nct(X,y,K,tau, threshold="soft"))
  #beta.nct <- beta.nct/sqrt( sum( beta.nct^2 ))#normaliser
  T2<-Sys.time()
  Time[3,iteration]<-round(difftime(T2, T1, units = "secs"),3)
  abscos[3,iteration]=round(abs(costeta(Beta,beta.nct)),3)
  #Dist[3,iteration]=round(Dproj(Beta,beta.nct),3)
  if(var(beta.nct)!=0) {Cor[3,iteration]=round(abs(cor(X%*%beta.nct,t)),3)} else {Cor[3,iteration]=0 };
  TPR[3,iteration]=round(TPR.funct(Beta, round(beta.nct,4)),3)
  FPR[3,iteration]=round(FPR.funct(Beta, round(beta.nct,4)),3)
  FNR[3,iteration]=round(FNR.funct(Beta, round(beta.nct,4)),3)
  
  
  ###SIR+WLS  
  T1<-Sys.time() 
  beta.wls<-as.vector(wls.sir(X,y))
  T2<-Sys.time()
  Time[4,iteration]<-round(difftime(T2, T1, units = "secs"),3)
  abscos[4,iteration]=round(abs(costeta(Beta, beta.wls)),3)
  #Dist[4,iteration]=round(Dproj(Beta, beta.wls),4)
  if(var(beta.nct)!=0) {Cor[4,iteration]=round(abs(cor(X%*% beta.wls,t)),3)} else {Cor[4,iteration]=0 };
  TPR[4,iteration]=round(TPR.funct(Beta, round(beta.wls,4)),3)
  FPR[4,iteration]=round(FPR.funct(Beta, round(beta.wls,4)),3)
  FNR[4,iteration]=round(FNR.funct(Beta, round(beta.wls,4)),3)
  
  
  
  
  ###lasso---
  T1<-Sys.time() 
  beta.lasso<-as.vector(Predict(X,y,method="LASSO"))
  T2<-Sys.time()
  Time[5,iteration]<-round(difftime(T2, T1, units = "secs"),3)
  abscos[5,iteration]=round(abs(costeta(Beta, beta.lasso)),3)
  #   Dist[5,iteration]=round(Dproj(Beta, beta.lasso),5)
  if(var(beta.nct)!=0) {Cor[5,iteration]=round(abs(cor(X%*% beta.lasso,t)),3)} else {Cor[5,iteration]=0 };
  TPR[5,iteration]=round(TPR.funct(Beta, round(beta.lasso,4)),3)
  FPR[5,iteration]=round(FPR.funct(Beta, round(beta.lasso,4)),3)
  FNR[5,iteration]=round(FNR.funct(Beta, round(beta.lasso,4)),3)
  
  
  iteration = iteration + 1
}

result$abscos["lassoSIR"]=sum(abscos[1,]) ;result$Cor["lassoSIR"]=sum(Cor[1,]) ;#result$Dist["lassoSIR"]=sum(Dist[1,])
result$Time["lassoSIR"]=sum(Time[1,]);result$FNR["lassoSIR"]=sum(FNR[1,]) 
result$TPR["lassoSIR"]=sum(TPR[1,]);result$FPR["lassoSIR"]=sum(FPR[1,]) 

result$abscos["sirSVD"]=sum(abscos[2,]) ;result$Cor["sirSVD"]=sum(Cor[2,]) #result$Dist["sirSVD"]=sum(Dist[2,]);
result$Time["sirSVD"]=sum(Time[2,]); result$FNR["sirSVD"]=sum(FNR[2,]) 
result$TPR["sirSVD"]=sum(TPR[2,]);result$FPR["sirSVD"]=sum(FPR[2,]) 


result$abscos["nctsir"]=sum(abscos[3,]) ;result$Cor["nctsir"]=sum(Cor[3,]) #result$Dist["nctsir"]=sum(Dist[3,]);
result$Time["nctsir"]=sum(Time[3,]);result$FNR["nctsir"]=sum(FNR[3,]) 
result$TPR["nctsir"]=sum(TPR[3,]);result$FPR["nctsir"]=sum(FPR[3,]) 


result$abscos["wls.sir"]=sum(abscos[4,]) ;result$Cor["wls.sir"]=sum(Cor[4,]) #result$Dist["wls.sir"]=sum(Dist[4,]);
result$Time["wls.sir"]=sum(Time[4,]) ;result$FNR["wls.sir"]=sum(FNR[4,]) 
result$TPR["wls.sir"]=sum(TPR[4,]);result$FPR["wls.sir"]=sum(FPR[4,]) 

result$abscos["LASSO"]=sum(abscos[5,]) ;result$Cor["LASSO"]=sum(Cor[5,]) #result$Dist["LASSO"]=sum(Dist[5,]);
result$Time["LASSO"]=sum(Time[5,]);result$FNR["LASSO"]=sum(FNR[5,]) 
result$TPR["LASSO"]=sum(TPR[5,]);result$FPR["LASSO"]=sum(FPR[5,]) 



res=lapply(result,  function(x) x/repetition)

res



