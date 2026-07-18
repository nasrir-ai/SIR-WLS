#### source------
rm(list = ls())
source("~/Desktop/phd /R codes /myfunctions.R")
setwd("~/Desktop/phd /R codes ")

# An example of fitting CHOMP and Adaptive CHOMP using PIC to select the tuning parameterlibrary(MASS)-------
n = 500; p = 100
# Generating X ~ N_p(0, Sigma)
set.seed(2000)
cov_sigma <- 0.5^abs(outer(1:p, 1:p, "-"))
omega <- diag(sqrt(runif(p, 0.5, 4)))
SIGMA <- omega %*% cov_sigma %*%  omega
X = mvrnorm(n, mu = rep(0,p), Sigma = SIGMA)
##### 
omega<-rho^abs(outer(1:p,1:p,"-"))
D <- diag(runif(p, 0.5, 2))
Sigma <- D %*% omega %*% D
X = mvrnorm(n, mu = rep(0,p), Sigma = SIGMA)
# Example 1: Single index model-------
s <- 5
beta <- rep(0, p)
beta[1:s] <- sample(c(-1,1), s, replace=TRUE)*runif(s, 1, 1.5)
y <- exp(X %*% beta) + rnorm(n)
beta.chomp<-as.numeric(fit.AdaptCHOMP_WeightPower1 <- AdaptCHOMPwithPIC(X, y, d = 1, gamma.pow = 1))
fit.AdaptCHOMP_WeightPower2 <- as.numeric(AdaptCHOMPwithPIC(X, y, d = 1, gamma.pow = 2))

select.chomp<-which(beta.chomp!=0)

sir.lasso <- LassoSIR( X,y,solution.path=FALSE, categorical=TRUE, nfolds=10,screening=TRUE,no.dim=1)
beta.lassosir <-as.vector(sir.lasso$beta)
select.lassosir<-which(beta.lassosir!=0)


#wls<-wls.sir( x=X.arcene_train_with_intercept,y=Y.arcene_train,categorical=TRUE )

wls<-wls.sir( X,y,categorical=TRUE )
beta.wls<-wls$betahat
#select.wls<-wls$select
beta.slect<-which(beta.wls!=0)

###
round(abs(costeta(beta,beta.lassosir)),3)
round(abs(costeta(beta,beta.wls)),3)
round(abs(costeta(beta,beta.chomp)),3)

### cor  round(abs(cor(X%*%beta.lassosir,X%*%beta)
round(abs(cor(X%*%beta.lassosir,X%*%beta)),3)
round(abs(cor(X%*%beta.wls,X%*%beta)),3)
round(abs(cor(X%*%beta.chomp,X%*%beta)),3)
### FPR , TPR 
round(TPR.funct(beta,beta.lassosir),3);round(FPR.funct(beta,beta.lassosir),3)
round(TPR.funct(beta,beta.wls),3);round(FPR.funct(beta,beta.wls),3)
round(TPR.funct(beta,beta.chomp),3);round(FPR.funct(beta,beta.chomp),3)



# Example 2: Double index model----------
# Generating sparse beta and outcome from a double index model-------
set.seed(200)
s <- list(s1 = 1:4, s2 = 3:6) # each element of s corresponds to the index of non-zero coefficient in each dimension
beta = matrix(0, p, length(s))
for (j in 1:length(s)){
  beta[s[[j]], j] = sample(c(-1,1), length(s[[j]]), replace=TRUE)*runif(length(s[[j]]), 1, 1.5)
}
Z = X %*% beta
y = (Z[, 1]) * exp(Z[,2] + rnorm(n))   

beta.chomp<- as.matrix(AdaptCHOMPwithPIC(X, y, d = 2))


select.chomp<-list(which(beta.chomp[,1]!=0),which(beta.chomp[,2]!=0))

sir.lasso <- LassoSIR( X,y,solution.path=FALSE, categorical=TRUE, nfolds=10,screening=TRUE,no.dim=2)
beta.lassosir <-as.matrix(sir.lasso$beta)
select.lassosir<-list(which(beta.lassosir[,1]!=0),which(beta.lassosir[,2]!=0))


#wls<-wls.sir( x=X.arcene_train_with_intercept,y=Y.arcene_train,categorical=TRUE )

wls<-wls.sir( X,y,categorical=TRUE,ndim = 2 )

beta.wls<-wls$betahat
#select.wls<-wls$select
select.wls<-list(which(beta.lassosir[,1]!=0),which(beta.lassosir[,2]!=0))

round(abs(cor(X%*%beta.lassosir,X%*%beta)),3)
round(abs(cor(X%*%beta.wls,X%*%beta)),3)
round(abs(cor(X%*%beta.chomp,X%*%beta)),3)

#### new examples FOR WLS -------------

#Continus example---- 
n=500;K=50;H=10 ;rho = 0.5 #(0.3, 0.4, 0.8)
tau=0.1
p=300
###seting (1)----
s <-c(1,10,15,20,25,30) # each element of s corresponds to the index of non-zero coefficient in each dimension
Beta <- rep(0, p) ;q<- length(s)
Beta[s] <- 1
landa<-floor(p/sqrt(n))
landa_sum<-80:0+landa
one<-rep(1,p-length(landa_sum))
spik<-diag(c(landa_sum,one)) 
sig=1
#### stting (2)  ------
s <-  c(1,10,15,20,25,30,35,40,45,50) 
#we let G =diag(50+⌈p/√n⌉,49 +⌈p/√n⌉,..,⌈p/√n⌉,1,..,1)  51 spike value 
Beta <- rep(0, p) ; q <- length(s)
Beta[s] <- 1
landa<-floor(p/sqrt(n))
landa_sum<-50:0+landa
one<-rep(1,p-length(landa_sum))
spik<-diag(c(landa_sum,one)) 
sig=1

#### setting 3 -----
s <-1:10 #10 #  each element of s corresponds to the index of non-zero coefficient in each dimension
Beta <- rep(0, p) ; q<- length(s)
Beta[s] <- 1

# setting (4)-------
s <-  c(1,10,15,20,25,30) ;q<- length(s)
Beta <- rep(0, p); 
Beta[s] <- 1

# setting 5-- block diagonal ----- 
s <-1:30 # or 10 each element of s corresponds to the index of non-zero coefficient in each dimension
Beta <- rep(0, p) ; q<- length(s)
Beta[s] <- 1
alpha=c(0.8, 0.2, 0.6)  #(0.5, 0.2, 0.4) ###c(0.8, 0.2, 0.4) #(0.8, 0.2, 0.6)   

######. boucle ----
r = c("Lasso-SIR"=0,"SIR-WLS"=0, "Lasso"=0)#,"nctsir"=0, "SirCHOMP"=0)
result = list("abscos"=r,"Cor"=r ,"Time"=r, "TPR"=r, "FPR"=r ,"FDR"=r) ###, "FNR")
se_result = list("abscos"=r,"Cor"=r , "TPR"=r, "FPR"=r ,"FDR"=r) ###, "FNR")
repetition = 1
abscos=matrix(nrow=3,ncol=repetition)
TPR=matrix(nrow=3,ncol=repetition)
FPR=matrix(nrow=3,ncol=repetition)
#FNR=matrix(nrow=5,ncol=repetition)
Cor=matrix(nrow=3,ncol=repetition)
Time=matrix(nrow=3,ncol=repetition)
FDR=matrix(nrow=3,ncol=repetition)

iteration = 1
set.seed(1)
#omega<-rho^abs(outer(1:p,1:p,"-"))
#set.seed(1234)
rm(X)
while (iteration <= repetition) {
#   X=X_rand(n=n,p=p,q=q,type = "homog",rho)
 X=X_rand(n=n,p=p,q=q,type = "autoR",rho)
#X=X_rand(n=n,p=p,q=q,type = "block",alpha=alpha)
#  X<- mvtnorm::rmvnorm(n, rep(0,p), diag(p))
  #X=X_rand(n=n,p=p,q=q,type = "homog2",rho,S=s)
  #D <- diag(runif(p, 0.5, 2))
  #Sigma <- D %*% omega %*% D
  #X = mvrnorm(n,  rep(0,p), Sigma )
  #t=X%*%Beta
   #Beta[s] <-1 # rnorm(s,mean=0,sd=1) 
  epsilon=rnorm(n,mean=0,sd=1)
#  epsilon=rt(n,10)
  
    #y=as.vector(1+exp(t/sqrt(3))+epsilon )  ## expo
#U<-mvrnorm(n, mu = rep(0,p), Sigma = diag(1,p,p))
#X<-U%*%spik;# dim(X)
#   X=mvrnorm(n, rep(0,p), diag(1,p))
  #  y=X[,1] + X[,2] + X[,3] + X[,4] + X[,5] + X[,6] + rnorm(n)
  t=X%*%Beta
#y <- as.numeric(t + epsilon > 0)
#y<- as.numeric(((t)^3)/2 + epsilon > 0)   ### model 2

  # y=as.vector(1+exp(t/sqrt(3))+epsilon )  ## expo
  #y= as.vector((X%*%Beta[,1]) /(0.5+(X%*%Beta[,2]+1)^2) + sig*epsilon)
  #y= as.vector((X[,s[1]]+X[,s[2]]+1.5*X[,s[3]]+1.2*X[,s[2]]) /(0.5+(X[,s[5]]+1.2*X[,s[6]]+1)^2) + sig*epsilon)
  # X=mvrnorm(n, rep(0,p), diag(1,p))

   #  y=as.vector(exp(t)+epsilon )  ## expo
  y=as.vector(t+epsilon) # linear cas 
# y=as.vector(t+ sig*epsilon  )       # Y<-X[,s[1]]+X[,s[2]]+X[,s[3]]+X[,s[4]]+X[,s[5]]+X[,s[6]] +epsilon
       
 # y=t+ epsilon 
  
#y=as.vector(t+ epsilon)            ### model 1 
 #y= as.vector(((t)^3)/2 + epsilon )    ### model 2
 # y=as.vector(exp(t)+epsilon)         ### model 3 
  #y=as.vector(as.numeric(t+ epsilon >0))
  # y=as.vector(sin(pi*t/2)+epsilon )### model 4
  
  #y=as.vector(sin(t)*exp(t)+epsilon )### model 4
  # y=t/sqrt(3)+ epsilon
  
  ####### Lasso-SIR-
  T1<-Sys.time() 
 sir.lasso <- LassoSIR( X, y, H, no.dim=1,solution.path=FALSE, categorical=FALSE, nfolds=10,screening=FALSE)
 # sir.lasso <- LassoSIR( X, y, H, no.dim=1,solution.path=FALSE, categorical=TRUE, nfolds=10,screening=FALSE)
  beta.lassosir <- sir.lasso$beta
  #beta.lassosir<-normalize_vector(beta.lassosir)
  T2<-Sys.time()
  Time[1,iteration]<-round(difftime(T2, T1, units = "secs"),3)
  abscos[1,iteration]=round(abs(costeta(Beta,beta.lassosir)),3)
  #Dist[1,iteration]=round(Dproj(Beta,beta.lassosir),3)
  if(var(beta.lassosir)!=0) {Cor[1,iteration]=round(abs(cor(X%*%beta.lassosir,t)),3)} else {Cor[1,iteration]=0 };
  TPR[1,iteration]=round(TPR.funct(Beta,beta.lassosir),3)
  FPR[1,iteration]=round(FPR.funct(Beta,beta.lassosir),3)
  FDR[1,iteration]=round((p-q)*FPR[1,iteration]/((p-q)*FPR[1,iteration]+ q*TPR[1,iteration]),3) 
  rm(T1,T2)
  
 
  ###SIR CHOMP --
#  T1<-Sys.time() 
 # beta.chomp<-as.numeric(AdaptCHOMPwithPIC(X, y, d = 1, gamma.pow = 2))
  
  #T2<-Sys.time()
  #Time[2,iteration]<-round(difftime(T2, T1, units = "secs"),3)
  #abscos[2,iteration]=round(abs(costeta(Beta,beta.chomp)),3)
  #Dist[2,iteration]=round(Dproj(Beta,beta.svd),3)
  #if(var(beta.chomp)!=0) {Cor[2,iteration]=round(abs(cor(X%*%beta.chomp,t)),3)} else {Cor[2,iteration]=0 };
  #FPR[2,iteration]=round(FPR.funct(Beta,beta.chomp),3)
  #FNR[2,iteration]=round(FNR.funct(Beta,beta.chomp),3)
  #TPR[2,iteration]=round(TPR.funct(Beta,beta.chomp),3)
  #rm(T1,T2)
  
  ###SIR+NCT natural cannonical tresholding direct --
  T1<-Sys.time() 
  
  # tau_min<-nct_cv_opt(x,y,k,tau_num =10, threshold="soft")
 # beta.nct<-as.vector(sir.nct(X,y,K,tau, threshold="soft"))
 # beta.nct<-normalize_vector(as.vector(sir.nct(X,y,K,tau, threshold="soft")))
  #beta.nct <- beta.nct/sqrt( sum( beta.nct^2 ))#normaliser
  #T2<-Sys.time()
  #Time[3,iteration]<-round(difftime(T2, T1, units = "secs"),3)
  #abscos[3,iteration]=round(abs(costeta(Beta,beta.nct)),3)
  #Dist[3,iteration]=round(Dproj(Beta,beta.nct),3)
  #if(var(beta.nct)!=0) {Cor[3,iteration]=round(abs(cor(X%*%beta.nct,t)),3)} else {Cor[3,iteration]=0 };
  #TPR[3,iteration]=round(TPR.funct(Beta,beta.nct),3)
  #FPR[3,iteration]=round(FPR.funct(Beta, beta.nct),3)
  #FNR[3,iteration]=round(FNR.funct(Beta, beta.nct),3)
  #rm(T1,T2)
  
  ###SIR+WLS  -
  T1<-Sys.time() 
  wlsir<-wls.sir(X,as.vector(y))
#  wlsir<-wls.sir(X,as.vector(y),categorical=TRUE)
  T2<-Sys.time()
  beta.wls<-Re(wlsir$betahat)
  Time[2,iteration]<-round(difftime(T2, T1, units = "secs"),3)
  
  abscos[2,iteration]=round(abs(costeta(Beta, beta.wls)),3)
  # Dist[2,iteration]=round(Dproj(Beta, beta.wls),4)
  if(var(beta.wls)!=0) {Cor[2,iteration]=round(abs(cor(X%*%beta.wls,t)),3)} else {Cor[2,iteration]=0 };
  TPR[2,iteration]=round(TPR.funct(Beta,beta.wls),3)
  FPR[2,iteration]=round(FPR.funct(Beta,beta.wls),3)
  #FNR[2,iteration]=round(FNR.funct(Beta,beta.wls),3)
  FDR[2,iteration]=round((p-q)*FPR[2,iteration]/((p-q)*FPR[2,iteration]+ q*TPR[2,iteration]),3) 
  
  ###lasso  -
  T1<-Sys.time() 
  beta.lasso<-as.vector(Predict(X,y,method="Lasso"))
#  cv.fit<-cv.glmnet(X, y, family = "binomial", alpha = 1,intercept=FALSE)
  # Fit logistic regression with LASSO regularization
 # fit <- glmnet(X, y, family = "binomial", alpha = 1, lambda = cv.fit$lambda.min, intercept = FALSE)
  #beta.lasso<-as.vector(fit$beta) ; 
  T2<-Sys.time()
  Time[3,iteration]<-round(difftime(T2, T1, units = "secs"),3)
  abscos[3,iteration]=round(abs(costeta(Beta, beta.lasso)),3)
  #Dist[3,iteration]=round(Dproj(Beta, beta.lasso),5)
  if(var(beta.lasso)!=0) {Cor[3,iteration]=round(abs(cor(X%*% beta.lasso,t)),3)} else {Cor[3,iteration]=0 };
  TPR[3,iteration]=round(TPR.funct(Beta,beta.lasso),3)
  FPR[3,iteration]=round(FPR.funct(Beta, beta.lasso),3)
  #FNR[3,iteration]=round(FNR.funct(Beta, beta.lasso),3)
  FDR[3,iteration]=round((p-q)*FPR[3,iteration]/((p-q)*FPR[3,iteration]+ q*TPR[1,iteration]),3) 
  iteration = iteration + 1
}
result$abscos["Lasso-SIR"]=sum(abscos[1,]) ;result$Cor["Lasso-SIR"]=sum(Cor[1,]) ;#result$Dist["Lasso-SIR"]=sum(Dist[1,])
result$Time["Lasso-SIR"]=sum(Time[1,]);#result$FNR["Lasso-SIR"]=sum(FNR[1,]) 
result$TPR["Lasso-SIR"]=sum(TPR[1,]);result$FPR["Lasso-SIR"]=sum(FPR[1,]) 
result$FDR["Lasso-SIR"]=sum(FDR[1,])
result$abscos["SIR-WLS"]=sum(abscos[2,]) ;result$Cor["SIR-WLS"]=sum(Cor[2,]) #result$Dist["SIR-WLS"]=sum(Dist[2,]);
result$Time["SIR-WLS"]=sum(Time[2,]) ;#result$FNR["SIR-WLS"]=sum(FNR[2,]) 
result$TPR["SIR-WLS"]=sum(TPR[2,]);result$FPR["SIR-WLS"]=sum(FPR[2,]) 
result$FDR["SIR-WLS"]=sum(FDR[2,])
result$abscos["Lasso"]=sum(abscos[3,]) ;result$Cor["Lasso"]=sum(Cor[3,]) #result$Dist["Lasso"]=sum(Dist[3,]);
result$Time["Lasso"]=sum(Time[3,]);#result$FNR["Lasso"]=sum(FNR[3,]) 
result$TPR["Lasso"]=sum(TPR[3,]);result$FPR["Lasso"]=sum(FPR[3,]) 
result$FDR["Lasso"]=sum(FDR[3,]) 

res=lapply(result,  function(x) x/repetition)
res

# L'écart-type est calculé directement à partir des données de l'échantillon en utilisant la formule :
 # $$s = \sqrt{\frac{\sum_{i=1}^{n} (x_i - \bar{x})^2}{n-1}}$$ sur R par sd()

# L'erreur standard est calculée à partir de l'écart-type de l'échantillon et de la taille de l'échantillon :
 # $$SE = \frac{s}{\sqrt{n}}$$

### compute standar error  ------
se_abscos = apply(abscos, 1, function(x) ifelse(sd(x) == 0, 0, sd(x) / sqrt(ncol(abscos))))
se_Cor = apply(Cor, 1, function(x) ifelse(sd(x) == 0, 0, sd(x) / sqrt(ncol(Cor))))
se_TPR = apply(TPR, 1, function(x) ifelse(sd(x) == 0, 0, sd(x) / sqrt(ncol(TPR))))
se_FPR = apply(FPR, 1, function(x) ifelse(sd(x) == 0, 0, sd(x) / sqrt(ncol(FPR))))
se_FDR = apply(FDR, 1, function(x) ifelse(sd(x) == 0, 0, sd(x) / sqrt(ncol(FDR))))


# Create a list to store the standard errors
se_result$abscos = round(se_abscos,4)*100
se_result$Cor =round( se_Cor,4)*100
se_result$TPR = round(se_TPR,4)*100
se_result$FPR =round(se_FPR,4)*100
se_result$FDR = round(se_FDR,4)*100


# Print the standard errors
se_result


# Your data

# Define the settings
settings <- c(300, 700, 1000, 1500, 2000)

# Create a data frame to hold all values
table_data <- expand.grid(Setting = settings, Method = names(res$abscos))
table_data$PC <- rep(res$abscos, length(settings))
table_data$Cor <- rep(res$Cor, length(settings))
table_data$TPR <- rep(res$TPR, length(settings))
table_data$FPR <- rep(res$FPR, length(settings))
table_data$FDR <- rep(res$FDR, length(settings))

# Function to format the table in LaTeX
generate_latex_table <- function(data) {
  cat("\\begin{tabular}{ |c|c|c|c|c|c| }\n")
  cat("\\hline\n")
  cat("Setting & Method & PC & Cor & TPR & FPR & FDR \\\\\n")
  cat("\\noalign{\\smallskip}\\hline\\noalign{\\smallskip}\n")
  
  for (setting in unique(data$Setting)) {
    first_row <- TRUE
    for (method in unique(data$Method)) {
      row <- data[data$Setting == setting & data$Method == method, ]
      if (first_row) {
        cat(paste0("$p=", setting, "$ & ", method, " & $", 
                   sprintf("%.2f", row$PC), "$ & $", 
                   sprintf("%.2f", row$Cor), "$ & $", 
                   sprintf("%.1f", row$TPR), "$ & $", 
                   sprintf("%.2f", row$FPR), "$ & $", 
                   sprintf("%.2f", row$FDR), "$ \\\\\n"))
        first_row <- FALSE
      } else {
        cat(paste0(" & ", method, " & $", 
                   sprintf("%.2f", row$PC), "$ & $", 
                   sprintf("%.2f", row$Cor), "$ & $", 
                   sprintf("%.1f", row$TPR), "$ & $", 
                   sprintf("%.2f", row$FPR), "$ & $", 
                   sprintf("%.2f", row$FDR), "$ \\\\\n"))
      }
    }
    cat("\\noalign{\\smallskip}\\hline\\noalign{\\smallskip}\n")
  }
  
  cat("\\end{tabular}\n")
}

# Generate the LaTeX table
generate_latex_table(table_data)


 ### FDR <-FP/(FP+TP)------
### for me FDR <-(p-q)*FPR/((p-q)*FPR+q*TPR) 
# Given values
#p <- total number of variables
#q <- number of true variables
#FPR <- False Positive Rate
#TPR <- True Positive Rate

# Estimate FDR using available information
#FDR <- (p - q) * FPR / ((p - q) * FPR + q * TPR)

## round((p-q)*0.031/((p-q)*0.031+ q*1),3)


#### ROC comparison --------
#lasso 
select.lasso<-which(beta.lasso!=0)
linear_predictor_lasso <- X %*% beta.lasso
predicted_probabilities_lasso<- as.vector(1 / (1 + exp(-linear_predictor_lasso)))
#Lasso-SIR 
select.lassosir<-which(beta.lassosir!=0)
linear_predictor_lassoSIR <- X %*% beta.lassosir 
predicted_probabilities_lassoSIR<- as.vector(1 / (1 + exp(-linear_predictor_lassoSIR)))
select.wls<-which(beta.wls!=0)
linear_predictor_wls <-X%*% beta.wls 
predicted_probabilities_wls<- as.vector(1 / (1 + exp(-linear_predictor_wls)))


# Plot ROC curve
roc_curve_lasso <- roc(y,predicted_probabilities_lasso, plot=TRUE, legacy.axes=TRUE)
roc_curve_lassoSIR<- roc(y,predicted_probabilities_lassoSIR, plot=FALSE, legacy.axes=TRUE)
roc_curve_wls <- roc(y,predicted_probabilities_wls, plot=TRUE, legacy.axes=TRUE)
#plot(roc_curve_wls, main = "ROC Curve for wlsSIR Logistic Regression", col = "red")

#### compare the three 

par(mfrow = c(1, 1))
plot(roc_curve_wls, main = "ROC Curve ", col = "red")
lines(roc_curve_lassoSIR, col = "blue")
lines(roc_curve_lasso, col = "green")
legend("bottomright", legend = c("ROC Curve_wls","ROC Curve_lassoSIR", "ROC Curve Lasso"), 
       col = c( "red","blue","green"), lwd = 2,cex = 0.6)
# Calculate AUC
auc_value_wls<- auc(roc_curve_wls); auc_value_wls
auc_value_lassSIR<-auc(roc_curve_lassoSIR );auc_value_lassSIR
auc_value_lasso<-auc(roc_curve_lasso);auc_value_lasso


#### plot comparison ------
t1=X%*%beta.lasso
t2=X%*%beta.lassosir
t3=X%*%beta.wls
#t=x%*%b1+(x%*%b2)^2
#t1=x%*%beta[,1]+(x%*%beta[,2])^2
#t2=x%*%dr$evectors[,1]+(x%*%dr$evectors[,2])^2
#y1= as.vector(((t1)^3)/2)    ### model 2
#y4= as.vector(((t4)^3)/2)    ### model 2
#y5= as.vector(((t5)^3)/2)    ### model 2
#par(mgp=c(2.5, 1, 0))
######### Y=f(XBeta)--in same ggplot ----
plot_real <-ggplot() +
  geom_point(aes(x = t, y = y), color = "red", shape = 16, size = 2) +
  labs(title = " real direction",  x = TeX('$X \\beta $'), y = "y") +
  theme_minimal()
plot_lasso <-ggplot() +
  geom_point(aes(x = t1, y = y), color = "orange", shape = 16, size = 2) +
  labs(title = "Lasso", x = TeX('$X \\hat{\\beta}_{Lasso} $') , y = "y") +
  theme_minimal()
plot_lassosir <-ggplot() +
  geom_point(aes(x = t2, y = y), color = "blue", shape = 16, size = 2) +
  labs(title = "Lasso-SIR", y= "y", x = TeX('$X \\hat{\\beta}_{Lasso-SIR} $')) +
  theme_minimal()
plot_wlsSIR <-ggplot() +
  geom_point(aes(x = t3, y = y), color = "green", shape = 16, size = 2) +
  labs(title = "SIR-WLS", x = TeX('$X \\hat{\\beta}_{SIR-WLS} $'), y = "y") +
  theme_minimal()
grid.arrange(plot_real, plot_lasso,plot_lassosir, plot_wlsSIR, 
             ncol = 2)




####  Linearity test ------
#plot_real <-ggplot() +
 # geom_point(aes(x = t, y = y), color = "red", shape = 16, size = 2) +
#  labs(title = " real direction", x = expression(X~beta), y = "y") +
 # theme_minimal()
plot_lasso <-ggplot() +
  geom_point(aes(x = t, y = t1), color = "orange", shape = 16, size = 2) +
    labs(title = "Lasso",  x = TeX('$X\\beta $') ,  y = TeX('$X \\hat{\\beta}_{Lasso} $')) +
  theme_minimal()
plot_lassosir <-ggplot() +
  geom_point(aes(x = t, y =- t2), color = "blue", shape = 16, size = 2) +
  #geom_line(aes(x = t, y = t), color = "red",linetype = 4) +
  labs(title = "Lasso-SIR", y=TeX('$X \\hat{\\beta}_{Lasso-SIR}$') , x= TeX('$X\\beta $')) +
  theme_minimal()
plot_wlsSIR <-ggplot() +
  geom_point(aes(x =t, y =t3), color = "green", shape = 16, size = 2) +
  #geom_line(aes(x = t, y = t), color = "red",linetype = 2) +
  labs(title = "SIR-WLS", x = TeX('$X\\beta $'), y=TeX('$X \\hat{\\beta}_{SIR-WLS}$') ) +
  theme_minimal()
grid.arrange( plot_lasso,plot_lassosir, plot_wlsSIR, 
             ncol = 3)


##### categaritocal responese  ------

# Simulate some data----
set.seed(123)
n=500; p=700;K=50;H=10;q=20; alpha=c(0.2, 0.3, 0.7) ;rho = 0.5 #(0.3, 0.4, 0.8)
tau=0.1 
s <-c(1,10,15,20,25,30) # each element of s corresponds to the index of non-zero coefficient in each dimension
#s <-  c(1,10,20,30,40)# each element of s corresponds to the index of non-zero coefficient in each dimension
#s<-1:20
#s <-  c(1,10,15,20,25,30,35,40,45,50)
#sig=1
#Beta<-c(rep(0,p))
#s <- 5
Beta <- rep(0, p)
Beta [s]<-1
#Beta[1:s] <- sample(c(-1,1), s, replace=TRUE)*runif(s, 1, 1.5)
#Beta<-matrix(0,ncol=2,nrow=p)
#s1<-c(1,10,15,20);s2<-c(25,30)
#val1<-c(1,1,1.5,1.2)
x <- matrix(rnorm(n * p), ncol = p)
#beta_true <- c(1, 2, 0, 0, 0, 0, 0, 0, 0, 0)  # True coefficients
y <- rbinom(n, 1, plogis(x %*% Beta))
cv.fit<-cv.glmnet(x, y, family = "binomial", alpha = 1,intercept=FALSE)
# Fit logistic regression with LASSO regularization
fit <- glmnet(x, y, family = "binomial", alpha = 1, lambda = cv.fit$lambda.min, intercept = FALSE)
#beta.lasso<-as.vector(coef(fit))
beta.lasso<-as.vector(fit$beta) ; 
which(beta.lasso!=0)
probabilities <- predict(fit, newx = x, type = "response")
linear_predictor_lasso <-x%*%beta.lasso

linear_predictor_lasso <-cbind(1, x)%*%beta.lasso

predicted_probabilities_lasso<- as.vector(1 / (1 + exp(-linear_predictor_lasso)))

#Plot ROC curve
roc_curve_lasso <- roc(y,predicted_probabilities_lasso, plot=TRUE, legacy.axes=TRUE)
roc_curve_lasso2<-roc(y,predicted_probabilities_lasso, plot=TRUE, legacy.axes=TRUE, colors="red")
plot(roc_curve_lasso)
lines(roc_curve_lasso2,col="red")

## lassoSIR 
## without intersept 
sir.lasso <- LassoSIR(x,y,solution.path=FALSE, categorical=TRUE, nfolds=10,screening=TRUE,no.dim=1)
beta.lassosir <-as.vector(sir.lasso$beta)
select.lassosir<-which(beta.lassosir!=0)
linear_predictor_lassoSIR <- x %*% beta.lassosir 
predicted_probabilities_lassoSIR<- as.vector(1 / (1 + exp(-linear_predictor_lassoSIR)))

roc_curve_lassoSIR <- roc(y,predicted_probabilities_lassoSIR, plot=FALSE, legacy.axes=TRUE)
plot(roc_curve_lassoSIR)
## with intercept 
sir.lasso <- LassoSIR(cbind(1, x),y,solution.path=FALSE, categorical=TRUE, nfolds=10,screening=TRUE,no.dim=1)
beta.lassosir <-as.vector(sir.lasso$beta)
select.lassosir<-which(beta.lassosir!=0)
linear_predictor_lassoSIR <- cbind(1, x) %*% beta.lassosir 
predicted_probabilities_lassoSIR<- as.vector(1 / (1 + exp(-linear_predictor_lassoSIR)))

roc_curve_lassoSIR2 <- roc(y,predicted_probabilities_lassoSIR, plot=FALSE, legacy.axes=TRUE)
plot(roc_curve_lassoSIR2)
lines(roc_curve_lassoSIR,col="red")


### wls without intercept
wls<-wls.sir( x,y,categorical=TRUE )
beta.wls<-wls$betahat
which(beta.wls!=0)
linear_predictor_wls <-x %*% beta.wls 
predicted_probabilities_wls<- as.vector(1 / (1 + exp(-linear_predictor_wls)))
# Plot ROC curve
roc_curve_wls <- roc(y,predicted_probabilities_wls, plot=TRUE, legacy.axes=TRUE)
plot(roc_curve_wls, main = "ROC Curve for wlsSIR Logistic Regression", col = "red")

### with intercept :
wls<-wls.sir( cbind(1, x),y,categorical=TRUE )
beta.wls<-wls$betahat
which(beta.wls!=0)
linear_predictor_wls <-cbind(1, x) %*% beta.wls 
predicted_probabilities_wls<- as.vector(1 / (1 + exp(-linear_predictor_wls)))
# Plot ROC curve
roc_curve_wls2<- roc(y,predicted_probabilities_wls, plot=TRUE, legacy.axes=TRUE)
plot(roc_curve_wls, main = "ROC Curve for wlsSIR Logistic Regression", col = "red")
lines(roc_curve_wls2)

#### compare the tree 

plot(roc_curve_wls, main = "ROC Curve ", col = "red")
lines(roc_curve_lassoSIR, col = "blue")
lines(roc_curve_lasso, col = "green")

# Calculate AU
auc_value_wls<- auc(roc_curve_wls); auc_value_wls
auc_value_lassSIR<-auc(roc_curve_lassoSIR);auc_value_lassSIR
auc_value_lasso<-auc(roc_curve_lasso);auc_value_lasso

auc_value_wls<- auc(roc_curve_wls); auc_value_wls
Area under the curve: 0.9118
> auc_value_lassSIR<-auc(roc_curve_lassoSIR);auc_value_lassSIR
Area under the curve: 0.9363
> auc_value_lasso<-auc(roc_curve_lasso);auc_value_lasso
Area under the curve: 0.9242

#with intercept-------
plot(roc_curve_wls2,  col = "red")
lines(roc_curve_lassoSIR2, col = "blue")
lines(roc_curve_lasso2, col = "green")

legend("bottomright", legend = c("ROC Curve_lassoSIR", "ROC Curve_wls","ROC Curve Lasso"), 
       col = c("blue", "red","green"), lwd = 2,cex = 0.6)
# Calculate AUC
auc_value_wls<- auc(roc_curve_wls2); auc_value_wls
auc_value_lassSIR<-auc(roc_curve_lassoSIR2);auc_value_lassSIR
auc_value_lasso<-auc(roc_curve_lasso2);auc_value_lasso
par(mfrow = c(1, 3))
plot(x%*%Beta,linear_predictor_lasso)
plot(x%*%Beta,linear_predictor_lassoSIR)
plot(x%*%Beta,linear_predictor_wls)



auc_value_wls<- auc(roc_curve_wls2); auc_value_wls
Area under the curve: 0.9071
> auc_value_lassSIR<-auc(roc_curve_lassoSIR2);auc_value_lassSIR
Area under the curve: 0.9363
> auc_value_lasso<-auc(roc_curve_lasso2);auc_value_lasso
Area under the curve: 0.9242
############### double INDIX example :-----
n=500; p=700;K=50;H=10 ;rho = 0.5 #(0.3, 0.4, 0.8);tau=0.1
sig<-1
Beta<-matrix(0,ncol=2,nrow=p)
s1<-1:10;s2<-1:6

Beta[s1,1]<-1
Beta[s2,2]<-1
#Beta[s]<-1
landa<-floor(p/sqrt(n))
#landa_sum<-80:0+landa
landa_sum<-50:0+landa
one<-rep(1,p-length(landa_sum))
spik<-diag(c(landa_sum,one)) 
epsilon=rnorm(n,mean=0,sd=1)

 U<-mvrnorm(n, mu = rep(0,p), Sigma = diag(1,p,p))
X<-U%*%spik;# dim(X)
t=X%*%Beta
y= as.vector(X%*%Beta[,1] +(0.5+(X%*%Beta[,2]+1)^3) + sig*epsilon) 

T1<-Sys.time() 
sir.lasso <- LassoSIR( X, y, H, no.dim=2,solution.path=FALSE, categorical=FALSE, nfolds=10,screening=FALSE)
beta.lassosir <- sir.lasso$beta
wlsir<-wls.sir(X,as.vector(y), ndim=2)
T2<-Sys.time()
beta.wls<-Re(wlsir$betahat)
###
round(abs(costeta(Beta,beta.lassosir)),3)
round(abs(costeta(Beta,beta.wls)),3)

### cor  round(abs(cor(X%*%beta.lassosir,X%*%beta)
round(abs(cor(X%*%beta.lassosir,X%*%Beta)),3)
round(abs(cor(X%*%beta.wls,X%*%Beta)),3)
### FPR , TPR 
round(TPR.funct(Beta[,1],beta.lassosir[,1]),3);round(FPR.funct(Beta,beta.lassosir),3)
round(TPR.funct(Beta,beta.wls),3);round(FPR.funct(Beta,beta.wls),3)

wls.slect<-which(beta.wls!=0)

lassosir.slect<-which(beta.lassosir!=0)

#plot multiple indice ------
t2=X%*%beta.lassosir
t1=X%*%beta.wls



######### Y=f(XBeta)--in same ggplot ----
#fisr direc 
plot_real <-ggplot() +
  #geom_point(aes(x = t[,2], y = y), color = "red", shape = 16, size = 2) +
geom_point(aes(x = t[,2], y = y), color = "red", shape = 16, size = 2) +
  labs(title = " real direction",  x = TeX('$X \\beta_1$'), y = "y") +
  theme_minimal()

plot_lassosir <-ggplot() +
  geom_point(aes(x = t2[,1], y = y), color = "blue", shape = 16, size = 2) +
  #geom_point(aes(x = t2[,2], y = y), color = "blue", shape = 16, size = 2) +
  labs(title = "Lasso-SIR", y= "y", x = TeX('$X \\hat{\\beta}_{LASSOSIR}_1 $')) +
  theme_minimal()
plot_wlsSIR <-ggplot() +
 # geom_point(aes(x = t1[,2], y = y), color = "green", shape = 16, size = 2) +
  geom_point(aes(x = -t1[,1], y = y), color = "green", shape = 16, size = 2) +
  labs(title = "SIR-WLS", x = TeX('$X \\hat{\\beta}_{SIR-WLS_1} $'), y = "y") +
  theme_minimal()
grid.arrange(plot_real,plot_lassosir, plot_wlsSIR, 
             ncol = 3)

#second direc : 
plot_real <-ggplot() +
  #geom_point(aes(x = t[,2], y = y), color = "red", shape = 16, size = 2) +
  geom_point(aes(x = t[,1], y = y), color = "red", shape = 16, size = 2) +
  labs(title = " real direction",  x = TeX('$X \\beta_1$'), y = "y") +
  theme_minimal()

plot_lassosir <-ggplot() +
  geom_point(aes(x = t2[,2], y = y), color = "blue", shape = 16, size = 2) +
  #geom_point(aes(x = t2[,2], y = y), color = "blue", shape = 16, size = 2) +
  labs(title = "Lasso-SIR", y= "y", x = TeX('$X \\hat{\\beta}_{LASSOSIR}_1 $')) +
  theme_minimal()
plot_wlsSIR <-ggplot() +
  # geom_point(aes(x = t1[,2], y = y), color = "green", shape = 16, size = 2) +
  geom_point(aes(x = -t1[,2], y = y), color = "green", shape = 16, size = 2) +
  labs(title = "SIR-WLS", x = TeX('$X \\hat{\\beta}_{SIR-WLS_1} $'), y = "y") +
  theme_minimal()
grid.arrange(plot_real,plot_lassosir, plot_wlsSIR, 
             ncol = 3)



######### Y=f(XBeta)--in same ggplot ----
plot_real <-ggplot() +
  #geom_point(aes(x = t[,2], y = y), color = "red", shape = 16, size = 2) +
  geom_point(aes(x = t[,1], y = y), color = "red", shape = 16, size = 2) +
  labs(title = " real direction",  x = TeX('$X \\beta_1$'), y = "y") +
  theme_minimal()

plot_lassosir <-ggplot() +
  geom_point(aes(x = t2[,2], y = y), color = "blue", shape = 16, size = 2) +
  #geom_point(aes(x = t2[,2], y = y), color = "blue", shape = 16, size = 2) +
  labs(title = "Lasso-SIR", y= "y", x = TeX('$X \\hat{\\beta}_{LASSOSIR}_1 $')) +
  theme_minimal()
plot_wlsSIR <-ggplot() +
  # geom_point(aes(x = t1[,2], y = y), color = "green", shape = 16, size = 2) +
  geom_point(aes(x = -t1[,2], y = y), color = "green", shape = 16, size = 2) +
  labs(title = "SIR-WLS", x = TeX('$X \\hat{\\beta}_{SIR-WLS_1} $'), y = "y") +
  theme_minimal()
grid.arrange(plot_real,plot_lassosir, plot_wlsSIR, 
             ncol = 3)


library(plotly)

# Create data

library(plotly)

# Assuming t3 is your dataset containing the coordinates and y-values
data <- data.frame(x = t[,1], y = t[,2], z = y)

# Create the 3D scatter plot
plot_real <- plot_ly(data, x = ~x, y = ~y, z = ~z, type = "scatter3d", mode = "markers", 
                        marker = list(color = "green", size = 3)) %>%
  layout(scene = list(title = "3D Scatter Plot",
                      aspectratio = list(x = 2, y = 1, z = 1)), # Adjust aspect ratio here
         title = "SIR-WLS Results")

# Print the plot
plot_real



####  Linearity test ------
#plot_real <-ggplot() +
# geom_point(aes(x = t, y = y), color = "red", shape = 16, size = 2) +
#  labs(title = " real direction", x = expression(X~beta), y = "y") +
# theme_minimal()
plot_lasso <-ggplot() +
  geom_point(aes(x = t, y = t1), color = "orange", shape = 16, size = 2) +
  labs(title = "Lasso",  x = TeX('$X\\beta $') ,  y = TeX('$X \\hat{\\beta}_{LASSO} $')) +
  theme_minimal()
plot_lassosir <-ggplot() +
  geom_point(aes(x = t, y = t2), color = "blue", shape = 16, size = 2) +
  #geom_line(aes(x = t, y = t), color = "red",linetype = 4) +
  labs(title = "Lasso-SIR", y=TeX('$X \\hat{\\beta}_{LASSOSIR}$') , x= TeX('$X\\beta $')) +
  theme_minimal()
plot_wlsSIR <-ggplot() +
  geom_point(aes(x = t, y = t3), color = "green", shape = 16, size = 2) +
  #geom_line(aes(x = t, y = t), color = "red",linetype = 2) +
  labs(title = "SIR-WLS", x = TeX('$X\\beta $'), y=TeX('$X \\hat{\\beta}_{SIR-WLS}$') ) +
  theme_minimal()
grid.arrange( plot_lasso,plot_lassosir, plot_wlsSIR, 
              ncol = 3)




##beta.chomp<-as.numeric( AdaptCHOMPwithPIC(X, y, d = 1, gamma.pow = 1))--------
fit.AdaptCHOMP_WeightPower2 <- as.numeric(AdaptCHOMPwithPIC(X, y, d = 1, gamma.pow = 2))

select.chomp<-which(beta.chomp!=0)

sir.lasso <- LassoSIR( X,y,solution.path=FALSE, categorical=TRUE, nfolds=10,screening=TRUE,no.dim=1)
beta.lassosir <-as.vector(sir.lasso$beta)
select.lassosir<-which(beta.lassosir!=0)



#wls<-wls.sir( x=X.arcene_train_with_intercept,y=Y.arcene_train,categorical=TRUE )

wls<-wls.sir( X,y,categorical=TRUE )

beta.wls<-wls$betahat
#select.wls<-wls$select
beta.slect<-which(beta.wls!=0)

###
round(abs(costeta(beta,beta.lassosir)),3)
round(abs(costeta(beta,beta.wls)),3)
#round(abs(costeta(beta,beta.chomp)),3)

### cor  round(abs(cor(X%*%beta.lassosir,X%*%beta)
round(abs(cor(X%*%beta.lassosir,X%*%beta)),3)
round(abs(cor(X%*%beta.wls,X%*%beta)),3)
#round(abs(cor(X%*%beta.chomp,X%*%beta)),3)
### FPR , TPR 
round(TPR.funct(beta,beta.lassosir),3);round(FPR.funct(beta,beta.lassosir),3)
round(TPR.funct(beta,beta.wls),3);round(FPR.funct(beta,beta.wls),3)
#round(TPR.funct(beta,beta.chomp),3);round(FPR.funct(beta,beta.chomp),3)
