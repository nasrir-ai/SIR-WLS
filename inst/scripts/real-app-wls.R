source("~/Desktop/phd /R codes /myfunctions.R")
#### arcene ------
X.arcene_train<- as.matrix( read.table("~/Desktop/phd /R codes /arcene/ARCENE/arcene_train.data", quote="\"", comment.char=""))
Y.arcene_train <- as.matrix(read.table("~/Desktop/phd /R codes /arcene/ARCENE/arcene_train.labels", quote="\"", comment.char=""))
#data.arcene<-as.matrix(bind_cols(X.arcene_train,Y.arcene_train))
#X.arcene_test <- as.matrix(read.table("~/Downloads/arcene/ARCENE/arcene_train.data", quote="\"", comment.char=""))
Y.arcene_valid <- as.matrix(read.table("~/Desktop/phd /R codes /arcene/ARCENE/arcene_valid.labels", quote="\"", comment.char=""))
X.arcene_valid <-as.matrix( read.table("~/Desktop/phd /R codes /arcene/ARCENE/arcene_valid.data", quote="\"", comment.char=""))

X.arcene_train_with_intercept <- cbind(1, X.arcene_train)
X.arcene_valid_with_intercept <- cbind(1, X.arcene_valid)

#beta.chomp <- AdaptCHOMPwithPIC( X.arcene_train,Y.arcene_train, d = 1, gamma.pow = 1)
#fit.AdaptCHOMP_WeightPower2 <- AdaptCHOMPwithPIC(X, y, d = 1, gamma.pow = 2)


sir.lasso <- LassoSIR( X.arcene_train_with_intercept,Y.arcene_train,solution.path=FALSE, categorical=TRUE, nfolds=10,no.dim=1)
#sir.lasso <- LassoSIR( X.arcene_train,Y.arcene_train,solution.path=FALSE, categorical=TRUE, nfolds=10,screening=TRUE,no.dim=1)
beta.lassosir <-as.vector(sir.lasso$beta)
#select.lassosir<-which(beta.lassosir!=0)


wls<-wls.sir( x=X.arcene_train_with_intercept,y=Y.arcene_train,categorical=TRUE )
 
#wls<-wls.sir( x=X.arcene_train,y=Y.arcene_train,categorical=TRUE )

beta.wls<-wls$betahat
#select.wls<-wls$select
#beta.slect<-which(beta.wls!=0)
#beta.rsir<- do.rsir(X.arcene_train,Y.arcene_train)


### lasso 
lasso.cv<-cv.glmnet( x=X.arcene_train,y=Y.arcene_train,alpha =1, family="binomial", nfolds=10, intercept=TRUE)   ## cross validation to obtain optimal mu
best_lam<-lasso.cv$lambda.min
fit<-glmnet( x=X.arcene_train,y=Y.arcene_train, alpha =1, family="binomial", lambda = best_lam,itercept=TRUE)
#beta.lasso <- as.vector(fit$beta)
beta.lasso <- as.vector(coef(fit))

#beta.slect.lasso<-which(beta.lasso!=0)
##### ROC curve 
## for lassoSIR

linear_predictor_lassoSIR <- X.arcene_valid_with_intercept %*% beta.lassosir 
#linear_predictor_lassoSIR <- X.arcene_train %*% beta.lassosir 
predicted_probabilities_lassoSIR<- as.vector(1 / (1 + exp(-linear_predictor_lassoSIR)))
# Plot ROC curve
library(pROC)
roc_curve_LassoSIR <- roc(Y.arcene_train,predicted_probabilities_lassoSIR, plot=TRUE, legacy.axes=TRUE)
roc_curve_LassoSIR <- roc(Y.arcene_valid,predicted_probabilities_lassoSIR, plot=TRUE, legacy.axes=TRUE)
#plot(roc_curve_LassoSIR, main = "ROC Curve for LassoSIR Logistic Regression", col = "blue")
### for wls_SIR 
linear_predictor_wls <- X.arcene_valid_with_intercept %*% beta.wls 
#linear_predictor_wls <- X.arcene_train %*% beta.wls 
predicted_probabilities_wls<- as.vector(1 / (1 + exp(-linear_predictor_wls)))

# Plot ROC curve
roc_curve_wls <- roc(Y.arcene_valid,predicted_probabilities_wls, plot=TRUE, legacy.axes=TRUE)

plot(roc_curve_wls, main = "ROC Curve for wlsSIR Logistic Regression", col = "red")
lines(roc_curve_LassoSIR, col = "blue")
### for lasso 
linear_predictor_lasso <- X.arcene_valid_with_intercept %*% beta.lasso 
#linear_predictor_lasso <- X.arcene_train %*% beta.lasso 
predicted_probabilities_lasso<- as.vector(1 / (1 + exp(-linear_predictor_lasso)))

 #Plot ROC curve
roc_curve_lasso <- roc(Y.arcene_train,predicted_probabilities_lasso, plot=TRUE, legacy.axes=TRUE)
par(mfrow = c(1, 1))
plot(roc_curve_wls, main = "ROC Curve ", col = "red")
lines(roc_curve_LassoSIR, col = "blue")
lines(roc_curve_lasso, col = "green")
legend("bottomright", legend = c("ROC Curve_wls","ROC Curve_lassoSIR", "ROC Curve Lasso"), 
       col = c( "red","blue","green"), lwd = 2,cex = 0.6)
# Calculate AUC
auc_value_wls<- auc(roc_curve_wls); auc_value_wls
auc_value_lassSIR<-auc(roc_curve_LassoSIR);auc_value_lassSIR
auc_value_lasso<-auc(roc_curve_lasso);auc_value_lasso


#plot(linear_predictor_wls)
#> auc_value_wls<- auc(roc_curve_wls); auc_value_wls
#Area under the curve: 0.8141
#> auc_value_lassSIR<-auc(roc_curve_LassoSIR);auc_value_lassSIR
#Area under the curve: 0.7853



### curve comparison
# Plot the first ROC curve
plot(roc_curve_LassoSIR, col = "blue", lwd = 2, main = "ROC Curves Comparison", 
     xlab = "False Positive Rate", ylab = "True Positive Rate")
# Add the second ROC curve to the same plot
lines(roc_curve_wls, col = "red", lwd = 2,  )
# Add legend
legend("bottomright", legend = c("ROC Curve_lassoSIR", "ROC Curve_wls"), 
       col = c("blue", "red"), lwd = 2,cex = 0.6)



library(ggplot2)

# Create a data frame for the first ROC curve
df1 <- data.frame(Specificity = roc_curve_LassoSIR$specificities,
                  Sensitivity = roc_curve_LassoSIR$sensitivities)

# Create a data frame for the second ROC curve
df2 <- data.frame(Specificity = roc_curve_wls$specificities,
                  Sensitivity = roc_curve_wls$sensitivities)

# Plot ROC curves using ggplot2
ggplot() +
  geom_line(data = df1, aes(x = Specificity, y = Sensitivity), color = "blue") +
 geom_line(data = df2, aes(x = Specificity, y = Sensitivity), color = "red") +
  labs(title = "ROC Curves Comparison",
       x = "Specificity",
       y = "Sensitivity") +
  theme_minimal() +
  scale_color_manual(values = c("blue", "red")) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.key = element_rect(fill = "green", color = "green"))



common_elements1<- Reduce(intersect,list(select.lasso,select.lassosir))
common_elements2 <- Reduce(intersect,list(select.lasso,select.wls))
# Display common elements
common_elements1
common_elements2
### #############@@
###### APlied for wine data ------
#M. Lichman. UCI machine learning repository, 2013. URL http://archive.ics.uci.edu/ml.
#https://archive.ics.uci.edu/dataset/109/wine

wine <- as.matrix(read.csv("~/Desktop/phd /R codes /wine/wine.data", header=FALSE))
which(wine[,1]==1)
which(wine[,1]==2)
which(wine[,1]==3)
 x.wine<-wine[,2:14]; y.wine<-wine[,1]
## ndim=1
smp_size <- floor(0.70 * nrow(x.wine))
set.seed(36)

#### Deviser les données en deux echontillons  "training set (donnée d'entraînement)" et "test set"  
train_ind <- sample(seq_len(nrow(wine)), size = smp_size)
x.wine.train <- as.matrix(x.wine[train_ind, ])
x.wine.test <- as.matrix(x.wine[-train_ind, ] )
y.wine.train<-as.vector(y.wine[train_ind])
y.wine.test<-as.vector(y.wine[-train_ind])

#sir.lassosir <- LassoSIR( wine[,2:14],wine[,1],solution.path=FALSE, categorical=TRUE, nfolds=10,screening=FALSE,no.dim=1)
## ndim=2 
sir.lassosir <- LassoSIR( x.wine.train ,y.wine.train,solution.path=FALSE, categorical=TRUE, nfolds=10,no.dim=2)
beta.lassosir<-sir.lassosir$beta
beta.lassosir<-normalize_matrix(beta.lassosir)
#lasso.sir.slct<-which(beta.lassosir!=0)
pc.lassoSIR<-x.wine.test%*%beta.lassosir
 #wls<-wls.sir(x= wine[,2:14], y=wine[,1], categorical=TRUE,ndim=1) # ndim=1
 wls<-wls.sir(x= x.wine.train, y=y.wine.train, categorical=TRUE,ndim=2)
#beta.wls<-as.vector(wls$betahat)
beta.wls<-wls$betahat
beta.wls<-normalize_matrix(beta.wls)
#beta.wls
#wls.slct<-wls$select # wls.slct 
pc.wls.sir<-x.wine.test%*%beta.wls
plot(pc.wls.sir[,1],pc.wls.sir[,2])
plot(pc.lassoSIR[,1],pc.lassoSIR[,2])

out.sir<-do.sir( x.wine.train, y.wine.train)#, categorical=TRUE,no.dim=2)
beta.sir<-normalize_matrix(out.sir$projection)
pc.sir<-x.wine.test%*%beta.sir



# Plot for pc.lassoSIR
class_labels <- y.wine.test
pc.lassoSIR_df <- as.data.frame(pc.lassoSIR)
pc.wls.sir_df <- as.data.frame(pc.wls.sir)
pc.sir_df <- as.data.frame(pc.sir)

plot_LassoSIR <- ggplot() +
  geom_point(data = pc.lassoSIR_df, aes(x = pc.lassoSIR[, 1], y = pc.lassoSIR[, 2], color = factor(class_labels)),show.legend = FALSE) +
  labs(x = TeX('$PC_1$'), y = TeX('$PC_2$')) +
  ggtitle("Lasso-SIR") +
  theme_minimal()

# Specifying colors for each class
plot_LassoSIR <- plot_LassoSIR + scale_color_manual(values = c( "red", "green","blue"))

# Plot for pc.wls.sir
plot_wlsSIR <- ggplot() +
  geom_point(data = pc.wls.sir_df, aes(y= pc.wls.sir[, 1],x = pc.wls.sir[, 2], color = factor(class_labels)),show.legend = FALSE) +
  labs(x = TeX('$PC_1$'), y = TeX('$PC_2$')) +
  ggtitle("SIR-WLS") +
  theme_minimal()

# Specifying colors for each class
plot_wlsSIR <- plot_wlsSIR + scale_color_manual(values = c( "red", "green","blue"))

#### for SIR 

# Plot for pc.wls.sir
plot_sir <- ggplot() +
  geom_point(data = pc.sir_df, aes(y = pc.sir[, 1], x = pc.sir[, 2], color = factor(class_labels)),show.legend = FALSE) +
  labs(x = TeX('$PC_1 $'), y = TeX('$PC_2$')) +
  ggtitle("SIR") +
  theme_minimal()

# Specifying colors for each class
plot_sir <- plot_sir + scale_color_manual(values = c( "red", "green","blue"))

library(gridExtra)

grid.arrange(ggplotGrob(plot_LassoSIR + theme(legend.position = "none")),
             ggplotGrob(plot_wlsSIR + theme(legend.position = "none")),
             ggplotGrob(plot_sir + theme(legend.position = "none")) , ncol = 3)

####----- boxplots------
# Convert scatter plots to box plots
plot_LassoSIR_boxplot <- ggplot(pc.lassoSIR_df, aes(x = factor(class_labels), y = pc.lassoSIR[, 2], fill = factor(class_labels))) +
  geom_boxplot() +
  labs(y = TeX('$PC_2$'), x = TeX('Factors')) +
  ggtitle("Lasso-SIR") +
  theme_minimal() +
  scale_fill_manual(values = c("red", "green", "blue"))
plot_wlsSIR_boxplot <- ggplot(pc.wls.sir_df, aes(x = factor(class_labels), y = pc.wls.sir[, 2], fill = factor(class_labels))) +
  geom_boxplot() +
  labs(y = TeX('$PC_2$'), x = TeX('Factors')) +
  ggtitle("SIR-WLS") +
  theme_minimal() +
  scale_fill_manual(values = c("red", "green", "blue"))
plot_sir_boxplot <- ggplot(pc.sir_df, aes(x = factor(class_labels), y = pc.sir[, 2], fill = factor(class_labels))) +
  geom_boxplot() +
  labs(y = TeX('$PC_2$'), x = TeX('Factors')) +
  ggtitle("SIR") +
  theme_minimal() +
  scale_fill_manual(values = c("red", "green", "blue"))
# Import the necessary library for grid.arrange
library(gridExtra)
# Arrange box plots in a grid without legends
grid.arrange(ggplotGrob(plot_LassoSIR_boxplot + theme(legend.position = "none")),
             ggplotGrob(plot_wlsSIR_boxplot + theme(legend.position = "none")),
             ggplotGrob(plot_sir_boxplot + theme(legend.position = "none")) , ncol = 3)


## boxplot pour PC1 : 
# Convert scatter plots to box plots
plot_LassoSIR_boxplot <- ggplot(pc.lassoSIR_df, aes(x = factor(class_labels), y = pc.lassoSIR[, 1], fill = factor(class_labels))) +
  geom_boxplot() +
  labs(y = TeX('$PC_1$'), x = TeX('Factors')) +
  ggtitle("Lasso-SIR") +
  theme_minimal() +
  scale_fill_manual(values = c("red", "green", "blue"))
plot_wlsSIR_boxplot <- ggplot(pc.wls.sir_df, aes(x = factor(class_labels), y = pc.wls.sir[, 1], fill = factor(class_labels))) +
  geom_boxplot() +
  labs(y = TeX('$PC_1$'), x = TeX('Factors')) +
  ggtitle("SIR-WLS") +
  theme_minimal() +
  scale_fill_manual(values = c("red", "green", "blue"))
plot_sir_boxplot <- ggplot(pc.sir_df, aes(x = factor(class_labels), y = pc.sir[, 1], fill = factor(class_labels))) +
  geom_boxplot() +
  labs(y = TeX('$PC_1$'), x = TeX('Factors')) +
   ggtitle("SIR") +
  theme_minimal() +
  scale_fill_manual(values = c("red", "green", "blue"))
# Import the necessary library for grid.arrange
library(gridExtra)
# Arrange box plots in a grid without legends
grid.arrange(ggplotGrob(plot_LassoSIR_boxplot + theme(legend.position = "none")),
             ggplotGrob(plot_wlsSIR_boxplot + theme(legend.position = "none")),
             ggplotGrob(plot_sir_boxplot + theme(legend.position = "none")) , ncol = 3)

###logistic 

lasso.cv<-cv.glmnet(x= wine[,2:14], y=wine[,1],alpha =1, family="multinomial", nfolds=10)   ## cross validation to obtain optimal mu
best_lam<-lasso.cv$lambda.min
fit<-glmnet(x=wine[,2:14], y=wine[,1], alpha =1, family="multinomial", lambda = best_lam)
beta <- fit$beta

####


pred <- prediction(ROCR.simple$predictions,ROCR.simple$labels)




#### riboflavin ---------
## import data
#library(qut) 
#data("riboflavin", package = "qut")
#x.riboflavin<- as.matrix(riboflavin[["x"]])
#y.riboflavin<-as.vector(riboflavin[["y"]])

##@ or  from hdi 
data("riboflavin", package = "hdi")
data(riboflavin) 
x.riboflavin <- riboflavin$x
y.riboflavin <- riboflavin$y
#y.riboflavin <- exp(riboflavin$y) #@ going back to linear case 

#set.seed(37)
#fit.stab <- hdi(x, y, method = "stability", B = 500, EV = 1, q = 20)
#fit.stab

## choisire le pourcentage de l'chontillon test selon  les donnés 
smp_size <- floor(0.60 * nrow(x.riboflavin))
set.seed(36)
#set.seed(123)

#### Deviser les données en deux echontillons  "training set (donnée d'entraînement)" et "test set"  
train_ind <- sample(seq_len(nrow(x.riboflavin)), size = smp_size)
x.train <- as.matrix(x.riboflavin[train_ind, ])
x.test <- as.matrix(x.riboflavin[-train_ind, ] )
y.train<-as.vector(y.riboflavin[train_ind])
y.test<-as.vector(y.riboflavin[-train_ind])

#####glmnet
### as gived in  paper Buhlmann of PeterBühlmann 
beta.lasso<-normalize_vector(as.vector(Predict(x=x.train,y=y.train,method="Lasso")))
# Fit an Elastic Net model (alpha = 0.5 for equal combination of L1 and L2 penalties)
#beta.elasticnet<- normalize_vector(as.vector(Predict(x=x.train,y=y.train,method="Elastic_net")))
H=10
### my predict function 
sir.lasso <-LassoSIR( x.train,y.train, H, no.dim=1,solution.path=FALSE, categorical=FALSE, nfolds=10,screening=FALSE)
beta.lassosir <- normalize_vector(as.vector(sir.lasso$beta))
wls<-wls.sir(x.train,y.train)
beta.wls<-normalize_vector(as.vector(wls$betahat))
#### selected variables :-------
which(beta.lasso!= 0)
which(beta.lassosir!= 0)
which(beta.wls!= 0)
#which(beta.elasticnet!= 0)
## lambda such that the cv-error is within 1 standard error of the minimum
b.lasso <-as.matrix(beta.lasso); b.lassosir <-as.matrix(beta.lassosir) ; b.wls <-as.matrix(beta.wls)  ;# b.elasticnet <-as.matrix(beta.elasticnet) 
rownames(b.lasso)=rownames(b.lassosir)=rownames(b.wls)<-colnames(x.train)
rownames(b.lasso)[b.lasso != 0]
rownames(b.lassosir)[b.lassosir != 0]
rownames(b.wls)[b.wls != 0]
#rownames(b.elasticnet)[b.elasticnet!=0]
select.lasso <-which(beta.lasso!= 0)
select.lassosir<- which(beta.lassosir!= 0)
select.wls<-which(beta.wls!= 0)
#select<-list(select.lm,select.lasso,select.lassosir,select.wls)
common_elements1<- Reduce(intersect,list(select.lasso,select.lassosir))
common_elements2 <- Reduce(intersect,list(select.lasso,select.wls))
# Display common elements
common_elements1
common_elements2

rownames(b.lasso)[common_elements1]
rownames(b.lasso)[common_elements2]

#dev.off()
# Create a ggplot for  X.train hat(beta) by beta train data  ----
# Load the necessary library
library(ggplot2)

# Create the scatter plot with a regression line
plot_LASSO <- ggplot() +
  geom_point(aes(x = x.train %*% beta.lasso, y = y.train), color = "blue", size = 2) +
  geom_smooth(aes(x = x.train %*% beta.lasso, y = y.train), method = "lm", color = "red", linetype = "dashed") +
  labs(x = TeX('$X_{train} \\hat{\\beta}_{Lasso}$'), y = TeX('$Y_{train}$')) +
  theme_minimal()

# Print the plot
#print(plot_LASSO)

plot_LASSOSIR <-ggplot() +
  geom_point(aes(x =x.train%*%beta.lassosir , y = y.train), color = "orange", size = 2) +
  geom_smooth(aes(x = x.train %*% beta.lassosir, y = y.train), method = "lm", color = "red", linetype = "dashed") +
  labs( x= TeX('$X_{train} \\hat{\\beta}_{Lasso-SIR} $'), y= TeX('$Y_{train}$'))+
  theme_minimal()

plot_wlsSIR<-ggplot() +
  geom_point(aes(x = x.train%*%beta.wls , y = y.train), color = "green", size = 2) +
 geom_smooth(aes(x = x.train %*% beta.wls, y = y.train), method = "lm", color = "red", linetype = "dashed") +
  labs( x= TeX('$X_{train} \\hat{\\beta}_{SIR-WLS} $'), y= TeX('$Y_{train}$'))+
  theme_minimal()

grid.arrange( plot_LASSO, plot_LASSOSIR,plot_wlsSIR, ncol = 3)

### test.data --- #log(y.test)
plot_LASSO <-ggplot() +
  geom_point(aes(x = x.test %*% beta.lasso, y = y.test), color = "blue", size = 2) +
  geom_smooth(aes(x = x.test %*% beta.lasso, y = y.test), method = "lm", color = "red", linetype = "dashed") +
  labs(x = TeX('$X_{test} \\hat{\\beta}_{Lasso}$'), y = TeX('$Y_{test}$')) +
  theme_minimal()

plot_LASSOSIR <-ggplot() +
  geom_point(aes(x =x.test%*%beta.lassosir , y = y.test), color = "orange", size = 2) +
  geom_smooth(aes(x =x.test %*% beta.lassosir, y = y.test), method = "lm", color = "red", linetype = "dashed") +
  labs( x= TeX('$X_{test} \\hat{\\beta}_{Lasso-SIR} $'), y= TeX('$Y_{test}$'))+
  theme_minimal()

plot_wlsSIR<-ggplot() +
  geom_point(aes(x = x.test%*%beta.wls , y = y.test), color = "green", size = 2) +
 geom_smooth(aes(x = x.test %*% beta.wls, y = y.test), method = "lm", color = "red", linetype = "dashed") +
  labs( x= TeX('$X_{test} \\hat{\\beta}_{SIR-WLS} $'), y= TeX('$Y_{test}$'))+
  theme_minimal()

grid.arrange( plot_LASSO, plot_LASSOSIR,plot_wlsSIR, ncol = 3)
library(Metrics)
# Compute MSE for each method
mse_lasso <- mse(y.test , x.test %*% beta.lasso) ; rmse_lasso <- rmse(y.test , x.test %*% beta.lasso)
mse_LassoSIR <- mse(y.test ,x.test %*% beta.lassosir);rmse_LassoSIR <- rmse(y.test ,x.test %*% beta.lassosir)
mse_proposed <- mse(y.test ,x.test %*% beta.wls) ;rmse_proposed <- rmse(y.test ,x.test %*% beta.wls)

# Display MSE values
print(paste("MSE RMSE for Lasso:", mse_lasso,rmse_lasso ))
print(paste("MSE RMSE for LassoSIR:", mse_LassoSIR,rmse_LassoSIR))
print(paste("MSE for Proposed Method:", mse_proposed,rmse_proposed))
#grid.arrange( plot_LASSO, plot_LASSOSIR,plot_wlsSIR, plot_SIRNCT,ncol = 2,nrow=2)





#### ggplot for X.test hatbeta 
##### result obteined for wls 
which(beta.lasso!= 0)
[1]   34   73  189  244  282  431  589  624  786  897  974 1100 1282 1425 1478 1503 1524 1762 1802 1810 1827 1859 2027 2117
[25] 2242 2443 2462 2474 2529 2564 2874 2923 3001 3276 3379 3491 3514 3539 3938 3950 4003
 which(beta.lassosir!= 0)
[1]  329  403  485 1462 3248
 which(beta.wls!= 0)
[1]   34   44   73  270  283  337  415  490  491  585  603  604  681  694  695  710  711  712  733 1069 1290 1291 1478 1501
[25] 1502 1503 1509 1510 1511 1923 2032 2055 2095 2613 2642 2923 2928 2929 3038 3153 3206 3239 3321 3408 3667 3668 3669 3670
[49] 4002 4003 4004 4005 4006

 common_elements1 #( lasso et lassoSIR )
integer(0)
 common_elements2 #( lasso et WLS )
[1]   34   73 1478 1503 2923 4003



rownames(b.wls)[common_elements2]
rownames(b.lassosir)[common_elements1]



stable slecetioned gene in other paper : LYSC_at, YOAB_at, and YXLD_at:
  Les gene selectioné 

 rownames(b.wls)[common_elements2]
[1] "ALD_at"  "ARGF_at" "YCDH_at" "YCGO_at" "YPTA_at" "YXLD_at"

by spars SIR  Hilafu1_ and Sandra E. Safo  :
  XHLA_at, YCGO_at, YHDX_r_at, YRZI_r_at, YTGD_at, YCKE_at, YXLD_at, YCDH_at, GAPB_at

> rownames(b.wls)[b.wls != 0]
[1] "ALD_at"      "AMYC_at"     "ARGF_at"     "DACA_at"     "DEGQ_r_at"   "ETFB_at"     "GAPB_at"    
[8] "GSIB_at"     "GSPA_at"     "LACA_at"     "LICA_at"     "LICB_at"     "MRGA_at"     "MTLA_at"    
[15] "MTLD_at"     "NADA_at"     "NADB_at"     "NADC_at"     "NDK_at"      "SIGY_at"     "XKDK_at"    
[22] "XKDM_at"     "YCDH_at"     "YCGM_at"     "YCGN_at"     "YCGO_at"     "YCIA_at"     "YCIB_at"    
[29] "YCIC_at"     "YFMH_r_at"   "YHDX_r_at"   "YHFH_r_at"   "YHZA_at"     "YOCH_at"     "YOEB_at"    
[36] "YPTA_at"     "YPUF_at"     "YPUG_at"     "YQGN-P_i_at" "YRBA_at"     "YRPE_at"     "YRZI_r_at"  
[43] "YTIA_at"     "YTZD_at"     "YVFK_at"     "YVFL_at"     "YVFM_at"     "YVFO_at"     "YXLC_at"    
[50] "YXLD_at"     "YXLE_at"     "YXLF_at"     "YXLG_at" 





##### ionosphre----- c'est pas bon -----
ionosphere <- read.csv("~/Downloads/ionosphere/ionosphere.data", header = FALSE)
x.ionosphere<-as.matrix(ionosphere[,1:34]) ;
 y.ionosphere<-as.vector(y.ionosphere <- ifelse(ionosphere[,35] == "g", 1, 0)) 
dim(x.ionosphere)


smp_size <- floor(0.60 * nrow(x.ionosphere))
set.seed(136)

#### Deviser les données en deux echontillons  "training set (donnée d'entraînement)" et "test set"  
train_ind <- sample(seq_len(nrow(x.ionosphere)), size = smp_size)
x.ionosphere.train <- as.matrix(x.ionosphere[train_ind, ])
x.ionosphere.test <- as.matrix(x.ionosphere[-train_ind, ] )
y.ionosphere.train<-as.vector(y.ionosphere[train_ind])
y.ionosphere.test<-as.vector(y.ionosphere[-train_ind])


#sir.lassosir <- LassoSIR(x.ionosphere,y.ionosphere,solution.path=FALSE, categorical=TRUE, nfolds=10,screening=FALSE,no.dim=1)
sir.lassosir <- LassoSIR(cbind(1, x.ionosphere.train),y.ionosphere.train,solution.path=FALSE, categorical=TRUE, nfolds=10,screening=FALSE,no.dim=1)

beta.lassosir<-sir.lassosir$beta
beta.lassosir<-normalize_matrix(beta.lassosir)
#lasso.sir.slct<-which(beta.lassosir!=0)
#wls<-wls.sir(x.ionosphere, y.ionosphere, categorical=TRUE,ndim=1)
wls<-wls.sir(cbind(1, x.ionosphere.train), y.ionosphere.train, categorical=TRUE,ndim=1)

#beta.wls<-as.vector(wls$betahat)
beta.wls<-wls$betahat
beta.wls<-normalize_matrix(beta.wls)
#beta.wls
#wls.slct<-wls$select # wls.slct 

#out.sir<-do.sir(x.ionosphere,y.ionosphere, ndim = 1)#, categorical=TRUE,no.dim=2)
out.sir<-do.sir(cbind(1, x.ionosphere.train),y.ionosphere.train, ndim = 1)# categorical=TRUE,no.dim=2)

beta.sir<-normalize_matrix(out.sir$projection)


##### for lassoSIR
linear_predictor_lassoSIR <- cbind(1, x.ionosphere.test)%*% beta.lassosir 
#linear_predictor_lassoSIR <- x.ionosphere %*% beta.lassosir 
predicted_probabilities_lassoSIR<- as.vector(1 / (1 + exp(-linear_predictor_lassoSIR)))
predicted_lassoSIR <- ifelse(predicted_probabilities_lassoSIR > 0.5, 1, 0)
# Plot ROC curve
library(pROC)
roc_curve_LassoSIR <- roc(y.ionosphere.test,predicted_probabilities_lassoSIR, plot=TRUE, legacy.axes=TRUE)
#plot(roc_curve_LassoSIR, main = "ROC Curve for LassoSIR Logistic Regression", col = "blue")
### for wls_SIR 
linear_predictor_wls <- cbind(1, x.ionosphere.test) %*% beta.wls 
#linear_predictor_wls <- x.ionosphere %*% beta.wls 
predicted_probabilities_wls<- as.vector(1 / (1 + exp(-linear_predictor_wls)))
predicted_wls <- ifelse(predicted_probabilities_wls > 0.5, 1, 0)

# Plot ROC curve
roc_curve_wls <- roc(y.ionosphere.test,predicted_probabilities_wls, plot=TRUE, legacy.axes=TRUE)
plot(roc_curve_wls, main = "ROC Curve for wlsSIR Logistic Regression", col = "red")
lines(roc_curve_LassoSIR, col = "blue")
### for lasso 
lasso.cv<-cv.glmnet( x=cbind(1, x.ionosphere.train),y=y.ionosphere.train,alpha =1, family="binomial", nfolds=10, intercept=TRUE)   ## cross validation to obtain optimal mu
best_lam<-lasso.cv$lambda.min
fit<-glmnet( x=cbind(1, x.ionosphere.train),y=y.ionosphere.train, alpha =1, family="binomial", lambda = best_lam,itercept=TRUE)
beta.lasso <- as.vector(fit$beta)
#beta.lasso <- as.vector(coef(fit))
linear_predictor_lasso <-cbind(1, x.ionosphere.test) %*% beta.lasso 
#linear_predictor_lasso <- x.ionosphere %*% beta.lasso 
predicted_probabilities_lasso<- as.vector(1 / (1 + exp(-linear_predictor_lasso)))
predicted_lasso <- ifelse(predicted_probabilities_lasso > 0.5, 1, 0)

#Plot ROC curve
roc_curve_lasso <- roc(y.ionosphere.test,predicted_probabilities_lasso, plot=TRUE, legacy.axes=TRUE)
par(mfrow = c(1, 1))
plot(roc_curve_wls, main = "ROC Curve ", col = "red")
lines(roc_curve_LassoSIR, col = "blue")
lines(roc_curve_lasso, col = "green")
legend("bottomright", legend = c("ROC Curve_wls","ROC Curve_lassoSIR", "ROC Curve Lasso"), 
       col = c( "red","blue","green"), lwd = 2,cex = 0.6)
# Calculate AUC
auc_value_wls<- auc(roc_curve_wls); auc_value_wls
auc_value_lassSIR<-auc(roc_curve_LassoSIR);auc_value_lassSIR
#auc_value_lasso<-auc(roc_curve_lasso);auc_value_lasso



misclassification_rate_lassoSIR<- mean(y.ionosphere.test!= predicted_lassoSIR)
# Create confusion matrix
conf_matrix <- table(y.ionosphere.test, predicted_lassoSIR)
# Compute sensitivity (recall)
sensitivity <- conf_matrix[2, 2] / sum(conf_matrix[2, ])
# Compute specificity
specificity <- conf_matrix[1, 1] / sum(conf_matrix[1, ])

# Output results
print(paste("Misclassification Rate:for lassoSIR", round(misclassification_rate_lassoSIR, 4)))
print(paste("Sensitivity (Recall):", round(sensitivity, 4)))
print(paste("Specificity:", round(specificity, 4)))
misclassification_rate_wls<- mean(y.ionosphere.test != predicted_wls)
# Create confusion matrix
conf_matrix <- table(y.ionosphere.test, predicted_wls)
# Compute sensitivity (recall)
sensitivity <- conf_matrix[2, 2] / sum(conf_matrix[2, ])
# Compute specificity
specificity <- conf_matrix[1, 1] / sum(conf_matrix[1, ])

# Output results
print(paste("Misclassification Rate:for WLS_SIR", round(misclassification_rate_wls, 4)))
print(paste("Sensitivity (Recall) for wls :", round(sensitivity, 4)))
print(paste("Specificity:for wls ", round(specificity, 4)))
boxplot(misclassification_rate_lassoSIR,misclassification_rate_wls)

#set.seed(1234)

####### datamicroarray package -----
library(devtools)
install_github('ramhiser/datamicroarray')

describe_data()

####----- alon ----- good without intercept -----
data(Colon)
#library(datamicroarray)
#data('alon', package = 'datamicroarray')
x.alon<-Colon$X        #as.matrix(alon$x)
y.alon<-Colon$Y   #as.vector( ifelse(alon$y == "t", 1, 0)) 
colnames(x.alon)<-Colon$gene.names
# of 40 tumour and 22 normal ;  tumor=1 
sir.lassosir <- LassoSIR(cbind(1,x.alon),y.alon,solution.path=FALSE, categorical=TRUE, nfolds=10,screening=FALSE,no.dim=1)
#sir.lassosir <- LassoSIR(x.alon,y.alon,solution.path=FALSE, categorical=TRUE, nfolds=10,screening=FALSE,no.dim=1)
beta.lassosir<-sir.lassosir$beta
beta.lassosir<-normalize_matrix(beta.lassosir)
#lasso.sir.slct<-which(beta.lassosir!=0)
wls<-wls.sir(cbind(1,x.alon), y.alon, categorical=TRUE,ndim=1)
#wls<-wls.sir(x.alon, y.alon, categorical=TRUE,ndim=1)
#beta.wls<-as.vector(wls$betahat)
beta.wls<-Re(wls$betahat)
beta.wls<-normalize_matrix(beta.wls)


## lasso 
fit.cv<-cv.glmnet(x.alon,y.alon,alpha =1, family="binomial", nfolds=10, intercept=TRUE)   ## cross validation to obtain optimal mu
best_lam<-fit.cv$lambda.min
fit<-glmnet(x.alon,y.alon,alpha =1,  family="binomial", lambda = best_lam,itercept=TRUE)
beta.lasso <- as.vector(coef(fit))

#@@@@ compar lasso intercept true with lasso cbind(1,x)
fit.cv<-cv.glmnet(cbind(1,x.alon),y.alon,alpha =1, family="binomial", nfolds=10, intercept=FALSE)   ## cross validation to obtain optimal mu
best_lam<-fit.cv$lambda.min
fit<-glmnet(cbind(1,x.alon),y.alon,alpha =1,  family="binomial", lambda = best_lam,itercept=FALSE)
beta.lasso2 <- as.vector(fit$beta)
beta.lasso<- as.vector(fit$beta)

linear_predictor<-rep(1,length(y.alon)) ; predicted_probabilities<- as.vector(1 / (1 + exp(-linear_predictor)))
linear_predictor_lasso <- cbind(1,x.alon) %*% beta.lasso
#linear_predictor_lasso2<- cbind(1,x.alon) %*% beta.lasso2 
#linear_predictor_lasso <- x.alon %*% beta.lasso 
predicted_probabilities_lasso<- as.vector(1 / (1 + exp(-linear_predictor_lasso)))
##predicted_probabilities_lasso2<- as.vector(1 / (1 + exp(-linear_predictor_lasso2)))
#predicted_probabilities_lasso==predicted_probabilities_lasso2
predicted_lasso <- ifelse(predicted_probabilities_lasso > 0.5, 1, 0)
###

#Plot ROC curve
roc_curve_lasso <- roc(y.alon,predicted_probabilities_lasso, plot=TRUE, legacy.axes=TRUE)
roc_curve_lasso2 <- roc(y.alon,predicted_probabilities_lasso2, plot=TRUE, legacy.axes=TRUE)

par(mfrow = c(1, 1))
plot(roc_curve_lasso, main = "ROC Curve ", col = "red")
lines(roc_curve_lasso2, col = "blue")

##### for lassoSIR
linear_predictor_lassoSIR <- cbind(1,x.alon) %*% beta.lassosir 
#linear_predictor_lassoSIR <- x.alon %*% beta.lassosir
predicted_probabilities_lassoSIR<- as.vector(1 / (1 + exp(-linear_predictor_lassoSIR)))
predicted_lassoSIR <- ifelse(predicted_probabilities_lassoSIR > 0.5, 1, 0)
# Plot ROC curve
library(pROC)
roc_curve_LassoSIR <- roc(y.alon,predicted_probabilities_lassoSIR, plot=FALSE, legacy.axes=TRUE)
#plot(roc_curve_LassoSIR, main = "ROC Curve for LassoSIR Logistic Regression", col = "blue")
### for wls_SIR 
linear_predictor_wls <- cbind(1,x.alon) %*% beta.wls 
#linear_predictor_wls <- x.alon %*% beta.wls 
predicted_probabilities_wls<- as.vector(1 / (1 + exp(-linear_predictor_wls)))
predicted_wls <- ifelse(predicted_probabilities_wls > 0.5, 1, 0)

# Plot ROC curve
roc_curve_wls <- roc(y.alon,predicted_probabilities_wls, plot=TRUE, legacy.axes=TRUE)
plot(roc_curve_wls, main = "ROC Curve for wlsSIR Logistic Regression", col = "red")
lines(roc_curve_LassoSIR, col = "blue")
### for lasso 
linear_predictor_lasso <- cbind(1,x.alon) %*% beta.lasso 
#linear_predictor_lasso <- x.alon %*% beta.lasso 
predicted_probabilities_lasso<- as.vector(1 / (1 + exp(-linear_predictor_lasso)))
predicted_lasso <- ifelse(predicted_probabilities_lasso > 0.5, 1, 0)
#predicted_lasso2 <- ifelse(predicted_probabilities_lasso2 > 0.5, 1, 0)
#Plot ROC curve
roc_curve_lasso <- roc(y.alon,predicted_probabilities_lasso, plot=TRUE, legacy.axes=TRUE)
par(mfrow = c(1, 1))
plot(roc_curve_wls, main = "ROC Curve ", col = "red")
lines(roc_curve_LassoSIR, col = "blue")
lines(roc_curve_lasso, col = "green")
legend("bottomright", legend = c("SIR-WLS","Lasso-SIR", "Lasso"), 
       col = c( "red","blue","green"), lwd = 2,cex = 0.6)
# Calculate AUC
auc_value_wls<- auc(roc_curve_wls); auc_value_wls
auc_value_lassSIR<-auc(roc_curve_LassoSIR);auc_value_lassSIR
auc_value_lasso<-auc(roc_curve_lasso);auc_value_lasso





####
misclassification_rate_lassoSIR<- mean(y.alon != predicted_lassoSIR)
# Create confusion matrix
conf_matrix <- table(y.alon, predicted_lassoSIR)
# Compute sensitivity (recall)
sensitivity <- conf_matrix[2, 2] / sum(conf_matrix[2, ])
# Compute specificity
specificity <- conf_matrix[1, 1] / sum(conf_matrix[1, ])



boxplot(misclassification_rate_lassoSIR,misclassification_rate_wls)
# for wls 
misclassification_rate_wls<- mean(y.alon != predicted_wls)
# Create confusion matrix
conf_matrix_wls<- table(y.alon, predicted_wls)
# Compute sensitivity (recall)
sensitivity_wls<- conf_matrix_wls[2, 2] / sum(conf_matrix_wls[2, ])
# Compute specificity
specificity_wls <- conf_matrix_wls[1, 1] / sum(conf_matrix_wls[1, ])

# for lasso 
misclassification_rate_lasso<- mean(y.alon != predicted_lasso)
# Create confusion matrix
conf_matrix_lasso<- table(y.alon, predicted_lasso)
# Compute sensitivity (recall)
sensitivity_lasso <- conf_matrix_lasso[2, 2] / sum(conf_matrix_lasso[2, ])
# Compute specificity
specificity <- conf_matrix_lasso[1, 1] / sum(conf_matrix_lasso[1, ])
# Output results
print(paste("Misclassification Rate:for lassoSIR", round(misclassification_rate_lassoSIR, 4)))
print(paste("Misclassification Rate:for wls", round(misclassification_rate_wls, 4)))
print(paste("Misclassification Rate:for lasso", round(misclassification_rate_wls, 4)))
# Output results
print(paste("Sensitivity (Recall):for lassoSIR", round(sensitivity, 4)))
print(paste("Sensitivity (Recall) for wls :", round(sensitivity_wls, 4)))
print(paste("Sensitivity (Recall) for lasso :", round(sensitivity_lasso, 4)))
# Output results
print(paste("Specificity:for lassoSIR", round(specificity, 4)))
print(paste("Specificity:for wls", round(specificity, 4)))
print(paste("Specificity:for lasso ", round(specificity, 4)))
boxplot(misclassification_rate_lassoSIR,misclassification_rate_wls,misclassification_rate_lasso)



###--- chin data -------- good with or without intercept ------
#library(devtools)
#install_github('ramhiser/datamicroarray')
 
library(datamicroarray)
data('chin', package = 'datamicroarray');
 x.chin<-as.matrix(chin$x); dim(x.chin) #scale(as.matrix(chin$x))
 y.chin<-as.vector( ifelse(chin$y == "positive", 1, 0)) 
 smp_size <- floor(1 * nrow(x.chin))
 # Set the seed for reproducibility
 set.seed(123)
 sir.lassosir <- LassoSIR(cbind(1,x.chin),y.chin,solution.path=FALSE, categorical=TRUE, nfolds=10,screening=TRUE,no.dim=1)
 #sir.lassosir <- LassoSIR(x.chin.train,y.chin.train,solution.path=FALSE, categorical=TRUE, nfolds=10,screening=TRUE,no.dim=1)
 beta.lassosir<-sir.lassosir$beta
 beta.lassosir<-normalize_matrix(beta.lassosir)

 
 wls<-wls.sir(cbind(1,x.chin), y.chin, categorical=TRUE,ndim=1)
 #wls<-wls.sir(x.chin.train, y.chin.train, categorical=TRUE,ndim=1)
 #beta.wls<-as.vector(wls$betahat)
 beta.wls<-wls$betahat
 beta.wls<-normalize_matrix(beta.wls)
 
 
 ## lasso 
 fit.cv<-cv.glmnet(x.chin,y.chin,alpha =1, family="binomial", nfolds=10, intercept=FALSE)   ## cross validation to obtain optimal mu
 best_lam<-fit.cv$lambda.min
 fit<-glmnet(x.chin,y.chin,alpha =1,  family="binomial", lambda = best_lam,itercept=FALSE)
 beta.lasso <- as.vector(coef(fit))
 #linear_predictor_lasso <- cbind(1,x.chin.train) %*% beta.lasso
 linear_predictor_lasso <- cbind(1,x.chin)%*% beta.lasso
  #linear_predictor_lasso <- x.chin %*% beta.lasso 
 predicted_probabilities_lasso<- as.vector(1 / (1 + exp(-linear_predictor_lasso)))
 #predicted_probabilities_lasso2<- as.vector(1 / (1 + exp(-linear_predictor_lasso2)))
 #predicted_probabilities_lasso==predicted_probabilities_lasso2
 predicted_lasso <- ifelse(predicted_probabilities_lasso > 0.5, 1, 0)
 ### #Plot ROC curve
 roc_curve_lasso <- roc(y.chin,predicted_probabilities_lasso, plot=TRUE, legacy.axes=TRUE)
 ##### for lassoSIR
 linear_predictor_lassoSIR <- cbind(1,x.chin) %*% beta.lassosir 
 #linear_predictor_lassoSIR <- x.chin %*% beta.lassosir
 predicted_probabilities_lassoSIR<- as.vector(1 / (1 + exp(-linear_predictor_lassoSIR)))
 predicted_lassoSIR <- ifelse(predicted_probabilities_lassoSIR > 0.5, 1, 0)
 # Plot ROC curve
 roc_curve_LassoSIR <- roc(y.chin,predicted_probabilities_lassoSIR, plot=FALSE, legacy.axes=TRUE)
 plot(roc_curve_LassoSIR, main = "ROC Curve for LassoSIR Logistic Regression", col = "blue")
 ### for wls_SIR 
 linear_predictor_wls <- cbind(1,x.chin) %*% beta.wls 
 #linear_predictor_wls <- x.chin %*% beta.wls 
 predicted_probabilities_wls<- as.vector(1 / (1 + exp(-linear_predictor_wls)))
 predicted_wls <- ifelse(predicted_probabilities_wls > 0.5, 1, 0)
 
 # Plot ROC curve
 roc_curve_wls <- roc(y.chin,predicted_probabilities_wls, plot=TRUE, legacy.axes=TRUE)
 plot(roc_curve_wls, main = "ROC Curve for wlsSIR Logistic Regression", col = "red")
 lines(roc_curve_LassoSIR, col = "blue")
 
 
 #Plot ROC curve
 par(mfrow = c(1, 1))
 plot(roc_curve_wls, main = "ROC Curve ", col = "red")
 lines(roc_curve_LassoSIR, col = "blue")
 lines(roc_curve_lasso, col = "green")
 legend("bottomright", legend = c("ROC Curve_wls","ROC Curve_lassoSIR", "ROC Curve Lasso"), 
        col = c( "red","blue","green"), lwd = 2,cex = 0.6)
 
 # Create a data frame for the ROC curve data
 roc_data_wls <- data.frame(
   Sensitivity = roc_curve_wls$sensitivities,
   Specificity = 1 - roc_curve_wls$specificities
 )
 
 # Create a data frame for the LassoSIR ROC curve data
 roc_data_LassoSIR <- data.frame(
   Sensitivity = roc_curve_LassoSIR$sensitivities,
   Specificity = 1 - roc_curve_LassoSIR$specificities
 )
 
 # Create a data frame for the Lasso ROC curve data
 roc_data_lasso <- data.frame(
   Sensitivity = roc_curve_lasso$sensitivities,
   Specificity = 1 - roc_curve_lasso$specificities
 )
 
 
 ggplot() +
   geom_line(data = roc_data_wls, aes(x = Specificity, y = Sensitivity,  color = "SIR-WLS")) +
   geom_line(data = roc_data_LassoSIR, aes(x = Specificity, y = Sensitivity,  color = "Lasso-SIR")) +
   geom_line(data = roc_data_lasso, aes(x = Specificity, y = Sensitivity,  color = "Lasso")) +
   labs(title = "ROC Curve", x = "Specificity (1 - False Positive Rate)", y = "Sensitivity (True Positive Rate)",
        color = "Method") +
   scale_color_manual(values = c("SIR-WLS" = "blue", "Lasso-SIR" = "red", "Lasso" = "green")) +
   theme_minimal()
 
 # Calculate AUC
 auc_value_wls<- auc(roc_curve_wls); auc_value_wls
 auc_value_lassSIR<-auc(roc_curve_LassoSIR);auc_value_lassSIR
 auc_value_lasso<-auc(roc_curve_lasso);auc_value_lasso
 
 
 
 
 library(ggplot2)
 
 # Create a data frame for the ROC curve data
 roc_data_wls <- data.frame(
   FPR = roc_curve_wls$specificities,
   TPR = roc_curve_wls$sensitivities
 )
 
 # Create a data frame for the LassoSIR ROC curve data
 roc_data_LassoSIR <- data.frame(
   FPR = roc_curve_LassoSIR$specificities,
   TPR = roc_curve_LassoSIR$sensitivities
 )
 
 # Create a data frame for the Lasso ROC curve data
 roc_data_lasso <- data.frame(
   FPR = roc_curve_lasso$specificities,
   TPR = roc_curve_lasso$sensitivities
 )
 
 # Plot the ROC curves using ggplot
 ggplot() +
   geom_line(data = roc_data_wls, aes(x = 1-FPR, y = TPR), color = "red") +
   geom_line(data = roc_data_LassoSIR, aes(x = 1-FPR, y = TPR), color = "blue") +
   geom_line(data = roc_data_lasso, aes(x = 1-FPR, y = TPR), color = "green") +
   labs(title = "ROC Curve", x = "False Positive Rate", y = "True Positive Rate") +
   scale_color_manual(values = c("red", "blue", "green"), 
                      labels = c("SIR-WLS", "Lasso-SIR", "Lasso")) +
   theme_minimal() +
   theme(legend.position = "bottom", legend.title = element_blank())
 #### selected variables :
 which(beta.lasso!= 0)
 which(beta.lassosir!= 0)
 which(beta.wls!= 0)
 #which(beta.elasticnet!= 0)
 ## lambda such that the cv-error is within 1 standard error of the minimum
 b.lasso <-as.matrix(beta.lasso); b.lassosir <-as.matrix(beta.lassosir) ; b.wls <-as.matrix(beta.wls)  
 #rownames(b.elasticnet)[b.elasticnet!=0]
 select.lasso <-which(beta.lasso!= 0)
 select.lassosir<- which(beta.lassosir!= 0)
 select.wls<-which(beta.wls!= 0)
 #select<-list(select.lm,select.lasso,select.lassosir,select.wls)
 common_elements1<- Reduce(intersect,list(select.lasso,select.lassosir))
 common_elements2 <- Reduce(intersect,list(select.lasso,select.wls))
 common_elements3 <- Reduce(intersect,list(select.lassosir,select.wls))
 
 # Display common elements
 common_elements1
 common_elements2
 common_elements3
plot(y.chin,x.chin%*%b.lasso)
 
  
 ###@. Leukemia  golub data ------
data(leukemia)
 #data('golub', package = 'datamicroarray')
 x.golub<-leukemia$X         #as.matrix(golub$x)
 y.golub<- leukemia$Y               #as.vector( ifelse(golub$y == "ALL", 1, 0)) 
 smp_size <- floor(0.50* nrow(x.golub))
 # Set the seed for reproducibility
 set.seed(123)
 #### Deviser les données en deux echontillons  "training set (donnée d'entraînement)" et "test set"  
 train_ind <- sample(seq_len(nrow(x.golub)), size = smp_size)
 x.golub.train <- as.matrix(x.golub[train_ind, ])
 x.golub.test <- as.matrix(x.golub[-train_ind, ] )
 y.golub.train<-as.vector(y.golub[train_ind])
 y.golub.test<-as.vector(y.golub[-train_ind])
 dim(x.golub.train); dim(x.golub)
 #sir.lassosir <- LassoSIR(x.golub,y.golub,solution.path=FALSE, categorical=TRUE, nfolds=10,screening=TRUE,no.dim =1)
 sir.lassosir <- LassoSIR( cbind(1,x.golub.train),y.golub.train,solution.path=FALSE, categorical=TRUE, nfolds=10,screening=TRUE,no.dim =1)
  beta.lassosir<-sir.lassosir$beta
 beta.lassosir<-normalize_matrix(beta.lassosir)
 #lasso.sir.slct<-which(beta.lassosir!=0)
 pc.lassoSIR<- cbind(1,x.golub.train)%*%beta.lassosir
 wls<-wls.sir(x= cbind(1,x.golub.train),y=y.golub.train, categorical=TRUE,ndim=1)

 #wls<-wls.sir(x=x.golub,y=y.golub, categorical=TRUE,ndim=1)
 #beta.wls<-as.vector(wls$betahat)
 beta.wls<-wls$betahat
 beta.wls<-normalize_matrix(beta.wls)
 #beta.wls
 #wls.slct<-wls$select # wls.slct 
 pc.wls.sir<- cbind(1,x.golub.train)%*%beta.wls
 ## lasso 
 fit.cv<-cv.glmnet(x.golub.train,y.golub.train,alpha =1, family="binomial", nfolds=10, intercept=TRUE)   ## cross validation to obtain optimal mu
 best_lam<-fit.cv$lambda.min
 fit<-glmnet(x.golub.train,y.golub.train,alpha =1,  family="binomial", lambda = best_lam,itercept=FALSE)
 beta.lasso <- as.vector(coef(fit))
# beta.lasso <- as.vector(fit$beta)
 
 # 
 #linear_predictor_lasso <- cbind(1,x.golub) %*% beta.lasso
 linear_predictor_lasso <- cbind(1,x.golub.test) %*% beta.lasso
 
 #linear_predictor_lasso <- x.golub %*% beta.lasso 
 predicted_probabilities_lasso<- as.vector(1 / (1 + exp(-linear_predictor_lasso)))
 #predicted_probabilities_lasso2<- as.vector(1 / (1 + exp(-linear_predictor_lasso2)))
 #predicted_probabilities_lasso==predicted_probabilities_lasso2
 predicted_lasso <- ifelse(predicted_probabilities_lasso > 0.5, 1, 0)
 ###
 
 #Plot ROC curve
 roc_curve_lasso <- roc(y.golub.test,predicted_probabilities_lasso, plot=TRUE, legacy.axes=TRUE)
 ##### for lassoSIR
 linear_predictor_lassoSIR <- cbind(1,x.golub.test) %*% beta.lassosir 
 #linear_predictor_lassoSIR <- x.golub %*% beta.lassosir
 predicted_probabilities_lassoSIR<- as.vector(1 / (1 + exp(-linear_predictor_lassoSIR)))
 predicted_lassoSIR <- ifelse(predicted_probabilities_lassoSIR > 0.5, 1, 0)
 # Plot ROC curve
 roc_curve_LassoSIR <- roc(y.golub.test,predicted_probabilities_lassoSIR, plot=FALSE, legacy.axes=TRUE)
 plot(roc_curve_LassoSIR, main = "ROC Curve for LassoSIR Logistic Regression", col = "blue")
 ### for wls_SIR 
 linear_predictor_wls <- cbind(1,x.golub.test) %*% beta.wls 
 #linear_predictor_wls <- x.golub%*% beta.wls 
 predicted_probabilities_wls<- as.vector(1 / (1 + exp(-linear_predictor_wls)))
 predicted_wls <- ifelse(predicted_probabilities_wls > 0.5, 1, 0)
 
 # Plot ROC curve
 roc_curve_wls <- roc(y.golub.test,predicted_probabilities_wls, plot=TRUE, legacy.axes=TRUE)
 plot(roc_curve_wls, main = "ROC Curve for wlsSIR Logistic Regression", col = "red")
 lines(roc_curve_LassoSIR, col = "blue")
 
 
 #Plot ROC curve
 par(mfrow = c(1, 1))
 plot(roc_curve_wls, main = "ROC Curve ", col = "red")
 lines(roc_curve_LassoSIR, col = "blue")
 lines(roc_curve_lasso, col = "green")
 legend("bottomright", legend = c("ROC Curve_wls","ROC Curve_lassoSIR", "ROC Curve Lasso"), 
        col = c( "red","blue","green"), lwd = 2,cex = 0.6)
 
 # Create a data frame for the ROC curve data
 roc_data_wls <- data.frame(
   Sensitivity = roc_curve_wls$sensitivities,
   Specificity = 1 - roc_curve_wls$specificities
 )
 
 # Create a data frame for the LassoSIR ROC curve data
 roc_data_LassoSIR <- data.frame(
   Sensitivity = roc_curve_LassoSIR$sensitivities,
   Specificity = 1 - roc_curve_LassoSIR$specificities
 )
 
 # Create a data frame for the Lasso ROC curve data
 roc_data_lasso <- data.frame(
   Sensitivity = roc_curve_lasso$sensitivities,
   Specificity = 1 - roc_curve_lasso$specificities
 )
 
 ggplot() +
   geom_line(data = roc_data_wls, aes(x = Specificity, y = Sensitivity,  color = "SIR-WLS")) +
   geom_line(data = roc_data_LassoSIR, aes(x = Specificity, y = Sensitivity,  color = "Lasso-SIR")) +
   geom_line(data = roc_data_lasso, aes(x = Specificity, y = Sensitivity,  color = "Lasso")) +
   labs(title = "ROC Curve", x = "Specificity (1 - False Positive Rate)", y = "Sensitivity (True Positive Rate)",
        color = "Method") +
   scale_color_manual(values = c("SIR-WLS" = "blue", "Lasso-SIR" = "red", "Lasso" = "green")) +
   theme_minimal()
 
 # Calculate AUC
 auc_value_wls<- auc(roc_curve_wls); auc_value_wls
 auc_value_lassSIR<-auc(roc_curve_LassoSIR);auc_value_lassSIR
 auc_value_lasso<-auc(roc_curve_lasso);auc_value_lasso
 
 
 #### selected variables :
 which(beta.lasso!= 0)
 which(beta.lassosir!= 0)
 which(beta.wls!= 0)
 #which(beta.elasticnet!= 0)
 colnames(x.golub.train)<-leukemia$gene.names
 ## lambda such that the cv-error is within 1 standard error of the minimum
 b.lasso <-as.matrix(beta.lasso); b.lassosir <-as.matrix(beta.lassosir) ; b.wls <-as.matrix(beta.wls)  
 rownames(b.lasso)=rownames(b.lassosir)=rownames(b.wls)<-colnames( cbind(1,x.golub.train))
 rownames(b.lasso)[b.lasso != 0]
 rownames(b.lassosir)[b.lassosir != 0]
 rownames(b.wls)[b.wls != 0]
 #rownames(b.elasticnet)[b.elasticnet!=0]
 select.lasso <-which(beta.lasso!= 0)
 select.lassosir<- which(beta.lassosir!= 0)
 select.wls<-which(beta.wls!= 0)
 #select<-list(select.lm,select.lasso,select.lassosir,select.wls)
 common_elements1<- Reduce(intersect,list(select.lasso,select.lassosir))
 common_elements2 <- Reduce(intersect,list(select.lasso,select.wls))
 common_elements3 <- Reduce(intersect,list(select.lassosir,select.wls))
 
  # Display common elements
 common_elements1
 common_elements2
 
 rownames(b.lasso)[common_elements1]
 rownames(b.lasso)[common_elements2]
 rownames(b.lasso)[common_elements3]
 
 ####------ khan --- 4 class -------
 library(plsgenomics)
 data(SRBCT)
 #library(datamicroarray)
  #data('khan', package = 'datamicroarray')
 #x.khan<-as.matrix(khan$x)
 #y.khan<-as.vector(khan$y)
#y.class<-unique(y.khan); length(y.class) 
#class_mapping <- c("EWS" = 1, "RMS" = 2, "NB" = 3, "BL" = 4)
#y.khan<-as.vector(as.numeric(class_mapping[y.khan]))
 #y.alon<-as.vector( ifelse(alon$y == "t", 1, 0)) 
 #dim(x.khan)

 #smp_size <- floor(0.60 * nrow(x.khan))
 #set.seed(36)
#### Deviser les données en deux echontillons  "training set (donnée d'entraînement)" et "test set"  
 #train_ind <- sample(seq_len(nrow(x.khan)), size = smp_size)
 #x.khan.train <- as.matrix(x.khan[train_ind, ])
#x.khan.test <- as.matrix(x.khan[-train_ind, ] )
#y.khan.train<-as.vector(y.khan[train_ind])
#y.khan.test<-as.vector(y.khan[-train_ind])
 x.khan<-as.matrix(SRBCT$X);  y.khan<-as.vector(SRBCT$Y)
 x.khan.train <- as.matrix(x.khan[1:65, ])
 x.khan.test <- as.matrix(x.khan[66:83, ] )
 y.khan.train<-as.vector(y.khan[1:65])
 y.khan.test<-as.vector(y.khan[66:83])
#sir.lassosir <- LassoSIR( khan[,2:14],khan[,1],solution.path=FALSE, categorical=TRUE, nfolds=10,screening=FALSE,no.dim=1)
## ndim=2 
sir.lassosir <- LassoSIR( x.khan.train ,y.khan.train,solution.path=FALSE, categorical=TRUE, nfolds=10,no.dim=3)
beta.lassosir<-sir.lassosir$beta
beta.lassosir<-normalize_matrix(beta.lassosir)
#lasso.sir.slct<-which(beta.lassosir!=0)
pc.lassoSIR<-x.khan.test%*%beta.lassosir
#wls<-wls.sir(x= khan[,2:14], y=khan[,1], categorical=TRUE,ndim=1) # ndim=1
wls<-wls.sir(x= x.khan.train, y=y.khan.train, categorical=TRUE,ndim=3)
#beta.wls<-as.vector(wls$betahat)
beta.wls<-wls$betahat
beta.wls<-normalize_matrix(beta.wls)
pc.wls.sir<-x.khan.test%*%beta.wls
 
#install.packages("scatterplot3d")
library(scatterplot3d)


# Plot for pc.lassoSIR
class_labels <- y.khan.test
pc.lassoSIR_df <- as.data.frame(pc.lassoSIR)
pc.wls.sir_df <- as.data.frame(pc.wls.sir)
#pc.sir_df <- as.data.frame(pc.sir)

library(scatterplot3d)

# Assuming pc.lassoSIR and pc.wls.sir are matrices with three principal components each
# Assuming class_labels is a vector containing the class labels for each sample
# Create a 3D scatter plot for LASSOSIR
scatterplot3d(pc.lassoSIR[, 1], pc.lassoSIR[, 2], pc.lassoSIR[, 3],pch=20,xlab = TeX('$PC_1$'), ylab = TeX('$PC_2$'),zlab = TeX('$PC_3$'), 
              color = c("red", "green", "blue", "orange")[as.factor(class_labels)],
              main = "Lasso-SIR")

# Create a 3D scatter plot for SIR-WLS  
scatterplot3d(pc.wls.sir[, 1], pc.wls.sir[, 2], pc.wls.sir[, 3],pch=20, 
              color = c("red", "green", "blue", "orange")[as.factor(class_labels)], xlab = TeX('$PC_1$'), ylab = TeX('$PC_2$'),zlab = TeX('$PC_3$'),
              main = "SIR-WLS")
### For LASSOSIR-----
plot_ly(x = pc.lassoSIR[, 1], y = pc.lassoSIR[, 2], z = pc.lassoSIR[, 3], 
                         type = "scatter3d", mode = "markers",
                         marker = list(color = c("red", "green", "blue", "orange")[as.factor(class_labels)],
                                       size = 5),
                         text = paste("PC1:", pc.lassoSIR[, 1], "<br>PC2:", pc.lassoSIR[, 2], "<br>PC3:", pc.lassoSIR[, 3]),
                         hoverinfo = "text",
                         name = "Lasso-SIR") %>%
  layout(scene = list(xaxis = list(title = "PC1"), 
                      yaxis = list(title = "PC2"), 
                      zaxis = list(title = "PC3")),
         title = "Lasso-SIR")

# For SIR-WLS
plot_ly(x = pc.wls.sir[, 1], y = pc.wls.sir[, 2], z = pc.wls.sir[, 3], 
                        type = "scatter3d", mode = "markers",
                        marker = list(color = c("red", "green", "blue", "orange")[as.factor(class_labels)],
                                      size = 5),
                        text = paste("PC1:", pc.wls.sir[, 1], "<br>PC2:", pc.wls.sir[, 2], "<br>PC3:", pc.wls.sir[, 3]),
                        hoverinfo = "text",
                        name = "SIR-WLS") %>%
  layout(scene = list(xaxis = list(title = "PC1"), 
                      yaxis = list(title = "PC2"), 
                      zaxis = list(title = "PC3")),
         title = "SIR-WLS")
#####
dev.off()
# Load necessary libraries
library(ggplot2)

library(ggplot2)
library(gridExtra)

# Plot 1: PC1 vs PC2
plot_LassoSIR1 <- ggplot(data = pc.lassoSIR_df, aes(x = pc.lassoSIR[, 1], y = pc.lassoSIR[, 2], color = factor(class_labels))) +
  geom_point() +
  labs(x = expression(PC[1]), y = expression(PC[2])) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("red", "green", "blue", "orange")) +
  ggtitle("Lasso-SIR")

plot_wlsSIR1 <- ggplot(data = pc.wls.sir_df, aes(x = pc.wls.sir[, 1], y = pc.wls.sir[, 2], color = factor(class_labels))) +
  geom_point() +
  labs(x = expression(PC[1]), y = expression(PC[2])) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("red", "green", "blue", "orange")) +
  ggtitle("SIR-WLS")

# Plot 2: PC2 vs PC3
plot_LassoSIR2 <- ggplot(data = pc.lassoSIR_df, aes(x = pc.lassoSIR[, 2], y = pc.lassoSIR[, 3], color = factor(class_labels))) +
  geom_point() +
  labs(x = expression(PC[2]), y = expression(PC[3])) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("red", "green", "blue", "orange")) +
  ggtitle("Lasso-SIR")

plot_wlsSIR2 <- ggplot(data = pc.wls.sir_df, aes(x = pc.wls.sir[, 2], y = pc.wls.sir[, 3], color = factor(class_labels))) +
  geom_point() +
  labs(x = expression(PC[2]), y = expression(PC[3])) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("red", "green", "blue", "orange")) +
  ggtitle(" SIR-WLS")

# Plot 3: PC3 vs PC1
plot_LassoSIR3 <- ggplot(data = pc.lassoSIR_df, aes(x = pc.lassoSIR[, 3], y = pc.lassoSIR[, 1], color = factor(class_labels))) +
  geom_point() +
  labs(x = expression(PC[3]), y = expression(PC[1])) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("red", "green", "blue", "orange")) +
  ggtitle("Lasso-SIR")

plot_wlsSIR3 <- ggplot(data = pc.wls.sir_df, aes(x = pc.wls.sir[, 3], y = pc.wls.sir[, 1], color = factor(class_labels))) +
  geom_point() +
  labs(x = expression(PC[3]), y = expression(PC[1])) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("red", "green", "blue", "orange")) +
  ggtitle(" SIR-WLS")

# Arrange the plots
grid.arrange(plot_LassoSIR1, plot_wlsSIR1, ncol = 2)
grid.arrange(plot_LassoSIR2, plot_wlsSIR2, ncol = 2)
grid.arrange(plot_LassoSIR3, plot_wlsSIR3, ncol = 2)
#### selected variables :
which(beta.lassosir!= 0)
which(beta.wls!= 0)
#which(beta.elasticnet!= 0)
## lambda such that the cv-error is within 1 standard error of the minimum
 b.lassosir <-as.matrix(beta.lassosir) ; b.wls <-as.matrix(beta.wls)  
rownames(b.lassosir)=rownames(b.wls)<-colnames( cbind(1,x.khan.train))
rownames(b.lassosir)[b.lassosir != 0]
rownames(b.wls)[b.wls != 0]
#rownames(b.elasticnet)[b.elasticnet!=0]
select.lassosir<- which(beta.lassosir!= 0)
select.wls<-which(beta.wls!= 0)
#select<-list(select.lm,select.lasso,select.lassosir,select.wls)
common_elements3 <- Reduce(intersect,list(select.lassosir,select.wls))

# Display common elements

rownames(b.lassosir)[common_elements3]




####------ Yeoh --- 6 class -------
library(datamicroarray)

data('yeoh', package = 'datamicroarray')
x.yeoh<-as.matrix(yeoh$x)
y.yeoh<-as.vector(yeoh$y)
y.class<-unique(y.yeoh); length(y.class) 
class_mapping <- c("BCR" = 1, "E2A" = 2, "Hyperdip" = 3, "MLL" = 4, "T" =5, "TEL"=6    )
y.yeoh<-as.vector(as.numeric(class_mapping[y.yeoh]))
#y.alon<-as.vector( ifelse(alon$y == "t", 1, 0)) 
dim(x.yeoh)

smp_size <- floor(0.70 * nrow(x.yeoh))
set.seed(1)
#### Deviser les données en deux echontillons  "training set (donnée d'entraînement)" et "test set"  
train_ind <- sample(seq_len(nrow(x.yeoh)), size = smp_size)
x.yeoh.train <- as.matrix(x.yeoh[train_ind, ])
x.yeoh.test <- as.matrix(x.yeoh[-train_ind, ] )
y.yeoh.train<-as.vector(y.yeoh[train_ind])
y.yeoh.test<-as.vector(y.yeoh[-train_ind])

#sir.lassosir <- LassoSIR( yeoh[,2:14],yeoh[,1],solution.path=FALSE, categorical=TRUE, nfolds=10,screening=FALSE,no.dim=1)
## ndim=2 
sir.lassosir <- LassoSIR( x.yeoh.train ,y.yeoh.train,solution.path=FALSE, categorical=TRUE, nfolds=10,no.dim=2)
beta.lassosir<-sir.lassosir$beta
beta.lassosir<-normalize_matrix(beta.lassosir)
#lasso.sir.slct<-which(beta.lassosir!=0)
pc.lassoSIR<-x.yeoh.test%*%beta.lassosir
#wls<-wls.sir(x= yeoh[,2:14], y=yeoh[,1], categorical=TRUE,ndim=1) # ndim=1
wls<-wls.sir(x= x.yeoh.train, y=y.yeoh.train, categorical=TRUE,ndim=2)
#beta.wls<-as.vector(wls$betahat)
beta.wls<-wls$betahat
beta.wls<-normalize_matrix(beta.wls)
#beta.wls
#wls.slct<-wls$select # wls.slct 
pc.wls.sir<-x.yeoh.test%*%beta.wls
plot(pc.wls.sir[,1],pc.wls.sir[,2])
plot(pc.lassoSIR[,1],pc.lassoSIR[,2])

#out.sir<-do.sir( x.yeoh.train, y.yeoh.train)#, categorical=TRUE,no.dim=2)
#beta.sir<-normalize_matrix(out.sir$projection)
#pc.sir<-x.yeoh.test%*%beta.sir



# Plot for pc.lassoSIR
class_labels <- y.yeoh.test
pc.lassoSIR_df <- as.data.frame(pc.lassoSIR)
pc.wls.sir_df <- as.data.frame(pc.wls.sir)
#pc.sir_df <- as.data.frame(pc.sir)

plot_LassoSIR <- ggplot() +
  geom_point(data = pc.lassoSIR_df, aes(x = pc.lassoSIR[, 1], y = pc.lassoSIR[, 2], color = factor(class_labels)),show.legend = FALSE) +
  labs(x = TeX('$PC_1$'), y = TeX('$PC_2$')) +
  ggtitle("Lasso-SIR") +
  theme_minimal()
# Specifying colors for each class
plot_LassoSIR <- plot_LassoSIR + scale_color_manual(values = c( "red", "green","blue","orange","violet","brown"))

# Plot for pc.wls.sir
plot_wlsSIR <- ggplot() +
  geom_point(data = pc.wls.sir_df, aes(y= pc.wls.sir[, 1],x = pc.wls.sir[, 2], color = factor(class_labels)),show.legend = FALSE) +
  labs(x = TeX('$PC_1$'), y = TeX('$PC_2$')) +
  ggtitle("SIR-WLS") +
  theme_minimal()

# Specifying colors for each class
plot_wlsSIR <- plot_wlsSIR + scale_color_manual(values = c( "red", "green","blue","orange","violet","brown"))

library(gridExtra)

grid.arrange(ggplotGrob(plot_LassoSIR + theme(legend.position = "none")),
             ggplotGrob(plot_wlsSIR + theme(legend.position = "none")) , ncol = 2)


### http://ico2s.org/datasets/microarray.html -----
dlbcl_preprocessed <- read.table("~/Downloads/dlbcl_preprocessed.txt", quote="\"", comment.char="")

breast_preprocessed <-t( read.table("~/Downloads/breast_preprocessed.txt", quote="\"", comment.char=""))

#Breast cancer dataset (47293 genes, 128 samples)
dim(breast_preprocessed)

x.golub<-as.matrix(golub$x)
y.golub<-as.vector( ifelse(golub$y == "ALL", 1, 0)) 

smp_size <- floor(0.60 * nrow(x.golub))
# Set the seed for reproducibility
set.seed(123)
#### Deviser les données en deux echontillons  "training set (donnée d'entraînement)" et "test set"  
train_ind <- sample(seq_len(nrow(x.golub)), size = smp_size)
x.golub.train <- as.matrix(x.golub[train_ind, ])
x.golub.test <- as.matrix(x.golub[-train_ind, ] )
y.golub.train<-as.vector(y.golub[train_ind])
y.golub.test<-as.vector(y.golub[-train_ind])
dim(x.golub.train); dim(x.golub)
#sir.lassosir <- LassoSIR(x.golub,y.golub,solution.path=FALSE, categorical=TRUE, nfolds=10,screening=TRUE,no.dim =1)
sir.lassosir <- LassoSIR( cbind(1,x.golub.train),y.golub.train,solution.path=FALSE, categorical=TRUE, nfolds=10,screening=TRUE,no.dim =1)
beta.lassosir<-sir.lassosir$beta
beta.lassosir<-normalize_matrix(beta.lassosir)
#lasso.sir.slct<-which(beta.lassosir!=0)
pc.lassoSIR<- cbind(1,x.golub.train)%*%beta.lassosir
wls<-wls.sir(x= cbind(1,x.golub.train),y=y.golub.train, categorical=TRUE,ndim=1)

#wls<-wls.sir(x=x.golub,y=y.golub, categorical=TRUE,ndim=1)
#beta.wls<-as.vector(wls$betahat)
beta.wls<-wls$betahat
beta.wls<-normalize_matrix(beta.wls)
#beta.wls
#wls.slct<-wls$select # wls.slct 
pc.wls.sir<- cbind(1,x.golub.train)%*%beta.wls
## lasso 
fit.cv<-cv.glmnet(x.golub.train,y.golub.train,alpha =1, family="binomial", nfolds=10, intercept=TRUE)   ## cross validation to obtain optimal mu
best_lam<-fit.cv$lambda.min
fit<-glmnet(x.golub.train,y.golub.train,alpha =1,  family="binomial", lambda = best_lam,itercept=FALSE)
beta.lasso <- as.vector(coef(fit))
# beta.lasso <- as.vector(fit$beta)

# 
#linear_predictor_lasso <- cbind(1,x.golub) %*% beta.lasso
linear_predictor_lasso <- cbind(1,x.golub.test) %*% beta.lasso

#linear_predictor_lasso <- x.golub %*% beta.lasso 
predicted_probabilities_lasso<- as.vector(1 / (1 + exp(-linear_predictor_lasso)))
#predicted_probabilities_lasso2<- as.vector(1 / (1 + exp(-linear_predictor_lasso2)))
#predicted_probabilities_lasso==predicted_probabilities_lasso2
predicted_lasso <- ifelse(predicted_probabilities_lasso > 0.5, 1, 0)
###

#Plot ROC curve
roc_curve_lasso <- roc(y.golub.test,predicted_probabilities_lasso, plot=TRUE, legacy.axes=TRUE)
##### for lassoSIR
linear_predictor_lassoSIR <- cbind(1,x.golub.test) %*% beta.lassosir 
#linear_predictor_lassoSIR <- x.golub %*% beta.lassosir
predicted_probabilities_lassoSIR<- as.vector(1 / (1 + exp(-linear_predictor_lassoSIR)))
predicted_lassoSIR <- ifelse(predicted_probabilities_lassoSIR > 0.5, 1, 0)
# Plot ROC curve
roc_curve_LassoSIR <- roc(y.golub.test,predicted_probabilities_lassoSIR, plot=FALSE, legacy.axes=TRUE)
plot(roc_curve_LassoSIR, main = "ROC Curve for LassoSIR Logistic Regression", col = "blue")
### for wls_SIR 
linear_predictor_wls <- cbind(1,x.golub.test) %*% beta.wls 
#linear_predictor_wls <- x.golub%*% beta.wls 
predicted_probabilities_wls<- as.vector(1 / (1 + exp(-linear_predictor_wls)))
predicted_wls <- ifelse(predicted_probabilities_wls > 0.5, 1, 0)

# Plot ROC curve
roc_curve_wls <- roc(y.golub.test,predicted_probabilities_wls, plot=TRUE, legacy.axes=TRUE)
plot(roc_curve_wls, main = "ROC Curve for wlsSIR Logistic Regression", col = "red")
lines(roc_curve_LassoSIR, col = "blue")


#Plot ROC curve
par(mfrow = c(1, 1))
plot(roc_curve_wls, main = "ROC Curve ", col = "red")
lines(roc_curve_LassoSIR, col = "blue")
lines(roc_curve_lasso, col = "green")
legend("bottomright", legend = c("ROC Curve_wls","ROC Curve_lassoSIR", "ROC Curve Lasso"), 
       col = c( "red","blue","green"), lwd = 2,cex = 0.6)

# Create a data frame for the ROC curve data
roc_data_wls <- data.frame(
  Sensitivity = roc_curve_wls$sensitivities,
  Specificity = 1 - roc_curve_wls$specificities
)

# Create a data frame for the LassoSIR ROC curve data
roc_data_LassoSIR <- data.frame(
  Sensitivity = roc_curve_LassoSIR$sensitivities,
  Specificity = 1 - roc_curve_LassoSIR$specificities
)

# Create a data frame for the Lasso ROC curve data
roc_data_lasso <- data.frame(
  Sensitivity = roc_curve_lasso$sensitivities,
  Specificity = 1 - roc_curve_lasso$specificities
)

ggplot() +
  geom_line(data = roc_data_wls, aes(x = Specificity, y = Sensitivity,  color = "SIR-WLS")) +
  geom_line(data = roc_data_LassoSIR, aes(x = Specificity, y = Sensitivity,  color = "Lasso-SIR")) +
  geom_line(data = roc_data_lasso, aes(x = Specificity, y = Sensitivity,  color = "Lasso")) +
  labs(title = "ROC Curve", x = "Specificity (1 - False Positive Rate)", y = "Sensitivity (True Positive Rate)",
       color = "Method") +
  scale_color_manual(values = c("SIR-WLS" = "blue", "Lasso-SIR" = "red", "Lasso" = "green")) +
  theme_minimal()

# Calculate AUC
auc_value_wls<- auc(roc_curve_wls); auc_value_wls
auc_value_lassSIR<-auc(roc_curve_LassoSIR);auc_value_lassSIR
auc_value_lasso<-auc(roc_curve_lasso);auc_value_lasso


#### selected variables :
which(beta.lasso!= 0)
which(beta.lassosir!= 0)
which(beta.wls!= 0)
#which(beta.elasticnet!= 0)
## lambda such that the cv-error is within 1 standard error of the minimum
b.lasso <-as.matrix(beta.lasso); b.lassosir <-as.matrix(beta.lassosir) ; b.wls <-as.matrix(beta.wls)  
rownames(b.lasso)=rownames(b.lassosir)=rownames(b.wls)<-colnames( cbind(1,x.golub.train))
rownames(b.lasso)[b.lasso != 0]
rownames(b.lassosir)[b.lassosir != 0]
rownames(b.wls)[b.wls != 0]
#rownames(b.elasticnet)[b.elasticnet!=0]
select.lasso <-which(beta.lasso!= 0)
select.lassosir<- which(beta.lassosir!= 0)
select.wls<-which(beta.wls!= 0)
#select<-list(select.lm,select.lasso,select.lassosir,select.wls)
common_elements1<- Reduce(intersect,list(select.lasso,select.lassosir))
common_elements2 <- Reduce(intersect,list(select.lasso,select.wls))
common_elements3 <- Reduce(intersect,list(select.lassosir,select.wls))

# Display common elements
common_elements1
common_elements2

rownames(b.lasso)[common_elements1]
rownames(b.lasso)[common_elements2]
rownames(b.lasso)[common_elements3]


####@
prostate_preprocessed <- read.table("~/Downloads/prostate_preprocessed.txt", quote="\"", comment.char="")

####  https://data.world/nrippner/cancer-trials ------
from this paper : Sparse Sliced Inverse Regression via Cholesky Matrix Penalization 
d=3 



### comaparesation with this paper :--------
#Feature Selection for Optimized High-Dimensional Biomedical Data Using an Improved Shuffled Frog Leaping Algorithm 


### try this data https://github.com/ayanban011/HandsonML/tree/main/regression%20analysis----



###wdbc <- read.csv("~/Downloads/breast+cancer+wisconsin+diagnostic/wdbc.data", header=FALSE)
dim(wdbc)
names.wbdc<-read.csv("~/Downloads/breast+cancer+wisconsin+diagnostic/wdbc.names", header=FALSE)

install.packages("survival")
library(survival)
# Installation et chargement du package curatedBreastData
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("curatedBreastData")
library(curatedBreastData)
