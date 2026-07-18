###setwd("~/Desktop/phd /R codes ")
# if u use my file in share  use this to put it in source : 
 #source("G:/Autres ordinateurs/Mon MacBook Pro/R codes/myfunctions.R")
####### ALL library  I need ----
libr<-(function(){
  #library(suffpcr)f
  library(orthoDr)
  library(MASS)
  #install.packages("xtable")
  library(xtable)
  library(dr)
  library(survival)
  library(ggplot2)
  library(gridExtra)
  library(genlasso)
  library(SIRthresholded)
  library(Matrix)
  library(expm)
  library(glmnet)
  library(devtools)
  # library(rgl)
  library(hdi) 
  #library(qut)
  # library(datamicroarray) ### collection of dataset for clasification 
  library(latex2exp)
  library(Rdimtools)
  library(minpack.lm)
  library(LassoSIR)
  library(Metrics)
  library(ddpca)
  library(mvtnorm)
  library(spec)
  library(dplyr)
  library(knitr)
  #library(SISIR)
  library(mvtnorm)
  #library(kableExtra)
  #library(gt)
  # library(xlsx)
  #library(nlshrink)
  library(POET)
  library(colorspace)
  library(igraph)
  library(igraphdata)
  #library(corpcor)
  #library(FarmSelect)
  #library(lattice)
  #library(fastclime)
  library(sparsepca)
  #library("FactoMineR")
  #install.packages("scales")  # Install the scales package if not already installed
  library(scales)
  library(caTools)
  #library(fpp2)
 # library(devtools)
  #library(corrplot)   
  library(pROC) 
 # library(corrplot)   
  #library(datamicroarray)
  library(plsgenomics)
  library(scatterplot3d)
  library("glmpath")
  library("glmnet")
  library("penalized") 
  library(ncvreg)
  library(Ball)
  #library(BeSS)
  #library(SIS)
  #library(parallel)
  #library(foreach)
  #library(iterators)
  #install_github('ramhiser/datamicroarray')
  #library(datamicroarray)
  library(Matrix)
})()
##### ---- graphclust function ----
##### this function use igraph package  to get the blocks of X then aply SIR in each block to get estimation of beta 
### x: matrix of data n*p 
### y: vector of dimension n 
### l : number of blocks 
  ## library needed  "igraph" 
graphclust<-function(x,y,l) {
  g <- graph.adjacency(abs(cor(x)),weighted = TRUE,mode="lower")      #### get the graph 
  adjM <- as_adjacency_matrix (delete.edges(g, E(g)[weight<l]))     #### adjacency matrix 
  adjM <- as.matrix(adjM)
  #adjM
  g2 <- graph_from_adjacency_matrix(adjM - diag(p), weighted=TRUE, mode = "undirected")  ### get the graph from adjacency matrix 
  #plot(g2)
  block <- groups(cluster_label_prop(g2))        #### extract the blocks 
  #block 
  C<-dim(block)          ## number of blocks  
  #classifier la matrice en c sous matrice selon la distance entre les variables Xi
  class=list()
  i=1
  for(i in 1:C){
    class[[i]]=as.matrix(x[,block[[i]]])
    colnames(class[[i]])=block[[i]]
    i=i+1
  }
  
  B=list()
  for (j in 1:C) {
    S <- crossprod(class[[j]])/n
    k <- computeRank(class[[j]])
    D <- ddpca::DDPCA_nonconvex(S,k)
    hatTheta.ddpca <- MASS::ginv(D$L + D$A)
    R.ddpca=sqrtm(hatTheta.ddpca)
    pj=ncol(class[[j]])    ##### column number in each class 
    H=pmax(3, floor(log2(n/sqrt(pj))))  #### number of slice for each class 
    B[[j]] = mysir((class[[j]]),y,H,R.ddpca)
  }
  # Storage in the same vecetor (of the estimated directions in each class ) as estimation of beta 
 
  P=c(rep(0,C))
  for (i in 1:C) {
    P[i]=ncol(class[[i]])
  }
  #P
  Pcum=c(P[1],rep(0,C-1)) 
  for (i in 2:C) {
    Pcum[i]=Pcum[i-1]+P[i]
    i=i+1
  }
  #Pcum
  V=c(rep(0,p))
  Init=1;ii=1
  for (ii in 1:C) {
    V[Init:Pcum[ii]]=block[[ii]]
    Init=Pcum[ii]+1
  }
  # V
  I=1
  B.p=c(rep(0,p))
  for (i in 1:C) {
    B.p[I:Pcum[i]]=B[[i]]
    I=Pcum[i]+1
  }
  # B.p
  B.final=c(rep(0,p))
  for (j in 1:p)    {
    for (i in 1:p)    {
      if(V[i]==j)  B.final[j]=B.p[i]
    }
  }
  #B.final
  return(B.final)
}


###### a function to calculate the cosinus between two  vector 
#*******   x and y tow vectors 

costeta<-function(x,y){
  t= crossprod(x,y)/(norm(as.matrix(x) ,"f")*norm(as.matrix(y) ,"f"))
  return(t)
}
###### normalise a vector -------
normalize_vector <- function(v) {
  magnitude <- sqrt(sum(v^2))
  if (magnitude == 0) {
    return(v)  # To avoid division by zero
  } else {
    return(v / magnitude)
  }
}
normalize_matrix <- function(mat) {
  apply(mat, MARGIN = 2, FUN = function(x) x / sqrt(sum(x^2)))
}
##### function to calculate distance between projection space
   ### Projection matrix : 
    projec<-function(u){
      P=u%*%t(u)/norm(u,"2")^2
      return(P)
    }
  #### distance with frobinuse norm betweene projection space of the two vector u and v .
  Dproj<-function(u,v){
  P1=projec(u)
  P2=projec(v)
  D=norm(as.matrix(P1-P2) ,"f")
  return(D)
}


####### SIR function for Slised Inverse regression 
mysir<-function(x,y,H=5,R){
  library(dr)  
  n=nrow(x);p=ncol(x)
  xbar=colMeans(x) 
  #########standarisé les X(i) 
  z<-matrix(0,n,p)                              
  for(i in 1:n)
  {
    z[i,]= R%*%(x[i,]-xbar)
    
  }
  ######### couper Y en tranches et estimer la  proportion p chapeau pour chaque tranche sh
  f=dr.slices(y,H)    ### tranchage de y 
  a=f$slice.indicator    #### les indices de chaque slice 
  sizes=f$slice.sizes
  pcha=sizes/n
  ############## estimer m(h) l'eperance conditionelle de X sachant que y est ds la tranche h 
  m<-matrix(0,H,p) 
  m1<-matrix(0,H,p) 
  for(h in 1:H)   
  {   
    for (i in 1:n) 
    { if (a[i]==h ) {   m1[h,]= m1[h,]+z[i,] }
      else {i=i+1}
    }
    m[h,]= m1[h,]/sizes[h]
  } 

  ######### estimer gama la covariance de la curbe de regression inverse standarisé  ########  
  gama<-matrix(0,p,p)
  gama1<-matrix(0,p,p)
  for (j in 1:H) {
    gama1=pcha[j]* m[j,]%*%t(m[j,])
    gama=gama+gama1
  }
  ######### determination des valeurs propres et  des vecteurs propres de gama  
  V= eigen(gama)
  val=V$values # valeurs propres
  vect= V$vectors # vecteurs propres
  ########******  estimer eta pui beta ******
  eta=vect[,1]
  beta=R%*%eta
  return( beta)
}

#### # mysir funtion  (SIR without standardisation )----
#*********************
              # x: matrix of data n*p 
              # y: vector od dimension p 
              # H : number of slices 
              # package needed :  library(dr) ,
######## value : estimation of beta 
mysir2<-function(x,y,H){
       n=nrow(x);p=ncol(x)
       
  ######### couper Y en tranches et éstimer la  proportion p chapeau pour chaque tranche sh
       library(dr)
        f=dr.slices(y,H)   #### slicing 
        a=f$slice.indicator
        sizes=f$slice.sizes
        pcha=sizes/n ### P(y in sh)
  ############## estimer m(h) l'eperance conditionelle de X sachant que y est ds tranche h 
  m<-matrix(0,H,p) 
  m1<-matrix(0,H,p) 
  for(h in 1:H)   
  {   
    for (i in 1:n) 
    { if (a[i]==h ) {   m1[h,]= m1[h,]+x[i,] }
      else {i=i+1}
    }
    m[h,]= m1[h,]/sizes[h]
  } 
  ######### estimer gama la covariance de la courbe de regression inverse standarisé  ########  
  gama<-matrix(0,p,p)
  gama1<-matrix(0,p,p)
  for (j in 1:H) {
    gama1=pcha[j]* m[j,]%*%t(m[j,])
    gama=gama+gama1
  }
  ######### determination des valeurs propres et  des vecteurs propres de gama  
  V= eigen(gama)
  val=V$values # valeurs propres
  vect= V$vectors # vecteurs propres
  ########******  estimer  beta ****** 
  beta=vect[,1]
  rm(V,val, gama, gama1 ,m1, m,f,a,sizes,pcha)
  return( beta)
}

#################@la fonction sirclust singl model with c detereminat  ########################@@
sirclust1<-function(x,y,c){ 
  #decouper X en c classes
  xt=t(x)
  #clust = hclust(dist(xt), "complete")
  clust = hclust(dist(xt), "ward.D")
  plot(clust,xlab="  ")
  #par1=cutree(clust,c)  
  U=rect.hclust(clust,c)
  
  #classifier la matrice en c sous matrice selon la distance entre les variables Xi
  class=list()
  for(i in 1:c){
    class[[i]]=xt[U[[i]],]
    rownames(class[[i]])=U[[i]]
  }
  #class
  #appliquer SIR sur chaque classe i
  sir=list()
  B=list()
  library(dr)
  for (j in 1:c) {
    sir[[j]] = dr(y~t(class[[j]]))
    B[[j]]=sir[[j]][["evectors"]][,1]
  }
  #sir
  #B
  #stoquer les directions estimées sur chaque classe dans un meme vecteur comme estimation de beta 
  P=c(rep(0,c))
  for (i in 1:c) {
    P[i]=nrow(class[[i]])
  }
  #P
  Pcum=c(P[1],rep(0,c-1)) 
  for (i in 2:c) {
    Pcum[i]=Pcum[i-1]+P[i]
    i=i+1
  }
  #Pcum
  V=c(rep(0,p))
  Init=1
  for (ii in 1:c) {
    V[Init:Pcum[ii]]=U[[ii]]
    Init=Pcum[ii]+1
  }
  #V
  I=1
  B.p=c(rep(0,p))
  for (i in 1:c) {
    B.p[I:Pcum[i]]=B[[i]]
    I=Pcum[i]+1
  }
  #B.p
  B.final=c(rep(0,p))
  for (j in 1:p)    {
    for (i in 1:p)    {
      if(V[i]==j)  B.final[j]=B.p[i]
    }
  }
  #B.final
  return(B.final)
}
#################@la fonction sirclust different de sirclust1  ########################@@

### the difrrence here is that we precise R=sqrt(Sigma inverse ) and  c the clust's numbre  is not detreminat before
sirclust<-function(x,y,R) { 
              #decouper X en c classes
  xt=t(x)
  #clust = hclust(dist(xt), "complete")
  clust = hclust(dist(xt), "ward.D")
  # plot(clust,xlab="  ")
  #par1=cutree(clust,c)  
  c=floor(p/n)+1
  U=rect.hclust(clust,c)
  
  #classifier la matrice en c sous matrice selon la distance entre les variables Xi
  class=list()
  for(i in 1:c){
    class[[i]]=xt[U[[i]],]
    rownames(class[[i]])=U[[i]]
  }
  #class
  #appliquer SIR sur chaque classe i
  B=list()
  library(dr)
  for (j in 1:c) {
    B[[j]]=mysir(y,class[[j]],H,R)
  }
  #sir
  #B
  #stoquer les directions estimées sur chaque classe dans un meme vecteur comme estimation de beta 
  P=c(rep(0,c))
  for (i in 1:c) {
    P[i]=nrow(class[[i]])
  }
  #P
  Pcum=c(P[1],rep(0,c-1)) 
  for (i in 2:c) {
    Pcum[i]=Pcum[i-1]+P[i]
    i=i+1
  }
  #Pcum
  V=c(rep(0,p))
  Init=1
  for (ii in 1:c) {
    V[Init:Pcum[ii]]=U[[ii]]
    Init=Pcum[ii]+1
  }
  #V
  I=1
  B.p=c(rep(0,p))
  for (i in 1:c) {
    B.p[I:Pcum[i]]=B[[i]]
    I=Pcum[i]+1
  }
  #B.p
  B.final=c(rep(0,p))
  for (j in 1:p)    {
    for (i in 1:p)    {
      if(V[i]==j)  B.final[j]=B.p[i]
    }
  }
  #B.final
  return(B.final)
}

#################@la fonction sirclust for mulitipl index model  ndim= nobr de direction  ########################@@
sirclust2<-function(x,y,c,ndim) {  #decouper X en c classes
  xt=t(x)
  clust = hclust(dist(xt), "complete")
  plot(clust)
  par1=cutree(clust,c)
  U=rect.hclust(clust,c)
  
  #classifier la matrice en c sous matrice selon la distance entre les variables Xi
  class=list()
  for(i in 1:c){
    class[[i]]=xt[U[[i]],]
    rownames(class[[i]])=U[[i]]
  }
  #class
  #appliquer SIR sur chaque classe i
  sir=list()
  B=list()
  library(dr)
  for (j in 1:c) {
    sir[[j]] = dr(y~t(class[[j]]))
    B[[j]]=sir[[j]][["evectors"]][,1:ndim]
  }
  #sir
  #B
  #stoquer les directions estimées sur chaque classe dans un meme vecteur comme estimation de beta 
  P=c(rep(0,c))
  for (i in 1:c) {
    P[i]=nrow(class[[i]])
  }
  #P
  Pcum=c(P[1],rep(0,c-1)) 
  for (i in 2:c) {
    Pcum[i]=Pcum[i-1]+P[i]
    i=i+1
  }
  #Pcum
  V=c(rep(0,p))
  Init=1
  for (ii in 1:c) {
    V[Init:Pcum[ii]]=U[[ii]]
    Init=Pcum[ii]+1
  }
  #V
  B.p=matrix(0,nrow=p,ncol=ndim)
  In=1
  for (i in 1:c) {
    for (j in 1:ndim) {
      # betach[E:Pcum[si],j]=new.betach[[i]]
      B.p[In:Pcum[i],j]=B[[i]][,j]
    }
    In=Pcum[i]+1
  }
  # B.p
  #Les beta_i sont mal ordonée il faut les reeordonées sur selon les names des variables : 
  Bfinal=matrix(0,nrow=p,ncol=ndim)
  for (j in 1:ndim) {
    for (i in 1:p)    {
      for (s in 1:p)    {
        if(V[s]==i)  Bfinal[i,j]=B.p[s,j]
      }
    }
  }
  #Bfinal
  
  return(Bfinal)
}


########+++++++++++++++function for matrix lasso--------
##-------- x:matrix data 
#--------- y: vector of dimension p 
#---------- H : number of slices 
#----------- sigma used here is the estimation empirique 
#---- this fuunction use the function "mysir2" and glmnet package
mlsir<- function(x,y,H){
  eta=mysir2(x,y,H) ## see algorithm 2 of the reference 
  n=nrow(x)
  hatsigma= crossprod(x)/n      #### Estimation of sample covariance 
  lasso.cv<-cv.glmnet(hatsigma,eta,alpha =1, family="gaussian", nfolds=10)   ## cross validation to obtain optimal mu
  best_lam<-lasso.cv$lambda.min
  best_lasso<-glmnet(hatsigma, eta, alpha =1, family="gaussian", lambda = best_lam)  # Reconstruction of the model with the best identified lambda value
  hatbeta=coef(best_lasso)[1:p]
  return(hatbeta)
}
# for more details see : LassoSIR ; Matrix lasso :Sparse Sliced Inverse Regression Via Lasso Qian Lin,Zhao, Zhigen Liu, Jun S (arxiv:1611.06655)

#################@la fonction sirclust pour CAH avec DDpca sur chaque bloc ########################@@
sirclustdd<-function(x,y)
{  #decouper X en c classes 
  clust = hclust(as.dist(cor(x)), "ward.D")
  plot(clust,xlab="  ")
  c=floor(p/n)+1
  U=rect.hclust(clust,c)
  #U
  #classifier la matrice en c sous matrice selon la distance entre les variables Xi
  class=list()
  for(i in 1:c){
    class[[i]]=t(x)[U[[i]],]
    rownames(class[[i]])=U[[i]]
  }
  #appliquer ddpaSIR sur chaque classe i
  B=list()
  for (j in 1:c) {
    S <- crossprod(t(class[[j]]))/n
    k <- computeRank(t(class[[j]]))
    D <- ddpca::DDPCA_nonconvex(S,k)
    hatTheta.ddpca <- MASS::ginv(D$L + D$A)
    R.ddpca=sqrtm(hatTheta.ddpca)
    pj=nrow(t(class[[j]]))
    H=pmax(3, floor(log2(n/sqrt(pj))))
    B[[j]] = mysir(t(class[[j]]),y,H,R.ddpca)
  }
  #stoquer les directions estimées sur chaque classe dans un meme vecteur comme estimation de beta 
  P=c(rep(0,c))
  for (i in 1:c) {
    P[i]=nrow(class[[i]])
  }
  # P
  Pcum=c(P[1],rep(0,c-1)) 
  for (i in 2:c) {
    Pcum[i]=Pcum[i-1]+P[i]
    i=i+1
  }
  #Pcum
  V=c(rep(0,p))
  Init=1
  for (ii in 1:c) {
    V[Init:Pcum[ii]]=U[[ii]]
    Init=Pcum[ii]+1
  }
  # V
  I=1
  B.p=c(rep(0,p))
  for (i in 1:c) {
    B.p[I:Pcum[i]]=B[[i]]
    I=Pcum[i]+1
  }
  #B.p
  B.final=c(rep(0,p))
  for (j in 1:p)    {
    for (i in 1:p)    {
      if(V[i]==j)  B.final[j]=B.p[i]
    }
  }
  # B.final
  return(B.final)
}

#*****************************
# LassoSIR ----
#*****************************

# X : 
# Y : 
# 

LassoSIR.v2 <- function (X, Y, H = 0, choosing.d = "automatic", solution.path = FALSE, nfolds = 10, screening = TRUE, no.dim = 0) 
{
  
  # H, p, and n ----
  if (no.dim != 0) 
    choosing.d = "given"
  if (H == 0) {
    H <- (function() {
      H <- readline("For the continuous response, please choose the number of slices:   ")
      H <- as.numeric(unlist(strsplit(H, ",")))
      return(dim)
    })()
  }
  
  p <- dim(X)[2]
  n <- dim(X)[1]
  
  # X.ord, Y.ord, and M ----
  ORD <- order(Y)
  X <- X[ORD, ]
  Y <- Y[ORD]
  ms <- array(0, n)
  m <- floor(n/H)
  c <- n%%H
  M <- matrix(0, nrow = H, ncol = n)
  if (c == 0) {
    M <- diag(H) %x% matrix(1, nrow = 1, ncol = m)/m
    ms <- m + ms
  } else {
    for (i in 1:c) {
      M[i, ((m + 1) * (i - 1) + 1):((m + 1) * i)] <- 1/(m + 
                                                          1)
      ms[((m + 1) * (i - 1) + 1):((m + 1) * i)] <- m
    }
    for (i in (c + 1):H) {
      M[i, ((m + 1) * c + (i - c - 1) * m + 1):((m + 
                                                   1) * c + (i - c) * m)] <- 1/m
      ms[((m + 1) * c + (i - c - 1) * m + 1):((m + 
                                                 1) * c + (i - c) * m)] <- m - 1
    }
  }
  
  # Screening or not ----
  if (screening == TRUE) {
    x.sliced.mean <- M %*% X
    sliced.variance <- apply(x.sliced.mean, 2, var)
    keep.ind <- sort(order(sliced.variance, decreasing = TRUE)[1:n])
  } else {
    keep.ind <- c(1:p)
  }
  
  # Compute X.H ----
  X <- X[, keep.ind]
  X.H <- matrix(0, nrow = H, ncol = dim(X)[2])
  grand.mean <- matrix(apply(X, 2, mean), nrow = 1, ncol = dim(X)[2])
  X.stand.ord <- X - grand.mean %x% matrix(1, nrow = dim(X)[1], 
                                           ncol = 1)
  X.H <- M %*% X.stand.ord
  
  # SVD of X.H ----
  svd.XH <- svd(X.H, nv = p)           #nu = H, nv = n or p (see keep.ind)
  res.eigen.value <- array(0, p)       #the later p-H are always equal to 0
  res.eigen.value[1:dim(X.H)[1]] <- (svd.XH$d)^2/H
  
  # Setting no.dim ----
  # (cette partie n'est pas intégrée dans le travail (no.dim = 1 alors choosing.d = "given"))
  if (choosing.d == "manual") {
    plot(c(1:p), res.eigen.value, ylab = "eigen values")
    no.dim <- (function() {
      dim <- readline("Choose the number of directions:   ")
      dim <- as.numeric(unlist(strsplit(dim, ",")))
      return(dim)
    })()
  }
  if (choosing.d == "automatic") {
    beta.hat <- array(0, c(p, min(p, H)))
    Y.tilde <- array(0, c(n, min(p, H)))
    for (ii in 1:min(p, H)) {
      eii <- matrix(0, nrow = dim(svd.XH$v)[2], ncol = 1)
      eii[ii] <- 1
      eigen.vec <- solve(t(svd.XH$v), eii)
      Y.tilde[, ii] <- t(M) %*% M %*% X.stand.ord %*% 
        eigen.vec/(res.eigen.value[ii]) * matrix(1/ms, 
                                                 nrow = n, ncol = 1)
    }
    mus <- array(0, min(p, H))
    for (ii in 1:min(p, H)) {
      lars.fit.cv <- cv.glmnet(X.stand.ord, Y.tilde[, 
                                                    ii], nfolds = nfolds)
      ind <- max(which(lars.fit.cv$cvm == min(lars.fit.cv$cvm)))
      if (ind == 1) 
        ind <- 2
      lambda <- lars.fit.cv$lambda[ind]
      mus[ii] <- lambda
      lars.fit <- glmnet(X.stand.ord, Y.tilde[, ii], lambda = lambda)
      beta.hat[keep.ind, ii] <- as.double(lars.fit$beta)
    }
    temp.2 <- sqrt(apply(beta.hat^2, 2, sum)) * res.eigen.value[1:H]
    temp <- temp.2/temp.2[1]
    res.kmeans <- kmeans(temp, centers = 2)
    no.dim <- min(sum(res.kmeans$cluster == 1), sum(res.kmeans$cluster == 
                                                      2))
  }
  
  # Compute Y.tilde ----
  beta.hat <- array(0, c(p, no.dim))
  Y.tilde <- array(0, c(n, no.dim))
  for (ii in 1:no.dim) {
    eii <- matrix(0, nrow = dim(svd.XH$v)[2], ncol = 1)
    eii[ii] <- 1
    eigen.vec <- solve(t(svd.XH$v), eii)
    Y.tilde[, ii] <- t(M) %*% M %*% X.stand.ord %*% eigen.vec/(res.eigen.value[ii]) * 
      matrix(1/ms, nrow = n, ncol = 1)
  }
  
  # Estimate Beta ----
  if (solution.path == FALSE) {
    mus <- array(0, no.dim)
    if (no.dim == 1) {
      lars.fit <- glmnet(X.stand.ord, Y.tilde)
      lars.fit.cv <- cv.glmnet(X.stand.ord, Y.tilde, nfolds = nfolds)
      ind <- max(which(lars.fit.cv$cvm == min(lars.fit.cv$cvm)))
      if (ind == 1) ind <- 2
      lambda <- lars.fit.cv$lambda[ind]
      beta.keep.covariates <- predict(lars.fit,type="coef",s=lambda)
      beta.hat[keep.ind] <- as.double(beta.keep.covariates)
      
      beta.hat
    }
    else {
      for (ii in 1:no.dim) {
        lars.fit <- glmnet(X.stand.ord, Y.tilde[,ii])
        lars.fit.cv <- cv.glmnet(X.stand.ord, Y.tilde[,ii], nfolds = nfolds)
        ind <- max(which(lars.fit.cv$cvm == min(lars.fit.cv$cvm)))
        if (ind == 1) ind <- 2
        lambda <- lars.fit.cv$lambda[ind]
        mus[ii] <- lambda
        beta.keep.covariates <- predict(lars.fit,type="coef",s=lambda)
        beta.hat[keep.ind, ii] <- as.double(beta.keep.covariates)
      }
    }
    list(beta = beta.hat, eigen.value = res.eigen.value, 
         no.dim = no.dim, H = H, categorical = categorical)
  }
  else {
    lars.fit.all <- list()
    for (ii in 1:no.dim) {
      lars.fit.all[[ii]] <- glmnet(X.stand.ord, Y.tilde[, ii])
    }
    lars.fit.all
  }
}

#*********************
# X_rand ----
#*********************
# n : sample size
# p : number of covariates
# q : number of non null coefficients (relevant variables)
# type : the structure of Sigma (block or autoR)
# alpha : vector of correlations (to be used with type = "block")
# rho : correlation in auto-regressive mode or homogeneous mode  (to be used with type = "autoR")
## S : indices of activ variables (needed in case those indics were note the frst q variables ) 
X_rand <- function(n,p,q,type = "block", alpha = c(0.5,0.7,0.9), rho = 0,S=S){
  
  if(type == "block"){
    C11 <- matrix(alpha[1], nrow = q, ncol = q)
    diag(C11) <- 1
    
    C22 <- matrix(alpha[3], nrow = p-q, ncol = p-q)
    diag(C22) <- 1
    
    C12 <- matrix(alpha[2], nrow = q, ncol = p-q)
    
    sigma <- rbind(cbind(C11, C12), cbind(t(C12),C22))
   
  }
  ###Auto-regressive stricture 
  if(type == "autoR"){
    sigma <- rho^abs(outer(1:p,1:p,"-"))
  }
  
  ###  homogeneous stricture 
  if(type == "homog"){
    sigma <- diag(1,p,p)
    for (i in 1:p ){
     for ( j in 1:p) { 
      if (i!= j) sigma[i,j]=rho 
    }}
    }
  ### this type is spicified in th article : sparse sliced inverse regression via lasso 2019 
  if(type == "homog2"){
  # Initialize covariance matrix
  sigma <- matrix(0, nrow = p, ncol = p)  # Initialize 
  sigma[S, S] <- rho  # Set diagonal elements for relevant variables to 1
  sigma[-S, -S] <- rho  # Set diagonal elements for non-relevant variables to 1
  sigma[S, -S] <- 0.1  # Set correlation between relevant and non-relevant variables
  sigma[-S, S] <- 0.1  # Set correlation between non-relevant and relevant variables
  diag(sigma)<-1
  # Print the covariance matrix
  #print(Sigma)
  }
  
  
    ### symetrisation 
  if(isSymmetric(sigma)==FALSE){
    sigma <- 0.5*(sigma+t(sigma))
  } 

  x <- mvtnorm::rmvnorm(n, rep(0,p), sigma)
  return(x)
}


# A function to compute the rank of the low rank component----

# ++++++++++++++++++
# computeRank
# ++++++++++++++++++
# x : the data matrix
# standardise : logical. standardises variables, Default is TRUE
computeRank <- function(x, standardise = TRUE){
  n = nrow(x)
  p = ncol(x)
  
  means = colMeans(x)
  sds = sqrt((n-1)/n) * apply(x,2,sd)
  
  sx = x - tcrossprod(rep(1,n),means)
  if (standardise) sx = sweep(sx,2L,sds,"/",check.margin = FALSE)
  sx[is.nan(sx)] <- 0
  
  sigma_sample <- ((n-1)/n) * cov(sx)
  eigen_object <- eigen(sigma_sample, symmetric = TRUE, only.values = TRUE)
  threshold <- 1 + 2*sqrt((p-1)/(n-1)) #see 2nd reference
  principal <- eigen_object$values[which(eigen_object$values > threshold)]
  rank = length(principal)
  
  return(rank)
  
  # For more details:
  # 1- https://arxiv.org/abs/1906.00051 (Diagonally-Dominant Principal Component Analysis) [section 4.2]
  # 2- Probabilit?s, analyse des donn?es et Statistique (3e edition) [Gilbert SAPORTA] (threshold: page 172)
}



# A function to estimate the Precision matrix using ROPE method

# ++++++++++++++++++
# rope
# ++++++++++++++++++
# S : p x p sample covariance matrix (zero centered data)
# rho : The penalty parameter
# target : A symmetric and positive definite matrix. The estimate approaches target matrix as the penalty parameter (rho) increases
#          (see formula 16, page 11, ROPE paper). Default is the null matrix [ROPE (0)]
rope <- function(S, rho, target = matrix(0,ncol(S),ncol(S))){
  
  p = ncol(S)
  
  # the matrix S^star (page 11, ROPE paper)
  # if Tar is the null matrix, Sigma will be the sample covariance matrix.
  S = S - 2*rho*target
  
  # compute eigenvalues and eigenvectors of Sigma
  eigen_object = eigen(x = S,
                       symmetric = TRUE)
  
  eigenvalues = eigen_object$values
  eigenvectors = eigen_object$vectors
  
  # calculate the eigenvalues of the precision matrix
  # (using property 1, page 10, ROPE paper)
  Lambda = 2 / (eigenvalues + sqrt(eigenvalues^2 + 8*rho))
  
  # sort Lambda in ascending order
  # Note: I think the expression of the Lambda matrix in the article is wrong!
  # it should be Lambda = diag(lambda_1,....,lambda_p) where lambda_1 <= .... <= lambda_p
  # i.e the eigenvalues are sorted in ascending order (not descending order).
  ind = order(Lambda, decreasing = TRUE)
  Lambda = diag(Lambda[ind],p)
  
  # rearrange the eigenvectors in the same order as Lambda
  M = eigenvectors[,ind]
  
  # calculate the estimated precision matrix
  # (using formula (15), page 10, ROPE paper)
  PrecisionMatrix = M %*% tcrossprod(Lambda,M)
  
  return(PrecisionMatrix)
}


#  Farm.select.sir A function to estimate the beta coefficient  using Farm select  method with Lasso SIR -----

# ++++++++++++++++++
#------------------------ Farm.select.sir-------------------------------------
# ++++++++++++++++++ Arguments 
# x : n x p data  matrix (n number of observation and  p variable's number )
# y: vector of p dimension 
# H : nuber of slices 
Farm.select.sir<- function(x, y,H,bic=FALSE,alpha=1){
  ##### SIR to get eta : ----
  ORD <- order(y)
  X.ord <- x[ORD, ]
  Y.ord<- y[ORD]
  rm(ORD)
  c <- floor(n/H)
  r<- n%%H
  M <- matrix(0, nrow = n, ncol = H)
  if (r== 0) {
    M <- diag(H) %x% matrix(1, nrow = c, ncol = 1)
  } 
  else stop(" n non divisble par H ") #### else : latter 
  ##### computing Y.tild  ---- see the LASSO SIR algorithm (algo 3 of reference 1 )
  X.H <-t(X.ord) %*% M /c
  Gama<- X.H %*% t(X.H) / H
  V<-eigen(Gama)
  lamda<-V$values[1]   # largest eigenvalue of Gama
  eta<-V$vectors[1,]   # associated eigenvector  
  Y.tild<- M %*% t(M) %*%X.ord %*% eta/(c*lamda) 
  rm(r,c,M, X.H, V,Gama,lamda,eta) 
  ##### Farm.select package ----
  FAD <- FarmSelect::farm.res(X.ord, robust = FALSE)
  F.FAD <- FAD$factors
  M.projection <- F.FAD %*% tcrossprod(solve(crossprod(F.FAD)), F.FAD)
  X.res <- (diag(n) - M.projection) %*% FAD$X.res
  Y.res <- (diag(n) - M.projection) %*% Y.tild
  rm(FAD, F.FAD, M.projection) 
  
  ##### LASSO ##
  if(bic==FALSE) {
  fit.farmselect.lasso <- ncvreg::cv.ncvreg(X = X.res,
                                            y = Y.res,
                                            penalty = "Lasso",
                                            nlambda = 100,
                                            eps = 1e-04)
  beta.coef <- coef(fit.farmselect.lasso, s = fit.original.farm$beta.chosen, exact = TRUE)
  betach <- beta.coef[-1]
                     }
  
   if (bic==TRUE)  {
    fit.lasso<-bic2.select.lamda(X.res,Y.res, gama=0.25 )
    best_lam<-fit.lasso$lambda.EBIC
    best_lasso<-glmnet( X.res, Y.res, alpha =1, family="gaussian", lambda = best_lam)  # Reconstruction of the model with the best identified lambda value
    betach=coef(best_lasso)[1:p]
  }
  
  return(betach)
}
### for mor details :
#++++ 1) Sparse Sliced Inverse Regression Via Lasso: Qian Lin,Zhao, Zhigen Liu, Jun S  (arxiv:1611.06655)
#++++ 2) Factor-Adjusted Regularized Model Selection :Jianqing Fan, Yuan Ke, Kaizheng Wang (arXiv:1612.08490v2)


##### selection of the lasso or ridge   parameter lambda by  BIC creteria ######
# ++++++++++++++++++ Arguments 
# x : n x p data  matrix (n number of observation and  p variable's number )
#  y: vetor of p dimension ,
# lambda : a sequence of value for lambda defaut value is seq(0.001,1,length=50)
# alpha = to use in glmnet equal 1 for lasso or 0 for ridge default value is 1 
bic.select.lamda<-function(x,y,lambdas = seq(0.001,1,length=50),alpha=1)
{
  BIC = c() 
  for (i in 1:length(lambdas)) {   
    fit<- glmnet(x, y, alpha = 1, lambda = lambdas[i])   
    tLL <- fit$nulldev - deviance(fit)
    k <- fit$df
    nobs <- fit$nobs
    BIC = c(BIC, log(nobs)*k - tLL ) 
    i=i+1
  } 
  opt.bic<-BIC[which.min(BIC)]
  optlambda=lambdas[which.min(BIC)]
  result=result = list(BIC=BIC,bic.min=opt.bic,lamda.min=optlambda)
  return(result)
}




###### BIC selection using the formula BIC=2log(sigma^2)+klogn#####
bic2.select.lamda<-function(x,y,lambdas = seq(0.001,1,length=50),alpha=1,gama=0.5)
{
  BIC = c() ;  EBIC = c() ;

  for (i in 1:length(lambdas)) {   
    fit<- glmnet(x, y, alpha = 1, lambda = lambdas[i],family="gaussian")   
    k <- fit$df
    nobs <- fit$nobs
    hatbeta=coef(fit)[1:p]
    coef.chosen<-which(hatbeta != 0)
    hatbetaS<-hatbeta[coef.chosen]
    Xs<-as.matrix(x[,coef.chosen])
    if(nobs<k)
    {
      res<-crossprod(y-Xs%*%hatbetaS)/nobs
      
    }
    if(nobs>k) res<-crossprod(y-Xs%*%hatbetaS)/(nobs-k) 
    
    BIC = c(BIC, log(res^2)-k*log(n) ) 
    EBIC=c(EBIC, -2*log(res^2)-k*log(n)+2*gama* log(p) ) 
    i=i+1
  } 
  opt.bic<-BIC[which.min(BIC)]
  opt.Ebic<-EBIC[which.min(EBIC)]
  lambda.BIC=lambdas[which.min(BIC)]
  lambda.EBIC=lambdas[which.min(EBIC)]
  
  result = list(lambda.BIC=lambda.BIC,lambda.EBIC=lambda.EBIC )
  return(result)
}










# ++++++++++++++++++
#------------------------ SPAC.LASSO.sir-------------------------------------
# ++++++++++++++++++ Arguments 
# x : n x p data  matrix (n number of observation and  p variable's number )
# y: vector of p dimension 
# H : nuber of slices 
spac.sir<- function(x, y,H,alpha=1,method){
  ##### SIR to get eta : ----
  ORD <- order(y)
  X.ord <- x[ORD, ]
  Y.ord<- y[ORD]
  rm(ORD)
  c <- floor(n/H)
  r<- n%%H
  M <- matrix(0, nrow = n, ncol = H)
  if (r== 0) {
    M <- diag(H) %x% matrix(1, nrow = c, ncol = 1)
  } 
  else stop(" n non divisble par H ") #### else : latter 
  ##### computing Y.tild  ---- see the LASSO SIR algorithm (algo 3 of reference 1 )
  X.H <-t(X.ord) %*% M /c
  Gama<- X.H %*% t(X.H) / H
  V<-eigen(Gama)
  lamda<-V$values[1]   # largest eigenvalue of Gama
  eta<-V$vectors[1,]   # associated eigenvector  
  Y.tild<- M %*% t(M) %*%X.ord %*% eta/(c*lamda) 
  rm(r,c,M, X.H, V,Gama,lamda,eta) 
 # if(method=="fastclime"){
    if(method=="clime"){
      
  #####  fastclime package ----
 # 
#fast1<-fastclime::fastclime(X.ord, lambda.min = 0.1, nlambda = 50)
clim<-clime::clime(X.ord, lambda.min = 0.1, nlambda = 50)

# fastclime gives lambdamtx	The sequence of regularization parameters for each column,a list of p by p  precision matrices  corresponding to lambdamtx
icovlist<-fast1$icovlist 
lambdamtx<-fast1$lambdamtx
fast2<-fastclime.selector(lambdamtx, icovlist, lambda=0.5)
# fastclime.selector gives the estimated precision matrix corresponding to lambda.
d<-sqrt(diag(fast2$icov)) 
         }
  
if( method=="DESP"){
  DESP<-estimateDESP (X.ord)
  d<-diag(sqrt(DESP))
}
  if( method=="DDPCA"){
    S <- crossprod(X.ord)/n
    k <- computeRank(X.ord)
    D <- ddpca::DDPCA_nonconvex(S,k)
    hatTheta.ddpca <- MASS::ginv(D$L + D$A)
    R.ddpca=sqrtm(hatTheta.ddpca)
    d<-sqrt(diag(R.ddpca))
  }
  
  V<-diag(1/d)  # matrice qui contient l'inverse des racine  des $d-jj$ 
  X.spac<-X.ord%*%V  ## changement sur la matrice x par la relation 15 page 9 da la version 2 de l'articl 
  ##### LASSO pour X.tild de spac et Y.tild de LASSO  ###### 
  lasso.cv<-cv.glmnet(X.spac,Y.tild,alpha =1, family="gaussian", nfolds=10)   ## cross validation to obtain optimal mu
  best_lam<-lasso.cv$lambda.min
  best_lasso<-glmnet(X.spac,Y.tild, alpha =1, family="gaussian", lambda = best_lam)  # Reconstruction of the model with the best identified lambda value
  beta.coef=coef(best_lasso,exact = TRUE) 
  beta.spac <- as.matrix(beta.coef[-1]) ### les gamma-j de spac 
  hatbeta<-beta.spac/d
  #costeta(beta,hatbeta)
  return(hatbeta)
}

#### spac SIR with scad and mcp and ALasso---- 
spac2.sir<- function(x, y,H,method){
  ##### SIR to get eta : ----
  ORD <- order(y)
  X.ord <- x[ORD, ]
  Y.ord<- y[ORD]
  rm(ORD)
  c <- floor(n/H)
  r<- n%%H
  M <- matrix(0, nrow = n, ncol = H)
  if (r== 0) {
    M <- diag(H) %x% matrix(1, nrow = c, ncol = 1)
  } 
  #else stop(" n non divisble par H ") #### else : latter 
  ##### computing Y.tild  ---- see the LASSO SIR algorithm (algo 3 of reference 1 )
  X.H <-t(X.ord) %*% M /c
  Gama<- X.H %*% t(X.H) / H
  V<-eigen(Gama)
  lamda<-V$values[1]   # largest eigenvalue of Gama
  eta<-V$vectors[1,]   # associated eigenvector  
  Y.tild<- M %*% t(M) %*%X.ord %*% eta/(c*lamda) 
  rm(r,c,M, X.H, V,Gama,lamda,eta) 
    S <- crossprod(X.ord)/n
    k <- computeRank(X.ord)
    D <- ddpca::DDPCA_nonconvex(S,k)
    hatTheta.ddpca <- MASS::ginv(D$L + D$A)
    R.ddpca=sqrtm(hatTheta.ddpca)
    d<-sqrt(diag(R.ddpca))
  V<-diag(1/d)  # matrice qui contient l'inverse des racine  des $d-jj$ 
  X.spac<-X.ord%*%V  ## changement sur la matrice x par la relation 15 page 9 da la version 2 de l'articl 
  ##### LASSO pour X.tild de spac et Y.tild de LASSO  ###### 
  if(method=="SCAD"){
    
    fit <- ncvreg::cv.ncvreg( X= X.spac, y = Y.tild,
                              penalty = "SCAD",
                              nlambda = 100,
                              eps = 1e-04)
  }                                          
  
  if(method=="MCP"){
    fit <- ncvreg::cv.ncvreg( X= X.spac, y = Y.tild,
                              penalty = "MCP",
                              nlambda = 100,
                              eps = 1e-04)
  }                                          
  if(method=="Alasso"){
    # Sélectionnez lambda qui minimise la formation MSE
    # Effectuer une validation croisée par 10 pour sélectionner lambda
    ridge.cv<-cv.glmnet(x=X.spac , y=Y.tild,alpha =0, family="gaussian", type.measure="mse", nfolds=10)
    best_lambda<-ridge.cv$lambda.min
    # Reconstruction du modèle avec la meilleure valeur lambda identifiée
    best_ridge<-glmnet(x=X.spac , y=Y.tild,alpha =0, family="gaussian", lambda = best_lambda)
    
    gamma=0.5
    first.step.coef=coef(best_ridge)[-1]
    poids.lasso=(abs(first.step.coef))^(-gamma)
    ada_lasso=glmnet(x=X.spac , y=Y.tild,alpha =1, family="gaussian",penalty.factor=poids.lasso)
    # Validation croisée
    # Effectuer une validation croisée pour sélectionner lambda
    set.seed(87) #Pour la reproductibilité
    ada_lasso_cv=cv.glmnet(x=X.spac , y=Y.tild, alpha =1, family="gaussian", penalty.factor=poids.lasso)
    # Meilleur lambda par validation croisée
    best_lam_ala<-ada_lasso_cv$lambda.min
    # Reconstruction du modèle avec la meilleure valeur lambda identifiée
    fit<-glmnet(x=X.spac , y=Y.tild,alpha =1, family="gaussian",lambda = best_lam_ala)
  }

  beta.coef <- coef(fit, exact = TRUE)
  betach <- as.matrix(beta.coef[-1])
  hatbeta<-betach/d 
  
  return(hatbeta)
}


# A function to estimate the diagonal elements of a sparse precision matrix via DESP#####

# ++++++++++++++++++++++++++++
# estimateDESP
# ++++++++++++++++++++++++++++
# x : input matrix
   estimateDESP <- function(x){
  
  n <- nrow(x)
  p <- ncol(x)
  
  #compute the sample mean
  barX <- apply(x, 2, mean)
  
  #subtract the mean from all the rows
  muX.m <- tcrossprod(rep(1, n), barX)
  cx <- x - muX.m[1:n,]
  
  if (n < p){
    # estimate the coefficient matrix (square-root LASSO)
    B <- DESP::DESP_SRL_B(cx, lambda = sqrt(2*log(p)/n))
    # compute the squared partial correlations
    SPC <- DESP::DESP_SqPartCorr(B,n)
    # re-estimate the coefficient matrix by OLS
    B_OLS <- DESP::DESP_OLS_B(cx,SPC)
    # estimate the diagonal of the precision matrix
    Dhat <- DESP::DESP_RV(cx, B_OLS)
  } else{
    # estimate the diagonal of the inverse of S = (t(X) %*% X)/n
    tmp <- crossprod(x)
    Dhat <- diag(solve(tmp/n))
  }
  
  return(Dhat)
  
  # For more details:
  # 1- https://arxiv.org/abs/1504.04696v4 (On estimation of the 
   }
## SPAC_LASSO.glm function
SPAC_LASSO.glm <- function (x, y, DESP){
  
  ## New explanatory matrix
  Vhat <- diag(sqrt(1/DESP))
  SPACx <- x%*%Vhat
  
  ## Resolving Lasso problem
  model = glmnet(SPACx, y, nlambda = 100, intercept = FALSE)
  
  return(model)
}
##### FPR.funct : false positive rate fucntion and FNR.funct false negative rate -----
# ++++++++++++++++++++++++++++
# estimate the false positive rate for the estimation hatbeta of beta 
# ++++++++++++++++++++++++++++
# beta : a vector of length p 
# hatbeta : the estimation of beta .
# p : the length of the vector beta 
# q: the number of non zero element of beta 
TPR.funct <- function(beta, hatbeta) {
  p <- length(beta)
  S <- which(beta != 0)  # Indices of true active features
  q <- length(S)          # Number of true active features
  
  if (q == 0) {
    return(0)  # Avoid division by zero
  }
  
  # Vectorized computation: Count True Positives (TP)
  TP <- sum(beta != 0 & hatbeta != 0)
  
  # Compute True Positive Rate (TPR)
  TPR <- TP / q
  
  return(TPR)
}



FPR.funct <- function(beta, hatbeta) {
  p <- length(beta)
  q <- sum(beta != 0)  # Number of true active features
  if ((p - q) == 0) return(0)  # Avoid division by zero
  
  FP <- sum(hatbeta != 0 & beta == 0)  # Count of False Positives
  FPR <- FP / (p - q)  # Normalize by total negatives
  
  return(FPR)
}

FNR.funct <- function(beta, hatbeta) {
  q <- sum(beta != 0)  # Number of true active features
  if (q == 0) return(0)  # Avoid division by zero
  
  FN <- sum(hatbeta == 0 & beta != 0)  # Count of False Negatives
  FNR <- FN / q  # Normalize by total positives
  
  return(FNR)
}

TNR.funct <- function(beta, hatbeta) {
  p <- length(beta)
  q <- sum(beta != 0)  # Number of true active features
  if ((p - q) == 0) return(0)  # Avoid division by zero
  
  TN <- sum(hatbeta == 0 & beta == 0)  # Count of True Negatives
  TNR <- TN / (p - q)  # Normalize by total negatives
  
  return(TNR)
}



# ++++++++++++++++++
#------------------------ SCAD.LASSO.sir or MCP LASSO  SIR -------------------------------------
# ++++++++++++++++++ Arguments 
# x : n x p data  matrix (n number of observation and  p variable's number )
# y: vector of p dimension 
# H : nuber of slices 
N.sir<- function(x, y,H,method="SCAD"){
  ##### SIR to get eta : ----
  ORD <- order(y)
  X.ord <- x[ORD, ]
  Y.ord<- y[ORD]
  rm(ORD)
  c <- floor(n/H) 
  r<- n%%H
  M <- matrix(0, nrow = n, ncol = H)
  if (r== 0) {
    M <- diag(H) %x% matrix(1, nrow = c, ncol = 1)
  } 
  else stop(" n non divisble par H ") #### else : latter 
  ##### computing Y.tild  ---- see the LASSO SIR algorithm (algo 3 of reference 1 )
  X.H <-t(X.ord) %*% M /c
  Gama<- X.H %*% t(X.H) / H
  V<-eigen(Gama)
  lamda<-V$values[1]   # largest eigenvalue of Gama
  eta<-V$vectors[1,]   # associated eigenvector  
  Y.tild<- M %*% t(M) %*%X.ord %*% eta/(c*lamda) 
  rm(r,c,M, X.H, V,Gama,lamda,eta) 
  
  if(method=="SCAD"){
    
  fit <- ncvreg::cv.ncvreg( X= X.ord, y = Y.tild,
                                        penalty = "SCAD",
                                            nlambda = 100,
                                                   eps = 1e-04)
  }                                          
  
  if(method=="MCP"){
    fit <- ncvreg::cv.ncvreg( X= X.ord, y = Y.tild,
                              penalty = "MCP",
                              nlambda = 100,
                              eps = 1e-04)
  }                                          
  if(method=="Alasso"){
  # Sélectionnez lambda qui minimise la formation MSE
  # Effectuer une validation croisée par 10 pour sélectionner lambda
  ridge.cv<-cv.glmnet(x=X.ord , y=Y.tild,alpha =0, family="gaussian", type.measure="mse", nfolds=10)
  best_lambda<-ridge.cv$lambda.min
  # Reconstruction du modèle avec la meilleure valeur lambda identifiée
  best_ridge<-glmnet(x=X.ord , y=Y.tild,alpha =0, family="gaussian", lambda = best_lambda)
  gamma=0.5
  first.step.coef=coef(best_ridge)[-1]
  poids.lasso=(abs(first.step.coef))^(-gamma)
  ada_lasso=glmnet(x=X.ord , y=Y.tild,alpha =1, family="gaussian",penalty.factor=poids.lasso)
  # Validation croisée
  # Effectuer une validation croisée pour sélectionner lambda
  set.seed(87) #Pour la reproductibilité
  ada_lasso_cv=cv.glmnet(x=X.ord , y=Y.tild, alpha =1, family="gaussian", penalty.factor=poids.lasso)
  # Meilleur lambda par validation croisée
  best_lam_ala<-ada_lasso_cv$lambda.min
  # Reconstruction du modèle avec la meilleure valeur lambda identifiée
  fit<-glmnet(x=X.ord , y=Y.tild,alpha =1, family="gaussian",lambda = best_lam_ala)
  }
  
  beta.coef <- coef(fit, exact = TRUE)
  betach <- as.matrix(beta.coef[-1])
  #costeta(beta,betach)
  return(betach)
}


####### LASSOSIR with  sparsACP modification to get eta ---
SIR.SPCA<-function(x, y,H){
  ### SIR to get eta : 
ORD <- order(y)
X.ord <- x[ORD, ]
Y.ord<- y[ORD]
rm(ORD)
c <- floor(n/H)
r<- n%%H
M <- matrix(0, nrow = n, ncol = H)
if (r== 0) {
  M <- diag(H) %x% matrix(1, nrow = c, ncol = 1)
} 
else stop(" n non divisble par H ") #### else : latter 
##### computing Y.tild  ---- see the LASSO SIR algorithm (algo 3 of reference 1 )
X.H <-t(X.ord) %*% M /c
#if(method==PCA){
#library(sparsepca)
#pca<- PCA(t(X.H), scale.unit = TRUE, ncp = 5, ind.sup = NULL, 
 #         quanti.sup = NULL, quali.sup = NULL, row.w = NULL, 
  #        col.w = NULL, graph = TRUE, axes = c(1,2))
#eta<-pca$loading[,1]
#}
#if(method==SPCA){
  library(sparsepca)
  spca<- spca(t(X.H),k=1 ,alpha = 1e-04, beta = 1e-04, center = TRUE,
             scale = FALSE, max_iter = 1000, tol = 1e-05, verbose = FALSE)
  
  eta<-spca$loading[,1]
#}
  lamda<-spca[["eigenvalues"]][1]
Y.tild<- M %*% t(M) %*%X.ord %*% eta/(c*lamda) 
  lasso.cv<-cv.glmnet(X.ord,Y.tild,alpha =1, family="gaussian", nfolds=10)   ## cross validation to obtain optimal mu
  best_lam<-lasso.cv$lambda.min
  best_lasso<-glmnet(X.ord,Y.tild, alpha =1, family="gaussian", lambda = best_lam)  # Reconstruction of the model with the best identified lambda value
  beta.coef=coef(best_lasso,exact = TRUE) 
  betach <- as.matrix(beta.coef[-1]) ### les gamma-j de spac 
#costeta(Beta,betach)
return(betach)
}

#### predict function to estimate beta by the usual regression method

Predict<-function(x,y,method="Lasso"){

  
  if(method=="Lasso"){
    lasso.cv<-cv.glmnet(x,y,alpha =1, family="gaussian", nfolds=10)   ## cross validation to obtain optimal mu
    best_lam<-lasso.cv$lambda.min
    fit<-glmnet(x,y, alpha =1, family="gaussian", lambda = best_lam)  # Reconstruction of the model with the best identified lambda value
  }  
  
  if(method=="ALASSO-lasso")
  {
    lasso.cv<-cv.glmnet(x , y,alpha =1, family="gaussian", type.measure="mse", nfolds=10)
    best_lambda<-lasso.cv$lambda.min
    # Reconstruction du modèle avec la meilleure valeur lambda identifiée
    best_lasso<-glmnet(x, y,alpha =1, family="gaussian", lambda = best_lambda)
    gamma=0.5
    first.step.coef=coef(best_lasso)[-1]
    poids.lasso=(abs(first.step.coef))^(-gamma)
    ada_lasso=glmnet(x, y,alpha =1, family="gaussian",penalty.factor=poids.lasso)
    # Validation croisée
    # Effectuer une validation croisée pour sélectionner lambda
    set.seed(87) #Pour la reproductibilité
    ada_lasso_cv=cv.glmnet(x , y, alpha =1, family="gaussian", penalty.factor=poids.lasso)
    # Meilleur lambda par validation croisée
    best_lam_ala<-ada_lasso_cv$lambda.min
    # Reconstruction du modèle avec la meilleure valeur lambda identifiée
    fit<-glmnet(x, y,alpha =1, family="gaussian",lambda = best_lam_ala)
  }
  
  if(method=="ALASSO-ridge")
  {
    lasso.cv<-cv.glmnet(x , y,alpha =0, family="gaussian", type.measure="mse", nfolds=10)
    best_lambda<-lasso.cv$lambda.min
    # Reconstruction du modèle avec la meilleure valeur lambda identifiée
    best_lasso<-glmnet(x, y,alpha =0, family="gaussian", lambda = best_lambda)
    gamma=0.5
    first.step.coef=coef(best_lasso)[-1]
    poids.lasso=(abs(first.step.coef))^(-gamma)
    ada_lasso=glmnet(x, y,alpha =1, family="gaussian",penalty.factor=poids.lasso)
    # Validation croisée
    # Effectuer une validation croisée pour sélectionner lambda
    set.seed(87) #Pour la reproductibilité
    ada_lasso_cv=cv.glmnet(x , y, alpha =1, family="gaussian", penalty.factor=poids.lasso)
    # Meilleur lambda par validation croisée
    best_lam_ala<-ada_lasso_cv$lambda.min
    # Reconstruction du modèle avec la meilleure valeur lambda identifiée
    fit<-glmnet(x, y,alpha =1, family="gaussian",lambda = best_lam_ala)
  }
  
  if(method=="Elastic_net"){
    Enet.cv<-cv.glmnet(x,y,alpha =0.5, family="gaussian", nfolds=10)   ## cross validation to obtain optimal mu
    best_lam<-Enet.cv$lambda.min
    fit<-glmnet(x,y, alpha =0.5, family="gaussian", lambda = best_lam)  # Reconstruction of the model with the best identified lambda value
  }  
 
   
  
  
  if(method=="SCAD"){
    
    fit <- ncvreg::cv.ncvreg( x, y ,
                              penalty = "SCAD",
                              nlambda = 100,
                              eps = 1e-04)
  }                                          
  
  if(method=="MCP"){
    fit <- ncvreg::cv.ncvreg( x, y ,
                              penalty = "MCP",
                              nlambda = 100,
                              eps = 1e-04)
  }
  
  
  
  beta.coef <- coef(fit, exact = TRUE)
  betach <- as.matrix(beta.coef[-1])
  
  return(betach)
}

##### sparse LASSO vs sparse ALASSO for usual linear regression 
sl.sal<-function(x, y,k=pmin(n,p,), method="ALASSO-ridge"){
  
  
  #spca<- spca(x,k ,alpha = 1e-04, beta = 1e-04, center = TRUE,
            #  scale = FALSE, max_iter = 1000, tol = 1e-05, verbose = FALSE)
  
  #Q<-spca$loadings
  pca<-prcomp(x,rank=k,center = FALSE, scale. = FALSE)  
  Q<-pca$rotation  
  Z<-pca$x      # Z<-x%*%Q
  if(method=="Lasso"){
    beta.tild<-  predict(Z,y,method="Lasso")  # Reconstruction of the model with the best identified lambda value
  }
  
  if(method=="ALASSO-lasso"){
    beta.tild<- predict(Z,y,method="ALASSO-lasso")
  }
  
  if(method=="ALASSO-ridge"){
    beta.tild<-predict(Z,y,method="ALASSO-ridge")
  }
  
  if(method=="scad"){
    beta.tild<-predict(Z,y,method="scad")
  }
  
  betach <- Q%*% beta.tild 
  
  #costeta(Beta,betach)
  return(betach)
}
### La focntion SIR.SPCA2 without pridict function------
# SIR.SPCA2<-function(x, y,H,k=pmin(n,p,), method="ALASSO-ridge"){
#   ### SIR to get eta : 
#   ORD <- order(y)
#   X.ord <- x[ORD, ]
#   Y.ord<- y[ORD]
#   rm(ORD)
#   c <- floor(n/H)
#   r<- n%%H
#   M <- matrix(0, nrow = n, ncol = H)
#   if (r== 0) {
#     M <- diag(H) %x% matrix(1, nrow = c, ncol = 1)
#   } 
#   else stop(" n non divisble par H ") #### else : latter 
#   ##### computing Y.tild  ---- see the LASSO SIR algorithm (algo 3 of reference 1 )
#   X.H <-t(X.ord) %*% M /c
#  Gama<- X.H %*% t(X.H) / H
#   V<-eigen(Gama)
#   lamda<-V$values[1]   # largest eigenvalue of Gama
#   eta<-V$vectors[1,]   # associated eigenvector  
#   Y.tild<- M %*% t(M) %*%X.ord %*% eta/(c*lamda) 
#   rm(r,c,M, X.H, V,Gama,lamda,eta) 
# 
#   spca<- spca(X.ord,k ,alpha = 1e-04, beta = 1e-04, center = TRUE,
#               scale = FALSE, max_iter = 1000, tol = 1e-05, verbose = FALSE)
#   
# Q<-spca$loadings
# Z<-X.ord%*%Q
# if(method=="Lasso"){
#   lasso.cv<-cv.glmnet(Z,Y.tild,alpha =1, family="gaussian", nfolds=10)   ## cross validation to obtain optimal mu
#   best_lam<-lasso.cv$lambda.min
#  fit<-glmnet(Z,Y.tild, alpha =1, family="gaussian", lambda = best_lam)  # Reconstruction of the model with the best identified lambda value
# }
#  
# if(method=="ALASSO-lasso"){
#   ## steps for ADaptive 
#   lasso.cv<-cv.glmnet(x=Z , y=Y.tild,alpha =1, family="gaussian", type.measure="mse", nfolds=10)
#   best_lambda<-lasso.cv$lambda.min
#   # Reconstruction du modèle avec la meilleure valeur lambda identifiée
#   best_lasso<-glmnet(x=Z, y=Y.tild,alpha =1, family="gaussian", lambda = best_lambda)
#   gamma=0.5
#   first.step.coef=coef(best_lasso)[-1]
#   poids.lasso=(abs(first.step.coef))^(-gamma)
#   ada_lasso=glmnet(x=Z, y=Y.tild,alpha =1, family="gaussian",penalty.factor=poids.lasso)
#   # Validation croisée
#   # Effectuer une validation croisée pour sélectionner lambda
#   set.seed(87) #Pour la reproductibilité
#   ada_lasso_cv=cv.glmnet(x=Z , y=Y.tild, alpha =1, family="gaussian", penalty.factor=poids.lasso)
#   # Meilleur lambda par validation croisée
#   best_lam_ala<-ada_lasso_cv$lambda.min
#   # Reconstruction du modèle avec la meilleure valeur lambda identifiée
#   fit<-glmnet(x=Z , y=Y.tild,alpha =1, family="gaussian",lambda = best_lam_ala)
# }
# 
# if(method=="ALASSO-ridge"){
#   ## steps for ADaptive 
#     ridge.cv<-cv.glmnet(x=Z , y=Y.tild,alpha =0, family="gaussian", type.measure="mse", nfolds=10)
#     best_lambda<-ridge.cv$lambda.min
#     # Reconstruction du modèle avec la meilleure valeur lambda identifiée
#     best_ridge<-glmnet(x=Z, y=Y.tild,alpha =0, family="gaussian", lambda = best_lambda)
#     
#     gamma=0.5
#     first.step.coef=coef(best_ridge)[-1]
#     poids.lasso=(abs(first.step.coef))^(-gamma)
#     ada_lasso=glmnet(x=Z, y=Y.tild,alpha =1, family="gaussian",penalty.factor=poids.lasso)
#     # Validation croisée
#     # Effectuer une validation croisée pour sélectionner lambda
#     set.seed(87) #Pour la reproductibilité
#     ada_lasso_cv=cv.glmnet(x=Z , y=Y.tild, alpha =1, family="gaussian", penalty.factor=poids.lasso)
#     # Meilleur lambda par validation croisée
#     best_lam_ala<-ada_lasso_cv$lambda.min
#     # Reconstruction du modèle avec la meilleure valeur lambda identifiée
#     fit<-glmnet(x=Z , y=Y.tild,alpha =1, family="gaussian",lambda = best_lam_ala)
# }
# 
# 
# beta.coef <- coef(fit, exact = TRUE)
# beta.tild <- as.matrix(beta.coef[-1])
# 
#   betach <- Q%*% beta.tild  
#   #betach <-sl.sal(X.ord,Y.tild,method)
#   #costeta(Beta,betach)
#   return(betach)
# }

### La focntion SIR.SPCA2 with  pridict function------

SIR.SPCA2<-function(x, y,H,k=pmin(n,p,), method="ALASSO-ridge"){
  ### SIR to get eta : 
  ORD <- order(y)
  X.ord <- x[ORD, ]
  Y.ord<- y[ORD]
  rm(ORD)
  c <- floor(n/H)
  r<- n%%H
  M <- matrix(0, nrow = n, ncol = H)
  if (r== 0) {
    M <- diag(H) %x% matrix(1, nrow = c, ncol = 1)
  } 
  else stop(" n non divisble par H ") #### else : latter 
  ##### computing Y.tild  ---- see the LASSO SIR algorithm (algo 3 of reference 1 )
  X.H <-t(X.ord) %*% M /c
  Gama<- X.H %*% t(X.H) / H
  V<-eigen(Gama)
  lamda<-V$values[1]   # largest eigenvalue of Gama
  eta<-V$vectors[1,]   # associated eigenvector  
  Y.tild<- M %*% t(M) %*%X.ord %*% eta/(c*lamda) 
  rm(r,c,M, X.H, V,Gama,lamda,eta) 
 betach <-sl.sal(X.ord,Y.tild,k,method)
  #costeta(Beta,betach)
  return(betach)
}

###### spac  function regression 

spac<- function(x, y,Inv="DESP",method ){

  if(Inv=="fastclime"){
    #####  fastclime package ----
    fast1<-fastclime::fastclime(x, lambda.min = 0.1, nlambda = 50)
    # fastclime gives lambdamtx	The sequence of regularization parameters for each column,a list of p by p  precision matrices  corresponding to lambdamtx
    icovlist<-fast1$icovlist 
    lambdamtx<-fast1$lambdamtx
    fast2<-fastclime.selector(lambdamtx, icovlist, lambda=0.5)
    # fastclime.selector gives the estimated precision matrix corresponding to lambda.
    d<-sqrt(diag(fast2$icov)) 
  }
  
  if( Inv=="DESP"){
    DESP<-estimateDESP (x)
    d<-diag(sqrt(DESP))
  }
  if(Inv=="DDPCA"){
    S <- crossprod(x)/n
    k <- computeRank(x)
    D <- ddpca::DDPCA_nonconvex(S,k)
    hatTheta.ddpca <- MASS::ginv(D$L + D$A)
    R.ddpca=sqrtm(hatTheta.ddpca)
    d<-sqrt(diag(R.ddpca))
  }

  V<-diag(1/d)  # matrice qui contient l'inverse des racine  des $d-jj$ 
  x.spac<-x%*%V  ## changement sur la matrice x par la relation 15 page 9 da la version 2 de l'articl 

  betach <- predict(x.spac,y,method="Lasso")
  hatbeta<-betach/d 
  
  return(hatbeta)
}

#### SPac---

##### spca(sparse pca added to spac function )
spca.spac<-function(x,y,k,method){

  spca<- spca(x,k ,alpha = 1e-04, beta = 1e-04, center = TRUE,
              scale = FALSE, max_iter = 1000, tol = 1e-05, verbose = FALSE)
  
  Q<-spca$loadings
  Z<-x%*%Q
  if(method=="spac.LASSO"){
    beta.tild<-  spac(Z,y,method="Lasso")  # Reconstruction of the model with the best identified lambda value
  }
  
  if(method=="spac.scad"){
    beta.tild<- spac(Z,y,method="scad")
  }
  
  betach <- Q%*% beta.tild
}

fast.spac<- function(x, y,method){
  #####  fastclime package ----
  # library(fastclime)
  fast1<-fastclime::fastclime(x, lambda.min = 0.1, nlambda = 50)
  # fastclime gives lambdamtx	The sequence of regularization parameters for each column,a list of p by p  precision matrices  corresponding to lambdamtx
  icovlist<-fast1$icovlist 
  lambdamtx<-fast1$lambdamtx
  fast2<-fastclime.selector(lambdamtx, icovlist, lambda=0.5)
  # fastclime.selector gives the estimated precision matrix corresponding to lambda.
  d<-sqrt(diag(fast2$icov)) 
  V<-diag(1/d)  # matrice qui contient l'inverse des racine  des $d-jj$ 
  x.spac<-x%*%V  ## changement sur la matrice x par la relation 15 page 9 da la version 2 de l'articl 
  
  betach <- predict(x.spac,y,method)
  hatbeta<-betach/d 
  
  return(hatbeta)
}


soft<- function(teta, tho){
  p<-nrow(teta)
  #tetas<-c(rep(0,p))
  for(i in 1:p){ 
  teta[i]<-teta[i] * max(1-tho/abs(teta[i]),0)  }
  return(teta)
}

soft2<- function(teta, tho, d ){
  p<-nrow(teta)
  #tetas<-c(rep(0,p))
  for(i in 1:p){ 
    teta[i]<-teta[i] * max(1-tho/(2*d[i]*abs(teta[i])),0)  }
  return(teta)
}

softthresh<-function (object, lambda, gamma) 
{
  if (!is.numeric(gamma) || gamma < 0) {
    stop("gamma must be nonnegative.")
  }
  if (!any(class(object) == "fusedlasso") || !is.null(object$X) || 
      (!is.null(object$gamma) && object$gamma != 0)) {
    warning(paste("Soft-thresholding only gives a valid primal solution when applied", 
                  "to a fused lasso problem with pure fusion (gamma=0) and identity predictor matrix X."))
  }
  if (missing(lambda)) 
    lambda = object$lambda
  beta = coef(object, lambda = lambda)$beta
  lams = matrix(gamma * lambda, nrow(beta), ncol(beta), byrow = TRUE)
  beta = sign(beta) * pmax(abs(beta) - lams, 0)
  return(beta)

#<bytecode: 0x000002243609aca8>
 # <environment: namespace:genlasso>
}
#### SIR.SVD function-----
sir.svd<-function(x,y,k){
  svdd<-svd(x,nv=k,nu=k)
 # pca<-prcomp(x,rank=k,center = FALSE, scale. = FALSE)  
  Q<-svdd$v #Q<-pca$rotation  
  Z<-x%*%Q #Z<-pca$x   
  ####### Application de SIR  ----
  outsir = do.sir(Z, y,ndim=1)    ### use the function do.sir from the package "rdimtools" we can also use package "dr"
  theta<-outsir$projection 
  #### mutiplication par Q pour revenir a la dimension p 
  beta.hat<-Q%*%theta
  return(beta.hat)
}
 
sir.svd2<-function(X,y,k,H, screening=TRUE){
  beta<-c(rep(0,dim(X)[2]))
 n<-dim(X)[1]
 if(screening==TRUE){
  ## Construct the  matrix M-----
  ms <- array(0, n)
  m <- floor( n/H )
  c <- n%%H
  M <- matrix(0, nrow=H, ncol=n )
  if( c==0 )
  {
    M <- diag( H ) %x% matrix( 1, nrow=1, ncol= m )/m
    ms <- m+ms
  }
  
  else{
    for(i in 1:c){
      M[i, ( (m+1)*(i-1)+1):( (m+1)*i )] <- 1/(m+1)
      ms[ ( (m+1)*(i-1)+1):( (m+1)*i) ] <- m
    }
    for( i in (c+1): H ){
      M[i, ( (m+1)*c + (i-c-1)*m +1):( (m+1)*c+(i-c)*m)] <- 1/m
      ms[ ( (m+1)*c +(i-c-1)*m+1):( (m+1)*c+(i-c)*m) ] <- m-1
    }
  }
  ### screninng step -----
  x.sliced.mean <- M%*%X
  sliced.variance <- apply( x.sliced.mean, 2, var )
  keep.ind <- sort( order( sliced.variance, decreasing=TRUE)[1:n] )

 }
 else keep.ind<-rep(1:p) 
 
 X<- X[, keep.ind]
 pca<-prcomp(X,rank=k,center = FALSE, scale. = FALSE)  
 Q<-pca$rotation  
 Z<-pca$x   #Z<-x%*%Q
 ####### Application de SIR  ----
 outsir = do.sir(Z, y,ndim=1)    ### use the function do.sir from the package "rdimtools" we can also use package "dr"
 theta<-outsir$projection 
 #### mutiplication par Q pour revenir a la dimension p 
 beta[keep.ind]<-Q%*%theta
beta[keep.ind ]<-outsir$projection 
  return(beta)
}





####SIR.NCT------
sir.nct<-function(x,y,k,tho, threshold="soft"){
  svdd<-svd(x,nv=k,nu=k)
  # pca<-prcomp(x,rank=k,center = FALSE, scale. = FALSE)  
  Q<-svdd$v #Q<-pca$rotation  
  Z<-x%*%Q #Z<-pca$x   

  ####### Application de SIR  ----
  outsir = do.sir(Z, y,ndim=1)    ### use the function do.sir from the package "rdimtools" we can also use package "dr"
  theta<-outsir$projection 
  ###thresholding 
  
  if( threshold=="soft"){
    theta<-soft(theta,tho)
  }
  
  if( threshold=="soft2"){
    d<-pca$sdev[1:k]
    theta<-soft2(theta,tho,d) 
  }
  
  
  #### mutiplication par Q pour revenir a la dimension p 
 
  beta.hat<-Q%*%theta
  return(beta.hat)
} 

#### 2eme version of SIR.NCT  -----
sir.nct2<-function(x,y,k,tho, threshold="soft"){
  pca<-prcomp(x,rank=k,center = FALSE, scale. = FALSE)  
  Q<-pca$rotation  
  Z<-pca$x   #Z<-x%*%Q
  ####### Application de SIR  ----
  outsir = do.sir(Z, y,ndim=1)    ### use the function do.sir from the package "rdimtools" we can also use package "dr"
  theta<-outsir$projection 
  #### mutiplication par Q pour revenir a la dimension p 
  beta.hat<-Q%*%theta
  ###thresholding 
  
  if( threshold=="soft"){
    beta.hat<-soft(beta.hat,tho)
  }
  
  if( threshold=="soft2"){
    d<-pca$sdev[1:k]
    beta.hat<-soft2(beta.hat,tho,d) 
  }
  
  
  return(beta.hat)
} 


###### cross validation function; folds conception-------
folds_conception <- function(n, nfolds){
  obs_per_fold <- n %/% nfolds
  remain <- n %% nfolds
  obs_folds <- c(rep(1:nfolds, each=obs_per_fold), c(1:nfolds)[1:remain])
  if(remain==0){
    obs_folds <- obs_folds[-(n+1)]
  }
  obs_folds <- sample(obs_folds)
  obs_folds
}

##### ridge fit model function. Output beta----
Ridge_fit <- function(X_data, y_data, alpha, intercept = FALSE){
  
  ############
  # solve the following  :   argmin Beta ||X*Beta - y||^2 + alpha * ||Beta||^2
  # a closed forme solution  :  Beta = (X'X + G)^-1 X'y,  with G = alpha * Identity_matrix
  ############
  
  #    n = num_lines,     p = num_columns
  n <- dim(X_data)[1]
  p <- dim(X_data)[2]
  
  G = diag(p) * alpha
  
  if(intercept){
    # add intercept to X
    X_data = cbind(rep(1,n), X_data)
    # do not penalize intercept
    G = diag(p+1) * alpha
    G[1,1] = 0
  }
  
  Beta <- solve( t(X_data)%*% X_data  +  G ) %*% t(X_data) %*% y_data
  Beta
}

########## ridge prediction function. output y_hat-----
Ridge_predict <- function(X_test, Beta_predict, intercept = FALSE){
  
  ## add intercept to data
  if(intercept){
    X_test <- cbind(1, X_test)
  }
  y_pred <- X_test %*% Beta_predict
  y_pred
}

########## ridge crossvalidation function
Ridge_cv <- function(X_data, y_data, alpha, intercept = FALSE, nfolds){
  
  n <- dim(X_data)[1] ## nrows
  p <- dim(y_data)[2] ## ncols
  error_folds <- numeric(nfolds)
  
  folds_index <- folds_conception(n, nfolds)
  
  for(fold in 1:nfolds){
    
    # partition observations using fold
    X_train <- X_data[folds_index!=fold,]
    y_train <- y_data[folds_index!=fold]
    X_test <- X_data[folds_index==fold,]
    y_test <- y_data[folds_index==fold]
    
    # fit model with current fold
    beta_pred <- Ridge_fit(X_train, y_train, alpha, intercept = intercept)
    y_data_pred <- Ridge_predict(X_test, beta_pred, intercept = intercept)
    error_folds[fold] <- mean((y_test - y_data_pred)**2)
  }
  mean(error_folds)
}


nw <- function(x, X, Y, h=0.5, K = dnorm) {
  # Arguments
  # x: evaluation points
  # X: vector (size n) with the predictors
  # Y: vector (size n) with the response variable
  # h: bandwidth
  # K: kernel
  # Matrix of size n x length(x) (rbind() is called for ensuring a matrix
  # output if x is a scalar)
  Kx <- rbind(sapply(X, function(Xi) K((x - Xi) / h) / h))
  # Weights
  W <- Kx / rowSums(Kx) # Column recycling!
  # Means at x ("drop" to drop the matrix attributes)
  drop(W %*% Y)
}

########## SIR-Nct  crossvalidation function-----
#### ntc-cv return the  vector of mean prediction  error for nflods E1,....E10
nct_cv <- function(x, y, k,tau, threshold,nfolds=10){
  error_folds <- numeric(nfolds)
  folds_index <- folds_conception(n, nfolds)
  for(fold in 1:nfolds){
    # partition observations using fold
    X_train <- x[folds_index!=fold,]
    y_train <- y[folds_index!=fold]
    X_test <- x[folds_index==fold,]
    y_test <- y[folds_index==fold]
    
    # fit model with current fold
    beta_hat <- sir.nct(X_train, y_train,k, tau, threshold )
    
   x_beta_hat <- X_test %*% beta_hat
   
   bw0 <-np::npregbw(xdat =  x_beta_hat, ydat = y_test, regtype = "lc")
   kre0 <- np::npreg(bws = bw0)
   #x_grid<-seq(-10, 10, l = 500)
   #nw(x=x_grid,X= x_beta_hat,Y=y)
   #y_data_pred<-nw(x=x.grl, X=x_beta_hat,y=y_test)  
   y_data_pred<-kre0$mean
   
    error_folds[fold] <- mean((y_test - y_data_pred)**2)
  }
  mean(error_folds)
}


nct_cv_opt<-function(X_data , y_data, k, tau_num = 10 , threshold,nfolds=10){
tau_sequene <- seq(0.1,1, length.out =tau_num)
error_cv <- numeric(tau_num)
for(iter in 1:tau_num){
  error_cv[iter] <- nct_cv(X_data, y_data,k,tau= tau_sequene[iter],threshold="soft",nfolds )
}     

#error_cv
tau_min <-tau_sequene[which.min(error_cv)]
return(tau_min)
}

#####@---- ALL VCPCR-function-------- 
####

#### Random Cluser 
GenRandomClusters <- function(K, p){
  
  cand <- 1:K
  clusters <- sample(cand, size = p, replace = T)
  num.unique <- length(unique(clusters))
  iter <- 0
  while(num.unique != K & iter < 200){
    iter <- iter + 1
    clusters <- sample(cand, size = p, replace = T)
    num.unique <- unique(clusters)
  }
  
  if (iter == 200){
    clusters <- rep(NA, p)
    rand.indexes <- sample(1:p, size = K, replace = F)
    clusters[rand.indexes] <- 1:K
    clusters[-rand.indexes] <- sample(cand, size = p - K, replace = T)
  }
  
  return(clusters)
}
## process_fold_data----
ProcessFoldData <- function(X, y, test.id, scale.X = T, scale.y = T, bin.y = F){
  
  if (is.null(test.id)){
    train.id <- 1:nrow(X)
    test.id <- 1:nrow(X)
  } else {
    train.id <- (1:nrow(X))[-test.id]
  }
  
  y.test <- y[test.id]
  X.test <- X[test.id, ]
  
  # Get mean and sd of features in training set
  
  X.mean <- apply(X[train.id, ], 2, mean) 
  X.sd <- apply(X[train.id, ], 2, sd)
  X.sd.all <- X.sd
  wh.zero.sd <- which(X.sd == 0)
  
  if (length(wh.zero.sd) > 0){
    X.sd[wh.zero.sd] <- 1 # replace with 1 to avoid dividing by 0
  }
  
  if (scale.X == F){
    X.sd <- rep(1, ncol(X))
  }
  
  # Get mean of y in training set
  y.mean <- mean(y[train.id])
  y.sd <- sd(y[train.id])
  
  # Scale and center training set
  X.norm <- scale(X[train.id, ], center = X.mean, scale = X.sd)
  if (bin.y){
    y.norm <- y[train.id]
    y.test <- y.test
  } else {
    y.norm <- scale(y[train.id], center = y.mean, scale = ifelse(scale.y, y.sd, F))[, 1]
    y.test <- scale(y.test, center = y.mean, scale = ifelse(scale.y, y.sd, F))[, 1]
  }
  
  
  # Scale and center test set using means and sds from training set
  X.test <- scale(X.test, center = X.mean, scale = X.sd)
  
  
  return(list(X.norm = X.norm, y.norm = y.norm, 
              X.test = X.test, y.test = y.test,
              X.sd = X.sd.all,
              y.sd = y.sd))
}
GetVtVInv <- function(V){
  VtV.diag <- diag(t(V)%*%V)
  VtV.diag[VtV.diag < (.Machine)$double.eps] <- 1
  VtV.inv <- diag(1/VtV.diag)
  
  return(VtV.inv)
}
##get_fold_data-----
GetFoldData <- function(exp.data, wh.rep, cv.loop.input, inner.loop = F){
  fold.data.list <- list()
  fold.data.params <- NULL
  for(data.index in 1:length(exp.data)){ # datasets
    for (rep.index in wh.rep){ # repetitions of each dataset
      for (outer.fold.index in 1:length(cv.loop.input[[data.index]][[rep.index]]$fold.ids.outer)){ # outer folds
        
        if (inner.loop){
          inner.fold.indexes <- 1:length(cv.loop.input[[data.index]][[rep.index]]$fold.ids.inner[[outer.fold.index]])
        } else {
          inner.fold.indexes <- 1
        }
        
        for (inner.fold.index in inner.fold.indexes){
          
          for (init.index in 1:nrow(cv.loop.input[[data.index]][[rep.index]]$init.params)){
            
            if (init.index == 1){
              X <- cv.loop.input[[data.index]][[rep.index]]$X
              y <- cv.loop.input[[data.index]][[rep.index]]$y
              if (inner.loop == T){
                fold.ids <- cv.loop.input[[data.index]][[rep.index]]$fold.ids.inner[[outer.fold.index]]
                fold.index <- inner.fold.index
              } else {
                fold.ids <- cv.loop.input[[data.index]][[rep.index]]$fold.ids.outer
                fold.index <- outer.fold.index
              }
              
              fold.data <- ProcessFoldData(X = X, 
                                           y = y, 
                                           test.id = fold.ids[[fold.index]], 
                                           scale.X = scale.X,
                                           bin.y = bin.y)
              fold.data.list <- c(fold.data.list, list(fold.data))
            }
            
            
            params <- cbind(data.index = data.index,
                            rep.index = rep.index,
                            outer.fold = outer.fold.index,
                            inner.fold = ifelse(inner.loop, inner.fold.index, NA),
                            init.index = init.index,
                            K.init = cv.loop.input[[data.index]][[rep.index]]$init.params[init.index, "K"],
                            wh.init = cv.loop.input[[data.index]][[rep.index]]$init.params[init.index, "init"],
                            wh.fold.data = length(fold.data.list))
            fold.data.params <- rbind.data.frame(fold.data.params, params)
          }
          
        }
        
      }
    }
  }
  
  return(list(fold.data.list = fold.data.list, fold.data.params = fold.data.params))
}
###VCPCR----
VCPCR <- function(fold.data, lambda, delta, Cs.init, X.sd = NULL, alpha = NA, bin.y = F){
  
  
  X <- fold.data$X.norm
  y <- fold.data$y.norm
  var.y <- var(y)
  
  if (!is.na(alpha)){
    if (alpha == 1){
      if (bin.y){
        lambda.lasso <- max( abs(t(y - mean(y)*(1-mean(y))) %*% X ) )/ ( nrow(X))
      } else {
        lambda.lasso <- max((1/nrow(X))*abs(t(X)%*%y), na.rm = T)
      }
      
    } else {
      lambda.lasso <- 1
    }
  }
  
  
  if (is.null(X.sd)){
    X.sd <- rep(1, ncol(X))
  }
  
  if (bin.y){
    fam <- "binomial"
  } else {
    fam <- "gaussian"
  }
  
  if (is.na(delta)){
    weights <- rep(1, ncol(X))
  } else {
    weights <- glmnet(x = X, y = y, family = fam,
                      alpha = alpha, lambda = delta*lambda.lasso, standardize = F)$beta[, 1]
  }
  
  
  
  # Step 1: W-SOS-NMF
  
  res <- WSOSNMF(X = X, weights = weights, lambda = lambda, Cs.init = Cs.init, X.sd = X.sd)
  
  # Step 2: Prediction
  
  V <- res$V
  M <- X%*%diag(X.sd)%*%V
  
  if (sum(abs(V)) > 0){
    temp <- GetBetasStep2(M = M, V = V, y = y, is.lasso = F, lambda = 0, bin.y = bin.y)
    U <- temp$M
    V <- temp$V
    betas <- temp$betas
    fit <- temp$fit
    Cs <- apply(V, 1, which.max)
    wh.non.zero <- which(apply(V, 1, sum) != 0)
    Cs[-wh.non.zero] <- max(Cs) + 1
  } else {
    U <- as.matrix(rep(0, nrow(X)))
    V <- V[, 1, drop = F]
    Cs <- rep(1, ncol(X))
    betas <- 0
    fit <- NULL
  }
  
  y.pred <- GetPredictionsTwoStep(V = V, Cs = Cs, betas = betas, X.train = X, 
                                  X.test = fold.data$X.test, fit = fit, bin.y = bin.y)
  
  return(list(V = V, U = U, Cs = Cs, betas = betas, weights = weights, y.pred = y.pred, y.true = fold.data$y.test))
  
}

###VCPCRAllHP-----
VCPCRAllHP <- function(fold.data, hps, Cs.init, bin.y = F){
  
  
  
  X.sd <- rep(1, ncol(fold.data$X.norm))
  
  
  all.hp <- lapply(1:nrow(hps), function(hp.index) 
    VCPCR(fold.data = fold.data,
          lambda = hps$lambda[hp.index], 
          delta = hps$delta[hp.index],
          X.sd = X.sd,
          Cs.init = Cs.init,
          alpha = hps$alpha[hp.index],
          bin.y = bin.y))
  
  
  return(all.hp)
}
### GetPredictionsTwoStep----
GetPredictionsTwoStep <- function(V, Cs, betas, X.train, X.test, fit, bin.y = F){
  require(glmnet)
  require(brglm2)
  if (is.null(fit)){
    y.pred <- rep(0, nrow(X.test))
    if (bin.y){
      y.pred <- cbind.data.frame(y.pred = y.pred, y.probs = y.pred)
    }
    
  } else {
    
    M.test <- as.matrix(X.test%*%V)
    if (ncol(M.test) < length(betas)){
      M.test <- cbind(M.test, 
                      matrix(0, 
                             nrow = nrow(M.test), 
                             ncol = length(betas) - ncol(M.test)))
    }
    if (bin.y){
      # preds = log(1-p)/log(p)
      # if p > 0.5 --> 1 
      # preds > 1 if p > 0.5
      # https://glmnet.stanford.edu/articles/glmnet.html#logistic-regression-family-binomial-
      if (class(fit)[1] == "glm"){
        preds <- as.matrix(predict.glm(object = fit, newdata = cbind.data.frame(M.test, y = rep(1, nrow(M.test))), type = "response"))
        
      } 
      if (class(fit)[1] == "brglmFit"){
        preds <- as.matrix(predict.glm(object = fit, newdata = cbind.data.frame(M.test, y = rep(1, nrow(M.test))), type = "response"))
        
      }
      if ("glmnet" %in% class(fit)){
        preds <- as.matrix(predict(fit, newx = M.test, type = "response", family = "binomial"))
        
      }
      
      
      preds <- preds[, ncol(preds)]
      y.pred <- as.numeric(preds > 0.5)
      y.pred <- cbind.data.frame(y.pred = y.pred, y.probs = preds)
      
    } else {
      
      y.pred <- as.numeric(M.test%*%betas)
      
    }
  }
  
  return(y.pred)
}
### Getfold data lasso -----
GetFoldDataLasso <- function(exp.data, wh.rep, cv.loop.input, inner.loop = F){
  fold.data.list <- list()
  fold.data.params <- NULL
  for(data.index in 1:length(exp.data)){ # datasets
    for (rep.index in wh.rep){ # repetitions of each dataset
      for (outer.fold.index in 1:length(cv.loop.input[[data.index]][[rep.index]]$fold.ids.outer)){ # outer folds
        
        if (inner.loop){
          inner.fold.indexes <- 1:length(cv.loop.input[[data.index]][[rep.index]]$fold.ids.inner[[outer.fold.index]])
        } else {
          inner.fold.indexes <- 1
        }
        
        for (inner.fold.index in inner.fold.indexes){
          
          
          X <- cv.loop.input[[data.index]][[rep.index]]$X
          y <- cv.loop.input[[data.index]][[rep.index]]$y
          if (inner.loop == T){
            fold.ids <- cv.loop.input[[data.index]][[rep.index]]$fold.ids.inner[[outer.fold.index]]
            fold.index <- inner.fold.index
          } else {
            fold.ids <- cv.loop.input[[data.index]][[rep.index]]$fold.ids.outer
            fold.index <- outer.fold.index
          }
          
          fold.data <- ProcessFoldData(X = X, 
                                       y = y, 
                                       test.id = fold.ids[[fold.index]], 
                                       scale.X = scale.X,
                                       bin.y = bin.y)
          fold.data.list <- c(fold.data.list, list(fold.data))
          
          
          
          params <- cbind(data.index = data.index,
                          rep.index = rep.index,
                          outer.fold = outer.fold.index,
                          inner.fold = ifelse(inner.loop, inner.fold.index, NA),
                          
                          wh.init = cv.loop.input[[data.index]][[rep.index]]$init.params[init.index, "init"],
                          wh.fold.data = length(fold.data.list))
          fold.data.params <- rbind.data.frame(fold.data.params, params)
        }
        
      }
      
    }
  }
  
  
  return(list(fold.data.list = fold.data.list, fold.data.params = fold.data.params))
}

### function SOSNMF ------
SOSNMF <- function(X.tilde, X.tilde.n, X.sd, V, gamma, max.iter = 500){
  
  
  iter <- 1
  
  
  VtV.inv <- GetVtVInv(V)
  U.tilde <- crossprod(t(X.tilde), V)%*%VtV.inv
  U <- ScaleNonzeroVar(U.tilde, center = T)
  
  err.old <- Inf
  err <- ObjectiveSOSNMF(X.tilde, U, V, X.sd, gamma)
  
  
  while (iter < max.iter && (abs(err.old - err)/err) > 1e-4){
    
    # Update V
    V <- SparseOrthogNNLS(M = X.tilde, U = U, Mn = X.tilde.n, gamma = gamma, M.sd = X.sd)
    
    # Update U
    VtV.inv <- GetVtVInv(V)
    U.tilde <- crossprod(t(X.tilde), V)%*%VtV.inv
    U <- ScaleNonzeroVar(U.tilde, center = T)
    
    
    err.old <- err
    err <- ObjectiveSOSNMF(X.tilde, U, V, X.sd, gamma)
    iter <- iter + 1
  }
  
  return(list(U = U, V = V))
}

####------SparseOrthogNNLS function----
SparseOrthogNNLS <- function(M, U, Mn, gamma, M.sd){
  
  #[(1/(n-1))*u_k^t m_j - \gamma]_+
  
  n <- nrow(M) - 1
  
  A <- crossprod((M), U) # Covariance of M and U
  wh.max <- apply(A, 1, which.max)
  R <- matrix(0, nrow = nrow(A), ncol = ncol(A))
  for (r in unique(wh.max)){
    wh.r <- which(wh.max == r)
    
    R[wh.r, r] <- pmax((((1/(n))*crossprod(M[, wh.r, drop = F], U[, r, drop = F])) - as.matrix((gamma)*M.sd[wh.r])), 0)
  }
  
  wh.zero <- which(R < (.Machine)$double.eps)
  R[wh.zero] <- 0
  Rn <- R
  return(Rn)
}


####----ScaleNonzeroVar----
ScaleNonzeroVar <- function(X, center = TRUE){
  scale.vec <- apply(X, 2, sd)
  scale.vec[scale.vec == 0] <- 1
  X.norm <- scale(X, center = center, scale = scale.vec)
  return(X.norm)
}

##---NonzeroCoeffs----
NonzeroCoeffs <- function(betas){
  s <- sum(sign(betas) != 0)
  return(s)
}

##----NormalizeX -----
NormalizeX <- function(X, scale.vars = F){
  
  X.norm <- scale(X, center = T, scale = F)
  if (scale.vars == T){
    X.norm <- scale(X.norm, center = F, scale = apply(X, 2, sd))
  } 
  
  return(X.norm)
}

#---L2Norm -----
L2Norm <- function(x, nonzero = F) {
  
  norm.val <- sqrt(sum(x*x))
  
  if (nonzero == T){
    norm.val <- ifelse(norm.val == 0, 1, norm.val)
  }
  
  return(norm.val)
  
}

### MakeBinaryVMatrix -----
MakeBinaryVMatrix <- function(clusters){
  
  nonzero <- clusters != 0
  
  if (sum(nonzero) == 0){
    V <- matrix(0, nrow = length(clusters), ncol = 1)
  } else {
    clusters[nonzero] <- as.numeric(factor(clusters[nonzero]))
    clust.indexes <- cbind(1:length(clusters), clusters)[nonzero, ]
    V <- matrix(0, nrow = length(clusters), ncol = max(clusters[nonzero]))
    V[clust.indexes] <- 1
  }
  
  
  return(V)
  
}


###------ScaleNonzeroL2Norm
ScaleNonzeroL2Norm <- function(X, center = FALSE){
  scale.vec <- apply(X, 2, L2Norm, nonzero = TRUE)
  X.norm <- scale(X, center = center, scale = scale.vec)
  return(X.norm)
}


###----ObjectiveSOSNMF -----
ObjectiveSOSNMF <- function(X.tilde, U, V, X.sd, gamma){
  
  n <- nrow(X.tilde)
  diff.PCA <- X.tilde - U%*%t(V)
  LS.PCA <- (1/(2*(n-1)))*sum(diff.PCA^2)
  lasso.term.clust <- sum(abs(t(V) %*% diag(gamma*X.sd)))
  obj <- LS.PCA + lasso.term.clust
  
  return(obj)
}





#--GetBetasStep2------
GetBetasStep2 <- function(M, V, y, is.lasso = T, lambda, bin.y = F){
  nonzero.var <- which(apply(M, 2, sd) != 0)
  if (length(nonzero.var) == 0){
    
    betas <- rep(0, ncol(M))
    fit <- NULL
    
  } else {
    M <- M[, nonzero.var, drop = F]
    V <- V[, nonzero.var, drop = F]
    
    if (bin.y){
      fam <- "binomial"
    } else {
      fam <- "gaussian"
      fit <- NULL
    }
    
    added.col <- F
    if (ncol(M) == 1){
      M <- cbind(M, 0)
      added.col <- T
    }
    
    
    if (bin.y & lambda == 0){
      
      fit <- glm(y ~ ., data = cbind.data.frame(M, y = y), family=binomial(link = "logit"),
                 method = "brglmFit",
                 type = "AS_mean", maxit = 1000)
      fit$betas <- fit$coefficients[-1]
      
    } else {
      fit <- glmnet(x = M, y = y, family = fam,
                    alpha = as.numeric(is.lasso), lambda = lambda, standardize = F)
    }
    
    
    
    if (added.col){
      M <- M[, -ncol(M), drop = F]
    }
    
    betas <- as.matrix(fit$beta)
    betas <- betas[, ncol(betas)]
  }
  return(list(M = M, V= V, betas = betas, fit = fit))
}
#### function WSOSNMF-----
WSOSNMF <- function(X, weights, lambda, Cs.init, X.sd){
  if (NonzeroCoeffs(weights) > 0){
    
    
    # Add weights to variables
    W <- diag(weights)
    XW <- X%*%W
    XW.n <- ScaleNonzeroL2Norm(XW)
    
    # initialize V
    V <- MakeBinaryVMatrix(Cs.init)
    
    # calculate max gamma val
    VtV.inv <- GetVtVInv(V)
    U.tilde <- crossprod(t(XW), V)%*%VtV.inv
    U <- ScaleNonzeroVar(U.tilde, center = T)
    cor.vals <- (1/(nrow(X)-1))*(t(U)%*%XW)
    gamma.max <- max(cor.vals)
    
    
    # Optimize U and V
    res <- SOSNMF(X.tilde = XW, 
                  X.tilde.n = XW.n, 
                  X.sd = X.sd, 
                  V = V, 
                  gamma = lambda*gamma.max)
  } else {
    
    V <- matrix(rep(0, ncol(X)))
    U <- matrix(rep(0, nrow(X)))
    res <- list(U = U, V = V)
    
  }
  
  return(res)
}

#### Update_Cs_CEnet-----
#update clusters with betas fixed
UpdateCs <- function(Xb, K, Cs.init, nstart=1){ 
  V <- MakeBinaryVMatrix(Cs.init)
  V <- scale(V, center = F, scale = apply(V, 2, sum))
  centers <- t(Xb%*%V)
  twonorm <- sqrt(apply(Xb^2,2,sum))
  if(sum(twonorm!=0)==0) {
    Cs <- rep(1,ncol(Xb))
  } else {
    wcss.old <- GetWCSS(Xb, Cs.init)  
    km.out <- kmeans(t(Xb), centers = unique(centers), algorithm = "Lloyd", iter.max = 100)
    
    if(sum(km.out$withinss) < wcss.old) {
      Cs <- km.out$cluster
    } else {
      Cs <- Cs.init
    }
    
  }  
  
  return(Cs)
}


##---GenInactiveVarExpData -----
GenInactiveVarExpData <- function(n, p, pk.inactive, pk.active, betas, SNR.model, n.rep, cor.x, inactive.type, order.types, wh.grouped = NULL){
  
  require(MASS)
  
  set.seed(92837)
  
  var.types <- AssignVarType(p = p, pk.active = pk.active, pk.inactive = pk.inactive, order.types = order.types)
  block.ids <- AssignBlock(p = p, pk.active = pk.active, pk.inactive = pk.inactive, order.types = order.types, inactive.type = inactive.type)
  var.info <- cbind.data.frame(index = 1:p,
                               var.type = var.types,
                               block.id = block.ids)
  if (!is.null(wh.grouped)){
    if (sum(wh.grouped) == 0){
      var.info$block.id <- 0
    } else {
      var.info$block.id[-wh.grouped] <- 0
    }
  }
  
  n.blocks <- length(unique(var.info$block.id[var.info$block.id != 0]))
  
  
  # Define Sigma_X
  Sigma <- diag(1, nrow = p, ncol = p)
  
  if (n.blocks > 0){
    pairs.list <- lapply(1:n.blocks, function(k) t(combn(var.info[which(var.info[, "block.id"] == k), "index"], 2)))
    for (block.index in 1:n.blocks){
      pairs.cor <- pairs.list[[block.index]]
      for (i in 1:nrow(pairs.cor)){
        row <- pairs.cor[i, 1]
        col <- pairs.cor[i, 2]
        Sigma[row, col] <- Sigma[col, row] <- cor.x
      }
    }
  }
  
  
  # Generate X
  X <- mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  
  # true betas
  betas.x <- rep(0, p)
  K <- length(pk.active)
  active.blocks <- unique(var.info$block.id[var.info$var.type == "active"])
  for (k in 1:K){
    wh.k <- which(var.info$block.id == active.blocks[k])
    betas.x[wh.k[1:(pk.active[k])]] <- betas[k]
  }
  
  # model error
  sd.y <- sqrt(c((t(betas.x)%*%Sigma%*%betas.x)/SNR.model))
  y <- X%*%betas.x + rnorm(n = n, sd = sd.y)
  
  ## Ground-truth clusters
  
  # Definition 1: variable blocks only
  Cs.true1 <- var.info$block.id
  # Definition 2: active in clusters + inactive cluster
  Cs.true2 <- Cs.true1
  Cs.true2[which(var.info$var.type == "inactive")] <- 0
  # Definition 3: active in clusters + inactive in clusters
  Cs.true3 <- Cs.true2
  inactive.blocks <- unique(var.info$block.id[var.info$var.type == "inactive" & var.info$block.id != 0])
  for (k in inactive.blocks){
    wh <- which(var.info$block.id == k & var.info$var.type == "inactive")
    Cs.true3[wh] <- max(Cs.true3) + 1
  }
  
  train.index <- split(1:n, seq(1, n, by=(n/(n.rep))))
  
  train <- lapply(train.index, function(x) list(X = X[x,],y = y[x], betas = betas, support = betas, sigma.x = Sigma, betas.x = betas.x))
  
  clusters <- list(def1 = Cs.true1,
                   def2 = Cs.true2,
                   def3 = Cs.true3)
  output <- list(train = train, clusters = clusters, sd.y = sd.y)
  return(output)
}


####AssignVarType -----
AssignVarType <- function(p, pk.active, pk.inactive, order.types){
  
  types <- NULL
  active.index <- 1
  inactive.index <- 1
  for (group.index in 1:length(order.types)){
    order.type <- order.types[group.index]
    
    if (order.type == 1){
      types <- c(types, rep("active", pk.active[active.index]))
      active.index <- active.index + 1
    }
    
    if (order.type == 0){
      types <- c(types, rep("inactive", pk.inactive[inactive.index]))
      inactive.index <- inactive.index + 1
    }
    
  }
  
  types <- c(types, rep("inactive", p - length(types)))
  return(types)
}


### AssignBlock -
AssignBlock <- function(p, pk.active, pk.inactive, order.types, inactive.type){
  
  block <- NULL
  active.index <- 1
  inactive.index <- 1
  block.index <- 1
  for (group.index in 1:length(order.types)){
    order.type <- order.types[group.index]
    
    if (order.type == 1){
      block <- c(block, rep(block.index, pk.active[active.index]))
      active.index <- active.index + 1
      block.index <- block.index + 1
    }
    
    if (order.type == 0){
      
      if (inactive.type == "corr active and inactive"){
        
        if (order.types[group.index - 1] == 1){
          block.index <- block.index - 1
        } 
        
        block <- c(block, rep(block.index, pk.inactive[inactive.index]))
        inactive.index <- inactive.index + 1
        block.index <- block.index + 1
        
      }
      
      if (inactive.type == "corr inactive"){
        
        block <- c(block, rep(block.index, pk.inactive[inactive.index]))
        inactive.index <- inactive.index + 1
        block.index <- block.index + 1
        
      }
      
      if (inactive.type == "uncorrelated"){
        block <- c(block, rep(0, pk.inactive[inactive.index]))
        inactive.index <- inactive.index + 1
        
      }
      
      
    }
    
  }
  
  block <- c(block, rep(0, p - length(block)))
  return(block)
}
###########################################@@@
#########----------------@VCPCR_SIR-----
VCPCR.SIR<-function(X,y,K,p){
  #source("~/Desktop/phd /R codes /VC-pcr codes /all VCPCR function.R")
Cs.init<-GenRandomClusters(K,p)
Prosfold<-ProcessFoldData(X, y, test.id=NULL , scale.X = T , scale.y = T , bin.y = F)
#source("~/Desktop/phd /R codes /VC-pcr codes /all VCPCR function.R")
VCPCR<-VCPCR(fold.data=Prosfold , lambda=0.5 , delta=1, Cs.init , X.sd = NULL, alpha = 0, bin.y = F)
V<-VCPCR$V
Z<-X%*%V
outsir = do.sir(Z, y,ndim=1)    ### use the function do.sir from the package "rdimtools" we can also use package "dr"
theta<-outsir$projection 
Beta.hat<-V%*%theta
Beta.hat <- Beta.hat/sqrt( sum( Beta.hat^2 ) )
return(Beta.hat )
}

##### FDR control function with different method ------
     # #### X a data  matrix  with (2*n)  row and P column 
        ### Y a 2*n vector   
                # H : number of slices 
                    # alf : the treshold coef 

FDR<-function(X,y,n,p,q,H,alf, metho="usual",K){

set.seed(123)
#### Deviser les données en deux echontillons  "training set (donnée d'entraînement)" et "test set"
ind<- sample(seq_len(nrow(X)), size = n)
X1<-X[ind,];Y1<-y[ind]
X2<-X[-ind,];Y2<-y[-ind]
      if(p<n){
  #### Lm in D1
  ######### couper Y en tranches et estimer la  proportion p chapeau pour chaque tranche sh
  set.seed(123)
  f=dr.slices(Y1,H)    ### tranchage de y 
  a=f$slice.indicator    #### les indices de chaque slice 
  ############## estimer m(h) l'eperance conditionelle de X sachant que y est ds la tranche h 
  beta.hat1<-matrix(nrow=H,ncol=p)
  for(h in 1:H)   
  {   
    indicator <-c(rep(0,n))
    for (i in 1:n) 
    { if (a[i]==h ) { indicator[i]=1 }
      else {i=i+1}
    }
    
    lm<-lm(indicator ~ X1)
    beta.hat1[h,]=lm$coefficients[-1]
    
    h=h+1
  } 
  
  ##### estimed Beta ind D2
  set.seed(123)
  f=dr.slices(Y2,H)    ### tranchage de y 
  a=f$slice.indicator    #### les indices de chaque slice 
  beta.hat2<-matrix(nrow=H,ncol=p)
  for(h in 1:H)   
  {   
    indicator <-c(rep(0,n))
    for (i in 1:n) 
    { if (a[i]==h ) { indicator[i]=1 }
      else {i=i+1}
    }
    
    lm<-lm(indicator ~ X2)
    beta.hat2[h,]=lm$coefficients[-1]
    h=h+1
  } 
  #### calcul des Sjk 
  I<-diag(p)
  sm2=sm1<-matrix(0,p,p)
  for (i in 1:n) {
    sm1= sm1+X1[i,]%*%t(X1[i,])
    sm2= sm2+X2[i,]%*%t(X2[i,])
  }
  S2=S1<-vector(mode="numeric",p)
  for (j in 1:p) {
    S1[j]<-t(I[,j])%*% solve(sm1/n)%*%I[,j]
    S2[j]<-t(I[,j])%*% solve(sm2/n)%*%I[,j]
    
  }
  
  #### calcul of W 
  W<-vector(mode="numeric",p)
  for( j in 1:p) W[j] <- t(beta.hat1[,j])%*%beta.hat2[,j]/ (sqrt(S1[j]*S2[j]))
  ### calcul of Aplus 
  t = sort(abs(W))
  Ta = sapply(t,function(x){(length(W[W<=(-x)]))/max(1,length(W[W>=x]))})   #for SD

}

        if(n<p){
          if(metho=="usual"){
            #### Lasso in D1
            ######### couper Y en tranches et estimer la  proportion p chapeau pour chaque tranche sh
            set.seed(123)
            f=dr.slices(Y1,H)    ### tranchage de y 
            a=f$slice.indicator 
            beta.hat1<-matrix(nrow=H,ncol=p)
            for(h in 1:H)   
            {   
              indicator <-c(rep(0,n))
              for (i in 1:n) 
              { if (a[i]==h ) { indicator[i]=1 }
                else {i=i+1}
              }
              
              lasso.cv<-cv.glmnet(indicator,x=X1,alpha =1,family="gaussian",  nfolds=10)   ## cross validation to obtain optimal mu
              best_lam<-lasso.cv$lambda.min
              best_lasso<-glmnet(indicator,x=X1, alpha =1, family="gaussian", lambda = best_lam)  # Reconstruction of the model with the best identified lambda value
              beta.hat1[h,]=coef(best_lasso)[-1]
              h=h+1
            } 
            #### la liste S of sparsity 
            betaS<-colSums(abs(beta.hat1)) 
            
            S<-which(betaS!=0) #### My S is indice of null coef 
            Sbar<-which(betaS==0)
            ##### estimed Beta ind D2
            set.seed(123)
            f=dr.slices(Y2,H)    ### tranchage de y 
            a=f$slice.indicator    #### les indices de chaque slice 
            beta.hat2<-matrix(nrow=H,ncol=p)
            for(h in 1:H)   
            {   
              indicator <-c(rep(0,n))
              for (i in 1:n) 
              { if (a[i]==h ) { indicator[i]=1 }
                else {i=i+1}
              }
              
              lasso.cv<-cv.glmnet(indicator,x=X2,alpha =1,family="gaussian",  nfolds=10)   ## cross validation to obtain optimal mu
              best_lam<-lasso.cv$lambda.min
              best_lasso<-glmnet(indicator,x=X2, alpha =1, family="gaussian", lambda = best_lam)  # Reconstruction of the model with the best identified lambda value
              beta.hat2[h,]=coef(best_lasso)[-1]
              h=h+1
            } 
            
            for(j in Sbar)   
            {beta.hat2[,j]<-rep(0,H)}
            
            
            
            #### calcul des Sjk 
            S2=S1<-vector(mode="numeric",p)
            s<-length(S)
            I<-diag(s)
            XS1<-X1[,S]
            XS2<-X2[,S]
            sm2=sm1<-matrix(0,s,s)
            for (i in 1:n) {
              sm1= sm1+XS1[i,]%*%t(XS1[i,])
              sm2= sm2+XS2[i,]%*%t(XS2[i,])
            }
            
            for (j in 1:s)
            {
              S1[j]<-t(I[,j])%*% solve(sm1/n)%*%I[,j]
              S2[j]<-t(I[,j])%*% solve(sm2/n)%*%I[,j]
              
            }
            
          
          }
        
          if(metho=="spac"){
            #choix <- menu(c("DDPCA", "fastclime"), title = "Choisissez une methode pour Spac: ")
            # Utilisation du choix de l'utilisateur
            #if (choix == 1) {
              #Inv <- "DDPCA"
            #} #else if (choix == 2) {Inv <- "fastclime"
            #} else {
              # Gérer le cas où l'utilisateur fait un choix invalide
             # stop("Choix invalide.")  }
            
            
            
            ######### couper Y en tranches et estimer la  proportion p chapeau pour chaque tranche sh
            set.seed(123)
            f=dr.slices(Y1,H)    ### tranchage de y 
            a=f$slice.indicator 
            beta.hat1<-matrix(nrow=H,ncol=p)
            for(h in 1:H)   
            {   
              indicator <-c(rep(0,n))
              for (i in 1:n) 
              { if (a[i]==h ) { indicator[i]=1 }
                else {i=i+1}
              }
              beta.hat1[h,]=spac(x=X1,y=indicator,Inv="DDPCA",method = "Lasso")

              h=h+1
            } 
            #### la liste S of sparsity 
            betaS<-colSums(abs(beta.hat1)) 
            S<-which(betaS!=0) #### My S is indice of null coef 
            Sbar<-which(betaS==0)
            ##### estimed Beta ind D2
            set.seed(123)
            f=dr.slices(Y2,H)    ### tranchage de y 
            a=f$slice.indicator    #### les indices de chaque slice 
            beta.hat2<-matrix(nrow=H,ncol=p)
            for(h in 1:H)   
            {   
              indicator <-c(rep(0,n))
              for (i in 1:n) 
              { if (a[i]==h ) { indicator[i]=1 }
                else {i=i+1}
              }
              beta.hat2[h,]=spac(x=X2,y=indicator,Inv="DDPCA",method="Lasso")
              rm(indicator)
              h=h+1
            } 
            
            for(j in Sbar)   
            {beta.hat2[,j]<-rep(0,H)}
            
            #### calcul des Sjk 
            S2=S1<-vector(mode="numeric",p)
            s<-length(S)
            I<-diag(s)
            XS1<-X1[,S]
            XS2<-X2[,S]
            sm2=sm1<-matrix(0,s,s)
            for (i in 1:n) {
              sm1= sm1+XS1[i,]%*%t(XS1[i,])
              sm2= sm2+XS2[i,]%*%t(XS2[i,])
            }
            
            for (j in 1:s)
                {
              S1[j]<-t(I[,j])%*% solve(sm1/n)%*%I[,j]
              S2[j]<-t(I[,j])%*% solve(sm2/n)%*%I[,j]
               }
          
          }
          if(metho=="VCPCR"){
#K<-as.numeric(readline("Please enter a value for K, the number of classes: "))
  #### split data ---
            source("~/Desktop/phd /R codes /VC-pcr codes /all VCPCR function.R")
            Cs.init<-GenRandomClusters(K,p)
            Prosfold1<-ProcessFoldData(X1, Y1, test.id=NULL , scale.X = T , scale.y = T , bin.y = F)
            #### transform X to Z  in D1
            VCPCR1<-VCPCR(fold.data=Prosfold1 , lambda=0.5 , delta=0, Cs.init=Cs.init , X.sd = NULL, alpha = 1, bin.y = F)
            V1<-VCPCR1$V
            Z1<-X1%*%V1
    
            ######### couper Y en tranches et estimer la  proportion p chapeau pour chaque tranche sh
            set.seed(123)
            f=dr.slices(Y1,H)    ### tranchage de y 
            a=f$slice.indicator  #### les indices de chaque slice
            teta.hat1<-matrix(nrow=H,ncol=ncol(V1))
            beta.hat1<-matrix(nrow=H,ncol=p)
            
            for(h in 1:H) {   
              indicator <-c(rep(0,n))
              for (i in 1:n) 
              { if (a[i]==h ) { indicator[i]=1 }
                else {i=i+1}
              }
              
              lm<-lm(indicator ~ Z1)
              teta.hat1[h,]=lm$coefficients[-1]
              beta.hat1[h,]<-V1%*%teta.hat1[h,]
              
              h=h+1
            } 
            
            
            ######## estimed Beta ind D2#
          Prosfold2<-ProcessFoldData(X2, Y2, test.id=NULL , scale.X = T , scale.y = T , bin.y = F)
          VCPCR2<-VCPCR(fold.data=Prosfold1 , lambda=0.5 , delta=1, Cs.init=Cs.init , X.sd = NULL, alpha = 0, bin.y = F)
          V2<-VCPCR2$V
          Z2<-X2%*%V2
          set.seed(123)
          f=dr.slices(Y2,H)    ### tranchage de y 
          a=f$slice.indicator    #### les indices de chaque slice 
          teta.hat2<-matrix(nrow=H,ncol=ncol(V2))
          beta.hat2<-matrix(nrow=H,ncol=p)
          
          for(h in 1:H) {   
            indicator <-c(rep(0,n))
            for (i in 1:n) 
            { if (a[i]==h ) { indicator[i]=1 }
              else {i=i+1}
            }
            lm2<-lm(indicator ~ Z2)
            teta.hat2[h,]=lm2$coefficients[-1]
            beta.hat2[h,]<-V2%*%teta.hat2[h,]
            
            h=h+1
          } 
          #### calcul des Sjk 
          E=S2=S1<-vector(mode="numeric",p)
          I<-diag(p)
          s<-K
          sm1<-matrix(0,ncol(Z1),ncol(Z1))
          for (i in 1:n) {
            sm1= sm1+Z1[i,]%*%t(Z1[i,])
          }
          
          sm2<-matrix(0,ncol(Z2),ncol(Z2))
          for (i in 1:n) {
            sm2= sm2+Z2[i,]%*%t(Z2[i,])
          }
          for (j in 1:p)
          {
            S1[j]<-t(I[,j])%*% V1%*%solve(sm1/n)%*%t(V1)%*%I[,j]
            S2[j]<-t(I[,j])%*%V2%*% solve(sm2/n)%*%t(V2)%*%I[,j]
            
          }
          
          
          }
         
          if(metho=="wls.sir"){
            #### Lasso in D1
            ######### couper Y en tranches et estimer la  proportion p chapeau pour chaque tranche sh
            set.seed(123)
            f=dr.slices(Y1,H)    ### tranchage de y 
            a=f$slice.indicator 
            beta.hat1<-matrix(nrow=H,ncol=p)
            for(h in 1:H)   
            {   
              indicator <-c(rep(0,n))
              for (i in 1:n) 
              { if (a[i]==h ) { indicator[i]=1 }
                else {i=i+1}
              }
             wlss1<-wls.sir(x=X1,y=as.vector(indicator),categorical=TRUE, ndim=1 )
              # Reconstruction of the model with the best identified lambda value
              beta.hat1[h,]=wlss1$betahat
              h=h+1
            } 
            #### la liste S of sparsity 
            betaS<-colSums(abs(beta.hat1)) 
            
            S<-which(betaS!=0) #### My S is indice of null coef 
            Sbar<-which(betaS==0)
            ##### estimed Beta ind D2
            set.seed(123)
            f=dr.slices(Y2,H)    ### tranchage de y 
            a=f$slice.indicator    #### les indices de chaque slice 
            beta.hat2<-matrix(nrow=H,ncol=p)
            for(h in 1:H)   
            {   
              indicator <-c(rep(0,n))
              for (i in 1:n) 
              { if (a[i]==h ) { indicator[i]=1 }
                else {i=i+1}
              }
              
              wlss2<-wls.sir(x=X2,y=as.vector(indicator),categorical=TRUE, ndim=1 )
              # Reconstruction of the model with the best identified lambda value
              beta.hat2[h,]=wlss2$betahat
              h=h+1
            } 
            
            for(j in Sbar)   
            {beta.hat2[,j]<-rep(0,H)}
            
            
            
            #### calcul des Sjk 
            S2=S1<-vector(mode="numeric",p)
            s<-length(S)
            I<-diag(s)
            XS1<-X1[,S]
            XS2<-X2[,S]
            sm2=sm1<-matrix(0,s,s)
            for (i in 1:n) {
              sm1= sm1+XS1[i,]%*%t(XS1[i,])
              sm2= sm2+XS2[i,]%*%t(XS2[i,])
            }
            
            for (j in 1:s)
            {
              S1[j]<-t(I[,j])%*% solve(sm1/n)%*%I[,j]
              S2[j]<-t(I[,j])%*% solve(sm2/n)%*%I[,j]
              
            }
            
            
          }
          

          #### calcul of W 
          W<-vector(mode="numeric",s)
          for( j in 1:s) W[j] <- t(beta.hat1[,j])%*%beta.hat2[,j]/ (sqrt(S1[j]*S2[j]))
          ### calcul of Aplus 
          t = sort(abs(W))
          Ta = sapply(t,function(x){(1+length(W[W<=(-x)]))/max(1,length(W[W>=x]))})  #for HD 
        }


   tmin<- min(t[which(Ta<=alf)])
  Aplus<-which(W>=tmin)
  betaselct = rep(0, p)
  betaselct [Aplus] = 1


FDP<-(length(which(betaselct[(q+1):p]!=0)))/max(1,length(Aplus)) #for FDP=E(FDR) 

AP<-(length(which(betaselct[1:q]!=0)))/max(1,q) # AP= E(TDP) 


result = list(betaselct=betaselct , FDP=FDP , AP=AP)
 
return(result)

}



inv_square_root<-function(S){
SVD_sigma <- svd(S)
U_sigma <- SVD_sigma$u
D_sigma <- SVD_sigma$d
square_root_sigma <- U_sigma%*%diag(sqrt(D_sigma))%*%t(U_sigma)
inv_diag <- ifelse(D_sigma<0.000001, 0, 1/sqrt(D_sigma))
inv_square_root_Sigma <- U_sigma%*%diag(inv_diag)%*%t(U_sigma)
inv_square_root_Sigma <- U_sigma%*%diag(1/sqrt(D_sigma))%*%t(U_sigma)
}


suffSIR <- function(X, Y, family = c("gaussian", "binomial"),
                    d = 3, n_lambda = 10, maxnvar = ncol(X),
                    lambda = NULL, lambda_max = NULL,
                    lambda_min = NULL,
                    lambda_seq = c("loglinear","linear"),
                    screening = TRUE){
  
  lambda_seq <- match.arg(lambda_seq, c("loglinear","linear"))
  assertthat::assert_that(is.null(lambda) || all(lambda >= 0),
                          msg = "lambda vector must be non-negative.")
  if (!is.null(lambda_max) && !is.null(lambda_min)) {
    assertthat::assert_that(lambda_min >= 0 && lambda_min < lambda_max,
                            msg = paste("lambda_min must be non-negative and",
                                        "strictly less than lambda_max"))
  }
  n = nrow(X)
  p = ncol(X)
  family <- match.arg(family, c("gaussian", "binomial"))
  n_d <- length(d)
  
  
  # center input data X
  xmeans = colMeans(X)
  xsd = apply(X, 2, sd)
  X = scale(X, center = xmeans, scale = xsd)
  
  # output holder
  #intercept <- matrix(0, nrow = n_lambda, ncol = n_d)
  Betahat <- matrix(0, nrow = p, ncol = n_lambda * n_d)
  Vhat.norm = matrix(0, nrow = p, ncol = n_lambda * n_d)
  
  # generate sample covariance matrix
  S = crossprod(X) / n
  
  if (is.null(lambda)) {
    lambda = compute_lambda(S, maxnvar, n_lambda,
                            lambda_max, lambda_min, lambda_seq)
  } else {
    lambda <- sort(lambda)
  }
  
  for(k in seq_len(n_d)) {
    # run fps
    approximate = fps::fps(S, ndim = d[k], ncomp = d[k],
                           lambda = lambda,
                           maxnvar = maxnvar, maxiter = 100)
    
    for (i in seq_len(n_lambda)) {
      # compute row-wise l2 norm of vhat
      vhat.norm = diag(approximate$projection[[i]])
      Vhat.norm[,(k - 1) * n_lambda - i] <- vhat.norm
      if (screening == TRUE) {
        l = sort(vhat.norm)
        t <- ifelse(i == 1, max(l), findt(l))
        # set rows and columns that has small diagonal values to be 0
        small <- (vhat.norm < t)
        approximate$projection[[i]][small, ] <- 0
        approximate$projection[[i]][ ,small] <- 0
      }
      if (all(approximate$projection[[i]] == 0)) {
        betahat = rep(0, p) # use ymean to predict
        #inter <- 0
      } else {
        # decompose projection matrix
        decomp = svd(approximate$projection[[i]], nu = d[k], nv = 0)
        vhat = decomp$u
        vhat[abs(vhat) < 1e-10] <-  0
        pchat = X %*% vhat
        # fit linear model
        # model = stats::glm(Y ~ pchat, family = family)
        outsir = do.sir(pchat, y)    ### use the function do.sir from the package "rdimtools" we can also use package "dr"
        gamhat<-outsir$projection 
        #inter <- coef(model)[1]
        betahat = vhat %*% gamhat
      }
      # store outputs
      Betahat[,(k - 1) * n_lambda - i] <- betahat
      #intercept[i,k] <- inter
    }
  }
  
  # Betahat <- Matrix::drop0(Betahat)
  #Vhat.norm <- Matrix::drop0(Vhat.norm)
  #nn <- expand.grid(seq_len(n_lambda), d)
  #nn <- stringr::str_glue_data(nn, "lamidx{Var1}_d{Var2}")
  #colnames(Betahat) <- nn
  # colnames(Vhat.norm) <- nn
  
  #out <- list(betahat = Betahat,
  # X = X,
  #vhat_norm = Vhat.norm,
  #    family = family,
  #   xmeans = xmeans,
  #  xsd = xsd,
  # d = d, n_lambda = n_lambda,
  #lambda = lambda)
  return(Betahat)
}


compute_lambda <- function(S, maxnvar, nsol,
                           lambda_max = NULL,
                           lambda_min = NULL,
                           lambda_seq = "loglinear"){
  maxoffdiag <- sort(compute_maxoffdiag(S), decreasing = T)
  if (is.null(lambda_max)) lambda_max <- maxoffdiag[1]
  if (is.null(lambda_min)) lambda_min <- maxoffdiag[maxnvar]
  # lambda_min = 0
  lambda <- switch(
    lambda_seq,
    linear = seq(from = lambda_max, to = lambda_min, length.out = nsol),
    loglinear = lambda_min * log10(seq(1, 10, length.out = nsol)) +
      lambda_max * (1 - log10(seq(1, 10, length.out = nsol)))
  )
  lambda
}

findt <- function(x){
  # this function returns the best threshold t for a sorted x
  # input x is a sorted x in ascending order
  
  if (max(x) == 0) return(0)
  p <- length(x)
  # compute total variance when divide the x in position i
  total_var <- rep(NA, p-1)
  total_var[1] <- (p - 1) * var(x[2:p])
  total_var[p-1] <- (p - 1) * var(x[1:(p - 1)])
  for(i in 2:(p-2)){
    var1 = var(x[1:i])
    var2 = var(x[(i + 1):p])
    total_var[i] = i * var1 + (p - i) * var2
  }
  # compute first derivative of total variance
  first_Derivative <- diff(total_var)
  # find the first point that moving it from the 2nd group to 1 group
  # increase the first derivative
  i <- 2
  index <- 1
  while (i < p - 2) {
    stopper <- first_Derivative[i] - first_Derivative[i-1] >
      mean(abs(first_Derivative[1:(i - 1)]))
    if (stopper) {
      index <- i
      break
    } else {
      i <- i + 1
    }
  }
  return(x[index])
}


# return the maximum off diagnal absolute value in each row
compute_maxoffdiag = function(S) {
  # S is the input matrix which is symmetric
  S = abs(S)
  n = nrow(S)
  out = rep(NA, n)
  out[1] = max(S[1, 2:n])
  for (i in 2:(n-1)) {
    out[i] = max(max(S[i, 1:(i-1)]), max(S[i, (i+1):n]))
  }
  out[n] = max(S[n, 1:(n-1)])
  return(out)
}


# ============================================================
# simple_ranking : classe les variables par corrélation
# marginale avec y (via SIR ou corrélation simple)
# ============================================================
simple_ranking <- function(X, y, H = 10) {
  
  n <- nrow(X)
  p <- ncol(X)
  
  # Découpage de y en H tranches (slices SIR)
  breaks <- quantile(y, probs = seq(0, 1, length.out = H + 1))
  breaks[1]     <- -Inf
  breaks[H + 1] <-  Inf
  slices <- cut(y, breaks = breaks, labels = FALSE, include.lowest = TRUE)
  
  # Moyenne de X dans chaque tranche
  slice_means <- matrix(0, nrow = H, ncol = p)
  for (h in 1:H) {
    idx <- which(slices == h)
    if (length(idx) > 0)
      slice_means[h, ] <- colMeans(X[idx, , drop = FALSE])
  }
  
  # Score SIR marginal pour chaque variable : variance inter-tranches pondérée
  nh     <- table(slices) / n
  grand  <- colMeans(X)
  scores <- numeric(p)
  
  for (j in 1:p) {
    scores[j] <- sum(as.numeric(nh) * (slice_means[, j] - grand[j])^2)
  }
  
  # Rang : rank 1 = variable la moins importante (poids fort dans glmnet)
  #        rank p = variable la plus importante  (poids faible)
  ranks <- rank(scores)   # ordre croissant du score
  
  return(ranks)
}
#wlsFunc
wls.sir<- function(x, y, nslice=10, cn1=0.1, cn2=1, choose.dir=FALSE, categorical=FALSE, ndim=1 ){
  ## require package MASS and dr
  ## basic information about x and y
  n <- dim(x)[1]
  p <- dim(x)[2]
  h <- nslice
  x = scale(x,scale=FALSE)
  
  ## slice
  if( categorical==TRUE ){
    index <- as.numeric(factor(y))
    nh <- summary(factor(y))
    h <- length(nh)
  }
  if( categorical == FALSE ){
    slice <- dr.slices.arc(y,h)
    index <- slice$slice.indicator		# Slice Index
    nh <- slice$slice.sizes				# Observations per Slice
  }
  ph <- nh/n
  
  ####################################
  ## calculate wls
  ####################################
  wls <- c()								# Weighted Leverage Score
  
  ####################################
  ## scenario 1
  ####################################
  svdx <- svd(x)
  u <- svdx$u
  d <- svdx$d
  v <- svdx$v
  
  if( choose.dir==FALSE) {
    dir = min(n, p)
  } else {
    ## selection of d   
    theta = d^2/(d[1])^2 + 1
    loglik <- penalty <- rep(0, length(d) )
    for( i in 1:length(d) ){
      if(i < length(d)) {
        loglik[i] <- sum( log(theta[(i+1):length(d)]) + 1 - theta[(i+1):length(d)] )
      } else loglik[i] = 0
      penalty[i] <- i*cn1/sqrt(n)
    }
    BIC = -loglik + penalty
    (dir <- which.min(BIC))
    print(dir)
  }
  
  ## calculate WLS
  w <- matrix(ncol = dir, nrow=length(nh))
  for(j in 1:dir){
    for(i in 1:length(nh)) w[i,j] <- sum(u[,j] * (index == i))/nh[i]
  }
  
  uut <- array(dim = c(length(nh), dir, dir))			# UUT Array
  for(i in 1:length(nh)) uut[i,,] <- (nh[i]) * ( w[i,] %*% t(w[i,]) )
  
  ## LEVERAGE SCORES
  for(j in 1:p) wls[j] <- t(v[j,1:dir]) %*% colSums(uut) %*% v[j,1:dir]
  
  
  ####################################
  ## BIC
  ####################################
  wls.sort = sort(wls, decreasing = TRUE)
  
  loglik <- penalty <- rep(0, min(n,p))
  for(k in 1:min(n,p)){
    temp_loglik = sum(wls.sort[1:k])
    (loglik[k] = -log(temp_loglik))
    penalty[k] = (log(n) + cn2*log(p))*k/max(n,p)
  }
  BIC <- loglik + penalty
  (sel.k <- which.min(BIC))
  
  select <- order(wls, decreasing = T)[1:sel.k]
  
  z<-x[,select]
  #outsir <- SIR(y, z, H = 10)
  #outsir = do.sir(z, y,ndim=ndim)    ### use the function do.sir from the package "rdimtools" we can also use package "dr"
  outsir = do.rsir(z, y,ndim=ndim,regmethod="Ridge")
  ### use the function do.sir from the package "rdimtools" we can also use package "dr"
  
    #theta<-outsir$b
  betahat<-matrix(0,nrow=ncol(x),ncol=ndim)
  for(j in 1:ndim){
  theta<-Re(outsir$projection[,j])
  betahat[select,j]<-theta 
}
  return( list(wls=wls, select=select,betahat=betahat) )
  
}
#wls.sir<-function(x,y,categorical){
wls<-function(x, y, nslice=10, cn1=0.1, cn2=1, choose.dir=FALSE, categorical){
  ## require package MASS and dr
  ## basic information about x and y
  n <- dim(x)[1]
  p <- dim(x)[2]
  h <- nslice
  x = scale(x,scale=FALSE)
  
  ## slice
  if( categorical==TRUE ){
    index <- as.numeric(factor(y))
    nh <- summary(factor(y))
    h <- length(nh)
  }
  if( categorical == FALSE ){
    slice <- dr.slices.arc(y,h)
    index <- slice$slice.indicator		# Slice Index
    nh <- slice$slice.sizes				# Observations per Slice
  }
  ph <- nh/n
  
  ####################################
  ## calculate wls
  ####################################
  wls <- c()								# Weighted Leverage Score
  
  ####################################
  ## scenario 1
  ####################################
  svdx <- svd(x)
  u <- svdx$u
  d <- svdx$d
  v <- svdx$v
  
  if( choose.dir==FALSE) {
    dir = min(n, p)
  } else {
    ## selection of d   
    theta = d^2/(d[1])^2 + 1
    loglik <- penalty <- rep(0, length(d) )
    for( i in 1:length(d) ){
      if(i < length(d)) {
        loglik[i] <- sum( log(theta[(i+1):length(d)]) + 1 - theta[(i+1):length(d)] )
      } else loglik[i] = 0
      penalty[i] <- i*cn1/sqrt(n)
    }
    BIC = -loglik + penalty
    (dir <- which.min(BIC))
    print(dir)
  }
  
  ## calculate WLS
  w <- matrix(ncol = dir, nrow=length(nh))
  for(j in 1:dir){
    for(i in 1:length(nh)) w[i,j] <- sum(u[,j] * (index == i))/nh[i]
  }
  
  uut <- array(dim = c(length(nh), dir, dir))			# UUT Array
  for(i in 1:length(nh)) uut[i,,] <- (nh[i]) * ( w[i,] %*% t(w[i,]) )
  
  ## LEVERAGE SCORES
  for(j in 1:p) wls[j] <- t(v[j,1:dir]) %*% colSums(uut) %*% v[j,1:dir]
  
  
  ####################################
  ## BIC
  ####################################
  wls.sort = sort(wls, decreasing = TRUE)
  
  loglik <- penalty <- rep(0, min(n,p))
  for(k in 1:min(n,p)){
    temp_loglik = sum(wls.sort[1:k])
    (loglik[k] = -log(temp_loglik))
    penalty[k] = (log(n) + cn2*log(p))*k/max(n,p)
  }
  BIC <- loglik + penalty
  (sel.k <- which.min(BIC))
  
  select <- order(wls, decreasing = T)[1:sel.k]
  

  return(select )
  
}




cen.wls <- function(x, y, delta, c = 0.05, n.slice = 10, h = NULL,
                    cn1 = 0.1, cn2 = 1, choose.dir = FALSE) {
  
  n <- dim(x)[1]; p <- dim(x)[2]
  x <- scale(x, scale = FALSE)
  
  svdx <- svd(x, nu = min(n, p), nv = min(n, p))
  u <- svdx$u; d <- svdx$d; v <- svdx$v   # d = valeurs singulières (papier: λ_i)
  
  if (is.null(h)) {
    beta_init <- prcomp(x, rank. = 1)$rotation[, 1]
    h <- bw.nrd(as.vector(x %*% beta_init))
  }
  
  # ── ÉTAPE 1 : choix de d_hat (Section 4.1, Eq. 9 du papier) ─────────────────
  if (choose.dir == FALSE) {
    d_hat<- min(n, p)
  } else {
    theta <- d^2 / (d[1])^2 + 1
    loglik <- penalty <- rep(0, length(d))
    for (i in 1:length(d)) {
      if (i < length(d)) {
        loglik[i] <- sum(log(theta[(i+1):length(d)]) + 1 - theta[(i+1):length(d)])
      } else loglik[i] <- 0
      penalty[i] <- i * cn1 / sqrt(n)
    }
    BIC_d <- -loglik + penalty # d_hat not dir its de spiekd value not the direction number 
    d_hat <- which.min(BIC_d)
  }
  
  # ── TRONCATURE IMMÉDIATE de u et v à dir colonnes (comme le fait wls()) ─────
  u <- u[, 1:d_hat, drop = FALSE]
  v <- v[, 1:d_hat, drop = FALSE]
  u.bar <- colMeans(u, na.rm = TRUE)
  
  e.hat.u <- matrix(0, nrow = (n.slice + 1), ncol = d_hat)   # dimension réduite
  
  ordr <- order(y, -delta)
  y <- y[ordr]; u <- u[ordr, , drop = FALSE]; delta <- delta[ordr]
  
  kernel.mtr <- kernel.est(u, h)
  s.hat <- cond.sur.est(kernel.mtr, y, c)
  
  t.slice <- slice.time(y, n.slice)
  m <- matrix(0, nrow = n.slice, ncol = d_hat)               # dimension réduite
  p.hat <- rep(1, (n.slice + 1)); P <- rep(1, n.slice)
  nobs <- length(y)
  
  for (j in 1:(n.slice - 1)) {
    weight.vect <- rep(1, nobs)
    idx <- (1:nobs)
    idx.case1 <- idx[y < t.slice[(j+1)] & delta == 0 & s.hat > c]
    idx.case2 <- idx[y >= t.slice[(j+1)]]
    for (i in idx.case1) {
      weight.vect[i] <- weight.est(y[i], t.slice[(j+1)], i, y, delta, kernel.mtr, s.hat)
    }
    p.hat[(j+1)] <- sum(weight.vect[idx.case1], na.rm=TRUE)/nobs +
      sum(weight.vect[idx.case2], na.rm=TRUE)/nobs
    
    if (length(idx.case1) > 1) {
      ave.case1 <- apply(u[idx.case1,,drop=FALSE] * weight.vect[idx.case1], 2, sum, na.rm=TRUE)/nobs
    } else if (length(idx.case1) == 1) {
      ave.case1 <- u[idx.case1,] * weight.vect[idx.case1] / nobs
    } else ave.case1 <- 0
    
    if (length(idx.case2) > 1) {
      ave.case2 <- apply(u[idx.case2,,drop=FALSE], 2, sum, na.rm=TRUE)/nobs
    } else if (length(idx.case2) == 1) {
      ave.case2 <- u[idx.case2,] / nobs
    } else ave.case2 <- 0
    
    e.hat.u[(j+1),] <- ave.case1 + ave.case2
    diff_p <- p.hat[j] - p.hat[(j+1)]
    if (abs(diff_p) < 1e-10) {
      m[j,] <- rep(0,  d_hat)
    } else {
      m[j,] <- (e.hat.u[j,] - e.hat.u[(j+1),]) / diff_p - u.bar
    }
    P[j] <- max(0, p.hat[j] - p.hat[(j+1)])
  }
  
  if (abs(p.hat[n.slice]) < 1e-10) {
    m[n.slice,] <- rep(0,  d_hat)
  } else {
    m[n.slice,] <- e.hat.u[n.slice,] / p.hat[n.slice] - u.bar
  }
  P[n.slice] <- p.hat[n.slice]
  
  mp <- m * sqrt(P)
  if (!is.matrix(mp)) mp <- matrix(mp, nrow = length(P))
  
  # sigma.eta est maintenant d x d (pas min(n,p) x min(n,p))
  sigma.eta <- t(mp) %*% mp
  if (!is.matrix(sigma.eta)) sigma.eta <- matrix(sigma.eta, nrow=
cen.wls <- function(x, y, delta, c = 0.05, n.slice = 10, h = NULL,
                    cn1 = 0.1, cn2 = 1, choose.dir = FALSE) {
  
  n <- dim(x)[1]; p <- dim(x)[2]
  x <- scale(x, scale = FALSE)
  
  svdx <- svd(x, nu = min(n, p), nv = min(n, p))
  u <- svdx$u; d <- svdx$d; v <- svdx$v   # d = valeurs singulières (papier: λ_i)
  
  if (is.null(h)) {
    beta_init <- prcomp(x, rank. = 1)$rotation[, 1]
    h <- bw.nrd(as.vector(x %*% beta_init))
  }
  
  # ── ÉTAPE 1 : choix de d_hat (Section 4.1, Eq. 9 du papier) ─────────────────
  if (choose.dir == FALSE) {
    d_hat<- min(n, p)
  } else {
    theta <- d^2 / (d[1])^2 + 1
    loglik <- penalty <- rep(0, length(d))
    for (i in 1:length(d)) {
      if (i < length(d)) {
        loglik[i] <- sum(log(theta[(i+1):length(d)]) + 1 - theta[(i+1):length(d)])
      } else loglik[i] <- 0
      penalty[i] <- i * cn1 / sqrt(n)
    }
    BIC_d <- -loglik + penalty # d_hat not dir its de spiekd value not the direction number 
    d_hat <- which.min(BIC_d)
  }
  
  # ── TRONCATURE IMMÉDIATE de u et v à dir colonnes (comme le fait wls()) ─────
  u <- u[, 1:d_hat, drop = FALSE]
  v <- v[, 1:d_hat, drop = FALSE]
  u.bar <- colMeans(u, na.rm = TRUE)
  
  e.hat.u <- matrix(0, nrow = (n.slice + 1), ncol = d_hat)   # dimension réduite
  
  ordr <- order(y, -delta)
  y <- y[ordr]; u <- u[ordr, , drop = FALSE]; delta <- delta[ordr]
  
  kernel.mtr <- kernel.est(u, h)
  s.hat <- cond.sur.est(kernel.mtr, y, c)
  
  t.slice <- slice.time(y, n.slice)
  m <- matrix(0, nrow = n.slice, ncol = d_hat)               # dimension réduite
  p.hat <- rep(1, (n.slice + 1)); P <- rep(1, n.slice)
  nobs <- length(y)
  
  for (j in 1:(n.slice - 1)) {
    weight.vect <- rep(1, nobs)
    idx <- (1:nobs)
    idx.case1 <- idx[y < t.slice[(j+1)] & delta == 0 & s.hat > c]
    idx.case2 <- idx[y >= t.slice[(j+1)]]
    for (i in idx.case1) {
      weight.vect[i] <- weight.est(y[i], t.slice[(j+1)], i, y, delta, kernel.mtr, s.hat)
    }
    p.hat[(j+1)] <- sum(weight.vect[idx.case1], na.rm=TRUE)/nobs +
      sum(weight.vect[idx.case2], na.rm=TRUE)/nobs
    
    if (length(idx.case1) > 1) {
      ave.case1 <- apply(u[idx.case1,,drop=FALSE] * weight.vect[idx.case1], 2, sum, na.rm=TRUE)/nobs
    } else if (length(idx.case1) == 1) {
      ave.case1 <- u[idx.case1,] * weight.vect[idx.case1] / nobs
    } else ave.case1 <- 0
    
    if (length(idx.case2) > 1) {
      ave.case2 <- apply(u[idx.case2,,drop=FALSE], 2, sum, na.rm=TRUE)/nobs
    } else if (length(idx.case2) == 1) {
      ave.case2 <- u[idx.case2,] / nobs
    } else ave.case2 <- 0
    
    e.hat.u[(j+1),] <- ave.case1 + ave.case2
    diff_p <- p.hat[j] - p.hat[(j+1)]
    if (abs(diff_p) < 1e-10) {
      m[j,] <- rep(0,  d_hat)
    } else {
      m[j,] <- (e.hat.u[j,] - e.hat.u[(j+1),]) / diff_p - u.bar
    }
    P[j] <- max(0, p.hat[j] - p.hat[(j+1)])
  }
  
  if (abs(p.hat[n.slice]) < 1e-10) {
    m[n.slice,] <- rep(0,  d_hat)
  } else {
    m[n.slice,] <- e.hat.u[n.slice,] / p.hat[n.slice] - u.bar
  }
  P[n.slice] <- p.hat[n.slice]
  
  mp <- m * sqrt(P)
  if (!is.matrix(mp)) mp <- matrix(mp, nrow = length(P))
  
  # sigma.eta est maintenant d x d (pas min(n,p) x min(n,p))
  sigma.eta <- t(mp) %*% mp
  if (!is.matrix(sigma.eta)) sigma.eta <- matrix(sigma.eta, nrow=d_hat, ncol=d_hat)
  
  # ── ÉTAPE 2 : calcul de omega_j sur les dir directions retenues ─────────────
  wls <- numeric(p)
  for (j in 1:p) {
    v_j <- matrix(v[j, ], ncol = 1)   # v déjà tronqué à dir colonnes
    wls[j] <- as.numeric(t(v_j) %*% sigma.eta %*% v_j)
  }
  
  # ── ÉTAPE 3 : choix de p0_hat (Section 4.2, Eq. 10 du papier) ───────────────
  wls.sort <- sort(wls, decreasing = TRUE)
  loglik <- penalty <- rep(0, min(n, p))
  for (k in 1:min(n, p)) {
    temp_loglik <- sum(wls.sort[1:k])
    loglik[k] <- -log(temp_loglik)
    penalty[k] <- (log(n) + cn2 * log(p)) * k / max(n, p)
  }
  BIC <- loglik + penalty
  sel.k <- which.min(BIC)
  select <- order(wls, decreasing = TRUE)[1:sel.k]
  
  return(list(n.sel = sel.k, select = select, dir_used =  d_hat))
}, ncol=d_hat)
  
  # ── ÉTAPE 2 : calcul de omega_j sur les dir directions retenues ─────────────
  wls <- numeric(p)
  for (j in 1:p) {
    v_j <- matrix(v[j, ], ncol = 1)   # v déjà tronqué à dir colonnes
    wls[j] <- as.numeric(t(v_j) %*% sigma.eta %*% v_j)
  }
  
  # ── ÉTAPE 3 : choix de p0_hat (Section 4.2, Eq. 10 du papier) ───────────────
  wls.sort <- sort(wls, decreasing = TRUE)
  loglik <- penalty <- rep(0, min(n, p))
  for (k in 1:min(n, p)) {
    temp_loglik <- sum(wls.sort[1:k])
    loglik[k] <- -log(temp_loglik)
    penalty[k] <- (log(n) + cn2 * log(p)) * k / max(n, p)
  }
  BIC <- loglik + penalty
  sel.k <- which.min(BIC)
  select <- order(wls, decreasing = TRUE)[1:sel.k]
  
  return(list(n.sel = sel.k, select = select, dir_used =  d_hat))
}

computeProjection <- function(one_beta) {
  # compute the projection matrix corresponding to one vector or one matrix
  one_beta = as.matrix(one_beta)
  a = min(sqrt(apply(one_beta,2,function(j) sum(j^2))))
  if (rcond(one_beta)<1e-10) 
    return(NA)
  else {
    # None of the columns equal to 0
    Pr = one_beta %*% solve(crossprod(one_beta)) %*% t(one_beta)
    return(Pr)
  }
}

getKernelMatrix <- function(X, y, method = c("sir", "save", "phdy"), nslices = 10){
  # get Kernel Matrix that corresponds to inverse regression methods
  n <- length(y)
  p <- ncol(X)
  order_y <- order(y) 
  slice_size <- floor(n/nslices);
  
  # Divide the feature matrix into slices
  slicedData <- array(0, dim = c(nslices, slice_size, p))
  for (h in 1:nslices){
    begin_index = slice_size*(h-1)+1
    end_index = slice_size*h
    slicedData[h, , ] = X[order_y[begin_index:end_index], ] 
  }
  
  M <- list()
  
  if ("sir" %in% method){
    tmp <- apply(slicedData, c(1, 3), mean, na.rm = TRUE)
    Msir <- cov(tmp)
    M <- c(M, Msir = list(Msir))
    rm(tmp)
  } 
  
  if ("save" %in% method){
    tmp <- lapply(1:nslices, function(oneslice) {
      Q <- cov(X) - cov(slicedData[oneslice, ,])
      Q %*% Q
    })
    Msave <- Reduce("+", tmp)/nslices
    M <- c(M, Msave = list(Msave))
    rm(tmp)
  } 
  
  if ('phdy' %in% method){
    y = y - mean(y)
    Mphdy <- t(X*y) %*% X/n
    M <- c(M, Mphdy = list(Mphdy))
  } 
  return(M)
}

ICProjection = function(fit, beta_ref_PIC, tau){
  # compute the PIC. 
  # fit: CHOMP/Adaptive CHOMP fit object, on which we use PIC to select the tuning parameter and final estimator
  # beta_ref_PIC: an estimator, typically the unpenalized estimator, that is consistent
  # tau: model complexity term
  
  p <- nrow(fit$beta)
  G = beta_ref_PIC %*% solve(crossprod(beta_ref_PIC)) %*% t(beta_ref_PIC)
  
  frobnorm <- apply(fit$beta, 2, function(x){
    tmp <- computeProjection(x)
    if (is.na(sum(tmp))) return(NA) else return(norm(tmp - G, "F")^2)
  })
  
  maxbetafit <- apply(abs(fit$beta),2,max)
  ICseq <- p*frobnorm + tau*fit$df
  ICseq[fit$df == 0 | maxbetafit < 1e-10] <- max(ICseq) 
  betaIC1 = coef(fit, s = fit$lambda[which.min(ICseq)])[-1]
  
  return(list(betaIC1 = betaIC1,
              IC = min(ICseq), 
              frobnorm = frobnorm,
              ICseq = ICseq))
}

AdaptCHOMP.fit = function(R, kappa, beta_ref_weight = NULL, beta_ref_PIC, tau){
  # R = t(L) transpose of the Cholesky factor; 
  # kappa: pseudo-response computed corresponding to each inverse regression method
  # beta_ref_weight: the initial estimator that is used for defining adaptive weight; if not supplied, then (unadaptive) CHOMP is fit
  # beta_ref_PIC: the initial estimator that is used to compute PIC 
  # tau: penalty for PIC
  
  p <- dim(R)[1]
  R <- as.matrix(R)
  kappa <- as.numeric(kappa)
  
  if (is.null(beta_ref_weight)){
    beta_ref_weight <- rep(1, p)
  }
  
  beta_ref_weight[beta_ref_weight == 0] <- 1e-8
  m_adapt_old = glmnet(x = R, y = kappa, penalty.factor = 1/abs(beta_ref_weight),
                       intercept = FALSE)
  max_lambda = m_adapt_old$lambda[1]
  # Grid of tuning parameters 
  lambda_seq = exp(seq(max(log(max_lambda),1e-6), log(max(1e-5*max_lambda, 1e-5)),
                       len=100))
  # Fit the adaptive CHOMP estimates
  m_adapt <- glmnet(x = R, y = kappa, penalty.factor = 1/abs(beta_ref_weight),
                    intercept = FALSE,
                    lambda = lambda_seq)
  
  ICnew <- ICProjection(fit = m_adapt, beta_ref_PIC = beta_ref_PIC, tau = tau)
  betaICnew <- ICnew$betaIC1
  betas <- betaICnew
  return(list(betas = betas, IC = ICnew$IC))
}

AdaptCHOMPwithPIC <- function(X, y, method = "sir", nslices = 10, d, gamma.pow = 2, 
                              adaptive = TRUE){
  # X : a n \times p design matrix; y : n \times 1 outcome
  # method: base method, either "sir", "save", or "phd"
  # nslices: require for sir and save
  # d: number of dimensions
  # gamma.pow: the power for defining the Adaptive CHOMP weight estimator, default is 2
  # adaptive: whether adaptive CHOMP should be fitted, if FALSE, unweighted CHOMP is fit
  
  SIR.kernel <- getKernelMatrix(X, y, method = method, nslices = nslices) 
  eta <- eigen(SIR.kernel$Msir)$vectors[, 1:d, drop = FALSE]
  unpen.fit <- solve(cov(X), eta)
  
  # ChOMP and Adaptive CHOMP fit
  R <- chol(cov(X)) # R is an upper triangular matrix
  L <- t(R)
  kappa <- forwardsolve(L, eta)
  # We fit CHOMP for each dimension separately
  if (adaptive){
    fit <- sapply(1:d, function(one_dim){
      AdaptCHOMP.fit(R, kappa[, one_dim, drop = FALSE], beta_ref_weight = NULL, beta_ref_PIC = unpen.fit[, one_dim, drop = FALSE], tau = log(p))
    })
  }
  else{
    fit <- sapply(1:d, function(one_dim){
      AdaptCHOMP.fit(R, kappa[, one_dim, drop = FALSE], beta_ref_weight = unpen.fit[, one_dim, drop = FALSE]^gamma.pow, beta_ref_PIC = unpen.fit[, one_dim, drop = FALSE], tau = log(p))  
    })
  }
  # Adaptive CHOMP with weight defined based on SIR.fit  
  
  estBeta <- as(Reduce(cbind, fit['betas', ]), 'sparseMatrix')
  colnames(estBeta) <- paste("Dim", 1:d)
  return(estBeta)
}


#### wls.SIR ------
wls.sir<- function(x, y, nslice=10, cn1=0.1, cn2=1, choose.dir=FALSE, categorical=FALSE, ndim=1 ){
  ## require package MASS and dr
  ## basic information about x and y
  n <- dim(x)[1]
  p <- dim(x)[2]
  h <- nslice
  x = scale(x,scale=FALSE)
  
  ## slice
  if( categorical==TRUE ){
    index <- as.numeric(factor(y))
    nh <- summary(factor(y))
    h <- length(nh)
  }
  if( categorical == FALSE ){
    slice <- dr.slices.arc(y,h)
    index <- slice$slice.indicator		# Slice Index
    nh <- slice$slice.sizes				# Observations per Slice
  }
  ph <- nh/n
  
  ####################################
  ## calculate wls
  ####################################
  wls <- c()								# Weighted Leverage Score
  
  ####################################
  ## scenario 1
  ####################################
  svdx <- svd(x)
  u <- svdx$u
  d <- svdx$d
  v <- svdx$v
  
  if( choose.dir==FALSE) {
    dir = min(n, p)
  } else {
    ## selection of d   
    theta = d^2/(d[1])^2 + 1
    loglik <- penalty <- rep(0, length(d) )
    for( i in 1:length(d) ){
      if(i < length(d)) {
        loglik[i] <- sum( log(theta[(i+1):length(d)]) + 1 - theta[(i+1):length(d)] )
      } else loglik[i] = 0
      penalty[i] <- i*cn1/sqrt(n)
    }
    BIC = -loglik + penalty
    (dir <- which.min(BIC))
    print(dir)
  }
  
  ## calculate WLS
  w <- matrix(ncol = dir, nrow=length(nh))
  for(j in 1:dir){
    for(i in 1:length(nh)) w[i,j] <- sum(u[,j] * (index == i))/nh[i]
  }
  
  uut <- array(dim = c(length(nh), dir, dir))			# UUT Array
  for(i in 1:length(nh)) uut[i,,] <- (nh[i]) * ( w[i,] %*% t(w[i,]) )
  
  ## LEVERAGE SCORES
  for(j in 1:p) wls[j] <- t(v[j,1:dir]) %*% colSums(uut) %*% v[j,1:dir]
  
  
  ####################################
  ## BIC
  ####################################
  wls.sort = sort(wls, decreasing = TRUE)
  
  loglik <- penalty <- rep(0, min(n,p))
  for(k in 1:min(n,p)){
    temp_loglik = sum(wls.sort[1:k])
    (loglik[k] = -log(temp_loglik))
    penalty[k] = (log(n) + cn2*log(p))*k/max(n,p)
  }
  BIC <- loglik + penalty
  (sel.k <- which.min(BIC))
  
  select <- order(wls, decreasing = T)[1:sel.k]
  
  z<-x[,select]
  #outsir <- SIR(y, z, H = 10)
  #outsir = do.sir(z, y,ndim=ndim)    ### use the function do.sir from the package "rdimtools" we can also use package "dr"
  outsir = do.rsir(z, y,ndim=ndim,regmethod="Ridge")
  ### use the function do.sir from the package "rdimtools" we can also use package "dr"
  
  #theta<-outsir$b
  betahat<-matrix(0,nrow=ncol(x),ncol=ndim)
  for(j in 1:ndim){
    theta<-Re(outsir$projection[,j])
    betahat[select,j]<-theta 
  }
  return( list(wls=wls, select=select,betahat=betahat) )
  
}




# compute the unpenalized SIR, assuming d = 2 being known
"cov1" <-
  function(x, use="pair"){
    n = nrow(x);
    c = cov(x, use=use)*(n-1)/n;
    c
  }

"kernel.est" <-
  function(Bx,h){
    n = nrow(Bx);
    r = ncol(Bx);
    b = rep(h, r);
    
    kernel.x = matrix(rep(0,n*n), nrow=n);
    for (i in 1:n){
      Bx.temp = t(matrix(rep(Bx[i,],n), nrow=r, ncol=n));
      d2 = ((Bx - Bx.temp)/b)^2;
      ds = apply(d2, 1, sum);
      fh = exp(-0.5*(ds));
      kernel.x[i,] = fh;
    }
    kernel.x
  }
"cond.sur.est" <-
  function(kernel.mtr, y, c){
    n = length(y);
    s.hat = rep(0, n)
    for (i in 1:n){
      yi = (1:n)[y>y[i]];
      fh = kernel.mtr[i,];
      if(length(yi)>0){
        s.hat[i] = sum(fh[yi], na.rm=TRUE)/sum(fh, na.rm=TRUE);
        s.hat[i] = max(s.hat[i], c, na.rm=TRUE);
      }else{
        s.hat[i] = c;
      }
    }
    s.hat
  }

"weight.est" <-
  function(t1, t2, index.x, y, delta, kernel.mtr, s.hat){
    if (t1 >= t2){
      cat(paste(t1, "\t", t2, "\n"));
      stop("t1 >= t2!");
    }
    
    # equation (5.2)
    fhatx0 = mean(kernel.mtr[index.x,], na.rm=TRUE)
    lambdahat.ttx.temp = kernel.mtr[index.x,] * (t1<y)*(y<t2) *  delta 
    lambdahat.ttx = mean((lambdahat.ttx.temp / s.hat), na.rm=TRUE) / fhatx0
    
    # equation (5.1)
    weight.ttx = exp(-lambdahat.ttx)
    weight.ttx
  }

"slice.time" <-
  function(t, h){
    lt=length(t);
    n1=ceiling(lt/h);  
    n2=floor(lt/h);  
    if(n1==n2){
      n = lt/h;
      index = rep(n, h);
    }else if(n1 > n2){
      h1 = lt%%h;
      h2 = h - h1;
      index = rep(c(n1, n2), c(h1, h2));
    }
    t.slice=numeric(h+1);  # (h+1) time points
    t.slice[1]=min(t) - 0.000001;
    end = 0;
    for(i in 1:h){
      bgn = end + 1;
      end = end + index[i];
      t.slice[i+1] = t[end];
    }
    t.slice 
  }

#### wls.sir with FDR control : -------
wls.sir.fdr<- function(x, y, nslice=10, cn1=0.1, cn2=1, choose.dir=FALSE, categorical=FALSE, ndim=1 ){
  ## require package MASS and dr
  ## basic information about x and y
  n <- dim(x)[1]
  p <- dim(x)[2]
  h <- nslice
  x = scale(x,scale=FALSE)
  
  ## slice
  if( categorical==TRUE ){
    index <- as.numeric(factor(y))
    nh <- summary(factor(y))
    h <- length(nh)
  }
  if( categorical == FALSE ){
    slice <- dr.slices.arc(y,h)
    index <- slice$slice.indicator		# Slice Index
    nh <- slice$slice.sizes				# Observations per Slice
  }
  ph <- nh/n
  
  ####################################
  ## calculate wls
  ####################################
  wls <- c()								# Weighted Leverage Score
  
  ####################################
  ## scenario 1
  ####################################
  svdx <- svd(x)
  u <- svdx$u
  d <- svdx$d
  v <- svdx$v
  
  if( choose.dir==FALSE) {
    dir = min(n, p)
  } else {
    ## selection of d   
    theta = d^2/(d[1])^2 + 1
    loglik <- penalty <- rep(0, length(d) )
    for( i in 1:length(d) ){
      if(i < length(d)) {
        loglik[i] <- sum( log(theta[(i+1):length(d)]) + 1 - theta[(i+1):length(d)] )
      } else loglik[i] = 0
      penalty[i] <- i*cn1/sqrt(n)
    }
    BIC = -loglik + penalty
    (dir <- which.min(BIC))
    print(dir)
  }
  
  ## calculate WLS
  w <- matrix(ncol = dir, nrow=length(nh))
  for(j in 1:dir){
    for(i in 1:length(nh)) w[i,j] <- sum(u[,j] * (index == i))/nh[i]
  }
  
  uut <- array(dim = c(length(nh), dir, dir))			# UUT Array
  for(i in 1:length(nh)) uut[i,,] <- (nh[i]) * ( w[i,] %*% t(w[i,]) )
  
  ## LEVERAGE SCORES
  for(j in 1:p) wls[j] <- t(v[j,1:dir]) %*% colSums(uut) %*% v[j,1:dir]
  
  
  ####################################
  ## BIC
  ####################################
  wls.sort = sort(wls, decreasing = TRUE)
  
  loglik <- penalty <- rep(0, min(n,p))
  for(k in 1:min(n,p)){
    temp_loglik = sum(wls.sort[1:k])
    (loglik[k] = -log(temp_loglik))
    penalty[k] = (log(n) + cn2*log(p))*k/max(n,p)
  }
  BIC <- loglik + penalty
  (sel.k <- which.min(BIC))
  
  select <- order(wls, decreasing = T)[1:sel.k]
  
  z<-x[,select]
  #outsir <- SIR(y, z, H = 10)
  #outsir = do.sir(z, y,ndim=ndim)    ### use the function do.sir from the package "rdimtools" we can also use package "dr"
  outsir = do.rsir(z, y,ndim=ndim,regmethod="Ridge")
  ### use the function do.sir from the package "rdimtools" we can also use package "dr"
  
  #theta<-outsir$b
  betahat<-matrix(0,nrow=ncol(x),ncol=ndim)
  for(j in 1:ndim){
    theta<-Re(outsir$projection[,j])
    betahat[select,j]<-theta 
  }
  
  
  tmin<- min(t[which(Ta<=alf)])
  Aplus<-which(W>=tmin)
  betaselct = rep(0, p)
  betaselct [Aplus] = 1
  
  
  FDP<-(length(which(betaselct[(q+1):p]!=0)))/max(1,length(Aplus)) #for FDP=E(FDR) 
  
  AP<-(length(which(betaselct[1:q]!=0)))/max(1,q) # AP= E(TDP) 
  return( list(wls=wls, select=select,betahat=betahat) )
  
}

##@@@ those function are from cenSIR method : --------
"cen.sir" <-
  function(y, delta, x, n.slice, joint.edrs, h, c){
    # order the data by survival time first
    ordr = order(y, -delta);
    y = y[ordr];
    x = x[ordr,];
    delta = delta[ordr];
    
    # first, estimate the kernel matrix and conditional survial function
    Bx = x %*% joint.edrs;
    kernel.mtr = kernel.est(Bx, h);
    s.hat = cond.sur.est(kernel.mtr, y, c);
    
    # sencondly, calculate the slice mean m_h
    t.slice = slice.time(y, n.slice);
    m = matrix(0, nrow=n.slice, ncol=ncol(x));
    p.hat = rep(1,(n.slice+1));
    p = rep(1,n.slice);
    e.hat.x = matrix(0, nrow=(n.slice+1), ncol=ncol(x));
    x.bar = apply(x, 2, mean, na.rm=TRUE);
    e.hat.x[1,] = x.bar;
    nobs = length(y);
    
    for(j in 1:(n.slice-1)){
      weight.vect = rep(1,nobs);
      idx = (1:nobs);
      idx.case1 = idx[y<t.slice[(j+1)] & delta==0 & s.hat>c];
      idx.case2 = idx[y>=t.slice[(j+1)]];
      
      for(i in idx.case1){
        weight.vect[i] = weight.est(y[i], t.slice[(j+1)], i, y, delta, 
                                    kernel.mtr, s.hat);
      }
      
      p.hat[(j+1)] = sum(weight.vect[idx.case1], na.rm=TRUE)/nobs;
      p.hat[(j+1)] = p.hat[(j+1)]+sum(weight.vect[idx.case2],na.rm=TRUE)/nobs;
      
      if(length(idx.case1)>1){
        x.case1 = x[idx.case1,]*weight.vect[idx.case1]
        ave.case1 = apply(x.case1, 2, sum, na.rm=TRUE)/nobs;
      }else if(length(idx.case1)==1){
        ave.case1 = x[idx.case1,]*weight.vect[idx.case1]/nobs;
      }else if(length(idx.case1)==0){
        ave.case1 = 0;
      }
      
      if(length(idx.case2)>1){
        ave.case2 = apply(x[idx.case2,], 2, sum, na.rm=TRUE)/nobs;
      }else if(length(idx.case2)==1){
        ave.case2 = x[idx.case2,]/nobs;
      }else if(length(idx.case2)==0){
        ave.case2 = 0;
      }
      
      e.hat.x[(j+1),] = ave.case1 + ave.case2;
      
      m[j,] = (e.hat.x[j,]-e.hat.x[(j+1),])/(p.hat[j]-p.hat[(j+1)]) - x.bar;
      p[j]  = max(0,p.hat[j] - p.hat[(j+1)]);
    }
    
    m[n.slice,] = e.hat.x[n.slice,]/p.hat[n.slice] - x.bar;
    p[n.slice]  = p.hat[n.slice];
    
    # at last, compute the SIR direction
    mp = m*sqrt(p);
    sigma.eta =  t(mp) %*% mp;#crossprod(mp)
    sigma.x   = cov1(x, use="pair");
    eign = eigen.decomp(sigma.eta, sigma.x);
    ordr  = order(eign$values, decreasing=TRUE);
    e.val = eign$values[ordr];
    e.vec = eign$vectors[,ordr];
    
    if(is.null(colnames(x))){
      x.names = paste("X", 1:ncol(x), sep="");
    }else{
      x.names = colnames(x);
    }
    dimnames(e.vec) = list(x.names, paste("Dir", 1:ncol(e.vec), sep=""));
    
    csir = list(y=y, x=x, delta=delta, eval=e.val, evec=e.vec, nslice=n.slice);
    csir$c = c;
    csir$h = h;
    class(csir) = "csir";
    csir
  }

"eigen.decomp" <-
  function(m1, m2, tol=1e-5){
    # cholesky decomposition of m2 = t(cm2) %*% cm2
    # and inverse(t(cm2)) = t(inverse(cm2)) 
    cm2 = chol(m2); 
    inv.cm2 = cm2;
    idx.cm2 = diag(cm2)>tol;
    inv.cm2[idx.cm2, idx.cm2] = solve(cm2[idx.cm2, idx.cm2]);
    m3 = t(inv.cm2) %*% m1 %*% (inv.cm2);
    svd.m3 = svd(m3);
    eigen.value = svd.m3$d;
    eigen.vect = inv.cm2 %*% svd.m3$v;
    result = list(values=eigen.value, vectors=eigen.vect);
    result
  }

"print.ds" <-
  function(object, digits = max(3, getOption("digits") - 3), ...)
  {
    cat("Double Slicing of Survival time\n\n")
    cat(paste("Number of observation:           ", length(object$y), "\n"))
    n.censor = length(object$y[object$delta==0]);
    cat(paste("Number of censored observation:  ", n.censor, "\n"))
    cat(paste("Number of predictors:            ", ncol(object$x), "\n"))
    cat(paste("Number of slices (uncensored):   ", object$nslice1, "\n"))
    cat(paste("Number of slices (censored):     ", object$nslice0, "\n\n"))
    cat("Eigenvalues:\n")
    cum.e = cumsum(object$eval/sum(object$eval))
    evals = rbind(object$eval, cum.e)
    rownames(evals) = c("Eigenvalues", "Cum.Sum.R^2");
    colnames(evals) = colnames(object$evec);
    print(evals, digits=digits)
    cat("\nAsym Chi-square test of SIR directions:\n")
    print(sir.test(object), digits=digits)
    cat("\nEigenvectors:\n")
    print(object$evec, digits=digits)
    invisible(object)
  }

"print.csir" <-
  function(object, digits = max(3, getOption("digits") - 3), ...)
  {
    cat("Censored Sliced Inverse Regression Model\n\n")
    cat(paste("Number of observation:           ", length(object$y), "\n"))
    n.censor = length(object$y[object$delta==0]);
    cat(paste("Number of censored observation:  ", n.censor, "\n"))
    cat(paste("Number of predictors:            ", ncol(object$x), "\n"))
    cat(paste("Number of slices:                ", (object$nslice), "\n"))
    cat(paste("Kernel width:                    ", (object$h), "\n\n"))
    cum.e = cumsum(object$eval/sum(object$eval))
    evals = rbind(object$eval, cum.e)
    rownames(evals) = c("Eigenvalues", "Cum.Sum.R^2");
    colnames(evals) = colnames(object$evec);
    cat("Eigenvalues:\n")
    print(evals, digits=digits)
    cat("\nAsym Chi-square test of SIR directions:\n")
    print(sir.test(object), digits=digits)
    cat("\nEigenvectors:\n")
    print(object$evec, digits=digits)
    invisible(object)
  }
"plot.sir" <-
  function(object, which=1:ncol(object$evec), logY=FALSE) {
    if(length(which)>6){
      stop("no more than 6 directions can be plot at one time");
    }
    
    dirs = object$x %*% object$evec[,which];
    if(logY){
      surv = log(object$y);
      labY = "log(survival time)";
    }else{
      surv = object$y;
      labY = "survival time";
    }
    labX = colnames(object$evec)[which];
    
    if(length(which)==2){ par(mfrow=c(1,2)); }
    if(length(which)==3 || length(which)==4){ par(mfrow=c(2,2)); }
    if(length(which)==5 || length(which)==6){ par(mfrow=c(2,3)); }
    
    for(i in 1:length(which)){
      plot(dirs[,i][object$delta==1], surv[object$delta==1], 
           col="blue", type="p", pch=16, 
           xlab=labX[i], ylab=labY, main="", 
           xlim=range(dirs[,i]), ylim=range(surv))
      points(dirs[,i][object$delta==0], surv[object$delta==0], 
             col="red", type="p", pch=10);
      legend(min(dirs[,i]) + 0.05*(max(dirs[,i]) - min(dirs[,i])), 
             min(surv) + 0.95*(max(surv) - min(surv)), 
             legend=c("censored", "uncensored"),  
             col=c("red", "blue"), pch=c(10, 16))
    }
  }

"plot.csir" <-
  function(object, which=1:ncol(object$evec), logY=FALSE) {
    plot.sir(object, which=which, logY=logY)
  }

"plot.ds" <-
  function(object, which=1:ncol(object$evec), logY=FALSE) {
    plot.sir(object, which=which, logY=logY)
  }

"edr.n" <-
  function(object, n){
    object$evec[,1:n]
  }

"edr.percent" <-
  function(object, percent){
    if(percent > 1 || percent < 0){
      stop("percent should be between 0 and 1\n");
    }
    
    values  = object$eval;
    vectors = object$evec;
    
    j=1;
    while(sum(values[1:j])/sum(values) < percent){
      j=j+1; # take the first j eigenvalues
    }
    pct = round(sum(values[1:j])/sum(values), 2);
    cat(paste("eigen.value =", "\n"));
    cat(signif(values,3));
    cat("\n");
    cat(paste("First", j, "eigen values explain",  pct, "variance", "\n"));
    cat("\n");
    
    vectors.mtr = vectors[,1:j];
    vectors.mtr
  }


"plot.3d.sir" <-
  function(object, which=1:2, angles=c(60, 120), z.plane=NULL,logY=FALSE) {
    library(scatterplot3d)
    
    if(length(which)!=2){
      stop("Please indicate 2 and only 2 directions for plot.3d\n");
    }
    if(length(angles)>6){
      stop("no more than 6 directions can be plot at one time");
    }
    
    dirs = object$x %*% object$evec[,which];
    if(logY){
      surv = log(object$y);
      labZ = "log(survival time)";
    }else{
      surv = object$y;
      labZ = "survival time";
    }
    labX = colnames(object$evec[,which]);
    
    if(length(angles)==2){ par(mfrow=c(1,2)); }
    if(length(angles)==3 || length(which)==4){ par(mfrow=c(2,2)); }
    if(length(angles)==5 || length(which)==6){ par(mfrow=c(2,3)); }
    
    d = object$delta;
    
    for(i in 1:length(angles)){
      s3d = scatterplot3d(dirs[,1], dirs[,2], xlab=labX[1], grid=TRUE, 
                          ylab=labX[2], zlab=labZ, surv, angle=angles[i], type="n");
      s3d$points3d(dirs[d==1,1], dirs[d==1,2], surv[d==1], 
                   col="blue", type="p", pch=16);
      s3d$points3d(dirs[d==0,1], dirs[d==0,2], surv[d==0], 
                   col="red", type="p", pch=10);
      if(!is.null(z.plane)){
        s3d$plane3d(z.plane, 0 ,0);
      }
    }
  }


"sir.test" <-
  function(object, nd=length(object$eval)) {
    if(class(object) != "csir" && class(object) != "ds") {
      stop("wrong class of object in sir.test\n"); 
    }
    e = sort(object$eval);
    p = length(object$eval);
    n = length(object$y);
    l = min(p, nd)
    chisq = numeric(l);
    df    = numeric(l);
    p.val = numeric(l);
    for (i in 0:(l-1)){
      j = i+1;
      chisq[j] = n*sum(e[1:(p-i)])
      df[j] = (p-i)*(object$nslice-i-1)
      p.val[j] = 1-pchisq(chisq[j], df[j])
    }
    test = data.frame(Chisq=chisq, df=df, p.value=p.val);
    rnames = character(l);
    for(i in 1:l){
      rnames[i] = paste("D=", i-1, " vs. ", "D>=", i, sep="");
    }
    rownames(test) = rnames;
    test
  }


"double.slice" <-
  function(y, delta, x, n.slice1, n.slice0){
    "cov1" <-
      function(x, use="pair"){
        n = nrow(x);
        c = cov(x, use=use)*(n-1)/n;
        c
      }
    "slice.idx" <-
      function(y, h, m){
        ly = length(y);
        if(m < max(y) ) {stop("y is not a valid index vector");}
        n1 = ceiling(ly/h);
        n2 = floor(ly/h);
        split = numeric(m);
        if(n1==n2){
          n = ly/h;
          index = rep(n, h);
        }else if(n1 > n2){
          h1 = ly%%h;
          h2 = h - h1;
          index = rep(c(n1, n2), c(h1, h2));
        }
        end = 0;
        for(i in 1:h){
          bgn = end + 1;
          end = end + index[i];
          split[y[bgn:end]] = i;
        }
        split
      }
    
    # order the data by survival time first
    ordr = order(y, -delta);
    y = y[ordr];
    x = x[ordr,];
    delta = delta[ordr];
    
    ds.idx = numeric(nrow(x));
    
    y.idx0 = (1:nrow(x))[delta==0];
    y.idx0 = y.idx0[order(y[delta==0])];
    ds.idx0 = slice.idx(y.idx0, n.slice0, nrow(x));
    
    y.idx1 = (1:nrow(x))[delta==1];
    y.idx1 = y.idx1[order(y[delta==1])];
    ds.idx1 = slice.idx(y.idx1, n.slice1, nrow(x));
    
    ds.idx[ds.idx0 > 0] = ds.idx0[ds.idx0 > 0];
    ds.idx[ds.idx1 > 0] = ds.idx1[ds.idx1 > 0] + n.slice0;
    
    n.per.slice = as.vector(table(ds.idx));
    ave.covx = matrix(0, ncol(x), ncol(x))  # Initialisation avec une matrice de zéros
    probs = n.per.slice/nrow(x);
    for(i in 1:length(n.per.slice)){
      x.slice = as.matrix(x[(1:nrow(x))[ds.idx==i],]);
      cov.slice = cov1(x.slice, use="pair");
      if(i==1){
        ave.covx = probs[i]*cov.slice;
      }else{
        ave.covx = ave.covx + probs[i]*cov.slice;
      }
    }
    
    sigma.x = cov1(x, use="pair");
    sigma.eta = sigma.x - ave.covx;
    
    eign = eigen.decomp(sigma.eta, sigma.x);
    ordr = order(eign$values, decreasing=TRUE);
    e.val = eign$values[ordr];
    e.vec = eign$vectors[,ordr];
    
    if(is.null(colnames(x))){
      x.names = paste("X", 1:ncol(x), sep="");
    }else{
      x.names = colnames(x);
    }
    dimnames(e.vec) = list(x.names, paste("Dir", 1:ncol(e.vec), sep=""));
    
    n.s = n.slice1 + n.slice0;
    ds = list(y=y, x=x, delta=delta, eval=e.val, evec=e.vec, nslice=n.s);
    ds$nslice0 = n.slice0;
    ds$nslice1 = n.slice1;
    
    class(ds) = "ds";
    ds
  }
### version corrected for the problemes of emty matrix 
double.slice2 <- function(y, delta, x, n.slice1, n.slice0) {
  cov1 <- function(x, use = "pair") {
    if (nrow(x) == 0) {
      return(matrix(0, ncol(x), ncol(x)))  # Return a zero matrix if x is empty
    }
    n = nrow(x)
    c = cov(x, use = use) * (n - 1) / n
    c
  }
  
  slice.idx <- function(y, h, m) {
    ly = length(y)
    if (m < max(y)) {
      stop("y is not a valid index vector")
    }
    n1 = ceiling(ly / h)
    n2 = floor(ly / h)
    split = numeric(m)
    if (n1 == n2) {
      n = ly / h
      index = rep(n, h)
    } else if (n1 > n2) {
      h1 = ly %% h
      h2 = h - h1
      index = rep(c(n1, n2), c(h1, h2))
    }
    end = 0
    for (i in 1:h) {
      bgn = end + 1
      end = end + index[i]
      split[y[bgn:end]] = i
    }
    split
  }
  
  # Order data by survival time
  ordr = order(y, -delta)
  y = y[ordr]
  x = x[ordr, ]
  delta = delta[ordr]
  
  ds.idx = numeric(nrow(x))
  
  y.idx0 = (1:nrow(x))[delta == 0]
  y.idx0 = y.idx0[order(y[delta == 0])]
  ds.idx0 = slice.idx(y.idx0, n.slice0, nrow(x))
  
  y.idx1 = (1:nrow(x))[delta == 1]
  y.idx1 = y.idx1[order(y[delta == 1])]
  ds.idx1 = slice.idx(y.idx1, n.slice1, nrow(x))
  
  ds.idx[ds.idx0 > 0] = ds.idx0[ds.idx0 > 0]
  ds.idx[ds.idx1 > 0] = ds.idx1[ds.idx1 > 0] + n.slice0
  
  n.per.slice = as.vector(table(ds.idx))
  probs = n.per.slice / nrow(x)
  
  # Initialize ave.covx
  ave.covx = matrix(0, ncol(x), ncol(x))
  
  for (i in 1:length(n.per.slice)) {
    idx = (1:nrow(x))[ds.idx == i]
    if (length(idx) > 1) {
      x.slice = as.matrix(x[idx, ])
      cov.slice = cov1(x.slice, use = "pair")
    } else {
      cov.slice = matrix(0, ncol(x), ncol(x))  # Use a zero matrix if the slice is empty or has only 1 row
    }
    ave.covx = ave.covx + probs[i] * cov.slice
  }
  
  sigma.x = cov1(x, use = "pair")
  sigma.eta = sigma.x - ave.covx
  
  eign = eigen.decomp(sigma.eta, sigma.x)
  ordr = order(eign$values, decreasing = TRUE)
  e.val = eign$values[ordr]
  e.vec = eign$vectors[, ordr]
  
  if (is.null(colnames(x))) {
    x.names = paste("X", 1:ncol(x), sep = "")
  } else {
    x.names = colnames(x)
  }
  dimnames(e.vec) = list(x.names, paste("Dir", 1:ncol(e.vec), sep = ""))
  
  n.s = n.slice1 + n.slice0
  ds = list(y = y, x = x, delta = delta, eval = e.val, evec = e.vec, nslice = n.s)
  ds$nslice0 = n.slice0
  ds$nslice1 = n.slice1
  
  class(ds) = "ds"
  ds
}

##### cen.wls.sir
cen.wls.sir <-function(x, y, delta, n.slice1, n.slice0, n.slice=10, cn1=0.1, cn2=1, choose.dir=FALSE, categorical=FALSE, ndim=1 ){
  ## require package MASS and dr
  ## basic information about x and y
  n <- dim(x)[1]
  p <- dim(x)[2]
  h <- n.slice
  x = scale(x,scale=FALSE)
  
  ## slice
  if( categorical==TRUE ){
    index <- as.numeric(factor(y))
    nh <- summary(factor(y))
    h <- length(nh)
  }
  if( categorical == FALSE ){
    slice <- dr.slices.arc(y,h)
    index <- slice$slice.indicator		# Slice Index
    nh <- slice$slice.sizes				# Observations per Slice
  }
  ph <- nh/n
  
  ####################################
  ## calculate wls
  ####################################
  wls <- c()								# Weighted Leverage Score
  
  ####################################
  ## scenario 1
  ####################################
  svdx <- svd(x)
  u <- svdx$u
  d <- svdx$d
  v <- svdx$v
  
  if( choose.dir==FALSE) {
    dir = min(n, p)
  } else {
    ## selection of d   
    theta = d^2/(d[1])^2 + 1
    loglik <- penalty <- rep(0, length(d) )
    for( i in 1:length(d) ){
      if(i < length(d)) {
        loglik[i] <- sum( log(theta[(i+1):length(d)]) + 1 - theta[(i+1):length(d)] )
      } else loglik[i] = 0
      penalty[i] <- i*cn1/sqrt(n)
    }
    BIC = -loglik + penalty
    (dir <- which.min(BIC))
    print(dir)
  }
  
  ## calculate WLS
  w <- matrix(ncol = dir, nrow=length(nh))
  for(j in 1:dir){
    for(i in 1:length(nh)) w[i,j] <- sum(u[,j] * (index == i))/nh[i]
  }
  
  uut <- array(dim = c(length(nh), dir, dir))			# UUT Array
  for(i in 1:length(nh)) uut[i,,] <- (nh[i]) * ( w[i,] %*% t(w[i,]) )
  
  ## LEVERAGE SCORES
  for(j in 1:p) wls[j] <- t(v[j,1:dir]) %*% colSums(uut) %*% v[j,1:dir]
  
  
  ####################################
  ## BIC
  ####################################
  wls.sort = sort(wls, decreasing = TRUE)
  
  loglik <- penalty <- rep(0, min(n,p))
  for(k in 1:min(n,p)){
    temp_loglik = sum(wls.sort[1:k])
    (loglik[k] = -log(temp_loglik))
    penalty[k] = (log(n) + cn2*log(p))*k/max(n,p)
  }
  BIC <- loglik + penalty
  (sel.k <- which.min(BIC))
  
  select <- order(wls, decreasing = T)[1:sel.k]
  z<-x[,select]
  ds = double.slice2(y, delta, z, n.slice1, n.slice0);
  joint.edrs  = edr.n(ds, 2);
  # --- find sir direction ---
  sir<- cen.sir(y, delta, z, n.slice, joint.edrs, h, c )
  # teta.sir<-sir$evec[,1]; #beta.sir
  #outsir = do.rsir(z, y,ndim=ndim,regmethod="Ridge")
  ### use the function do.sir from the package "rdimtools" we can also use package "dr"
  
  #theta<-outsir$b
  betahat<-matrix(0,nrow=ncol(x),ncol=ndim)
  for(j in 1:ndim){
    theta<-sir$evec[,j]
    betahat[select,j]<-theta 
  }
  return( list(wls=wls, select=select,betahat=betahat,n.sel=sel.k) )
  
}

#### uses cen.wls directly 
cen.wls.sir2<- function(x, y, delta, n.slice1, n.slice0, n.slice=10, cn1=0.1, cn2=1, choose.dir=FALSE, ndim=1 ,c=0.05){
  #CENSORsir+wls
  cwlsir<-cen.wls(x,as.vector(y),delta, c , n.slice , cn1 , cn2, choose.dir) 
  select<-cwlsir$select;n.sel<-cwlsir$n.sel
  Z<-x[,select]
  ds = double.slice2(y, delta, Z, n.slice1, n.slice0);
  joint.edrs  = edr.n(ds, 2);
  beta_init <- prcomp(x, rank. = 1)$rotation[, 1]
  # Rapide et justifié théoriquement → Silverman
  h<- bw.nrd(as.vector(x%*% beta_init))
  # Plus précis si vous voulez → Sheather-Jones
  #h <- bw.SJ(as.vector(X %*% beta_init))
  censir = cen.sir(y, delta, Z, n.slice, joint.edrs, h, c )
  beta.hat<-matrix(0,nrow=ncol(x),ncol=ndim)
  eval.hat<-rep(0,ncol(x))
  eval.hat[select]<-censir$eval
  for(j in 1:ndim){
    theta<-censir$evec[,j];
    beta.hat[select,j]<-theta 
  }
  rownames(beta.hat)=colnames(x)
  return(list(beta.hat = beta.hat, eval.hat=eval.hat,select = select,n.sel=n.sel))
}


cen.wls.cox <-function(x, y, delta, n.slice1, n.slice0, nslice=10, cn1=0.1, cn2=1, choose.dir=FALSE, categorical=FALSE, ndim=1 ){
  ## require package MASS and dr
  ## basic information about x and y
  n <- dim(x)[1]
  p <- dim(x)[2]
  h <- nslice
  x = scale(x,scale=FALSE)
  
  ## slice
  if( categorical==TRUE ){
    index <- as.numeric(factor(y))
    nh <- summary(factor(y))
    h <- length(nh)
  }
  if( categorical == FALSE ){
    slice <- dr.slices.arc(y,h)
    index <- slice$slice.indicator		# Slice Index
    nh <- slice$slice.sizes				# Observations per Slice
  }
  ph <- nh/n
  
  ####################################
  ## calculate wls
  ####################################
  wls <- c()								# Weighted Leverage Score
  
  ####################################
  ## scenario 1
  ####################################
  svdx <- svd(x)
  u <- svdx$u
  d <- svdx$d
  v <- svdx$v
  
  if( choose.dir==FALSE) {
    dir = min(n, p)
  } else {
    ## selection of d   
    theta = d^2/(d[1])^2 + 1
    loglik <- penalty <- rep(0, length(d) )
    for( i in 1:length(d) ){
      if(i < length(d)) {
        loglik[i] <- sum( log(theta[(i+1):length(d)]) + 1 - theta[(i+1):length(d)] )
      } else loglik[i] = 0
      penalty[i] <- i*cn1/sqrt(n)
    }
    BIC = -loglik + penalty
    (dir <- which.min(BIC))
    print(dir)
  }
  
  ## calculate WLS
  w <- matrix(ncol = dir, nrow=length(nh))
  for(j in 1:dir){
    for(i in 1:length(nh)) w[i,j] <- sum(u[,j] * (index == i))/nh[i]
  }
  
  uut <- array(dim = c(length(nh), dir, dir))			# UUT Array
  for(i in 1:length(nh)) uut[i,,] <- (nh[i]) * ( w[i,] %*% t(w[i,]) )
  
  ## LEVERAGE SCORES
  for(j in 1:p) wls[j] <- t(v[j,1:dir]) %*% colSums(uut) %*% v[j,1:dir]
  
  
  ####################################
  ## BIC
  ####################################
  wls.sort = sort(wls, decreasing = TRUE)
  
  loglik <- penalty <- rep(0, min(n,p))
  for(k in 1:min(n,p)){
    temp_loglik = sum(wls.sort[1:k])
    (loglik[k] = -log(temp_loglik))
    penalty[k] = (log(n) + cn2*log(p))*k/max(n,p)
  }
  BIC <- loglik + penalty
  (sel.k <- which.min(BIC))
  
  select <- order(wls, decreasing = T)[1:sel.k]
  z<-x[,select]
  ds = double.slice(y, delta, z, n.slice1, n.slice0);
  joint.edrs  = edr.n(ds, 2);
  # --- find sir direction ---
  cox<-coxph(Surv(y)~z)
  theta<-cox$coef
  
  # sir<- cen.sir(y, delta, z, n.slice, joint.edrs, h, c )
  # teta.sir<-sir$evec[,1]; #beta.sir
  #outsir = do.rsir(z, y,ndim=ndim,regmethod="Ridge")
  ### use the function do.sir from the package "rdimtools" we can also use package "dr"
  
  #theta<-outsir$b
  betahat<-matrix(0,nrow=ncol(x),ncol=ndim)
  for(j in 1:ndim){
    betahat[select,j]<-theta 
  }
  return( list(wls=wls, select=select,betahat=betahat) )
  
}




compute_TPR_FPR <- function(estimated_active, true_active, total_features) {
  # Convert to unique sets to avoid duplicates
  estimated_active <- unique(estimated_active)
  true_active <- unique(true_active)
  
  # Calculate True Positives (TP), False Positives (FP), and False Negatives (FN)
  TP <- length(intersect(estimated_active, true_active))  # True Positives
  FP <- length(setdiff(estimated_active, true_active))    # False Positives
  FN <- length(setdiff(true_active, estimated_active))    # False Negatives
  
  # Calculate True Negatives (TN)
  TN <- total_features - (TP + FP + FN)                   # True Negatives
  
  # Check if TP + FP + FN exceeds total_features (this should not happen)
  if (TP + FP + FN > total_features) {
    stop("The sum of TP, FP, and FN exceeds the total number of features (p). Check your data!")
  }
  
  # Compute TPR and FPR with checks to avoid division by zero
  TPR <- ifelse((TP + FN) > 0, TP / (TP + FN), 0)
  FPR <- ifelse((FP + TN) > 0, FP / (FP + TN), 0)
  FDR <- ifelse((FP + TP) > 0, FP / (FP + TP), 0)  # False Discovery Rate
  
  # Return the metrics as a list
  return(list(TPR = TPR, FPR = FPR, FDR = FDR, TP = TP, FP = FP, FN = FN, TN = TN))
}

#Model-Free Feature screening Based on Concordance Index Statistic--------
#'
#' A model-free and data-adaptive feature screening method for
#' ultrahigh-dimensional data and even survival data. The proposed method is based
#' on the concordance index which measures concordance between random vectors even
#' if one of the vectors is a survival object Surv. This rank correlation based
#' method does not require specifying a regression model, and applies robustly to data
#' in the presence of censoring and heavy tails. It enjoys both sure screening and rank
#' consistency properties under weak assumptions.
#'
#' @param X The design matrix of dimensions n * p. Each row is an observation vector.
#' @param Y The response vector of dimension n * 1. For survival models,
#' Y should be an object of class Surv, as provided by the function
#' Surv() in the package survival.
#' @param nsis Number of predictors recruited by CSIS. The default is n/log(n).
#'
#' @return the labels of first nsis largest active set of all predictors
#'
#' @importFrom survival concordancefit
#' @importFrom survival Surv
#' @import foreach
#' @import parallel
#' @import doParallel
#'
#' @export
#' @author Xuewei Cheng \email{xwcheng@hunnu.edu.cn}
#' @examples
#'
#' Cheng X, Li G, Wang H. The concordance filter: an adaptive model-free feature screening procedure[J]. Computational Statistics, 2023: 1-24.
CSIS <-function(X, Y, nsis = (dim(X)[1]) / log(dim(X)[1])) {
  if (dim(X)[1] != length(Y)) {
    stop("X and Y should have same number of rows!")
  }
  if (missing(X) | missing(Y)) {
    stop("The data is missing!")
  }
  if (TRUE %in% (is.na(X) | is.na(Y) | is.na(nsis))) {
    stop("The input vector or matrix cannot have NA!")
  }
  n <- dim(X)[1] ## sample size
  p <- dim(X)[2] ## dimension
  B <- vector(mode = "numeric", length = p)
  Cindex <- vector(mode = "numeric", length = p)
  if (n * p <= 2000000) {
    for (i in 1:p) {
      Cindex[i] <- concordancefit(Y, X[, i])$concordance
    }
  } else {
    # Real physical cores in the computer
    cores <- detectCores(logical = FALSE)
    cl <- makeCluster(cores)
    registerDoParallel(cl, cores = cores)
    j <- NULL
    Cindex <- foreach::foreach(
      j = 1:p, .combine = "c",
      .packages = c("survival")
    ) %dopar%
      concordancefit(Y, X[, j])$concordance
    stopImplicitCluster()
    stopCluster(cl)
  }
  num <- which(Cindex < 0.5)
  B[num] <- 1 - Cindex[num]
  B[-num] <- Cindex[-num]
  A <- order(B, decreasing = TRUE)
  return(A[1:nsis])
}

#' Fan, J. and J. Lv (2008). Sure independence screening for ultrahigh dimensional feature space. Journal of the Royal Statistical Society: Series B (Statistical Methodology) 70(5),849–911.
#'
#' Li, R., W. Zhong, and L. Zhu (2012). Feature screening via distance correlation learning. Journal of the American Statistical Association 107(499), 1129–1139.
DCSIS=function(X,Y,nsis=(dim(X)[1])/log(dim(X)[1])){
  if (dim(X)[1]!=length(Y)) {
    stop("X and Y should have same number of rows!")
  }
  if (missing(X)|missing(Y)) {
    stop("The data is missing!")
  }
  if (TRUE%in%(is.na(X)|is.na(Y)|is.na(nsis))) {
    stop("The input vector or matrix cannot have NA!")
  }
  if (inherits(Y,"Surv")) {
    stop("DCSIS can not implemented with object  of Surv")
  }
  n=dim(X)[1]; ##sample size
  p=dim(X)[2]; ##dimension
  B=matrix(1,n,1);
  C=matrix(1,1,p);
  sxy1=matrix(0,n,p);
  sxy2=matrix(0,n,p);
  sxy3=matrix(0,n,1);
  sxx1=matrix(0,n,p);
  syy1=matrix(0,n,1);
  for (i in 1:n){
    XX1=abs(X-B%*%X[i,]);
    YY1=sqrt(apply((Y-B%*%Y[i])^2,1,sum))
    sxy1[i,]=apply(XX1*(YY1%*%C),2,mean);
    sxy2[i,]=apply(XX1,2,mean);
    sxy3[i,]=mean(YY1);
    XX2=XX1^2;
    sxx1[i,]=apply(XX2,2,mean);
    YY2=YY1^2;
    syy1[i,]=mean(YY2);
  }
  SXY1=apply(sxy1,2,mean);
  SXY2=apply(sxy2,2,mean)*apply(sxy3,2,mean);
  SXY3=apply(sxy2*(sxy3%*%C),2,mean);
  SXX1=apply(sxx1,2,mean);
  SXX2=apply(sxy2,2,mean)^2;
  SXX3=apply(sxy2^2,2,mean);
  SYY1=apply(syy1,2,mean);
  SYY2=apply(sxy3,2,mean)^2;
  SYY3=apply(sxy3^2,2,mean);
  dcovXY=sqrt(SXY1+SXY2-2*SXY3);
  dvarXX=sqrt(SXX1+SXX2-2*SXX3);
  dvarYY=sqrt(SYY1+SYY2-2*SYY3);
  dcorrXY=dcovXY/sqrt(dvarXX*dvarYY);
  A=order(dcorrXY,decreasing=TRUE)
  return (A[1:nsis])
}



  #@references
#'
#' Zhu, L.-P., L. Li, R. Li, and L.-X. Zhu (2011). Model-free feature screening for ultrahigh-dimensional data. Journal of the American Statistical Association 106(496), 1464–1475.
#'
SIRS=function(X,Y,nsis=(dim(X)[1])/log(dim(X)[1])){
  if (dim(X)[1]!=length(Y)) {
    stop("X and Y should have same number of rows!")
  }
  if (missing(X)|missing(Y)) {
    stop("The data is missing!")
  }
  if (TRUE%in%(is.na(X)|is.na(Y)|is.na(nsis))) {
    stop("The input vector or matrix cannot have NA!")
  }
  if (inherits(Y,"Surv")) {
    stop("SIRS can not implemented with object  of Surv")
  }
  n=dim(X)[1]; ##sample size
  posit=order(Y,decreasing=FALSE)
  Y=Y[posit];
  X=X[posit,];
  xx=(X-as.matrix(rep(1,n))%*%apply(X,2,mean))/(as.matrix(rep(1,n))%*%apply(X,2,sd))
  B=matrix(1,n,n)
  B[!upper.tri(B,diag = TRUE)]=0
  USIRS=apply((t(xx)%*%B/n)^2,1,mean);
  A=order(USIRS,decreasing=TRUE)
  B=A[1:nsis]
  return (B)
}


#' @references
#'
#' Fan, J. and J. Lv (2008). Sure independence screening for ultrahigh dimensional feature space. Journal of the Royal Statistical Society: Series B (Statistical Methodology) 70(5),849–911.
SIS=function(X,Y,nsis=(dim(X)[1])/log(dim(X)[1])){
  if (dim(X)[1]!=length(Y)) {
    stop("X and Y should have same number of rows!")
  }
  if (missing(X)|missing(Y)) {
    stop("The data is missing!")
  }
  if (TRUE%in%(is.na(X)|is.na(Y)|is.na(nsis))) {
    stop("The input vector or matrix cannot have NA!")
  }
  if (inherits(Y,"Surv")) {
    stop("SIS can not implemented with object  of Surv")
  }
  A=order(abs(t(X)%*%Y),decreasing=TRUE)
  return(A[1:nsis])
}

compute_censoring_rate <- function(c, Y0, errc) {
  C <- exp(errc + c)
  y <- pmin(Y0, C)
  delta <- as.numeric(Y0 <= C)
  censoring_rate <- 1 - mean(delta)
  return(censoring_rate)
}

# Function to find c for a given censoring rate
find_c_for_censoring_rate <- function(target_rate, Y0, errc) {
  objective_function <- function(c) {
    current_rate <- compute_censoring_rate(c, Y0, errc)
    return(current_rate - target_rate)
  }
  
  # Use root-finding to solve for c
  result <- uniroot(objective_function, interval = c(-10, 10), extendInt = "yes")
  return(result$root)
}

