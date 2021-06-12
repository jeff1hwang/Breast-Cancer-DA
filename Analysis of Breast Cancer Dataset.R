#########################
## Jeff Hwang
# Copyright (c) 2021, Jeff Hwang
# All rights reserved.
# 
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree. 
#########################

#########################
### (a) Cluster Analysis using EM algorithm
##########################

## Read the data
data <- read.csv("Data_Breast_Cancer.csv", header = T)
head(data)
dim(data) # check the dimension
datanew <- data[, -31] # exclude column 31
Y <- data.matrix(datanew)
Y[c(1:3),] # double check

## Scale  and center the data
X <- scale(Y, scale = T)

### Likelihood ignoring the constant term 
### and assuming each row of X, is a 2-vector x 
## from a multivariate normal with big Sigma 
### equal to sigma^2I  
f <- function(x,mu,sigma2){
  exp((-1/(2*sigma2))*sum((x-mu)^2))
}

EM_algorithm <- function(X, sigma2=1) {
  n <- nrow(X)
  
  # Initialization 
  ## Setting initial alpha and mu values
  alpha <- matrix(rep(1/2,2),ncol=2) # probability of being located in each cluster is the same
  ### Setting 2 distinct observations as the starting cluster centers
  for (i in 2:nrow(X)){
    if(sum(X[1,]-X[i,])!=0){
      mu <- matrix(c(X[1,],X[i,]),ncol = 2)
      break
    }else{
      if (i==nrow(X)) stop("All observations are the same")
    }
  }
  ### Setting different indicators for testing convergence
  indicator <- 1:nrow(X)
  pastIndicator <- nrow(X):1
  t <- 1 ## for storing values in the alpha matrix
  while(sum(pastIndicator!=indicator)!=0){
    pastIndicator <- indicator
    
    ## E-step - finding the weights given alpha_k and mu_k (1pt)
    w1 <- apply(X,1,function(x){
      alpha[t,1]*f(x,as.numeric(mu[,1]),sigma2)
    })
    w2 <- apply(X,1,function(x){
      alpha[t,2]*f(x,as.numeric(mu[,2]),sigma2)
    })
    W <- cbind(w1,w2)
    W <- W/rowSums(W)
    
    ## Assigning observations to clusters based on weights
    indicator <- max.col(W)
    
    ## M step - Updating alphas (1pt)
    nk <- colSums(W)
    alpha <- rbind(alpha,nk/n)
    t <- t+1
    
    ## M step - updating the mus (1pt)
    temp1 <- (W[,1]*X)/nk[1]
    temp2 <- (W[,2]*X)/nk[2]
    
    mu[,1] <- colSums(temp1)
    mu[,2] <- colSums(temp2)
  }
  cat("The mu for the first cluster is (",as.numeric(mu[1,1]),",",as.numeric(mu[2,1]),")\n",sep = "")
  cat("The mu for the second cluster is (",as.numeric(mu[1,2]),",",as.numeric(mu[2,2]),")\n\n",sep = "")
  
  cat("The probability of being in the first cluster is ",alpha[t,1],"\n",sep = "")
  cat("The probability of being in the second cluster is ",alpha[t,2],"\n\n",sep = "")
  
  cat("The final allocation vector is given by\n",indicator,"\n\n",sep = "")
  
  cat("This algorithm took ",t," iterations\n",sep = "")
  
  indicator
}
X_EM <- EM_algorithm(X)
X_EM
cluster=X_EM


###############################
###### Summary for the result of clustering
###############################
## Z vector
Z <- cluster

#### Merging the cluster allocation with the original data set
data.cluster <- cbind(Z, data)

## group of cluster 1
data.cluster.1 <- data.cluster[Z==1,]
nrow(data.cluster.1) # number of observations in cluster 1

## group of cluster 2
data.cluster.2 <- data.cluster[Z==2,]
nrow(data.cluster.2) # number of observations in cluster 2

######## Summary for Cluster 1
mean.1 <- colMeans(data.cluster.1)

median.1 <- apply(data.cluster.1, 2, median)

sd.1 <- apply(data.cluster.1, 2, sd)

## Combine the above data and print it out
summary.1 <- rbind(mean.1, median.1, sd.1)
summary.1

######### Summary for cluster 2
mean.2 <- colMeans(data.cluster.2)

median.2 <- apply(data.cluster.2, 2, median)

sd.2 <- apply(data.cluster.2, 2, sd)

### combine the above data and print it out
summary.2 <- rbind(mean.2, median.2, sd.2)
summary.2

# Scatter plot that has mean area and mean concavity
plot(data[, "mean.area"], data[, "mean.concavity"],
     col=c("red", "blue")[unclass(Z)],
     main = "(A) Mean Area vs. Mean Concavity", xlab = "Mean Area",
     ylab = "Mean Concavity")
legend("bottomright", c("Cluster 1", "Cluster 2"), pch = c("R", "B"), col = c("red", "blue"))

# Scatter plot that has mean area and mean compactness
plot(data[, "mean.area"], data[, "mean.compactness"],
     col=c("red", "blue")[unclass(Z)],
     main = "(B) Mean Area vs. Mean Compactness", xlab = "Mean Area",
     ylab = "Mean Compactness")
legend("bottomright", c("Cluster 1", "Cluster 2"), pch = c("R", "B"), col = c("red", "blue"))

# Scatter plot that has mean symmetric and mean compactness
plot(data[, "mean.area"], data[, "mean.concave.points"],
     col=c("red", "blue")[unclass(Z)],
     main = "(C) Mean Area vs. Mean Concave Points", xlab = "Mean Area",
     ylab = "Mean Concave points")
legend("bottomright", c("Cluster 1", "Cluster 2"), pch = c("R", "B"), col = c("red", "blue"))

###############################
#### (b) PC analysis
###############################

## different units, so we use the centered and scaled data X from above
X[1:5, 1:5] # double check
## check that correlation is present between variables
correlation <- cor(X)
# Because too many variables, we only print the first few variables to check
print(correlation[1:5, 1:5])

### Calculate the variance covariance matrix of the Xcs
Sx <- var(X)
print(Sx[1:5, 1:5])

### For PC, you need the eigenvalues and eigenvectors 
## of the variance-covariance matrix

EP <- eigen(Sx)
V <- EP$vectors
print(V[1:5, 1:5]) # view the first few data because it't too big

## Check eigenvectors are orthonormal
temp1 <- V%*%t(V)
temp2 <- t(V)%*%V
print(temp1[1:5, 1:5]) # view the first few variables
print(temp2[1:5, 1:5]) # view the first few variables

Lambda <- EP$values
Lambda

# see what percent of variance 
# is contained in each PC
## Since eigenvalues are variances of the PC variables, 
## and sum of PC variance- sum of data variance, the 
## eigenvalues tell us how much variance we are preserving 
## and the following calculation help us do dimension reduction.

cumsum(Lambda) / sum(Lambda) # if we want to preserve close to 90% variance,
## choose first 7 variables is what this calculation is saying

# Calculate the matrix of PC values for each observation in the 
# original data set. 
PC <- X%*%V
print(PC[1:5, 1:5])
## Check that the PC variables are uncorrelated
PC_cor <- cor(PC)
print(PC_cor[1:5, 1:5])
PC
# Plot the PC
plot(cumsum(Lambda) / sum(Lambda), xaxt = "n", col = "red", ylab = "cumsum of the lambda")
axis(side = 1, at = seq(1, 30, by = 1))
lines(cumsum(Lambda) / sum(Lambda))


##############################
### Find MLE and Confidence Interval with Newton
##############################

### Finding initial values from contour plot

## Estimation of variable distribution 
hist(datanew[,4], probability = T, ylim = c(0,0.0025))
x <- 0:max(datanew$mean.area)
# Simulate the Gamma Distribution
curve(dgamma(x, shape = 4, rate = 0.005), add = TRUE, col = "red")

## start calculation
y <- datanew$mean.area
n <- length(datanew$mean.area)

# The log-likelihood function
fn <- function(p, y) {
  n * log(p[1]) + (p[1]-1) * sum(log(y)) - n * p[1] * log(p[2]) - ((p[2])^(-p[1])) * sum(y^p[1])
}

x1 <- seq(0,20, by=0.01) # values of alpha

x2 <- seq(0,0.05,by=0.01) #values of lambda

f <- matrix(0,nrow=length(x1),ncol=length(x2))

for(i in 1:length(x1)){
  for(j in 1:length(x2)) {
    f[i,j]=(x1[i]-1)*sum(log(y)) - sum(y)*x2[j] - (length(y))*(lgamma(x1[i])) +(length(y))*x1[i]*log(x2[j])
  }}
# contour plot takes time 
contour(x1,x2,f,nlevels=60,xlab="alpha",ylab="lambda")

##### Now numberical optimization
xt <- c(100,0)
eps <- 0.00000000001
xtp1 <- c(5,0.01)
xHist <- matrix(xtp1,2,1)
fHist=c() # history of log likelihood 
xHist  # history of parameter vectors
# objective function=log likelihood 
f <- (xtp1[1]-1)*sum(log(y)) - sum(y)*xtp1[2] - n*(lgamma(xtp1[1])) +n*xtp1[1]*log(xtp1[2])
while(sum((xtp1-xt)^2)>eps){
  xt=xtp1  # the last x^(t+1) will be used to compute Newton's
  xt
  gradient=as.vector(c(n*log(xt[2])+sum(log(y))-n*digamma(xt[1]) ,
                       n*xt[1]/xt[2]-sum(y) ))
  gradient
  hessian=matrix(c(-n*trigamma(xt[1]),n/xt[2], n/xt[2], -n*xt[1]/xt[2]^2 ),ncol=2,nrow=2)
  hessian
  ###
  # compute xtp1 solve(hessian*gradient=hessian^{-1}*gradient)
  ###
  xtp1=xt-solve(hessian)%*%gradient  # Newton iteration
  xtp1
  ###
  #save history
  xHist=matrix(c(xHist,xtp1),2)
  xHist
  f=(xtp1[1]-1)*sum(log(y)) - sum(y)*xtp1[2] - n*(lgamma(xtp1[1])) +n*xtp1[1]*log(xtp1[2])
  fHist=c(fHist,f)
}
xHist
fHist

### Compute, after all this, the hessian
### at the stationary point
gradient
hessian

### Now find the standard errors of the estimators
se_newton <- sqrt(diag(solve(-hessian)))
diag(solve(-hessian))
### Calculate the Confidence Interval
## define a function for calculating Confidence Interval
confidenceInterval <- function(level=0.95, stdErr, mle) {
  CI <- mle + c(-1,1)*qnorm((1+level)/2)*stdErr
  return(CI)
}

## The CI for alpha
confidenceInterval(stdErr = se_newton, mle = xHist[1, 11])

## The CI for lambda
confidenceInterval(stdErr = se_newton, mle = xHist[2, 11])

### Finally we plot again to double check
hist(datanew[,4], probability = T, ylim = c(0,0.0025))
x <- 0:max(datanew$mean.area)
curve(dgamma(x, shape = 4.282259349, rate = 0.006538908), add = TRUE, col = "red")




# Copyright (c) <year>, <copyright holder>
#   All rights reserved.
# 
# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree. 
# 










