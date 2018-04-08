rm(list=ls())
#graphics.off()
library(lattice)
library(ggplot2)
library(car) 
#install.packages("olsrr")
library(olsrr)
#=======================================================================================#
# load data
load(file="5factors_2008_2018.RData")
#mydata2 <- read.table("F-3factors_2017_2018.csv",header=TRUE,sep=",")
dim(mydata2)

# Extract Fama-French Factors and Fund Returns
mydate <- mydata2$Date
range(mydate)
rmrf <- mydata2$Mkt.RF
smb <- mydata2$SMB
hml <- mydata2$HML
rf <- mydata2$RF
rmw <- mydata2$RMW
cma <- mydata2$CMA
# Calculate Excess Returns for Target fund
lo_30.xcess <- mydata2$Lo.30 - rf
#==================================================
# Run Fama-French Regression
lo_30.ffregression <- lm(lo_30.xcess ~ rmrf + smb + hml + rmw + cma)
# Print summary of regression results
lo_30.sum <- summary(lo_30.ffregression)
lo_30.sum
#===========================================
## Construct X data frame
# X data: rmrf, smb, hml
X <- data.matrix(cbind(rep(1, dim(mydata2)[1]), rmrf, smb, hml, rmw, cma))
is.matrix(X)
dim(X)
# # obervations
M <- dim(X)[1]
# True beta
beta <- lo_30.sum$coefficients[,1]
is.vector(beta)
class(beta)
length(beta)
#
mu <- 0
sig2e <- var(lo_30.xcess)
# True sigma: std of distribution of error term
sde <- sqrt(sig2e)
sde
#=========================
### Boxplots
#library(lattice)
bwplot(~lo_30.xcess, xlab = "Portfolio excess return") # Both side outliers
bwplot(~rmrf, xlab = expression(K[m] - R[f])) 
bwplot(~smb, xlab = "SMB") 
bwplot(~hml, xlab = "HML") 
bwplot(~rmw, xlab = "RMW") 
bwplot(~cma, xlab = "CMA") 
# density plot
densityplot(lo_30.xcess, xlab = "Portfolio excess return")

#library(car)     # red line affected by outlier, green ignoring the outlier (robust)
scatterplot(x = lo_30.xcess, y = rmrf)       #(dependent)
scatterplot(x = lo_30.xcess, y = smb)       #(dependent)
scatterplot(x = lo_30.xcess, y = hml)       #(dependent)
scatterplot(x = lo_30.xcess, y = rmw)       #(dependent)
scatterplot(x = lo_30.xcess, y = cma)       #(dependent)

#library(ggplot2)
# scatterplot using ggplot
ggplot(mapping = aes(x = rmrf, y = lo_30.xcess)) + geom_point(colour = 'skyblue') + geom_smooth(method = 'lm')
ggplot(mapping = aes(x = smb, y = lo_30.xcess)) + geom_point(colour = 'skyblue') + geom_smooth(method = 'lm')
ggplot(mapping = aes(x = hml, y = lo_30.xcess)) + geom_point(colour = 'skyblue') + geom_smooth(method = 'lm')
ggplot(mapping = aes(x = rmw, y = lo_30.xcess)) + geom_point(colour = 'skyblue') + geom_smooth(method = 'lm')
ggplot(mapping = aes(x = cma, y = lo_30.xcess)) + geom_point(colour = 'skyblue') + geom_smooth(method = 'lm')
#=======================================================================================#
# One example
#=======================================================================================#
set.seed(12345)
Y <- X %*% beta + rnorm(M, mean=0, sd=sde)
class(Y)
is.vector(Y)
is.matrix(Y)
length(Y)

# Full Model
ff<-lm(Y~X[,-1])
sum <- summary(ff)
sum
rsq.adj <- sum$adj.r.squared
rsq.adj
aic <- AIC(ff)
aic
bic <- AIC(ff, k = log(M))
bic

#
ff2<-lm(Y~X[,c(-1,-2)])
sum2 <- summary(ff2)
sum2
rsq.adj2 <- sum2$adj.r.squared
rsq.adj2
aic2 <- AIC(ff)
aic2
bic2 <- AIC(ff, k = log(M))
bic2
#install.packages("olsrr")
#library(olsrr)
mlcp2 <- ols_mallows_cp(ff2, ff)
mlcp2

#==============================
# Confidence interval
cfi<-confint(ff)
cfi
class(cfi)
dim(cfi)
# If true beta is in confidence interval of estimated beta
flag.beta <- logical(length(beta))
for (k in 1:length(beta)){
flag.beta[k] <- cfi[k,1]<=beta[k]&&cfi[k,2]>=beta[k]
}
flag.beta
class(flag.beta)
length(flag.beta)

#==============================
# estimated error variance
err.var <- deviance(ff)/df.residual(ff)
err.var

#==============================
# \hat{\beta}
ff$coefficients

#==============================
# T-Test: Null Hypothesis: betahat = beta
# Significant level of t-test
alp <- 0.05
# T statistics = (betahat-beta)/s.e.
t.stat <- (sum$coefficients[,1]-beta)/sum$coefficients[,2]
# t_{1-\alpha/2, n-k}
thr <- qt(1-alp/2, df = df.residual(ff))
# whether reject null hypothesis: False - not reject; True - reject (Type I)
# We should not reject the null hypothesis of beta = betahat
# But if we reject it (<alp), we are making Type I error, right?
flag.type1 <- t.stat>thr | t.stat < - thr
flag.type1
class(flag.type1)
length(flag.type1)

#=======================================================================================#
# Running Monte Carlo Simulation
#=======================================================================================#
# Vector of N (Replication number)
#vecN <- c(10,50,100)
# vecN <- c(10,50,100,500,1e3,5e3,1e4,5e4,1e5,5e5)
# level of t-test
alpha <- 0.05
N <- 1e5
#
temp <- 1:length(beta)
nbv2 <- temp[-2]
nbv3 <- temp[-3]
nbv4 <- temp[-4]
nbv5 <- temp[-5]
nbv6 <- temp[-6]
nbv7 <- temp[1:4]
nbv8 <- temp[1:2]

mat.beta <- matrix(0,nrow = N, ncol = length(beta))
mat.beta2 <- matrix(0,nrow = N, ncol = length(nbv2))
mat.beta3 <- matrix(0,nrow = N, ncol = length(nbv3))
mat.beta4 <- matrix(0,nrow = N, ncol = length(nbv4))
mat.beta5 <- matrix(0,nrow = N, ncol = length(nbv5))
mat.beta6 <- matrix(0,nrow = N, ncol = length(nbv6))
mat.beta7 <- matrix(0,nrow = N, ncol = length(nbv7))
mat.beta8 <- matrix(0,nrow = N, ncol = length(nbv8))

vec.flag.beta <- logical(N*length(beta))

vec.err.var <- numeric(N)
vec.err.var2 <- numeric(N)
vec.err.var3 <- numeric(N)
vec.err.var4 <- numeric(N)
vec.err.var5 <- numeric(N)
vec.err.var6 <- numeric(N)
vec.err.var7 <- numeric(N)
vec.err.var8 <- numeric(N)

rsq.adj <- numeric(N)
rsq.adj2 <- numeric(N)
rsq.adj3 <- numeric(N)
rsq.adj4 <- numeric(N)
rsq.adj5 <- numeric(N)
rsq.adj6 <- numeric(N)
rsq.adj7 <- numeric(N)
rsq.adj8 <- numeric(N)

aic <- numeric(N)
aic2 <- numeric(N)
aic3 <- numeric(N)
aic4 <- numeric(N)
aic5 <- numeric(N)
aic6 <- numeric(N)
aic7 <- numeric(N)
aic8 <- numeric(N)

bic <- numeric(N)
bic2 <- numeric(N)
bic3 <- numeric(N)
bic4 <- numeric(N)
bic5 <- numeric(N)
bic6 <- numeric(N)
bic7 <- numeric(N)
bic8 <- numeric(N)

mlcp2 <- numeric(N)
mlcp3 <- numeric(N)
mlcp4 <- numeric(N)
mlcp5 <- numeric(N)
mlcp6 <- numeric(N)
mlcp7 <- numeric(N)
mlcp8 <- numeric(N)

vec.flag.type1 <- logical(N*length(beta))
vec.flag.type1.2 <- logical(N*(length(nbv2)))
vec.flag.type1.3 <- logical(N*(length(nbv3)))
vec.flag.type1.4 <- logical(N*(length(nbv4)))
vec.flag.type1.5 <- logical(N*(length(nbv5)))
vec.flag.type1.6 <- logical(N*(length(nbv6)))
vec.flag.type1.7 <- logical(N*(length(nbv7)))
vec.flag.type1.8 <- logical(N*(length(nbv8)))

mat.t.stat <- matrix(0,nrow = N, ncol = length(beta))
mat.t.stat2 <- matrix(0,nrow = N, ncol = length(nbv2))
mat.t.stat3 <- matrix(0,nrow = N, ncol = length(nbv3))
mat.t.stat4 <- matrix(0,nrow = N, ncol = length(nbv4))
mat.t.stat5 <- matrix(0,nrow = N, ncol = length(nbv5))
mat.t.stat6 <- matrix(0,nrow = N, ncol = length(nbv6))
mat.t.stat7 <- matrix(0,nrow = N, ncol = length(nbv7))
mat.t.stat8 <- matrix(0,nrow = N, ncol = length(nbv8))
# 
for (i in 1:N) {
  set.seed(27+6*i)
  Y <- X %*% beta + rnorm(M, mean=0, sd=sde)
  # Full Linear Regression
  ff<-lm(Y~X[,-1])
  # Omitt Rm-Rf variable
  ff2<-lm(Y~X[,c(-1,-2)])
  # Omitt SMB variable
  ff3<-lm(Y~X[,c(-1,-3)])
  # Omitt HML variable
  ff4<-lm(Y~X[,c(-1,-4)])
  # Omitt RMW variable
  ff5<-lm(Y~X[,c(-1,-5)])
  # Omitt CMA variable
  ff6<-lm(Y~X[,c(-1,-6)])
  # Omitt 2 variables, RMW&CMA variable, Fama French 3 factor model
  ff7<-lm(Y~X[,2:4])
  # Omitting 4 variables, CAPM model
  ff8<-lm(Y~X[,2])
  
  
  # Confidence interval
  cfi<-confint(ff)
  # If true beta is in confidence interval of estimated beta
  for (k in 1:length(beta)){
    vec.flag.beta[(i-1)*length(beta)+k] <- cfi[k,1]<=beta[k]&&cfi[k,2]>=beta[k]
  }
  # estimated error variance
  vec.err.var[i] <- deviance(ff)/df.residual(ff)
  vec.err.var2[i] <- deviance(ff2)/df.residual(ff2)
  vec.err.var3[i] <- deviance(ff3)/df.residual(ff3)
  vec.err.var4[i] <- deviance(ff4)/df.residual(ff4)
  vec.err.var5[i] <- deviance(ff5)/df.residual(ff5)
  vec.err.var6[i] <- deviance(ff6)/df.residual(ff6)
  vec.err.var7[i] <- deviance(ff7)/df.residual(ff7)
  vec.err.var8[i] <- deviance(ff8)/df.residual(ff8)
  
  # \hat{\beta}
  mat.beta[i,] <- as.matrix(ff$coefficients)
  mat.beta2[i,] <- as.matrix(ff2$coefficients)
  mat.beta3[i,] <- as.matrix(ff3$coefficients)
  mat.beta4[i,] <- as.matrix(ff4$coefficients)
  mat.beta5[i,] <- as.matrix(ff5$coefficients)
  mat.beta6[i,] <- as.matrix(ff6$coefficients)
  mat.beta7[i,] <- as.matrix(ff7$coefficients)
  mat.beta8[i,] <- as.matrix(ff8$coefficients)
  
  # T statistics = (betahat-beta)/s.e.
  sum <- summary(ff)
  rsq.adj[i] <- sum$adj.r.squared
  aic[i] <- AIC(ff)
  bic[i] <- AIC(ff, k = log(M))
  mat.t.stat[i,] <- (sum$coefficients[,1]-beta)/sum$coefficients[,2]
  
  sum2 <- summary(ff2)
  rsq.adj2[i] <- sum2$adj.r.squared
  aic2[i] <- AIC(ff2)
  bic2[i] <- AIC(ff2, k = log(M))
  mlcp2[i] <- ols_mallows_cp(ff2, ff)
  mat.t.stat2[i,] <- (sum2$coefficients[,1]-beta[nbv2])/sum2$coefficients[,2]
  
  sum3 <- summary(ff3)
  rsq.adj3[i] <- sum3$adj.r.squared
  aic3[i] <- AIC(ff3)
  bic3[i] <- AIC(ff3, k = log(M))
  mlcp3[i] <- ols_mallows_cp(ff3, ff)
  mat.t.stat3[i,] <- (sum3$coefficients[,1]-beta[nbv3])/sum3$coefficients[,2]
  
  sum4 <- summary(ff4)
  rsq.adj4[i] <- sum4$adj.r.squared
  aic4[i] <- AIC(ff4)
  bic4[i] <- AIC(ff4, k = log(M))
  mlcp4[i] <- ols_mallows_cp(ff4, ff)
  mat.t.stat4[i,] <- (sum4$coefficients[,1]-beta[nbv4])/sum4$coefficients[,2]
  
  sum5 <- summary(ff5)
  rsq.adj5[i] <- sum5$adj.r.squared
  aic5[i] <- AIC(ff5)
  bic5[i] <- AIC(ff5, k = log(M))
  mlcp5[i] <- ols_mallows_cp(ff5, ff)
  mat.t.stat5[i,] <- (sum5$coefficients[,1]-beta[nbv5])/sum5$coefficients[,2]
  
  sum6 <- summary(ff6)
  rsq.adj6[i] <- sum2$adj.r.squared
  aic6[i] <- AIC(ff6)
  bic6[i] <- AIC(ff6, k = log(M))
  mlcp6[i] <- ols_mallows_cp(ff6, ff)
  mat.t.stat6[i,] <- (sum6$coefficients[,1]-beta[nbv6])/sum6$coefficients[,2]
  
  sum7 <- summary(ff7)
  rsq.adj7[i] <- sum7$adj.r.squared
  aic7[i] <- AIC(ff7)
  bic7[i] <- AIC(ff7, k = log(M))
  mlcp7[i] <- ols_mallows_cp(ff7, ff)
  mat.t.stat7[i,] <- (sum7$coefficients[,1]-beta[nbv7])/sum7$coefficients[,2]
  
  sum8 <- summary(ff8)
  rsq.adj8[i] <- sum8$adj.r.squared
  aic8[i] <- AIC(ff8)
  bic8[i] <- AIC(ff8, k = log(M))
  mlcp8[i] <- ols_mallows_cp(ff8, ff)
  mat.t.stat8[i,] <- (sum8$coefficients[,1]-beta[nbv8])/sum8$coefficients[,2]
  
  # t_{1-\alpha/2, n-k}
  threshold <- qt(1-alpha/2, df = df.residual(ff))
  # whether reject null hypothesis: False - not reject; True - reject (Type I)
  # We should not reject the null hypothesis of beta = betahat
  # But if we reject it (<alp), we are making Type I error, right?
  for (k in 1:length(beta)){
    vec.flag.type1[(i-1)*length(beta)+k] <- mat.t.stat[i,k] > threshold | mat.t.stat[i,k] < -threshold
  }
  
  for (k in 1:length(nbv2)){
    vec.flag.type1.2[(i-1)*length(nbv2)+k] <- mat.t.stat2[i,k] > threshold | mat.t.stat2[i,k] < -threshold
  }
  for (k in 1:length(nbv3)){
    vec.flag.type1.3[(i-1)*length(nbv3)+k] <- mat.t.stat3[i,k] > threshold | mat.t.stat3[i,k] < -threshold
  }
  for (k in 1:length(nbv4)){
    vec.flag.type1.4[(i-1)*length(nbv4)+k] <- mat.t.stat4[i,k] > threshold | mat.t.stat4[i,k] < -threshold
  }
  for (k in 1:length(nbv5)){
    vec.flag.type1.5[(i-1)*length(nbv5)+k] <- mat.t.stat5[i,k] > threshold | mat.t.stat5[i,k] < -threshold
  }
  for (k in 1:length(nbv6)){
    vec.flag.type1.6[(i-1)*length(nbv6)+k] <- mat.t.stat6[i,k] > threshold | mat.t.stat6[i,k] < -threshold
  }
  for (k in 1:length(nbv7)){
    vec.flag.type1.7[(i-1)*length(nbv7)+k] <- mat.t.stat7[i,k] > threshold | mat.t.stat7[i,k] < -threshold
  }
  for (k in 1:length(nbv8)){
    vec.flag.type1.8[(i-1)*length(nbv8)+k] <- mat.t.stat8[i,k] > threshold | mat.t.stat8[i,k] < -threshold
  }
}
mat.flag.beta <- matrix(vec.flag.beta, ncol = length(beta), byrow = T)
mat.flag.type1 <- matrix(vec.flag.type1, ncol = length(beta), byrow = T)
mat.flag.type1.2 <- matrix(vec.flag.type1.2, ncol = length(nbv2), byrow = T)
mat.flag.type1.3 <- matrix(vec.flag.type1.3, ncol = length(nbv3), byrow = T)
mat.flag.type1.4 <- matrix(vec.flag.type1.4, ncol = length(nbv4), byrow = T)
mat.flag.type1.5 <- matrix(vec.flag.type1.5, ncol = length(nbv5), byrow = T)
mat.flag.type1.6 <- matrix(vec.flag.type1.6, ncol = length(nbv6), byrow = T)
mat.flag.type1.7 <- matrix(vec.flag.type1.7, ncol = length(nbv7), byrow = T)
mat.flag.type1.8 <- matrix(vec.flag.type1.8, ncol = length(nbv8), byrow = T)
#=======
object <- list(toN=1:N, 
               #
               mat.beta=mat.beta,
               mat.beta2 = mat.beta2,
               mat.beta3 = mat.beta3,
               mat.beta4 = mat.beta4,
               mat.beta5 = mat.beta5,
               mat.beta6 = mat.beta6,
               mat.beta7 = mat.beta7,
               mat.beta8 = mat.beta8,
               #
               vec.err.var = vec.err.var, 
               vec.err.var2 = vec.err.var2, 
               vec.err.var3 = vec.err.var3, 
               vec.err.var4 = vec.err.var4, 
               vec.err.var5 = vec.err.var5, 
               vec.err.var6 = vec.err.var6, 
               vec.err.var7 = vec.err.var7, 
               vec.err.var8 = vec.err.var8, 
               #
               mat.flag.beta = mat.flag.beta, 
               #
               mat.flag.type1 = mat.flag.type1, 
               mat.flag.type1.2 = mat.flag.type1.2,
               mat.flag.type1.3 = mat.flag.type1.3,
               mat.flag.type1.4 = mat.flag.type1.4,
               mat.flag.type1.5 = mat.flag.type1.5,
               mat.flag.type1.6 = mat.flag.type1.6,
               mat.flag.type1.7 = mat.flag.type1.7,
               mat.flag.type1.8 = mat.flag.type1.8,
               #
               mat.t.stat = mat.t.stat,
               #
               rsq.adj = rsq.adj,
               rsq.adj2 = rsq.adj2,
               rsq.adj3 = rsq.adj3,
               rsq.adj4 = rsq.adj4,
               rsq.adj5 = rsq.adj5,
               rsq.adj6 = rsq.adj6,
               rsq.adj7 = rsq.adj7,
               rsq.adj8 = rsq.adj8,
               #
               aic = aic,
               aic2 = aic2,
               aic3 = aic3,
               aic4 = aic4,
               aic5 = aic5,
               aic6 = aic6,
               aic7 = aic7,
               aic8 = aic8,
               #
               bic = bic,
               bic2 = bic2,
               bic3 = bic3,
               bic4 = bic4,
               bic5 = bic5,
               bic6 = bic6,
               bic7 = bic7,
               bic8 = bic8,
               #
               mlcp2 = mlcp2,
               mlcp3 = mlcp3,
               mlcp4 = mlcp4,
               mlcp5 = mlcp5,
               mlcp6 = mlcp6,
               mlcp7 = mlcp7,
               mlcp8 = mlcp8
               )
# names(object) <- c("vecN", "betahat1", "betahat2", "betahat3", "betahat4",  "betahat1_2", "betahat2_2", "betahat4_2", "vebh1", "vebh2", "vebh3", "vebh4", "err.varh", "err.varh2", "prob.flag.beta1", "prob.flag.beta2", "prob.flag.beta3", "prob.flag.beta4", "prob.flag.type1.bh1", "prob.flag.type1.bh2", "prob.flag.type1.bh3", "prob.flag.type1.bh4","prob.flag.type1.bh1_2", "prob.flag.type1.bh2_2", "prob.flag.type1.bh4_2", "t.stathat1", "t.stathat2", "t.stathat3", "t.stathat4")
# Save an object to a RData
class(object)
length(object)
save(object, file = "mcff5_2008_2018_mutiplemodel.RData")
#save(object, file = "mcff5_2008_2018.RData")
# Restore the object
#load(file = "mcff.RData")
#========

