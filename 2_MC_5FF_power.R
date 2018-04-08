rm(list=ls())
#graphics.off()
#=======================================================================================#
# load data
load(file="5factors_2008_2018.RData")
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
#===========================================
## Construct X data frame
# X data: rmrf, smb, hml
X <- data.matrix(cbind(rep(1, dim(mydata2)[1]), rmrf, smb, hml, rmw, cma))
# # obervations
M <- dim(X)[1]
# True beta
beta <- lo_30.sum$coefficients[,1]
#
mu <- 0
sig2e <- var(lo_30.xcess)
# True sigma: std of distribution of error term
sde <- sqrt(sig2e)
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
# Number of sequence for each true beta for the power simulation
Ns = 3
step = 0.1
seq.beta <- matrix(0, nrow = (Ns-1)/step+1, ncol = length(beta))
for (i in (1:length(beta))){
  seq.beta[,i] <- seq(from = beta[i]-(Ns-1)/2, to = beta[i] + (Ns-1)/2, by = step)
}
# beta1
mat.t.stat.power1 <- matrix(0,nrow = N, ncol = dim(seq.beta)[1])
mat.flag.reject1 <- matrix(0,nrow = N, ncol = dim(seq.beta)[1])
# beta2
mat.t.stat.power2 <- matrix(0,nrow = N, ncol = dim(seq.beta)[1])
mat.flag.reject2 <- matrix(0,nrow = N, ncol = dim(seq.beta)[1])
# beta3
mat.t.stat.power3 <- matrix(0,nrow = N, ncol = dim(seq.beta)[1])
mat.flag.reject3 <- matrix(0,nrow = N, ncol = dim(seq.beta)[1])
# beta4
mat.t.stat.power4 <- matrix(0,nrow = N, ncol = dim(seq.beta)[1])
mat.flag.reject4 <- matrix(0,nrow = N, ncol = dim(seq.beta)[1])
# beta5
mat.t.stat.power5 <- matrix(0,nrow = N, ncol = dim(seq.beta)[1])
mat.flag.reject5 <- matrix(0,nrow = N, ncol = dim(seq.beta)[1])
# beta6
mat.t.stat.power6 <- matrix(0,nrow = N, ncol = dim(seq.beta)[1])
mat.flag.reject6 <- matrix(0,nrow = N, ncol = dim(seq.beta)[1])

for (i in 1:N) {
  set.seed(27+6*i)
  Y <- X %*% beta + rnorm(M, mean=0, sd=sde)
  # Full Linear Regression
  ff<-lm(Y~X[,-1])
  # summary
  sum <- summary(ff)
  # For different t test of null hypothesis
  for (s in 1:dim(seq.beta)[1]){
    # beta 1
    mat.t.stat.power1[i,s] <- (sum$coefficients[1,1] - 
                                 seq.beta[s,1]) / sum$coefficients[1,2]
    # beta 2
    mat.t.stat.power2[i,s] <- (sum$coefficients[2,1] - 
                                 seq.beta[s,2]) / sum$coefficients[2,2]
    # beta 3
    mat.t.stat.power3[i,s] <- (sum$coefficients[3,1] - 
                                 seq.beta[s,3]) / sum$coefficients[3,2]
    # beta 4
    mat.t.stat.power4[i,s] <- (sum$coefficients[4,1] - 
                                 seq.beta[s,4]) / sum$coefficients[4,2]
    # beta 5
    mat.t.stat.power5[i,s] <- (sum$coefficients[5,1] - 
                                 seq.beta[s,5]) / sum$coefficients[5,2]
    # beta 6
    mat.t.stat.power6[i,s] <- (sum$coefficients[6,1] - 
                                 seq.beta[s,6]) / sum$coefficients[6,2]
  }
  
  # t_{1-\alpha/2, n-k}
  threshold <- qt(1-alpha/2, df = df.residual(ff))
  # whether reject null hypothesis: False - not reject; True - reject (Type I)
  # For different t test of null hypothesis
  for (s in 1:dim(seq.beta)[1]){
    # beta 1
    mat.flag.reject1[i,s] <- 
      mat.t.stat.power1[i,s] > threshold | mat.t.stat.power1[i,s] < -threshold
    # beta 2
    mat.flag.reject2[i,s] <- 
      mat.t.stat.power2[i,s] > threshold | mat.t.stat.power2[i,s] < -threshold
    # beta 3
    mat.flag.reject3[i,s] <- 
      mat.t.stat.power3[i,s] > threshold | mat.t.stat.power3[i,s] < -threshold
    # beta 4
    mat.flag.reject4[i,s] <- 
      mat.t.stat.power4[i,s] > threshold | mat.t.stat.power4[i,s] < -threshold
    # beta 5
    mat.flag.reject5[i,s] <- 
      mat.t.stat.power5[i,s] > threshold | mat.t.stat.power5[i,s] < -threshold
    # beta 6
    mat.flag.reject6[i,s] <- 
      mat.t.stat.power6[i,s] > threshold | mat.t.stat.power6[i,s] < -threshold
  }
}
#=======
power <- list(toN=1:N,
              mat.flag.reject1 = mat.flag.reject1,
              mat.flag.reject2 = mat.flag.reject2,
              mat.flag.reject3 = mat.flag.reject3,
              mat.flag.reject4 = mat.flag.reject4,
              mat.flag.reject5 = mat.flag.reject5,
              mat.flag.reject6 = mat.flag.reject6,
              Ns = Ns,
              step = step,
              seq.beta = seq.beta
)
# names(object) <- c("vecN", "betahat1", "betahat2", "betahat3", "betahat4",  "betahat1_2", "betahat2_2", "betahat4_2", "vebh1", "vebh2", "vebh3", "vebh4", "err.varh", "err.varh2", "prob.flag.beta1", "prob.flag.beta2", "prob.flag.beta3", "prob.flag.beta4", "prob.flag.type1.bh1", "prob.flag.type1.bh2", "prob.flag.type1.bh3", "prob.flag.type1.bh4","prob.flag.type1.bh1_2", "prob.flag.type1.bh2_2", "prob.flag.type1.bh4_2", "t.stathat1", "t.stathat2", "t.stathat3", "t.stathat4")
# Save an object to a RData
class(power)
length(power)
save(power, file = "mcff5_2008_2018_power.RData")
#========

