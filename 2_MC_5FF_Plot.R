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
lo_30.ffregression <- lm(lo_30.xcess ~ 
                           rmrf + smb + hml + rmw + cma)
# Print summary of regression results
lo_30.sum <- summary(lo_30.ffregression)
#
sink("output_lo_30_sum.txt")
print(lo_30.sum)
sink()
#===========================================
## Construct X data frame
# X data: rmrf, smb, hml
X <- data.matrix(cbind(rep(1, dim(mydata2)[1]), 
                       rmrf, smb, hml, rmw, cma))

# Coloum Summary
summaryfun <- function(x){
  nr = dim(x)[2]
  nc = 6
  temp <- rep(NA, nr*nc)
  out <- matrix(temp, nrow = nr, ncol = nc) 
  for (i in 1:dim(X)[2]){
    out[i,] <- c(N = length(x[,i]),
                 Mean = mean(x[,i]),
                 Median = median(x[,i]),
                 StdDev = sd(x[,i]),
                 Min = min(x[,i]),
                 Max = max(x[,i]))
  }
  out1 <- as.data.frame(out)
  return(out1)
}
des.stat <- summaryfun(X)
names(des.stat) <- c("N",
                     "Mean", 
                     "Median",
                     "St.Dev",
                     "Min",
                     "Max")
row.names(des.stat) <- c("1",
                             "RmRf", 
                             "SMB",
                             "HML",
                             "RMW",
                             "CMA")


# stargazer
#install.packages("stargazer") #Use this to install it, do this only once
library(stargazer)
stargazer(as.data.frame(X), type = "text", title="Descriptive statistics", digits=4, out="output_star_X.txt")
#
sink("output_X.txt")
print(des.stat)
sink()
# # obervations
M <- dim(X)[1]
# True beta
beta <- lo_30.sum$coefficients[,1]
#
mu <- 0
sig2e <- var(lo_30.xcess)
# True sigma: std of distribution of error term
sde <- sqrt(sig2e)
#==================================
# Load runned data
#load(file = "mcff5_2008_2018.RData")
load(file = "mcff5_2008_2018_mutiplemodel.RData")
attach(object)
class(object)
length(object)
#===========================
# Significance level of t-test
alpha <- 0.05
# Ordinal number for second model with omitted variable
temp <- 1:length(beta)
nbv2 <- temp[-2]
nbv3 <- temp[-3]
nbv4 <- temp[-4]
nbv5 <- temp[-5]
nbv6 <- temp[-6]
nbv7 <- temp[1:4]
nbv8 <- temp[1:2]
#=======================
#vecN = c(10,50,100)
vecN = c(10,50,100,500,1e3,5e3,seq(from = 1e4, to = 1e5, by = 1e4))
mat.betahat <- matrix(0, 
                      nrow = length(vecN), 
                      ncol = length(beta))
mat.betahat2 <- matrix(0, 
                       nrow = length(vecN), 
                       ncol = length(nbv2))

mat.vebh <- matrix(0, 
                   nrow = length(vecN), 
                   ncol = length(beta))
vec.err.varh <- numeric(length(vecN))
vec.prob.err.varh <- numeric(length(vecN))
vec.err.varh2 <- numeric(length(vecN))
vec.prob.err.varh2 <- numeric(length(vecN))
prob.flag.beta <- matrix(0, 
                         nrow = length(vecN), 
                         ncol = length(beta))
prob.flag.type1 <- matrix(0, 
                          nrow = length(vecN), 
                          ncol = length(beta))
prob.flag.type1.2 <- matrix(0, 
                            nrow = length(vecN), 
                            ncol = length(beta)-1)
mat.t.stathat <- matrix(0, 
                        nrow = length(vecN), 
                        ncol = length(beta))
for (j in 1:length(vecN)){
  N <- vecN[j]
  # tolerance
  epsilon <- 0.05
  # \hat{\beta}
  mat.betahat[j,] <- colMeans(mat.beta[1:N,])
  #
  mat.betahat2[j,] <- colMeans(mat.beta2[1:N,])
  # \variance of betahat
  mat.vebh[j,] <- apply(mat.beta[1:N,], 2, var)
  # error variance
  vec.err.varh[j] <- mean(vec.err.var[1:N])
  # P(abs(err.varh - sig2e) >= epsilon)
  vec.prob.err.varh[j] <- sum(abs(vec.err.var[1:N] - sig2e) 
                              >= epsilon)/N
  # error variance for model 2
  vec.err.varh2[j] <- mean(vec.err.var2[1:N])
  # P(abs(err.varh - sig2e) >= epsilon)
  vec.prob.err.varh2[j] <- sum(abs(vec.err.var2[1:N] - sig2e) 
                               >= epsilon)/N
  # Probability of true beta is in confidence interval of estimated beta
  prob.flag.beta[j,] <- colSums(mat.flag.beta[1:N,])/N
  # Prob of Type I
  prob.flag.type1[j,] <- colSums(mat.flag.type1[1:N,])/N
  # Prob of Type I
  prob.flag.type1.2[j,] <- colSums(mat.flag.type1.2[1:N,])/N
  # 
  mat.t.stathat[j,] <- colMeans(mat.t.stat[1:N,])
}
#==============================================
# Model Selection Table
v.rsq.adj <- c(mean(rsq.adj),mean(rsq.adj2),mean(rsq.adj3),mean(rsq.adj4),
           mean(rsq.adj5),mean(rsq.adj6),mean(rsq.adj7),mean(rsq.adj8))
v.aic <- c(mean(aic),mean(aic2),mean(aic3),mean(aic4),
           mean(aic5),mean(aic6),mean(aic7),mean(aic8))
v.bic <- c(mean(bic),mean(bic2),mean(bic3),mean(bic4),
           mean(bic5),mean(bic6),mean(bic7),mean(bic8))
v.mlcp <- c(NaN,mean(mlcp2),mean(mlcp3),mean(mlcp4),
           mean(mlcp5),mean(mlcp6),mean(mlcp7),mean(mlcp8))
model_select <- as.data.frame(cbind(v.rsq.adj,
                              v.aic,
                              v.bic,
                              v.mlcp))
names(model_select) <- c("Adjusted R-squared", 
                       "AIC",
                       "BIC",
                       "Mallow's Cp")
row.names(model_select) <- c("Full Model",
                             "Model 2",
                             "Model 3",
                             "Model 4",
                             "Model 5",
                             "Model 6",
                             "Model 7",
                             "Model 8")
sink("output_model_select.txt")
print(model_select)
sink()
#==============================================
# Table of Probability of Type I error for different model
temp <- rep(NA, 8*length(beta))
prob.type1 <- matrix(temp,nrow = 8, ncol = length(beta)) 
prob.type1[1,]<- colSums(mat.flag.type1)/tail(toN,1)
prob.type1[2,nbv2] <- colSums(mat.flag.type1.2)/tail(toN,1)
prob.type1[3,nbv3] <- colSums(mat.flag.type1.3)/tail(toN,1)
prob.type1[4,nbv4] <- colSums(mat.flag.type1.4)/tail(toN,1)
prob.type1[5,nbv5] <- colSums(mat.flag.type1.5)/tail(toN,1)
prob.type1[6,nbv6] <- colSums(mat.flag.type1.6)/tail(toN,1)
prob.type1[7,nbv7] <- colSums(mat.flag.type1.7)/tail(toN,1)
prob.type1[8,nbv8] <- colSums(mat.flag.type1.8)/tail(toN,1)
prob.type1 <- as.data.frame(prob.type1)
names(prob.type1) <- c("beta1", 
                         "beta2", 
                         "beta3", 
                         "beta4", 
                         "beta5", 
                         "beta6")
row.names(prob.type1) <- c("Full Model",
                             "Model 2",
                             "Model 3",
                             "Model 4",
                             "Model 5",
                             "Model 6",
                             "Model 7",
                             "Model 8")
sink("output_prob_typeI.txt")
print(prob.type1)
sink()
#==============================================
# Difference between simulated value and true value
temp <- rep(NA, 
            8*(length(beta) + 1))
tab.bias <- matrix(temp,nrow = 8, 
                   ncol = length(beta) + 1) 
tab.bias[1,] <- c(
  # Bias of OLS estmator of beta
  abs(beta-colMeans(mat.beta)),
  # Bias of OLS estimator of Error Variance
  abs(sig2e-mean(vec.err.var)))
tab.bias[2,c(nbv2,7)] <- c(
  # Bias of OLS estmator of beta
  abs(beta[nbv2]-colMeans(mat.beta2)),
  # Bias of OLS estimator of Error Variance
  abs(sig2e-mean(vec.err.var2)))
tab.bias[3,c(nbv3,7)] <- c(
  # Bias of OLS estmator of beta
  abs(beta[nbv3]-colMeans(mat.beta3)),
  # Bias of OLS estimator of Error Variance
  abs(sig2e-mean(vec.err.var3)))
tab.bias[4, c(nbv4,7)] <- c(
  # Bias of OLS estmator of beta
  abs(beta[nbv4]-colMeans(mat.beta4)),
  # Bias of OLS estimator of Error Variance                    
  abs(sig2e-mean(vec.err.var4)))
tab.bias[5, c(nbv5,7)] <- c(
  # Bias of OLS estmator of beta
  abs(beta[nbv5]-colMeans(mat.beta5)),
  # Bias of OLS estimator of Error Variance                    
  abs(sig2e-mean(vec.err.var5)))
tab.bias[6, c(nbv6,7)] <- c(
  # Bias of OLS estmator of beta
  abs(beta[nbv6]-colMeans(mat.beta6)),
  # Bias of OLS estimator of Error Variance                    
  abs(sig2e-mean(vec.err.var6)))
tab.bias[7, c(nbv7,7)] <- c(
  # Bias of OLS estmator of beta
  abs(beta[nbv7]-colMeans(mat.beta7)),
  # Bias of OLS estimator of Error Variance                    
  abs(sig2e-mean(vec.err.var7)))
tab.bias[8, c(nbv8,7)] <- c(
  # Bias of OLS estmator of beta
  abs(beta[nbv8]-colMeans(mat.beta8)),
  # Bias of OLS estimator of Error Variance                    
  abs(sig2e-mean(vec.err.var8)))
tab.bias <- as.data.frame(tab.bias)
names(tab.bias) <- c("beta1", 
                       "beta2", 
                       "beta3", 
                       "beta4", 
                       "beta5", 
                       "beta6",
                     "Error Variance")
row.names(tab.bias) <- c("Full Model",
                           "Model 2",
                           "Model 3",
                           "Model 4",
                           "Model 5",
                           "Model 6",
                           "Model 7",
                           "Model 8")
sink("output_table_bias.txt")
print(tab.bias)
sink()
#=======================================================================================#
# Plot of N = 10^4
#=======================================================================================#
# In Distribution pic, bulue line is true beta, red line is betahat
foo = expression(hat(beta)[1], hat(beta)[2],hat(beta)[3],
                 hat(beta)[4],hat(beta)[5],hat(beta)[6])
for (i in 1:length(beta)){
  mypath <- file.path(getwd(),"Figure",
                      paste("density_beta_",i, ".jpeg", sep = ""))
  jpeg(mypath, height = 8, width = 8,
       units = 'in', res = 300)
  # beta_{i}
  plot(density(mat.beta[,i]), 
       xlab = foo[[i]], 
       main = "Distribution")
  curve(dnorm(x, mat.betahat[length(vecN),i], 
              sqrt(mat.vebh[length(vecN),i])), 
        col="red", add=TRUE)
  abline(v=mat.betahat[length(vecN),i], 
         col = "red", lty = 2)
  abline(v=beta[i], col = "blue", lty = 2)
  # "T" for "True", "S" for "Simulated"
  legend("topright", legend=c("T", "S"), 
         lty=1, col=c("red", "black"))
  dev.off()
  ##
  mypath <- file.path(getwd(),"Figure",
                      paste("hist_beta_",i, ".jpeg", sep = ""))
  jpeg(mypath, height = 8, width = 8, 
       units = 'in', res = 300)
  hist(mat.beta[,i], prob=TRUE, 
       xlab = foo[[i]], main = "Histogram")
  curve(dnorm(x, mat.betahat[length(vecN),i], 
              sd=sqrt(mat.vebh[length(vecN),i])), 
        col="red", add=TRUE)
  # "T" for "True", "S" for "Simulated"
  legend("topright", legend=c("T", "S"), 
         lty=1, col=c("red", "black"))
  dev.off()
}

#=======================================================================================#
# Answer six questions
#=======================================================================================#
foo = expression(hat(beta)[1], hat(beta)[2],hat(beta)[3],
                 hat(beta)[4],hat(beta)[5],hat(beta)[6])
# i) The O.L.S. estimators of the unknown coefficients and error variance are unbiased
for (i in 1:length(beta)){
  mypath <- file.path(getwd(),"Figure",
                      paste("i_OLS_betahat_",i, ".jpeg", sep = ""))
  jpeg(mypath, height = 8, width = 8, 
       units = 'in', res = 300)
  plot(vecN, mat.betahat[,i], type="b",
       xlab = "#Replication", ylab = foo[[i]], 
       main="Convergence Graph")
  abline(h=beta[i], col = "blue", lty = 2)
  dev.off()
}
# In one graph
mypath <- file.path(getwd(),"Figure")
mydimnames1 <- list(vecN, 
                    name=c("beta_1", "beta_2", "beta_3", 
                           "beta_4", "beta_5", "beta_6"))
dimnames(mat.betahat)=mydimnames1
#
mydimnames2 <- list(vecN, 
                    name=c("beta_1", "beta_2", "beta_3", 
                           "beta_4", "beta_5", "beta_6"))
trueb <- t(replicate(length(vecN),beta))
dimnames(trueb)=mydimnames2
#
method = gl(2,(length(beta))*length(vecN),
            (length(beta))*length(vecN)*2, 
            labels = c("Simulated", "True"))
tempdata <- as.data.frame(as.table(cbind(mat.betahat,trueb)))
colnames(tempdata) <- c("N","beta","Value")
tempdata$Method <- method
#
library(ggplot2)
ggplot(data = tempdata, aes(x = N, 
                            y = Value,
                            linetype = Method,
                            colour = beta,
                            group = interaction(beta,Method))) + 
  geom_line(size = 1) + 
  xlab("#Replication") + 
  ylab(expression(hat(beta))) + 
  # color ref: http://sape.inf.usi.ch/quick-reference/ggplot2/colour
  # linetype ref: http://sape.inf.usi.ch/quick-reference/ggplot2/linetype
  theme_light() +
  theme(axis.text.x = element_text(face="bold", color="#993333", 
                                   size=10, angle=45),
        axis.text.y = element_text(face="bold", color="#993333", 
                                   size=10, angle=45))
ggsave("i_prob_OLS_betah.jpeg", 
       plot = last_plot(), 
       path = mypath, 
       width = 8, 
       height = 8)

# In one graph of betahat - beta
mypath <- file.path(getwd(),"Figure")
mydimnames1 <- list(vecN, 
                    name=c("beta_1", "beta_2", "beta_3", 
                           "beta_4", "beta_5", "beta_6"))
dimnames(mat.betahat)=mydimnames1
#
mydimnames2 <- list(vecN, 
                    name=c("beta_1", "beta_2", "beta_3", 
                           "beta_4", "beta_5", "beta_6"))
trueb <- t(replicate(length(vecN),beta))
dimnames(trueb)=mydimnames2
#
tempdata <- as.data.frame(as.table(mat.betahat-trueb))
colnames(tempdata) <- c("N","betahat_min_beta","Value")
#
library(ggplot2)
ggplot(data = tempdata, aes(x = N, 
                            y = Value,
                            colour = betahat_min_beta,
                            group = betahat_min_beta)) + 
  geom_line(size = 1) + 
  geom_hline(yintercept = 0, color = "blue", linetype="dotted") +
  xlab("#Replication") + 
  ylab(expression(hat(beta)~"-"~beta)) + 
  # color ref: http://sape.inf.usi.ch/quick-reference/ggplot2/colour
  # linetype ref: http://sape.inf.usi.ch/quick-reference/ggplot2/linetype
  theme_light() +
  theme(axis.text.x = element_text(face="bold", color="#993333", 
                                   size=10, angle=45),
        axis.text.y = element_text(face="bold", color="#993333", 
                                   size=10, angle=45),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle(expression("Plot of "~hat(beta)~"-"~beta))
  
ggsave("i_prob_OLS_betah2.jpeg", 
       plot = last_plot(), 
       path = mypath, 
       width = 8, 
       height = 8)
#================
# beta - betahat ~ N(0, sigma^2 (X^T X)^-1)
var.beta <- sig2e*solve(t(X)%*%X)
foo2 = expression("var of "~hat(beta)[1], 
                  "var of "~hat(beta)[2],
                  "var of "~hat(beta)[3],
                  "var of "~hat(beta)[4],
                  "var of "~hat(beta)[5],
                  "var of "~hat(beta)[6])
for (i in 1:length(beta)){
  mypath <- file.path(getwd(),"Figure",
                      paste("i_var_betahat_",i, ".jpeg", sep = ""))
  jpeg(mypath, height = 8, width = 8, 
       units = 'in', res = 300)
  plot(vecN, mat.vebh[,i], 
       type="b", 
       xlab = "#Replication", 
       ylab = foo2[[i]], 
       main="Convergence Graph")
  abline(h=var.beta[i,i], col = "blue", lty = 2)
  dev.off()
}

# Error Variance
# Unbised estimater of sigma^2
mypath <- file.path(getwd(),
                    "Figure","i_var_err.jpeg")
jpeg(mypath, height = 8, width = 8, 
     units = 'in', res = 300)
plot(vecN, vec.err.varh, 
     type="b", 
     xlab = "#Replication", 
     ylab = "Error Var.")
abline(h=sig2e, col = "blue", lty = 2)
dev.off()
#=======================================================================================
# ii) The “correct” meaning of a 100(1‐α)% confidence interval of an unknown coefficient
for (i in 1:length(beta)){
  mypath <- file.path(getwd(),"Figure",
                      paste("ii_prob_beta_",i, ".jpeg", sep = ""))
  jpeg(mypath, height = 8, width = 8, 
       units = 'in', res = 300)
  plot(vecN, prob.flag.beta[,i], 
       type="b", 
       xlab = "#Replication", 
       ylab = "Prob(CI contains beta)", 
       main=foo[[i]], 
       ylim = c(0.6,1))
  abline(h=0.95, col = "blue", lty = 2)
  dev.off()
}
# in one graph
mypath <- file.path(getwd(),"Figure")
mydimnames <- list(vecN, name=c("beta_1", "beta_2", "beta_3", 
                                "beta_4", "beta_5", "beta_6"))
dimnames(prob.flag.beta)=mydimnames
tempdata <- as.data.frame(as.table(prob.flag.beta))
library(ggplot2)
ggplot(data = tempdata, 
       aes(x = Var1, y=Freq, group = name)) + 
  geom_line(aes(linetype=name, color = name), size=1) + 
  xlab("#Replication") + 
  ylab(expression("Pr(CI contains "~beta~")")) + 
  ylim(0.9,1) + 
  scale_y_continuous(labels = scales::percent) +
  geom_hline(yintercept = .95, color = "blue", linetype="dotted") +
  theme_light() +
  theme(axis.text.x = element_text(face="bold", color="#993333", 
                                   size=10, angle=45),
        axis.text.y = element_text(face="bold", color="#993333", 
                                   size=10, angle=45))
ggsave("ii_prob_beta.jpeg", plot = last_plot(), path = mypath, width = 8, height = 8)
#=======================================================================================
# iii) The significance level of the t test for testing a linear hypothesis concerning one or more coefficients is the probability of committing a Type I error
for (i in 1:length(beta)){
  mypath <- file.path(getwd(),"Figure",
                      paste("iii_prob_type1err_",i, ".jpeg", sep = ""))
  jpeg(mypath, height = 8, width = 8, 
       units = 'in', res = 300)
  plot(vecN, prob.flag.type1[,i], 
       type="b", xlab = "#Replication",
       ylab = "Prob(Type 1 Error)", 
       main=foo[[i]], ylim = c(0,1))
  abline(h=0.05, col = "blue", lty = 2)
  dev.off()
}
# in one graph
mypath <- file.path(getwd(),"Figure")
mydimnames <- list(vecN, name=c("beta_1", "beta_2", "beta_3", 
                                "beta_4", "beta_5", "beta_6"))
dimnames(prob.flag.type1)=mydimnames
tempdata <- as.data.frame(as.table(prob.flag.type1))
library(ggplot2)
ggplot(data = tempdata, aes(x = Var1, y=Freq, group = name)) + 
  geom_line(aes(linetype = name, color = name), size = 1) + 
  xlab("#Replication") + 
  ylab("Prob(Type 1 Error)") + 
  ylim(0,0.1) + 
  geom_hline(yintercept = .05, color = "blue", linetype = "dotted")+
  theme_light() +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_text(face="bold", color="#993333", 
                                   size=10, angle=45),
        axis.text.y = element_text(face="bold", color="#993333", 
                                   size=10, angle=45))
ggsave("iii_prob_type1err.jpeg", plot = last_plot(), path = mypath, width = 8, height = 8)
#=============================
# iv) The t test is unbiased
for (i in 1:length(beta)){
  mypath <- file.path(getwd(),"Figure",paste("iv_tstat_",i, ".jpeg", sep = ""))
  jpeg(mypath, height = 8, width = 8, units = 'in', res = 300)
  plot(vecN, mat.t.stathat[,i], type="b", xlab = "#Replication", ylab = "t-stat", main=foo[[i]])
  abline(h=0, col = "blue", lty = 2)
  dev.off()
}
# in one graph
mypath <- file.path(getwd(),"Figure")
mydimnames <- list(vecN, name=c("beta_1", "beta_2", "beta_3", "beta_4", "beta_5", "beta_6"))
dimnames(mat.t.stathat)=mydimnames
tempdata <- as.data.frame(as.table(mat.t.stathat))
library(ggplot2)
ggplot(data = tempdata, aes(x = Var1, y=Freq, group = name)) + 
  geom_line(aes(linetype = name, color = name), size = 1) + 
  xlab("#Replication") + 
  ylab("t-stat") + 
  geom_hline(yintercept = 0, color = "blue", linetype = "dotted") +
  theme_light() +
  theme(axis.text.x = element_text(face="bold", color="#993333", 
                                   size=10, angle=45),
        axis.text.y = element_text(face="bold", color="#993333", 
                                   size=10, angle=45))
ggsave("iv_tstat.jpeg", plot = last_plot(), path = mypath, width = 8, height = 8)
#=======================================================================================
# v) The result in Part iii) no longer holds if some relevant explanatory variables have been omitted from the model
foo3 = expression(hat(beta)[12], hat(beta)[32],hat(beta)[42],hat(beta)[52],hat(beta)[62])
for (i in 1:(length(beta)-1)){
  mypath <- file.path(getwd(),"Figure",paste("v_prob_type1err_",i, ".jpeg", sep = ""))
  jpeg(mypath, height = 8, width = 8, units = 'in', res = 300)
  plot(vecN, prob.flag.type1.2[,i], type="b", xlab = "#Replication", ylab = "Prob(Type 1 Error)", main=foo3[[i]], ylim = c(0,1))
  abline(h=0.05, col = "blue", lty = 2)
  dev.off()
}
# in one graph
mypath <- file.path(getwd(),"Figure")
mydimnames <- list(vecN, name=c("beta_12", "beta_32", "beta_42", "beta_52", "beta_62"))
dimnames(prob.flag.type1.2)=mydimnames
tempdata <- as.data.frame(as.table(prob.flag.type1.2))
library(ggplot2)
ggplot(data = tempdata, aes(x = Var1, y=Freq, group = name)) + 
  geom_line(aes(linetype = name, color = name), size = 1) + 
  xlab("#Replication") + 
  ylab("Prob(Type 1 Error)")+ 
  ylim(0,1) + 
  geom_hline(yintercept = 0.05, color = "blue", linetype = "dotted") +
  ggtitle("Model with omitting variable") +
  theme_light() +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_text(face="bold", color="#993333", 
                                   size=10, angle=45),
        axis.text.y = element_text(face="bold", color="#993333", 
                                   size=10, angle=45))
ggsave("v_prob_type1err.jpeg", plot = last_plot(), path = mypath, width = 8, height = 8)
#=============================
# vi) The estimator of say, the coefficient of X2, is no longer unbiased if the decision of whether to include X1 in the model is dependent on the outcome of a t test. Based on your findings, discuss the wider implications of “model selection” for statistical modeling and the lessons to be learnt for practitioners.
temp <- 1:length(beta)
nbv <- temp[-2]
for (i in 1:(length(beta)-1)){
  mypath <- file.path(getwd(),"Figure",paste("vi_batahat_",i, ".jpeg", sep = ""))
  jpeg(mypath, height = 8, width = 8, units = 'in', res = 300)
  plot(vecN, mat.betahat2[,i], type="b", xlab = "#Replication", ylab = foo3[[i]], main="Convergence Graph")
  abline(h=beta[nbv[i]], col = "blue", lty = 2)
  dev.off()
}
# One graph
mypath <- file.path(getwd(),"Figure")
#
mydimnames1 <- list(vecN, 
                    name=c("beta_12", "beta_32", 
                           "beta_42", "beta_52", "beta_62"))
dimnames(mat.betahat2)=mydimnames1
#
mydimnames2 <- list(vecN, 
                    name=c("beta_12", "beta_32", 
                           "beta_42", "beta_52", "beta_62"))
trueb <- t(replicate(length(vecN),beta[nbv]))
dimnames(trueb)=mydimnames2
#
method = gl(2,(length(beta)-1)*length(vecN),
            (length(beta)-1)*length(vecN)*2, 
            labels = c("Simulated", "True")) 
tempdata <- as.data.frame(as.table(cbind(mat.betahat2,trueb)))
colnames(tempdata) <- c("N","beta","Value")
tempdata$Method <- method

library(ggplot2)
ggplot(data = tempdata, aes(x = N, y = Value, 
                            linetype = Method, colour = beta, 
                            group = interaction(beta, Method))) +
  geom_line(size = 1) + 
  xlab("#Replication") + 
  ylab(expression(hat(beta))) + 
  ggtitle("Model with omitting variable") +
  theme_light() +
  theme(axis.text.x = element_text(face="bold", color="#993333", 
                                   size=10, angle=45),
        axis.text.y = element_text(face="bold", color="#993333", 
                                   size=10, angle=45))
ggsave("vi_prob_batahat.jpeg", plot = last_plot(), path = mypath, width = 8, height = 8)

# In one graph of betahat - beta
mypath <- file.path(getwd(),"Figure")
mydimnames1 <- list(vecN, 
                    name=c("beta_12", "beta_32", 
                           "beta_42", "beta_52", "beta_62"))
dimnames(mat.betahat2)=mydimnames1
#
mydimnames2 <- list(vecN, 
                    name=c("beta_12", "beta_32", 
                           "beta_42", "beta_52", "beta_62"))
trueb <- t(replicate(length(vecN),beta[nbv]))
dimnames(trueb)=mydimnames2
#
tempdata <- as.data.frame(as.table(mat.betahat2-trueb))
colnames(tempdata) <- c("N","betahat_min_beta","Value")
#
library(ggplot2)
ggplot(data = tempdata, aes(x = N, 
                            y = Value,
                            colour = betahat_min_beta,
                            group = betahat_min_beta)) + 
  geom_line(size = 1) + 
  geom_hline(yintercept = 0, color = "blue", linetype="dotted") +
  xlab("#Replication") + 
  ylab(expression(hat(beta)~"-"~beta)) + 
  # color ref: http://sape.inf.usi.ch/quick-reference/ggplot2/colour
  # linetype ref: http://sape.inf.usi.ch/quick-reference/ggplot2/linetype
  theme_light() +
  theme(axis.text.x = element_text(face="bold", color="#993333", 
                                   size=10, angle=45),
        axis.text.y = element_text(face="bold", color="#993333", 
                                   size=10, angle=45),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle(expression("Plot of "~hat(beta)~"-"~beta~"of omitting model"))

ggsave("vi_prob_batahat2.jpeg", 
       plot = last_plot(), 
       path = mypath, 
       width = 8, 
       height = 8)

# Error Variance
# Unbised estimater of sigma^2
mypath <- file.path(getwd(),"Figure","vi_var_err.jpeg")
jpeg(mypath, height = 8, width = 8, units = 'in', res = 300)
# bty for Box Style
plot(vecN, vec.err.varh2, type="b", xlab = "#Replication", ylab = "Error Var.", main = "Model of omitting")
abline(h=sig2e, col = "blue", lty = 2)
dev.off()



