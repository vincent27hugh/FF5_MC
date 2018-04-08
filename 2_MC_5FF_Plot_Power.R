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
#===========================================
## Construct X data frame
# X data: rmrf, smb, hml
X <- data.matrix(cbind(rep(1, dim(mydata2)[1]), 
                       rmrf, smb, hml, rmw, cma))
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
load(file = "mcff5_2008_2018_power.RData")
attach(power)
class(power)
length(power)
#===========================
# Significance level of t-test
alpha <- 0.05
N=tail(toN,1)
#=======================
temp <- rep(NA, dim(seq.beta)[1]*length(beta))
prob.reject <- matrix(temp, nrow = dim(seq.beta)[1], ncol = length(beta)) 
prob.reject[,1] <- colSums(mat.flag.reject1)/N
prob.reject[,2] <- colSums(mat.flag.reject2)/N
prob.reject[,3] <- colSums(mat.flag.reject3)/N
prob.reject[,4] <- colSums(mat.flag.reject4)/N
prob.reject[,5] <- colSums(mat.flag.reject5)/N
prob.reject[,6] <- colSums(mat.flag.reject6)/N
#=======================================================================================
foo = expression(hat(beta)[1], hat(beta)[2],hat(beta)[3],
                 hat(beta)[4],hat(beta)[5],hat(beta)[6])
# Power
for (i in 1:length(beta)){
  mypath <- file.path(getwd(),"Figure",
                      paste("iv_power_",i, ".jpeg", sep = ""))
  jpeg(mypath, height = 8, width = 8, 
       units = 'in', res = 300)
  plot(seq.beta[,i]-beta[i], prob.reject[,i], 
       type="b", xlab = "#Replication",
       ylab = "Power of Test", 
       main=foo[[i]], ylim = c(0,1))
  abline(h=0.05, col = "blue", lty = 2)
  dev.off()
}
# in one graph
mypath <- file.path(getwd(),"Figure")
mydimnames <- list(seq(from = -(Ns-1)/2, to = (Ns-1)/2, by = step), name=c("beta_1", "beta_2", "beta_3", 
                                "beta_4", "beta_5", "beta_6"))
dimnames(prob.reject)=mydimnames
tempdata <- as.data.frame(as.table(prob.reject))
library(ggplot2)
ggplot(data = tempdata, aes(x = Var1, y=Freq, group = name)) + 
  geom_line(aes(linetype = name, color = name), size = 1) + 
  xlab(expression(hat(beta)~"-"~beta^"*")) + 
  ylab("Power of test") + 
  ylim(0,1) + 
  geom_hline(yintercept = .05, color = "blue", linetype = "dotted")+
  theme_light() +
  theme(axis.text.x = element_text(face="bold", color="#993333", 
                                   size=10, angle=45),
        axis.text.y = element_text(face="bold", color="#993333", 
                                   size=10, angle=45))
ggsave("iv_power.jpeg", plot = last_plot(), path = mypath, width = 8, height = 8)
#=============================




