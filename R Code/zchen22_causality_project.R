# Load the libraries 
#install.packages("vars")
#install.packages("urca")
#install.packages("pcalg")
library(vars)
library(urca)
library(pcalg)

# Read the input data 
data = read.csv(file="../Input Data/data.csv",head=TRUE,sep=",")

# Build a VAR model 
var <- VAR(data, type="const")

# Extract the residuals from the VAR model 
var_residuals <- resid(var)

# Check for stationarity using the Augmented Dickey-Fuller test 
ur.df(var_residuals[,1],type="none",lags=0)
ur.df(var_residuals[,2],type="none",lags=0)
ur.df(var_residuals[,3],type="none",lags=0)
# The interpretation of the ur.df test value can be done using critical value 0.05 which is qnorm(0.05/2) = -1.960
# Since the values obtains from the test do not exceed -1.960, the residuals are stationary

# Check whether the variables follow a Gaussian distribution
ks.test(var_residuals[,1],"pnorm",mean=mean(var_residuals[,1]),sd = sd(var_residuals[,1]))
ks.test(var_residuals[,2],"pnorm",mean=mean(var_residuals[,2]),sd = sd(var_residuals[,2]))
ks.test(var_residuals[,3],"pnorm",mean=mean(var_residuals[,3]),sd = sd(var_residuals[,3]))
# from the p-value, they are not normal distributions

# PC algorithm
suffStat=list(C=cor(var_residuals), n=1000)
pc_fit <- pc(suffStat, indepTest=gaussCItest, alpha=0.1, labels=colnames(data), skel.method="original")
plot(pc_fit, main="PC Output")

# LiNGAM algorithm
lingam_fit <- LINGAM(var_residuals)
show(lingam_fit)

# Write the residuals to a csv file to build causal graphs using Tetrad software
write.csv(var_residuals, file = "residuals.csv")
