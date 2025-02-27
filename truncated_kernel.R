# d. m. chakrabarti
# applied econometrics for macro
# Problem set 1
# problem number 6a

########## loading packagaes #############
library(readxl)
library(lmtest)
library(sandwich)

#############################

data <- read_excel('/Users/dwaipayanchakrabarti/Downloads/yen.xlsx')

# compute logarithms of rates

data$S <- log(data[[2]]) # log spot rate
data$F <- log(data[[3]]) # log 30-day forward rate
data$S30 <- log(data[[4]]) # log spot rate 30 days later

# annualized differences (times 1200)
data$Y <- 1200 * (data$S30 - data$S) # dependent variable
data$X <- 1200 * (data$F - data$S) # independent variable

# OLS regression
model <- lm(Y ~ X, data = data)
summary(model)

# robust standard errors using truncated kernel (qt = 4)
# residuals and design matrix
residVec <- model$residuals
Xmat <- model.matrix(model)
T_obs <- nrow(Xmat)

gamma0 <- t(Xmat) %*% ( (residVec^2) * Xmat ) / T_obs # gamma 0
k <- ncol(Xmat)
gammaSum <- matrix(0, ncol = k, nrow = k)

qT <- 4

for(j in 1:qT) {
  gamma_j <- matrix(0, ncol = k, nrow = k)
  for(t in (j+1):T_obs) {
    gamma_j <- gamma_j + t(residVec[t] * Xmat[t, , drop = FALSE]) %*% (residVec[t-j] * Xmat[t-j, , drop = FALSE])
  }
  gamma_j <- gamma_j / T_obs
  gammaSum <- gammaSum + gamma_j + t(gamma_j)
}

bread <- solve(t(Xmat) %*% Xmat)

robustVar <- bread %*% gammaSum %*% bread

robustSE <- sqrt(diag(robustVar))
print("Robust standard errors (truncated kernel, qT=4):")
print(robustSE)

# wald test
betaHat <- coef(model)
betaNull <- c(0,1)

diffBeta <- betaHat - betaNull

waldStat <- as.numeric(t(diffBeta) %*% solve(robustVar) %*% diffBeta)
cat("Wald statistic =", waldStat, "\n")

# p-value
p_value <- 1 - pchisq(waldStat, df = 2)
cat("p-value =", p_value, "\n")
