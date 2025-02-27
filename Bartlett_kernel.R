# d. m. chakrabarti
# applied econometrics for macro
# Problem set 1
# problem number 6b

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

qT <- 12

for(j in 1:qT) {
  gamma_j <- matrix(0, ncol = k, nrow = k)
  for(t in (j+1):T_obs) {
    gamma_j <- gamma_j + t(residVec[t] * Xmat[t, , drop = FALSE]) %*% (residVec[t-j] * Xmat[t-j, , drop = FALSE])
  }
  gamma_j <- gamma_j / T_obs
  
  # Bartlett weight: weight = 1 - j/qT
  weight <- 1 - (j / qT)
  
  gammaSum <- gammaSum + weight * (gamma_j + t(gamma_j))
}

# (X'X)^{-1}
bread <- solve(t(Xmat) %*% Xmat)

# covariance matrix using Bartlett weights:
robustVar_Bartlett <- bread %*% (gamma0 + gammaSum) %*% bread

# robust standard errors:
robustSE_Bartlett <- sqrt(diag(robustVar_Bartlett))
print("Robust standard errors (Bartlett kernel, qT=12):")
print(robustSE_Bartlett)

# wald Test
# null hypothesis
betaHat <- coef(model)
betaNull <- c(0, 1)

diffBeta <- betaHat - betaNull

waldStat <- as.numeric(t(diffBeta) %*% solve(robustVar_Bartlett) %*% diffBeta)
cat("Wald statistic =", waldStat, "\n")

# p-value from chi-squared distribution with 2 degrees of freedom:
p_value <- 1 - pchisq(waldStat, df = 2)
cat("p-value =", p_value, "\n")
