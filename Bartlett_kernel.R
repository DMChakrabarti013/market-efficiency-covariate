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

nw_se <- coeftest(model, vcov = NeweyWest(model, lag = 12, prewhite = FALSE))

print(nw_se)

vcov_newWest <- NeweyWest(model, lag = 12, prewhite = FALSE)

# wald Test
# null hypothesis
betaHat <- coef(model)
betaNull <- c(0, 1)

diffBeta <- betaHat - betaNull

waldStat <- as.numeric(t(diffBeta) %*% solve(vcov_newWest) %*% diffBeta)
cat("Wald statistic =", waldStat, "\n")

# p-value from chi-squared distribution with 2 degrees of freedom:
p_value <- 1 - pchisq(waldStat, df = 2)
cat("p-value =", p_value, "\n")
