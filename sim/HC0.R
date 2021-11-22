set.seed(12)

# adapted from https://stats.stackexchange.com/questions/110091/how-to-calculate-the-robust-standard-error-of-predicted-y-from-a-linear-regressi

# test data
d <- fread("../test/data/data-outlier.csv")
x <- d$x
y = d$y
X = cbind(1,x)

#Linear model
m = lm(y~X-1)
summary(m)
betahat = as.matrix(coef(m))

#Non-HC standard errors of fitted values
se.yhat = sqrt(diag(X%*%vcov(m)%*%t(X)))
se = predict(m,se.fit=TRUE)$se.fit
all.equal(se,se.yhat)#Matrix formula gives same result as se.fit option of predict method for lm's

#Now getting se's of fitted values based on HC-robust paramter covariance matrix
ehat = residuals(m)
vcov.HC = solve(t(X)%*%X) %*% t(X)%*%diag(ehat^2)%*%X %*% solve(t(X)%*%X)

se.yhat.HC = sqrt(diag(X%*%vcov.HC%*%t(X)))

#Showing that this method replicates regression parameter standard errors as given by the coeftest function
library(lmtest)
library(sandwich)
coeftest(m,vcov=vcovHC(m,type="HC0"))
sqrt(diag(vcov.HC))