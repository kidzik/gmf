require(assertthat)

n = 1009
d = 5
X = matrix(rnorm(n*d),n,d)
trueBeta = rnorm(d,mean = 1)

trueBeta
set.seed(0)

# test regressions exact
Y = X%*%trueBeta
estBeta = gmf:::glm.basic(X,Y,NULL,family=gaussian())
assert_that(mean((trueBeta - estBeta)**2) < 1e-10)

# test regressions Gaussian
Y = apply(X%*%trueBeta,MARGIN = 1,FUN = function(lambda) rnorm(1,lambda))
estBeta = gmf:::glm.basic(X,Y,family=gaussian(),steps = 100,stepsize = 0.1,tol = 0)
assert_that(mean((trueBeta - estBeta)**2) / mean(trueBeta**2) < 0.05)

# test regression Poisson
Y = apply(X%*%trueBeta,MARGIN = 1,FUN = function(lambda) rpois(1,exp(lambda)))
estBeta = gmf:::glm.basic(X,Y,family=poisson(),steps = 1000,stepsize = 0.01)
assert_that(mean((trueBeta - estBeta)**2) / mean(trueBeta**2) < 0.05)

# comapre with a standard GLM
fit = glm.fit(X,Y,family=poisson())
assert_that(mean((fit$coefficients - estBeta)**2) / mean(fit$coefficients**2) < 0.05)

# test offset
offset = 5
Y = apply(X%*%trueBeta + offset,MARGIN = 1,FUN = function(lambda) rpois(1,exp(lambda)))
estBeta = gmf:::glm.basic(X,Y,family=poisson(),offset = offset,steps = 1000,stepsize = 0.01)
assert_that(mean((trueBeta - estBeta)**2) / mean(trueBeta**2) < 0.05)

# comapre with a standard GLM
fit = glm.fit(X,Y,family=poisson(),offset = 5)
assert_that(mean((fit$coefficients - estBeta)**2) / mean(fit$coefficients**2) < 0.05)

