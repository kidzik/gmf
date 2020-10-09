library("gllvm")
library("gmf")

n = 100
m = 100
d = 3
tol = 1e-3

family = poisson()

sim.data = gmf.simulate(n = n, m = m, d = d, family=family)
Y = sim.data$Y

npoints = length(c(Y))
mask = matrix(FALSE,nrow(Y),ncol(Y))
mask[sample(npoints)[1:floor(npoints*0.1)]] = TRUE
Ytest = Y[mask]

if (family$family=="poisson")
  Y[mask] = rpois(sum(mask), family$linkinv(mean(sim.data$M)))
if (family$family=="binomial")
  Y[mask] = rbinom(sum(mask), 1, family$linkinv(mean(sim.data$M)))
X = sim.data$X

ptm <- proc.time()
model.newton = gmf(Y = Y, X = sim.data$X, d = sim.data$d, gamma= 5e-1, maxIter = 1000,
                          family=family, method = "quasi", tol = tol)
time.newton = proc.time() - ptm

ptm <- proc.time()
model.gmf.airwls = gmf(Y = Y, X = sim.data$X, d = sim.data$d, gamma= 1e-1, maxIter = 1000,
                                family=family, method = "airwls", tol = tol)
time.gmf.airwls = proc.time() - ptm

ptm <- proc.time()
model.gllvm = gllvm(y = Y, X = sim.data$X,
                    formula = ~ ., family=family,
                    num.lv = sim.data$d,
                    reltol = tol, trace = TRUE)
time.gllvm = proc.time() - ptm

print(time.gllvm)
print(time.gmf.airwls)
print(time.newton)

M.gllvm = family$linkinv(residuals(model.gllvm)$linpred)
M.gmf.airwls = model.gmf.airwls$fit
M.gmf.newton = model.newton$fit

# goodness of fit
dev.null = matrix.deviance(mean(Y), Y, family)
1 - matrix.deviance(M.gllvm, Y, family) / dev.null
1 - matrix.deviance(M.gmf.airwls, Y, family) / dev.null
1 - matrix.deviance(M.gmf.newton, Y, family) / dev.null
