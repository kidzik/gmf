#library("roxygen2") ; roxygenize()
#library("devtools") ; install(".")
library("gmf")
library(raster) 

n = 100
m = 100
d = 2
reps = 1
family_name = "poisson"

args = commandArgs(trailingOnly=TRUE)
if (length(args) > 1){
   n = as.numeric(args[1])
   m = as.numeric(args[2])
   d = as.numeric(args[3])
   family_name = args[4]
   reps = 50
}
print(args)

if (family_name == "poisson"){
   family = poisson()
}
if (family_name == "binomial"){
   family = binomial()
}

options(repr.plot.width=4, repr.plot.height=4) # set the size of the plots
theme_set(theme_minimal(base_size = 15)) # choose a theme and large font for presentation clarity

data(antTraits)
Y = as.matrix(antTraits$abund)
X = as.matrix(antTraits$env)
X = X[,-(3:5),drop=FALSE]
X = scale(X,center = TRUE, scale = TRUE)
names(antTraits$env)

#model.gllvm = gllvm(y = Y, X = X, formula = ~ . - 1,  family=family_name, num.lv = d, beta0com =FALSE)
model.gmf.newton = gmf(Y = Y, X = X, d = d, gamma= 0.1, maxIter = 1000, family=poisson(), penaltyU = 1, penaltyV = 0, init = "random")

cov.X = cov(X)
cov.U = cov(model.gmf.newton$u)
cov.V = cov(model.gmf.newton$v)

matrix.sqrt = function(M){
  S = svd(M)
  S$u %*% sqrt(diag(S$d))
}

simulate = function(n, m, cov.X, cov.U, cov.V, d = ncol(cov.U), family = binomial(), sd=1, sdnoise=0.3){
  p = ncol(cov.X)
  beta = matrix(rnorm(m*p,sd=sd), m, p)*10
  
  X = matrix(rnorm(p*n),n,p) %*% matrix.sqrt(cov.X)*3
  U = matrix(rnorm(d*n),n,d) #%*% matrix.sqrt(cov.U) * sd
  V = matrix(rnorm(d*m),m,d) %*% matrix.sqrt(cov.V)*10
  
  linpred = (U%*%t(V) + X%*%t(beta))*sd
  dims = dim(linpred)
  noise = matrix(rnorm(prod(dims))*sdnoise,dims[1],dims[2])
  linnoise = linpred + noise
  th = 10
  
  M = family$linkinv(linnoise)
  M[M>th] = log(M[M>th])
  
  if (family$family == "poisson")
    rr = rpois
  if (family$family == "binomial")
    rr = function(n, p) rbinom(n, 1, p)
  if (family$family == "Gamma")
    rr = function(n, p) rgamma(n, 1, scale = p)
  if (family$family == "gaussian")
    rr = rnorm
  if (family$family == "Tweedie")
    rr = function(n, p) rTweedie(p, p = var.power)
  if (substr(family$family,1,17) == "Negative Binomial"){
    theta = as.numeric(substr(family$family,19,nchar(family$family)-1))
    rr = function(n, p) { rnegbin(n, 1, p) }
  }
 
  list( X = X,
        U = U,
        V = V,
        M = M,
        noise = noise,
        d = d,
        beta = beta,
        linpred = linpred,
        family = family,
        Y = apply(M, 1:2, function(x){
          rr(1, x)
        }) 
        )
}

gamma = 0.1
for (i in 1:reps){
    tol = 1e-3
    simulation = simulate(n,m,cov.X,cov.U,cov.V,d=d,family=family,sd=0.1,sdnoise=0)
    
    # 
    ptm <- proc.time()
    model.newton = gmf(Y = simulation$Y, X = simulation$X, d = simulation$d,
                              gamma= gamma, maxIter = 1000, family=simulation$family,
                              penaltyU = 1, penaltyV = 0, penaltyBeta = 0, tol=tol*1e-2, #damping = norm(simulation$Y,"M"),
                              normalize_uv = TRUE, method="quasi")
    time.newton <- proc.time() - ptm
    
    ptm <- proc.time()
    model.gllvm = gllvm(y = simulation$Y, X = simulation$X, formula = ~ . - 1,
                        family=simulation$family, num.lv = simulation$d,
                        beta0com = FALSE, reltol = tol)
    time.gllvm <- proc.time() - ptm

    ptm <- proc.time()
    model.airwls = gmf(Y = simulation$Y, X = simulation$X, d = simulation$d,
                              gamma= gamma, maxIter = 100, family=simulation$family, tol=tol,
                              penaltyU = 1, penaltyV = 0, penaltyBeta = 0, normalize_uv = TRUE, 
                              method = "airwls")
    time.airwls <- proc.time() - ptm

     save(simulation,
          model.gllvm, model.newton, model.airwls, 
          time.gllvm, time.newton, time.airwls,
          file = paste0("output/",paste(i,n,m,d,simulation$family$family,format(Sys.time(), "%s"),sep="-"),".RData") )
}

## Present results of the last simulation
M.gllvm = simulation$family$linkinv(residuals(model.gllvm)$linpred)
M.gmf.newton = model.newton$fit
M.gmf.airwls = model.airwls$fit

# Mean deviance
matrix.deviance(mean(simulation$Y), simulation$Y, simulation$family)
matrix.deviance(M.gllvm, simulation$Y, simulation$family)
matrix.deviance(M.gmf.newton, simulation$Y, simulation$family)
matrix.deviance(M.gmf.airwls, simulation$Y, simulation$family)

# MSE of fixed effects
mean((simulation$beta)**2)
mean((simulation$beta - model.gllvm$params$Xcoef)**2)
mean((simulation$beta - t(model.newton$beta))**2)
mean((simulation$beta - t(model.airwls$beta))**2)

# Procrustes error
ss = sum((simulation$V)**2)
norm.procrustes(simulation$V, model.gllvm$params$theta)$ss / ss
norm.procrustes(simulation$V, model.newton$v)$ss / ss
norm.procrustes(simulation$V, model.airwls$v)$ss / ss

# Rasters
cuts = c(0,5,10,15,20,25,30)
title("Ants")
plot(raster(as.matrix(simulation$Y)),lab.breaks=cuts)
plot(raster(as.matrix(M.gllvm)),lab.breaks=cuts)
plot(raster(as.matrix(M.gmf.newton)),lab.breaks=cuts)
plot(raster(as.matrix(M.gmf.airwls)),lab.breaks=cuts)
title("")

# Process all simulations using "simulation.ants.combine.R"