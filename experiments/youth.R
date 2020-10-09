library(gmf)
library(psychotools)
library(vegan)
library(gllvm)
library(raster)
library(psychotools)

data("YouthGratitude")
X = YouthGratitude[,2,drop=FALSE]
Y = YouthGratitude[,-(1:3)]
X = as.matrix(X[1:100,,drop=FALSE])
Y = as.matrix(Y[1:100,,drop=FALSE])
p = 2

fm = gaussian()

options(warn = -1) 
ptm <- proc.time()
model.gllvm = gllvm(y = Y, X = X, formula = ~ ., family="gaussian", num.lv = p)
time.gllvm = proc.time() - ptm

ptm <- proc.time()
model.gmf.newton = gmf(Y = Y, X = X, p = p, gamma = 0.05, maxIter = 1000, tol = 1e-5, method = "quasi", family = fm)
time.gmf.newton = proc.time() - ptm

ptm <- proc.time()
model.gmf.airwls = gmf(Y = Y, X = X, p = p, gamma = 0.1, maxIter = 100, method = "airwls", parallel = 1, family = fm)
time.gmf.airwls = proc.time() - ptm

M.gllvm = fm$linkinv(residuals(model.gllvm)$linpred)
M.gmf.newton = model.gmf.newton$fit
M.gmf.airwls = model.gmf.airwls$fit

options(repr.plot.width=4, repr.plot.height=4)
cuts=c(0,5,10,15,20,25,30)
plot(raster(Y),lab.breaks=cuts,zlim=c(0,max(cuts)))
title("Abundance data")
plot(raster(M.gllvm),lab.breaks=cuts,zlim=c(0,max(cuts)))
title("gllvm")
plot(raster(M.gmf.newton),lab.breaks=cuts,zlim=c(0,max(cuts)))
title("Newton")
plot(raster(M.gmf.airwls),lab.breaks=cuts,zlim=c(0,max(cuts)))
title("Alternating GLMM")

denom = matrix.deviance(mean(Y), Y, family = fm)
1 - matrix.deviance(M.gllvm, Y, family = fm) / denom
1 - matrix.deviance(M.gmf.newton, Y, family = fm) / denom
1 - matrix.deviance(M.gmf.airwls, Y, family = fm) / denom