library(gmf)
library(MASS)

data(mtcars)

Y = scale(mtcars)
pca = prcomp(Y,retx = TRUE)
plot(pca$rotation[,1:2])

# PCA implementation from MASS
Ypred = pca$x[,1:2] %*% t(pca$rotation[,1:2])
mean((Ypred-Y)**2) / mean(Y**2)

# Minimal GMF
decomp = gmf(Y=Y,family=gaussian(), p=2)
mean((decomp$fit-Y)**2) / mean(Y**2)

# More problem specific and much faster
decomp = gmf(Y=Y,method="quasi",family=gaussian(),d=2,gamma = 0.1,maxIter=100,tol=1e-6,verbose=FALSE)
mean((decomp$fit-Y)**2) / mean(Y**2)

plot(norm.procrustes(pca$rotation[,1:2],decomp$v))
