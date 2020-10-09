#library("devtools") ; install(".")
library("gllvm") ; library("gmf")
library(pROC)

cutoff = 50000

### Experimental data
load("experiments/mydata/covariates.species_new.RData")

n = min(cutoff,nrow(covariates.species))
ns = min(cutoff,ncol(covariates.species))
Y = as.matrix(covariates.species[1:n,14:ns])
X = scale(as.matrix(covariates.species[1:n,5:13]))
sum(Y)/sum(Y>-1)

rn = 0
cn = 0
while ((mean(rn) < 1 - 1e-10) && (mean(cn) < 1 - 1e-10)){
  rn = (rowMeans(Y) > 1e-3) & (apply(Y, FUN=min, MARGIN=1) >= -1e-10)
  cn = (colMeans(Y) > 1e-3) & (apply(Y, FUN=min, MARGIN=2) >= -1e-10) 
  Y = Y[rn,cn]
  X = X[rn,]
}
family = binomial()
Yfull = Y
Xfull = X

shrink.size = 0.005

shrink.dataset = function(Xfull, Yfull, shrink.size){
  yrsums = rowSums(Yfull)
  leaveout_rows = order(yrsums)[1:floor(nrow(Yfull)*(1-shrink.size))]
  Y = Yfull[-leaveout_rows,]
  X = Xfull[-leaveout_rows,]
  
  ycsums = colSums(Y)
  leaveout_cols = order(ycsums)[1:floor(ncol(Y)*(1-shrink.size))]
  Y = Y[,-leaveout_cols]
  
  ycsums = colSums(Y)
  Y = Y[,ycsums>2]
  list(X=X,Y=Y)
}
