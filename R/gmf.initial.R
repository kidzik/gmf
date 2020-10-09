qrrange = function(X,q=c(0.05,0.95)){
  range = quantile(c(X),q)
  X[X>range[2]] = range[2]
  X[X<range[1]] = range[1]
  X
}

gmf.initial = function(Y,X,d=min(dim(Y)),family = poisson()){
  cf = c()
  res = c()
  for (c in 1:ncol(Y)){
    yy = Y[,c]
    if (family$family=="binomial"){
      yy[is.na(yy)] = rbinom(sum(is.na(yy)), 1, mean(yy,na.rm=TRUE))
      yy = as.factor(yy)
    }
    
    m = glm(yy ~ X-1, family = family)
    cf = cbind(cf, m$coefficients)
    res = cbind(res, m$residuals)
  }
  
  s = svd(res, nu = d, nv = d)
  u = s$u %*% diag(sqrt(s$d[1:d]))
  v = s$v %*% diag(sqrt(s$d[1:d]))

  list(beta = cf,
       u = u,
       v = v
  )
}


