norm_vec <- function(x) sqrt(sum(x^2))

glm.step = function(X,Y,beta,family,offset = 0, stepsize = 1, penalized = 0)
{
  eta = offset + X %*% beta
  mu = family$linkinv(eta) 

  mu.eta = family$mu.eta( eta )
  var.mu = family$variance(mu)
  
  Winv = var.mu / mu.eta**2
  W = mu.eta**2 / var.mu
  
  Z = (eta - offset) + (Y - mu) * Winv

  # Use only rows that give reasonable values  
  thresh = 1e20
  keep.rows = !is.nan(c(W)) & !is.infinite(c(W)) & !is.nan(c(Z)) & !is.infinite(c(Z)) & (c(W)>1/thresh) & (c(W)<thresh) & (c(abs(Z))>1/thresh) & (c(abs(Z))<thresh)
  
  if (sum(keep.rows) < ncol(X)+1)
  {
    stop("Too many rows with infinite values")
  }
  Xfull = X[keep.rows,,drop=FALSE]
  Yfull = Z[keep.rows,,drop=FALSE]
  Wfull = c(W)[keep.rows]
  
  if (penalized){
    Yfull = rbind(Yfull,matrix(0, ncol(Xfull), 1))
    Xfull = rbind(Xfull,diag(1, ncol(Xfull)))
    Wfull = c(Wfull,rep(penalized, ncol(Xfull)))
  }
  
  fit = suppressWarnings(lsfit(Xfull,Yfull,Wfull,intercept = FALSE))$coefficients
  fit*stepsize + (1-stepsize)*beta
}

glm.basic = function(X,Y,beta=NULL,family=gaussian(),tol=1e-5,offset=0, stepsize = 1, penalized = 0, steps=1){
  if (is.null(beta))
    beta = matrix(0,ncol(X),1)
  for (i in 1:steps){
    betaold = beta
    beta = glm.step(X,Y,beta,family,offset=offset,stepsize=stepsize,penalized=penalized)
    if (norm_vec(betaold - beta)/norm_vec(beta) < tol)
      break
  }
  beta
}

