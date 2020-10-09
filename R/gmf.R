#' Fit a GLLVM model, i.e. decompose a matrix of non-Gaussian responses to low-rank matrices. This is an alternative to PCA and SVD
#' in situations where observed data do not follow Gaussian distribution. This approach is particularly useful when we observe
#' responses following a distribution from the exponential family, such as
#' binomial (i.e. binary, e.g. species present or not at certain locations), nonnegative-binomial (e.g. count data
#' of genes observed in patients), Poisson, etc. In such cases PCA and similar techniques can lead to unreliable results
#' when normality assumptions are violated. If additional information for each observational unit are present,
#' these data can be incorporated into the model as a linear fixed effect. Moreover, users of this package
#' have freedom to chose the number of underlying latent factors to estimate.
#' 
#' To fit the model, the function optimizes a penalized quasi-likelihood function using the Newton method (gradient descent). Newton step
#' is performed using either a direct differntiation with diagonal approximation of Hessians (for fast computation), or using
#' Alternating Iterative Reweighted Least Squares (AIRWLS, default). AIRWLS leverages the fact that the problem for each row
#' and column can be seen as a regular GLM problem.
#' 
#' We assume that observed matrix of \eqn{n \times m}responses follow the distribution
#' \eqn{y_{ij} | \mu_{ij} \sim \mathcal{F}(\mu_{ij},\phi_j)}
#' where \eqn{g} is a link function, i.e. \eqn{g(\mu_{ij}) = \eta_{ij} = \beta_{0j} + x_i'\beta_j + u_i'\lambda_j} for a
#' set distribution \eqn{\mathcal{F}}. Quantities \eqn{\beta}, \eqn{\lambda} and \eqn{\phi} are model parameters. \eqn{\beta}
#' correspond to \eqn{d} fixed-effects or known factors. \eqn{\lambda} are  \eqn{p} random-effects or unobserved factors. \eqn{\phi}
#' are dispersion parameters. The parameter \eqn{p} is selected by the user.
#'
#' We find model parameters by optimizing
#' \deqn{L(\Psi) =-\sum_{i=1}^n \sum_{j=1}^m (y_{ij}\tilde\eta_{ij} - b(\tilde\eta_{ij})) + \frac{\lambda_U}{2} \sum_{i=1}^n u_i'u_i}
#' where \eqn{b} is the cumulant function of the given distribution in the exponential family and \eqn{\lambda_U} is the
#' weight for the L2 penalty on \eqn{U}. This is \eqn{1} and can be tuned with \code{penaltyU}, together with additional penalty
#' of the same form on V (\code{penaltyV}) or beta coefficients (\code{penaltyBeta}).
#' 
#' Each step of the Newton algorithm takes form
#' \deqn{\theta_{t+1} = \theta_{t} + s [-d^2L(\theta_t)]^{-1} \nabla L(\theta_t)}
#' where \eqn{\theta_t} is the parameter we are optimizing (\eqn{U},\eqn{V}, or \eqn{\beta}), \eqn{d^2L} stands for Hessians
#' while \eqn{\nabla L} is the gradient of L. In \code{"newton"} method these quantities are approximated directly, while
#' in \code{"newton"} an update step is transformed to a weighted linear regression problem as in GLMs.
#' 
#' We perform at most \code{maxIter} steps of the Newton algorithm with the stepsize \code{gamma}, stopping if the relative
#' change in deviance is smaller than \code{tol}. After performing optimization, the latent scores and variables are rotated
#' such that the matrix \eqn{U} has the covariance matrix equal to identity, while the matrix \eqn{V} is upper-triangular
#' with positive elements on the diagonal.
#' 
#' To initiate the Newton algorithm we either draw random variables (if \code{init} is \code{"random"}) or we regress out the
#' fixed effects and perform the SVD on the residuals. Matrices U and V derived from SVD are then used as initialization points
#' for the latent scores and variables.
#' 
#' Each of the two optimization methods takes additional parameters. In the AIRWLS method we can specify the number
#' of steps of each of the internal optimizations (i.e. IRWLS steps) of \eqn{U}, \eqn{V}, and \eqn{\beta}, before we move on
#' to optimizing the next parameter. By default \code{airwls.internal.steps=1} corresponding to only one step of the internal
#' optimization and then passing to the next iteration of the external AIRWLS loop. The parameter \code{parallel} is used
#' to distribute computations of coefficients for separate collumns and rows to different CPU cores, effectively speeding up
#' the computation. In the Newton method, damping parameter adds a diagonal regularization term to Hessians, making the 
#' inversion of Hessians more stable.
#' 
#' @title Factorize a matrix of non-Gaussian observations 
#'
#' @param Y matrix of responses
#' @param X matrix of fixed effects
#' @param family a family as in the \code{\link{glm}} interface
#' @param p number of random effects to estimate (default 2)
#' @param penaltyU penalty on latent scores
#' @param penaltyV penalty on latent variables
#' @param penaltyBeta penalty on fixed-effect coefficients
#' @param intercept should the model include the intercept
#' @param maxIter maximum number of Newton steps iterations 
#' @param gamma step size in the newton algorithm
#' @param tol threshold for the stopping criterium (a relative change in deviance) 
#' @param init initialization of model parameters: "random" (default) or "svd"
#' @param normalize_uv normalize UV to iid gaussian U and upper diagonal V with a positive diagonal
#' @param verbose output deviance every 10 iterations
#' @param method "airwls" (Alternating Iterative Reweighted Least Squares, default) or "newton" (Quasi-Newton)
#' @param damping regularization parameter of the hessian in the direct quasi Newton
#' @param airwls.internal.steps number of steps of the IRWLS in each inner loop for AIRWLS
#' @param parallel number of cores on which to run AIRWLS iterations
#' @return Returns model parameters and its deviance
#' \itemize{
#'   \item beta - fixed-effect coefficients
#'   \item U - latent scores
#'   \item V - latent vectors
#'   \item fit - fitted matrix of means
#'   \item deviance - model mean deviance
#' }
#' @references Lukasz Kidzinski, Francis K.C. Hui, David I. Warton, Trevor J. Hastie
#' \emph{Generalized Matrix Factorization}
#' arXiv, 2020
#' @import assertthat parallel
#' @importFrom whitening whiteningMatrix
#' @examples
#' # In this example we simulate synthetic matrix of resoponses from the Binomial
#' # family and we fit the GLLVM model 
#' data = gmf.simulate(family=binomial())
#' 
#' # Fit a GLLVM model using Newton method
#' model = gmf(data$Y)
#' 
#' # Plot scores
#' plot(model$u)
#' 
#' # Compute the mean deviance
#' print(matrix.deviance(model$fit, data$Y, model$family))
#' @export
gmf = function(Y,
               X=NULL,
               family=poisson(),
               p=2,
               penaltyU=1,
               penaltyV=0,
               penaltyBeta=0,
               intercept=FALSE,
               maxIter=NULL,
               gamma=NULL,
               tol=1e-3,
               damping=0,
               init="random",
               normalize_uv=TRUE,
               verbose=TRUE,
               method="airwls",
               airwls.internal.steps = 1,
               parallel=1){
  # Adjust parameters to defaults for given methods
  # AIRWLS requires less steps and takes longer steps, but is slower
  # quasi-Newton is fast but since gradients are less accurate it may need shorters steps
  if (is.null(maxIter)){
    if (method == "quasi"){
      maxIter = 1000
    }
    if (method == "airwls"){
      maxIter = 100
    }
  }
  if (is.null(gamma)){
    if (method == "quasi"){
      gamma = 0.01
    }
    if (method == "airwls"){
      gamma = 0.1
    }
  }
  
  # Derive dimensions
  n = dim(Y)[1]
  m = dim(Y)[2]
  
  # Add the intercept as a column of ones in X
  if (intercept)
    X = cbind(matrix(1,n,1),X)
  
  # X still may be NULL if no intercept and empty input X. Then d = 0
  d = 0
  if (!is.null(X))
    d = dim(X)[2]
  
  # Initialize U, V and beta using the selected method
  if (init == "svd"){
    initialization = gmf.initial(Y,X,p,family)
    U = initialization$u #[,1:p,drop=FALSE]
    V = initialization$v #[,1:p,drop=FALSE]
    beta = rbind(initialization$beta)
  }
  if (init=="random"){
    sd = 1e-1
    udim = c(n,p)
    vdim = c(m,p)
    betadim = c(d,m)
    U = array(rnorm(prod(udim))/prod(dim(udim))*sd, udim)
    V = array(rnorm(prod(vdim))/prod(dim(vdim))*sd, vdim)
    beta = array(rnorm(prod(betadim))/prod(betadim)*sd, betadim)
  }
  
  # Remember the last deviance
  llast = Inf
  
  # remember where is missing data so that we can keep replacing it
  # with more and more accurate models
  isna = is.na(Y)
  Y[isna] = mean(Y,na.rm=TRUE)
  
  # Since regression and latent predictions are used multiple times
  # throughout the computation, we will store them
  regpart = 0
  if (d)
    regpart = X %*% beta
  latentpart = U %*% t(V)
  
  for (i in 1:maxIter){
    # After the first iteration, replace NAs with model values
    if (i > 1){
      Y[isna] = family$linkinv(latentpart + regpart)[isna]
    }
    
    if (method == "airwls"){
      # Perform airwls.internal.steps internal steps of
      # the regularized IRWLS
      for (j in 1:airwls.internal.steps){
        # Get column coefficients
        coefs = slice.glm(cbind(X,U), Y, m, coefs = rbind(beta,t(V)),
                          offset = NULL, penalized=penaltyV, parallel = parallel,
                          family = family, method = method, stepsize = gamma)
        if (d)
          beta = coefs[1:d,,drop=FALSE]
        V = t(coefs[(d+1):(d+p),,drop=FALSE])
      }
      for (j in 1:airwls.internal.steps){
        # Get row coefficients, correct for regression using the offset
        offset = NULL
        if (d)
          offset = t(regpart)
        U = t(slice.glm(V, t(Y), n, coefs = t(U), offset=offset,
                        penalized=penaltyU, parallel = parallel, family = family,
                        method = method, stepsize = gamma))
      }
    }
    else {
      # Directly compute gradients and Hessians
      eta = latentpart + regpart
      M = family$linkinv(eta)
      
      # Set up helper matrices for computing differentials 
      dratio = family$mu.eta(eta) / family$variance(M)
      ddiff = ((Y - M) * dratio)
      ddratio = dratio * family$mu.eta(eta)
      
      # gaussion()$mu.eta returns a vector instead of a matrix
      # potentially, other families too, so we add the following umbrella conversion
      if (!is.matrix(ddratio)){ 
        ddratio = matrix(ddratio,nrow(eta),ncol(eta)) 
      }
      
      # Update model parameters, updating U and V at once
      newU = update.params(U,V,penaltyU,ddiff,ddratio,gamma,damping)
      V = update.params(V,U,penaltyV,t(ddiff),t(ddratio),gamma,damping)
      U = newU
      
      # Update beta if exists
      if (d){
        beta = t(update.params(t(beta),X,penaltyBeta,t(ddiff),t(ddratio),gamma,damping))
      }
    }
    
    # Update linear predictors
    if (d)
      regpart = X %*% beta
    latentpart = U %*% t(V)
    
    # Get predicted means
    predicted = family$linkinv(latentpart + regpart)
    
    # Compute penalties and mean deviance
    penalties = norm(U,"F")*penaltyU + norm(V,"F")*penaltyV + norm(beta,"F")*penaltyBeta
    l = (matrix.deviance(predicted, Y, family) + penalties)/prod(dim(Y))
    
    # Diagnostics output: iteration -> deviance    
    if (i %% 10 == 0 && verbose){
      cat("iteration ",i,", PQL = ",l,"\n",sep = "")
    }
    
    # Check the stopping criteria
    if (abs(llast - l) / l < tol){
      break
    }
    else{
      llast = l
    }
  }
  
  # Rotate U and V so that cov(U) is the identity matrix, and V is upper diagonal
  if (normalize_uv){
    decomp = correct.uv(U, V)
    U = decomp$u
    V = decomp$v
    #    assert_that(norm(latentpart - U%*%t(V) ) < 1e-2)
  }
  
  print(paste("Stopped after",i,"iterations with deviance",l))
  
  list(beta = beta,
       u = U,
       v = V,
       fit = family$linkinv(latentpart + regpart),
       family = family,
       deviance = l)
}

# Rotates U and V such that the covariance of U is diagonal, and V is upper-triangular with a positive diagonal
#' @importFrom whitening whiteningMatrix
correct.uv = function(U, V){
  S = cov(U)

  if (ncol(U)==1){
    return(list(u=U/sqrt(c(S)),v=V*sqrt(c(S))))
  }

  # Make cov of U identity
  W = whiteningMatrix(S)
  U = U %*% W
  V = V %*% t(solve(W))

  # Make V lower triangular
  V.qr = qr(t(V))
  U = U %*% qr.Q(V.qr)
  V = t(qr.R(V.qr))
  
  # Positive diagonal of V
  d = diag(V)
  V = t(sign(d)*t(V))
  U = t(sign(d)*t(U))
  
  list(u=U,v=V)
}

# compute gradients for the parameters. For stability, in estimation remove values
# that are too large
update.params = function(U,V,penalty,ddiff,ddratio,gamma,damping){
#  stop()
  th = 1e2
  
  # filter out rows with at least one large value of V 
  keep.rows = rowSums(abs(V)>th) < 0.5

  nfull = sum(keep.rows)
  
  # correct rows with large ddiff (remove those and reweight rows)
  corrections = nfull / rowSums(abs(ddiff)<th)
  ddiff[abs(ddiff)>th] = 0
  dU = - sweep(ddiff[,keep.rows] %*% V[keep.rows,], MARGIN=1, corrections, `*`) + U*penalty

  # correct rows with large ddratio (remove those and reweight rows)
  corrections = nfull / rowSums(abs(ddratio)<th)
  ddratio[abs(ddratio)>th] = 0
  ddU = sweep(ddratio[,keep.rows] %*% V[keep.rows,]**2, MARGIN=1, corrections, `*`) + penalty + damping

  gradU = dU / ddU
  U = matrix(U - gamma * gradU,nrow(U),ncol(U)) 
}