#' Generate model parameters of a GLLVM and draw a matrix from the model. Return the matrix, model parameters, and underlying
#' observed and unobserved factors. The drawn matrix is of size \eqn{n \times m}.
#' 
#' A matrix \eqn{[y_{ij}]} for \eqn{1 \leq i \leq n} and \eqn{1 \leq j \leq m} is drawn from
#' \eqn{y_{ij} | \mu_{ij} \sim \mathcal{F}(\mu_{ij},\phi_j)}
#' where \eqn{g} is a link function, i.e. \eqn{g(\mu_{ij}) = \eta_{ij} = \beta_{0j} + x_i'\beta_j + u_i'\lambda_j} for a
#' set distribution \eqn{\mathcal{F}}. Quantities \eqn{\beta}, \eqn{\lambda} and \eqn{\phi} are model parameters. \eqn{\beta}
#' correspond to fixed-effects or known factors. \eqn{\lambda} are random-effects or unobserved factors. \eqn{\phi} are dispersion parameters.
#'  
#' Model parameters as well as observed and unobserved parameters are simulated from a Gaussian distribution. These are tuned such that
#' \eqn{\eta_{ij}} are approximately on the scale \eqn{[-3,3]} enabling testing a wide variety of distributions and link functions.
#' 
#' To generate a matrix we first generate model parameters. Next, we generate observed \eqn{x_i} and latent scores \eqn{u_i}.
#' Together with the link function that identifies mean values \eqn{\mu_{ij}} which allow us to draw \eqn{y_{ij}}
#' from the defined distribution.
#' 
#' @title Simulate a GLLVM model and draw a matrix from the model
#'
#' @param n number of rows (observational units, subjects, etc.)
#' @param m number of columns (species, genes, etc.)
#' @param d number of fixed effects (observed)
#' @param p number of random effects (unobserved)
#' @param family a family as in the \code{\link{glm}} interface
#' @return Returns a list with model parameters and simulated values
#' \itemize{
#'   \item Y - simulated matrix of responses 
#'   \item X - simulated fixed effects
#'   \item M - matrix of means \eqn{\mu_{ij}}
#'   \item U - matrix of latent scores
#'   \item V - matrix of latent variables
#'   \item beta - matrix of coefficients of fixed effects
#' }
#' @export
#' @import stats
#' @import MASS
gmf.simulate = function (n = 30, m = 60, d = 1, p = 2, family = poisson()){
  V = matrix(rnorm(m*p,0,1),p,m) # site variables

  U = matrix(rnorm(n*p,0,1),p,n) # species variables
  S = matrix(rnorm(n*n),n,n)
  Sigma = S %*% t(S)
  Sigma = 0.8 * Sigma / norm(Sigma)
  U = matrix(mvrnorm(n = p, mu = rep(0,n), Sigma = Sigma),p,n)
  
  trueBeta = matrix(runif(m*d,-1,1),d,m) # species variables
  X = matrix(rnorm(d*n),n,d)
  
  regX = X %*% trueBeta
  
  latPart = t(U) %*% V #+ rnorm(m*N,0,1)
  M = latPart + regX
  
  if (family$family == "poisson")
    rr = rpois
  if (family$family == "binomial")
    rr = function(n, p) rbinom(n, 1, p)
  if (family$family == "Gamma")
    rr = function(n, p) rgamma(n, 1, scale = p)
  if (family$family == "gaussian")
    rr = rnorm
  if (substr(family$family,1,17) == "Negative Binomial"){
    theta = as.numeric(substr(family$family,19,nchar(family$family)-1))
    rr = function(n, p) { rnegbin(n, p, theta) }
  }

  Y = apply(latPart + regX,1:2,function(x){ rr(1, family$linkinv(x) ) } )
  list(Y=Y,X=X,M=M,latPart = latPart, regPart = regX, d=d,p=p,m=m,n=n,U=t(U),V=t(V),beta=trueBeta)
}