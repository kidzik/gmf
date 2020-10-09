#' Compute a deviance on all entries of a matrix, using \code{dev.resids} function
#' from the provided family.
#' 
#' We define deviance as
#' \deqn{D(Y,M)=2\sum_{i=1}^n\sum_{j=1}^m \log (p(Y_{ij} | {M}_{ij}))-\log (p(Y_{ij}| {M }_{0,ij}),}
#' where \eqn{M} are the parameters predicted from the model, \eqn{{M}_{0}}
#' are the linear predictors from the saturated model, and \eqn{p} is the
#' density function of the selected \code{family}.
#' 
#' Instead of using deviance as a stand-alone metric of the fit, it recommended to use it
#' for assessing the relative difference between two models or a fitted model and the null model.
#' Ratio of the fitted model and the null model can be interpreted as the unexplained variance.
#' 
#' @title Compute a deviance between predicted and observed matrices
#'
#' @param pred predicted matrix
#' @param obs observed matrix
#' @param family distribution as in the \code{\link{glm}} interface
#' 
#' @return Returns the mean deviance omitting NA values
#' @export
#' @examples
#' # In this example we compare a fitted model with a null model
#' data = gmf.simulate(family=binomial())
#' 
#' # Fit a GLLVM model using Newton method
#' model = gmf(data$Y)
#' 
#' # Compute mean deviance of the null model and the fitted model
#' dev.null = matrix.deviance(mean(data$Y), data$Y, model$family)
#' dev.model = matrix.deviance(model$fit, data$Y, model$family)
#' 
#' cat("Deviance explained by the model: ", 1-dev.model/dev.null)
matrix.deviance = function(pred, obs, family){
  if (length(pred) == 1){
    pred.matrix = obs
    pred.matrix[] = pred
    pred = pred.matrix
  }
  isna = is.na(obs)
  pred = pred[!isna]
  obs = obs[!isna]
  
  mean(family$dev.resids(obs, pred, 1),na.rm = TRUE)
}

