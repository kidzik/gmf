#' Procrustes error measures the difference between two linear spaces spanned by a set of vectors.
#' In this function, columns of A and B are interpreted as spanning vectors. The smaller the Procrustes
#' norm the more similar are the spaces. In the context of GMF this is used to compare spaces spanned
#' by the latent variables.
#' 
#' We follow Niku, et al. (2019) and for the purpose of this work we define Procrustes as
#' \deqn{P(A, B) = \min_{R'R = I}\|A - R B\|,}
#' where \eqn{\|\cdot\|} is the Frobenius norm, \code{A} and \code{B} are matrices to compare
#' and \eqn{R} is a rotation matrix. Here, we first scale \code{A} and \code{B} such that \eqn{\|A\|=\|B\|=1.}
#' 
#' @title Compute Procurstes error between normalized A and B
#'
#' @references Niku, Jenni, et al. "Efficient estimation of generalized linear latent variable models." PloS one 14.5 (2019): e0216129.
#' @param A first matrix 
#' @param B second matrix
#' @return Returns a result of the same form as \code{\link[vegan]{procrustes}}.
#' @export
#' @import vegan
#' @examples
#' # In this example we compare estimated latent space with the true
#' # latent space which is known in simulations
#' data = gmf.simulate(family=binomial())
#' 
#' # Fit a GLLVM model using Newton method
#' model = gmf(data$Y)
#' 
#' # Compute the procrustes error
#' pnorm = norm.procrustes(data$U, model$u)
#' 
norm.procrustes = function(A,B){
  UA = A / norm(A, type = "F")
  UB = B / norm(B, type = "F")
  procrustes(UA, UB)
}
