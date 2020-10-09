slice.glm = function(X, Y, slices, coefs, offset=NULL, penalized=0, parallel=1, family = poisson(), method="step", stepsize = 1){
  res = c()
  steps = 10
  if (method == "step")
    steps = 1

  res = mclapply(1:slices,function(i){
    if (is.null(offset)){
      offset.slice = 0
    } else {
      offset.slice = offset[,i,drop=FALSE]
    }
    ok = !is.na(Y[,i])
    glm.basic(X[ok,,drop=FALSE],
              Y[ok,i,drop=FALSE],
              beta=coefs[,i],
              family=family,
              offset = offset.slice,
              stepsize = stepsize,
              penalized = penalized,
              steps=steps)
  },mc.cores = parallel)
  arr = simplify2array(res)
  if (is.vector(arr))
    arr = rbind(arr, NULL)
  arr
}