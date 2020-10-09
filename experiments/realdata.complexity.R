source("experiments/load.data.R")

###################################
# Experiment 1: Compute time vs n,m
###################################
time.gllvm = c()
sizes = 0.005 * 1:2 # In paper 1:12
time.gmf.newton = c()
time.gmf.airwls = c()
time.gmf.airwls.glmnet = c()

dev.null = c()
dev.gllvm = c()
dev.gmf.newton = c()
dev.gmf.airwls = c()
dev.gmf.airwls.glmnet = c()

for (shrink.size in sizes){
  data = shrink.dataset(Xfull,Yfull,shrink.size)
  
  X = data$X
  Y = data$Y

  d = 3
  
  ptm <- proc.time()
  model.gllvm = gllvm(y = Y,
                      X = X,
                      formula = ~ .,
                      family=binomial(),
                      num.lv = d,
                      reltol = 1e-3,
                      trace = TRUE,
                      method="VA")
  tm = proc.time() - ptm
  time.gllvm = c(time.gllvm, tm[3])
  
  ptm <- proc.time()
  model.gmf.newton = gmf.newton(Y = Y,
                                X = X,
                                d = d,
                                family = binomial(),
                                maxIter=5000,
                                gamma=0.1,
                                penaltyU = 1,
                                penaltyV = 0,
                                penaltyBeta = 0,
                                method = "quasi")
  tm = proc.time() - ptm
  time.gmf.newton = c(time.gmf.newton, tm[3])
  
  ptm <- proc.time()
  model.gmf.airwls = gmf.newton(Y = Y,
                            X = X,
                            d = d,
                            gamma = 0.01,
                            maxIter = 10000,
                            family = family,
                            penaltyU = 1,
                            penaltyV = 0,
                            method = "airwls")
  time.gmf.airwls <- c(time.gmf.airwls, (proc.time() - ptm)[3])
  
  family = binomial()
  M.gllvm = family$linkinv(residuals(model.gllvm)$linpred)
  M.gmf.newton = model.gmf.newton$fit
  M.gmf.airwls = model.gmf.airwls$fit

  dev.null = c(dev.null, matrix.deviance(mean(Y), Y, family))

  dev.gllvm = c(dev.gllvm, matrix.deviance(M.gllvm, Y, family))
  dev.gmf.newton = c(dev.gmf.newton, matrix.deviance(M.gmf.newton, Y, family))
  dev.gmf.airwls = c(dev.gmf.airwls, matrix.deviance(M.gmf.airwls, Y, family))

  df = data.frame(dev.null, dev.gllvm, dev.gmf.newton, dev.gmf.airwls, #dev.gmf.airwls.glmnet,
       time.gllvm, time.gmf.newton, time.gmf.airwls #, time.gmf.airwls.glmnet
       )
  save(model.gllvm, model.gmf.newton, model.gmf.airwls, Y, df, file=paste0("output/iteration-", shrink.size, ".Rdata"))
}

# Run postprocess-realdata.R after experiments are done