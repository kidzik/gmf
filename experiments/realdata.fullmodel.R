source("experiments/load.data.R")

#####################################
# Experiment 2: Missing at random AUC
#####################################
load("old/models/model.real.big.newton-5.Rda")
load("old/model.real.big.newton-20-leaveout-rows-samples-1599440413.1167.Rda")

# Parameters:
nsample = 500
train.sample = 0.95
d = 20

leaveout_rows = sample(nrow(X))[1:floor(nrow(X)*(1-train.sample))]
ones = sample(which(model.gmf.nuc$fit[-leaveout_rows,] > 0),nsample)
zeros = sample(which(model.gmf.nuc$fit[-leaveout_rows,] < 0),nsample)
leaveout_indexes = c(ones, zeros)

Yleftout = Y[-leaveout_rows,]
Yleftout[leaveout_indexes] = NA
Y[-leaveout_rows,][leaveout_indexes]
mean(Y[-leaveout_rows,][leaveout_indexes])
Xint=cbind(1,X)

model.gmf.newton = gmf(Y = Yleftout, X = Xint[-leaveout_rows,],
                              d = d, family = binomial(),maxIter=500, gamma=1e-2,
                              penaltyU = 1, penaltyV = 0,
                              penaltyBeta = 0, method = "quasi")
save(model.gmf.newton, leaveout_rows, ones, zeros, file=paste0("model.real.big.newton-",d,"-leaveout-rows-samples-",as.numeric(Sys.time()),".Rda"))

model.gmf.airwls = gmf(Y = Yleftout, X = X[-leaveout_rows,],
                              d = d, family = family,
                              gamma = 0.01, maxIter = 500,
                              penaltyU = 1, penaltyV = 0, method = "airwls")
save(model.gmf.newton, leaveout_rows, ones, zeros, file=paste0("model.real.big.newton-",d,"-leaveout-rows-samples-",as.numeric(Sys.time()),".Rda"))

model.gmf.newton = model.gmf.nuc

regpart = Xint[-leaveout_rows,] %*% model.gmf.newton$beta
S = svd(t(model.gmf.newton$v))
df = data.frame(factor.number=1:length(S$d),singular.value=S$d)
g = ggplot(df,aes(x=factor.number,y=singular.value)) +
  geom_point() + ggtitle("Scree plot of latent factors") + 
  xlab("Ordered factor number") + ylab("Singular value") +
  ylim(0,100)
g
ggsave("figures/scree.pdf",g,width = 3.5,height=5)

dd=3 
latentpart = model.gmf.newton$u %*% S$u[,1:dd] %*% diag(S$d[1:dd]) %*% t(S$v[,1:dd])
model.gmf.newton$fit = latentpart + regpart


#save(model.gmf.nuc, leaveout_rows, leaveout_indexes, file="model.real.big.newton-5-leaveout-rows-samples.Rda")
#load("model.real.big.newton-1.Rda")
time.gmf.nuc = proc.time() - ptm
time.gmf.nuc

#init = gmf.initial(Y[-leaveout,],X[-leaveout,],fm)

f = function(x) {matrix.deviance(family$linkinv(model.gmf.newton$fit*x)[leaveout_indexes], Y[-leaveout_rows,][leaveout_indexes], family)}
g = function(x) {matrix.deviance(family$linkinv(X[-leaveout_rows,] %*% model.gmf.newton$beta*x)[leaveout_indexes], Y[-leaveout_rows,][leaveout_indexes], family)}
h = function(x) {matrix.deviance(model.gmf.newton$fit*x, Y[-leaveout_rows,], family)}
args = 1:10/10
plot(args, simplify2array(lapply(args, f)),t='l',ylim=c(0,2))
lines(args, simplify2array(lapply(args, g)),col='red')
lines(args, simplify2array(lapply(args, h)),col='blue')

matrix.deviance(family$linkinv(Xint[-leaveout_rows,] %*% model.gmf.newton$beta)[leaveout_indexes], Y[-leaveout_rows,][leaveout_indexes], family)  

xx = Xint[-leaveout_rows,]
#xx = cbind(1,xx)
#beta = rbind(model.gmf.newton$beta0, model.gmf.newton$beta)
beta = model.gmf.newton$beta

calibration = 1
dev.null = matrix.deviance(mean(Y[-leaveout_indexes]), Y[-leaveout_rows,][leaveout_indexes], family)
1 - matrix.deviance(calibration*model.gmf.newton$fit[leaveout_indexes], Y[-leaveout_rows,][leaveout_indexes], family) / dev.null
1 - matrix.deviance(family$linkinv(calibration*model.gmf.newton$u %*% (model.gmf.newton$v))[leaveout_indexes], Y[-leaveout_rows,][leaveout_indexes], family) / dev.null
1 - matrix.deviance(family$linkinv(calibration*xx %*% beta)[leaveout_indexes], Y[-leaveout_rows,][leaveout_indexes], family) / dev.null

library(ggplot2)

roc.full = roc(Y[-leaveout_rows,][leaveout_indexes],family$linkinv(model.gmf.newton$fit)[leaveout_indexes])
roc.fixed = roc(Y[-leaveout_rows,][leaveout_indexes],family$linkinv(xx %*% beta)[leaveout_indexes])

g <- ggroc(list('Full Model'=roc.full, 'Fixed effects'=roc.fixed), size = 1.25) +
  ggtitle("ROC of fitted models") + labs(color = "Model")
g
ggsave("figures/roc-d10.pdf",g,width = 6,height=5)


  ## ANALYSIS IN-SAMPLE
reg.part.0 = t(t(rep(1,nrow(Y[leaveout_rows,])))) %*% model.gmf.nuc$beta0
reg.part = X[leaveout_rows,] %*% model.gmf.nuc$beta
#lat.part = model.gmf.nuc$u %*% model.gmf.nuc$v

matrix.deviance(fm$linkinv(reg.part.0 + reg.part),Y[leaveout,],fm)

df = data.frame(x = c(Y[rows,]), y = c(ypred))
boxplot(fm$linkinv(df$y) ~ df$x)

## DEVIANCE OUT-SAMPLE
reg.part.0 = t(t(rep(1,nrow(Y[leaveout,])))) %*% model.gmf.nuc$beta0
reg.part = X[leaveout,] %*% model.gmf.nuc$beta
M.gmf.nuc = fm$linkinv(reg.part + reg.part.0)

matrix.deviance(M.gmf.nuc, Y[leaveout,], fm)

## DEVIANCE OUT-SAMPLE RANDOM POINTS
ones = sample(which(Y==1),1000)
zeros = sample(which(Y==0),1000)
Yleftout = Y[]
Y[c(ones,zeros)] = NA
topredict = Yleftout[is.na(Y)]
preds = fm$linkinv(model.gmf.newton$fit[is.na(Y)])

matrix.deviance(preds, topredict, fm)
roc.curve = roc(c(topredict),c(preds))
plot(roc.curve)

png("roc-d1-random.png",width=3.25,height=3.25,units="in",res=1200)
plot(roc.curve)
dev.off()

## PLOT ROC FOR OUT-SAMPLE
roc.curve = roc(c(Y[leaveout,]),c(M.gmf.nuc))

png("roc-d1-5.png",width=3.25,height=3.25,units="in",res=1200)
plot(roc.curve)
dev.off()

auc(roc.curve)
