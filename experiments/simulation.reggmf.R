library("ggplot2")
library("gllvm")
library("gmf")

scaleFUN <- function(x) sprintf("%.1f", x)

thm = theme_minimal(base_size = 15)
theme_set(thm) # for presentation clarity

## Simulate data
d = 2
n = 100
m = 100
set.seed(20)
family = poisson()
sim.data = gmf.simulate(n=n,m=m,d=d,p=2,family = family)

## Present data
heatmap(sim.data$Y, Colv = NA, Rowv = NA, scale="column")
heatmap(family$linkinv(t(sim.data$U)%*%sim.data$V), Colv = NA, Rowv = NA, scale="column")

# Cross validate deviance for a given model training function
cv.deviance = function(X,Y,func,family,nfolds = 10,test.fraction=0.05){
  res = c()
  for (folds in 1:nfolds){
    set.seed(folds)
    smpl = sample(length(Y), test.fraction * length(Y))
    Y.obs = Y
    Y.obs[smpl] = NA
    Y.heldout = Y[smpl]
    model = func(X = X, Y = Y.obs, family = family)
    res = c(res, matrix.deviance(model$fit[smpl],
                              Y.heldout,
                              family = family))
  }
  model = func(X = X, Y = Y, family = family)
  list(model = model, cv.dev = median(res), cv.dev.sd = sd(res), dev = model$deviance)
}

# Cross validate across different lambdas
cv.gmf = function(X,Y,d,trueM,lambdas = 1:30, nfolds = 10, method = "newton"){
  family = poisson()
  paths = c()
  trueScores = c()
  estScores = c()
  insampleDeviance = c()
  models = list()
  for (i in 1:length(lambdas)){
    print(i)
    penaltyU = penaltyV = lambdas[i]
    
    gamma = 0.1
    print(gamma)
    
    func = function(X,Y,family){
      gmf(Y = Y,
                 d = min(d,ceiling(sqrt(ncol(Y)))),
                 family = family,
                 penaltyU = penaltyU,
                 penaltyV = penaltyV,
                 maxIter = 2000,
                 method = method,
                 gamma = gamma,
                 normalize_uv = TRUE,
                 tol = 1e-3,
                 init = "random")
    }
    res = cv.deviance(NULL,Y,func,family,nfolds = nfolds,test.fraction = 0.05)
    models[[i]] = res$model
    
    print(mean(family$dev.resids(Y,family$linkinv(res$model$fit),1)))
    print(norm(res$model$v,"F"))
    print(norm(res$model$u,"F"))
    
    svals = svd(res$model$v)$d
    paths = rbind(paths, svals)
      
    plot(svals,ylim=c(0,max(svals)))
    estScores = c(estScores, res$cv.dev)
  }
  list(model = models[[which.min(estScores)]],
       estScores = estScores,
       paths = paths)
}

# Generate 20 lambdas between 0.5 and e^4
nlambdas = 20
lambdas = seq(0.5,exp(4),length.out=nlambdas)

# Run cross-validation
res.newton = cv.gmf(NULL,sim.data$Y,10,sim.data$M, lambdas = lambdas, nfolds = 10, method = "quasi")
plot(svd(res.newton$model$v)$d)

# Out-sample deviance as a function of lambda
bestLambda = lambdas[which.min(res.newton$estScores)]
g = ggplot( data=data.frame(gamma=lambdas, deviance=res.newton$estScores), aes(x=gamma, y=deviance)) +
  ylim(1.0,2.5) +
  xlab("regularization parameter") + 
  geom_line( size=1.5) +
  geom_vline(xintercept = bestLambda, linetype="dashed", 
             color = "red", size=0.75)
g
ggsave("figures/gamma-deviance.pdf",g,width = 4, height = 3)

# Paths of singular values vs lambda
paths = c()
for (i in 1:ncol(res.newton$paths))
  paths = rbind(paths, cbind(lambdas,i,res.newton$paths[,i]))
colnames(paths) = c("gamma","num","singval")
paths = data.frame(paths)
paths$num = as.factor(paths$num)

g = ggplot( data=paths, aes(x=gamma, y=singval, group=num, color=num)) +
  geom_line( size=1) +
  xlab("regularization parameter") + 
  ylab("singular value") + 
  scale_y_continuous(labels=scaleFUN) +
  guides(color=FALSE) +
  geom_vline(xintercept = bestLambda, linetype="dashed", 
             color = "red", size=0.75)
g
ggsave("figures/sv-paths.pdf",g,width = 4, height = 3)

# Compared with the null model
func.null = function(X,Y,family){
  M = Y
  M[] = mean(M,na.rm=TRUE)
  list(fit = M, deviance = matrix.deviance(M,Y,family = family))
}
res.null = cv.deviance(sim.data$X,sim.data$Y,func.null,family)
dev.null = matrix.deviance(mean(sim.data$Y),sim.data$Y,family = poisson())
print(dev.null)

# Cross validate d
cv.devs = c()
ds = 1:50
paths = matrix(NA, length(ds), length(ds))
dd = c()

for (d in ds){
  func = function(X,Y,family){
    gmf(X = X,
               Y = Y,
             d = d,
             family = family,
             penaltyU = 1,
             penaltyV = 0,
             maxIter = 2000,
             method = "quasi",
             gamma = 0.01,
             normalize_uv = TRUE,
             tol = 1e-3,
             init = "random")
  }
  res = cv.deviance(NULL,sim.data$Y,func,family,nfolds = 10,test.fraction = 0.05)
  svals = svd(res$model$v)$d
  dd = rbind(dd, data.frame(d=d, i=1:d, value=svals))
  
  plot(svals,ylim=c(0,max(svals)))
  cv.devs = c(cv.devs, res$cv.dev)
}

# Plot off-sample deviance as a function of d
bestD = which(min(df$dev) == df$dev)
df = data.frame(num=ds,dev=cv.devs)
g = ggplot( data=df, aes(x=num, y=dev)) +
  geom_point( size=1) +
  geom_line( size=0.2) +
  ylab("deviance") + 
  xlab("number of latent factors") + 
  ylim(1,2.5) +
  geom_vline(xintercept = bestD, linetype="dashed", 
             color = "dark green", size=0.75) +
  guides(color=FALSE)
g
ggsave("figures/dev-by-d.pdf",g,width = 4, height = 3)

# Plot paths of singular values as a function of d
dd$i = factor(dd$i)
g = ggplot( data=dd, aes(x=d, y=value, group=i, color=i)) +
  geom_line( size=1) +
  ylab("singular value") + 
  guides(color=FALSE) +
  xlab("number of latent factors") + 
  scale_y_continuous(labels=scaleFUN) +
  geom_vline(xintercept = bestD, linetype="dashed", 
             color = "dark green", size=0.75)
g
ggsave("figures/d-paths.pdf",g,width = 4, height = 3)
