library("gmf")

trials = c()
ns = c()
ms = c()
ds = c()
dists = c()
res = c()
dates = c()

get_results = function(simulation, model.gllvm, model.newton, model.airwls, model.airwls.glmnet){
    M.gmf.newton = simulation$family$linkinv(model.newton$fit)
    M.gmf.airwls = simulation$family$linkinv(model.airwls$fit)
    M.gmf.airwls.glmnet = simulation$family$linkinv(model.airwls.glmnet$fit)
    M.gllvm = tryCatch({
       simulation$family$linkinv(residuals(model.gllvm)$linpred)
    }, error = function(e) {
       dims = dim(M.gmf.airwls.glmnet)
       matrix(0, dims[1], dims[2])
    })

    dev.null = gmf.deviance(mean(simulation$Y), simulation$Y, simulation$family)

    cors = c(    
    cor(c(M.gllvm), c(simulation$linpred)),
    cor(c(M.gmf.newton), c(simulation$linpred)),
    cor(c(M.gmf.airwls), c(simulation$linpred)),
    cor(c(M.gmf.airwls.glmnet), c(simulation$linpred)))

    devs = c(
    gmf.deviance(M.gllvm, simulation$Y, simulation$family),
    gmf.deviance(M.gmf.newton, simulation$Y, simulation$family),
    gmf.deviance(M.gmf.airwls, simulation$Y, simulation$family),
    gmf.deviance(M.gmf.airwls.glmnet, simulation$Y, simulation$family))

    procrustes = c(
        gmf.procrustes(simulation$V, model.gllvm$params$theta)$ss,
        gmf.procrustes(simulation$V, model.newton$v)$ss,
        gmf.procrustes(simulation$V, model.airwls$v)$ss,
        gmf.procrustes(simulation$V, model.airwls.glmnet$v)$ss
    )

    mse = c(
        mean((simulation$beta - model.gllvm$params$Xcoef)**2),
        mean((simulation$beta - t(model.newton$beta))**2),
        mean((simulation$beta - t(model.airwls$beta))**2),
        mean((simulation$beta - t(model.airwls.glmnet$beta))**2)
    )

    c(dev.null, cors,devs,procrustes,mse)
}

i=1
for (file in list.files("output")){
    try({
    load(paste0("output/",file))

    r = c()
    r = get_results(simulation, model.gllvm, model.newton, model.airwls, model.airwls.glmnet)
    r = c(r, time.gllvm[3], time.newton[3], time.airwls[3], time.airwls.glmnet[3])

    basename = strsplit(file,"[.]")[[1]][1]
    vars = strsplit(basename,"-")[[1]]
    print(vars)
    trials = c(trials, as.numeric(vars[1]))
    ns = c(ns, as.numeric(vars[2]))
    ms = c(ms, as.numeric(vars[3]))
    ds = c(ds, as.numeric(vars[4]))
    dists = c(dists, vars[5])
    dates = c(dates, vars[6])

    res = rbind(res, r)
    },
    silent = TRUE)
    i=i+1
    if (i==50)
        break
}

res = data.frame(res)
algs = c("gllvm", "gmf.newton", "gmf.airwls", "gmf.airwls.glmnet")
colnames(res) = c("dev.null", paste("cor",algs,sep="."), paste("devexp",algs,sep="."), paste("procrustes",algs,sep="."), paste("mse",algs,sep="."), paste("time",algs,sep="."))
df = data.frame(trial=trials,n=ns,m=ms,d=ds,date=dates,dist=dists,res=res)

save(df, file="output/simulations-ants-combine.Rda")