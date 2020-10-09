devtools::install(".")
library("gmf")
library(gllvm)
library(reshape2)
library(ggplot2)

set.seed(0)

# Plotting
options(repr.plot.width=4, repr.plot.height=4) # set the size of the plots
theme_set(theme_minimal(base_size = 15)) # choose a theme and large font for presentation clarity

# Setup
data(antTraits)
Y = as.matrix(antTraits$abund)
X = as.matrix(antTraits$env)
X = X[,-3,drop=FALSE]
p = 2
fm = poisson()
names(antTraits$env)
tol = 1e-3

## Newton Method
ptm <- proc.time()
model.gmf.newton = gmf(Y = Y, X = X, p = p,
                       gamma=0.2, maxIter = 1000,
                       family = fm,
                       method = "quasi",
                       tol = tol*1e-1)
time.gmf.newton = proc.time() - ptm

## AIRWLS Method
ptm <- proc.time()
model.gmf.airwls = gmf(Y = Y, X = X, p = p,
                       gamma=0.1, maxIter = 10,
                       family = fm, 
                       parallel = 1,
                       penaltyU = 0.01,
                       method = "airwls",
                       tol = tol)
time.gmf.airwls = proc.time() - ptm

## GLLVM Method
ptm <- proc.time()
model.gllvm = gllvm(y = Y, X = X, formula = ~ .,
                    family=fm$family, num.lv = p,
                    beta0com = TRUE,
                    reltol = tol)
time.gllvm = proc.time() - ptm

# Print time
print(paste("Time Newton:",time.gmf.newton[3],"seconds"))
print(paste("Time AIRWLS:",time.gmf.airwls[3],"seconds"))
print(paste("Time gllvm:",time.gllvm[3],"seconds"))

# Print deviances
M.gllvm = fm$linkinv(residuals(model.gllvm)$linpred)
M.gmf.newton = model.gmf.newton$fit
M.gmf.airwls = model.gmf.airwls$fit

devnull = matrix.deviance(pred = mean(Y), obs = Y, family = fm)
print(1 - matrix.deviance(pred = M.gmf.airwls, obs = Y, family = fm)/devnull)
print(1 - matrix.deviance(pred = M.gmf.newton, obs = Y, family = fm)/devnull)
print(1 - matrix.deviance(pred = M.gllvm, obs = Y, family = fm)/devnull)

options(repr.plot.width=4, repr.plot.height=4)
cuts=c(0,5,10,15,20,25,30)

plot_abundance = function(Y, title = "Observed abundance of species", legend = "Count"){
  longData<-melt(Y)
  longData$value = as.numeric(longData$value)

  b = seq(0,30,5)
  
  ggplot(longData, aes(x = Var2, y = Var1)) + 
    geom_raster(aes(fill=value)) + 
    scale_fill_gradient() +
    labs(x="Species", y="Sites", title=title) +
    scale_x_continuous(expand = c(0, 0), limits = range(longData$Var2)) +
    scale_y_continuous(expand = c(0, 0), limits = range(longData$Var1)) +
    guides(fill=guide_legend(title=legend,reverse = TRUE)) +
    scale_fill_gradientn(limits = c(0,30),
                         colours=c("navyblue", "gold"),
                         breaks=b, labels=format(b)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank()
          )
}

width = 5.5
height = 4

plot_abundance(Y)
ggsave("figures/abundance.pdf",width = width, height = height)

plot_abundance(M.gllvm, title = "GLLVM", legend = "Mean")
ggsave("figures/GLLVM.pdf",width = width, height = height)

plot_abundance(M.gmf.airwls, title = "AIRWLS", legend = "Mean")
ggsave("figures/AIRWLS.pdf",width = width, height = height)

plot_abundance(M.gmf.newton, title = "Newton", legend = "Mean")
ggsave("figures/Newton.pdf",width = width, height = height)


pdf("figures/Newton.pdf",width=4.8,4.7)
plot(raster(M.gmf.newton),lab.breaks=cuts,zlim=c(0,max(cuts)))
title("Newton")
dev.off()

pdf("figures/AIRWLS.pdf",width=4.8,4.7)
plot(raster(M.gmf.airwls),lab.breaks=cuts,zlim=c(0,max(cuts)))
title("AIRWLS")
dev.off()

# Plot PC1-PC2 plane
U = model.gmf.newton$u
colnames(U) = c("PC1","PC2")
rownames(U) = paste("Site",rownames(antTraits$env))
U = data.frame(U)
U$Shrub.cover = antTraits$env[,3]

p<-ggplot(U,aes(x=-PC1,y=PC2,color=Shrub.cover,label=row.names(U) ))+
  geom_point(size=2) +
  scale_colour_gradientn(colours=c("navyblue", "gold"))
p
ggsave("figures/pc-plane.pdf", width= 5.9, height = 6)

# Test the correlation with Shrub.cover
cor.test(U$Shrub.cover,U$PC2)

# Plot correlation with Shrub.cover
options(repr.plot.width=3, repr.plot.height=4)
p = ggplot(U, aes(y = PC2, x = Shrub.cover)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "blue")
p
ggsave("figures/pc-shrub.pdf", width= 4, height = 6)

  
