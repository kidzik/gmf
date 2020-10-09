library(latex2exp)
library(gmf)
library(ggplot2)
load("output/simulations-ants-combine.Rda")

deviance.df = data.frame(method = factor(rep(c("GLLVM","AIRWLS","Newton"), each=nrow(df)), levels = c("GLLVM","AIRWLS","Newton")),
                         mse = c(df$res.mse.gllvm, df$res.mse.gmf.airwls, df$res.mse.gmf.newton),
                         procrustes = c(df$res.procrustes.gllvm, df$res.procrustes.gmf.airwls, df$res.procrustes.gmf.newton),
                         deviance = c(df$res.devexp.gllvm, df$res.devexp.gmf.airwls, df$res.devexp.gmf.newton),
                         time = (c(df$res.time.gllvm, df$res.time.gmf.airwls, df$res.time.gmf.newton)))
g = ggplot(deviance.df, aes(x=method, y=deviance, color=method)) + 
  ylab("Mean deviance") + xlab("") +
  ggtitle("Mean deviance") + 
  guides(color=FALSE) +
  geom_boxplot(lwd=1)
g
ggsave("figures/sim-deviance.pdf",width = 4, height = 6)

g = ggplot(deviance.df, aes(x=method, y=procrustes, color=method)) + 
  ylab("Procrustes error") + xlab("") +
  ggtitle("Procrustes error") + 
  guides(color=FALSE) +
  geom_boxplot(lwd=1)
g
ggsave("figures/sim-procrustes.pdf",width = 4, height = 6)

g = ggplot(deviance.df, aes(x=method, y=mse, color=method)) + 
  ylab("MSE") + xlab("") +
  ggtitle("MSE of fixed effects") + 
  guides(color=FALSE) +
  geom_boxplot(lwd=1)
g
ggsave("figures/sim-mse.pdf",width = 4, height = 6)

g = ggplot(deviance.df, aes(x=method, y=time, color=method)) + 
  ylab(TeX("Computation time \\[seconds\\]")) + xlab("") +
  ggtitle("Computation time") + 
  scale_y_log10() +
  guides(color=FALSE) +
  geom_boxplot(lwd=1)
g
ggsave("figures/sim-time.pdf",width = 4, height = 6)

apply(cbind(df$res.time.gllvm/60, df$res.time.gmf.airwls/60, df$res.time.gmf.newton/60),2,mean)

library("plyr")
ddply(df, c("n", "m", "d"), summarise,
      mean = mean(res.devexp.gmf.airwls, na.rm=TRUE)
)