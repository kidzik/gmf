library(tidyr)
library(latex2exp)

load("output/iteration-0.01.Rdata")

df$size = 0.005 * 1:nrow(df)

df.dev = df[,c(2:4,8)]
colnames(df.dev) = c("GLLVM","Newton","AIRWLS","size")
df.dev = df.dev %>% pivot_longer(cols = GLLVM:AIRWLS, names_to = "Method")
df.dev$Method = factor(df.dev$Method, levels = c("GLLVM", "AIRWLS", "Newton"))

g = ggplot(df.dev, aes(x=size, y=value, color=Method)) + 
  ylab("Mean deviance") + xlab("Fraction of the dataset") +
  ggtitle("Mean deviance") + 
  geom_line(size=1.5)  + theme(legend.position = "none")
g
ggsave("figures/deviance.pdf",width = 3.5, height = 5)

df.time = df[,c(5:7,8)]
colnames(df.time) = c("GLLVM","Newton","AIRWLS","size")
df.time = df.time %>% pivot_longer(cols = GLLVM:AIRWLS, names_to = "Method")
df.time$Method = factor(df.time$Method, levels = c("GLLVM", "AIRWLS", "Newton"))
#df.time$value = df.time$value*60**2

g = ggplot(df.time, aes(x=size, y=value, color=Method)) + 
  ylab("Time [sec]") + xlab("Fraction of the dataset") +
  ggtitle("Computation time") + 
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", function(x) format(10^x,big.mark=",",scientific=FALSE))
  ) +
  #  ylim(-3,12) + 
  geom_line(size=1.5)  #+ theme(legend.position = "none")
g
ggsave("figures/compute.pdf",width = 5.5, height = 5)
