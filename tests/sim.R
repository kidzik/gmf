library(gmf)

data = gmf.simulate()

model = gmf(data$Y)
plot(model$u)
plot(model$v)
print(matrix.deviance(model$fit, data$Y, model$family))
