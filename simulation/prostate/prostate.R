rm(list = ls())
library(diffnet)
library(pracma)

prostmat <- read.csv(file = "../comparision/prostate/prostmat.csv", header = T)
prostz <- read.csv(file = "../comparision/prostate/prostz.txt")
control <- t(prostmat[,1:50])
cancer <- t(prostmat[, 51:102])



## ---------------------------- setting -------------------------------------
nlambda = 2
lambda.min.ratio = 0.8
stop.tol = 1e-3
perturb = T
correlation = T
tuning = "aic"

start = proc.time()
fit.lasso = diffnet(cancer, control, verbose = T, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, stop.tol=stop.tol, method = "lasso", perturb = perturb, correlation = correlation)
proc.time() - start

fit.lasso = diffnet.select(fit.lasso, tuning = tuning)


