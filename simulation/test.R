rm(list = ls())
library(diffnet)

setwd("/home/johntan/Desktop/Dropbox/Diffnet/diffnet/")

iter = 1
ii = 1

infity.lasso = matrix(0, nrow = iter, ncol = 3,  dimnames = list(c(),c("TP", "TN", "TD")) )

infity.dantzig = infity.admm = infity.asadmm = infity.aslass   = infity.lasso
fro.dantzig    = fro.admm    = fro.asadmm    = fro.aslasso    = fro.lasso = infity.lasso
time.dantzig   = time.admm   = time.asadmm   = time.aslasso   = time.lasso = rep(0, iter)

## ---------------------------- preprocess -----------------------------------------
n_X = 200; n_Y = n_X;
p_X = 100; p_Y = p_X;
nlambda = 50
tuning = "none"
case = "case1"
stop.tol = 1e-5
perturb = F
correlation = F
max.iter = 100
lambda.min.ratio = 0.5


data = diffnet.case(n=n_X, p=p_X, method = case)

X = data$X
Y = data$Y
diff.Omega = data$diff.Omega
print(sum(diff.Omega!=0))



lasso = NULL
start = proc.time()[3]
lasso = diffnet(X, Y, verbose = F, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, max.iter = max.iter, stop.tol=stop.tol, method = "lasso", perturb =perturb, correlation = correlation)
time.lasso[ii] = proc.time()[3] - start
print("symmetric apg done!")


aslasso = NULL
dyn.load("../comparision/assymmetric_apg/assymmetric_apg.so")
source("../comparision/assymmetric_apg/assymmetric_apg.R")
start = proc.time()[3]
aslasso = assymmetric_apg(X, Y, verbose = F, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, max.iter = max.iter, stop.tol=stop.tol, method = "lasso", perturb =perturb, correlation = correlation)
time.aslasso[ii] = proc.time()[3] - start
print("assymmetric apg done!")


# ---------------------------------- new 2 (symmetric apg) method ----------------------------------



for(i in 1:nlambda){
  print(
    sprintf("%2i %3g %3g", i,
            nnzero(lasso$path[[i]]), nnzero(aslasso$path[[i]])
    )
  )
}



