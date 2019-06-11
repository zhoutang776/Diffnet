rm(list = ls())
library(diffnet)

setwd("/home/johntan/Desktop/Dropbox/Diffnet/diffnet/")
source("../comparision/refinement.R")
iter = 5
ii = 1
infity.lasso = matrix(0, nrow = iter, ncol = 3,  dimnames = list(c(),c("TP", "TN", "TD")) )

infity.dantzig = infity.admm = infity.asadmm = infity.aslass   = infity.lasso
fro.dantzig    = fro.admm    = fro.asadmm    = fro.aslasso    = fro.lasso = infity.lasso
time.dantzig   = time.admm   = time.asadmm   = time.aslasso   = time.lasso = rep(0, iter)

## ---------------------------- preprocess -----------------------------------------
n_X = 100; n_Y = n_X;
p_X = 800; p_Y = p_X;
nlambda = 50
tuning = "none"
case = "case1"
stop.tol = 1e-5
perturb = F
correlation = F
max.iter = 800
lambda.min.ratio = 0.5

# ------------ data generating -----------------
data = diffnet.case(n=n_X, p=p_X, method = case)

X = data$X
Y = data$Y
diff.Omega = data$diff.Omega
print(sum(diff.Omega!=0))


# ---------------------------------- new 1 (assymmetric apg) method ----------------------------------
aslasso = NULL
dyn.load("../comparision/assymmetric_apg/assymmetric_apg.so")
source("../comparision/assymmetric_apg/assymmetric_apg.R")
start = proc.time()[3]
aslasso = assymmetric_apg(X, Y, verbose = F, nlambda = nlambda, max.iter = max.iter, lambda.min.ratio = lambda.min.ratio, stop.tol=stop.tol, method = "lasso", perturb =perturb, correlation = correlation)
time.aslasso[ii] = proc.time()[3] - start
print("assymmetric apg done!")
# aslasso = diffnet.select(aslasso, tuning = tuning)

# ---------------------------------- new 2 (symmetric apg) method ----------------------------------
lasso = NULL
start = proc.time()[3]
lasso = diffnet(X, Y, verbose = F, nlambda = nlambda, max.iter = max.iter, lambda.min.ratio = lambda.min.ratio, stop.tol=stop.tol, method = "lasso", perturb =perturb, correlation = correlation)
time.lasso[ii] = proc.time()[3] - start
print("symmetric apg done!")
# lasso = diffnet.select(lasso, tuning = tuning)

# ---------------------------------- ADMM 1 (assymmetric)  method ----------------------------------
asadmm = NULL
dyn.load("../comparision/assymmetric_admm/assymmetric_admm.so")
source("../comparision/assymmetric_admm/assymmetric_admm.R")

start = proc.time()[3]
asadmm = assymmetric_admm(X, Y, nlambda=nlambda, max.iter = max.iter, lambda.min.ratio=lambda.min.ratio, stop.tol=stop.tol, perturb=perturb, correlation = correlation)
time.asadmm[ii] = proc.time()[3] - start
print("assymmetric admm done!")
# asadmm = diffnet.select(asadmm, tuning = tuning)

# ---------------------------------- ADMM 2 (symmetric Xi)  method ----------------------------------
admm = NULL
dyn.load("../comparision/symmetric_admm/symmetric_admm.so")
source("../comparision/symmetric_admm/Functionsnew.R")

start = proc.time()[3]
admm = symmetric_admm(X, Y, nlambda=nlambda, max.iter = max.iter, lambda.min.ratio=lambda.min.ratio, stop.tol=stop.tol, perturb=perturb, correlation = correlation)
time.admm[ii] = proc.time()[3] - start
print("symmetric admm done!")
# admm = diffnet.select(admm, tuning = tuning)


# # ---------------------------------- Dantzig (zhao) method ----------------------------------
dantzig = NULL
# dyn.load("../comparision/dantzig/dpm.so")
# source("../comparision/dantzig/dpm.R")
# start = proc.time()[3]
# dantzig = dpm(X, Y, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, perturb =perturb, correlation = correlation)
# time.dantzig[ii] = proc.time()[3] - start
# dantzig = diffnet.select(dantzig, tuning = tuning)


time.aslasso
time.lasso
time.asadmm
time.admm




mean(time.asadmm)
sd(time.asadmm)
mean(time.admm)
sd(time.admm)
mean(time.aslasso)
sd(time.aslasso)
mean(time.lasso)
sd(time.lasso)




for(i in 1:nlambda){
    print(
        sprintf("%2i %3g %3g %3g %3g", i,
                nnzero(lasso$path[[i]]), nnzero(aslasso$path[[i]]),
                nnzero(admm$path[[i]]), nnzero(asadmm$path[[i]])
        )
    )
}






diffnet.ic(data$Sigma.X, data$Sigma.Y, lasso$path, 200, 2)
diffnet.ic(lasso$Sigma.X, lasso$Sigma.Y, lasso$path, 200, 2)

diffnet.ic(data$Sigma.X, data$Sigma.Y, asadmm$path, 200, 2)
diffnet.ic(asadmm$Sigma.X, asadmm$Sigma.Y, asadmm$path, 200, 2)



compute.table(lasso$path[[24]], diff.Omega)

compute.table(asadmm$path[[24]], diff.Omega)


# diffnet.ic(data$Sigma.X, data$Sigma.Y, list(diff.Omega), 200, 2)
id = 24
aslasso$path[[id]][1:5,1:5]
lasso$path[[id]][1:5,1:5]
asadmm$path[[id]][1:5,1:5]
admm$path[[id]][1:5,1:5]



















# # -------------------------------- print nnz -----------------------------------------
























































