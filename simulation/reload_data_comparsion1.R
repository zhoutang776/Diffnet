rm(list = ls())
library(diffnet)
# set.seed(1)
setwd("/home/zhoutang/paper/Dropbox/Diffnet/comparision")
source("../comparision/refinement.R")
iter = 1
ii = 1
infity.lasso = matrix(0, nrow = iter, ncol = 3,  dimnames = list(c(),c("TP", "TN", "TD")) )

infity.dantzig = infity.admm = infity.asadmm = infity.aslass   = infity.lasso
fro.dantzig    = fro.admm    = fro.asadmm    = fro.aslasso    = fro.lasso = infity.lasso
time.dantzig   = time.admm   = time.asadmm   = time.aslasso   = time.lasso = rep(0, iter)
iter.dantzig   = iter.admm   = iter.asadmm   = iter.aslasso   = iter.lasso = rep(0, iter)

## ---------------------------- preprocess -----------------------------------------
n_X = 200; n_Y = n_X;
p_X = 100; p_Y = p_X;
nlambda = 50
tuning = "none"
case = "case1"
stop.tol = 1e-5
perturb = F
correlation = F
max.iter = 800
lambda.min.ratio = 0.5
for(ii in 1:iter){
    # ------------ data generating -----------------
    data = diffnet.case(n=n_X, p=p_X, method = case)

    X = data$X
    Y = data$Y
    diff.Omega = data$diff.Omega
    print(sum(diff.Omega!=0))


    # ---------------------------------- new 1 (assymmetric apg) method ----------------------------------
    aslasso = NULL
    # dyn.load("../comparision/assymmetric_apg/assymmetric_apg.so")
    # source("../comparision/assymmetric_apg/assymmetric_apg.R")
    # start = proc.time()[3]
    # aslasso = assymmetric_apg(X, Y, verbose = F, nlambda = nlambda, max.iter = max.iter, lambda.min.ratio = lambda.min.ratio, stop.tol=stop.tol, method = "lasso", perturb =perturb, correlation = correlation)
    # time.aslasso[ii] = proc.time()[3] - start
    # iter.aslasso[ii] = sum(aslasso$iter)
    # print("assymmetric apg done!")
    # # aslasso = diffnet.select(aslasso, tuning = tuning)

    # ---------------------------------- new 2 (symmetric apg) method ----------------------------------
    lasso = NULL
    # start = proc.time()[3]
    # lasso = diffnet(X, Y, verbose = F, nlambda = nlambda, max.iter = max.iter, lambda.min.ratio = lambda.min.ratio, stop.tol=stop.tol, method = "lasso", perturb =perturb, correlation = correlation)
    # time.lasso[ii] = proc.time()[3] - start
    # iter.lasso[ii] = sum(lasso$iter)
    # print("symmetric apg done!")
    # # lasso = diffnet.select(lasso, tuning = tuning)


    # ---------------------------------- ADMM 1 (assymmetric)  method ----------------------------------
    asadmm = NULL
    dyn.load("../comparision/assymmetric_admm/assymmetric_admm.so")
    source("../comparision/assymmetric_admm/assymmetric_admm.R")

    # start = proc.time()[3]
    # asadmm = assymmetric_admm(X, Y, nlambda=nlambda, max.iter = max.iter, lambda.min.ratio=lambda.min.ratio, stop.tol=stop.tol, perturb=perturb, correlation = correlation)
    # time.asadmm[ii] = proc.time()[3] - start
    # iter.asadmm[ii] = sum(asadmm$iter)
    # print("assymmetric admm done!")
    # # asadmm = diffnet.select(asadmm, tuning = tuning)

    # ---------------------------------- ADMM 2 (symmetric Xi)  method ----------------------------------
    admm = NULL
    # dyn.load("../comparision/symmetric_admm/symmetric_admm.so")
    # source("../comparision/symmetric_admm/Functionsnew.R")
    #
    # start = proc.time()[3]
    # admm = symmetric_admm(X, Y, nlambda=nlambda, max.iter = max.iter, lambda.min.ratio=lambda.min.ratio, stop.tol=stop.tol, perturb=perturb, correlation = correlation)
    # time.admm[ii] = proc.time()[3] - start
    # iter.admm[ii] = sum(admm$iter)
    # print("symmetric admm done!")
    # # admm = diffnet.select(admm, tuning = tuning)


    # # ---------------------------------- Dantzig (zhao) method ----------------------------------
    dantzig = NULL
    dyn.load("../comparision/dantzig/dpm.so")
    source("../comparision/dantzig/dpm.R")
    start = proc.time()[3]
    dantzig = dpm(X, Y, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, perturb =perturb, correlation = correlation)
    time.dantzig[ii] = proc.time()[3] - start
    dantzig = diffnet.select(dantzig, tuning = tuning)

}

time.dantzig
dantzig$iter





print("======================================")
time.asadmm
time.admm
time.aslasso
time.lasso

iter.asadmm
iter.admm
iter.aslasso
iter.lasso



mean(time.asadmm)
mean(time.admm)
mean(time.aslasso)
mean(time.lasso)

sd(time.asadmm)
sd(time.admm)
sd(time.aslasso)
sd(time.lasso)

print("======================================")


mean(iter.asadmm)
mean(iter.admm)
mean(iter.aslasso)
mean(iter.lasso)

sd(iter.asadmm)
sd(iter.admm)
sd(iter.aslasso)
sd(iter.lasso)








print("======================================")
sum(aslasso$iter)/time.aslasso
sum(lasso$iter)/time.lasso
sum(asadmm$iter)/time.asadmm
sum(admm$iter)/time.admm



























p100case1 = c(533.8, 1899.6, 791.6, 784.2)
p100case2 = c(818.4, 2980.2, 968, 917.6)
p200case1 = c(832.4, 2598.6, 757.2, 755.4)
p200case2 = c(1413.8, 4488.6, 1113, 1040.8)
p400case1 = c(1422.4, 4562.6, 1024.8, 957.8)
p400case2 = c(1417.7, 6268, 1194, 1089.2)
p600case1 = c(1754.4, 3913, 1317.2, 1274.6)
p600case2 = c(1675, 7128.2, 1228.2, 1158)
p800case1 = c(1125, 4601.6, 1147.4, 1092)
p800case2 = c(2463.4, 6538.2, 1468.4, 1388.2)
data_summary = rbind(p100case1, p100case2, p200case1, p200case2, p400case1, p400case2, p600case1, p600case2, p800case1, p800case2)
case1_iter = data_summary[c(1,3,5,7,9),]
case2_iter = data_summary[c(2,4,6,8,10),]
colnames(case1_iter) = c("asadmm", "admm", "aslasso", "lasso")
colnames(case2_iter) = c("asadmm", "admm", "aslasso", "lasso")


