rm(list = ls())
library(diffnet)
source("../comparision/refinement.R")

dyn.load("../comparision/assymmetric_apg/assymmetric_apg.so")
source("../comparision/assymmetric_apg/assymmetric_apg.R")
dyn.load("../comparision/symmetric_admm/symmetric_admm.so")
source("../comparision/symmetric_admm/Functionsnew.R")

iter = 10
ii = 1

table2 = data.frame(matrix(0, iter, 5))
colnames(table2) = c("fro", "infinity", "size", "rate", "time")
aslasso.table = lasso.table = table2
refined.table = table2
admm.table = table2

## ---------------------------- preprocess -----------------------------------------
n_X = 100; n_Y = n_X;
p_X = 1000; p_Y = p_X;
nlambda = 50
tuning = "aic"
case = "case2"
stop.tol = 1e-3
perturb = F
correlation = F
lambda.min.ratio = 0.5
max.iter = 50
savefile = T

for(ii in 1:iter){
    # ------------ data generating -----------------
    data = diffnet.case(n=n_X, p=p_X, method = case)

    X = data$X
    Y = data$Y
    diff.Omega = data$diff.Omega

    print(sum(diff.Omega!=0))
    select = function(fit){
        return(diffnet.select(fit, data$Sigma.X, data$Sigma.Y, tuning))
    }

    lasso = NULL
    # ---------------------------------- symmetric  method ---------------------------------

    start = proc.time()[3]
    lasso = diffnet(X, Y, verbose=F, nlambda=nlambda, max.iter=max.iter, lambda.min.ratio=lambda.min.ratio, stop.tol=stop.tol, method = "lasso", perturb=perturb, correlation=correlation)
    lasso.table[ii,5] = proc.time()[3] - start
    lasso = select(lasso)
    id = max(lasso$opt[[1]][1], lasso$opt[[1]][2])
    res = compute.table(lasso$path[[id]], diff.Omega)
    lasso.table[ii,1:4] = res
    print("lasso done!")
    print(lasso.table[ii,])

    refined_Delta = symmetric_refinement(lasso$path[[id]], lasso$Sigma.X, lasso$Sigma.Y)
    refined.table[ii,1:4] = compute.table(refined_Delta, diff.Omega)
    print("refine done!")
    print(refined.table[ii,])


    aslasso = NULL
    # # ---------------------------------- assymmetric  method ----------------------------------

    start = proc.time()[3]
    aslasso = assymmetric_apg(X, Y, verbose=F, nlambda=nlambda, max.iter=max.iter, lambda.min.ratio=lambda.min.ratio, stop.tol=stop.tol, method = "lasso", perturb=perturb, correlation=correlation)
    aslasso.table[ii,5] = proc.time()[3] - start
    aslasso = select(aslasso)
    id = max(lasso$opt[[1]][1], lasso$opt[[1]][2])
    res = compute.table(aslasso$path[[id]], diff.Omega)
    aslasso.table[ii,1:4] = res
    print("aslasso done!")
    print(aslasso.table[ii,])



    admm = NULL
    # ---------------------------------- symmetric admm method ----------------------------------
    # start = proc.time()[3]
    # admm = symmetric_admm(X, Y, nlambda=nlambda, lambda.min.ratio=lambda.min.ratio, max.iter = max.iter, perturb=perturb, correlation = correlation, stop.tol=stop.tol)
    # admm.table[ii,5] = proc.time()[3] - start
    # print("symmetric admm done!")
    # admm = select(admm)
    # id = max(lasso$opt[[1]][1], lasso$opt[[1]][2])
    # res = compute.table(admm$path[[id]], diff.Omega)
    # admm.table[ii,1:4] = res
    # print("admm done!")
    # print(admm.table[ii,])


    print("------------------------------------------------------------")


    #------------------------ save something ---------------------------------
    if(savefile){
        for(i in 1:nlambda){
            aslasso$path[[i]] = Matrix(aslasso$path[[i]], sparse = T)
            lasso$path[[i]] = Matrix(lasso$path[[i]], sparse = T)
        }
        fileaddress = "/home/johntan/Desktop"
        filename = paste(case, paste("p=", p_X, sep=""), paste("iter=", ii, sep=""), sep="_")
        filepath = file.path(fileaddress, paste(filename, "RData", sep = "."))

        save(list = c("aslasso", "lasso", "data", "nlambda", "max.iter", "lambda.min.ratio", "stop.tol", "perturb", "correlation"), file = filepath)
    }

}



apply(aslasso.table, 2, mean)
apply(aslasso.table, 2, sd)


apply(lasso.table, 2, mean)
apply(lasso.table, 2, sd)


apply(refined.table, 2, mean)
apply(refined.table, 2, sd)


apply(admm.table, 2, mean)
apply(admm.table, 2, sd)































for(i in 1:nlambda){print(    c(i, nnzero(admm$path[[i]]), nnzero(lasso$path[[i]]))   )}


## ---------------------------------- using admm to check algorithm -----------------------------------
dyn.load("../comparision/symmetric_admm/symmetric_admm.so")
source("../comparision/symmetric_admm/Functionsnew.R")
ii = 1
admm.table = table2
# ---------------------------------- symmetric admm method ----------------------------------
start = proc.time()[3]
admm = symmetric_admm(X, Y, nlambda=nlambda, lambda.min.ratio=lambda.min.ratio, max.iter = max.iter, perturb=perturb, correlation = correlation, stop.tol=stop.tol)
admm.table[ii,5] = proc.time()[3] - start
print("symmetric admm done!")
admm = select(admm)
id = admm$opt[[1]][1]
res = compute.table(admm$path[[id]], diff.Omega)
admm.table[ii,1:4] = res
print("admm done!")


apply(admm.table, 2, mean)
apply(admm.table, 2, sd)


for(i in 1:nlambda){print(   c(i,nnzero(aslasso$path[[i]]), nnzero(lasso$path[[i]]), nnzero(admm$path[[i]]))    )}

id = 50
lasso$path[[id]][1:5, 1:5]
## ------------------------------------------ end  -----------------------------------
