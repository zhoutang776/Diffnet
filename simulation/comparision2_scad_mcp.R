library(diffnet)
source("../comparision/refinement.R")

dyn.load("../comparision/assymmetric_apg/assymmetric_apg.so")
source("../comparision/assymmetric_apg/assymmetric_apg.R")

iter = 2
ii = 1

table2 = data.frame(matrix(0, iter, 5))
colnames(table2) = c("fro", "infinity", "size", "rate", "time")
aslasso.table = asscad.table = asmcp.table = lasso.table = scad.table = mcp.table = table2
refined.table = table2

## ---------------------------- preprocess -----------------------------------------
p_X = 1000
tuning = "aic"
case = "case1"
fileaddress = "/home/johntan/Desktop"
savefile = F

ii=9
for(ii in 1:iter){
    filename = paste(case, paste("p=", p_X, sep=""), paste("iter=", ii, sep=""), sep="_")
    filepath = file.path(fileaddress, paste(filename, "RData", sep = "."))
    load(file = filepath)

    diff.Omega = data$diff.Omega
    print(sum(diff.Omega!=0))
    select = function(fit){
        return(diffnet.select(fit, data$Sigma.X, data$Sigma.Y, tuning))
    }

    scad = mcp = NULL
    # ---------------------------------- symmetric  method ----------------------------------

    id = lasso$opt[[1]][1]
    res = compute.table(aslasso$path[[id]], diff.Omega)
    lasso.table[ii,1:4] = res
    print("lasso done!")


    refined_Delta = symmetric_refinement(lasso$path[[id]], lasso$Sigma.X, lasso$Sigma.Y)
    refined.table[ii,1:4] = compute.table(refined_Delta, diff.Omega)
    print("refine done!")


    start = proc.time()[3]
    scad = diffnet(data$X, data$Y, verbose=F, nlambda=nlambda, max.iter=max.iter, lambda.min.ratio=lambda.min.ratio, stop.tol=stop.tol, method="scad", perturb=perturb, correlation=correlation, Delta.init=lasso$path)
    scad.table[ii,5] = proc.time()[3] - start
    scad = select(scad)
    id = scad$opt[[1]][1]
    res = compute.table(scad$path[[id]], diff.Omega)
    scad.table[ii,1:4] = res
    print("scad done!")

    start = proc.time()[3]
    mcp = diffnet(data$X, data$Y, verbose=F, nlambda=nlambda, max.iter=max.iter, lambda.min.ratio=lambda.min.ratio, stop.tol=stop.tol, method = "mcp", perturb =perturb, correlation=correlation, Delta.init=lasso$path)
    mcp.table[ii,5] = proc.time()[3] - start
    mcp = select(mcp)
    id = mcp$opt[[1]][1]
    res = compute.table(mcp$path[[id]], diff.Omega)
    mcp.table[ii,1:4] = res
    print("mcp done!")

   asscad = asmcp = NULL
    # ---------------------------------- assymmetric  method ----------------------------------

    id = aslasso$opt[[1]][1]
    res = compute.table(aslasso$path[[id]], diff.Omega)
    aslasso.table[ii,1:4] = res
    print("aslasso done!")


    start = proc.time()[3]
    asscad = assymmetric_apg(data$X, data$Y, verbose=F, nlambda=nlambda, max.iter=max.iter, lambda.min.ratio=lambda.min.ratio, stop.tol=stop.tol, method = "scad", perturb =perturb, correlation=correlation, Delta.init=aslasso$path)
    asscad.table[ii, 5] = proc.time()[3] - start
    asscad = select(asscad)
    id = asscad$opt[[1]][1]
    res = compute.table(asscad$path[[id]], diff.Omega)
    asscad.table[ii,1:4] = res
    print("asscad done!")

    start = proc.time()[3]
    asmcp = assymmetric_apg(data$X, data$Y, verbose=F, nlambda=nlambda, max.iter=max.iter, lambda.min.ratio=lambda.min.ratio, stop.tol=stop.tol, method = "mcp", perturb=perturb, correlation=correlation, Delta.init=aslasso$path)
    asmcp.table[ii,5] = proc.time()[3] - start
    asmcp = select(asmcp)
    id = asmcp$opt[[1]][1]
    res = compute.table(asmcp$path[[id]], diff.Omega)
    asmcp.table[ii,1:4] = res
    print("asmcp done!")

    #------------------------ save something ---------------------------------
    if(savefile){
        for(i in 1:nlambda){
            aslasso$path[[i]] = Matrix(aslasso$path[[i]], sparse = T)
            asscad$path[[i]] = Matrix(asscad$path[[i]], sparse = T)
            asmcp$path[[i]] = Matrix(asmcp$path[[i]], sparse = T)

            lasso$path[[i]] = Matrix(lasso$path[[i]], sparse = T)
            scad$path[[i]] = Matrix(scad$path[[i]], sparse = T)
            mcp$path[[i]] = Matrix(mcp$path[[i]], sparse = T)
        }
        fileaddress = "/home/johntan/Desktop"
        filename = paste(case, paste("p=", p_X, sep=""), paste("iter=", ii, sep=""), sep="_")
        filepath = file.path(fileaddress, paste(filename, "RData", sep = "."))
        save(list = c("aslasso", "asscad", "asmcp", "lasso", "scad", "mcp", "data"), file = filepath)
    }

}


apply(aslasso.table, 2, mean)
apply(aslasso.table, 2, sd)

apply(asscad.table, 2, mean)
apply(asscad.table, 2, sd)

apply(asmcp.table, 2, mean)
apply(asmcp.table, 2, sd)

apply(lasso.table, 2, mean)
apply(lasso.table, 2, sd)

apply(scad.table, 2, mean)
apply(scad.table, 2, sd)

apply(mcp.table, 2, mean)
apply(mcp.table, 2, sd)

apply(refined.table, 2, mean)
apply(refined.table, 2, sd)

for(i in 1:nlambda){
    print(c(i,
            nnzero(aslasso$path[[i]]), nnzero(asscad$path[[i]]), nnzero(asmcp$path[[i]]),
            nnzero(lasso$path[[i]]), nnzero(scad$path[[i]]), nnzero(mcp$path[[i]])
    )
    )
}



































