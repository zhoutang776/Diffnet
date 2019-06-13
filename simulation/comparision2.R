rm(list = ls())
library(diffnet)
source("../comparision/refinement.R")

iter = 1
ii = 1

table2 = data.frame(matrix(0, iter, 5))
colnames(table2) = c("fro", "infinity", "size", "rate", "time")
aslasso.table = asscad.table = asmcp.table = lasso.table = scad.table = mcp.table = table2
refined.table = table2

## ---------------------------- preprocess -----------------------------------------
n_X = 100; n_Y = n_X;
p_X = 100; p_Y = p_X;
nlambda = 50
tuning = "bic"
case = "case1"
stop.tol = 5e-8
perturb = F
correlation = F
lambda.min.ratio = 0.5
max.iter = 600
savefile = F

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

    lasso = scad = mcp = NULL
    # ---------------------------------- symmetric  method ----------------------------------


    start = proc.time()[3]
    lasso = diffnet(X, Y, verbose=F, nlambda=nlambda, max.iter=max.iter, lambda.min.ratio=lambda.min.ratio, stop.tol=stop.tol, method = "lasso", perturb=perturb, correlation=correlation)
    lasso.table[ii,5] = proc.time()[3] - start
    lasso = select(lasso)
    id = lasso$opt[[1]][1]
    res = compute.table(lasso$path[[id]], diff.Omega)
    lasso.table[ii,1:4] = res
    print("lasso done!")

    refined_Delta = symmetric_refinement(lasso$path[[id]], lasso$Sigma.X, lasso$Sigma.Y)
    refined.table[ii,1:4] = compute.table(refined_Delta, diff.Omega)
    print("refine done!")



    start = proc.time()[3]
    scad = diffnet(X, Y, verbose=F, nlambda=nlambda, max.iter=max.iter, lambda.min.ratio=lambda.min.ratio, stop.tol=stop.tol, method="scad", perturb=perturb, correlation=correlation, Delta.init=lasso$path)
    scad.table[ii,5] = proc.time()[3] - start
    scad = select(scad)
    id = scad$opt[[1]][1]
    res = compute.table(scad$path[[id]], diff.Omega)
    scad.table[ii,1:4] = res
    print("scad done!")

    start = proc.time()[3]
    mcp = diffnet(X, Y, verbose=F, nlambda=nlambda, max.iter=max.iter, lambda.min.ratio=lambda.min.ratio, stop.tol=stop.tol, method = "mcp", perturb =perturb, correlation=correlation, Delta.init=lasso$path)
    mcp.table[ii,5] = proc.time()[3] - start
    mcp = select(mcp)
    id = mcp$opt[[1]][1]
    res = compute.table(mcp$path[[id]], diff.Omega)
    mcp.table[ii,1:4] = res
    print("mcp done!")

    aslasso = asscad = asmcp = NULL
    # ---------------------------------- assymmetric  method ----------------------------------

    dyn.load("../comparision/assymmetric_apg/assymmetric_apg.so")
    source("../comparision/assymmetric_apg/assymmetric_apg.R")

    start = proc.time()[3]
    aslasso = assymmetric_apg(X, Y, verbose=F, nlambda=nlambda, max.iter=max.iter, lambda.min.ratio=lambda.min.ratio, stop.tol=stop.tol, method = "lasso", perturb=perturb, correlation=correlation)
    aslasso.table[ii,5] = proc.time()[3] - start
    aslasso = select(aslasso)
    id = aslasso$opt[[1]][1]
    res = compute.table(aslasso$path[[id]], diff.Omega)
    aslasso.table[ii,1:4] = res
    print("aslasso done!")

    start = proc.time()[3]
    asscad = assymmetric_apg(X, Y, verbose=F, nlambda=nlambda, max.iter=max.iter, lambda.min.ratio=lambda.min.ratio, stop.tol=stop.tol, method = "scad", perturb =perturb, correlation=correlation, Delta.init=aslasso$path)
    asscad.table[ii, 5] = proc.time()[3] - start
    asscad = select(asscad)
    id = asscad$opt[[1]][1]
    res = compute.table(asscad$path[[id]], diff.Omega)
    asscad.table[ii,1:4] = res
    print("asscad done!")

    start = proc.time()[3]
    asmcp = assymmetric_apg(X, Y, verbose=F, nlambda=nlambda, max.iter=max.iter, lambda.min.ratio=lambda.min.ratio, stop.tol=stop.tol, method = "mcp", perturb=perturb, correlation=correlation, Delta.init=aslasso$path)
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































lasso.table
diffnet.ic(aslasso$Sigma.X, aslasso$Sigma.Y, aslasso$path, 200, 2)
diffnet.ic(data$Sigma.X, data$Sigma.Y, aslasso$path, 200, 2)


diffnet.ic(data$Sigma.X, data$Sigma.Y, list(lasso$path[[12]]), 200, 2)
# diffnet.ic(data$Sigma.X, data$Sigma.Y, list(diff.Omega), 200, 2)

#











































for(i in 1:nlambda){
    print(c(i,
        nnzero(aslasso$path[[i]]), nnzero(asscad$path[[i]]), nnzero(asmcp$path[[i]]),
        nnzero(lasso$path[[i]]), nnzero(scad$path[[i]]), nnzero(mcp$path[[i]])
        )
    )
}












# # -------------------------------- plot loss -----------------------------------------

# ------------------------------- Frobenius norm --------------------------------
loss.aslasso = loss.asmcp = loss.asscad = loss.mcp = loss.scad = loss.lasso = rep(0, nlambda)

for(i in 1:nlambda){
    loss.lasso[i] = norm(lasso$path[[i]] - diff.Omega, "F")
    loss.scad[i] = norm(scad$path[[i]] - diff.Omega, "F")
    loss.mcp[i] = norm(mcp$path[[i]] - diff.Omega, "F")
    loss.aslasso[i] = norm(aslasso$path[[i]] - diff.Omega, "F")
    loss.asscad[i] = norm(asscad$path[[i]] - diff.Omega, "F")
    loss.asmcp[i] = norm(asmcp$path[[i]] - diff.Omega, "F")
}

lambdas = lasso$lambdas
ymax = max(loss.aslasso, loss.asmcp, loss.asscad, loss.mcp, loss.scad, loss.lasso)
ymin = min(loss.aslasso, loss.asmcp, loss.asscad, loss.mcp, loss.scad, loss.lasso)

plot(lambdas, loss.lasso, ylab = "Frobenius loss", xlab = "lambdas", ylim = c(ymin, ymax))
lines(lambdas, loss.lasso)
# points(lambdas, loss.scad,  col=c("black"))
lines(lambdas, loss.scad,  col=c("black"), lty = 2)
# points(lambdas, loss.mcp,  col=c("black"))
lines(lambdas, loss.mcp,  col=c("black"), lty = 3)

# points(lambdas, loss.aslasso,  col=c("green"))
# lines(lambdas, loss.aslasso,  col=c("green"))
# points(lambdas, loss.asscad,  col=c("green"))
# lines(lambdas, loss.asscad,  col=c("green"), lty = 2)
# points(lambdas, loss.asmcp,  col=c("green"))
# lines(lambdas, loss.asmcp,  col=c("green"), lty = 3)

legend("bottomright", legend=c("lasso","scad", "mcp"), lty = c(1,2,3))












# ----------------------------------------- max norm --------------------------------
lasso_loss = rep(0, nlambda)
assymm_loss = rep(0, nlambda)

for(i in 1:nlambda){
    lasso_loss[i] = max(abs(lasso$path[[i]] - diff.Omega))
    assymm_loss[i] = max(abs(assymm$path[[i]] - diff.Omega))
}

lambdas = lasso$lambdas
# lambdas = max(lambdas) - lambdas
ymax = max(lasso_loss, assymm_loss)
ymin = min(lasso_loss, assymm_loss)

plot(lambdas, lasso_loss, ylab = "infinity loss", xlab = "lambdas", ylim = c(ymin, ymax))
lines(lambdas, lasso_loss)

points(lambdas, assymm_loss,  col=c("green"))
lines(lambdas, assymm_loss,  col=c("green"))

legend("topleft", legend=c("lasso","dtrace"), col=c("black","green"), lty = 1)




# ## ---------------------------- print loss ------------------------------------
#
# for(i in 1:nlambda){
#     print(c(
#         norm(lasso$path[[i]]-diff.Omega, "F"),
#         norm(dtrace$path[[i]]-diff.Omega, "F")
#     )
#     )
# }
#
# for(i in 1:nlambda){
#     print(c(
#         max(abs(lasso$path[[i]]-diff.Omega)),
#         max(abs(dtrace$path[[i]]-diff.Omega))
#     )
#     )
# }



























































