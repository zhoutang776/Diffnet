diffnet = function(X, Y, lambdas=NULL, lambda.min.ratio = NULL, nlambda = NULL, a = NULL,
                         method = c("lasso", "scad", "mcp"), tuning=c("none","aic","bic"),
                         perturb = FALSE, stop.tol=1e-5, max.iter=500,
                         correlation=FALSE, verbose=FALSE, Delta.init = NULL){
    gcinfo(FALSE)
    fit = list()
    fit$method = method[1]

    n.X = nrow(X)
    n.Y = nrow(Y)
    fit$n.X = n.X
    fit$n.Y = n.Y

    if(ncol(X) == ncol(Y)){
        p = ncol(X)
    }else{
        cat("the dimension of the matrices are inconsistent!\n")
        return(NULL)
    }
    if(correlation){
        fit$Sigma.X = cor(X); fit$Sigma.Y = cor(Y);
        X = scale(X, center = T, scale = T) / sqrt(n.X-1) # Sigma.X = t(X)%*%X
        Y = scale(Y, center = T, scale = T) / sqrt(n.Y-1) # Sigma.Y = t(Y)%*%Y
    }else{
        fit$Sigma.X = cov(X)*(1 - 1/n.X); fit$Sigma.Y = cov(Y)*(1 - 1/n.Y)
        X = scale(X, center = T, scale = F) / sqrt(n.X)  # Sigma.X = t(X)%*%X
        Y = scale(Y, center = T, scale = F) / sqrt(n.Y)  # Sigma.Y = t(Y)%*%Y
    }

    if(perturb){
        ## same perturbation as the clime software
        cat("perturbing...")
        eigvals.X = eigen(fit$Sigma.X, only.values=TRUE)$values;
        eigvals.Y = eigen(fit$Sigma.Y, only.values=TRUE)$values;
        epsilon.X = max(max(eigvals.X)-p*min(eigvals.X), 0)/(p-1);
        epsilon.Y = max(max(eigvals.Y)-p*min(eigvals.Y), 0)/(p-1);
        fit$Sigma.X = fit$Sigma.X + epsilon.X*diag(p)
        fit$Sigma.Y = fit$Sigma.Y + epsilon.Y*diag(p)
        cat("done\n")
    }else{ epsilon.X = 0; epsilon.Y = 0; }

    lip = eigen(fit$Sigma.X, symmetric = T, only.values = T)$value[1] * eigen(fit$Sigma.Y, symmetric = T, only.values = T)$value[1]

    if(verbose){ cat(c("the lip constant is :", lip, "\n"))}

    fit$lip = lip

    if(!is.null(lambdas)) {
        nlambda = length(lambdas)
    }else{
        if(is.null(nlambda)) nlambda = 10
        if(is.null(lambda.min.ratio)) lambda.min.ratio = 0.3
        lambda.max = max(abs(fit$Sigma.X - fit$Sigma.Y))
        lambda.min = lambda.min.ratio*lambda.max
        lambdas = exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
    }
    lambdas = sort(lambdas, decreasing = T)
    fit$lambdas = lambdas

    if(is.null(a)) {
        if(method[1] == "scad") a = 3.7;
        if(method[1] == "mcp") a = 3;
    }


    fit$path = vector("list", nlambda)
    fit$sparsity = rep(0, nlambda)
    fit$iter = rep(0, nlambda)

    if(is.null(Delta.init)) Delta = matrix(0, p, p)
    else if(is.matrix(Delta.init)) Delta = Delta.init

    start = proc.time()[3]
    for(i in 1:length(lambdas)){
        lambda = lambdas[i]

        iter = 0
        Lambda = matrix(lambda,p,p)

        out = .C("diffnet_lasso", as.double(fit$Sigma.X), as.double(fit$Sigma.Y), as.integer(p),
                 as.double(Lambda), as.double(X), as.double(Y), as.integer(n.X), as.integer(n.Y),
                 as.double(epsilon.X), as.double(epsilon.Y),
                 as.double(lip), as.double(stop.tol), as.integer(max.iter), iter = as.integer(iter), as.integer(verbose),
                 Delta = as.double(Delta))
        Delta = matrix(out$Delta, ncol=p)
        fit$iter[i] = out$iter

        if(method[1] == "scad"){
            out = .C("diffnet_scad", as.double(fit$Sigma.X), as.double(fit$Sigma.Y), as.integer(p),
                     as.double(lambda), as.double(a), as.double(X), as.double(Y), as.integer(n.X), as.integer(n.Y),
                     as.double(epsilon.X), as.double(epsilon.Y),
                     as.double(lip), as.double(stop.tol), as.integer(max.iter), as.integer(verbose),
                     Delta = as.double(Delta))
            Delta = matrix(out$Delta, ncol=p)
        }
        if(method[1] == "mcp"){
            out = .C("diffnet_mcp", as.double(fit$Sigma.X), as.double(fit$Sigma.Y), as.integer(p),
                     as.double(lambda), as.double(a), as.double(X), as.double(Y), as.integer(n.X), as.integer(n.Y),
                     as.double(epsilon.X), as.double(epsilon.Y),
                     as.double(lip), as.double(stop.tol), as.integer(max.iter), as.integer(verbose),
                     Delta = as.double(Delta))
            Delta = matrix(out$Delta, ncol=p)
        }

        fit$path[[i]] = Matrix::Matrix(Delta, sparse = T)
        fit$sparsity[i] = sum(Delta != 0)/p/(p-1)
    }
    fit$elapse =  proc.time()[3] - start

    if(!is.null(tuning)){
        fit = diffnet.select(fit, tuning = tuning[1])
    }

    rm(X, Y, Delta)
    gc()
    class(fit) = "diffnet"
    return(fit)
}

print.diffnet = function(x, ...)
{
    if(x$method == "lasso")
        cat("Model: graph estimation via lasso\n")
    if(x$method == "scad")
        cat("Model: graph estimation via SCAD\n")
    if(x$method == "mcp")
        cat("Model: graph estimation via MCP\n")

    cat("Path length:",length(x$lambda),"\n")
    cat("Graph dimension:",ncol(x$Sigma.X),"\n")
    cat("lambda:", min(x$lambdas), "----->", max(x$lambdas),"\n")
    cat("Sparsity level:", min(x$sparsity), "----->", max(x$sparsity),"\n")
}

