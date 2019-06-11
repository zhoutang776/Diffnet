assymmetric_admm<- function(X, Y, lambda=NULL,nlambda=10, lambda.min.ratio=NULL,
                rho=1, correlation=FALSE, perturb=FALSE, stop.tol=1e-5, max.iter = 500){
    fit = list()

    n.X = nrow(X);
    n.Y = nrow(Y);
    fit$n.X = n.X
    fit$n.Y = n.Y

    if(ncol(X) == ncol(Y)){
        p = ncol(X)
    }else{
        cat("the dimension of the matrices are inconsistent!\n")
        return(NULL)
    }

    if(correlation){
        fit$Sigma.X = cor(X);
        fit$Sigma.Y = cor(Y);
    }else{
        fit$Sigma.X = cov(X)*(1-1/n.X);
        fit$Sigma.Y = cov(Y)*(1-1/n.Y);
    }

    if(perturb){
        cat("perturbing...")
        eigvals.X = eigen(fit$Sigma.X, only.values=TRUE)$values;
        eigvals.Y = eigen(fit$Sigma.Y, only.values=TRUE)$values;
        ## same perturbation as the clime software
        epsilon.X = max(max(eigvals.X)-p*min(eigvals.X), 0)/(p-1);
        epsilon.Y = max(max(eigvals.Y)-p*min(eigvals.Y), 0)/(p-1);
        fit$Sigma.X = fit$Sigma.X + epsilon.X*diag(p)
        fit$Sigma.Y = fit$Sigma.Y + epsilon.Y*diag(p)
        cat("done\n")
    }
    if(!is.null(lambda)) {
        nlambda = length(lambda)
    }else{
        if(is.null(nlambda)) nlambda = 10
        if(is.null(lambda.min.ratio)) lambda.min.ratio = 0.3
        lambda.max = max(abs(fit$Sigma.X - fit$Sigma.Y))
        lambda.min = lambda.min.ratio*lambda.max
        lambda = exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
    }

    fit$lambda = lambda


    M1 = eigen(fit$Sigma.X)
    M2 = eigen(fit$Sigma.Y)
    Ux = M1$vectors
    Dx = M1$values
    Uy = M2$vectors
    Dy = M2$values
    Ham = 1/(Dx %*% t(Dy) + rho)

    Delta = matrix(0, p, p)
    Lambda = matrix(0, p, p)

    fit$path = vector("list", nlambda);
    fit$iter = rep(0, nlambda);

    start = proc.time()[3]
    for(i in 1:length(lambda)){
        iter = 0
        lambda = fit$lambda[i]
        out = .C("assymmetric_admm", Delta = as.double(Delta), Lambda = as.double(Lambda),
                 as.double(fit$Sigma.X), as.double(fit$Sigma.Y), as.double(rho), as.double(Ham),
                 as.double(Ux),as.double(Uy),as.double(stop.tol),as.integer(p),as.double(lambda), as.integer(max.iter), iter = as.integer(iter));
        Delta = matrix(out$Delta, p, p)
        fit$iter[i] = out$iter

        fit$path[[i]] = Matrix(Delta, sparse = T);
        fit$sparsity[i] = sum(Delta != 0)/p/(p-1)

    }
    fit$elapse =  proc.time()[3] - start

    rm(X, Y, Delta)
    gc()

    return(fit);

}



