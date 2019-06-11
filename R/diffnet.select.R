diffnet.select = function(fit, Sigma.X = NULL, Sigma.Y = NULL, tuning=c("none","aic","bic")){
    gcinfo(FALSE)
    fit$tuning = tuning[1]
    n.X = fit$n.X
    n.Y = fit$n.Y
    if(is.null(Sigma.X)) Sigma.X = fit$Sigma.X
    if(is.null(Sigma.Y)) Sigma.Y = fit$Sigma.Y

    res <- switch(tuning[1], ## default is "none"
                  none=NA,
                  aic=diffnet.ic(Sigma.X, Sigma.Y, fit$path, n.X+n.Y, 2),
                  bic=diffnet.ic(Sigma.X, Sigma.Y, fit$path, n.X+n.Y, log(n.X+n.Y))
                  );
    if(tuning[1] != "none"){
        fit$opt = res$opt
        fit$ic = res$ic
    }

    class(fit) = "diffnet"
    rm(n.X, n.Y)
    gc()
    return(fit)
}

diffnet.ic = function(Sigma.X, Sigma.Y, path, n, penalty){
    lowtri = which(lower.tri(path[[1]], diag=TRUE))
    df = sapply(path, function(x){ sum(x[lowtri]!=0); })
    ic = scale(n*sapply(path, lossFunction, Sigma.X, Sigma.Y), center = -penalty*df, scale = FALSE)
    opt = vector("list", 2)

    opt.ind = apply(ic, 1, which.min)
    names(opt.ind) <- c("max","F");
    opt[[1]] = opt.ind

    opt.loss= apply(ic, 1, min)
    names(opt.loss) <- c("max","F");
    opt[[2]] = opt.loss
    return(list(opt=opt, ic=ic))
}


lossFunction = function(Delta, Sigma.X, Sigma.Y){
    err = Sigma.X%*%Delta%*%Sigma.Y - Sigma.X + Sigma.Y
    return(
        c(max(abs(err)), ## max
          sqrt(sum(err^2))## frobenius
        )
    )
}
