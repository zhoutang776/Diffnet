L1_dts <- function(SigmaX, SigmaY, rho, lambda, Delta0 = NULL, Lambda0 = NULL, Ux = NULL, Dx = NULL, Uy = NULL, Dy = NULL, C1 = NULL, C2 = NULL)
{
    if(is.null(Delta0)) Delta0<-solve(SigmaY+diag(nrow(SigmaY)))-solve(SigmaX+diag(nrow(SigmaX)))
    if(is.null(Lambda0)) Lambda0<-matrix(1,nrow(SigmaX),ncol(SigmaX))

    tol = 1e-3
    p = dim(Delta0)[1]

    Delta1 = Delta0
    Delta2 = Delta0
    Delta3 = Delta0
    Lambda1 = Lambda0
    Lambda2 = Lambda0
    Lambda3 = Lambda0


    #while (k<10000){
        #z1 = kro(SigmaX-SigmaY+rho*Delta2+rho*Delta3+Lambda3-Lambda1,C1,Ux,Uy)
        #z2 = kro(SigmaX-SigmaY+rho*z1+rho*Delta3+Lambda1-Lambda2,C2,Uy,Ux)
        #z3 = S((Lambda2/rho-Lambda3/rho+z1+z2)/2,(lambda/rho)/2)
        #if (F_dis(Delta1,z1)<tol&&F_dis(Delta2,z2)<tol&&F_dis(z1,z2)<tol)
        #{
         #   break
        #}
        #Delta1 = z1
        #Delta2 = z2
        #Delta3 = z3
        #Lambda1 = Lambda1+rho*(Delta1-Delta2)
        #Lambda2 = Lambda2+rho*(Delta2-Delta3)
        #Lambda3 = Lambda3+rho*(Delta3-Delta1)
        #k = k+1
    #}
    l=p^2
    iter = 0
    result = .C("symmetric_admm", as.double(Delta0),Delta3=as.double(Delta3), as.double(Lambda0), Lambda3=as.double(Lambda0), as.double(SigmaX),
        as.double(SigmaY), as.integer(rho), as.double(C1),as.double(C2),
        as.double(Ux),as.double(Uy),as.double(tol/1000),as.integer(p),as.double(lambda), max.iter = as.integer(1000), iter=as.integer(iter));
    Delta3<-matrix(result$Delta3, nrow=p,ncol=p) 
    Lambda3 = result$Lambda3
    
    # print(Delta3[1:5,1:5])
    iter = result$iter

    return (list(Delta3, Lambda3, iter))
}

symmetric_admm<- function(X,Y, lambda=NULL,nlambda=10,lambda.min.ratio=NULL,
                rho=NULL, correlation=FALSE, perturb=FALSE,
                tuning=c("none","aic","bic","nbic"), stop.tol=1e-4){

    n.X <- nrow(X); n.Y <- nrow(Y);
    ## the number of parameters is p(p+1)/2
    p <- ncol(X);

    if(correlation){ Sigma.X <- cor(X); Sigma.Y <- cor(Y); }
    else{ Sigma.X <- cov(X)*(1-1/n.X); Sigma.Y <- cov(Y)*(1-1/n.Y); }
    if(perturb){
        cat("perturbing...")
        eigvals.X = eigen(Sigma.X, only.values=TRUE)$values;
        eigvals.Y = eigen(Sigma.Y, only.values=TRUE)$values;
        ## same perturbation as the clime software
        epsilon.X = max(max(eigvals.X)-p*min(eigvals.X), 0)/(p-1);
        epsilon.Y = max(max(eigvals.Y)-p*min(eigvals.Y), 0)/(p-1);
        Sigma.X = Sigma.X + epsilon.X*diag(p)
        Sigma.Y = Sigma.Y + epsilon.Y*diag(p)
        cat("done\n")
    }

    if(!is.null(lambda)){ nlambda <- length(lambda); }
    if(is.null(lambda))
        {
            if(is.null(lambda.min.ratio)){ lambda.min.ratio <- 0.04; }
            lambda.max <- 2*max(abs(Sigma.X-Sigma.Y));
            lambda.min <- lambda.min.ratio*lambda.max;
            lambda <- exp(seq(log(lambda.max), log(lambda.min), length=nlambda));

        }

    if(is.null(rho)){ rho <- 1; }

    M1 = eigen(Sigma.X)
    M2 = eigen(Sigma.Y)
    Ux = M1$vectors
    Dx = M1$values
    Uy = M2$vectors
    Dy = M2$values
    C1 = matrix(0,p,p)
    C2 = matrix(0,p,p)
    for (i in 1:p){
        for (j in 1:p){
            C1[i,j] = 1/(Dy[j]*Dx[i]+2*rho)
            C2[i,j] = 1/(Dy[i]*Dx[j]+2*rho)
        }
    }
    Delta0 <-solve(Sigma.Y+diag(nrow(Sigma.Y)))-solve(Sigma.X+diag(nrow(Sigma.X)))
    Lambda0 <-matrix(1,nrow(Sigma.X),ncol(Sigma.X))
    ret = vector("list", nlambda);
    iter = rep(0, nlambda);

    for(i in 1:nlambda)
    {
        ret[[i]] <- matrix(NA,nrow=p,ncol=p);
        result <- L1_dts(Sigma.X,Sigma.Y,rho,lambda[i], Delta0, Lambda0, Ux, Dx, Uy, Dy, C1, C2);
        ret[[i]] = result[[1]]
        Delta0 = result[[1]]
        Lambda0 = result[[2]]
        iter[i] = result[[3]]
    }



    return(list(path=ret,Sigma.X=Sigma.X, Sigma.Y=Sigma.Y, n.X = n.X, n.Y = n.Y, lambda=lambda,nlambda=nlambda, iter=iter));

}
#