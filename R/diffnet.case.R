diffnet.case = function(n, p, method){
    result = switch(method[1],
                    case1 = case1(n, p),
                    case2 = case2(n, p)
                    )
    data = list()
    data$Omega.X = result$Omega.X
    data$Omega.Y = result$Omega.Y
    data$diff.Omega = result$Omega.Y - result$Omega.X
    data$Sigma.X = solve(data$Omega.X)
    data$Sigma.Y = solve(data$Omega.Y)

    data$X = mvrnorm(n = n, mu = rep(0, p), Sigma = data$Sigma.X )
    data$Y = mvrnorm(n = n, mu = rep(0, p), Sigma = data$Sigma.Y )
    return(data)
}
generate.Delta = function(p){
    Delta=matrix(0,p,p);
    Delta[1:2,1:2]=matrix(c(0, -1, -1, 2),2)
    # Delta[1:3,1:3]=matrix(c(0, -1, 0, -1, 2, -1, 0, -1, 0), 3)
    # Delta[1:4,1:4]=matrix(c(0, -1, 0, 0, -1, 2, -1, 0, 0, -1, 2, -1, 0, 0, -1, 0), 4)
    return(Delta)
}

## sparse case
case1 = function(n, p){
    Delta = generate.Delta(p)

    Omega.X = as.matrix(band(matrix(2/3, p, p) + diag(1, p, p), -1, 1))
    Omega.X[1,1] = 4/3; Omega.X[p,p] = 4/3;
    Omega.Y = Omega.X + Delta

    case = list()
    case$Omega.X = Omega.X
    case$Omega.Y = Omega.Y
    return(case)
}


## Asymptotic sparse case
case2 = function(n, p){
    Delta = generate.Delta(p)

    Omega.X = toeplitz(0.5^(0:(p-1)))
    Omega.Y = Omega.X + Delta

    case = list()
    case$Omega.X = Omega.X
    case$Omega.Y = Omega.Y
    return(case)
}
