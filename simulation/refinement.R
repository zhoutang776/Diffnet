symmetric_refinement = function(Delta, Sigma.X, Sigma.Y){
    if(!isSymmetric(Delta) || !isSymmetric(Sigma.X) || !isSymmetric(Sigma.Y)){
        print("somethiong wrong in here!")
    }
    p = nrow(Delta)
    index = apply(abs(Delta), 1, sum)
    index = as.logical(index)
    small.p = sum(index)

    if(small.p > 33) print("maybe it is a little dense!")
    if(small.p == 0) return(Delta)

    # first deduction
    small.Sigma.X = Sigma.X[index, index]
    small.Sigma.Y = Sigma.Y[index, index]
    small.Delta = Delta[index, index]
    A = (kronecker(small.Sigma.Y, small.Sigma.X) + kronecker(small.Sigma.X, small.Sigma.Y))/2
    x = as.double(small.Delta)
    b = as.double(small.Sigma.X - small.Sigma.Y)
    if(small.p == 1) {
        Delta[index, index] = matrix(b/A, 1, 1)
        return(Delta)
    }
    # now we want to have A * x = b, since we have already now some part of x must equal 0, hence we have second deduction

    # second deduction
    small.index = as.logical(x)
    small.A = A[small.index, small.index]
    small.b = b[small.index]

    ans.x = x
    ans.small.x = solve(small.A) %*% small.b
    ans.x[small.index] = ans.small.x
    small.Delta = matrix(ans.x, nrow(small.Delta), nrow(small.Delta))

    ans = matrix(0, p, p)
    ans[index, index] = small.Delta
    return(ans)
}


compute.table = function(Delta, diff.Omega){
    size = nnzero(Delta)
    rate = sum(Delta & diff.Omega)/sum(as.logical(diff.Omega)) * 100
    fro.loss = norm(Delta - diff.Omega, "F")
    infity.loss = max(abs(Delta - diff.Omega))
    return(c(fro.loss, infity.loss, size, rate))
}
