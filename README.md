# Introduction for R package Diffnet
Differential network is an important tool to capture the changes of conditional correlations under two sample cases. We develop an efficient fista algorithm via the symmetric quadratic loss for differential matrix estimation. The computation complexity of our algorithm is linear in the sample size and the number of parameters, which is optimal in the sense that it is of the same order as computing two sample covariance matrices.

# Getting Started
These instructions will give you a toy example for implementing the package.

## Prerequisites
What things you need to install the software and how to install them. The key functions of the package is writing in C++. So, make sure your OS can complies C++ code. For example, you should install Rtools under Windows and Xcode under MacOS. After that, the following R packages are also necessary.

    install.packages("devtools")
    install.packages("MASS")
    install.packages("Matrix")

## Install Diffnet

### Install from Github

    devtools::install_github("zhoutang776/Diffnet")
    
    
# Toy example
    library(diffnet)
    rm(list = ls())
    set.seed(123)
    library('MASS')
    library('Matrix')
    library('diffnet')
    ## ---------------------------- preprocess -----------------------------------------
    n_X = 100; n_Y = n_X;
    p_X = 100; p_Y = p_X;
    nlambda = 50
    tuning = "none"
    case = "case1"
    stop.tol = 1e-5
    perturb = FALSE
    correlation = FALSE
    max.iter = 800
    lambda.min.ratio = 0.5
    
    # ------------ data generating -----------------
    data = diffnet.case(n=n_X, p=p_X, method = case)
    X = data$X
    Y = data$Y
    diff.Omega = data$diff.Omega
    print(sum(diff.Omega!=0))
    
    start = proc.time()[3]
    result = diffnet(X, Y, verbose = FALSE, nlambda = nlambda, max.iter = max.iter, lambda.min.ratio = lambda.min.ratio, stop.tol=stop.tol, method = "lasso", perturb =perturb, correlation = correlation)
    elapse = proc.time()[3] - start
    
    print(result$path[[50]][1:5,1:5])
