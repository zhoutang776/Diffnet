\name{diffnet}
\alias{diffnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Differential Network
}
\description{
A Fast Iterative Algorithm for High-dimensional Differential Network
}
\usage{
diffnet = function(X, Y, lambdas=NULL, lambda.min.ratio = NULL, nlambda = NULL, a = NULL,
                         method = c("lasso", "scad", "mcp"), tuning=c("none","aic","bic"),
                         perturb = FALSE, stop.tol=1e-5, max.iter=800,
                         correlation=FALSE, verbose=FALSE, Delta.init = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{design matrix of group 1}
  \item{Y}{design matrix of group 2}
  \item{lambda}{regularization parameter, which control the sparisty of the solution}
  \item{lambda.min.ratio}{it is the smallest value for lambda, as a fraction of the uppperbound (MAX) of the regularization/thresholding parameter which makes all estimates equal to 0. The program can automatically generate lambda as a sequence of length = nlambda starting from MAX to lambda.min.ratio*MAX in log scale.}
  \item{nlambda}{The number of regularization/thresholding paramters. The default value is 10}
  \item{method}{Graph estimation methods with 3 penalty: "lasso", "mcp" and "scad". The defaulty value is "lasso".}
  \item{tuning}{Model selection criterion.}
  \item{perturb}{if perturb=TRUE, then perturb the sample covariance matrix to make it definitive positive}
  \item{stop.tol}{the stop condtion of the algorithm}
  \item{max.iter}{the maximum of the iteration}
  \item{correlation}{standarize the variables or not}
  \item{verbose}{if verbose equals TRUE, it will return the covergent information}
  \item{Delta.init}{the initial solution of the algorithm.}
}

\value{
An object with S3 class "diffnet" is returned:
\item{method}{
    The penalty of graph estimation
}
\item{n.X}{
the number of group 1
}
\item{n.Y}{
the number of group 2
}
  \item{Sigma.X}{
The sample covariance of group 1
}
  \item{Sigma.Y}{
The sample covariance of group 2
}
    \item{lip}{
The lipschitz constant of the problem
}

    \item{lambdas}{
The lambda of solution path
}
  \item{path}{
The solution path
}
  \item{sparsity}{
The sparsity of each solution
}
  \item{iter}{
The number of iterations of each solution
}
  \item{tuning}{
The model selection criterion
}
  \item{opt}{
The index of optimal solution
}
  \item{ic}{
The score of each solution under "aic" or "bic" selection criterion
}

}

\references{
Tang, Zhou, Zhangsheng Yu, and Cheng Wang. "A fast iterative algorithm for high-dimensional differential network." Computational Statistics 35, no. 1 (2020): 95-109.
}
\author{
Zhou Tang, Zhangsheng Yu, Cheng Wang
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
rm(list = ls())
library(diffnet)
set.seed(123)
library('MASS')
library('Matrix')
library('diffnet')
## ---------------------------- preprocess ------------------------------------
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
result = diffnet(
    X, Y, verbose = FALSE, nlambda = nlambda, max.iter = max.iter,
    lambda.min.ratio = lambda.min.ratio, stop.tol=stop.tol, method = "lasso",
    perturb =perturb, correlation = correlation
)
elapse = proc.time()[3] - start

print(result$path[[50]][1:5,1:5])

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~Differential }
\keyword{ ~Gaussian Graphical Model }% __ONLY ONE__ keyword per line
