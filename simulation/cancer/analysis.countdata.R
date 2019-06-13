rm(list = ls())
library(diffnet)

## ---------------------------- parameter setting --------------------------------------

npn = T
perturb = T
correlation = T
stop.tol = 1e-2
max.iter = 1000
tuning = "aic"

nlambda = 10
lambda.min.ratio = 0.7
lambdas = NULL

load("../comparision/cancer/selected.data.rda")


header = names(cancer)

## ---------------------------- npn or not -----------------------------------------

if(npn){
    cancer = diffnet.npn(cancer)
    normal = diffnet.npn(normal)
}


## ----------------------------  lasso  -----------------------------------------

start = proc.time()
fit.lasso = diffnet(cancer, normal, lambda.min.ratio = lambda.min.ratio, verbose = T, method = "lasso", nlambda = nlambda, stop.tol = stop.tol, perturb = perturb, correlation = correlation)
proc.time() - start


fit.lasso$lip
fit.lasso$lambdas



dyn.load("../comparision/Difdtl/src/dtl.so")
source("../comparision/Difdtl/R/Functionsnew.R")
start = proc.time()
fit.dtrace = Dpmdtl(cancer, normal, nlambda=nlambda, lambda.min.ratio=lambda.min.ratio/2, perturb=perturb, correlation = correlation)
proc.time() - start

fro.infinity = 1

fit.lasso = diffnet.select(fit.lasso, tuning = tuning)
index = lasso$opt[[1]][fro.infinity]
lasso.opt = fit.lasso$path[[index]]

fit.dtrace = diffnet.select(fit.dtrace, tuning = tuning)
index = dtrace$opt[[1]][fro.infinity]
dtrace.opt = fit.dtrace$path[[index]]


nnz(lasso.opt)
nnz(dtrace.opt)



Delta = fit.lasso$path[[1]]
sum(Delta[upper.tri(Delta, diag = F)]!=0)





# save.image(file = "fit.cancer.0.533881.RData")


## ---------------------------- select -----------------------------------------

Delta.lasso = fit.lasso$path[[1]]

threshold = stop.tol / fit.lasso$lip * 10000
Delta.lasso[abs(Delta.lasso) < threshold] = 0

ind = as.logical(apply(abs(Delta.lasso), 2, sum))
header.simplified = header[ind]
Delta.simplified = as.data.frame(Delta.lasso[ind, ind])
rownames(Delta.simplified) = header.simplified
colnames(Delta.simplified) = header.simplified




## ---------------------------- simplified plot -----------------------------------------

par(mfrow = c(1,1))
simplified.plot = function(Delta, names){
    g = graph.adjacency(as.matrix(Delta.simplified != 0) , mode="undirected", diag=FALSE) %>%
        set_vertex_attr("name", value = names) %>% set_edge_attr("color", value = "blue")

    g.neg = graph.adjacency(as.matrix(Delta.simplified < 0) , mode="undirected", diag=FALSE) %>%
        set_vertex_attr("name", value = names) %>% set_edge_attr("color", value = "red")

    # g = union(g.pos, g.neg)
    # g = graph.adjacency(Delta.simplified , mode="undirected", diag=FALSE, weighted = T) %>%
    #     set_vertex_attr("name", value = names)

    layout.grid = layout.fruchterman.reingold(g)
    plot(g, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=0, vertex.label.cex = 0.8)
    print(sum(Delta.simplified[upper.tri(Delta.simplified, diag = F)]!=0))
}

simplified.plot(Delta.simplified, header.simplified)

