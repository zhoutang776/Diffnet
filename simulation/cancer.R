rm(list = ls())
library(diffnet)
library(ggplot2)
load(file = "../comparision/cancer/two.pathway.rda")

header = colnames(cancer)
p = ncol(cancer)
## ---------------------------- setting -------------------------------------

nlambda = 100
lambda.min.ratio = 0.3
max.iter = 1000
stop.tol = 1e-5
npn = T
perturb = F
correlation = T
tuning = "aic"
p_X = p_Y = ncol(cancer)
n_X = nrow(cancer)
n_Y = nrow(normal)
## ---------------------------- npn or not -----------------------------------------
if(npn){
    cancer = diffnet.npn(cancer)
    normal = diffnet.npn(normal)
}

## ---------------------------- apg -------------------------------------
start = proc.time()
fit.lasso = diffnet(cancer, normal, verbose=T, nlambda=nlambda, lambda.min.ratio=lambda.min.ratio, stop.tol=stop.tol, max.iter=max.iter, method="lasso", perturb=perturb, correlation=correlation)
proc.time() - start
fit.lasso = diffnet.select(fit.lasso, tuning = tuning)

# ## ---------------------------- check  -----------------------------------------
print(fit.lasso$opt)
id = fit.lasso$opt[[1]][2]
print(nnzero(fit.lasso$path[[id]]))


## ----------------------------- compute solution path ------------------------------------
lasso = fit.lasso
solution = matrix(lasso$path[[50]], p_X, p_X)
lowtri = lower.tri(solution, diag = T)
nnz = nnzero(solution[lowtri])
solution.path = as.data.frame(matrix(0, nnz, nlambda))
colnames(solution.path) = lasso$lambdas
lasso.nnz.path = solution.path
rownames(lasso.nnz.path) = sapply(c(1:nnz), function(x){paste("lasso",as.character(x),sep = "_")})

ind = (solution != 0) & lowtri
for(i in 1:nlambda){
    res = matrix(lasso$path[[i]], p_X, p_X)
    lasso.nnz.path[,i] = res[ind]
}

plot.lassoPath <- function(coeffs, lambdas, opt.lambda, title=NULL, scipen=999) {
    options(scipen=scipen)
    coeffs = as.data.frame(coeffs)
    coeffs$vname <- rownames(coeffs)
    lambda.max = lambdas[1]
    melt.coeffs <- reshape2::melt(coeffs, id.var="vname")
    melt.coeffs$variable = as.numeric(levels(melt.coeffs$variable))[melt.coeffs$variable]
    colour = rep("black", nrow(coeffs))
    g <- ggplot(melt.coeffs, aes(x=variable, y=value, col=vname, group=vname)) + geom_line() +
        scale_x_continuous(breaks=round(seq(0, ceiling(lambda.max), length.out=10), 2)) +
        scale_colour_manual(values = colour) +
        labs(x="lambda", y="Coefficients", subtitle = title) +
        theme(panel.background = element_rect(fill = "transparent",colour = NA),
              legend.position="None", legend.title = element_blank(),
              panel.border = element_rect(linetype = "solid", fill = NA)) +
        geom_vline(xintercept=opt.lambda, colour="red")
    return(g)
}

p1 = plot.lassoPath(lasso.nnz.path, lasso$lambdas, lasso$lambdas[id])

plot(p1)


## ---------------------------- select -----------------------------------------
ic = 1 # infinty
# ic = 2 # fro

opt.ind.lasso = fit.lasso$opt[[1]][ic]
Delta.lasso = matrix(fit.lasso$path[[opt.ind.lasso]], p_X, p_X)

## ---------------------------- simplified plot -----------------------------------------

simplified.plot = function(Delta, header){
    ind = as.logical(apply(abs(Delta), 2, sum))
    header = header[ind]
    Delta = Delta[ind, ind]

    g = graph.adjacency(as.matrix(Delta != 0) , mode="undirected", diag=FALSE) %>%
        set_vertex_attr("name", value = header)

    V(g)$label.cex = 1.5

    layout.grid = layout.fruchterman.reingold(g)
    plot(g, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=0,
         margin=c(0, 0, 0, 0))
    print(sum(Delta[upper.tri(Delta, diag = F)]!=0))
}

## ------------ simplified  ----------------


ind = as.logical(apply(abs(Delta.lasso), 2, sum))
header.simplified = header[ind]
Delta.simplified = as.data.frame(Delta.lasso[ind, ind])
rownames(Delta.simplified) = header.simplified
colnames(Delta.simplified) = header.simplified

simplified.plot(Delta.simplified, header.simplified)



pdf("cancer_path_npn_correlation_iter1000.pdf")
plot(p1)
dev.off()

pdf("liver_cancer_npn_correlation_iter1000.pdf")
simplified.plot(Delta.simplified, header.simplified)
dev.off()

