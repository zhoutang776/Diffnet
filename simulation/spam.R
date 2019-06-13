rm(list = ls())
library(diffnet)
library(ggplot2)
data(spam)


nlambda = 100
tuning = "bic"
npn = T
perturb = F
correlation = T
max.iter = 200
stop.tol = 1e-5
lambda.min.ratio = 0.1


header = names(spam)
header = gsub("word_freq_", "", header)
header = gsub("char_freq_", "", header)
header = gsub("capital_run_length_average", "average_length", header)
header = gsub("capital_run_length_total", "total_length", header)

n_X = num.spam = nrow(spam)
n_Y = num.nonspam = nrow(nonspam)
p_X = p_Y = num.feature = ncol(spam)

spam.train = spam
nonspam.train = nonspam
## ---------------------------- npn or not -----------------------------------------

if(npn){
    spam.train = diffnet.npn(spam.train)
    nonspam.train = diffnet.npn(nonspam.train)
}
## ----------------------------  lasso  -----------------------------------------

fit.lasso = diffnet(spam.train, nonspam.train, method = "lasso", nlambda=nlambda, lambda.min.ratio=lambda.min.ratio, max.iter=max.iter, stop.tol=stop.tol, perturb=perturb, correlation=correlation)
fit.lasso = diffnet.select(fit.lasso, tuning = tuning)

# ## ---------------------------- check  -----------------------------------------
print(fit.lasso$opt)
id = fit.lasso$opt[[1]][1]
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
Delta.lasso = matrix(fit.lasso$path[[opt.ind.lasso]], num.feature, num.feature)

## ---------------------------- simplified plot -----------------------------------------

simplified.plot = function(Delta, title){
    p = ncol(Delta)
    ind = apply(apply(Delta, c(1,2), as.logical), 2, sum)>2
    # ind = as.logical(apply(abs(Delta), 2, sum))
    vertex_attr_names = header[ind]
    Delta.simplified = Delta[ind, ind]
    g = graph.adjacency(as.matrix(Delta.simplified != 0) , mode="undirected", diag=FALSE) %>%
        set_vertex_attr("name", value = vertex_attr_names)
    V(g)$label.cex = 2
    layout.grid = layout.fruchterman.reingold(g)
    plot(g, sub=title, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=0,
         margin=c(0, 0.1, 0, 0.1))

    print(sum(Delta.simplified[upper.tri(Delta.simplified, diag = F)]!=0))
}


simplified.plot(Delta.lasso, "")


pdf("spam.pdf")
simplified.plot(Delta.lasso, "")
dev.off()

pdf("spam_path.pdf")
plot(p1)
dev.off()






## ---------------------------- normal plot -----------------------------------------

normal.plot = function(Delta){
    g = graph.adjacency(as.matrix(Delta != 0), mode="undirected", diag=FALSE)

    # g = graph.adjacency(Delta, mode="undirected", diag=FALSE, weighted = T)

    layout.grid = layout.fruchterman.reingold(g)
    plot(g, layout=layout.grid, edge.color='gray50',vertex.color="red", vertex.size=2.5)
    print(sum(Delta[upper.tri(Delta, diag = F)]!=0))
}
normal.plot(Delta.lasso)
# normal.plot(Delta.dtrace)
