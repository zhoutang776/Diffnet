rm(list = ls())
library(diffnet)
library(ggplot2)
source(file = "../comparision/refinement.R")
iter = 1
ii = 1

table2 = data.frame(matrix(0, iter, 5))
colnames(table2) = c("fro", "infinity", "size", "rate", "time")
lasso.table = table2


## ---------------------------- preprocess -----------------------------------------
n_X = 200; n_Y = n_X;
p_X = 400; p_Y = p_X;
case = "case2"

nlambda = 50
stop.tol = 1e-5
perturb = F
correlation = F
lambda.min.ratio = 0.5
max.iter = 50
savefile = F


# ------------ data generating -----------------
data = diffnet.case(n=n_X, p=p_X, method = case)

X = data$X
Y = data$Y
diff.Omega = data$diff.Omega

print(sum(diff.Omega!=0))
select = function(fit){
    return(diffnet.select(fit, data$Sigma.X, data$Sigma.Y, tuning))
}

lasso = NULL
# ---------------------------------- lasso  ---------------------------------

start = proc.time()[3]
lasso = diffnet(X, Y, verbose=F, nlambda=nlambda, max.iter=max.iter, lambda.min.ratio=lambda.min.ratio, stop.tol=stop.tol, method = "lasso", perturb=perturb, correlation=correlation)
lasso.table[ii,5] = proc.time()[3] - start
lasso = select(lasso)
id = max(lasso$opt[[1]][1], lasso$opt[[1]][2])
res = compute.table(lasso$path[[id]], diff.Omega)
lasso.table[ii,1:4] = res
print("lasso done!")
print(lasso.table[ii,])


apply(lasso.table, 2, mean)
apply(lasso.table, 2, sd)


## ----------------------------- compute solution path ------------------------------------
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



#------------------------------ save something ------------------------------------
if(savefile){
    for(i in 1:nlambda){
        lasso$path[[i]] = Matrix(lasso$path[[i]], sparse = T)
    }
    fileaddress = "/home/johntan/Desktop"
    filename = paste(case, paste("n", n_X, sep=""), paste("p", p_X, sep=""), sep="_")
    filepath = file.path(fileaddress, paste(filename, "RData", sep = "."))

    save(list = c("lasso", "lasso.nnz.path", "data", "nlambda",
                  "max.iter", "lambda.min.ratio", "stop.tol", "perturb", "correlation"), file = filepath)
}


#------------------------------- plot  --------------------------------------
plot.lassoPath <- function(coeffs, lambdas, title=NULL, scipen=999) {
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
        ylim(-1,2)+
        theme(panel.background = element_rect(fill = "transparent",colour = NA),
              legend.position="None", legend.title = element_blank(),
              panel.border = element_rect(linetype = "solid", fill = NA))
    return(g)
}
solution.nnz.path = lasso.nnz.path
p1 = plot.lassoPath(solution.nnz.path, lasso$lambdas, "")
plot(p1)




## ---------------------------------------------------- load file ----------------------------------------

# ---------------- first fig ----------------

n_X = 400
p_X = 800
case = "case2"
fileaddress = "/home/johntan/Desktop"
filename = paste(case, paste("n", n_X, sep=""), paste("p", p_X, sep=""), sep="_")
filepath = file.path(fileaddress, paste(filename, "RData", sep = "."))
load(file = filepath)



solution.nnz.path = lasso.nnz.path
p1 = plot.lassoPath(solution.nnz.path, lasso$lambdas, "")

pdf(paste(filename, ".pdf", sep = ""))
plot(p1)
dev.off()
