# ---------------------------- generate countdata -------------------------------
rm(list = ls())
countdata = read.csv(file = "/home/johntan/Downloads/exp_seq.tsv", header = T, na.strings = "NA", sep = "\t")
countdata = countdata[, c(5,8,10)]

library(reshape2)
countdata = dcast(countdata, submitted_sample_id ~ gene_id)

row.names = as.character(countdata[,1])
col.names = colnames(countdata)[2:ncol(countdata)]
countdata = as.data.frame(countdata[, 2:ncol(countdata)], col.names = col.names, row.names = row.names)

save(countdata, file = "../comparision/cancer/countdata.rda")

rm(list = ls())
load(file = "../comparision/cancer/countdata.rda")
n = nrow(countdata)
gene.name = colnames(countdata)

## ---------------------------- delete missing value -------------------------------

cancer = as.data.frame(countdata[seq(1, n, 2),])
normal = as.data.frame(countdata[seq(2, n, 2),])


missing.count = apply(is.na(normal), 2, sum)
ind1 = missing.count / nrow(normal) <= 0.2

missing.count = apply(is.na(cancer), 2, sum)
ind2 = missing.count / nrow(cancer) <= 0.2

normal = normal[, ind1&ind2]
cancer = cancer[, ind1&ind2]

## ---------------------------- imputation -------------------------------

library(impute)
normal = t(impute.knn(t(normal))$data)
cancer = t(impute.knn(t(cancer))$data)
# round
normal = round(normal)
cancer = round(cancer)
# integer
normal = apply(normal, c(1,2), as.integer)
cancer = apply(cancer, c(1,2), as.integer)

normal = t(impute.knn(t(normal))$data)
cancer = t(impute.knn(t(cancer))$data)
# round
normal = round(normal)
cancer = round(cancer)
# integer
normal = apply(normal, c(1,2), as.integer)
cancer = apply(cancer, c(1,2), as.integer)

## ---------------------------- DESeq2 -------------------------------
library(DESeq2)
normal.T = t(normal)
cancer.T = t(cancer)

data.T = cbind.data.frame(normal.T, cancer.T)

condition = factor(c(rep("normal", ncol(normal.T)), rep("cancer", ncol(cancer.T))))

coldata = data.frame(row.names = colnames(data.T), condition)

# object construction
dds = DESeqDataSetFromMatrix(countData=data.T, colData=coldata, design=~condition)

# standard analysis
dds = DESeq(dds)

# save dds
save(dds, file = "dds.RData")
load(file = "../comparision/cancer/dds.RData")

# moderated log2 fold changes
res = results(dds)
res = res[order(res$pvalue),]

summary(res)
table(res$padj<0.05)

diff_gene_deseq2 = subset(res, abs(log2FoldChange) > 0.5)
dim(diff_gene_deseq2)

diff_gene_deseq2 = as.data.frame(diff_gene_deseq2)

diff.gene = sort(rownames(diff_gene_deseq2))

cancer = as.data.frame(cancer[, diff.gene])
normal = as.data.frame(normal[, diff.gene])

save(file = "../comparision/cancer/selected.data.rda", list = c("cancer", "normal"))

