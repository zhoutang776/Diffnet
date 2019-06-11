# # ---------------------------- generate countdata -------------------------------
# rm(list = ls())
# countdata = read.csv(file = "/home/johntan/Downloads/exp_seq.tsv", header = T, na.strings = "NA", sep = "\t")
# countdata = countdata[, c(5,8,9)]
#
# library(reshape2)
# countdata = dcast(countdata, submitted_sample_id ~ gene_id)
#
# row.names = as.character(countdata[,1])
# col.names = colnames(countdata)[2:ncol(countdata)]
# countdata = as.data.frame(countdata[, 2:ncol(countdata)], col.names = col.names, row.names = row.names)
#
# save(countdata, file = "../comparision/cancer/countdata.rda")

rm(list = ls())
load(file = "../comparision/cancer/countdata.rda")

n = nrow(countdata)
header = names(countdata)
## ---------------------------- choosing pathway -----------------------------------------

kegg2gene = read.csv(file = "../comparision/cancer/kegg2gene.txt", sep = '\t', header = F)
# 05200  Pathways in cancer
# 05202  Transcriptional misregulation in cancer
# 05203  Viral carcinogenesis
# 05204  Chemical carcinogenesis
# 05205  Proteoglycans in cancer
# 05206  MicroRNAs in cancer
# 05230  Central carbon metabolism in cancer
# 05231  Choline metabolism in cancer

# 05225  Hepatocellular carcinoma


pathway = c(05200, 05202, 05203, 05204, 05205, 05206, 05230, 05231, 05225)
ind = NULL
for(i in 1:length(pathway)){
    tmp = which(kegg2gene[,1] == pathway[i])
    ind = union(ind, tmp)
}
sub.EntrezID = unique(kegg2gene[ind,2])
length(sub.EntrezID)

## ---------------------------- convert entrez to symbol name -------------------------------------

library(annotables)
symbol = rep(0, length(sub.EntrezID))
for(i in 1:length(sub.EntrezID)){
    ind = which(grch38[,2] == sub.EntrezID[i])
    tmp = unlist(grch38[ind, 3])
    if(length(tmp) > 1){
        for(j in 1:length(tmp)){
            if(which(header == tmp[j])){
                tmp = tmp[j]
                break
            }
        }
    }
    if(length(tmp) == 1 && is.element(tmp, header)) symbol[i] = tmp
}
symbol = symbol[-which(symbol == "0")]

## ---------------------------- convert entrez to symbol name -------------------------------------
cancer = as.data.frame(countdata[seq(1, n, 2), symbol])
colnames(cancer) = symbol
normal = as.data.frame(countdata[seq(2, n, 2), symbol])
colnames(normal) = symbol
## ---------------------------- imputation -------------------------------

library(impute)
normal = t(impute.knn(t(normal))$data)
cancer = t(impute.knn(t(cancer))$data)





## ---------------------------- save -------------------------------

normal = as.data.frame(normal)
cancer = as.data.frame(cancer)

# save(file = "../comparision/cancer/selected.data.rda", list = c("cancer", "normal"))

save(file = "../comparision/cancer/two.pathway.rda", list = c("cancer", "normal"))

# ## ---------------------------- delete missing value -------------------------------
#
# missing.count = apply(is.na(normal), 2, sum)
# ind1 = missing.count / nrow(normal) <= 0.2
#
# missing.count = apply(is.na(cancer), 2, sum)
# ind2 = missing.count / nrow(cancer) <= 0.2
#
# normal = normal[, ind1&ind2]
# cancer = cancer[, ind1&ind2]



