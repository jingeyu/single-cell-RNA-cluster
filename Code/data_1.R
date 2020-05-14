#rm(list = ls(all = TRUE))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(factoextra))
suppressMessages(library(plotly))
library(cluster)

scRNA <- read.table('GSM3405531_PDAC-B-indrop1.tsv', header=TRUE, sep="\t")
ST <- read.table('GSM3405534_PDAC-B-ST1.tsv', header = TRUE, sep = '\t', quote = "")
rownames(scRNA) <- scRNA[,1]
scRNA <- scRNA[,-1]
rownames(ST) <- ST[,1]
ST <- ST[,-1]
dim(ST)
dim(scRNA)

head(ST[,1:5])
head(scRNA[,1:5])

t_f <- apply(scRNA != 0, 2, sum)
index <- which(t_f == 0)
length(index)

# index = 0, all gens are expressed.
#zero-inflated effects
effect = colSums(scRNA == 0)/dim(scRNA)[1]
head(effect)

effect = rowSums(ST == 0)/dim(ST)[2]
head(effect)

#表达值在100以上的
sum(ST>10)
sum(ST > 100) #232
sum(ST > 200) #60
sum(ST > 500) #7
sum(ST > 1000) #0

#Data normalization
tmp = median(colSums(scRNA))/colSums(scRNA)
RNA = floor(sweep(scRNA,2,tmp,'*'))
tmp = median(rowSums(ST))/rowSums(ST)
ST_n = floor(sweep(ST,1,tmp,'*'))

#Heat map.
sd_rna = apply(RNA, 1, sd)
sd_st = apply(ST_n, 2, sd)
rna_varible = scRNA[order(sd_rna, decreasing = T),] %>% head(500)
st_varible = t(ST[,order(sd_st, decreasing = T)]) %>% head(500)
max(rna_varible) #723
max(st_varible) #775


res1 = get_clust_tendency(df_1, 50, graph = TRUE, gradient = list(low = "steelblue", high = "white"))
res1$hopkins_stat
# 0.9073394

res2 = get_clust_tendency(df_2, 50, graph = TRUE, gradient = list(low = "steelblue", high = "white"))
res2$hopkins_stat
# 0.9044891

#continuous 
# RNAC <- log(2*(RNA +1))
# STC <- log(2*(ST_n +1))





  
  
  
  