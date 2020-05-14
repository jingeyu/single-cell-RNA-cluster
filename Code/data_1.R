rm(list = ls(all = TRUE))
library(ggplot2)
library(dplyr)
library(factoextra)
library(plotly)
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
rna_varible = scRNA[order(sd_rna, decreasing = T),] %>% head(100)
st_varible = t(ST[,order(sd_st, decreasing = T)]) %>% head(1000)
max(rna_varible) #723
max(st_varible) #5


vals <- unique(scales::rescale(c(volcano)))
o <- order(vals, decreasing = FALSE)
cols <- scales::col_numeric("Blues", domain = NULL)(vals)
colz <- setNames(data.frame(vals[o], cols[o]), NULL)
plot_ly(z = as.matrix(rna_varible)[,1:100], zmax= 500 ,zmin=0, colorscale = colz, type = "heatmap")%>%
  layout(title = "Heatmap of scRNA-seq data",
         xaxis = list(title = "Cells"),
         yaxis = list(title = "Gene"))

cols2 <- scales::col_numeric("Reds", domain = NULL)(vals)
colz2 <- setNames(data.frame(vals[o], cols2[o]), NULL)
plot_ly(z = as.matrix(st_varible), zmax= 10 ,zmin=0, colorscale = colz2, type = "heatmap")%>%
  layout(title = "Heatmap of spots data",
         xaxis = list(title = "spots"),
         yaxis = list(title = "Gene"))

set.seed(1996)
df_1 <- t(RNA[order(sd_rna, decreasing = T),] %>% head(100))
fviz_nbclust(df_1, kmeans, method = "wss",k.max=10) + geom_vline(xintercept = 7, linetype = 2)
km_result <- kmeans(df_1, 7)
fviz_cluster(km_result, data = df_1,
             ellipse = F,
             show.clust.cent = F,
             shape = 1,
             repel = F,
             ggtheme = theme_grey()
)


df_2 <- ST_n[,order(sd_st, decreasing = T)][,1:100]
fviz_nbclust(df_2, kmeans, method = "wss",k.max=20) 
km_result <- kmeans(df_2, 6)
fviz_cluster(km_result, data = df_2,
             ellipse = F,
             show.clust.cent = F,
             shape = 1,
             repel = F,
             ggtheme = theme_grey()
)


res1 = get_clust_tendency(df_1, 50, graph = TRUE, gradient = list(low = "steelblue", high = "white"))
res1$hopkins_stat
# 0.9073394
res1$plot

res2 = get_clust_tendency(df_2, 50, graph = TRUE, gradient = list(low = "steelblue", high = "white"))
res2$hopkins_stat
# 0.9831106
res2$plot




suppressMessages(library(factoextra))
df<-t(RNA)
fviz_nbclust(df, kmeans, method = "wss",k.max=10) + geom_vline(xintercept = 6, linetype = 2)

#continuous 
# RNAC <- log(2*(RNA +1))
# STC <- log(2*(ST_n +1))





  
  
  
  