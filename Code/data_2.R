rm(list = ls(all = TRUE))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(factoextra))
suppressMessages(library(plotly))
library(cluster)
library(stringr)
library(miscTools)

#-----data preprocessing----
scRNA <- read.table('GSM3405531_PDAC-B-indrop1.tsv', header=TRUE, sep="\t")
ST <- read.table('GSM3405534_PDAC-B-ST1.tsv', header = TRUE, sep = '\t', quote = "")
rownames(scRNA) <- scRNA[,1]
scRNA <- scRNA[,-1]
rownames(ST) <- ST[,1]
ST <- ST[,-1]
ST <- as.data.frame(t(ST))

# #find the same genes of two data sets(index as well)
# gen_1 <- rownames(scRNA)
# gen_2 <- rownames(ST)
# same_gen <- intersect(gen_1,gen_2)
# G = length(same_gen)
# re1 <- c()
# re2 <- c()
# for (i in 1:G) {
#   re1[i] <- which(gen_1 == same_gen[i])
#   re2[i] <- which(gen_2 == same_gen[i])
# }
# 
# scRNA_same <- scRNA[re1, ]
# ST_same <- ST[re2,]

#Data normalization
tmp = median(colSums(scRNA))/colSums(scRNA)
RNA_norm = floor(sweep(scRNA,2,tmp,'*'))
tmp = median(colSums(ST))/colSums(ST)
ST_norm = floor(sweep(ST,2,tmp,'*'))

#find the 1000 most variable genes
sd_rna = apply(RNA_norm, 1, sd)
sd_st = apply(ST_norm, 1, sd)
# RNA_1000 = RNA_norm[order(sd_rna, decreasing = T),] %>% head(1000)
# ST_1000 = ST_norm[order(sd_st, decreasing = T),] %>% head(1000)
RNA_500 = RNA_norm[order(sd_rna, decreasing = T),] %>% head(500)
ST_500 = ST_norm[order(sd_st, decreasing = T),] %>% head(500)
# write.csv(RNA_500, file = "scRNA_processed.csv")
# write.csv(ST_500, file = "ST_processed.csv")

#----heat map of both data----
# rna_1000 <- scRNA[order(sd_rna, decreasing = T),] %>% head(1000)
# st_1000 <- ST[order(sd_st, decreasing = T),] %>% head(1000)
rna_500 <- scRNA[order(sd_rna, decreasing = T),] %>% head(500)
st_500 <- ST[order(sd_st, decreasing = T),] %>% head(500)

vals <- unique(scales::rescale(c(volcano)))
o <- order(vals, decreasing = FALSE)
cols <- scales::col_numeric("Blues", domain = NULL)(vals)
colz <- setNames(data.frame(vals[o], cols[o]), NULL)
plot_ly(z = as.matrix(rna_500), zmax= 100 ,zmin=0, colorscale = colz, type = "heatmap")%>%
  layout(title = "Heatmap of scRNA-seq data",
         xaxis = list(title = "Cells"),
         yaxis = list(title = "Gene"))

cols2 <- scales::col_numeric("Reds", domain = NULL)(vals)
colz2 <- setNames(data.frame(vals[o], cols2[o]), NULL)
plot_ly(z = as.matrix(st_500), zmax= 100 ,zmin=0, colorscale = colz2, type = "heatmap")%>%
  layout(title = "Heatmap of spots data",
         xaxis = list(title = "spots"),
         yaxis = list(title = "Gene"))

#----transform spots coordinates to numeric ones----
#get the coordinates of each spot.
sp_1 <- colnames(ST_500)
sp_2 <- str_split(sp_1, "x")
index_get <- function(x){
  x <- as.numeric(x)
  x <- as.vector(x)
}
ind <- as.data.frame(t(sapply(sp_2, index_get)))
names(ind) <- c("row_ind", "col_ind")
#ind <- arrange(ind, col_ind, row_ind)
ind$col_ind <- ind$col_ind - 1
ind$row_ind <- ind$row_ind - 1
#number of rows and colums
L <- length(table(ind$row_ind))
W <- length(table(ind$col_ind))

#----null/NA process----
#NULL：31*33-996 = 27
ST_2 <- as.data.frame(t(ST_500))
ST_2$row_ind <- ind$row_ind
ST_2$col_ind <- ind$col_ind
ST_2$coor_ind <- colnames(ST_500)
ST_2 <- arrange(ST_2, col_ind, row_ind)
ind <- arrange(ind, col_ind, row_ind)
ggplot(ind, aes(col_ind, row_ind)) + geom_point(alpha = 0.6) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) 
#from plot: 16 + 11, NULL.
#rearrange ST data by the order of coordinates of spots
ST_500 <- t(ST_2[,-c(501,502,503)])
#find NA coordinates:
find_in <- function(x, coord_ind){
  tmp <- coord_ind == x
  if(sum(tmp) == 2)
    flag <- 1
  else
    flag <- 0
  return(flag)
}

NULL_find <- function(){
  Null <- NULL
  s = 0
  for(i in 1:L){
    for(j in 1:W){
      tmp <- apply(ind, 1, find_in, coord_ind = c(i, j))
      temp <- sum(tmp)
      if(temp == 0){
        Null <- rbind(Null, c(i,j))
      }
    }
  }
  return(Null)
}

#finding neighbors
neigh_ind <- function(ell, w){
  tmp <- matrix(c(ell, w-1, ell-1, w,
                  ell+1, w, ell, w+1), 4, 2, byrow=TRUE) 
  a <- rowSums((tmp <= 0)|(tmp > L)|(tmp > W))
  ind_r <- tmp[a==0,]
  ind_x <- NULL
  for(i in 1:nrow(ind_r)){
    ind_x <- c(ind_x, ind_r[i,1]+(ind_r[i,2]-1)*L)
  }#findex in x is rom smaller to bigger
  ret_list <- list(ind_r, ind_x)
  names(ret_list) <- c("ind_r", "ind_x")
  return(ret_list)
}
#null coordinate index
null_na <- as.data.frame(NULL_find())
null_na <- arrange(null_na, null_na[,2], null_na[,1])
null_na
write.csv(null_na, file = 'null_index.csv')
#plot a more clear figure
data_plus <- null_na[17:27,]
names(data_plus) <- c("row_ind", "col_ind")
tmp<- rbind(ind, data_plus)
tmp$I <- c(rep(2,996),rep(13,11))
ggplot(tmp, aes(col_ind, row_ind, color = letters[I])) + geom_point(alpha = 0.8) +
theme(axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank()) +
  geom_hline(yintercept = c(4.5,0.5), linetype="dotted") +
  geom_vline(xintercept = c(0.5,4.5), linetype="dotted") +
  scale_fill_discrete(guide=FALSE)

#get neighbor's index of null data(11 x 1)
nei_null <- list()
for(i in 17:27){
  nei_null[[i-16]] <- neigh_ind(null_na[i,1], null_na[i,2])[[1]]
}

#11个缺失值用周围的spots数值平均值填充。
find_nei_ind <- function(a){
  nr <- nrow(a)
  ind_ori <- c()
  for(i in 1:nr){
    tmp <- apply(ind, 1, find_in, coord_ind = a[i,])
    ind_ori[i] <- which(tmp == 1)
  }
  return(ind_ori)
}
null_nei_ori <- lapply(nei_null, find_nei_ind)

# continuous
RNA_c <- log(2*(RNA_500 +1))
ST_c <- log(2*(ST_500 +1))
null_com <- matrix(NA, 500, 11)
for(i in 1:length(null_nei_ori)){
  s = 0
  for(j in null_nei_ori[[i]]){
    s = s + ST_c[,j]
  }
  null_com[,i] <- s/length(null_nei_ori[[i]])
}

tmp <- ST_c
j = 1
for(i in 1:length(null_nei_ori)){
  tmp <- insertCol(tmp, null_nei_ori[[i]][2] + j , null_com[,i])
  j = j+1
}

#first 16 null/na values assumed to be 0 in all genes,
# and initilized as one cluster different from others
c1 <- log(2 * matrix(1, 4, 500))
tmp <- matrix(as.numeric(tmp), dim(tmp)[1], dim(tmp)[2])
for(i in 1:4){
  tmp <- insertCol(tmp, L*(i-1) + 1, c1)
}
dim(tmp)
ST_complete <- tmp
rownames(ST_complete) <- rownames(ST_c)
write.csv(ST_complete, file = 'ST_complete.csv')
write.csv(RNA_c, file = 'RNA_c.csv')
