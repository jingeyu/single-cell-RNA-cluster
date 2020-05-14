
# ------------------------------------------------------------
# Single cell RNA clustering
# Yu jinge
# Thu May 14 11:27:25 2020
# ------------------------------------------------------------


# ------------------------------------------------------------
# Reproducibility: Reproducible_Single cell RNA clustering
# ------------------------------------------------------------

# This file contains instructionss for reproducing data,
# figures, tables and all analysis in the report.
# As denoted below, Steps 3 and 4 are computationally 
# intensive and take a *very* long time to run. 

# ------------------------------------------------------------
# Step 1: Exploratory Data Analysis
# ------------------------------------------------------------
source('~/Desktop/Code/data_1.R', echo=TRUE)
#heat maps of scRNA data and spots data
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
plot_ly(z = as.matrix(st_varible), zmax= 100 ,zmin=0, colorscale = colz2, type = "heatmap")%>%
  layout(title = "Heatmap of spots data",
         xaxis = list(title = "spots"),
         yaxis = list(title = "Gene"))

#kmeans cluster of both data
set.seed(1996)
df_1 <- t(RNA[order(sd_rna, decreasing = T),] %>% head(500))
fviz_nbclust(df_1, kmeans, method = "wss",k.max=10) + geom_vline(xintercept = 7, linetype = 2)
km_result <- kmeans(df_1, 7)
fviz_cluster(km_result, data = df_1,
             ellipse = F,
             show.clust.cent = F,
             shape = 1,
             repel = F,
             ggtheme = theme_grey()
)

df_2 <- ST_n[,order(sd_st, decreasing = T)][,1:500]
fviz_nbclust(df_2, kmeans, method = "wss",k.max=20) 
km_result <- kmeans(df_2, 6)
fviz_cluster(km_result, data = df_2,
             ellipse = F,
             show.clust.cent = F,
             shape = 1,
             repel = F,
             ggtheme = theme_grey()
)

#Dissimilarity matrix of scRNA data
res1$plot
res2$plot

# ------------------------------------------------------------
# Step 2: Data processing 
# ------------------------------------------------------------
source('~/Desktop/Code/data_2.R', echo=TRUE)


# ------------------------------------------------------------
# Step 3: Single cell RNA data cluster
# ------------------------------------------------------------
Rcpp::sourceCpp('~/Desktop/Code/sim_pj.cpp')
source('~/Desktop/Code/data_3.R', echo=TRUE)

# ------------------------------------------------------------
# Step 4: Spot data cluster
# ------------------------------------------------------------
source('~/Desktop/Code/data_4.R', echo=TRUE)
