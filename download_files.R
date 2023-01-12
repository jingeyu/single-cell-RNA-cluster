n <- 20000
GG <- matrix(rnorm(n*n),n,n)



dir.create("STATMap")
download <- function(name) {
  url <- "https://www.math.ttu.edu/~atrindad/tsdata/index.html/"
  download.file(paste0(url, name), name, quiet = TRUE)
}
download("q-gdp-ukcaus.txt")
download("pathways.rda")
download("lrpairs.rda")

for(i in 1:6){
  download(paste0("mouse2_sample",i,".txt"))
}

download <- function(name) {
  url <- "https://raw.githubusercontent.com/jingeyu/data/main/"
  download.file(paste0(url, name), name, quiet = TRUE)
}
download("insurance2.csv")
data = "https://raw.githubusercontent.com/jingeyu/data/main/insurance2.csv"
dat <- read.csv(data,header = T)
dim(dat)
dattt = read.csv("insurance2.csv",header = TRUE)