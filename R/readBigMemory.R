library(data.table)
library(bigmemory)

x <- read.big.matrix(filename = "../plink.tped", sep=" ", type = "char")


x2 <- sapply(1:nrow(x), function(row.ind) {
  a <- cps:::recodeAD(x[row.ind,-(1:4)])
  c(a[1,], a[2,])
})

dim(x)

x[1,]
