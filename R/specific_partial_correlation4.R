if(!require(ppcor)){
  install.packages("ppcor", repos='http://lib.stat.cmu.edu/R/CRAN/')
}
library(ppcor)

args = commandArgs(trailingOnly = TRUE)

df <- na.omit(read.csv(args[1]))

res <- pcor.test(x = df$x, y = df$y, z = matrix(c(df$z1, df$z2, df$z3, df$z4), ncol = 4), method = args[2])

write.csv(c(res$estimate, res$p.value), args[3], row.names=FALSE)