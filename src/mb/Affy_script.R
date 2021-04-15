library(affy)
data <- ReadAffy()
eset <- rma(data)
write.exprs(eset,file="data.txt")