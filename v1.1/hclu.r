#Programmed by guofeifei, used for hierarchical clustering
#usage: Rscript hclu.r simMatrix-test.txt clu-test.pdf
library(ape)
params<-commandArgs(trailingOnly=TRUE)
inputfile=params[1]
outfile=params[2] #输出文件
mat<-read.table(inputfile) #read data matrix
hc <- hclust(as.dist(mat),method="complete") # saperate data into xxx cluster
#hcd = as.dendrogram(hc)
hcd = as.phylo(hc)
ad=paste(getwd(),outfile,sep='/')
pdf(ad,width=8, height=8)#设置画布大小
plot(hcd)
dev.off()