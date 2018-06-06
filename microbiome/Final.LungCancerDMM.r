#Load Libraries
library(DirichletMultinomial)
library(lattice)
library(xtable)
library(parallel)
library(phyloseq)

#load data table
dmm.table = read.table('otu_table.MSQ34.Airway.rare.biom.txt', header=T, sep="\t", as.is=TRUE)

#convert table to matrix
count <- t(as.matrix(dmm.table, row.names=1, col.names=1))

#fix matrix with right colnames and convert to numeric from string
colnames(count) <- count[1,]
count<-count[-1,]
class(count) <-"numeric"

#display matrix
count[1:5, 1:3]

#fit model for 1-7 clusters
fit.all <- mclapply(1:7, dmn, count=count, verbose=TRUE)

#save models to file
save(fit.all, file="fit.Lung.Cancer.Airways.rda")
save.image(file="fit.Lung.Cancer.Airways.RData")


#save all the laplace goodness of fit
lplc <- sapply(fit.all, laplace)

#plot the Lapplce based on number of clusters
pdf("Which_k_laplace.Airways.pdf", height = 10, width = 10)
plot(lplc, type="b", xlab="Number of Dirichlet Components",
     ylab="Model Fit")
dev.off()

#display the number of clusters with the lowest laplace
(best <- fit.all[[which.min(lplc)]])

mixturewt(best)

#create table with grouping of clusters
grouping <- mixture(best, assign=TRUE)
write.table(grouping, file="Lung.Cancer.Airways.grouping.txt", sep="\t")

#save models to file
save(fit.all, file="fit.Lung.Cancer.Airways.rda")
save.image(file="fit.Lung.Cancer.Airways.RData")
