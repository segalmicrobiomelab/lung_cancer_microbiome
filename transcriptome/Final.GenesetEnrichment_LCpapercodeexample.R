## Michelle Badri
## Gene Set Enrichment Analysis using fgsea package 


library(fgsea)

## Dummy example

# # create .rnk file from LogFoldChange .csv files # Dataset1
# x <- read.table("DATASET1_foldchange.csv",sep=",",header=T)
# colnames(x) <- c("Name","basemean","logFC", "log2FCunshrunk","pvalue","padj")
# x <- x[!is.na(x$padj),]
# x <- x[x$padj <= 0.1,] # set pdj threshold
# x$fcSign <- sign(x$logFC)
# x$logP   <- -log10(x$padj)
# x$metric <- x$logP/x$fcSign
# y <- x[,c("Name", "metric")]
# write.table(y,file="DATASET1.rnk",quote=F,sep="\t",row.names=F) # write .rnk file
# 
# # create gmt list from LogFoldChange .csv files # Dataset2
# a<-read.table("DATASET2_foldchange.csv",sep=",",header=T)
# colnames(a) <- c("Name","basemean","logFC", "log2FCunshrunk","pvalue","padj")
# a <- a[!is.na(a$padj),]
# a <- a[a$padj <= 0.1,] # set pdj threshold
# a$fcSign <- sign(a$logFC)
# gmt.file<-c()
# gmt.file$GenesetUp<-a[a$fcSign==1,]$Name
# gmt.file$GenesetDown<-a[a$fcSign==-1,]$Name
# 
# 
# ## Compare .rnk file to GMT list 
# rnk.file <- "DATASET1.rnk"
# ranks    <- read.table(rnk.file,header=TRUE, colClasses = c("character", "numeric"))
# ranks    <- setNames(ranks$metric,ranks$Name)
# ranks    <- ranks[!is.infinite(ranks)] 
# #run fgsea 
# fgseaRes <- fgsea(gmt.file, ranks, nperm=1000)
# fgseaRes





### Working example comparison from LC paper

# create .rnk file from LogFoldChange .csv files # Dataset1
x <- read.table("A549 dilution buffer vs CSC + sup-Veil.csv",sep=",",header=T)
colnames(x) <- c("Name","basemean","logFC", "log2FCunshrunk","pvalue","padj")
x <- x[!is.na(x$padj),]
x <- x[x$padj <= 0.1,] # set pdj threshold
x$fcSign <- sign(x$logFC)
x$logP   <- -log10(x$padj)
x$metric <- x$logP/x$fcSign
y <- x[,c("Name", "metric")]
write.table(y,file="dilutionbuff_vs_CSCsupVeil.rnk",quote=F,sep="\t",row.names=F)

# create gmt list from LogFoldChange .csv files # Dataset2
a<-read.table("A549 dilution buffer vs CSC + Veil.csv",sep=",",header=T)
colnames(a) <- c("Name","basemean","logFC", "log2FCunshrunk","pvalue","padj")
a <- a[!is.na(a$padj),]
a <- a[a$padj <= 0.1,] # set pdj threshold
a$fcSign <- sign(a$logFC)
gmt.file<-c()
gmt.file$GenesetUp<-a[a$fcSign==1,]$Name
gmt.file$GenesetDown<-a[a$fcSign==-1,]$Name


## Compare .rnk file to GMT list 
rnk.file <- "dilutionbuff_vs_CSCsupVeil.rnk"
ranks    <- read.table(rnk.file,header=TRUE, colClasses = c("character", "numeric"))
ranks    <- setNames(ranks$metric,ranks$Name)
ranks    <- ranks[!is.infinite(ranks)] 
#run fgsea 
fgseaRes <- fgsea(gmt.file, ranks, nperm=1000)
fgseaRes


