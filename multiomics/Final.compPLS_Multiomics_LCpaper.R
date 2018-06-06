##Scripts for Multiomics:16S Transcriptome compPLS Association networks 
# Michelle Badri


#Load Phyloseq
library(phyloseq)
library(igraph)
library(compPLS) # code and install info github.com/zdk123/compPLS

##Load the files needed
file = "otu_table.Merged.New.biom"
map = "Mapping.lung.cancer.project.new.txt"



# Load the abundace table and mapping table 
abundance.table = import_biom(file, taxaPrefix=F)
mapping.table=sample_data(read.table(map, header=T, sep="\t", row.names=1))


lung.physeq=phyloseq(otu_table(abundance.table),tax_table(abundance.table), mapping.table)


colnames(tax_table(lung.physeq))=c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "OTU")

# Load the tree file (use the unannotated.tree)
treefile = "97_otus_unannotated.tree"
tree.obj = import_qiime(treefilename = treefile) 


# Now merge the three separate phyloseq objects into a single object
OTU.Table = merge_phyloseq(lung.physeq, mapping.table, tree.obj)


sample_data(OTU.Table)

rownames(sample_data(OTU.Table))
colnames(sample_data(OTU.Table))

# Remove taxa with 0 abundance
OTU.Table = subset_taxa(OTU.Table, rowSums(otu_table(OTU.Table)) != 0)

##If you want to nomalize OTU table before
## To normalize data you need to set a function
normalizeSample = function(x) {
  x/sum(x)
}
OTU.Rel.Table = transformSampleCounts(OTU.Table, normalizeSample)


sample_data(OTU.Table)$Sample_Type
sample_data(OTU.Table)$Description_LNS




# # Create phyllum and order tables (do it after normalization and out of the relative table)
Phylum.rel.table = tax_glom(OTU.Rel.Table, taxrank = "Phylum")
Class.rel.table = tax_glom(OTU.Rel.Table, taxrank = "Class")
Order.rel.table = tax_glom(OTU.Rel.Table, taxrank = "Order")
Family.rel.table = tax_glom(OTU.Rel.Table, taxrank = "Family")
Genus.rel.table = tax_glom(OTU.Rel.Table, taxrank = "Genus")
OTU.rel.table = tax_glom(OTU.Rel.Table, taxrank = "OTU")


Paired.Microbiome.TranscriptomeGenus.Rel.table = subset_samples(Genus.rel.table, Good_RNAseqdata =="1")


data <-data.frame(Paired.Microbiome.TranscriptomeGenus.Rel.table@tax_table@.Data)
data$id <- paste(data$Domain, data$Phylum, data$Class, data$Order,data$Family,data$Genus, sep=";")
data$id<- gsub("'","",data$id)

tab = "otu_table.MSQ34.Airway.Clus_L6_OTUgroup.csv"

# Load the clustering info and select only taxa from cluster 1
clustering.taxa=read.csv(tab,check.names=FALSE, header =1)

group1ID <- clustering.taxa[clustering.taxa$group==1,]
group1labels <- rownames(data[group1ID$OTUID,])
group1labels <- group1labels[1:32]


## Create new OTU table from group 1 taxa
otu.tab ="otutablegroup1.csv"
otu.table=read.csv(otu.tab,check.names=FALSE, header =1,row.names = 1)
otu.table<- otu_table(data.matrix(otu.table), taxa_are_rows = TRUE)
rownames(otu.table) <- group1labels

taxa <- group1ID[c(1:32),c(121:126)]
taxtab<-tax_table(taxa)
rownames(taxtab) <- group1labels

mapping.tab<-mapping.table[colnames(otu.table),]

lung.physeq=phyloseq(otu_table(otu.table, taxa_are_rows = TRUE),tax_table(taxtab), mapping.tab)
Paired.Microbiome.TranscriptomeGenus.Rel.table1B <- lung.physeq
#Paired.Microbiome.TranscriptomeGenus.Rel.table1B <- subset_samples(Paired.Microbiome.TranscriptomeGenus.Rel.table1B, Sample_Type=="Airway")
Paired.Microbiome.TranscriptomeGenus.Rel.table = subset_samples(Paired.Microbiome.TranscriptomeGenus.Rel.table1B, Good_RNAseqdata =="1")



#Paired.Microbiome.TranscriptomeGenus.Rel.table<- prune_taxa(group1labels,Paired.Microbiome.TranscriptomeGenus.Rel.table)
Paired.Microbiome.TranscriptomeGenus.Rel.wh1 = genefilter_sample(Paired.Microbiome.TranscriptomeGenus.Rel.table, filterfun_sample(function(x) x >0.002), A = 0.02 * nsamples(Paired.Microbiome.TranscriptomeGenus.Rel.table))
Paired.Microbiome.TranscriptomeGenus.Rel.table1B = prune_taxa(Paired.Microbiome.TranscriptomeGenus.Rel.wh1, Paired.Microbiome.TranscriptomeGenus.Rel.table)
colnames(sample_data(Paired.Microbiome.TranscriptomeGenus.Rel.table1B))
dim(Paired.Microbiome.TranscriptomeGenus.Rel.table1B@otu_table@.Data)
Paired.Microbiome.TranscriptomeGenus.Rel.table1B



## Need to make a table out of the phyloseq object pruned genus table.
Paired.Microbiome.Transcriptome.GenusData = otu_table(Paired.Microbiome.TranscriptomeGenus.Rel.table1B)


##Transpose the microbiome data
Paired.Microbiome.Transcriptome.GenusData.Transp = as.data.frame(t(Paired.Microbiome.Transcriptome.GenusData)) #as data frame makes it a readable matrix by R
rownames(Paired.Microbiome.Transcriptome.GenusData.Transp)  <- gsub(".64","",rownames(Paired.Microbiome.Transcriptome.GenusData.Transp))
rownames(Paired.Microbiome.Transcriptome.GenusData.Transp)  <- gsub(".65","",rownames(Paired.Microbiome.Transcriptome.GenusData.Transp))

paired.microbe.SAMdata<-data.frame(sample_data(Paired.Microbiome.TranscriptomeGenus.Rel.table1B))
rownames(paired.microbe.SAMdata)  <- gsub(".64","",rownames(paired.microbe.SAMdata))
rownames(paired.microbe.SAMdata)  <- gsub(".65","",rownames(paired.microbe.SAMdata))



table <- "transcriptome.simple.MichelleID.csv"
newlev  <-read.table(table,check.names=FALSE, header =1,sep=",")
newlev$michelleID <- as.character(newlev$michelleID)

paired.microbe.SAMdata<- sample_data(paired.microbe.SAMdata)[order(newlev$michelleID)]
Paired.Microbiome.Transcriptome.GenusData.Transp<- Paired.Microbiome.Transcriptome.GenusData.Transp[order(newlev$michelleID),]

paired.microbe.SAMdata<-cbind(paired.microbe.SAMdata, newlev$control3groups, newlev$`threegrroups-c`)

paired.microbe.TAXTAB<-as.data.frame(tax_table(Paired.Microbiome.TranscriptomeGenus.Rel.table1B))
##############################################   TRANCRIPTOME     ###################


#import the transcriptome data
Transcriptome.Data = read.csv('transcriptome.simple.Michelle.csv', row.names = 6, header = T, check.names=F)
Transcriptome.Data.Transp = as.data.frame(t(Transcriptome.Data)) #as data frame makes it a readable matrix by R


### Remove the non-sig pathways
Transcriptome.Data.Transp<-t(Transcriptome.Data[Transcriptome.Data$All_Sig==1,])
#Transcriptome.Data.Transp<-t(Transcriptome.Data[Transcriptome.Data$Pre_Knowledge_Based_Sig==1,])
#Transcriptome.Data.Transp<-t(Transcriptome.Data[Transcriptome.Data$Data_Driven_Sig==1,])

### Remove the metadata and annotation from transcription pathways 
Transcriptome.Data.Transp<- Transcriptome.Data.Transp[6:93,]

#reorder GenusData to match the Trascriptome matrix (samples need to be in row)
Transcriptome.Data.Order <- Transcriptome.Data.Transp[c(rownames(Paired.Microbiome.Transcriptome.GenusData.Transp)),]


##Transcriptome data 
dim(t(Transcriptome.Data.Order))
##OTU tab
otu_genus <- Paired.Microbiome.Transcriptome.GenusData.Transp
colnames(Y)
airway<-paired.microbe.SAMdata$Location_sort_code

Y<- Transcriptome.Data.Order
Y <- as.matrix(Y)
class(Y) <- "numeric"
Y[is.na(Y)]<-0
Y[Y == "n.a"] <- NA


## Variance decomposition 
Y.split <- compPLS:::.Split.variation.one.level(Y,
                                                #paired.microbe.SAMdata$Lung_Cancer,
                                                paired.microbe.SAMdata$`newlev$control3groups`,
                                                paired.microbe.SAMdata$Subject_ID) 


Y <- Y.split$Xw


otu.pclr <-  clr(t(Paired.Microbiome.Transcriptome.GenusData.Transp), 2)
otu.pclr <- t(otu.pclr)

otu.pclr.split <- compPLS:::.Split.variation.one.level(otu.pclr,
                                               #paired.microbe.SAMdata$Lung_Cancer,
                                               paired.microbe.SAMdata$`newlev$control3groups`,
                                               paired.microbe.SAMdata$Subject_ID) 

Y <- scale(Y, center=TRUE, scale=TRUE)
X <- scale(otu.pclr.split$Xw, center=TRUE, scale=FALSE) 

## Filter features by number of zeros
#Y<-Y[,apply(Y, 2, function(x)  sum(sign(x)==1,na.rm=TRUE)>30), drop=FALSE]
#X<-X[,apply(X, 2, function(x)  sum(sign(x)==1,na.rm=TRUE)>5), drop=FALSE]
#X <- X[-colVars(X)>0.2,]
### This works 
# X<-X[,apply(X, 2, function(x)  sum(sign(x)==1,na.rm=TRUE)>5), drop=FALSE]

## Check cross-covariance matrix for value of K 
XY <- cov(X, Y)
dim(XY)
tmp <- svd(XY)
tmp$d
barplot(tmp$d) 

K <-22

## Subsampling the data using sPLS
out.stars <- compPLS:::spls.stars(X, Y, rep.num=50, K=K, eta=(seq((.799), (.999), length.out=10)), ncores=1 )
out.stars$opt.ind <- which.min((out.stars$variability[,1][out.stars$variability[,1] >= .1]))
out.stars$eta.opt <- as.numeric(row.names(out.stars$variability)[out.stars$opt.ind])
out.spls <- spls::spls(X, Y, eta=out.stars$eta.opt, K=K)

## Prediction error
plot.new()
par(xpd=FALSE)
Ypred <- scale(X[,out.spls$A] %*% out.spls$betahat[out.spls$A,,drop=FALSE])
correlation.lm = lm(Ypred[,1]~Y[,1])
summary(correlation.lm)$r.squared
plot(Ypred[,1], Y[,1], xlab="Predicted Expression Level", ylab="Actual Expression Level")
abline(correlation.lm)


## plot Biplot
scores <- X[,out.spls$A,drop=FALSE] %*% out.spls$projection
taxtab <- as.matrix(paired.microbe.TAXTAB)[colnames(X[,out.spls$A]),]
taxtab<- apply(taxtab, 2, function(x) gsub("\\w__", "", x)) ## Remove labels from unknown levels (we use OTU Numbers to track)
class <- paste(taxtab[,5],"_",taxtab[,6]) ## set levels for labels
dim(taxtab)


tol21rainbow= c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")


basecols <-tol21rainbow
ccols <-  basecols[as.numeric(as.factor(class))]

legmat <- unique(cbind(class[order(as.numeric(as.factor(class)))],
                       ccols[order(as.numeric(as.factor(class)))]))


legmat <- unique(cbind(class,
                       ccols))
plot.new()
legend('left', legmat[,1], col=legmat[,2], lty=1, lwd=5,cex=0.3)

gg <- ggbiplot.splsda(out.spls, 
                      shape=as.factor(paired.microbe.SAMdata$Location_sort_code),size=3,
                      col=as.factor(paired.microbe.SAMdata$Lung_Cancer), sel=c(1,2),
                      col.loadings=ccols,scale.loadings=25,xlab="PLS 1", ylab="PLS 2") + theme_light() 
plot(gg)

## Calculate p-values from model co-efficients
.bfun <- function(x, indices, ...) spls::spls(x[indices,], Y[indices,], ...)$betahat
.pfun <- function(x, indices, ...) {
  x.rand <- apply(x, 2, function(x) sample(x))[indices, ]
  y.rand <- apply(Y, 2, function(x) sample(x))[indices, ]
  spls::spls(x.rand, y.rand, ...)$betahat
}

## Set R to number of bootstraps
out.spls.boot <- compPLS:::splsdaboot(X[,out.spls$A], .bfun, .pfun, eta=0, K=out.spls$K, R=5000)
pv <- compPLS:::pval.splsdaboot(out.spls.boot)


## Plot network and filter by FDR apdj
mat<-pv
mat <- melt(mat)
mat$value <- p.adjust(mat$value, method = "fdr")
mat<-mat[mat$value <= 0.1,] 

## Change node names
ig.raw <- graph.data.frame(mat, directed = FALSE)
vnames<-cbind(rownames(taxtab),taxtab[,6])
data.frame(vnames)
vnames<- gsub("g__", "",vnames)

ig.mat <- graph.data.frame(mat, directed = FALSE)

## Set node size to median relative abundance and scale
medians<-apply(as.matrix(Paired.Microbiome.Transcriptome.GenusData.Transp[,only]), 2, mean) *1000
library(scales)
V(ig.mat)$size[1:len]<-as.numeric(rescale(medians,to=c(10,17)))
V(ig.mat)$shape[1:len]<- "circle"
V(ig.mat)$color[1:len] <- "lightblue2"

## Scale edge width to pvalue
E(ig.mat)$width <-(1/E(ig.mat)$value)

## Network edge list for formatting edges
net <- ig.raw
el1 <- apply(get.edgelist(net), 1, paste, collapse="-")
el <- get.edgelist(net)

### Set edge color to coefficient red for negative, blue for positive
posorneg <- list(rep(TRUE, length(el[,1])))
posorneg <-unlist(posorneg)
for (i in 1:length(el[,1])){
  first <- el[,1][i]
  second <- el[,2][i] 
  if(out.spls$betahat[first,second] < 0) posorneg[i] <- FALSE
}
#E(ig.mat)$color <-  ifelse(posorneg, "#b3b3b3", "#fb6a4a")
E(ig.mat)$color <-  ifelse(posorneg, adjustcolor("blue", .6), adjustcolor("#fb6a4a", .6))

## Choose layout and plot network
l = layout.fruchterman.reingold(ig.mat)
l = layout.graphopt(ig.mat)
plot(layout=l,ig.mat,vertex.shape= V(ig.mat)$shape,vertex.color= V(ig.mat)$color,vertex.label.font=V(ig.mat)$label.font, 
     vertex.label=V(ig.mat)$name, vertex.frame.color="white",vertex.label.color="black",rescale=TRUE, 
     vertex.label.family="Helvetica",ylim=c(-1,1),xlim=c(-1,1), asp = 0)


