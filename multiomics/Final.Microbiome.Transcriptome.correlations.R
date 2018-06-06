####### Codes for Microbiome analysis.  #####
#Load Phyloseq
library(phyloseq)

#Load other libraries: 
library("ade4")
library("vegan")




##########Correlation analysis between microbiome and trasncriptome
#Select only samples where transcriptomics are availabe
Paired.Microbiome.TranscriptomeGenus.Rel.table = subset_samples(Genus.rel.table, Rn.aseq=="1")
Paired.Microbiome.TranscriptomeGenus.Rel.table
otu_table()   OTU Table:         [ 1244 taxa and 87 samples ]
sample_data() Sample Data:       [ 87 samples by 102 sample variables ]
tax_table()   Taxonomy Table:    [ 1244 taxa by 7 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 1244 tips and 1243 internal nodes ]


# Prune the data at genus level: select for genera present in >2% relative abundance in 0.5% of the samples (this approach brings in > 70% of the data in almost all samples)
Paired.Microbiome.TranscriptomeGenus.Rel.wh1 = genefilter_sample(Paired.Microbiome.TranscriptomeGenus.Rel.table, filterfun_sample(function(x) x > 0.02), A = 0.005 * nsamples(Paired.Microbiome.TranscriptomeGenus.Rel.table))
Paired.Microbiome.TranscriptomeGenus.Rel.table1B = prune_taxa(Paired.Microbiome.TranscriptomeGenus.Rel.wh1, Paired.Microbiome.TranscriptomeGenus.Rel.table)
colnames(sample_data(Paired.Microbiome.TranscriptomeGenus.Rel.table1B))
plot_bar(Paired.Microbiome.TranscriptomeGenus.Rel.table1B, fill="Genus")



#SPEARMAN RHO AND p-VALUES TO NETWORK
library('Hmisc')
library("reshape2")

#Need to make a table out of the phyloseq object pruned genus table.
Paired.Microbiome.Transcriptome.GenusData = otu_table(Paired.Microbiome.TranscriptomeGenus.Rel.table1B)



##Now: Network for taxa vs. taxa
#spearman between taxa. In order to keep order, family, genus names through network, change the names of the OTU #s before calculating spearman 
#change names
#Here we are able to change the names for genuses that are labelled as "g__" 
Paired.New.Names = prune_taxa(tail(names(sort(taxa_sums(Paired.Microbiome.TranscriptomeGenus.Rel.table1B))), ntaxa(Paired.Microbiome.TranscriptomeGenus.Rel.table1B)), Paired.Microbiome.TranscriptomeGenus.Rel.table1B)
	tax_table(Paired.Microbiome.TranscriptomeGenus.Rel.table1B)
	# Add a new rank, Strain, with the Genus ids
	tax_table(Paired.New.Names) <- cbind(tax_table(Paired.New.Names), Strain=taxa_names(Paired.New.Names))
	# Define the ranks you want to include
	myranks = c("Order", "Family", "Genus", "Strain")
	mylabels = apply(tax_table(Paired.New.Names)[, myranks], 1, paste, sep="", collapse="|")
	# Add concatenated labels as a new rank after strain
	tax_table(Paired.New.Names) <- cbind(tax_table(Paired.New.Names), catglab=mylabels)
	#set new row names 
	rownames(Paired.Microbiome.Transcriptome.GenusData) <-as.matrix(tax_table(Paired.New.Names))[,"catglab"]

#Check names
rownames(Paired.Microbiome.Transcriptome.GenusData)


#spearman between samples 
Paired.Sample.Spear.dist = cor(Paired.Microbiome.Transcriptome.GenusData, method='spearman')

#View p values for samples 
Paired.Spear.Sammples.Sig = rcorr(Paired.Microbiome.Transcriptome.GenusData, type="spearman")$P


#spearman between taxa 
Paired.Taxa.Spear.dist = cor(t(Paired.Microbiome.Transcriptome.GenusData), method='spearman')

#View p values for taxa 
Paired.Taxa.Sammples.Sig = rcorr(t(Paired.Microbiome.Transcriptome.GenusData), type="spearman")$P


#Set correlations as a matrix for Taxa
Taxa.Spear.matrix = as.matrix(Paired.Taxa.Spear.dist)

#Set pval as a matrix 
Paired.Taxa.Sammples.Sig.matrix = as.matrix(Paired.Taxa.Sammples.Sig)

#set bottom half of all matrices to "dup" for elimination later 
Taxa.Spear.matrix.Final <- Taxa.Spear.matrix
Taxa.Spear.matrix.Final[upper.tri(Taxa.Spear.matrix.Final, diag=TRUE)] <- "dup"

Paired.Taxa.Sammples.Sig.matrix.Final <- Paired.Taxa.Sammples.Sig.matrix
Paired.Taxa.Sammples.Sig.matrix.Final[upper.tri(Paired.Taxa.Sammples.Sig.matrix.Final, diag=TRUE)] <- "dup"

#melt into columns and save as txt for excel 
write.table(melt(Taxa.Spear.matrix.Final), file = "SpearmanCorrelationMatrix.txt")
write.table(melt(Paired.Taxa.Sammples.Sig.matrix.Final), file = "SpearmanCorrelation.Significance.Matrix.txt")



#New code for Correlations with permutations (leaving one out)
#16S data: Paired.Microbiome.TranscriptomeGenus.Rel.table1B

#taxa by taxa
ntrys <- 100 #number of resampling

#Here we are able to change the names for genuses that are labelled as "g__" 
New.Names = prune_taxa(tail(names(sort(taxa_sums(Paired.Microbiome.TranscriptomeGenus.Rel.table1B))), ntaxa(Paired.Microbiome.TranscriptomeGenus.Rel.table1B)), Paired.Microbiome.TranscriptomeGenus.Rel.table1B)
tax_table(Paired.Microbiome.TranscriptomeGenus.Rel.table1B)
# Add a new rank, Strain, with the Genus ids
tax_table(New.Names) <- cbind(tax_table(New.Names), Strain=taxa_names(New.Names))
# Define the ranks you want to include
myranks = c("Order", "Family", "Genus", "Strain")
mylabels = apply(tax_table(New.Names)[, myranks], 1, paste, sep="", collapse="_")
# Add concatenated labels as a new rank after strain
tax_table(New.Names) <- cbind(tax_table(New.Names), catglab=mylabels)

#Make a table from the phyloseq object
GenusData.FullNames <- otu_table(Paired.Microbiome.TranscriptomeGenus.Rel.table1B)

#set new row names 
rownames(GenusData.FullNames) <-as.matrix(tax_table(New.Names))[,"catglab"]



#Start with Baseline taxa by taxa 
Genus.Microbiome.data <- as.data.frame(t(GenusData.FullNames))
resamples.b <- sapply(1: ntrys,  function(i) 
	Genus.Microbiome.data[sample(1:nrow(Genus.Microbiome.data), nrow(Genus.Microbiome.data)-1, replace=FALSE),])
# You can check resampling by calling resamples.b[5,5]
resamples.b[5,5]
resamples.b[5,99]
resamples.b[,99]

#Create a list out of the matrices of the RHO values 
Spearman.Resamples <- lapply(1: ntrys, function(i)
	rcorr(as.matrix(as.data.frame(resamples.b[,i])), type='spearman')$r)
#Create a list of the matrices of the p values 
Pval.Spearman.Resamples <- lapply(1: ntrys, function(i)
	rcorr(as.matrix(as.data.frame(resamples.b[,i])), type='spearman')$P)
	
#Median for RHO values 
#laply makes an array from a list 
#aaply splits an array, applies the function, and then returns a list
Median.Spear <- aaply(laply(Spearman.Resamples, as.matrix, na.rm=TRUE), c(2,3), median, na.rm=TRUE)
#Max for P values 
Max.Pval <- aaply(laply(Pval.Spearman.Resamples, as.matrix, na.rm=T), c(2,3), max, na.rm=T)

#make lists of values into a matrix 
#for Median RHO
Row.Col.Names <- rownames(Pval.Spearman.Resamples[[1]])
Median.Matrix.Spear <-matrix(unlist(Median.Spear), nrow=154, ncol=154)#number needs to change based on how many taxa there are in the table (Paired.Microbiome.TranscriptomeGenus.Rel.table1B)
rownames(Median.Matrix.Spear) <-Row.Col.Names
colnames(Median.Matrix.Spear) <- Row.Col.Names

#for Max P
Max.P.Matrix.Spear <-matrix(unlist(Max.Pval), nrow=154, ncol=154)
rownames(Max.P.Matrix.Spear) <-Row.Col.Names
colnames(Max.P.Matrix.Spear) <- Row.Col.Names

#set bottom half of all matrices to "dup" for elimination later 
Median.Spear.Final <- Median.Matrix.Spear
Median.Spear.Final[upper.tri(Median.Spear.Final, diag=TRUE)] <- "dup"

Max.P.Spear.Final <- Max.P.Matrix.Spear
Max.P.Spear.Final[upper.tri(Max.P.Spear.Final, diag=TRUE)] <- "dup"

#melt into columns and save as txt for excel 
write.table(melt(Median.Spear.Final), file = "Median_SpearmanCorrelationMatrix.100.txt")
write.table(melt(Max.P.Spear.Final), file = "MaxP_SpearmanPvalueMatrix.100.txt")


