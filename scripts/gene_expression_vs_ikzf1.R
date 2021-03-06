#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Set variables
comparison = args[1]
#comparison = "Tu.Female.Mutated-Tu.Female.Low_expression"

# Set the working directory
setwd("/Users/Alec/Documents/Bioinformatics/MDV_Project/gene_expression_vs_ikzf1")
plot_dir <- paste0("./analysis/plots/", comparison, "/")
dir.create(plot_dir, showWarnings = FALSE)

# Create a table of samples and their appropriate groups, and feed this into an object 
targets <- read.delim("./data/targets_IKZF1.txt", stringsAsFactors=FALSE)
targets

# Group your samples accordingly, and feed into another object
group <- paste(targets$group, targets$description, targets$ikzf1_status, sep=".")
group <- factor(group)
table(group)

# Create an object with sample names
samples <- c("017733", "017748", "017820", "017824", "017936", "017939", "017945", "017947", "017738-1", "017741-1", "017766-1", "017777-3", "017787-2", "017794-1", "017798-1_1", "017798-1_2", "017833-1", "017835-1", "017841-3", "017842-2_1", "017842-2_2", "017855-1_1", "017855-1_2", "017863-1", "017884-2", "017901-2_1", "017901-2_2", "017906-1", "017911-1_1", "017911-1_2", "017918-3", "017927-2")

# This is a function to read our gene count files that resulted from HTSeq
read.sample <- function(sample.name) {
        file.name <- paste("./data/counts/", sample.name, "_ensembl_gene_counts.txt", sep="")
        result <- read.delim(file.name, col.names=c("SYMBOL", "count"), sep="\t", colClasses=c("character", "numeric"), row.names=1)
}

# Read the samples into properly named variables
sample.1 <- read.sample(samples[1])
sample.2 <- read.sample(samples[2])
sample.3 <- read.sample(samples[3])
sample.4 <- read.sample(samples[4])
sample.5 <- read.sample(samples[5])
sample.6 <- read.sample(samples[6])
sample.7 <- read.sample(samples[7])
sample.8 <- read.sample(samples[8])
sample.9 <- read.sample(samples[9])
sample.10 <- read.sample(samples[10])
sample.11 <- read.sample(samples[11])
sample.12 <- read.sample(samples[12])
sample.13 <- read.sample(samples[13])
sample.14 <- read.sample(samples[14])
sample.15 <- read.sample(samples[15])
sample.16 <- read.sample(samples[16])
sample.17 <- read.sample(samples[17])
sample.18 <- read.sample(samples[18])
sample.19 <- read.sample(samples[19])
sample.20 <- read.sample(samples[20])
sample.21 <- read.sample(samples[21])
sample.22 <- read.sample(samples[22])
sample.23 <- read.sample(samples[23])
sample.24 <- read.sample(samples[24])
sample.25 <- read.sample(samples[25])
sample.26 <- read.sample(samples[26])
sample.27 <- read.sample(samples[27])
sample.28 <- read.sample(samples[28])
sample.29 <- read.sample(samples[29])
sample.30 <- read.sample(samples[30])
sample.31 <- read.sample(samples[31])
sample.32 <- read.sample(samples[32])

# Load all the sample variables into one data frame
all.data <- data.frame(sample.1, sample.2$count, sample.3$count, sample.4$count, sample.5$count, sample.6$count, sample.7$count, sample.8$count, sample.9$count, sample.10$count, sample.11$count, sample.12$count, sample.13$count, sample.14$count, sample.15$count, sample.16$count, sample.17$count, sample.18$count, sample.19$count, sample.20$count, sample.21$count, sample.22$count, sample.23$count, sample.24$count, sample.25$count, sample.26$count, sample.27$count, sample.28$count, sample.29$count, sample.30$count, sample.31$count, sample.32$count)
# Temporarily label the column names a certain way for the sake of an output file
# First, set row names to their own column
# library(data.table)
# setDT(all.data, keep.rownames = TRUE)[]
# colnames(all.data) <- c("ENSEMBL_GENE_ID", "017733", "017748", "017820", "017824", "017936", "017939", "017945", "017947", "017738-1", "017741-1", "017766-1", "017777-3", "017787-2", "017794-1", "017798-1_1", "017798-1_2", "017833-1", "017835-1", "017841-3", "017842-2_1", "017842-2_2", "017855-1_1", "017855-1_2", "017863-1", "017884-2", "017901-2_1", "017901-2_2", "017906-1", "017911-1_1", "017911-1_2", "017918-3", "017927-2")
# write.table(all.data, "/Users/Alec/Documents/Bioinformatics/MDV_Project/rna_gene_expression_analysis/data/gene_counts/expression_data_ensembl_chicken_gene_id.txt", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Label the columns of all.data
colnames(all.data)[1:ncol(all.data)] <- samples

# source("https://bioconductor.org/biocLite.R")
# biocLite("BiocUpgrade")
library("DESeq")
library("limma")
library("edgeR")

# Create a variable for different comparisons
if (comparison == "NC.Young.Wildtype-Tu.Female.Low_expression") {
        cohort <- c(1:4,9,19,22,23,25,28)
} else if (comparison == "NC.Mature.Wildtype-Tu.Female.Low_expression") {
        cohort <- c(5:8,9,19,22,23,25,28)
} else if (comparison == "NC.Young.Wildtype-Tu.Male.Low_expression") {
        cohort <- c(1:4,24)
} else if (comparison == "NC.Mature.Wildtype-Tu.Male.Low_expression") {
        cohort <- c(5:8,24)
} else if (comparison == "NC.Young.Wildtype-Tu.Female.Mutated") {
        cohort <- c(1:4,12,20,21,31)
} else if (comparison == "NC.Mature.Wildtype-Tu.Female.Mutated") {
        cohort <- c(5:8,12,20,21,31)
} else if (comparison == "NC.Young.Wildtype-Tu.Male.Mutated") {
        cohort <- c(1:4,26,27,29,30)
} else if (comparison == "NC.Mature.Wildtype-Tu.Male.Mutated") {
        cohort <- c(5:8,26,27,29,30)
} else if (comparison == "Tu.Female.Mutated-Tu.Female.Low_expression") {
        cohort <- c(12,20,21,31, 9,19,22,23,25,28)
} else if (comparison == "Tu.Male.Mutated-Tu.Male.Low_expression") {
        cohort <- c(26,27,29,30,24)
} else if (comparison == "Tu.Female.Wildtype-Tu.Female.Mutated") {
        cohort <- c(18,32,12,20,21,31)
} else if (comparison == "Tu.Male.Wildtype-Tu.Male.Mutated") {
        cohort <- c(10,11,13:17,26,27,29,30)
} else if (comparison == "Tu.Female.Wildtype-Tu.Female.Low_expression") {
        cohort <- c(18,32,9,19,22,23,25,28)
} else if (comparison == "Tu.Male.Wildtype-Tu.Male.Low_expression") {
        cohort <- c(10,11,13:17,24)
}

# Enter gene counts into a DGEList object with edgeR
group <- factor(group[cohort])
all.data <- all.data[,cohort]
dge = DGEList(counts=all.data, group=group, genes=row.names(all.data))

# Add gene annotation 
# Load a gallus gallus annotation bank to add additional annotation for downstream analyses
library(org.Gg.eg.db)
keytypes(org.Gg.eg.db)
columns(org.Gg.eg.db)
head(keys(org.Gg.eg.db, "GENENAME"))

# Add "ENTREZID" and "GENENAME" annotation to dge$genes
dim(dge)
#names(dge$genes)[names(dge$genes) == 'genes'] <- 'EnsemblID'
dge$genes$EntrezID <- mapIds(org.Gg.eg.db, rownames(all.data), "ENTREZID", "ENSEMBL")
dge$genes$GeneName <- mapIds(org.Gg.eg.db, rownames(all.data), "GENENAME", "ENSEMBL")
dge$genes$Symbol <- mapIds(org.Gg.eg.db, rownames(all.data), "SYMBOL", "ENSEMBL")
dge$counts[dge$genes$Symbol=="JUN",]
dge$counts[dge$genes$genes=="ENSGALG00000001826", ]
#all.data[rownames(all.data)=="ENSGALG00000009903", ]
#all.data[rownames(all.data)=="ENSGALG00000009903", "ENSGALG00000011849", "ENSGALG00000006827", "ENSGALG00000023822", "ENSGALG00000006801", "ENSGALG00000006329", "ENSGALG00000000892", "ENSGALG00000016677", "ENSGALG00000016678", "ENSGALG00000009904", "ENSGALG00000011844"]
#write.table(all.data, file = "/Users/Alec/Documents/Bioinformatics/MDV_Project/rna_gene_expression_analysis/data/raw_counts_cohort.txt")

# Remove any genes that did not have a corresponding ENTREZID
dge <- dge[!is.na(dge$genes$EntrezID), ]
dim(dge$genes)

# Remove genes that show low read counts (keep genes that have above 0.5 count per million in atleast 12 libraries; 12 males and 12 females)
# 0.5 comes from 10/L where L is the minimum library size in millions. Lowest library size: 20925538 (017901-2_1).
#keep <- rowSums(cpm(dge) > 0.5) >= 2
#table(keep)

# Subset the DGEList object to retain only the non-filtered genes
# kepp.lib.sizes=FALSE recalculates the library sizes post filtering 
#dge <- dge[keep, keep.lib.sizes=FALSE]

# Normalize the composition bias
# Normalization by trimmed mean of M values (TMM) is performed by calcNormFactors, 
# which calculates a set of normalization factors for each library to elimintae composition bias.
dge <- calcNormFactors(dge, method="TMM")
dge$samples

# Cluster the RNA samples in two dimensions using a multi-dimensional scaling (MDS) plot. 
# Looks at overall differences in expression profiles
pch <- c(0,1)
pdf(paste0(plot_dir, comparison, "_MDS.pdf"))
plot <- plotMDS(dge, pch=pch[group])  
# To determine the samples that cluster together. Note, this is not a very attractive graph; need to edit for future use
# plotMDS(dge, labels = rownames(dge$samples))
legend("topleft", legend=levels(group), pch=pch)
dev.off()

# Explore expression profiles of individual samples with mean-difference (MD) plots. 
# "The MD plot visualizes the library size-adjusted log-fold change between two libraries (the difference) against the
# average log-expression across those libraries (the mean)" ~Chen, Lun, Smyth 2016
#plotMD(dge, column=1)
#abline(h=0, col="red", lty=2, lwd=2)

# "The design matrix, which is required for DE analysis, records which treatment conditions were applied to 
# each samples, and it also defines how the experimental effects are parametrized in the linear models."
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
design

# "The dispersion parameter of the NB distribution accounts for variability between biological replicates. 
# edgeR estimates an empirical Bayes moderated dispersion for each individual gene. It also estimates a 
# common dispersion, which is a global dispersion estimate averaged over all genes, and a trended dispersion 
# where the dispersion of a gene is predicted from its abundance."
dge <- estimateDisp(dge, design, robust=TRUE)

# Plot the dispersion estimates with plotBCV
# "The vertical axis of the plotBCV plot shows square-root dispersion, also known as biological coefficient 
# of variation (BCV)"
# "the asymptotic value for the BCV tends to be in range from 0.05 to 0.2 for genetically identical mice 
# or cell lines, whereas somewhat larger values (> 0.3) are observed for human subjects."
#plotBCV(dge)

# "The NB model can be extended with quasi-likelihood (QL) methods to account for gene-specific variability 
# from both biological and technical sources7,12. Under the QL framework, the NB dispersion trend is used 
# to describe the overall biological variability across all genes, and gene-specific variability above and 
# below the overall level is picked up by the QL dispersion. In the QL approach, the individual (tagwise) 
# NB dispersions are not used."
# "Setting robust=TRUE in glmQLFit is usually recommended21. This allows gene-specific prior df estimates, 
# with lower values for outlier genes and higher values for the main body of genes. This reduces the 
# Chance of getting false positives from genes with extremely high or low raw dispersions, while at the 
# same time increasing statistical power to detect differential expression for the main body of genes."
fit <- glmQLFit(dge, design, robust=TRUE)
head(fit$coefficients)

# Visualize the QL dispersions with plotQLDisp
#plotQLDisp(fit)
summary(fit$df.prior)

# Test for differential gene expression between appropriate groups
# Group1 vs Group2
# A positive log2-fold-change (logFC) indicates up-regulation in Group1 relative to Group2 
# A negative log2-fold-change (logFC) indicates up-regulation in Group1 relative to Group2
#NMvsTM <- makeContrasts(NC.Mature-Tu.Male, levels=design)
vs <- makeContrasts(comparison, levels=design)

res <- glmQLFTest(fit, contrast=vs)
# View the top DE genes
topTags(res)

# The total number of DE genes identified at an FDR of 5% can be shown with decideTestsDGE
# Use decideTestsDGE to determine the number of DE genes with an FDR at or below 5%.
is.de <- decideTestsDGE(res)
summary(is.de)
topDE200 <- topTags(res, n=200)



# glmQLFTest shows differential expression no matter how small that difference is
# "A better and more rigorous approach is to modify the statistical test so as to 
# detect expression changes greater than a specified threshold."
# "Here we test whether the differential expression fold changes are significantly 
# greater than 1.5, that is, whether the logFCs are significantly greater than log2(1.5)"
# IMPORTANT: "The p-values from glmTreat are larger than those from glmQLFTest, and the 
# number of significantly DE genes is fewer, because it is testing an interval null 
# hypothesis and requires stronger evidence for differential expression than does a 
# conventional test. It provides greater specificity for identifying the most important 
# genes with large fold changes."
tr <- glmTreat(fit, contrast=vs, lfc=log2(1.5))
topDE200 <- topTags(tr, n=200)
topDE200
write.table(topDE200, file = paste0(plot_dir, comparison, "_topDE200.txt"), sep = "\t")

# Use decideTestsDGE to determine the number of DE genes with an FDR at or below 5%.
is.de <- decideTestsDGE(tr)
summary(is.de)

# Visualize the DE genes with an MD plot
pdf(paste0(plot_dir, comparison, "_MD.pdf"))
plot <- plotMD(res, status=is.de, values=c(1,-1), col=c("red","blue"), legend="topright")
dev.off()

# Generate a heat map
# Covert read counts to log2 counts per million (logCPM)
logCPM <- cpm(dge, prior.count=2, log=TRUE)
rownames(logCPM) <- dge$genes$Symbol
colnames(logCPM) <- paste(rownames(dge$samples))

# Choose the top 40 DE genes
o <- order(tr$table$PValue)
logCPM <- logCPM[o[1:40],]
logCPM

logCPM <- t(scale(t(logCPM)))

library(gplots)
col.pan <- colorpanel(100, "blue", "white", "red")
pdf(paste0(plot_dir, comparison, "_heatmap.pdf"))
plot <- heatmap.2(logCPM, col=col.pan, Rowv=TRUE, scale="none", trace="none", dendrogram="both", cexRow=1, cexCol=1.4, density.info="none", margin=c(10,9), lhei=c(2,10), lwid=c(2,6))
dev.off()
# "Because we have pre-standardized the rows of the logCPM matrix, 
# the Euclidean distance between each pair of genes is proportional to (1 − r)2, 
# where r is the Pearson correlation coefficient between the two genes. This 
# shows that the heatmap will cluster together genes that have positively 
# correlated logCPM values, because large positive correlations correspond to small distances."

# Load the appropriate chicken annotation databases
# biocLite("GO.db")
library(GO.db)
library(org.Gg.eg.db)

# Use goana to find the top differentiated or enriched gene ontology pathways
go <- goana(tr, geneid = tr$genes$EntrezID, species="Gg")
go50 <- topGO(go, n=50)
go50
write.table(go50, file = paste0(plot_dir, comparison, "_topgo50.txt"), sep = "\t")

# Use kegga to find the top differentiated or enriched KEGG pathways
keg <- kegga(tr, geneid = tr$genes$EntrezID, species="Gg")
kegg50 <- topKEGG(keg, n=50)
write.table(kegg50, file = paste0(plot_dir, comparison, "_topkegg50.txt"), sep = "\t")

# Install pathview to map RNA dif expressed genes and DNA variants to pathways
library(pathview)

# Obtain mapping data between compound or gene IDs and KEGG accessions
data(cpd.accs)
data(cpd.names)
data(kegg.met)
data(ko.ids)
data(rn.list)
data(gene.idtype.list)
data(gene.idtype.bods)
data(cpd.simtypes)

# Replace a dataframe that works with pathview (unknown why it doesn't work) with the column of data wanted
# Covert read counts to log2 counts per million (logCPM)
logCPM <- cpm(dge, prior.count=2, log=TRUE)

# Replace a dataframe that works with pathview (unknown why it doesn't work) with the column of data wanted
logCPM[, 1] <- tr$table[, "logFC"]
rownames(logCPM) <- dge$genes$genes
colnames(logCPM) <- paste(rownames(dge$samples))
#mka <- df[df$PATH=="04010", ]
#mkg <- mka[!is.na(mka$PATH), ]

# Need to change directory for purpose of saving info
setwd(plot_dir)
# Iterate through the top 50 kegg pathways and visualize the differences in gene expression
for (path in rownames(kegg50)) {
        fields <- strsplit(path, ":")[[1]]
        path_id <- (fields[2])
        print(path_id)
        path_name <- (gsub("\\s", "_", kegg50$Pathway[rownames(kegg50)==path]))
        print(path_name)
        pv.out <- pathview(gene.data = logCPM[, 1], gene.idtype = "ensembl", pathway.id = path_id, species = "gga", out.suffix = paste0(comparison, "_", path_name), kegg.native = T)
}

# Return to working directory
setwd("../../../")