#===============================================================================
#
#         FILE: /Users/Alec/Documents/Bioinformatics/MDV_Project/gene_expression_vs_ikzf1/gene_expression_vs_ikzf1_main_documentation.sh
#
#        USAGE: for documentation purposes, scripts inside
#
#  DESCRIPTION:  This script serves as a step by step documentation script and development script for determining if IKZF1 perterbation affects gene expression
#                
# REQUIREMENTS:  ---
#        NOTES:  ---
#       AUTHOR:  Alec Steep, steepale@msu.edu
#  AFFILIATION:  Michigan State University (MSU), East Lansing, MI, United States
#				         USDA ARS Avian Disease and Oncology Lab (ADOL), East Lansing, MI, United States
#				         Technical University of Munich (TUM), Weihenstephan, Germany
#      VERSION:  1.0
#      CREATED:  2017.06.22
#     REVISION:  
#===============================================================================

# PROJECT DIRECTORY (MacBook)
proj_dir="/Users/Alec/Documents/Bioinformatics/MDV_Project/gene_expression_vs_ikzf1"
cd $proj_dir

# Make proper directories
mkdir -p ./{data,scripts,analysis,jobs}
mkdir ./data/counts

# Sources to cite:
# http://f1000researchdata.s3.amazonaws.com/manuscripts/9996/16888710-6725-433b-ab2a-13c04bfe6cb5_8987_-_gordon_smyth_v2.pdf
# http://2014-msu-rnaseq.readthedocs.io/en/latest/model-orgs.html
# https://wikis.utexas.edu/display/bioiteam/Differential+gene+expression+analysis#Differentialgeneexpressionanalysis-Optional:edgeR

# Create a "targets_IKZF1.txt" file that categorizes samples for analysis 
vi ./data/targets_IKZF1.txt

#files	group	description	ikzf1_status
#017733_ensembl_gene_counts.txt	NC	Young	Wildtype
#017748_ensembl_gene_counts.txt	NC	Young	Wildtype
#017820_ensembl_gene_counts.txt	NC	Young	Wildtype
#017824_ensembl_gene_counts.txt	NC	Young	Wildtype
#017936_ensembl_gene_counts.txt	NC	Mature	Wildtype
#017939_ensembl_gene_counts.txt	NC	Mature	Wildtype
#017945_ensembl_gene_counts.txt	NC	Mature	Wildtype
#017947_ensembl_gene_counts.txt	NC	Mature	Wildtype
#017738-1_ensembl_gene_counts.txt	Tu	Female	Low_expression
#017741-1_ensembl_gene_counts.txt	Tu	Male	Wildtype
#017766-1_ensembl_gene_counts.txt	Tu	Male	Wildtype
#017777-3_ensembl_gene_counts.txt	Tu	Female	Mutated
#017787-2_ensembl_gene_counts.txt	Tu	Male	Wildtype
#017794-1_ensembl_gene_counts.txt	Tu	Male	Wildtype
#017798-1_1_ensembl_gene_counts.txt	Tu	Male	Wildtype
#017798-1_2_ensembl_gene_counts.txt	Tu	Male	Wildtype
#017833-1_ensembl_gene_counts.txt	Tu	Male	Wildtype
#017835-1_ensembl_gene_counts.txt	Tu	Female	Wildtype
#017841-3_ensembl_gene_counts.txt	Tu	Female	Low_expression
#017842-2_1_ensembl_gene_counts.txt	Tu	Female	Mutated
#017842-2_2_ensembl_gene_counts.txt	Tu	Female	Mutated
#017855-1_1_ensembl_gene_counts.txt	Tu	Female	Low_expression
#017855-1_2_ensembl_gene_counts.txt	Tu	Female	Low_expression
#017863-1_ensembl_gene_counts.txt	Tu	Male	Low_expression
#017884-2_ensembl_gene_counts.txt	Tu	Female	Low_expression
#017901-2_1_ensembl_gene_counts.txt	Tu	Male	Mutated
#017901-2_2_ensembl_gene_counts.txt	Tu	Male	Mutated
#017906-1_ensembl_gene_counts.txt	Tu	Female	Low_expression
#017911-1_1_ensembl_gene_counts.txt	Tu	Male	Mutated
#017911-1_2_ensembl_gene_counts.txt	Tu	Male	Mutated
#017918-3_ensembl_gene_counts.txt	Tu	Female	Mutated
#017927-2_ensembl_gene_counts.txt	Tu	Female	Wildtype
## Explanation of terms:
# NC: Normal Control
# Tu: Tumor
# Young: Uninfected bird, from which CD4+ RNA was extracted at 2 weeks of age
# Mature: Uninfected bird, from which CD4+ RNA was extracted at 6 weeks of age
# Female: Female Sex
# Male: Male Sex
# Wildtype: Seemingly normal IKZF1 function and expression
# Low_expression: IKZF1 with low gene expresison that is not mutated
# Mutated: Mutated IKZF1 with either a somatic snv or a somatic indel (usually in zinc finger binding domain)

# Transfer the gene counts data from the rna_gene_expression_analysis folder
rsync -avp \
/Users/Alec/Documents/Bioinformatics/MDV_Project/rna_gene_expression_analysis/data/gene_counts/*_ensembl_gene_counts.txt \
./data/counts/

# Generate comparisons to be made with these groups
# NC.Mature.Wildtype
# NC.Young.Wildtype
# Tu.Female.Low_expression
# Tu.Female.Mutated
# Tu.Female.Wildtype
# Tu.Male.Low_expression
# Tu.Male.Mutated
# Tu.Male.Wildtype

# Comparison 1: Compare normal samples to tumors with low IKZF1 expression
# Capture common DE genes to act as true positives
declare -a C1
C1[1]="NC.Mature.Wildtype-Tu.Female.Low_expression"
C1[2]="NC.Young.Wildtype-Tu.Female.Low_expression"
C1[3]="NC.Mature.Wildtype-Tu.Male.Low_expression"
C1[4]="NC.Young.Wildtype-Tu.Male.Low_expression"

# Comparison 2: Compare normal samples to tumors with mutated IKZF1
# Capture common DE genes
declare -a C2
C2[1]="NC.Mature.Wildtype-Tu.Female.Mutated"
C2[2]="NC.Young.Wildtype-Tu.Female.Mutated"
C2[3]="NC.Mature.Wildtype-Tu.Male.Mutated"
C2[4]="NC.Young.Wildtype-Tu.Male.Mutated"

# Comparison 3: Compare tumors without IKZF1 mutations with tumors with IKZF1 mutations
declare -a C3
C3[1]="Tu.Female.Wildtype-Tu.Female.Mutated"
C3[2]="Tu.Male.Wildtype-Tu.Male.Mutated"

# Comparison 4: Compare tumors with normal IKZF1 expression with tumors with low IKZF1 expression
declare -a C4
C4[1]="Tu.Female.Wildtype-Tu.Female.Low_expression"
C4[2]="Tu.Male.Wildtype-Tu.Male.Low_expression"

# Comparison 5: Compare IKZF1 mutated tumors to tumors with low IKZF1 expression
declare -a C5
C5[1]="Tu.Female.Mutated-Tu.Female.Low_expression"
C5[2]="Tu.Male.Mutated-Tu.Male.Low_expression"


### Run the majority of the analysis in R (interactive can be done in R studio)
for comp in ${C1[@]} ${C2[@]} ${C3[@]} ${C4[@]} ${C5[@]}
do
Rscript --vanilla ./scripts/gene_expression_vs_ikzf1.R ${comp}
done

# ./scripts/gene_expression_vs_ikzf1.R
#######################################
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
################################################

# Comparison 1: Compare normal samples to tumors with low IKZF1 expression
# Capture common DE genes to act as true positives
declare -a C1
C1[1]="NC.Mature.Wildtype-Tu.Female.Low_expression"
C1[2]="NC.Young.Wildtype-Tu.Female.Low_expression"
C1[3]="NC.Mature.Wildtype-Tu.Male.Low_expression"
C1[4]="NC.Young.Wildtype-Tu.Male.Low_expression"

# Comparison 2: Compare normal samples to tumors with mutated IKZF1
# Capture common DE genes
declare -a C2
C2[1]="NC.Mature.Wildtype-Tu.Female.Mutated"
C2[2]="NC.Young.Wildtype-Tu.Female.Mutated"
C2[3]="NC.Mature.Wildtype-Tu.Male.Mutated"
C2[4]="NC.Young.Wildtype-Tu.Male.Mutated"

# Comparison 3: Compare tumors without IKZF1 mutations with tumors with IKZF1 mutations
declare -a C3
C3[1]="Tu.Female.Wildtype-Tu.Female.Mutated"
C3[2]="Tu.Male.Wildtype-Tu.Male.Mutated"

# Comparison 4: Compare tumors with normal IKZF1 expression with tumors with low IKZF1 expression
declare -a C4
C4[1]="Tu.Female.Wildtype-Tu.Female.Low_expression"
C4[2]="Tu.Male.Wildtype-Tu.Male.Low_expression"

# Comparison 5: Compare IKZF1 mutated tumors to tumors with low IKZF1 expression
declare -a C5
C5[1]="Tu.Female.Mutated-Tu.Female.Low_expression"
C5[2]="Tu.Male.Mutated-Tu.Male.Low_expression"

# Determine the commonly differentially expressed genes from comparison 1
for comp in ${C1[@]}
do
# seperate gene symbols
cut -f5 ./analysis/plots/${comp}/${comp}_topDE200.txt | sort > ./analysis/plots/${comp}/${comp}_topDE200_gene_symbols.txt
done

# Compare the lists
# males
comm -1 -2 \
./analysis/plots/${C1[1]}/${C1[1]}_topDE200_gene_symbols.txt \
./analysis/plots/${C1[3]}/${C1[3]}_topDE200_gene_symbols.txt > ./data/${C1[1]}_vs_${C1[3]}_shared.txt
# females
comm -1 -2 \
./analysis/plots/${C1[2]}/${C1[2]}_topDE200_gene_symbols.txt \
./analysis/plots/${C1[4]}/${C1[4]}_topDE200_gene_symbols.txt > ./data/${C1[2]}_vs_${C1[4]}_shared.txt

# Grab the differentially mutated genes shared by all 4 comparisons
comm -1 -2 \
./data/${C1[1]}_vs_${C1[3]}_shared.txt \
./data/${C1[1]}_vs_${C1[3]}_shared.txt > \
./data/${C1}_shared.txt

# Determine the commonly differentially expressed genes from comparison 2
for comp in ${C2[@]}
do
# seperate gene symbols
cut -f5 ./analysis/plots/${comp}/${comp}_topDE200.txt | sort > ./analysis/plots/${comp}/${comp}_topDE200_gene_symbols.txt
done

# Determine commonly differentially expressed genes from comparisons 3
for comp in ${C3[@]}
do
# seperate gene symbols
cut -f5 ./analysis/plots/${comp}/${comp}_topDE200.txt | sort > ./analysis/plots/${comp}/${comp}_topDE200_gene_symbols.txt
done

# Compare the lists
comm -1 -2 \
./analysis/plots/${C3[1]}/${C3[1]}_topDE200_gene_symbols.txt \
./analysis/plots/${C3[2]}/${C3[2]}_topDE200_gene_symbols.txt > ./data/${C3[1]}_vs_${C3[2]}_shared.txt

# Determine commonly differentially expressed genes from comparisons 4
# Comparison 4: Compare tumors with normal IKZF1 expression with tumors with low IKZF1 expression
for comp in ${C4[@]}
do
# seperate gene symbols
cut -f5 ./analysis/plots/${comp}/${comp}_topDE200.txt | sort > ./analysis/plots/${comp}/${comp}_topDE200_gene_symbols.txt
done
# Compare the lists
comm -1 -2 \
./analysis/plots/${C4[1]}/${C4[1]}_topDE200_gene_symbols.txt \
./analysis/plots/${C4[2]}/${C4[2]}_topDE200_gene_symbols.txt > ./data/${C4[1]}_vs_${C4[2]}_shared.txt

# There are no shared DE genes between comparisons 3 and 4, suggesting that IKZF1 perturbation
# does not work in the same manner to regulate gene expression at the gene level.
# Perhaps there are not enough samples in each comparison and my comparisons are too stringent
# Instead compare the females from comparison 3 to the males in comparison 4 and vice versa
comm -1 -2 \
./analysis/plots/${C3[1]}/${C3[1]}_topDE200_gene_symbols.txt \
./analysis/plots/${C4[2]}/${C4[2]}_topDE200_gene_symbols.txt > \
./data/C3vsC4_possible_IKZF1_gene_targets.int

comm -1 -2 \
./analysis/plots/${C3[2]}/${C3[2]}_topDE200_gene_symbols.txt \
./analysis/plots/${C4[1]}/${C4[1]}_topDE200_gene_symbols.txt >> \
./data/C3vsC4_possible_IKZF1_gene_targets.int
# Sort the file and remove a header that got in there
sort ./data/C3vsC4_possible_IKZF1_gene_targets.int | \
grep -v "logFC" > \
./data/C3vsC4_possible_IKZF1_gene_targets.txt

rm ./data/C3vsC4_possible_IKZF1_gene_targets.int

# Determine commonly differentially expressed genes from comparisons 5
for comp in ${C5[@]}
do
# seperate gene symbols
cut -f5 ./analysis/plots/${comp}/${comp}_topDE200.txt | sort > ./analysis/plots/${comp}/${comp}_topDE200_gene_symbols.txt
done

# Compare the lists
comm -1 -2 \
./analysis/plots/${C5[1]}/${C5[1]}_topDE200_gene_symbols.txt \
./analysis/plots/${C5[2]}/${C5[2]}_topDE200_gene_symbols.txt > ./data/${C5[1]}_vs_${C5[2]}_shared.txt

# Males comparisons 1, 2, 3, 4
cat \
./analysis/plots/${C1[1]}/${C1[1]}_topDE200_gene_symbols.txt \
./analysis/plots/${C1[3]}/${C1[3]}_topDE200_gene_symbols.txt \
./analysis/plots/${C3[2]}/${C3[2]}_topDE200_gene_symbols.txt \
./analysis/plots/${C4[2]}/${C4[2]}_topDE200_gene_symbols.txt | sort | uniq -c | sort -n
# Female comparison 1, 2, 3, 4
cat \
./analysis/plots/${C1[2]}/${C1[2]}_topDE200_gene_symbols.txt \
./analysis/plots/${C1[4]}/${C1[4]}_topDE200_gene_symbols.txt \
./analysis/plots/${C3[1]}/${C3[1]}_topDE200_gene_symbols.txt \
./analysis/plots/${C4[1]}/${C4[1]}_topDE200_gene_symbols.txt | sort | uniq -c | sort -n

# Comparison 5: Compare IKZF1 mutated tumors to tumors with low IKZF1 expression
comm -1 -2 \
./analysis/plots/${C5[1]}/${C5[1]}_topDE200_gene_symbols.txt \
./analysis/plots/${C5[2]}/${C5[2]}_topDE200_gene_symbols.txt > ./data/${C5[1]}_vs_${C5[2]}_shared.txt

# It seems that the prior process of elimination may have been too stringent witht he limited amount of samples
# we possess. Instead, I will sum all comparisons and take the top genes that are differentially expressed
# in each across comparisons. Then, I'll perform linear regression analyses on these genes to see if there is a correlation
# between their mRNA counts and those of IKZF1.

# Create a list of comparison targets with ensembl gene id and gene symbol
for comp in ${C1[@]} ${C2[@]} ${C3[@]} ${C4[@]} ${C5[@]} 
do
# seperate gene symbols
cut -f2,5 ./analysis/plots/${comp}/${comp}_topDE200.txt | grep -v "logFC" | sort > ./analysis/plots/${comp}/${comp}_topDE200_gene_ids_symbols.txt
done

# Genes at the top are more common
printf "ENSEMBL_GENE_ID\tENSEMBL_GENE_SYMBOL\n" > ./data/topDE_genes_across_comparisons_IKZF1_targets.txt
cat \
./analysis/plots/${C1[1]}/${C1[1]}_topDE200_gene_ids_symbols.txt \
./analysis/plots/${C1[3]}/${C1[3]}_topDE200_gene_ids_symbols.txt \
./analysis/plots/${C2[1]}/${C2[1]}_topDE200_gene_ids_symbols.txt \
./analysis/plots/${C2[2]}/${C2[2]}_topDE200_gene_ids_symbols.txt \
./analysis/plots/${C2[3]}/${C2[3]}_topDE200_gene_ids_symbols.txt \
./analysis/plots/${C2[4]}/${C2[4]}_topDE200_gene_ids_symbols.txt \
./analysis/plots/${C3[2]}/${C3[2]}_topDE200_gene_ids_symbols.txt \
./analysis/plots/${C4[2]}/${C4[2]}_topDE200_gene_ids_symbols.txt \
./analysis/plots/${C1[2]}/${C1[2]}_topDE200_gene_ids_symbols.txt \
./analysis/plots/${C1[4]}/${C1[4]}_topDE200_gene_ids_symbols.txt \
./analysis/plots/${C3[1]}/${C3[1]}_topDE200_gene_ids_symbols.txt \
./analysis/plots/${C4[1]}/${C4[1]}_topDE200_gene_ids_symbols.txt | \
sort | uniq -c | sort -nr | cut -d" " -f 5 | sed 's/"//g' >> \
./data/topDE_genes_across_comparisons_IKZF1_targets.txt

# Perform a linear regression anlysis on all potentially differentially expressed genes to determine
# which genes are most associated with IKZF1 expression levels

# ./scripts/gene_expression_vs_ikzf1_linear_regression.R
##########################################
# Determine if there is a correlation via linear regression between quantities of functioning IKZF1 transcripts vs other genes in the genome

# Set the working directory
setwd("/Users/Alec/Documents/Bioinformatics/MDV_Project/gene_expression_vs_ikzf1")
plot_dir <- "./analysis/plots/lm_mrna_counts/"
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

# Enter gene counts into a DGEList object with edgeR
# Remove birds with mutated IKZF1 from df
cohort <- c(1:11, 13:19, 22:25, 32)
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
dge$counts[dge$genes$Symbol=="RFK", ]
dge$counts[dge$genes$genes=="ENSGALG00000013086", ]
#all.data[rownames(all.data)=="ENSGALG00000009903", ]
#all.data[rownames(all.data)=="ENSGALG00000009903", "ENSGALG00000011849", "ENSGALG00000006827", "ENSGALG00000023822", "ENSGALG00000006801", "ENSGALG00000006329", "ENSGALG00000000892", "ENSGALG00000016677", "ENSGALG00000016678", "ENSGALG00000009904", "ENSGALG00000011844"]
#write.table(all.data, file = "/Users/Alec/Documents/Bioinformatics/MDV_Project/rna_gene_expression_analysis/data/raw_counts_cohort.txt")

# Normalize the composition bias
# Normalization by trimmed mean of M values (TMM) is performed by calcNormFactors, 
# which calculates a set of normalization factors for each library to elimintae composition bias.
dge <- calcNormFactors(dge, method="TMM")
dge$samples

# Load the most commonly affected genes from IKZF1 alteration
IKZF1_targets <- read.table("./data/topDE_genes_across_comparisons_IKZF1_targets.txt", sep="\t", header=TRUE)

# Create a mock data.frame
pr_df <- data.frame(5000, "FILLER", 0, 0)
names(pr_df) <- c("RANK", "SYMBOL", "R_SQUARED", "P-VALUE")

# Loop through all the samples
for (id in IKZF1_targets[,1]) {
        symbol_int <- IKZF1_targets[IKZF1_targets$ENSEMBL_GENE_ID==id, ]
        symbol <- symbol_int$ENSEMBL_GENE_SYMBOL
        rank <- as.numeric(rownames(IKZF1_targets[IKZF1_targets$ENSEMBL_GENE_ID==id, ]))
        print(IKZF1_targets[IKZF1_targets$ENSEMBL_GENE_ID==id, ])
        # Extract gene counts of target (response variable)
        r_df <- dge$counts[dge$genes$genes==id, ]
        # Extract Gene counts of IKZF1 (explanatory variable)
        e_df <- as.vector(dge$counts[dge$genes$genes=="ENSGALG00000013086", ])
        # Create a dataframe with these two gene's counts
        lin_reg = data.frame(
                response = r_df,
                IKZF1 = e_df
        )
        # Generate sample labels for dots
        names <- rownames(lin_reg)
        # Generate a scatterplot for linear regression analysis (make sure to label the samples)
        pdf(paste0(plot_dir, "IKZF1_vs_", symbol, "_lin_reg_transcripts.pdf"))
        plot <- xyplot(response ~ IKZF1, data = lin_reg,
               xlab = "mRNA Transcript Count for IKZF1",
               ylab = paste0("mRNA Transcript Count for ", symbol),
               main = "Linear Regression Analysis IKZF1-vs-Target Expression",
               panel=function(x, y, ...) {
                       panel.xyplot(x, y, ...);
                       ltext(x=x, y=y, labels=names, pos=1, offset=1, cex=0.3)
               }
        )
        print(plot)
        dev.off()
        # Fit a linear model to the data
        lm_IKZF1 <- lm(response ~ IKZF1, data = lin_reg)
        
        # Capture p-value and r-squared
        p_val <- summary(lm_IKZF1)$coefficients[2, 4]
        r_sqr <- summary(lm_IKZF1)$r.squared
        # Create temp dataframe
        pr_df_int <- data.frame(rank, symbol, r_sqr, p_val)
        names(pr_df_int) <- c("RANK", "SYMBOL", "R_SQUARED", "P-VALUE")
        # Combine dataframes
        pr_df <- rbind(pr_df_int, pr_df)
}

# Re-order the dataframe by rank
sapply(pr_df, mode)
pr_df_sort <- pr_df[order(pr_df$R_SQUARED, decreasing=TRUE), ]
View(pr_df_sort)

# Sample by random
pr_df_random <- pr_df_sort[sample(nrow(pr_df_sort), 12), ]


# Test counts for certain genes
# WWC1
r_df <- as.vector(dge$counts[dge$genes$genes=="ENSGALG00000001826", ])
# CHD4
r_df <- dge$counts[dge$genes$genes=="ENSGALG00000013565", ]
# IKZF1
e_df <- as.vector(dge$counts[dge$genes$genes=="ENSGALG00000013086", ])

# Create a dataframe with these two gene's counts
lin_reg = data.frame(
        response = r_df,
        IKZF1 = e_df
)
# Remove outliers
#lin_reg <- lin_reg[1:21,]
#lin_reg <- log(lin_reg)

# Generate sample labels for dots
names <- rownames(lin_reg)

# Generate a scatterplot for linear regression analysis
xyplot(response ~ IKZF1, data = lin_reg,
       xlab = "mRNA Transcript Count for IKZF1",
       ylab = "mRNA Transcript Count for CHD4",
       main = "Linear Regression Analysis IKZF1-vs-Target Expression",
       panel=function(x, y, ...) {
               panel.xyplot(x, y, ...);
               ltext(x=x, y=y, labels=names, pos=1, offset=2, cex=0.3)
       }
)

# Fit a linear model to the data
lm_IKZF1 <- lm(response ~ IKZF1, data = lin_reg)

# Generate a summary of the linear model
summary(lm_IKZF1)
##########################################






