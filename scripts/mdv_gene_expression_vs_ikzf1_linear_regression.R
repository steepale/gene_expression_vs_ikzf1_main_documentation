# Determine if there is a correlation via linear regression between quantities of functioning IKZF1 transcripts vs other genes in the mdv genome

# Set the working directory
setwd("/Users/Alec/Documents/Bioinformatics/MDV_Project/gene_expression_vs_ikzf1")
plot_dir <- "./analysis/plots/lm_mrna_counts_mdv/"
dir.create(plot_dir, showWarnings = FALSE)

# Create a table of samples and their appropriate groups, and feed this into an object 
targets <- read.delim("./data/targets_IKZF1.txt", stringsAsFactors=FALSE)
targets

# Group your samples accordingly, and feed into another object
group <- paste(targets$group, targets$description, targets$ikzf1_status, sep=".")
group <- factor(group)
table(group)

# Create an object with sample names
## Note: 017834-2 and 017756-3 removed because of low quality
samples <- c("017738-1", "017741-1", "017766-1", "017777-3", "017787-2", "017794-1", "017798-1_1", "017798-1_2", "017833-1", "017835-1", "017841-3", "017842-2_1", "017842-2_2", "017855-1_1", "017855-1_2", "017863-1", "017884-2", "017901-2_1", "017901-2_2", "017906-1", "017911-1_1", "017911-1_2", "017918-3", "017927-2")

# This is a function to read our gene count files that resulted from HTSeq
read.sample <- function(sample.name) {
        file.name <- paste("./data/gene_counts_mdv/", sample.name, "_ensembl_gene_counts.txt", sep="")
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

# Load all the sample variables into one data frame
all.data <- data.frame(sample.1, sample.2$count, sample.3$count, sample.4$count, sample.5$count, sample.6$count, sample.7$count, sample.8$count, sample.9$count, sample.10$count, sample.11$count, sample.12$count, sample.13$count, sample.14$count, sample.15$count, sample.16$count, sample.17$count, sample.18$count, sample.19$count, sample.20$count, sample.21$count, sample.22$count, sample.23$count, sample.24$count)
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
# Remove birds with mutated IKZF1 from df and low quality samples
#group <- group[!grepl("NC|Mutated", group)]
cohort_group = c(10:11, 13:17, 24)
group <- factor(group[cohort_group])
group <- factor(group)
group
# Remove those values from all.data
#cohort = c(1:3, 5:11, 14:17, 20, 24)
cohort_alldata = c(2:3, 5:9, 16)
all.data <- all.data[, cohort_alldata]

dge = DGEList(counts=all.data, group=group, genes=row.names(all.data))

# Add gene annotation 
# Load a gallus gallus annotation bank to add additional annotation for downstream analyses
library(org.Gg.eg.db)
keytypes(org.Gg.eg.db)
columns(org.Gg.eg.db)
head(keys(org.Gg.eg.db, "GENENAME"))

# Add "ENTREZID" and "GENENAME" annotation to dge$genes
dim(dge)

#dge$genes$EntrezID <- mapIds(org.Gg.eg.db, rownames(all.data), "ENTREZID", "ENSEMBL")
#dge$genes$GeneName <- mapIds(org.Gg.eg.db, rownames(all.data), "GENENAME", "ENSEMBL")
#dge$genes$Symbol <- mapIds(org.Gg.eg.db, rownames(all.data), "SYMBOL", "ENSEMBL")
#dge$counts[dge$genes$Symbol=="RFK", ]
#dge$counts[dge$genes$genes=="ENSGALG00000013086", ]
#write.table(all.data, file = "/Users/Alec/Documents/Bioinformatics/MDV_Project/rna_gene_expression_analysis/data/raw_counts_cohort.txt")

# Normalize the composition bias
# Normalization by trimmed mean of M values (TMM) is performed by calcNormFactors, 
# which calculates a set of normalization factors for each library to elimintae composition bias.
dge <- calcNormFactors(dge, method="TMM")
dge$samples

# Generate a list of samples and their IKAROS transcript amounts
e_df <- read.delim("/Users/Alec/Documents/Bioinformatics/MDV_Project/gene_expression_vs_ikzf1/data/IKZF1_counts_vs_mdv.txt", col.names=c("SAMPLE", "count"), sep=" ", colClasses=c("character", "numeric"))

row.names(dge$counts)

# Create a mock data.frame
pr_df <- data.frame("FILLER", 0, 0)
names(pr_df) <- c("GENE", "R_SQUARED", "P-VALUE")

# Iterate through each gene and compare it to IKZF1
for (gene in row.names(dge$counts)) {
        # Extract gene counts of target (response variable)
        r_df <- dge$counts[row.names(dge$counts)==gene, ]
        # Create a dataframe with these two gene's counts
        lin_reg = data.frame(
                response = r_df,
                IKZF1 = e_df[,2]
        )
        # Generate sample labels for dots
        names <- rownames(lin_reg)
        # Generate a scatterplot for linear regression analysis (make sure to label the samples)
        pdf(paste0(plot_dir, "IKZF1_vs_", gene, "_lin_reg_transcripts.pdf"))
        plot <- xyplot(response ~ IKZF1, data = lin_reg,
                       xlab = "mRNA Transcript Count for IKZF1",
                       ylab = paste0("mRNA Transcript Count for ", gene),
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
        pr_df_int <- data.frame(gene, r_sqr, p_val)
        names(pr_df_int) <- c("GENE", "R_SQUARED", "P-VALUE")
        # Combine dataframes
        pr_df <- rbind(pr_df_int, pr_df)
}

# Re-order the dataframe by R squared
sapply(pr_df, mode)
pr_df_sort <- pr_df[order(pr_df$R_SQUARED, decreasing=TRUE), ]
View(pr_df_sort)
