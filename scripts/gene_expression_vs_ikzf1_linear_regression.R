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


