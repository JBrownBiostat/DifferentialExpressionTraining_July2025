##############################################################################
## Tutorial script for single-cell RNA-seq differential expression analysis ##
##############################################################################


###############
## Libraries ##
###############
library(DESeq2)
library(ggplot2)
library(TENxPBMCData)
library(Matrix)
library(irlba)
library(matrixStats)
library(sparseMatrixStats)
library(scran)
library(BiocParallel)
library(scater)


##########
## Data ##
##########
set.seed(93051421)
pbmc68k <- TENxPBMCData(dataset = "pbmc68k")
pbmc68k <- pbmc68k[, sample.int(ncol(pbmc68k), 5e3)]
counts(pbmc68k) <- as(counts(pbmc68k), "sparseMatrix")

dim(pbmc68k)
pbmc68k <- pbmc68k[rowSums(counts(pbmc68k) > 0) >= 1e-2 * ncol(pbmc68k), ]
dim(pbmc68k)

print(pbmc68k)
counts(pbmc68k)[1:8, 1:8]

prllWorkers <- 6 # number of cores (strictly threads) to use for parallel computation


################################
## Normalization size factors ##
################################
summary(rowMeans(counts(pbmc68k) > 0))
?scuttle::computePooledFactors
pbmc68k <- computePooledFactors(
    pbmc68k, BPPARAM = SnowParam(workers = prllWorkers)
)
?scuttle::logNormCounts
pbmc68k <- scuttle::logNormCounts(pbmc68k)


###################
## Quick cluster ##
###################
# Top PCs #
top_pbmc68k <- getTopHVGs(pbmc68k, n = 2500)
pbmc68k <- fixedPCA(pbmc68k, subset.row = top_pbmc68k, rank = 20) 
plotReducedDim(pbmc68k, dimred="PCA")

# umap dimension reduction #
pbmc68k <- runUMAP(pbmc68k, dimred="PCA")
plotReducedDim(pbmc68k, dimred="UMAP")

# clustering #
graphClusters <- clusterCells(pbmc68k, use.dimred="PCA")
graphClusters <- LETTERS[graphClusters]
table(graphClusters)
colLabels(pbmc68k) <- graphClusters
plotReducedDim(pbmc68k, "UMAP", colour_by="label")


#############################
## Differential expression ##
#############################
# Build DESeq2 model #
metadat <- data.frame(
    cluster = factor(graphClusters, levels = sort(unique(graphClusters)))
)
dds <- DESeqDataSetFromMatrix(
    countData = counts(pbmc68k),
    colData = metadat,
    design = ~ cluster - 1
)
sizeFactors(dds) <- pbmc68k$sizeFactor
?DESeq
dds <- DESeq(
    dds,
    minReplicatesForReplace = Inf,
    useT = TRUE,
    minmu = 1e-6,
    parallel = TRUE,
    BPPARAM = SnowParam(workers = prllWorkers)
)

# Calculate A-B test results #
resultsNames(dds)
shrunk_AvsB <- lfcShrink(
    dds,
    contrast = c("cluster", "A", "B"),
    type = "ashr",
    svalue = TRUE
)
head(shrunk_AvsB)
table(shrunk_AvsB$svalue <= 1e-3)

# Visualize model fit and results #
plotDispEsts(dds)
plot(
    x = log2(shrunk_AvsB$baseMean),
    y = shrunk_AvsB$log2FoldChange,
    col = ifelse(shrunk_AvsB$svalue <= 1e-3, "blue", "black"),
    pch = 19, cex = 0.25
)

# Split out cluster - A #
dds_A <- DESeqDataSetFromMatrix(
    countData = counts(pbmc68k)[, graphClusters == "A"],
    colData = metadat[graphClusters == "A", ],
    design = ~ 1
)
sizeFactors(dds_A) <- pbmc68k$sizeFactor[graphClusters == "A"]
dds_A <- DESeq(
    dds_A,
    minReplicatesForReplace = Inf,
    useT = TRUE,
    minmu = 1e-6,
    parallel = TRUE,
    BPPARAM = SnowParam(workers = prllWorkers)
)

# Split out cluster - B #
dds_B <- DESeqDataSetFromMatrix(
    countData = counts(pbmc68k)[, graphClusters == "B"],
    colData = metadat[graphClusters == "B", ],
    design = ~ 1
)
sizeFactors(dds_B) <- pbmc68k$sizeFactor[graphClusters == "B"]
dds_B <- DESeq(
    dds_B,
    minReplicatesForReplace = Inf,
    useT = TRUE,
    minmu = 1e-6,
    parallel = TRUE,
    BPPARAM = SnowParam(workers = prllWorkers)
)

# Visualize dispersion estimates #
par(mfrow = c(1, 3))
plotDispEsts(dds, main = "All cluster model")
plotDispEsts(dds_A, main = "Only cluster 'A'")
points(
    x = mcols(dds)$baseMean,
    y = mcols(dds)$dispFit,
    pch = 19, cex = 0.25, col = 3
)
plotDispEsts(dds_B, main = "Only cluster 'B'")
points(
    x = mcols(dds)$baseMean,
    y = mcols(dds)$dispFit,
    pch = 19, cex = 0.25, col = 3
)
par(mfrow = c(1, 1))

# Recalculating test results on gene-subset #
subCounts_A <- counts(pbmc68k)[, graphClusters == "A"]
subCounts_B <- counts(pbmc68k)[, graphClusters == "B"]
subGenes <- (rowSums(subCounts_A > 0) >= max(4, 5e-2 * ncol(subCounts_A))) | 
    (rowSums(subCounts_B > 0) >= max(4, 5e-2 * ncol(subCounts_B))) |
    (
        rowSums(subCounts_B > 0) >= max(2, 2.5e-2 * ncol(subCounts_B)) & 
        rowSums(subCounts_A > 0) >= max(2, 2.5e-2 * ncol(subCounts_A))
    )
table(subGenes)

shrunk_AvsB_subset <- lfcShrink(
    dds[subGenes, ],
    contrast = c("cluster", "A", "B"),
    type = "ashr",
    svalue = TRUE
)
head(shrunk_AvsB_subset)
table(shrunk_AvsB_subset$svalue <= 1e-3)

par(mfrow = c(1, 2))
plot(
    x = log2(shrunk_AvsB$baseMean),
    y = shrunk_AvsB$log2FoldChange,
    col = ifelse(shrunk_AvsB$svalue <= 1e-3, "blue", "black"),
    pch = 19, cex = 0.25,
    main = "Original all-cluster results"
)
plot(
    x = log2(shrunk_AvsB_subset$baseMean),
    y = shrunk_AvsB_subset$log2FoldChange,
    col = ifelse(shrunk_AvsB_subset$svalue <= 1e-3, "blue", "black"),
    pch = 19, cex = 0.25,
    main = "Subset all-cluster results"
)
par(mfrow = c(1, 1))

# Rebuild full model on gene/cluster subset #
dds_Targeted <- DESeqDataSetFromMatrix(
    countData = counts(pbmc68k)[subGenes, graphClusters %in% c("A", "B")],
    colData = metadat[graphClusters %in% c("A", "B"), , drop = FALSE],
    design = ~ cluster - 1
)
sizeFactors(dds_Targeted) <- pbmc68k$sizeFactor[graphClusters %in% c("A", "B")]
dds_Targeted <- DESeq(
    dds_Targeted,
    minReplicatesForReplace = Inf,
    useT = TRUE,
    minmu = 1e-6,
    parallel = TRUE,
    BPPARAM = SnowParam(workers = prllWorkers)
)
shrunk_AvsB_Targeted <- lfcShrink(
    dds_Targeted,
    contrast = c("cluster", "A", "B"),
    type = "ashr",
    svalue = TRUE
)
head(shrunk_AvsB_Targeted)
table(shrunk_AvsB_Targeted$svalue <= 1e-3)

# Compare dispersion fit #
par(mfrow = c(1, 2))
plotDispEsts(dds, main = "Original all-cluster dispersion")
plotDispEsts(dds_Targeted, main = "Gene/Cluster subset dispersion")
par(mfrow = c(1, 1))

# Compare testing results #
par(mfrow = c(1, 3))
plot(
    x = log2(shrunk_AvsB$baseMean),
    y = shrunk_AvsB$log2FoldChange,
    col = ifelse(shrunk_AvsB$svalue <= 1e-3, "blue", "black"),
    pch = 19, cex = 0.25,
    main = "Original all-cluster results"
)
plot(
    x = log2(shrunk_AvsB_subset$baseMean),
    y = shrunk_AvsB_subset$log2FoldChange,
    col = ifelse(shrunk_AvsB_subset$svalue <= 1e-3, "blue", "black"),
    pch = 19, cex = 0.25,
    main = "Subset all-cluster results"
)
plot(
    x = log2(shrunk_AvsB_Targeted$baseMean),
    y = shrunk_AvsB_Targeted$log2FoldChange,
    col = ifelse(shrunk_AvsB_Targeted$svalue <= 1e-3, "blue", "black"),
    pch = 19, cex = 0.25,
    main = "Gene/Cluster subset results"
)
par(mfrow = c(1, 1))






