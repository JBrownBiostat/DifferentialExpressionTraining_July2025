#######################################################################
## Tutorial script for bulk RNA-seq differential expression analysis ##
#######################################################################


###############
## Libraries ##
###############
library(pasilla)
library(DESeq2)
library(ggplot2)
library(matrixStats)


####################################################
## Data: scripts adapted from the DESeq2 vignette ##
####################################################
# Import Pasilla Data #
cntFilePath <- system.file(
    "extdata",
    "pasilla_gene_counts.tsv",
    package="pasilla", 
    mustWork=TRUE
)
metaFilePath <- system.file(
    "extdata",
    "pasilla_sample_annotation.csv",
    package="pasilla", 
    mustWork=TRUE
)
cntMat <- as.matrix(read.csv(cntFilePath, sep = "\t", row.names = "gene_id"))
cntMat <- cntMat[rowSums(cntMat > 0) >= 3, ]
metaDat <- read.csv(metaFilePath, row.names = 1)
metaDat$condition <- factor(metaDat$condition)
metaDat$type <- factor(metaDat$type)

# Check sample structure #
print(metaDat)

# Import ZBTB3 Data #
load(paste(
    dirname(rstudioapi::getSourceEditorContext()$path),
    "/ZBTB3_TestData.RData",
    sep = ""
))
# Note: if the above fails, manually change
# 'dirname(rstudioapi::getSourceEditorContext()$path)' to the folder containing
# the data 'ZBTB3_TestData.RData', by default this is the same folder as 
# contains this script
print(metaDat_zbtb3)
# ZBTB3 data are a cleaned subset of the data from Paul Park, et al., Molecular Cell, 2024


#################################
## Quick-start simple analysis ##
#################################
# Match rows of meta data against columns of cntMat #
rownames(metaDat) <- gsub("fb", "", rownames(metaDat))
all(rownames(metaDat) %in% colnames(cntMat))
all(colnames(cntMat) %in% rownames(metaDat))
cntMat <- cntMat[, rownames(metaDat)]
identical(colnames(cntMat), rownames(metaDat))

# Setup DESeqDataSet #
dds_quickStart <- DESeqDataSetFromMatrix(
    countData = round(cntMat),
    colData = metaDat,
    design = ~ condition
)

# Run analysis and view results #
dds_quickStart <- DESeq(dds_quickStart)
head(results(dds_quickStart))
plotMA(dds_quickStart)


##################################
## Custom normalization factors ##
##################################
# Factors by library size: Reads per million #
rpmSF <- colSums(cntMat) / 1e6
normMat_rpm <- t(t(cntMat) / rpmSF)

# Factors by Median-of-Ratios #
filterCnts <- cntMat[rowSums(cntMat <= 5) == 0, ]
dim(filterCnts)
filterGeoMeans <- exp(rowMeans(log(filterCnts)))
mrSF <- as.numeric(apply(filterCnts / filterGeoMeans, 2, median))
normMat_mr <- t(t(cntMat) / mrSF)

# Factors calculated by default #
defaultSF <- sizeFactors(dds_quickStart)

# Compare custom to default size factors #
plotOrd <- order(defaultSF)
plot(
    x = log2(defaultSF[plotOrd]),
    y = log2(rpmSF[plotOrd]),
    type = "b", col = "red", pch = 19, 
    ylim = range(log2(c(rpmSF, mrSF))) + c(-0.1, 0.1),
    xlim = range(log2(defaultSF)) + c(-0.1, 0.1),
    xlab = "Default Factors (log2)", ylab = "Custom Factors (log2)"
)
lines(
    x = log2(defaultSF[plotOrd]),
    y = log2(mrSF[plotOrd]),
    type = "b", col = "blue", pch = 19, lty = 2
)
abline(0, 1, col = "black", lty = 2)
legend("topleft", legend = c("RPM", "MR"), col = c("red", "blue"), lty = 1:2)

# Visualize normalization quality #
normDat <- data.frame(
    normExpr = c(as.numeric(normMat_rpm), as.numeric(normMat_mr)),
    normMethod = rep(
        c("RPM", "MR"),
        rep(prod(dim(cntMat)), 2)
    ),
    sample = rep(
        rep(
            colnames(cntMat),
            rep(nrow(cntMat), ncol(cntMat))
        ),
        2
    )
)
normQualityPlot <- ggplot(
    normDat, 
    aes(x = sample, y = normExpr, color = sample, fill = sample)
) + 
    theme_classic() + 
    facet_grid(. ~ normMethod) + 
    geom_boxplot(alpha = 0.5) + 
    scale_color_brewer(palette = "Dark2") + 
    scale_fill_brewer(palette = "Dark2") + 
    scale_y_continuous(
        trans = "log1p",
        breaks = 10 ^ (0:6)
    ) + 
    labs(
        x = "Sample", y = "Normalized expression", color = "", fill = ""
    ) + 
    theme(
        axis.text.x = element_text(angle = -25, hjust = 0)
    )
normQualityPlot

# DESeq Analysis: RPM #
dds_RPM <- DESeqDataSetFromMatrix(
    countData = round(cntMat),
    colData = metaDat,
    design = ~ condition
)
sizeFactors(dds_RPM) <- rpmSF
dds_RPM <- DESeq(dds_RPM)
head(results(dds_RPM))

# DESeq Analysis: MR #
dds_MR <- DESeqDataSetFromMatrix(
    countData = round(cntMat),
    colData = metaDat,
    design = ~ condition
)
sizeFactors(dds_MR) <- mrSF
dds_MR <- DESeq(dds_MR)
head(results(dds_MR))

# Visuaize #
par(mfrow = c(1, 2))
plotMA(dds_RPM, main = "RPM - norm.")
plotMA(dds_MR, main = "MR - norm.")
par(mfrow = c(1, 1))


######################################################################
## Exploring estimated dispersion parameters - normalization effect ##
######################################################################
par(mfrow = c(1, 2))
# Dispersion plot: MR #
plotDispEsts(dds_MR, main = "MR - norm.")

# Dispersion plot: RPM #
plotDispEsts(dds_RPM, main = "RPM - norm.")
par(mfrow = c(1, 1))

# Visualize mean-variance relationship #
meanVarDat <- data.frame(
    meanVals = c(
        rowMeans(cntMat),
        rowMeans(normMat_mr),
        rowMeans(normMat_rpm)
    ),
    varVals = c(
        rowVars(cntMat),
        rowVars(normMat_mr),
        rowVars(normMat_rpm)
    ),
    normType = factor(
        rep(c("unnormalized", "mr", "rpm"), each = nrow(cntMat)),
        levels = c("unnormalized", "mr", "rpm"), ordered = TRUE
    )
)
fitDat <- data.frame(
    x = exp(seq(0, log(max(meanVarDat$meanVals)), length = 400))
)
fitDat$y1 <- fitDat$x + 0.01 * fitDat$x ^ 2
fitDat$y2 <- fitDat$x + 0.1 * fitDat$x ^ 2
meanVarPlot <- ggplot(
    meanVarDat, aes(x = meanVals, y = varVals, col = normType)
) + 
    theme_classic() + 
    facet_grid(. ~ normType) + 
    geom_point(size = 0.5) + 
    geom_abline(slope = 1, intercept = 0, lty = 1) + 
    geom_line(data = fitDat, aes(x = x, y = y1), lty = 2, inherit.aes = FALSE) +
    geom_line(data = fitDat, aes(x = x, y = y2), lty = 3, inherit.aes = FALSE) +
    scale_x_continuous(
        trans = "sqrt", 
        breaks = round(((0:10) * 50) ^ 2),
        limits = c(0, quantile(meanVarDat$meanVals, 0.995))
    ) + 
    scale_y_continuous(
        trans = "sqrt", 
        breaks = round(((0:10) * 1e3) ^ 2),
        limits = c(0, quantile(meanVarDat$varVals, 0.995))
    ) + 
    labs(
        x = "Gene means", y = "Gene vars.",
        title = "Mean-variance relationship",
        subtitle = "Solid line: mean = var (phi = 0) -- Dashed line: phi = 0.01 -- dotted line: phi = 0.1"
    ) +
    theme(
        legend.position = "none",
        axis.text.x = element_text(angle = -25)
    )
meanVarPlot

meanVarPlot_Log <- ggplot(
    meanVarDat, aes(x = meanVals, y = varVals, col = normType)
) + 
    theme_classic() + 
    facet_grid(. ~ normType) + 
    geom_point(size = 0.5) + 
    geom_abline(slope = 1, intercept = 0, lty = 1) + 
    geom_line(data = fitDat, aes(x = x, y = y1), lty = 2, inherit.aes = FALSE) +
    geom_line(data = fitDat, aes(x = x, y = y2), lty = 3, inherit.aes = FALSE) +
    scale_x_continuous(
        trans = "log1p", 
        breaks = c(1, 1e2, 1e4, 1e6, 1e8),
        limits = c(0, quantile(meanVarDat$meanVals, 0.995))
    ) + 
    scale_y_continuous(
        trans = "log1p", 
        breaks = 10 ^ (0:10),
        limits = c(0, quantile(meanVarDat$varVals, 0.995))
    ) + 
    labs(
        x = "Gene means", y = "Gene vars.",
        title = "Mean-variance relationship",
        subtitle = "Solid line: mean = var (phi = 0) -- Dashed line: phi = 0.01 -- dotted line: phi = 0.1"
    ) +
    theme(
        legend.position = "none",
        axis.text.x = element_text(angle = -25)
    )
meanVarPlot_Log

plot(
    x = log1p(rowMeans(normMat_rpm)),
    y = log1p(rowVars(normMat_rpm)) - log1p(rowMeans(normMat_rpm)),
    pch = 19, cex = 0.25,
    xlab = "gene means (log1p)", ylab = "gene var - gene mean (log1p)"
)
abline(h = 0, lty = 2, col = 2)
mean(log1p(rowVars(normMat_rpm)) - log1p(rowMeans(normMat_rpm)) < 0)

# Re-normalize relative to library size #
rpmScale <- rpmSF / exp(mean(log(rpmSF)))
rpmScale
scaleNorm <- t(t(cntMat) / rpmScale)

# DESeq Analysis: re-scaled RPM #
dds_scaleRPM <- DESeqDataSetFromMatrix(
    countData = round(cntMat),
    colData = metaDat,
    design = ~ condition
)
sizeFactors(dds_scaleRPM) <- rpmScale
dds_scaleRPM <- DESeq(dds_scaleRPM)
head(results(dds_scaleRPM))
plotMA(dds_scaleRPM)
plotDispEsts(dds_scaleRPM)


##############################################################
## Exploring estimated dispersion parameters - batch effect ##
##############################################################
# Build multimodal RNAseq/ChIPseq analysis #
print(metaDat_zbtb3)
dds_Merged <- DESeqDataSetFromMatrix(
    countData = countDat_zbtb3,
    colData = metaDat_zbtb3, 
    design = ~ Treatment - 1
)
sizeFactors(dds_Merged) <- metaDat_zbtb3$SF
dds_Merged <- DESeq(dds_Merged)
plotDispEsts(dds_Merged)
par(mfrow = c(1, 2))
plotMA(
    results(
        dds_Merged, contrast = list(c("TreatmentRNA_DOX_p402"), c("TreatmentRNA_DOX_p401"))
    ),
    main = "RNA-seq"
)
plotMA(
    results(
        dds_Merged, contrast = list(c("TreatmentChIP_DOX_p402"), c("TreatmentChIP_DOX_p401"))
    ),
    main = "ChIP-seq"
)
par(mfrow = c(1, 1))

# Separate RNA-seq only model #
dds_RNA <- DESeqDataSetFromMatrix(
    countData = countDat_zbtb3[, metaDat_zbtb3$Method == "RNAseq"],
    colData = metaDat_zbtb3[metaDat_zbtb3$Method == "RNAseq", ], 
    design = ~ Treatment - 1
)
sizeFactors(dds_RNA) <- metaDat_zbtb3$SF[metaDat_zbtb3$Method == "RNAseq"]
dds_RNA <- DESeq(dds_RNA)

# Separate ChIP-seq only model #
dds_ChIP <- DESeqDataSetFromMatrix(
    countData = countDat_zbtb3[, metaDat_zbtb3$Method == "ChIPseq"],
    colData = metaDat_zbtb3[metaDat_zbtb3$Method == "ChIPseq", ], 
    design = ~ Treatment - 1
)
sizeFactors(dds_ChIP) <- metaDat_zbtb3$SF[metaDat_zbtb3$Method == "ChIPseq"]
dds_ChIP <- DESeq(dds_ChIP)

# Visualize new dispersion plots #
par(mfrow = c(1, 3))
plotDispEsts(dds_Merged, main = "Merged model")
plotDispEsts(dds_RNA, main = "RNA only")
plotDispEsts(dds_ChIP, main = "ChIP only")
par(mfrow = c(1, 3))

# Visualize new MA plots #
par(mfrow = c(2, 2))
plotMA(
    results(
        dds_Merged, contrast = list(c("TreatmentRNA_DOX_p402"), c("TreatmentRNA_DOX_p401"))
    ),
    main = "Merged: RNA-seq"
)
plotMA(
    results(
        dds_Merged, contrast = list(c("TreatmentChIP_DOX_p402"), c("TreatmentChIP_DOX_p401"))
    ),
    main = "Merged: ChIP-seq"
)
plotMA(
    results(
        dds_RNA, contrast = list(c("TreatmentRNA_DOX_p402"), c("TreatmentRNA_DOX_p401"))
    ),
    main = "Seprate: RNA-seq"
)
plotMA(
    results(
        dds_ChIP, contrast = list(c("TreatmentChIP_DOX_p402"), c("TreatmentChIP_DOX_p401"))
    ),
    main = "Seprate: ChIP-seq"
)
par(mfrow = c(1, 1))


###########################################
## Design matrix and contrasts - Pasilla ##
###########################################
# Default model matrix #
model.matrix(~ condition, metaDat)

# Default coefficients #
resultsNames(dds_quickStart)
quickStartRes <- results(dds_quickStart)
head(quickStartRes)

# Specifying comparison through contrasts #
dds_Simple <- DESeqDataSetFromMatrix(
    countData = round(cntMat),
    colData = metaDat,
    design = ~ condition - 1
)
dds_Simple <- DESeq(dds_Simple)
resultsNames(dds_Simple)
?results
simpleRes1 <- results(
    dds_Simple, 
    contrast = c("condition", "treated", "untreated")
)
head(simpleRes1)
simpleRes2 <- results(
    dds_Simple, 
    contrast = list(
        c("conditiontreated"),
        c("conditionuntreated")
    )
)
head(simpleRes2)
simpleRes3 <- results(
    dds_Simple, 
    contrast = c(1,-1)
)
head(simpleRes3)

# Full model #
model.matrix(~ condition + type - 1, data = metaDat)
dds_Full <- DESeqDataSetFromMatrix(
    countData = round(cntMat),
    colData = metaDat,
    design = ~ condition + type - 1
)
dds_Full <- DESeq(dds_Full)
resultsNames(dds_Full)
resTrtUnt <- results(
    dds_Full,
    contrast = list(
        c("conditiontreated"),
        c("conditionuntreated")
    )
)
head(resTrtUnt)
resEndBias <- results(
    dds_Full,
    name = "typesingle.read"
)
head(resEndBias[order(resEndBias$padj), ])


#########################################
## Design matrix and contrasts - ZBTB4 ##
#########################################
# Ratio-of-ratios test for RNA-seq response #
resultsNames(dds_RNA)
resResponseFC <- results(
    dds_RNA,
    contrast = list(
        c("TreatmentRNA_DOX_p402", "TreatmentRNA_DMSO_p401"),
        c("TreatmentRNA_DMSO_p402", "TreatmentRNA_DOX_p401")
    )
)
head(resResponseFC)
plotMA(
    resResponseFC,
    main = "ZBTB3 - RNA relative treatment response"
)


##########################
## Shrinkage correction ##
##########################
# ashr contrast treated against untreated #
shrunkTrtUnt <- lfcShrink(
    dds_Full,
    contrast = list(
        c("conditiontreated"),
        c("conditionuntreated")
    ),
    type = "ashr",
    svalue = TRUE
)
head(shrunkTrtUnt)

# Check genes with zeros in one group #
zeroGenes <- which(rowSums(cntMat[, 1:3]) == 0 | rowSums(cntMat[, 4:7]) == 0)
shrunkTrtUnt[head(zeroGenes), ]
cntMat[head(zeroGenes), ]

# apeglm untreated against treated #
apeUntTrt <- lfcShrink(
    dds_quickStart,
    coef = "condition_untreated_vs_treated",
    type = "apeglm",
    svalue = TRUE
)
head(apeUntTrt)

# apeglm treated against untreated #
metaDat$condition <- relevel(metaDat$condition, "untreated")
dds_Relevel <- DESeqDataSetFromMatrix(
    countData = round(cntMat),
    colData = metaDat,
    design = ~ condition
)
dds_Relevel <- DESeq(dds_Relevel)
resultsNames(dds_Relevel)
apeTrtUnt <- lfcShrink(
    dds_Relevel,
    coef = "condition_treated_vs_untreated",
    type = "apeglm",
    svalue = TRUE
)
head(apeTrtUnt)

# ashr vs apeglm #
plot(
    x = -log10(shrunkTrtUnt$svalue),
    y = -log10(apeTrtUnt$svalue),
    pch = 19, cex = 0.25, 
    xlab = "apeglm s-value (-log10)",
    ylab = "ashr s-value (-log10)",
    xlim = c(0, 50), ylim = c(0, 50)
)
abline(0, 1, col = 2, lty = 2)


