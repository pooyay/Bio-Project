setwd("D:/University/Term 7/Bio/Project/Bio-Project/")
library(GEOquery)
library(limma)
library(Biobase)
library(umap)
library(pheatmap)
library(gplots)
library(ggplot2)
library(reshape2)
library(plyr)
library(Rtsne)

Sys.setenv("VROOM_CONNECTION_SIZE" = 100 * 1000 * 1000)
GEOseries <- "GSE48558"

gset <- getGEO(GEOseries, GSEMatrix =TRUE, AnnotGPL=TRUE, destdir = "Data/")
if (length(gset) > 1) idx <- grep("GPL6244", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

#gsms <- paste0("1111111111111XXXXXXXXXXXXXXXXXXXXXXXXXXX0XXX0XXXXX",
#               "XXXXXXXXXXXXXXXXXX0X0XXX0X0000X0XX00XX00X0X0X0X0X0",
#               "XXX0XXX0XXXXXXXXXXXXXXXXXXXXXXXXXXXXX0000000110!111",
#               "00000000000000000000")
#sml <- strsplit(gsms, split="")[[1]]

#sel <- which(sml != "X")
#sml <- sml[sel]
#gset <- gset[,sel]

gr <- c(rep("aml",13), rep("X",27), "healthy", rep("X",3), "healthy",
        rep("X",23), "healthy", "X", "healthy", rep("X",3), "healthy",
        "X", rep("healthy",4), "X", "healthy", rep("X",2), rep("healthy",2),
        rep("X",2), rep("healthy",2), "X", "healthy", "X", "healthy", "X",
        "healthy", "X", "healthy", "X", "healthy", rep("X",3), "healthy",
        rep("X",3), "healthy", rep("X",29), rep("healthy",7), rep("aml",2),
        "healthy", rep("aml",3), rep("healthy", 20))


ex <- exprs(gset)
# Max of ex is not too large (13.76154). So logarithm is not necessary.
# Normalizing is not necessary for this dataset.

# **** Boxplot of expression matrix ****
pdf("Results/boxplot.pdf", width=170)
boxplot(ex)
dev.off()

# **** Correlation heatmap of expression matrix ****
pdf("Results/CorHeatmap.pdf", width=170, height=170)
pheatmap(cor(ex), labels_row=gr, labels_col=gr, color=bluered(256), border_color=NA)
dev.off()

# **** Dimensionality Reduction ****
# Principal Component Analysis:
pc <- prcomp(ex)
pdf("Results/PC.pdf")
plot(pc)
plot(pc$x[,1:2])
dev.off()
# As the codes above show, PC1 does not give us useful information.
# We try to ???? 

ex.scale <- t(scale(t(ex), scale=FALSE))
pc <- prcomp(ex.scale)
pdf("Results/PC_scaled.pdf")
plot(pc)
plot(pc$x[,1:2])
dev.off()
# Now more than one principal components have information and PC1 is not
# the only principal components that varies genes.

# tSNE:
tsne_results <- Rtsne(ex, perplexity=30, check_duplicates = FALSE)
pdf("Results/tSNE.pdf")
plot(tsne_results$Y)
dev.off()

pcr <- data.frame(pc$r[,1:3], Group=gr)
pdf("Results/PC_samples.pdf")
ggplot(pcr, aes(PC1, PC2, color=Group)) + geom_point(size=3) + theme_bw()
dev.off()

# **** Differential Expression Analysis ****
#gr_factor <- factor(gr)
gr <- factor(gr)
gset$description <- gr
design_matrix <- model.matrix()





