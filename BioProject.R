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

#real_gset <- gset
#gset <- real_gset

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

ex_norm <- normalizeQuantiles(ex)

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

pcr <- data.frame(pc$r[,1:3], Group=gr)
pdf("Results/PC_samples.pdf")
ggplot(pcr, aes(PC1, PC2, color=Group)) + geom_point(size=3) + theme_bw()
dev.off()


# MDS:
dist_matrix <- dist(t(ex))
ex_mds <- cmdscale(dist_matrix)
pdf("Results/MDS.pdf")
plot(ex_mds[,1], ex_mds[,2])
dev.off()

ex_mds_dataframe <- data.frame(ex_mds[,1:2], Group=gr)
pdf("Results/MDS_samples.pdf")
ggplot(ex_mds_dataframe, aes(X1, X2, color=Group)) + geom_point(size=3) + theme_bw()
dev.off()


# tSNE:
tsne_results <- Rtsne(ex, perplexity=30, check_duplicates = FALSE)
pdf("Results/tSNE.pdf")
plot(tsne_results$Y)
dev.off()

#tSNE_dataframe <- data.frame(tsne_results$Y[,1:2], Group=gr)
#pdf("Results/tSNE_samples.pdf")
#plot(tSNE_dataframe)
#dev.off()



# **** Differential Expression Analysis ****

#gr_factor <- factor(gr)
gr <- factor(gr)
gset$description <- gr
design <- model.matrix(~ description, gset)
colnames(design) <- levels(gr)


fit <- lmFit(gset, design)
cont_matrix <- makeContrasts(aml - healthy, levels=design)
fit2 <- contrasts.fit(fit, cont_matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
tT <- subset(tT, select=c("Gene.symbol", "Gene.ID", "adj.P.Val", "logFC"))
write.table(tT, "Results/AML_vs_Healthy.txt", row.names=FALSE, sep='\t', quote=FALSE)

aml.up <- subset(tT, logFC > 1 & adj.P.Val < 0.05)
aml.up.genes <- unique(as.character(strsplit2(aml.up$Gene.symbol, "///")))
write.table(aml.up.genes, file="Results/AML_vs_Healthy_UP.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
aml.down <- subset(tT, logFC < -1 & adj.P.Val < 0.05)
aml.down.genes <- unique(as.character(strsplit2(aml.down$Gene.symbol, "///")))
write.table(aml.down.genes, file="Results/AML_vs_Health_Down.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)












