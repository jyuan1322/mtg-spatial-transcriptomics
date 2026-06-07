library(Seurat)

# BN2003
path1 <- "/mnt/data/projects/ASAP/Xenium/xenium/20240530__205709__MTG_AddON_30may2024/output-XETG00322__0038488__Region_1__20240530__205717"
xenium.obj1 <- LoadXenium(path1, fov = "fov")
# remove cells with 0 counts
xenium.obj1 <- subset(xenium.obj1, subset = nCount_Xenium > 0)

# BN1730
path2 <- "/mnt/data/projects/ASAP/Xenium/xenium/20240530__205709__MTG_AddON_30may2024/output-XETG00322__0038489__Region_1__20240530__205717"
xenium.obj2 <- LoadXenium(path2, fov = "fov")
# remove cells with 0 counts
xenium.obj2 <- subset(xenium.obj2, subset = nCount_Xenium > 0)




# get gene names
# rownames(xenium.obj@assays$Xenium@counts)

fig4_genes1 <- c("RELN", "AQP4", "LAMP5", "NEFM", "TSHZ2")
fig4_genes2 <- c("HTR2C", "RASGRP1", "NR4A2", "THEMIS", "CTNNA3")
# fig4_genes <- c("RELN", "AQP4", "LAMP5", "NEFM", "HTR2C", "RASGRP1", "NR4A2", "THEMIS")

# fig4_genes_layer <- c("RELN", "AQP4", "LAMP5", "NEFM", "TSHZ2", "HTR2C", "RASGRP1", "NR4A2", "THEMIS", "CTNNA3")
# fig4_genes_layer <- c("AQP4", "LAMP5", "NEFM", "TSHZ2", "HTR2C", "RASGRP1", "CTNNA3")
fig4_genes_layer <- c("AQP4", "LAMP5", "NEFM", "TSHZ2", "HTR2C", "NR4A2", "MOBP")

pdf("xenium_test1_09242024.pdf")
ImageDimPlot(xenium.obj1, fov = "fov", molecules = fig4_genes1, nmols = 20000)
ImageFeaturePlot(xenium.obj1, features = fig4_genes1, size = 0.75, cols = c("white", "red"))
ImageDimPlot(xenium.obj1, fov = "fov", molecules = fig4_genes2, nmols = 20000)
ImageFeaturePlot(xenium.obj1, features = fig4_genes2, size = 0.75, cols = c("white", "red"))
dev.off()
pdf("xenium_test2_09242024.pdf")
ImageDimPlot(xenium.obj2, fov = "fov", molecules = fig4_genes1, nmols = 20000)
ImageFeaturePlot(xenium.obj2, features = fig4_genes1, size = 0.75, cols = c("white", "red"))
ImageDimPlot(xenium.obj2, fov = "fov", molecules = fig4_genes2, nmols = 20000)
ImageFeaturePlot(xenium.obj2, features = fig4_genes2, size = 0.75, cols = c("white", "red"))
dev.off()

pdf("xenium_test1_09242024.pdf")
ImageDimPlot(xenium.obj1, fov = "fov", molecules = fig4_genes1, nmols = 20000)
ImageFeaturePlot(xenium.obj1, features = fig4_genes1, size = 0.75, cols = c("white", "red"))
ImageDimPlot(xenium.obj1, fov = "fov", molecules = fig4_genes2, nmols = 20000)
ImageFeaturePlot(xenium.obj1, features = fig4_genes2, size = 0.75, cols = c("white", "red"))
dev.off()

random_seed <- 101
set.seed(101)
pdf("xenium1_MTG_layer_09242024.pdf")
ImageDimPlot(xenium.obj1, fov = "fov", molecules = fig4_genes_layer, nmols = 10000, mols.size=0.1, mols.alpha=1.0, coord.fixed=TRUE, axes=TRUE, border.size=NA)
dev.off()
pdf("xenium2_MTG_layer_09242024.pdf")
ImageDimPlot(xenium.obj2, fov = "fov", molecules = fig4_genes_layer, nmols = 10000, mols.size=0.1, mols.alpha=1.0, coord.fixed=TRUE, axes=TRUE, border.size=NA)
dev.off()
