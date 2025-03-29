# see "~/Downloads/process_layer_output_top5_heatmap_02012024.R

library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(matrixStats)
library(corrplot)
library(testit)
library(ggdendro)
library(DESeq2)
library(ggrepel)

set.seed(101)

# contsOnly <- FALSE
contsOnly <- TRUE

# lognormalization
lognorm_const <- 1
# lognorm_const <- 0.1
# scale_factor <- 10^4

st_layers <- c(paste0("Layer", 1:6), "WM")

agg_path <- "/Users/jy1008/Downloads/pseudobulk_DE_DESeq2_V4_12052023/agg_counts"
clinical_data <- read.csv("/Users/jy1008/Downloads/pseudobulk_DE_DESeq2_V4_12052023/all_samples_meta_data_V4_01042024.csv")
clinical_data$X <- NULL

# in nctx_temporal, there is only 1 sample with a score of 4: set this to 3 so as not to overly bias the regression
clinical_data$nctx_temporal_fix <- clinical_data$nctx_temporal
clinical_data$nctx_temporal_fix[clinical_data$nctx_temporal_fix > 3] <- 3

layer_path <- "/Users/jy1008/Downloads/pseudobulk_DE_DESeq2_V4_12052023/layer_vs_layer"
# pdf(sprintf("gsea_kegg_MAplot_layer_vs_layer_log%sp_lfcshrink_02082024.pdf", lognorm_const))
pdf(sprintf("gsea_kegg_MAplot_layer_vs_layer_log%sp_lfcshrink_04292024.pdf", lognorm_const))
for(st_layer in st_layers) {
    de_file <- sprintf("pseudobulk_layer_deseq2_is_%s_contsonly%s_12052023.csv", st_layer, contsOnly)
    de_genes <- read.csv(file.path(layer_path, de_file))
    colnames(de_genes)[colnames(de_genes) == "X"] <- "gene"
    
    # make sure there are no NA p-adj values
    # assert(sum(is.na(de_genes$padj)) == 0)
    de_genes[is.na(de_genes$padj), "padj"] <- 1.0
    
    # lfcshrink
    de_rds <- sprintf("pseudobulk_layer_deseq2_dds_lrt_is_%s_contsonly%s_12052023.rds", st_layer, contsOnly)
    dds_lrt <- readRDS(file.path(layer_path, de_rds))
    res <- results(dds_lrt)
    lfc_coef <- length(resultsNames(dds_lrt)) # take the last one
    resLFC.ape <- lfcShrink(dds_lrt, coef = lfc_coef, res = res, type = "apeglm")
    assert(all(de_genes$gene == resLFC.ape$gene))
    de_genes$log2FoldChange_shrink <- resLFC.ape$log2FoldChange
    
    # ribosomal category
    de_genes$ribo_cat <- "other"
    de_genes[startsWith(de_genes$gene, "RPL"), "ribo_cat"] <- "RPL"
    de_genes[startsWith(de_genes$gene, "RPS"), "ribo_cat"] <- "RPS"
    
    de_genes$logp_basemean <- log(de_genes$baseMean + lognorm_const)
    
    pt <- ggplot(de_genes, aes(x=logp_basemean, y=log2FoldChange_shrink, color=ribo_cat)) +
        geom_point(data=de_genes[de_genes$ribo_cat == "other", ], size=1, aes(color="other")) +
        geom_point(data=de_genes[de_genes$ribo_cat == "RPL", ], size=1, aes(color="RPL")) +
        # geom_text_repel(data=de_genes[de_genes$ribo_cat == "RPL", ], aes(x=logp_basemean, y=log2FoldChange_shrink, label=gene), max.overlaps=20) +
        geom_point(data=de_genes[de_genes$ribo_cat == "RPS", ], size=1, aes(color="RPS")) +
        # geom_text_repel(data=de_genes[de_genes$ribo_cat == "RPS", ], aes(x=logp_basemean, y=log2FoldChange_shrink, label=gene), max.overlaps=20) +
        geom_line(aes(y = 0), col = "black") + theme_classic() +
        scale_color_manual(values = c("gray90", "#377eb8", "#4daf4a")) + ggtitle(paste("Ribosome", st_layer))
    print(pt)
    
    # proteasome
    de_genes$proteasome <- "other"
    de_genes[startsWith(de_genes$gene, "PSM"), "proteasome"] <- "proteasome"
    
    pt <- ggplot(de_genes, aes(x=logp_basemean, y=log2FoldChange_shrink, color=proteasome)) +
        geom_point(data=de_genes[de_genes$proteasome == "other", ], size=1, aes(color="other")) +
        geom_point(data=de_genes[de_genes$proteasome == "proteasome", ], size=1, aes(color="proteasome")) +
        # geom_text_repel(data=de_genes[de_genes$proteasome == "proteasome", ], aes(x=logp_basemean, y=log2FoldChange_shrink, label=gene), max.overlaps=20) +
        geom_line(aes(y = 0), col = "black") + theme_classic() +
        scale_color_manual(values = c("gray90", "#4daf4a")) + ggtitle(paste("Proteasome", st_layer))
    print(pt)
    
    # translational proteins
    trans_protein_list <- c("ZBP1", "FMR1", "DKC1", "EIF4E", "STAU1", "STAU2")
    de_genes$trans_protein <- "other"
    de_genes[de_genes$gene %in% trans_protein_list, "trans_protein"] <- "translational_protein"

    pt <- ggplot(de_genes, aes(x=logp_basemean, y=log2FoldChange_shrink, color=trans_protein)) +
        geom_point(data=de_genes[de_genes$trans_protein == "other", ], size=1, aes(color="other")) +
        geom_point(data=de_genes[de_genes$trans_protein == "translational_protein", ], size=1, aes(color="translational_protein")) +
        # geom_text_repel(data=de_genes[de_genes$trans_protein == "translational_protein", ], aes(x=logp_basemean, y=log2FoldChange_shrink, label=gene), max.overlaps=20) +
        geom_line(aes(y = 0), col = "black") + theme_classic() +
        scale_color_manual(values = c("gray90", "#4daf4a")) + ggtitle("Translational Protein", st_layer)
    print(pt)

}
dev.off()
