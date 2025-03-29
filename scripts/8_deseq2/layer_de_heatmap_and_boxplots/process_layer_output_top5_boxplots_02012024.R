# see "~/Downloads/process_layer_output_top5_heatmap_02012024.R

library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(matrixStats)
library(corrplot)
library(testit)
library(ggdendro)

set.seed(101)

# contsOnly <- FALSE
contsOnly <- TRUE

# lognormalization
lognorm_const <- 1
# lognorm_const <- 0.1
scale_factor <- 10^4


agg_path <- "/Users/jy1008/Downloads/pseudobulk_DE_DESeq2_V4_12052023/agg_counts"
clinical_data <- read.csv("/Users/jy1008/Downloads/pseudobulk_DE_DESeq2_V4_12052023/all_samples_meta_data_V4_01042024.csv")
clinical_data$X <- NULL

# in nctx_temporal, there is only 1 sample with a score of 4: set this to 3 so as not to overly bias the regression
clinical_data$nctx_temporal_fix <- clinical_data$nctx_temporal
clinical_data$nctx_temporal_fix[clinical_data$nctx_temporal_fix > 3] <- 3


st_layers <- c(paste0("Layer", 1:6), "WM")
# sig_genes <- c("SNCA", "SNCG", "BIN3", "CARTPT", "RORB") # visual check
sig_genes <- c(
"FABP7", "MSX2", "RELN", "AL031056.1", "LINC01094", "CD74", "AHNAK", "DTNA", "AQP4", "SPARC",
"SOWAHA", "HPCAL1", "PDZD2", "LAMP5", "NECAB2", "LRATD1", "DGKB", "KCNIP2", "ITPKA", "SIPA1L1",
"SMPX", "NEFM", "TPGS2", "GPX3", "COL5A2", "GPR26", "NGB", "DCC", "SCN1B", "SCN1A",
"TSHZ2", "PLCXD3", "RORB", "TPBG", "FRMPD2", "HPCA", "PTPRT", "PKD2L1", "NREP", "PLCXD2",
"HTR2C", "GRIN3A", "CD24", "CAMK2D", "SMYD2", "TUBB2A", "STX1A", "KIF17", "IFITM10", "UHMK1",
"SLC35F1", "B3GALT2", "DIRAS2", "RASGRP1", "AK5", "TMEM233", "MREG", "ISLR", "NR4A2", "THEMIS",
"AQP1", "SLC7A14-AS1", "PRR18", "CERCAM", "CTNNA3", "LDB3", "MBP", "MOBP", "BCAS1", "TP53INP2")

layer_genes <- c()
for(st_layer in st_layers) {
    agg_file <- sprintf("agg_counts_filtered_smoothed_label_s5_%s_V4_12052023.csv", st_layer)
    layer_df <- read.csv(file.path(agg_path, agg_file), row.names = 1)
    if(contsOnly) {
        clin_samples <- clinical_data[clinical_data$diagnosis == "Control", "sample_name"]
        clin_samples <- intersect(clin_samples, colnames(layer_df))
        layer_df <- layer_df[, clin_samples]
    }
    print(paste("Layer:", st_layer, "Controls Only:", contsOnly, "Num Samples:", ncol(layer_df)))
    
    # convert counts to lognormalized expression
    for(missing_gene in setdiff(sig_genes, rownames(layer_df))) { # insert missing genes with 0 counts
        layer_df[missing_gene, ] <- 0
    }
    layer_df <- sweep(layer_df, 2, colSums(layer_df),`/`)
    layer_df <- log(layer_df * scale_factor + lognorm_const)
    matching_sig_genes <- intersect(sig_genes, rownames(layer_df))
    layer_df <- layer_df[matching_sig_genes, ]
    assert(nrow(layer_df) == length(sig_genes))

    layer_df <- as.data.frame(t(layer_df))
    layer_df$sample_name <- rownames(layer_df)
    layer_df <- layer_df %>% pivot_longer(!sample_name, names_to = "gene", values_to = "lognorm_expr") %>% as.data.frame()
    layer_df$layer <- st_layer
    layer_genes <- rbind(layer_genes, layer_df)
}

pdf(sprintf("pseudobulk_layer_deseq2_V4_top10_contsonly%s_log%sp_boxplots_02012024.pdf", contsOnly, lognorm_const))
for(sig_gene in sig_genes) {
    layer_gene <- layer_genes[layer_genes$gene == sig_gene, ]
    layer_gene <- merge(layer_gene, clinical_data[, c("sample_name", "diagnosis")], by="sample_name")
    layer_gene$layer <- factor(layer_gene$layer, levels=st_layers) # make all layers show up, even if empty
    pt <- ggplot(layer_gene, aes(x=layer, y=lognorm_expr, fill=as.factor(diagnosis))) + geom_boxplot() + theme_classic() + ylab(sig_gene)
    print(pt)
}
dev.off()
