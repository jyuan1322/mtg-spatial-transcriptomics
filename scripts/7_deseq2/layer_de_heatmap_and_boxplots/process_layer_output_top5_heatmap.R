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

# grab the top N DE genes from each layer
topN <- 10
# contsOnly <- FALSE
contsOnly <- TRUE

# lognormalization
lognorm_const <- 1
# lognorm_const <- 0.1
scale_factor <- 10^4

st_layers <- c(paste0("Layer", 1:6), "WM")

agg_path <- "/Users/jy1008/Downloads/pseudobulk_DE_DESeq2_V4_12052023/agg_counts"
clinical_data <- read.csv("/Users/jy1008/Downloads/pseudobulk_DE_DESeq2_V4_12052023/all_samples_meta_data_V4_01042024.csv")
clinical_data$X <- NULL

# in nctx_temporal, there is only 1 sample with a score of 4: set this to 3 so as not to overly bias the regression
clinical_data$nctx_temporal_fix <- clinical_data$nctx_temporal
clinical_data$nctx_temporal_fix[clinical_data$nctx_temporal_fix > 3] <- 3

layer_path <- "/Users/jy1008/Downloads/pseudobulk_DE_DESeq2_V4_12052023/layer_vs_layer"
sig_genes <- c()
for(st_layer in st_layers) {
    de_file <- sprintf("pseudobulk_layer_deseq2_is_%s_contsonly%s_12052023.csv", st_layer, contsOnly)
    de_genes <- read.csv(file.path(layer_path, de_file))
    colnames(de_genes)[colnames(de_genes) == "X"] <- "gene"
    
    # make sure there are no NA p-adj values
    # assert(sum(is.na(de_genes$padj)) == 0)
    de_genes[is.na(de_genes$padj), "padj"] <- 1.0

    # keep only positively-associated genes
    de_genes_sig <- de_genes[de_genes$log2FoldChange > 0, ]
    de_genes_sig <- de_genes_sig[order(de_genes_sig$padj), ]
    de_genes_sig$layer <- st_layer
    # get topN x 3 sig genes just so we have some substitutes when gene duplicates are removed
    sig_genes <- rbind(sig_genes, de_genes_sig[1:(topN*3), c("gene", "pvalue", "padj", "layer")])
}
# for genes which appear in multiple layers, only pick the layer with the most significant pvalue
# NOTE: sorting by padj can result in ties if small padjs are set to 0
sig_genes <- sig_genes %>% group_by(gene) %>% slice(which.min(pvalue)) %>% as.data.frame()
# arrange() smallest to largest
sig_genes <- sig_genes %>% arrange(pvalue) %>% group_by(layer) %>% slice(1:topN) %>% as.data.frame()

# NOTE: do we want to keep genes which are exclusively upregulated in one layer? This produces very few genes.

if(contsOnly) {
    clinical_labels <- c("Control")
} else {
    clinical_labels <- c("Control", "ILBD", "Case")
}
layer_de_heatmap <- c()
for(st_layer in st_layers) {
    agg_file <- sprintf("agg_counts_filtered_smoothed_label_s5_%s_V4_12052023.csv", st_layer)
    layer_df <- read.csv(file.path(agg_path, agg_file), row.names = 1)

    # convert counts to lognormalized expression
    layer_df <- sweep(layer_df, 2, colSums(layer_df),`/`)
    layer_df <- log(layer_df * scale_factor + lognorm_const)
    matching_sig_genes <- intersect(sig_genes$gene, rownames(layer_df))
    layer_df <- layer_df[matching_sig_genes, ]
    assert(nrow(layer_df) == length(matching_sig_genes))
    
    # calculate means according to clinical feature

    for(clinical_label in clinical_labels) {
        clin_samples <- clinical_data[clinical_data$diagnosis == clinical_label, "sample_name"]
        clin_samples <- intersect(clin_samples, colnames(layer_df))
        layer_df_sub <- layer_df[, clin_samples]
        layer_df_sub <- as.data.frame(rowMeans(layer_df_sub))
        layer_df_sub$gene <- rownames(layer_df_sub)
        layer_df_sub$group <- paste(st_layer, clinical_label, sep='_')
        rownames(layer_df_sub) <- 1:nrow(layer_df_sub)
        colnames(layer_df_sub)[1] <- "mean_lognorm_count"
        layer_de_heatmap <- rbind(layer_de_heatmap, layer_df_sub)
    }
}

# NOTE: if using lognormalization constants other than 1 (i.e. log0.1p), the min value is not zero but rather log(lognorm_const)
# layer_de_heatmap <- layer_de_heatmap %>% pivot_wider(names_from = group, values_from = mean_lognorm_count, values_fill = 0) %>% as.data.frame()
layer_de_heatmap <- layer_de_heatmap %>% pivot_wider(names_from = group, values_from = mean_lognorm_count, values_fill = log(lognorm_const)) %>% as.data.frame()
assert(length(unique(layer_de_heatmap$gene)) == nrow(layer_de_heatmap))
rownames(layer_de_heatmap) <- layer_de_heatmap$gene
layer_de_heatmap$gene <- NULL

# reorder columns based on diagnosis
col_order_diagnosis <- c()
for(st_layer in st_layers) {
    for(clinical_label in clinical_labels) {
        col_order_diagnosis <- c(col_order_diagnosis, paste(st_layer, clinical_label, sep='_'))
    }
}
layer_de_heatmap <- layer_de_heatmap[, col_order_diagnosis]

# standard normalize each row
# scale() standard normalizes by column
layer_de_heatmap_scaled <- t( scale( t(layer_de_heatmap) ))

# reorder rows (genes) according to dendrogram
gene_row_order <- c()
pdf(sprintf("pseudobulk_layer_deseq2_V4_top%s_contsonly%s_log%sp_dendros_02012024.pdf", topN, contsOnly, lognorm_const))
for(st_layer in sort(unique(sig_genes$layer))) {
    sig_genes_sub <- sig_genes[sig_genes$layer == st_layer, "gene"]
    layer_de_sub <- layer_de_heatmap_scaled[sig_genes_sub, ]
    distance.row <- dist(layer_de_sub, method = 'euclidean') # the default
    cluster.row <- hclust(distance.row, method = 'complete') # the default
    dendro_pt <- ggdendrogram(cluster.row, rotate=TRUE, size=20) + ggtitle(st_layer)
    print(dendro_pt)
    gene_row_order <- c( gene_row_order, rownames(layer_de_sub)[order.dendrogram(as.dendrogram(cluster.row))] )
}
dev.off()

layer_de_heatmap_scaled <- layer_de_heatmap_scaled[gene_row_order, ]

col1<-colorRampPalette(c('#0571b0', '#92c5de', '#f7f7f7', '#f4a582', '#ca0020'))
col_limits <- c( floor(min(layer_de_heatmap_scaled)*10)/10,
                 ceiling(max(layer_de_heatmap_scaled)*10)/10 )
num_increments <- 100
pdf(sprintf("pseudobulk_layer_deseq2_V4_top%s_contsonly%s_log%sp_02012024.pdf", topN, contsOnly, lognorm_const), width=20, height=50)
corrplot(as.matrix(layer_de_heatmap_scaled), is.corr=FALSE, method='color', col=col1(num_increments), col.lim=c(-max(col_limits), max(col_limits)) )
dev.off()

