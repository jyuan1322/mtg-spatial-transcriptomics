library(stringr)
library(dplyr)
library(tidyr)
library(testit)
library(ggplot2)

# de_genes <- read.csv("/mnt/data0/projects/ASAP/ST/clustering_pipeline/pseudobulk_DE_DESeq2_V4_12052023/pseudobulk_clinical_deseq2_all_summary_stats_12122023.csv", row.names=1)


# input table should be
# gene log2fc pval padj comparison layer

# "","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"
de_dir <- "/mnt/data0/projects/ASAP/ST/clustering_pipeline/pseudobulk_DE_DESeq2_V4_12052023"
layers <- c(paste("Layer", 1:6), "WM")
file_list <- NULL
# layer_vs_layer
# for(trait in c("layer_vs_layer")) {
for(trait in c("nctx_temporal_fix", "brain_stem_sn", "last_mmse_test_score", "motor_updrs_score")) {
    for(layer in layers) {
        layer_nospace <- gsub(' ', '', layer)
        # layer_vs_layer
        # filename <- file.path(de_dir, trait, sprintf("pseudobulk_layer_deseq2_is_%s_contsonlyFALSE_12052023.csv", layer_nospace))
        filename <- file.path(de_dir, trait, sprintf("pseudobulk_clinical_deseq2_%s_%s_12052023.csv", trait, layer_nospace))
        if(file.exists(filename)) {
            file_list <- rbind(file_list, data.frame(trait=trait, layer=layer, filename=filename))
        } else {
            print(paste("File doesn't exist:", filename))
        }
    }
}

dotplot_df <- NULL
for(i in rownames(file_list)) {
    row <- file_list[i, ]
    results <- read.csv(row$filename)
    colnames(results)[colnames(results) == "X"] <- "gene"
    results$trait <- row$trait
    results$layer <- row$layer
    # results$nlog10_pval <- -log10(results$pvalue)

    dotplot_df <- rbind(dotplot_df, results)
}
# add nlog10 and set padj=NULL to padj=1
dotplot_df$nlog10_pval <- -log10(dotplot_df$pvalue)
max_val <- max(dotplot_df[!is.infinite(dotplot_df$nlog10_pval), ]$nlog10_pval)
inf_subst <- round(max_val + 50, -2)
if ( any(is.infinite(dotplot_df$nlog10_pval)) ) {
    dotplot_df[is.infinite(dotplot_df$nlog10_pval), ]$nlog10_pval <- inf_subst
}
dotplot_df[is.na(dotplot_df$padj), "padj"] <- 1

# keep genes which are significant in at least one comparison
keep_genes <- unique(dotplot_df[dotplot_df$padj < 0.05, "gene"])
dotplot_df <- dotplot_df[dotplot_df$gene %in% keep_genes, ]




dotplot_gene_order <- function(dotplot_df, names_col, outname) {
    # set row similarity using only significant values
    dotplot_df$log2FoldChange_sigonly <- dotplot_df$log2FoldChange
    dotplot_df[dotplot_df$padj >= 0.05, "log2FoldChange_sigonly"] <- 0
    # gene_table_wide <- pivot_wider(dotplot_df, id_cols="gene", names_from=all_of(names_col), values_from="log2FoldChange")
    gene_table_wide <- pivot_wider(dotplot_df, id_cols="gene", names_from=all_of(names_col), values_from="log2FoldChange_sigonly")
    gene_table_wide <- as.data.frame(gene_table_wide)
    rownames(gene_table_wide) <- gene_table_wide$gene
    gene_table_wide$gene <- NULL
    gene_table_wide[is.na(gene_table_wide)] <- 0 # to handle distance function
    
    distance.row <- dist(gene_table_wide, method = 'euclidean') # the default
    cluster.row <- hclust(distance.row, method = 'complete') # the default
    
    pdf(sprintf("pseudobulk_clinical_deseq2_gene_dendro_%s.pdf", outname))
    # pdf(sprintf("pseudobulk_clinical_deseq2_gene_dendro_%s_V2_09202023.pdf", clinical_feat))
    order_score <- rowSums(gene_table_wide * c(10^6, 10^5, 10^4, 10^3, 10^2, 10, 1))
    reordered.cluster.row <- reorder(as.dendrogram(cluster.row), order_score)
    plot(reordered.cluster.row, horiz=TRUE)
    # plot(as.dendrogram(cluster.row), horiz=TRUE)
    dev.off()
    # gene_display_order <- cluster.row$label[cluster.row$order]
    gene_display_order <- rownames(gene_table_wide)[order.dendrogram(reordered.cluster.row)]
    return(gene_display_order)
}

degene_dotplot <- function(dotplot_df, outname, width=6, height=1000) {
    dotplot_df$comparison_vis <- paste(dotplot_df$trait, dotplot_df$layer)
    comparison_order <- sort(unique(dotplot_df$comparison_vis))
    
    
    gene_display_order <- dotplot_gene_order(dotplot_df, "comparison_vis", outname)
    
    # width <- 6
    # height <- 1000
    
    pdf(sprintf("pseudobulk_clinical_deseq2_dotplot_effects_%s.pdf", outname), width=width, height=height)
    dt_pt <- ggplot(dotplot_df, aes(x=factor(comparison_vis, level=comparison_order), y = factor(gene, level = gene_display_order), color = log2FoldChange, size = nlog10_pval)) +
    geom_point(shape = 21,  stroke = 1, aes(colour = padj < 0.05, fill = log2FoldChange)) +
    # scale_fill_gradient2(low = "#998ec3", mid = "white", high = "#f1a340", midpoint = 0, name = 'log2FoldChange') +
    scale_fill_gradient2(low = "#8f00ff", mid = "white", high = "#ffa500", midpoint = 0, name = 'log2FoldChange') +
    scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = rgb(0,0,0,0)), name = "p_adj < 0.05") +
    cowplot::theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_y_discrete(label = function(x) stringr::str_trunc(x, 25))
    print(dt_pt)
    dev.off()
}
outname <- "gene_alltraits_12192023"

lfc_thresh <- 0.8
max_lfc_per_gene <- dotplot_df %>%
  group_by(gene) %>%
  summarise(max = max(log2FoldChange, na.rm=TRUE))
df_gene_list <- max_lfc_per_gene[max_lfc_per_gene$max > lfc_thresh, "gene"] %>% pull(gene)
dotplot_df_filtered <- dotplot_df[dotplot_df$gene %in% df_gene_list, ]

# layer_vs_layer
# outname <- sprintf("gene_layer_vs_layer_lfcthresh%s_01112024", lfc_thresh)

degene_dotplot(dotplot_df_filtered, outname, height=1000)
