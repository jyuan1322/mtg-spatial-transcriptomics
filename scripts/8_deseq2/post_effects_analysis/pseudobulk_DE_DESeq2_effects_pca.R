# see "pseudobulk_DE_DESeq2_post_analysis_subcats_02142023.R

# https://lashlock.github.io/compbio/R_presentation.html
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#differential-expression-analysis

library(Seurat)
library(DESeq2)
library(stringr)
library(testit)
library(ggplot2)
library(dplyr)
library(tidyr)
library(umap)



# GSEA
# https://github.com/sterding/BRAINcode/blob/7b4c5e816ff1cf9af86041326b71cf3f3e2e4bf6/modules/_pathway.R

# run on all cases vs controls
library(dplyr)
library(tibble)
library(fgsea)
library(ggplot2)

source("/Users/jy1008/Downloads/pseudobulk_DE_DESeq2_GSEA_12292023_utils.R")

run_gsea <- function(de_files_list, width, height, clinical_feat, outname, input_gmt_file, nes_cutoff=0, Q_CUTOFF=NULL) {
    combined_gt <- c()
    for(st_layer in names(de_files_list)) {
        print(sprintf("%s %s", clinical_feat, st_layer))
        
        results <- read.csv(de_files_list[[st_layer]])
        results$gene <- results$X
        results$X <- NULL
        results$comparison <- sprintf("%s * %s", clinical_feat, gsub(' ', '_', st_layer))

        results$nlog10_pval <- -log10(results$pvalue)
        max_val <- max(results[!is.infinite(results$nlog10_pval), ]$nlog10_pval)
        inf_subst <- round(max_val + 50, -2)
        if ( any(is.infinite(results$nlog10_pval)) ) {
            results[is.infinite(results$nlog10_pval), ]$nlog10_pval <- inf_subst
        }
        results$gsea_rank <- results$nlog10_pval * sign(results$log2FoldChange)
        gsea_stats <- deframe(results[, c("gene", "gsea_rank")])
        
        # Q_CUTOFF <- 0.05
        gt <- data.frame()
        # gmt.file = sprintf("~/Downloads/gsea_msigdb_pathways/%s", i)
        pathways <- gmtPathways(input_gmt_file)
        # eps = 0 gives better p-value estimation below 1e-10
        fgseaRes <- fgsea(pathways, stats=gsea_stats, minSize=10, maxSize=500, nPermSimple=10000, eps=0)
        fgseaRes <- fgseaRes[order(fgseaRes$padj),]  # default is BH-adjusted p-value, same as FDR
        
        if (!is.null(Q_CUTOFF)) {
            fgseaRes <- subset(fgseaRes, padj < Q_CUTOFF)
        }
        # only create plots for sig pathways, even if writing all pathways to file
        fgseaRes_sig <- subset(fgseaRes, padj < 0.05)

        if(nrow(fgseaRes)>0) {
            # plot
            # for(x in fgseaRes_sig$pathway) {p=plotEnrichment(pathways[[x]], stats=gsea_stats) + labs(title=x); print(p);}
          
            # flatten the list of leadingEdge
            fgseaRes$leadingEdge <- vapply(fgseaRes$leadingEdge, paste, collapse = ", ", character(1L))

            gt <- rbind(gt, data.frame(gene_set=basename(input_gmt_file), fgseaRes))
        }
        # }
        # dev.off()
        # write.table(gt, file=file.path(paste0(output_file, '.fGSEA.xls')), sep="\t", quote=F, row.names=F, col.names=T)
        gt$comparison <- sprintf("%s * %s", clinical_feat, gsub(' ', '_', st_layer))
        combined_gt <- rbind(combined_gt, gt)
    }
    # split pathway db from pathway name
    combined_gt$pathway_id <- NA
    if(all(startsWith(combined_gt$pathway, "KEGG")) || all(startsWith(combined_gt$pathway, "REACTOME"))) {
        combined_gt[, c("pathway_id", "pathway")] <- str_split_fixed(combined_gt$pathway, "_", n=2)
    }
    col_order <- c("gene_set", "pathway_id", "pathway", "pval", "padj", "log2err", "ES", "NES", "size", "leadingEdge", "comparison")
    combined_gt <- combined_gt[, col_order]
    output_file <- sprintf("de_genes_DESeq2_pseudobulk_%s_%s", clinical_feat, outname)
    write.table(combined_gt, file=file.path(paste0(output_file, '.fGSEA.xls')), sep="\t", quote=F, row.names=F, col.names=T)

    # create dotplot across layers
    # ignore NA p-values: these are produced by fgsea when estimation is poor
    combined_gt <- combined_gt[!is.na(combined_gt$pval), ]
    
    # keep only genes which are significant in more than one category?
    num_sig_tests_cutoff <- 1
    combined_gt_filt <- combined_gt[combined_gt$padj < 0.05, ]
    pathway_sig_counts <- table(combined_gt_filt$pathway)
    sig_pathways <- names(pathway_sig_counts)[pathway_sig_counts >= num_sig_tests_cutoff]
    combined_gt <- combined_gt[combined_gt$pathway %in% sig_pathways, ]
    
    # keep only genes passing log fold change filter
    # nes_cutoff <- 0
    combined_gt_filt <- combined_gt[abs(combined_gt$NES) > nes_cutoff, ]
    sig_pathways <- unique(combined_gt_filt$pathway)
    combined_gt <- combined_gt[combined_gt$pathway %in% sig_pathways, ]

    # -log10(pval=0) = Inf, which doesn't show up in dotplot
    # set Inf values to equal the highest finite value
    combined_gt$nlog_pval <- -log10(combined_gt$pval)
    combined_gt$nlog_pval[is.infinite(combined_gt$nlog_pval)] <- max(combined_gt$nlog_pval[is.finite(combined_gt$nlog_pval)])
    
    layer_order <- unique(combined_gt$comparison)
    
    # if(nrow(combined_gt) == 0) {
    # if there are fewer than 2 unique pathways, clustering rows will fail
    if(length(unique(combined_gt$pathway)) < 2) {
        print("skip")
        return(NULL)
    }

    # NA padj get rendered with gray border otherwise
    # combined_results[is.na(combined_results$padj), "padj"] <- 1

    
    # set row similarity using only significant values
    pathway_table_wide <- pivot_wider(combined_gt, id_cols="pathway", names_from="comparison", values_from="NES")
    pathway_table_wide <- as.data.frame(pathway_table_wide)
    rownames(pathway_table_wide) <- pathway_table_wide$pathway
    pathway_table_wide$pathway <- NULL
    # set NA NES values in pathway_table_wide to 0
    pathway_table_wide[is.na(pathway_table_wide)] <- 0
    
    pathway_padjs_wide <- pivot_wider(combined_gt, id_cols="pathway", names_from="comparison", values_from="padj")
    pathway_padjs_wide <- as.data.frame(pathway_padjs_wide)
    rownames(pathway_padjs_wide) <- pathway_padjs_wide$pathway
    pathway_padjs_wide$pathway <- NULL
    pathway_padjs_wide <- pathway_padjs_wide < 0.05
    
    print(paste(clinical_feat, st_layer, input_gmt_file))
    distance.row <- dist(pathway_table_wide, method = 'euclidean') # the default
    cluster.row <- hclust(distance.row, method = 'complete') # the default
    
    if(is.na(height)) {
        height <- max(10, nrow(pathway_table_wide) / 5)
    }
    pdf(sprintf("gsea_deseq2_pseudobulk_dendro_%s_%s.pdf", clinical_feat, outname), width=width, height=height)
    # order_score <- rowSums(pathway_table_wide * c(10^6, 10^5, 10^4, 10^3, 10^2, 10, 1))
    order_score <- rowSums(pathway_padjs_wide * c(10^6, 10^5, 10^4, 10^3, 10^2, 10, 1)[1:ncol(pathway_padjs_wide)])
    reordered.cluster.row <- reorder(as.dendrogram(cluster.row), order_score)
    plot(reordered.cluster.row, horiz=TRUE)
    # plot(as.dendrogram(cluster.row), horiz=TRUE)
    dev.off()
    # gene_display_order <- cluster.row$label[cluster.row$order]
    pathway_display_order <- rownames(pathway_table_wide)[order.dendrogram(reordered.cluster.row)]

    print(paste0("number of rows in dotplot: ", nrow(pathway_table_wide)))
    pdf(sprintf("gsea_deseq2_pseudobulk_dotplot_%s_NES%s_%s.pdf", clinical_feat, nes_cutoff, outname), width=width, height=height)
    dt_pt <- ggplot(combined_gt, aes(x = factor(comparison, levels = layer_order), y = factor(pathway, levels = pathway_display_order), color = NES, size = nlog_pval)) +
    geom_point(shape = 21,  stroke = 1, aes(colour = padj < 0.05, fill = NES)) +
    # scale_fill_gradient2(low = "#998ec3", mid = "white", high = "#f1a340", midpoint = 0, name = 'log2FoldChange') +
    # purple, orange
    # scale_fill_gradient2(low = "#8f00ff", mid = "white", high = "#ffa500", midpoint = 0, name = 'NES') +
    scale_fill_gradient2(low = "#0571b0", mid = "#f7f7f7", high = "#ca0020", midpoint = 0, name = 'NES') +
    scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = rgb(0,0,0,0)), name = "p_adj < 0.05") +
    cowplot::theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          panel.grid = element_line(color = "#dddddd", size = 0.75, linetype = 2)) # +
    # scale_y_discrete(label = function(x) stringr::str_trunc(x, 25))
    print(dt_pt)
    dev.off()
}

random_seed <- 101
layers <- c("Layer 1", "Layer 2", "Layer 3", "Layer 4", "Layer 5", "Layer 6", "WM")
effects_path <- "/Users/jy1008/Downloads/pseudobulk_DE_DESeq2_V4_12052023"

# out_label <- 01102024
# out_label <- "01162024"
out_gsea_path <- "/Users/jy1008/Downloads/pseudobulk_DE_GSEA_V4_01102024"

# clinical_features <- c("diagnosis_cont_ilbd",
#                        "diagnosis_cont_case",
#                        "last_mmse_test_score",
#                        "motor_updrs_score",
#                        "nctx_temporal_fix",
#                        "brain_stem_sn")
clinical_features <- c("last_mmse_test_score",
                       "motor_updrs_score",
                       "nctx_temporal_fix")
include_layer_vs <- TRUE



# set.seed(random_seed)

combined_res <- c()
for(clinical_feat in clinical_features) {
    de_files_list <- gsea_clinical_per_layer_filelist(effects_path, clinical_feat, layers)
    # clinical_feat is only used to fill the "comparison" field in output summary files, i.e. "trait * Layer 1"

    for(st_layer in names(de_files_list)) {
        print(sprintf("%s %s", clinical_feat, st_layer))
        
        results <- read.csv(de_files_list[[st_layer]])
        results$gene <- results$X
        results$X <- NULL
        results$comparison <- sprintf("%s * %s", clinical_feat, gsub(' ', '_', st_layer))
        combined_res <- rbind(combined_res, results)
    }
}
if(include_layer_vs) {
    # layer vs layer
    de_files_list <- gsea_layer_vs_layer_filelist(effects_path, layers, contsonly=FALSE)
    for(st_layer in names(de_files_list)) {
        print(sprintf("%s %s", "layer_vs", st_layer))
        results <- read.csv(de_files_list[[st_layer]])
        results$gene <- results$X
        results$X <- NULL
        results$comparison <- sprintf("%s * %s", "layer_vs", gsub(' ', '_', st_layer))
        combined_res <- rbind(combined_res, results)
    }
}
# apply a loose filter for the p-value
pval_filter <- 0.1
combined_res <- combined_res[combined_res$pvalue < pval_filter, ]

# flip the sign of last mmse test score (lower value is later PD)
combined_res[startsWith(combined_res$comparison, "last_mmse_test_score"), "log2FoldChange"] <- -1 * combined_res[startsWith(combined_res$comparison, "last_mmse_test_score"), "log2FoldChange"]

combined_wide <- pivot_wider(combined_res, id_cols="comparison", names_from="gene", values_from="log2FoldChange") %>% as.data.frame()
rownames(combined_wide) <- combined_wide$comparison
combined_wide$comparison <- NULL

combined_wide[is.na(combined_wide)] <- 0
scale_pca <- TRUE
comb_pca <- prcomp(combined_wide, scale = scale_pca)

comb_pcs <- as.data.frame(comb_pca$x)
comb_pcs$comparison <- rownames(comb_pcs)
comb_pcs$clinfeature <- str_split_fixed(comb_pcs$comparison, " * ", n=2)[, 1]
comb_pcs$layer <- str_split_fixed(comb_pcs$comparison, " * ", n=2)[, 2]
comb_pcs$layer <- gsub("^\\* ", "", comb_pcs$layer)
comb_pcs$layer <- gsub("_", " ", comb_pcs$layer)
# pdf(sprintf("all_clinical_effects_pca_scale%s_03292024.pdf", scale_pca))
# ggplot(comb_pcs, aes(x= PC1, y= PC2, label=comparison)) +
#   geom_point() + geom_text(hjust="inward", vjust="inward") + theme_bw()
# dev.off()

# test UMAP
comb_umap <- umap(comb_pcs[, colnames(comb_pcs)[startsWith(colnames(comb_pcs), "PC")]])
comb_umap <- comb_umap$layout
colnames(comb_umap) <- c("UMAP1", "UMAP2")
assert(all(rownames(comb_umap) == rownames(comb_pcs)))
comb_pcs <- cbind(comb_pcs, comb_umap)

group.colors <- c("Layer 1" = "#f43545", "Layer 2"="#ff8901", "Layer 3"="#fad717", "Layer 4"="#00ba71", "Layer 5"="#00c2de", "Layer 6"="#ba55d3", "WM"="#652175")
# clin.colors <- c("layer_vs" = 16, "nctx_temporal_fix" = 110, "brain_stem_sn" = 98, "motor_updrs_score" = 111, "last_mmse_test_score" = 109, "diagnosis_cont_ilbd" = 105, "diagnosis_cont_case" = 99)
clin.colors <- c("layer_vs" = 16, "nctx_temporal_fix" = 110, "motor_updrs_score" = 111, "last_mmse_test_score" = 109)
# pdf(sprintf("all_clinical_effects_pca_scale%s_04042024.pdf", scale_pca))
# pdf(sprintf("all_clinical_effects_pca_scale%s_layervs%s_04182024.pdf", scale_pca, include_layer_vs))
pdf(sprintf("all_clinical_effects_pca_scale%s_layervs%s_04262024.pdf", scale_pca, include_layer_vs))
ggplot(comb_pcs, aes(x= PC1, y= PC2, group=clinfeature, shape=clinfeature, color=layer)) +
    geom_point(size=3) +
    scale_shape_manual(values=clin.colors) +
    scale_color_manual(values=group.colors) +
    theme_bw()
ggplot(comb_pcs, aes(x= PC1, y= PC2, group=clinfeature, shape=clinfeature, color=layer)) +
    geom_point(size=5) +
    scale_shape_manual(values=clin.colors) +
    scale_color_manual(values=group.colors) +
    theme_bw() + xlim(-50,25) + ylim(-50,25) + coord_fixed(ratio=1.2)
ggplot(comb_pcs, aes(x= PC1, y= PC3, group=clinfeature)) +
    geom_point(aes(shape=clinfeature, color=layer)) +
    scale_color_manual(values=group.colors) +
    theme_bw()
ggplot(comb_pcs, aes(x= PC2, y= PC3, group=clinfeature)) +
    geom_point(aes(shape=clinfeature, color=layer)) +
    scale_color_manual(values=group.colors) +
    theme_bw()
ggplot(comb_pcs, aes(x= PC3, y= PC4, group=clinfeature)) +
    geom_point(aes(shape=clinfeature, color=layer)) +
    scale_color_manual(values=group.colors) +
    theme_bw()
ggplot(comb_pcs, aes(x= UMAP1, y= UMAP2, group=clinfeature)) +
    geom_point(aes(shape=clinfeature, color=layer)) +
    scale_color_manual(values=group.colors) +
    theme_bw()
dev.off()

# results$nlog10_pval <- -log10(results$pvalue)
# max_val <- max(results[!is.infinite(results$nlog10_pval), ]$nlog10_pval)
# inf_subst <- round(max_val + 50, -2)
# if ( any(is.infinite(results$nlog10_pval)) ) {
#     results[is.infinite(results$nlog10_pval), ]$nlog10_pval <- inf_subst
# }







# Venn diagram of DE genes across traits
library(stringr)
library(testit)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggVennDiagram)

source("/Users/jy1008/Downloads/pseudobulk_DE_DESeq2_GSEA_12292023_utils.R")

random_seed <- 101
layers <- c("Layer 1", "Layer 2", "Layer 3", "Layer 4", "Layer 5", "Layer 6", "WM")
effects_path <- "/Users/jy1008/Downloads/pseudobulk_DE_DESeq2_V4_12052023"

# out_label <- 01102024
# out_label <- "01162024"
out_gsea_path <- "/Users/jy1008/Downloads/pseudobulk_DE_GSEA_V4_01102024"

clinical_features <- c("diagnosis_cont_ilbd",
                       "diagnosis_cont_case",
                       "last_mmse_test_score",
                       "motor_updrs_score",
                       "nctx_temporal_fix",
                       "diagnosis_contilbd_case")

combined_res <- c()
for(clinical_feat in clinical_features) {
    de_files_list <- gsea_clinical_per_layer_filelist(effects_path, clinical_feat, layers)
    # clinical_feat is only used to fill the "comparison" field in output summary files, i.e. "trait * Layer 1"

    for(st_layer in names(de_files_list)) {
        print(sprintf("%s %s", clinical_feat, st_layer))
        
        results <- read.csv(de_files_list[[st_layer]])
        results$gene <- results$X
        results$X <- NULL
        results$comparison <- sprintf("%s * %s", clinical_feat, gsub(' ', '_', st_layer))
        combined_res <- rbind(combined_res, results)
    }
}
# apply a loose filter for the p-value
pval_filter <- 0.1
combined_res <- combined_res[combined_res$pvalue < pval_filter, ]

# flip the sign of last mmse test score (lower value is later PD)
combined_res[startsWith(combined_res$comparison, "last_mmse_test_score"), "log2FoldChange"] <- -1 * combined_res[startsWith(combined_res$comparison, "last_mmse_test_score"), "log2FoldChange"]


# prepare list for venn diagram
combined_res_sig <- combined_res[!is.na(combined_res$padj) & combined_res$padj < 0.05, ]
combined_res_sig$clin_trait <- str_split_fixed(combined_res_sig$comparison, " * ", n=2)[, 1]
venn_allgenes <- list(
    "nctx_temporal_fix" = unique(combined_res_sig[combined_res_sig$clin_trait == "nctx_temporal_fix", "gene"]),
    "last_mmse_test_score" = unique(combined_res_sig[combined_res_sig$clin_trait == "last_mmse_test_score", "gene"]),
    "motor_updrs_score" = unique(combined_res_sig[combined_res_sig$clin_trait == "motor_updrs_score", "gene"]),
    "diagnosis_contilbd_case" = unique(combined_res_sig[combined_res_sig$clin_trait == "diagnosis_contilbd_case", "gene"])
)
pdf("st_de_venn_diagram_06132024.pdf")
ggVennDiagram(venn_allgenes)
dev.off()

# get list of genes shared across all
common_de_set <- venn_allgenes[[1]]
for(delist_idx in 2:length(venn_allgenes)) {
    common_de_set <- intersect(common_de_set, venn_allgenes[[delist_idx]])
}
