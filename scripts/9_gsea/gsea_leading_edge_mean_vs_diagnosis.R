# see ~/Downloads/pseudobulk_DE_DESeq2_V4_ordinal_plots_12132023.R


library(Seurat)
library(testit)
library(dplyr)
library(ggplot2)
library(ggrepel)

comb_seurat <- readRDS("/mnt/data/projects/ASAP/ST/clustering_pipeline/layer_smoothing_V4_seed101_11302023_replication_test/comb_seurat_post_smoothed_label_s5_11302023.rds")

comb_seurat@meta.data$nctx_temporal_fix <- comb_seurat@meta.data$nctx_temporal
comb_seurat@meta.data$nctx_temporal_fix[comb_seurat@meta.data$nctx_temporal_fix > 3] <- 3
assert(all(colnames(comb_seurat@assays$Spatial@data) == rownames(comb_seurat@meta.data)))



plot_gene_diagnosis_all_layers <- function(gene_list, comb_seurat, outname) {
    meta_data <- comb_seurat@meta.data[, c("sample_name", "nCount_Spatial", "diagnosis", "smoothed_label_s5", "nctx_temporal_fix")]
    meta_data$diagnosis <- factor(meta_data$diagnosis, levels=c("Control", "ILBD", "Case"))
    meta_data$nctx_temporal_fix <- factor(meta_data$nctx_temporal_fix)

    if(length(gene_list) == 1) {
        meta_data <- cbind(meta_data, comb_seurat@assays$Spatial@counts[gene_list[1], ])
        colnames(meta_data)[ncol(meta_data)] <- "gene"
    } else {
        meta_data <- cbind(meta_data, t(as.matrix(comb_seurat@assays$Spatial@counts[gene_list, ])))
        meta_data$gene <- rowSums(meta_data[, gene_list])
        meta_data <- meta_data[, c("sample_name", "nCount_Spatial", "diagnosis", "smoothed_label_s5", "nctx_temporal_fix", "gene")]
    }
    layers <- sort(unique(meta_data$smoothed_label_s5))
    meta_data <- meta_data[!is.na(meta_data$nctx_temporal_fix), ]

    # get aggregated sum counts
    meta_data_agg <- meta_data %>%
                    group_by(sample_name, diagnosis, smoothed_label_s5, nctx_temporal_fix) %>%
                    summarise(gene_sum = sum(gene), count_sum = sum(nCount_Spatial)) %>%
                    as.data.frame()
    scale_factor <- 10^-4
    lognorm_const <- 1.0
    meta_data_agg$lognorm_gene_sum <- log(meta_data_agg$gene_sum / meta_data_agg$count_sum / scale_factor + lognorm_const)

    group.colors <- c("Layer 1" = "#f43545", "Layer 2"="#ff8901", "Layer 3"="#fad717", "Layer 4"="#00ba71", "Layer 5"="#00c2de", "Layer 6"="#ba55d3", "WM"="#652175")
    pdf(sprintf("pseudobulk_V4_nctx_temporal_fix_vs_%s.pdf", outname))
    meta_sub_agg <- meta_data_agg %>% group_by(smoothed_label_s5, nctx_temporal_fix) %>%
                        summarise(diag_mean_sum = mean(lognorm_gene_sum),
                                  diag_sd_sum = sd(lognorm_gene_sum)) %>%
                        as.data.frame()
    sp_pt <- ggplot(meta_sub_agg, aes(x=nctx_temporal_fix, y=diag_mean_sum, color=smoothed_label_s5, group=smoothed_label_s5)) +
    geom_line(position=position_dodge(0.3)) +
    geom_point(position=position_dodge(0.3)) +
    geom_errorbar(aes(ymin=diag_mean_sum-diag_sd_sum,
                      ymax=diag_mean_sum+diag_sd_sum),
                  position=position_dodge(0.3),
                  width=.5, linewidth=1) +
    scale_color_manual(values=group.colors) + theme_bw()
    print(sp_pt)

    dev.off()
}

# gsea_file <- "/Users/jy1008/Downloads/pseudobulk_DE_GSEA_V4_01102024/GO-BP/de_genes_DESeq2_pseudobulk_nctx_temporal_fix_GO-BP_01162024.fGSEA.xls"
# on tiger server
gsea_file <- "/mnt/data/projects/ASAP/ST/GSEA_analysis_backup_09242024/pseudobulk_DE_GSEA_V4_01102024/GO-BP/de_genes_DESeq2_pseudobulk_nctx_temporal_fix_GO-BP_01162024.fGSEA.xls"
gsea_table <- read.csv(gsea_file, sep='\t')


pathway <- "neutrophil mediated immunity"
gene_list <- c()
# aggregate leading edges across all layers
for(le_set in gsea_table[gsea_table$pathway==pathway, "leadingEdge"]) {
# use only leading edges from WM
# for(le_set in gsea_table[gsea_table$pathway==pathway & grepl("WM", gsea_table$comparison), "leadingEdge"]) {
    le_set <- strsplit(le_set,split=", ",fixed=TRUE)[[1]]
    gene_list <- union(gene_list, le_set)
    print(length(gene_list))
}
# plot_gene_diagnosis(gene_list, comb_seurat, "GO-BP_neutrophil_mediated_immunity_05022024")
# plot_gene_diagnosis(gene_list, comb_seurat, "GO-BP_neutrophil_mediated_immunity_line_05022024", sd_bar=TRUE)
plot_gene_diagnosis_all_layers(gene_list, comb_seurat, "GO-BP_neutrophil_mediated_immunity_all_layers_09242024")


pathway <- "SRP-dependent cotranslational protein targeting to membrane"
gene_list <- c()
# aggregate leading edges across all layers
for(le_set in gsea_table[gsea_table$pathway==pathway, "leadingEdge"]) {
# use only leading edges from WM
# for(le_set in gsea_table[gsea_table$pathway==pathway & grepl("WM", gsea_table$comparison), "leadingEdge"]) {
    le_set <- strsplit(le_set,split=", ",fixed=TRUE)[[1]]
    gene_list <- union(gene_list, le_set)
    print(length(gene_list))
}
# plot_gene_diagnosis(gene_list, comb_seurat, "GO-BP_neutrophil_mediated_immunity_05022024")
# plot_gene_diagnosis(gene_list, comb_seurat, "GO-BP_neutrophil_mediated_immunity_line_05022024", sd_bar=TRUE)
plot_gene_diagnosis_all_layers(gene_list, comb_seurat, "GO-BP_SRP_cotranslational_protein_all_layers_09242024")


pathway <- "signal release from synapse"
gene_list <- c()
# aggregate leading edges across all layers
for(le_set in gsea_table[gsea_table$pathway==pathway, "leadingEdge"]) {
# use only leading edges from WM
# for(le_set in gsea_table[gsea_table$pathway==pathway & grepl("WM", gsea_table$comparison), "leadingEdge"]) {
    le_set <- strsplit(le_set,split=", ",fixed=TRUE)[[1]]
    gene_list <- union(gene_list, le_set)
    print(length(gene_list))
}
# plot_gene_diagnosis(gene_list, comb_seurat, "GO-BP_neutrophil_mediated_immunity_05022024")
# plot_gene_diagnosis(gene_list, comb_seurat, "GO-BP_neutrophil_mediated_immunity_line_05022024", sd_bar=TRUE)
plot_gene_diagnosis_all_layers(gene_list, comb_seurat, "GO-BP_signal_release_from_synapse_all_layers_09242024")




gene_list <- c("HBB", "TRAM1L1", "HCRT", "PSMB2", "POU3F3", "NR2F1", "NR3C1", "ARHGEF12", "KIF5A", "ARHGAP5", "JDP2", "APRT", "TSHZ1", "TMSB4X", "FAM199X", "CDC42", "NOTCH2", "RNF181", "SSH1", "STAT5B", "SERPINA3")
for(de_gene in gene_list) {
    print(de_gene)
    plot_gene_diagnosis_all_layers(c(de_gene), comb_seurat, sprintf("%s_all_layers_09242024", de_gene))
}





# read DE files from nctx_temporal_fix: check which ones are positively and negatively affected
de_genes <- c()
st_layers <- c(paste0("Layer", 1:6), "WM")
de_files_path <- "/mnt/data/projects/ASAP/ST/clustering_pipeline/pseudobulk_DE_DESeq2_V4_12052023"
clinical_feat <- "nctx_temporal_fix"
file_run_date <- "12052023"
for(st_layer in st_layers) {
    de_file <- sprintf("pseudobulk_clinical_deseq2_%s_%s_%s.csv", clinical_feat, st_layer, file_run_date)
    de_file <- file.path(de_files_path, clinical_feat, de_file)
    de_layer_list <- read.csv(de_file)
    de_layer_list$st_layer <- st_layer
    de_genes <- rbind(de_genes, de_layer_list)
}
colnames(de_genes)[colnames(de_genes) == "X"] <- "gene"
de_genes$nlog10_pvalue <- -log10(de_genes$pvalue)
de_genes[is.na(de_genes$padj), "padj"] <- 1

pathway <- "neutrophil mediated immunity"
gene_list <- c()
# aggregate leading edges across all layers
for(le_set in gsea_table[gsea_table$pathway==pathway, "leadingEdge"]) {
# use only leading edges from WM
# for(le_set in gsea_table[gsea_table$pathway==pathway & grepl("WM", gsea_table$comparison), "leadingEdge"]) {
    le_set <- strsplit(le_set,split=", ",fixed=TRUE)[[1]]
    gene_list <- union(gene_list, le_set)
    print(length(gene_list))
}
pdf(sprintf("de_gene_validation_%s_11122024.pdf", gsub(' ', '_', pathway)))
for(st_layer in st_layers) {
    # get the p-value corresponding to padj < 0.05
    de_genes_sub <- de_genes[de_genes$st_layer == st_layer, ]
    de_genes_sub_temp <- de_genes_sub[de_genes_sub$padj < 0.05, ]
    nom_pval_sig_thresh <- min(de_genes_sub_temp$nlog10_pvalue)

    # filter for genes in the pathway of interest
    de_genes_sub <- de_genes[de_genes$st_layer == st_layer & de_genes$gene %in% gene_list, ]
    # volcano plot
    pt <- ggplot(de_genes_sub, aes(x=log2FoldChange, y=nlog10_pvalue, label=gene)) +
          geom_point() +
          geom_text(hjust=0, vjust=0) +
          geom_hline(yintercept = nom_pval_sig_thresh, color = "red") +
          ggtitle(st_layer) +
          theme_bw()
    print(pt)
}
dev.off()

# compare to reference DE for PD

# see ~/Downloads/compare_de_effects_to_refs.R
# https://stackoverflow.com/a/7549819
lm_eqn <- function(df, x, y){
    m <- lm(as.formula(paste(y, '~', x)), df);
    eq <- sprintf("y == %s + %s * x, r^2 = %s \n SP corr = %s",
              format(unname(coef(m)[1]), digits = 2),
              format(unname(coef(m)[2]), digits = 2),
              format(summary(m)$r.squared, digits = 3),
              format(cor(df[, x], df[, y], method="spearman"), digits = 3))
    # as.character(as.expression(eq));
    return(eq)
}
count_point_quadrants <- function(x, y) {
    temp_xy <- data.frame(x=x, y=y)
    q1 <- nrow(temp_xy[temp_xy$x > 0 & temp_xy$y > 0, ])
    q2 <- nrow(temp_xy[temp_xy$x < 0 & temp_xy$y > 0, ])
    q3 <- nrow(temp_xy[temp_xy$x < 0 & temp_xy$y < 0, ])
    q4 <- nrow(temp_xy[temp_xy$x > 0 & temp_xy$y < 0, ])
    out_string <- sprintf("Q1:%s, Q2:%s, Q3:%s, Q4:%s", q1, q2, q3, q4)
    return(out_string)
}

de_vs_ref <- function(de_genes, ref_df, ref_df_name, pathway, gene_list, st_layers) {
    pdf(sprintf("de_gene_validation_vs_ref_%s_%s_11122024.pdf", gsub(' ', '_', pathway), ref_df_name))
    for(st_layer in st_layers) {
        # filter for genes in the pathway of interest
        de_genes_sub <- de_genes[de_genes$st_layer == st_layer & de_genes$gene %in% gene_list, ]
        temp_df <- merge(de_genes_sub[, c("gene", "log2FoldChange", "baseMean", "padj", "nlog10_pvalue")],
                         ref_df[, c("gene", "log2FoldChange")],
                         by="gene", suffixes=c(".ours", ".ref"))

        sp_pt_label <- paste(
            lm_eqn(temp_df, "log2FoldChange.ours", "log2FoldChange.ref"), '\n',
            count_point_quadrants(temp_df$log2FoldChange.ours, temp_df$log2FoldChange.ref))
        
        sp_pt <- ggplot(temp_df, aes(x=log2FoldChange.ours, y=log2FoldChange.ref, label=gene, size=baseMean, colour=nlog10_pvalue)) +
        geom_point() +
        scale_colour_viridis_c() +
        geom_smooth(method = "lm", formula = y ~ x) +
        geom_text_repel(hjust=0, vjust=0) +
        geom_hline(yintercept=0) +
        geom_vline(xintercept=0) +
        ggtitle(paste(st_layer, clinical_feat)) + theme_bw()
        
        sp_pt <- sp_pt + annotate(geom='text', label=sp_pt_label, x=-Inf, y=Inf, hjust=0, vjust=1)
        print(sp_pt)
    }
    dev.off()
}


ref_dir <- "/home/jy743/neutrophil_investigation_11122024/de_references"

ref_file <- file.path(ref_dir, "40478_2020_932_MOESM4_ESM.csv")
ref_df_name <- "40478_Nido_2020"
ref_df <- read.csv(ref_file)
colnames(ref_df)[colnames(ref_df) == "GeneSymbol"] <- "gene"
de_vs_ref(de_genes, ref_df, ref_df_name, pathway, gene_list, st_layers)

ref_file <- file.path(ref_dir, "41467_2024_47867_MOESM10_ESM.csv")
ref_df_name <- "41467_Flones_2024_CIPDm1"
ref_df <- read.csv(ref_file)
ref_df <- ref_df[, c("gene_name", "log2FoldChange_CI.PDmodel1", "pvalue_CI.PDmodel1", "padj_CI.PDmodel1")]
colnames(ref_df) <- c("gene", "log2FoldChange", "pvalue", "padj")
de_vs_ref(de_genes, ref_df, ref_df_name, pathway, gene_list, st_layers)

ref_file <- file.path(ref_dir, "41467_2024_47867_MOESM10_ESM.csv")
ref_df_name <- "41467_Flones_2024_CIPDm2"
ref_df <- read.csv(ref_file)
ref_df <- ref_df[, c("gene_name", "log2FoldChange_CI.PDmodel2", "pvalue_CI.PDmodel2", "padj_CI.PDmodel2")]
colnames(ref_df) <- c("gene", "log2FoldChange", "pvalue", "padj")
de_vs_ref(de_genes, ref_df, ref_df_name, pathway, gene_list, st_layers)

ref_file <- file.path(ref_dir, "Ruifeng_TableS4_neutrophil_genes.csv")
ref_df_name <- "Ruifeng_AMPPD_neutrophil"
ref_df <- read.csv(ref_file)
ref_df <- ref_df[, c("gene_name", "log2FoldChange_Rep", "pvalue_Rep", "padj_Rep")]
colnames(ref_df) <- c("gene", "log2FoldChange", "pvalue", "padj")
de_vs_ref(de_genes, ref_df, ref_df_name, pathway, ref_df$gene, st_layers)

ref_file <- file.path(ref_dir, "abo1997_Oligo_DEGs.csv")
ref_df_name <- "Zhu_2024_PD_SC_Science_TranslatMed"
ref_df <- read.csv(ref_file)
colnames(ref_df)[colnames(ref_df) == "X"] <- "gene"
ref_df <- ref_df[, c("gene", "avg_logFC", "p_val_adj")] # nominal p-value missing
colnames(ref_df) <- c("gene", "log2FoldChange", "padj")
de_vs_ref(de_genes, ref_df, ref_df_name, pathway, gene_list, st_layers)

ref_file <- file.path(ref_dir, "Keo_42003_2020.csv")
ref_df_name <- "Keo_2020_NatCommsBio_BrainRegVulPD"
ref_df <- read.csv(ref_file)
ref_df <- ref_df[, c("Gene.symbol", "Fold.change.between.R1.and.R6", "P.value.of.fold.change.between.R1.and.R6..BH.corrected.")] # nominal p-value missing
colnames(ref_df) <- c("gene", "log2FoldChange", "padj")
de_vs_ref(de_genes, ref_df, ref_df_name, pathway, gene_list, st_layers)





# It's unclear why synapse and neurotransmitter genes are enriched, especially in WM
# Test the possibility that these genes stay the same, while other genes (e.g. mitochondria?) are depleted
# pathway and gene_list defined above
temp_mat <- comb_seurat@assays$Spatial@counts[gene_list, ]
comb_seurat@meta.data$pathway_read_count <- colSums(temp_mat)
meta_data <- comb_seurat@meta.data[, c("sample_name", "smoothed_label_s5", "nctx_temporal_fix", "nCount_Spatial", "pathway_read_count")]
meta_agg <- meta_data %>%
                group_by(sample_name, smoothed_label_s5,
                         nctx_temporal_fix) %>%
                summarise(sum_total_count = sum(nCount_Spatial),
                          sum_pathway_count = sum(pathway_read_count)) %>%
                as.data.frame()
pdf("mod_chem_synapse_sum_idea_11132024.pdf")
for(st_layer in unique(meta_agg$smoothed_label_s5)) {
    meta_agg_sub <- meta_agg[meta_agg$smoothed_label_s5 == st_layer, ]
    pt <- ggplot(meta_agg_sub, aes(x=sum_pathway_count, y=sum_total_count, color=nctx_temporal_fix)) +
        geom_point() +
        ggtitle(st_layer) +
        theme_bw()
    print(pt)
    
    meta_agg_sub$pathway_prop <- meta_agg_sub$sum_pathway_count / meta_agg_sub$sum_total_count
    meta_agg_sub <- meta_agg_sub[order(meta_agg_sub$nctx_temporal_fix),]
    meta_agg_sub$sample_name <- factor(meta_agg_sub$sample_name, levels=meta_agg_sub$sample_name)
    pt <- ggplot(meta_agg_sub, aes(x=sample_name, y=pathway_prop, color=nctx_temporal_fix)) +
        geom_point() +
        ggtitle(st_layer) +
        theme_bw()
    print(pt)
}
dev.off()
