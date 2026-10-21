library(Seurat)
library(DESeq2)
library(stringr)
library(testit)

library(ggplot2)
library(GGally)

library(UpSetR)

dds_lrt <- readRDS("pseudobulk_clinical_deseq2_dds_lrt_diagnosis_Layer6_12122023.rds")

tempQ <- results(dds_lrt, name="sample_meta...feat_ordered..Q")
tempL <- results(dds_lrt, name="sample_meta...feat_ordered..L")
# mcols(dds_lrt)[, "sample_meta...feat_ordered..L"]

tempQ <- tempQ[!is.na(tempQ$padj) & tempQ$padj < 0.05, ]
tempL <- tempL[!is.na(tempL$padj) & tempL$padj < 0.05, ]


tempL <- as.data.frame(tempL)
colnames(tempL) <- c("baseMean", "L.beta", "L.SE", "stat", "pvalue", "padj")
tempQ <- as.data.frame(tempQ)
colnames(tempQ) <- c("baseMean", "Q.beta", "Q.SE", "stat", "pvalue", "padj")

temp_comb <- merge(tempL, tempQ[, c("Q.beta", "Q.SE")], by=0)
rownames(temp_comb) <- temp_comb$Row.names
temp_comb$Row.names <- NULL
temp_comb <- temp_comb[, c("baseMean", "stat", "pvalue", "padj", "L.beta", "L.SE", "Q.beta", "Q.SE")]


ggplot(temp_comb, aes(x=L.beta, y=Q.beta)) +
  geom_point() + # Show dots
  geom_text(
    label=rownames(temp_comb),
    nudge_x = 0.02, nudge_y = 0.02,
    check_overlap = T
  ) + theme_bw()







# layer-wise union of all DE genes
layers <- sort(unique(comb_seurat@meta.data$smoothed_label_s5))
clinical_feats <- c("nctx_temporal", "brain_stem_sn", "motor_updrs_score", "last_mmse_test_score") # note diagnosis added separately
summary_stat_dir <- "/mnt/data0/projects/ASAP/ST/clustering_pipeline/pseudobulk_DE_DESeq2_V4_12052023"
res_all_all_layers <- c()
for(layer in rev(layers)) { # rev so that columns missing from Layer 1 are in correct order
    
    res_all <- c()
    for(clinical_feat in clinical_feats) {
        csv_file <- file.path(summary_stat_dir, clinical_feat, sprintf("pseudobulk_clinical_deseq2_%s_%s_12052023.csv", clinical_feat, gsub(' ', '', layer)))
        if(!file.exists(csv_file)) {
            print(paste("Skipped file", csv_file))
            next
        }
        res <- read.csv(csv_file, row.names=1)
        colnames(res) <- paste(colnames(res), clinical_feat, sep='_')
        res_all <- merge(res_all, res, by=0, all=TRUE)
        rownames(res_all) <- res_all$Row.names
        res_all$Row.names <- NULL
        print(paste("Finished", layer, clinical_feat))
    }
    # add diagnosis ordinal factor
    dds_file <- file.path(summary_stat_dir, "diagnosis_ordinal_12122023", sprintf("pseudobulk_clinical_deseq2_dds_lrt_%s_%s_12122023.rds", "diagnosis", gsub(' ', '', layer)))
    dds_lrt <- readRDS(dds_file)

    tempQ <- results(dds_lrt, name="sample_meta...feat_ordered..Q")
    tempL <- results(dds_lrt, name="sample_meta...feat_ordered..L")
    # mcols(dds_lrt)[, "sample_meta...feat_ordered..L"]
    # tempQ <- tempQ[!is.na(tempQ$padj) & tempQ$padj < 0.05, ]
    # tempL <- tempL[!is.na(tempL$padj) & tempL$padj < 0.05, ]

    tempL <- as.data.frame(tempL)
    colnames(tempL) <- c("baseMean", "L.beta", "L.SE", "stat", "pvalue", "padj")
    tempQ <- as.data.frame(tempQ)
    colnames(tempQ) <- c("baseMean", "Q.beta", "Q.SE", "stat", "pvalue", "padj")

    temp_comb <- merge(tempL, tempQ[, c("Q.beta", "Q.SE")], by=0)
    rownames(temp_comb) <- temp_comb$Row.names
    temp_comb$Row.names <- NULL
    # rearrange columns
    temp_comb <- temp_comb[, c("baseMean", "L.beta", "L.SE", "Q.beta", "Q.SE", "stat", "pvalue", "padj")]
    colnames(temp_comb) <- paste(colnames(temp_comb), "diagnosis", sep='_')

    pdf(file.path(summary_stat_dir, "diagnosis_ordinal_12122023", sprintf("pseudobulk_clinical_deseq2_diagnosis_ordinal_%s_12122023.pdf", gsub(' ', '', layer))))
    temp_comb_filt <- temp_comb[!is.na(temp_comb$padj) & temp_comb$padj < 0.05, ]
    sp_pt <- ggplot(temp_comb_filt, aes(x=L.beta_diagnosis, y=Q.beta_diagnosis)) + geom_point() + geom_text( label=rownames(temp_comb_filt), nudge_x=0.02, nudge_y=0.02, check_overlap=T ) + theme_bw()
    print(sp_pt)
    dev.off()
    
    res_all <- merge(res_all, temp_comb, by=0, all=TRUE)
    rownames(res_all) <- res_all$Row.names
    res_all$Row.names <- NULL
    print(paste("Finished", layer, "diagnosis"))
    
    res_all$layer <- layer
    res_all$gene <- rownames(res_all)
    res_all_all_layers <- plyr::rbind.fill(res_all_all_layers, res_all)
}
# reorder columns
write.csv(res_all_all_layers, file.path(summary_stat_dir, "pseudobulk_clinical_deseq2_all_summary_stats_12122023.csv"))

# upset plot of DE genes across 5 clinical features
venn_sets <- list()
for(clinical_feat in c(clinical_feats, "diagnosis")) {
    venn_sets[[clinical_feat]] = unique(res_all_all_layers[!is.na(res_all_all_layers[, paste0("padj_", clinical_feat)]) & res_all_all_layers[, paste0("padj_", clinical_feat)]<0.05, "gene"])
}
pdf(file.path(summary_stat_dir, "clinical_outs", "pseudobulk_clinical_deseq2_all_summary_stats_venn_upset_12122023.pdf"))
upset(fromList(venn_sets), order.by = "freq")
dev.off()


# plot individual genes
agg_counts <- list()
agg_normcounts <- list()
for(layer in layers) {
    agg_file <- file.path(summary_stat_dir, "diagnosis_ordinal_12122023", sprintf("agg_counts_filtered_smoothed_label_s5_%s_V4_12052023.csv", gsub(' ', '', layer)))
    agg_counts[[layer]] <- read.csv(agg_file, row.names=1)
    
    agg_normcounts[[layer]] <- sweep(agg_counts[[layer]],2,colSums(agg_counts[[layer]]),`/`) * 10000
}
meta_clinical <- comb_seurat@meta.data[, c("sample_name", "diagnosis", "brain_stem_sn", "nctx_temporal_fix", "motor_updrs_score", "last_mmse_test_score")]
meta_clinical <- meta_clinical[!duplicated(meta_clinical$sample_name), ]
meta_clinical$diagnosis <- factor(meta_clinical$diagnosis, ordered=TRUE, levels=c("Control", "ILBD", "Case"))

plot_gene_diagnosis <- function(gene, agg_list) {
    # agg_list is either agg_counts or agg_normcounts
    pdf(sprintf("pseudobulk_V4_%s_vs_diagnosis_12122023.pdf", gene))
    for(layer in layers) {
        if(gene %in% rownames(agg_list[[layer]])) {
            agg_gene <- as.data.frame(t(agg_list[[layer]][gene, ]))
            colnames(agg_gene) <- c("gene_name")
            agg_gene$sample_name <- rownames(agg_gene)
            temp <- merge(agg_gene, meta_clinical[, c("sample_name", "diagnosis")])
            # pdf(sprintf("pseudobulk_V4_%s_%s_vs_diagnosis_12122023.pdf", gsub(' ', '', layer), gene))
            sp_pt <- ggplot(temp, aes(x=diagnosis, y=gene_name)) + geom_boxplot() + theme_bw()
            print(sp_pt)
            
            # for a closer look at the learned quadratic function
            # plot the function along with the standardized cont/ilbd/case inputs
            
            temp$diagnosisnum <- as.numeric(temp$diagnosis)
            temp$diagnosisnum <- scale(temp$diagnosisnum, center=TRUE, scale=FALSE)
            # temp$diagnosisnum <- as.factor(as.character(temp$diagnosisnum))
            
            
            dds_file <- file.path(summary_stat_dir, "diagnosis_ordinal_12122023", sprintf("pseudobulk_clinical_deseq2_dds_lrt_%s_%s_12122023.rds", "diagnosis", gsub(' ', '', layer)))
            dds_lrt <- readRDS(dds_file)
            poly_a <- mcols(dds_lrt)[gene, "sample_meta...feat_ordered..Q"]
            poly_b <- mcols(dds_lrt)[gene, "sample_meta...feat_ordered..L"]
            # poly_c <- mcols(dds_lrt)[gene, "Intercept"]
            poly_c <- mean(temp$gene_name)
            x = seq(-1.5, 1.5, by = 0.05)
            temp2 <- data.frame(x = x,
                                y = poly_a * x**2 + poly_b * x + poly_c)
            sp_pt <- ggplot() + geom_boxplot(data=temp, aes(x=diagnosisnum, y=gene_name, group=diagnosisnum), width=0.3) +
                     geom_line(data=temp2, aes(x=x, y=y), color="red") + theme_bw()
            print(sp_pt)
            
            # dds_lrt@colData
            #
            # table(as.data.frame(dds_lrt@design[, c("sample_meta[, feat_ordered].L", "sample_meta[, feat_ordered].Q")]))
            
            
            # temp <- as.data.frame(dds_lrt@design[, c("sample_meta[, feat_ordered].L", "sample_meta[, feat_ordered].Q")])
            # temp <- cbind(temp, dds_lrt@colData[, c("sample_name", "diagnosis")])
            
            # recover the contrasts using
            # contrasts(temp$diagnosis)
            
            # temp$diagnosis.num <- as.numeric(temp$diagnosis)
            # temp$diagnosis.num <- scale(temp$diagnosis.num, center=TRUE, scale=FALSE)
            # sp_pt <- ggplot(temp, aes(x=diagnosis.num, y=gene_name))
            
            # dev.off()
        }
    }
    dev.off()
}
plot_gene_diagnosis("MTSS2", agg_normcounts)
plot_gene_diagnosis("FMNL2", agg_normcounts)
plot_gene_diagnosis("NRSN2", agg_normcounts)

plot_gene_diagnosis("NEAT1", agg_normcounts)
plot_gene_diagnosis("DNAJB1", agg_normcounts)
plot_gene_diagnosis("VGF", agg_normcounts)
plot_gene_diagnosis("SST", agg_normcounts)
plot_gene_diagnosis("CIB2", agg_normcounts)
plot_gene_diagnosis("ADCYAP1", agg_normcounts)





# Comparison of clinical features
comb_seurat <- readRDS("/mnt/data0/projects/ASAP/ST/clustering_pipeline/layer_smoothing_V4_seed101_11302023_replication_test/comb_seurat_post_smoothed_label_s5_11302023.rds")

comb_seurat@meta.data$nctx_temporal_fix <- comb_seurat@meta.data$nctx_temporal
comb_seurat@meta.data$nctx_temporal_fix[comb_seurat@meta.data$nctx_temporal_fix > 3] <- 3

meta_clinical <- comb_seurat@meta.data[, c("sample_name", "diagnosis", "brain_stem_sn", "nctx_temporal_fix", "motor_updrs_score", "last_mmse_test_score")]
meta_clinical <- meta_clinical[!duplicated(meta_clinical$sample_name), ]
meta_clinical$diagnosis <- factor(meta_clinical$diagnosis, ordered=TRUE, levels=c("Control", "ILBD", "Case"))
meta_clinical$nctx_temporal_fix_factor <- factor(meta_clinical$nctx_temporal_fix, ordered=TRUE)
meta_clinical$brain_stem_sn_factor <- factor(meta_clinical$brain_stem_sn, ordered=TRUE)

pdf("clinical_features_V4_12142023.pdf")
sp_pt <- ggpairs(meta_clinical[, c("diagnosis", "brain_stem_sn", "nctx_temporal_fix", "motor_updrs_score", "last_mmse_test_score")]) + theme_bw()
print(sp_pt)
sp_pt <- ggpairs(meta_clinical[, c("diagnosis", "brain_stem_sn_factor", "nctx_temporal_fix_factor", "motor_updrs_score", "last_mmse_test_score")]) + theme_bw()
print(sp_pt)
dev.off()




# 5/2/2024 Plot single gene vs nctx_temporal_fix, and cont/ilbd/case
library(Seurat)
library(testit)
library(dplyr)
comb_seurat <- readRDS("/mnt/data0/projects/ASAP/ST/clustering_pipeline/layer_smoothing_V4_seed101_11302023_replication_test/comb_seurat_post_smoothed_label_s5_11302023.rds")

comb_seurat@meta.data$nctx_temporal_fix <- comb_seurat@meta.data$nctx_temporal
comb_seurat@meta.data$nctx_temporal_fix[comb_seurat@meta.data$nctx_temporal_fix > 3] <- 3
assert(all(colnames(comb_seurat@assays$Spatial@data) == rownames(comb_seurat@meta.data)))

plot_gene_diagnosis <- function(gene_list, comb_seurat, outname, sd_bar=FALSE) {
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
    
    pdf(sprintf("pseudobulk_V4_diagnosis_vs_%s.pdf", outname))
    for(layer in layers) {
        meta_sub <- meta_data_agg[meta_data_agg$smoothed_label_s5 == layer, ]
        if(sd_bar) {
            meta_sub_agg2 <- meta_sub %>% group_by(diagnosis) %>%
                                summarise(diag_mean_sum = mean(lognorm_gene_sum),
                                          diag_sd_sum = sd(lognorm_gene_sum)) %>%
                                as.data.frame()
            sp_pt <- ggplot(meta_sub_agg2, aes(x=diagnosis, y=diag_mean_sum, group=1)) +
            geom_line() + geom_point() +
            geom_errorbar(aes(ymin=diag_mean_sum-diag_sd_sum,
                              ymax=diag_mean_sum+diag_sd_sum), width=.4, linewidth=1) +
            theme_bw() + ggtitle(layer)
            sp_pt <- sp_pt + geom_jitter(data=meta_sub, aes(x=diagnosis, y=lognorm_gene_sum), width=0.1, height=0, size=0.5)
        } else {
            sp_pt <- ggplot(meta_sub, aes(x=diagnosis, y=lognorm_gene_sum)) +
                        geom_boxplot(outlier.color=NA) +
                        geom_jitter() +
                        theme_bw() + ggtitle(layer)
        }
        print(sp_pt)
    }
    dev.off()
    pdf(sprintf("pseudobulk_V4_nctx_temporal_fix_vs_%s.pdf", outname))
    for(layer in layers) {
        meta_sub <- meta_data_agg[meta_data_agg$smoothed_label_s5 == layer, ]
        if(sd_bar) {
            meta_sub_agg2 <- meta_sub %>% group_by(nctx_temporal_fix) %>%
                                summarise(diag_mean_sum = mean(lognorm_gene_sum),
                                          diag_sd_sum = sd(lognorm_gene_sum)) %>%
                                as.data.frame()
            sp_pt <- ggplot(meta_sub_agg2, aes(x=nctx_temporal_fix, y=diag_mean_sum, group=1)) +
            geom_line() + geom_point() +
            geom_errorbar(aes(ymin=diag_mean_sum-diag_sd_sum,
                              ymax=diag_mean_sum+diag_sd_sum), width=.4, linewidth=1) +
            theme_bw() + ggtitle(layer)
            sp_pt <- sp_pt + geom_jitter(data=meta_sub, aes(x=nctx_temporal_fix, y=lognorm_gene_sum), width=0.1, height=0, size=0.5)
        } else {
            sp_pt <- ggplot(meta_sub, aes(x=nctx_temporal_fix, y=lognorm_gene_sum)) +
                        geom_boxplot(outlier.color=NA) +
                        geom_jitter() +
                        theme_bw() + ggtitle(layer)
        }
        print(sp_pt)
    }
    dev.off()
}


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


plot_gene_diagnosis(c("SNCA"), comb_seurat, "SNCA_05022024")
plot_gene_diagnosis(c("RPL27"), comb_seurat, "RPL27_05022024")
plot_gene_diagnosis(c("ITGAM"), comb_seurat, "ITGAM_neutrophil_05022024")
plot_gene_diagnosis(c("PRDX6"), comb_seurat, "PRDX6_neutrophil_05022024")


# grab sets of leading edge genes (NOTE: this is one panda server - move from home to proper dir)
gsea_file <- "~/de_genes_DESeq2_pseudobulk_nctx_temporal_fix_GO-BP_01162024.fGSEA.xls"
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
plot_gene_diagnosis(gene_list, comb_seurat, "GO-BP_neutrophil_mediated_immunity_05022024")
plot_gene_diagnosis(gene_list, comb_seurat, "GO-BP_neutrophil_mediated_immunity_line_05022024", sd_bar=TRUE)
plot_gene_diagnosis_all_layers(gene_list, comb_seurat, "GO-BP_neutrophil_mediated_immunity_all_layers_05022024")

pathway <- "SRP-dependent cotranslational protein targeting to membrane"
gene_list <- c()
for(le_set in gsea_table[gsea_table$pathway==pathway, "leadingEdge"]) {
    le_set <- strsplit(le_set,split=", ",fixed=TRUE)[[1]]
    gene_list <- union(gene_list, le_set)
    print(length(gene_list))
}
plot_gene_diagnosis(gene_list, comb_seurat, "GO-BP_SRP_cotranslational_protein_05022024")
plot_gene_diagnosis(gene_list, comb_seurat, "GO-BP_SRP_cotranslational_protein_line_05022024", sd_bar=TRUE)

pathway <- "modulation of chemical synaptic transmission"
gene_list <- c()
for(le_set in gsea_table[gsea_table$pathway==pathway, "leadingEdge"]) {
    le_set <- strsplit(le_set,split=", ",fixed=TRUE)[[1]]
    gene_list <- union(gene_list, le_set)
    print(length(gene_list))
}
plot_gene_diagnosis(gene_list, comb_seurat, "GO-BP_modulation_synaptic_transmission_05022024")
plot_gene_diagnosis(gene_list, comb_seurat, "GO-BP_modulation_synaptic_transmission_line_05022024", sd_bar=TRUE)


pathway <- "signal release from synapse"
gene_list <- c()
for(le_set in gsea_table[gsea_table$pathway==pathway, "leadingEdge"]) {
    le_set <- strsplit(le_set,split=", ",fixed=TRUE)[[1]]
    gene_list <- union(gene_list, le_set)
    print(length(gene_list))
}
plot_gene_diagnosis(gene_list, comb_seurat, "GO-BP_modulation_signal_release_synapse_05022024")
plot_gene_diagnosis(gene_list, comb_seurat, "GO-BP_modulation_signal_release_synapse_line_05022024", sd_bar=TRUE)
plot_gene_diagnosis_all_layers(gene_list, comb_seurat, "GO-BP_modulation_signal_release_synapse_all_layers_05022024")














# for last MMSE test score (quantitative feature)
library(ggpubr)
plot_gene_quant <- function(gene_list, comb_seurat, clinical_feat, outname) {
    meta_data <- comb_seurat@meta.data[, c("sample_name", "nCount_Spatial", "diagnosis", "smoothed_label_s5", clinical_feat)]
    meta_data$diagnosis <- factor(meta_data$diagnosis, levels=c("Control", "ILBD", "Case"))
    meta_data$clinical_feat <- meta_data[, clinical_feat]
    # meta_data$clinical_feat <- factor(meta_data$clinical_feat)
    
    if(length(gene_list) == 1) {
        meta_data <- cbind(meta_data, comb_seurat@assays$Spatial@counts[gene_list[1], ])
        colnames(meta_data)[ncol(meta_data)] <- "gene"
    } else {
        meta_data <- cbind(meta_data, t(as.matrix(comb_seurat@assays$Spatial@counts[gene_list, ])))
        meta_data$gene <- rowSums(meta_data[, gene_list])
        meta_data <- meta_data[, c("sample_name", "nCount_Spatial", "diagnosis", "smoothed_label_s5", "clinical_feat", "gene")]
    }
    layers <- sort(unique(meta_data$smoothed_label_s5))
    meta_data <- meta_data[!is.na(meta_data$clinical_feat), ]
    
    # get aggregated sum counts
    meta_data_agg <- meta_data %>%
                    group_by(sample_name, diagnosis, smoothed_label_s5, clinical_feat) %>%
                    summarise(gene_sum = sum(gene), count_sum = sum(nCount_Spatial)) %>%
                    as.data.frame()
    scale_factor <- 10^-4
    lognorm_const <- 1.0
    meta_data_agg$lognorm_gene_sum <- log(meta_data_agg$gene_sum / meta_data_agg$count_sum / scale_factor + lognorm_const)

    pdf(sprintf("pseudobulk_V4_%s_vs_%s.pdf", clinical_feat, outname))
    for(layer in layers) {
        meta_sub <- meta_data_agg[meta_data_agg$smoothed_label_s5 == layer, ]
        sp_pt <- ggplot(meta_sub, aes(x=clinical_feat, y=lognorm_gene_sum)) +
                    geom_point(aes(color=diagnosis)) +
                    geom_smooth(method = "lm",
                                formula = y ~ x) +
                    stat_cor(label.y=3.2) +
                    stat_regline_equation(label.y=3.1) +
                    theme_bw() + ggtitle(layer)
        print(sp_pt)
    }
    dev.off()
}

# grab sets of leading edge genes (NOTE: this is one panda server - move from home to proper dir)
gsea_file <- "~/de_genes_DESeq2_pseudobulk_last_mmse_test_score_GO-BP_01162024.fGSEA.xls"
gsea_table <- read.csv(gsea_file, sep='\t')

deseq2_file <- "~/pseudobulk_clinical_deseq2_last_mmse_test_score_WM_12052023.csv"
deseq2_table <- read.csv(deseq2_file, sep=',')
colnames(deseq2_table)[1] <- "gene"
deseq2_table_sig <- deseq2_table[!is.na(deseq2_table$padj) & deseq2_table$padj < 0.05, ]

pathway <- "synaptic vesicle cycle"
gene_list <- c()
# aggregate leading edges across all layers
for(le_set in gsea_table[gsea_table$pathway==pathway, "leadingEdge"]) {
# use only leading edges from WM
# for(le_set in gsea_table[gsea_table$pathway==pathway & grepl("WM", gsea_table$comparison), "leadingEdge"]) {
    le_set <- strsplit(le_set,split=", ",fixed=TRUE)[[1]]
    gene_list <- union(gene_list, le_set)
    print(length(gene_list))
}
gene_list_sig <- intersect(gene_list, deseq2_table_sig$gene)
deseq2_table_sig[deseq2_table_sig$gene %in% gene_list_sig, ]
plot_gene_quant(gene_list, comb_seurat, "last_mmse_test_score", "GO-BP_synaptic_vesicle_cycle_05022024")
plot_gene_quant(gene_list_sig, comb_seurat, "last_mmse_test_score", "GO-BP_synaptic_vesicle_cycle_WMsig_05022024")






# ------
# 6/12/2024:

meta_clinical <- comb_seurat@meta.data[, c("sample_name", "diagnosis", "nctx_temporal", "nctx_temporal_fix", "motor_updrs_score", "last_mmse_test_score")]
meta_clinical <- meta_clinical[!duplicated(meta_clinical$sample_name), ]
write.csv(meta_clinical, "meta_clinical_neutrophil_debug_06122024.csv")

# grab sets of leading edge genes (NOTE: this is one panda server - move from home to proper dir)
gsea_file <- "~/de_genes_DESeq2_pseudobulk_nctx_temporal_fix_GO-BP_01162024.fGSEA.xls"
gsea_table <- read.csv(gsea_file, sep='\t')
head(gsea_table)
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


meta_data <- comb_seurat@meta.data[, c("sample_name", "nCount_Spatial", "diagnosis", "smoothed_label_s5", "nctx_temporal", "last_mmse_test_score", "motor_updrs_score")]
# meta_data <- cbind(meta_data, t(as.matrix(comb_seurat@assays$Spatial@counts[gene_list, ])))
meta_data <- cbind(meta_data, t(as.matrix(comb_seurat@assays$Spatial@data[gene_list, ])))
# meta_data$gene <- rowSums(meta_data[, gene_list])
# meta_data <- meta_data[, c("sample_name", "nCount_Spatial", "diagnosis", "smoothed_label_s5", "nctx_temporal", "nctx_temporal_fix", "last_mmse_test_score", "motor_updrs_score", "gene")]

layers <- sort(unique(meta_data$smoothed_label_s5))
# meta_data <- meta_data[!is.na(meta_data$nctx_temporal_fix), ]

# get aggregated sum counts
# meta_data_agg <- meta_data %>%
#                 group_by(sample_name, diagnosis, smoothed_label_s5, nctx_temporal_fix) %>%
#                 summarise(gene_sum = sum(gene), count_sum = sum(nCount_Spatial)) %>%
#                 as.data.frame()
meta_data_agg <- meta_data %>%
                group_by(sample_name, diagnosis, smoothed_label_s5, nctx_temporal, last_mmse_test_score, motor_updrs_score) %>%
                summarise_at(gene_list, mean, na.rm=TRUE)  %>%
                as.data.frame()
# scale_factor <- 10^-4
# lognorm_const <- 1.0
# meta_data_agg$lognorm_gene_sum <- log(meta_data_agg$gene_sum / meta_data_agg$count_sum / scale_factor + lognorm_const)
write.csv(meta_data_agg, "meta_clinical_neutrophil_debug2_06122024.csv")




plot_gene_neutrophil_investigation <- function(gene_list, outname) {
    # BN2003, BN1504, BN1412, BN1865, BN1812
    # NOTE: some of these are our worst performant samples, though not all case samples are like this
    # redo the boxplots, but group by control, ilbd, case, late-case
    meta_data <- comb_seurat@meta.data[, c("sample_name", "nCount_Spatial", "diagnosis", "smoothed_label_s5", "nctx_temporal")]
    meta_data$late_pd <- meta_data$nctx_temporal
    meta_data[meta_data$sample_name %in% c("BN1412", "BN1812", "BN1865", "BN2003"), "late_pd"] <- "late_case"

    meta_data$diagnosis <- factor(meta_data$diagnosis, levels=c("Control", "ILBD", "Case"))
    meta_data$nctx_temporal <- factor(meta_data$nctx_temporal)
    meta_data$late_pd <- factor(meta_data$late_pd, levels=c('0', '1', '2', '3', '4', 'late_case'))

    if(length(gene_list) == 1) {
        meta_data <- cbind(meta_data, comb_seurat@assays$Spatial@counts[gene_list[1], ])
        colnames(meta_data)[ncol(meta_data)] <- "gene"
    } else {
        meta_data <- cbind(meta_data, t(as.matrix(comb_seurat@assays$Spatial@counts[gene_list, ])))
        meta_data$gene <- rowSums(meta_data[, gene_list])
        meta_data <- meta_data[, c("sample_name", "nCount_Spatial", "diagnosis", "smoothed_label_s5", "nctx_temporal", "late_pd", "gene")]
    }
    layers <- sort(unique(meta_data$smoothed_label_s5))
    # meta_data <- meta_data[!is.na(meta_data$nctx_temporal_fix), ]
    meta_data <- meta_data[!is.na(meta_data$late_pd), ]

    # get aggregated sum counts
    meta_data_agg <- meta_data %>%
                    group_by(sample_name, diagnosis, smoothed_label_s5, late_pd, nctx_temporal) %>%
                    summarise(gene_sum = sum(gene), count_sum = sum(nCount_Spatial)) %>%
                    as.data.frame()
    scale_factor <- 10^-4
    lognorm_const <- 1.0
    meta_data_agg$lognorm_gene_sum <- log(meta_data_agg$gene_sum / meta_data_agg$count_sum / scale_factor + lognorm_const)

    group.colors <- c("Layer 1" = "#f43545", "Layer 2"="#ff8901", "Layer 3"="#fad717", "Layer 4"="#00ba71", "Layer 5"="#00c2de", "Layer 6"="#ba55d3", "WM"="#652175")
    pdf(sprintf("pseudobulk_V4_late_case_investigation_%s_06142024.pdf", outname))
    meta_sub_agg <- meta_data_agg %>% group_by(smoothed_label_s5, late_pd) %>%
                        summarise(diag_mean_sum = mean(lognorm_gene_sum),
                                  diag_sd_sum = sd(lognorm_gene_sum)) %>%
                        as.data.frame()
    sp_pt <- ggplot(meta_sub_agg, aes(x=late_pd, y=diag_mean_sum, color=smoothed_label_s5, group=smoothed_label_s5)) +
    geom_line(position=position_dodge(0.3)) +
    geom_point(position=position_dodge(0.3)) +
    geom_errorbar(aes(ymin=diag_mean_sum-diag_sd_sum,
                      ymax=diag_mean_sum+diag_sd_sum),
                  position=position_dodge(0.3),
                  width=.5, linewidth=1) +
    scale_color_manual(values=group.colors) + theme_bw()
    print(sp_pt)

    meta_sub_agg <- meta_data_agg %>% group_by(smoothed_label_s5, nctx_temporal) %>%
                        summarise(diag_mean_sum = mean(lognorm_gene_sum),
                                  diag_sd_sum = sd(lognorm_gene_sum)) %>%
                        as.data.frame()
    sp_pt <- ggplot(meta_sub_agg, aes(x=nctx_temporal, y=diag_mean_sum, color=smoothed_label_s5, group=smoothed_label_s5)) +
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
plot_gene_neutrophil_investigation(c("VEGFA"), "VEGFA")
plot_gene_neutrophil_investigation(c("VEGFB"), "VEGFB")
plot_gene_neutrophil_investigation(c("VEGFC"), "VEGFC")
plot_gene_neutrophil_investigation(c("VEGFD"), "VEGFD")
plot_gene_neutrophil_investigation(c("ICAM1"), "ICAM1")
plot_gene_neutrophil_investigation(c("HBB"), "HBB")

plot_gene_neutrophil_investigation(c("CD177"), "CD177")
plot_gene_neutrophil_investigation(c("CCL2"), "CCL2")
plot_gene_neutrophil_investigation(c("ARG2"), "ARG2")

# from neutrophil list
summary_stat_dir <- "/mnt/data0/projects/ASAP/ST/clustering_pipeline/pseudobulk_DE_DESeq2_V4_12052023"
res_all_all_layers <- read.csv(file.path(summary_stat_dir, "pseudobulk_clinical_deseq2_all_summary_stats_12122023.csv"))
res_all_all_layers$X <- NULL
res_all_all_layers <- res_all_all_layers[!is.na(res_all_all_layers$padj_nctx_temporal) & res_all_all_layers$padj_nctx_temporal < 0.05, ]
# NOTE this is on nctx_temporal with single sample with max score = 4, not nctx_temporal_fix with max score 3
# > sort(intersect(gene_list, unique(res_all_all_layers$gene)))
#  [1] "A1BG"     "ACAA1"    "ACLY"     "ACTR10"   "ACTR1B"   "ADAM10"
#  [7] "ADGRE5"   "AGA"      "AP1M1"    "APEH"     "APRT"     "ASAH1"
# [13] "ATP6V1D"  "COMMD3"   "CPPED1"   "CSNK2B"   "DDOST"    "DEGS1"
# [19] "DERA"     "DNASE1L1" "DPP7"     "DYNLT1"   "FUCA1"    "FUCA2"
# [25] "GAA"      "GSN"      "HBB"      "HMOX2"    "ILF2"     "IMPDH2"
# [31] "IST1"     "LAMTOR1"  "LAMTOR2"  "MOSPD2"   "MVP"      "NEU1"
# [37] "NIT2"     "OSTF1"    "PADI2"    "PDAP1"    "PFKL"     "PKM"
# [43] "PPIE"     "PSAP"     "PSMA5"    "PSMB1"    "PSMC3"    "PSMD13"
# [49] "PSMD6"    "PTGES2"   "PTPRN2"   "PYGB"     "RAB4B"    "RAB5C"
# [55] "SERPINB1" "TMEM179B" "TOM1"     "TRAPPC1"  "VPS35L"   "XRCC6"

plot_gene_neutrophil_investigation(c("SERPINB1"), "SERPINB1")
plot_gene_neutrophil_investigation(c("PSMA5"), "PSMA5")
plot_gene_neutrophil_investigation(c("PKM"), "PKM")

group.colors <- c("Layer 1" = "#f43545", "Layer 2"="#ff8901", "Layer 3"="#fad717", "Layer 4"="#00ba71", "Layer 5"="#00c2de", "Layer 6"="#ba55d3", "WM"="#652175")
pdf(sprintf("pseudobulk_V4_nctx_temporal_fix_vs_%s.pdf", outname))
meta_sub_agg <- meta_data_agg %>% group_by(smoothed_label_s5, nctx_temporal_fix) %>%
                    summarise(diag_mean_sum = mean(lognorm_gene_sum),
                              diag_sd_sum = sd(lognorm_gene_sum)) %>%
                    as.data.frame()
sp_pt <- ggplot(meta_sub_agg, aes(x=nctx_temporal_fix, y=diag_mean_sum, color=smoothed_label_s5, group=smoothed_label_s5))
