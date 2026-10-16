# see ~/Downloads/pseudobulk_DE_DESeq2_V4_ordinal_plots_12132023.R

library(Seurat)
library(DESeq2)
library(stringr)
library(testit)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(ggpmisc) # for regression trendline

# library(GGally)
# library(UpSetR)



# Comparison of clinical features
comb_seurat <- readRDS("/mnt/data/projects/ASAP/ST/clustering_pipeline/layer_smoothing_V4_seed101_11302023_replication_test/comb_seurat_post_smoothed_label_s5_11302023.rds")

comb_seurat@meta.data$nctx_temporal_fix <- comb_seurat@meta.data$nctx_temporal
comb_seurat@meta.data$nctx_temporal_fix[comb_seurat@meta.data$nctx_temporal_fix > 3] <- 3

# meta_clinical <- comb_seurat@meta.data[, c("sample_name", "diagnosis", "brain_stem_sn", "nctx_temporal_fix", "motor_updrs_score", "last_mmse_test_score")]
# meta_clinical <- meta_clinical[!duplicated(meta_clinical$sample_name), ]
meta_clinical <- comb_seurat@meta.data %>%
    group_by(sample_name, sex, expired_age_num,
             PMI_num, RIN_num, diagnosis,
             last_mmse_test_score, motor_updrs_score,
             nctx_temporal, brain_stem_sn) %>%
    summarise(nspots=n(),
              mean_read_count=mean(nCount_Spatial)) %>%
    as.data.frame()
meta_clinical$diagnosis <- as.character(meta_clinical$diagnosis)
meta_clinical$diagnosis[meta_clinical$diagnosis == "ILBD"] <- "Prodromal"
meta_clinical$diagnosis <- factor(meta_clinical$diagnosis, ordered=TRUE, levels=c("Control", "Prodromal", "Case"))
meta_clinical$nctx_temporal_factor <- factor(meta_clinical$nctx_temporal, ordered=TRUE)
meta_clinical$brain_stem_sn_factor <- factor(meta_clinical$brain_stem_sn, ordered=TRUE)

meta_clinical_sort <- meta_clinical %>% arrange(diagnosis, nctx_temporal_factor, brain_stem_sn_factor, desc(last_mmse_test_score))
meta_clinical_sort$sample_name <- factor(meta_clinical_sort$sample_name, levels=meta_clinical_sort$sample_name)

meta_clinical_sort$sex_char <- "M"
meta_clinical_sort$sex_char[as.numeric(meta_clinical_sort$sex) == 2] <- "F"


pt1 <- ggplot(data = meta_clinical_sort, aes(x=sample_name, y=1, fill=diagnosis)) + geom_tile() + scale_fill_manual(values=c("Control"="#fee0d2", "Prodromal"="#fc9272", "Case"="#de2d26")) + theme_bw() + theme(legend.direction="horizontal", legend.key.height=unit(0.25,"cm"),
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())
pt2 <- ggplot(data = meta_clinical_sort, aes(x=sample_name, y=1, fill=nctx_temporal_factor)) + geom_tile() + scale_fill_viridis_d(option = "plasma") + theme_bw() + theme(legend.direction="horizontal", legend.key.height=unit(0.25,"cm"),
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())
pt3 <- ggplot(data = meta_clinical_sort, aes(x=sample_name, y=1, fill=brain_stem_sn_factor)) + geom_tile() + scale_fill_viridis_d(option = "plasma") + theme_bw() + theme(legend.direction="horizontal", legend.key.height=unit(0.25,"cm"),
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())
pt4 <- ggplot(data = meta_clinical_sort, aes(x=sample_name, y=1, fill=last_mmse_test_score)) + geom_tile() +
    # scale_fill_viridis_c(trans = "reverse", option = "viridis") +
    scale_fill_gradient2(low = "#bd0026", midpoint = 24, mid = "#fd8d3c", high = "#ffffb2", space="Lab") +
    theme_bw() + theme(legend.direction="horizontal", legend.key.height=unit(0.25,"cm"),
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())
pt5 <- ggplot(data = meta_clinical_sort, aes(x=sample_name, y=1, fill=motor_updrs_score)) + geom_tile() +
    # scale_fill_viridis_c(option = "viridis") +
    scale_fill_gradient2(low = "#91cf60", midpoint = 45.5, mid = "#ffffbf", high = "#fc8d59", space="Lab") +
    theme_bw() + theme(legend.direction="horizontal", legend.key.height=unit(0.25,"cm"),
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())
pt6 <- ggplot(data = meta_clinical_sort, aes(x=sample_name, y=1, fill=sex_char)) + geom_tile() + theme_bw() + theme(legend.direction="horizontal", legend.key.height=unit(0.25,"cm"),
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())
pt7 <- ggplot(data = meta_clinical_sort, aes(x=sample_name, y=1, fill=expired_age_num)) + geom_tile() +
    # scale_fill_viridis_c(option = "mako") +
    scale_fill_distiller(palette="YlGn") +
    theme_bw() + theme(legend.direction="horizontal", legend.key.height=unit(0.25,"cm"),
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())
pt8 <- ggplot(data = meta_clinical_sort, aes(x=sample_name, y=1, fill=PMI_num)) + geom_tile() +
    # scale_fill_viridis_c(option = "mako") +
    scale_fill_gradient(low="#f2f0f7", high="#54278f", space="Lab") +
    theme_bw() + theme(legend.direction="horizontal", legend.key.height=unit(0.25,"cm"),
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())
pt9 <- ggplot(data = meta_clinical_sort, aes(x=sample_name, y=1, fill=RIN_num)) + geom_tile() +
    # scale_fill_viridis_c(option = "mako") +
    scale_fill_distiller(palette="RdPu") +
    theme_bw() + theme(legend.direction="horizontal", legend.key.height=unit(0.25,"cm"),
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())
pt10 <- ggplot(data = meta_clinical_sort, aes(x=sample_name, y=1, fill=nspots)) + geom_tile() +
    # scale_fill_viridis_c(option = "mako") +
    scale_fill_distiller(palette="GnBu") +
    theme_bw() + theme(legend.direction="horizontal", legend.key.height=unit(0.25,"cm"),
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())
pt11 <- ggplot(data = meta_clinical_sort, aes(x=sample_name, y=1, fill=mean_read_count)) + geom_tile() +
    # scale_fill_viridis_c(option = "mako") +
    # scale_fill_gradient(low="#f1a340", high="#998ec3", space="Lab") +
    # scale_fill_gradient2(midpoint = mean(meta_clinical_sort$mean_read_count), low = "#d7191c", mid = "#ffffbf", high = "#2b83ba") +
    scale_fill_gradient2(midpoint = 4000, low = "#d7191c", mid = "#ffffbf", high = "#2b83ba") +
    theme_bw() + theme(legend.direction="horizontal", legend.key.height=unit(0.25,"cm"),
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank())
pdf("st_clinical_overview_grid_09272024.pdf")
plot_grid(pt1, pt2, pt3, pt4, pt5, pt6, pt7, pt8, pt9, pt10, pt11, ncol=1, align='v')
dev.off()

# alternative 3D-like plot, showing nctx_temporal, brain_stem_sn, diagnosis
library(vcd)
temp <- meta_clinical_sort
temp$nctx_temporal_factor_disp <- as.character(temp$nctx_temporal_factor)
temp$nctx_temporal_factor_disp[is.na(temp$nctx_temporal_factor_disp)] <- "NoData"
temp$nctx_temporal_factor_disp[temp$nctx_temporal_factor_disp == "4"] <- "3"
temp$nctx_temporal_factor_disp[temp$nctx_temporal_factor_disp == "3"] <- ">=3"
temp$nctx_temporal_factor_disp <- factor(temp$nctx_temporal_factor_disp, ordered=TRUE, levels=c("0", "1", "2", ">=3", "NoData"))
pdf("st_clinical_3d_grid_10022024.pdf")
struct <- structable(~ nctx_temporal_factor_disp + brain_stem_sn_factor + diagnosis, data = temp)
mosaic(struct, data = temp, pop=F)
labeling_cells(text = as.table(struct), margin = 0)(as.table(struct))
# vcd::mosaic(diagnosis ~ nctx_temporal_factor + brain_stem_sn_factor, data = meta_clinical_sort)
dev.off()


# pdf("clinical_features_V4_12142023.pdf")
random_seed <- 101
set.seed(random_seed)
pdf("clinical_features_plotgrid_V4_09232024.pdf", width=12.5, height=10)

# MTG LB x SN LB
mtgxsn <- meta_clinical_sort %>%
            select(nctx_temporal_factor, brain_stem_sn_factor) %>%
            group_by(nctx_temporal_factor, brain_stem_sn_factor) %>%
            summarise(count=n()) %>%
            na.omit() %>%
            as.data.frame()
diagxsn <- meta_clinical_sort %>%
            select(diagnosis, brain_stem_sn_factor) %>%
            group_by(diagnosis, brain_stem_sn_factor) %>%
            summarise(count=n()) %>%
            na.omit() %>%
            as.data.frame()
diagxmtg <- meta_clinical_sort %>%
            select(diagnosis, nctx_temporal_factor) %>%
            group_by(diagnosis, nctx_temporal_factor) %>%
            summarise(count=n()) %>%
            na.omit() %>%
            as.data.frame()
mtgxmmse <- meta_clinical_sort %>%
            select(nctx_temporal_factor, last_mmse_test_score) %>%
            na.omit() %>%
            as.data.frame()
mtgxupdrs <- meta_clinical_sort %>%
            select(nctx_temporal_factor, motor_updrs_score) %>%
            na.omit() %>%
            as.data.frame()
pt11 <- ggplot(diagxsn, aes(x=brain_stem_sn_factor,
                           y=diagnosis, fill=count)) +
        geom_tile() +
        geom_text(aes(brain_stem_sn_factor,
                      diagnosis, label=count),
                      colour = "white", check_overlap = TRUE) +
        theme_bw() +
        theme(axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              legend.position = "none")
pt12 <- ggplot(meta_clinical_sort, aes(x=last_mmse_test_score,
                                       y=diagnosis)) +
        geom_boxplot(outlier.color=NA) +
        geom_jitter() +
        scale_x_reverse() +
        theme_bw() +
        theme(axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.text.y = element_blank(),
              axis.title.y = element_blank())
pt13 <- ggplot(diagxmtg, aes(x=nctx_temporal_factor,
                           y=diagnosis, fill=count)) +
        geom_tile() +
        geom_text(aes(nctx_temporal_factor,
                      diagnosis, label=count),
                      colour = "white", check_overlap = TRUE) +
        theme_bw() +
        theme(axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              legend.position = "none")
pt21 <- ggplot(meta_clinical_sort, aes(x=brain_stem_sn_factor,
                                       y=motor_updrs_score)) +
        geom_boxplot(outlier.color=NA) +
        geom_jitter() +
        theme_bw() +
        theme(axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              legend.position = "none")
pt22 <- ggplot(meta_clinical_sort, aes(x=last_mmse_test_score,
                                       y=motor_updrs_score)) +
        stat_poly_line() +
        stat_poly_eq(use_label(c("eq", "R2"))) +
        geom_point() +
        scale_x_reverse() +
        theme_bw() +
        theme(axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.text.y = element_blank(),
              axis.title.y = element_blank())
pt23 <- ggplot(mtgxupdrs, aes(x=nctx_temporal_factor,
                              y=motor_updrs_score)) +
        geom_boxplot(outlier.color=NA) +
        geom_jitter() +
        theme_bw() +
        theme(axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.text.y = element_blank(),
              axis.title.y = element_blank())
pt31 <- ggplot(mtgxsn, aes(x=brain_stem_sn_factor,
                           y=nctx_temporal_factor, fill=count)) +
        geom_tile() +
        geom_text(aes(brain_stem_sn_factor,
                      nctx_temporal_factor, label=count),
                      colour = "white", check_overlap = TRUE) +
        theme_bw() +
        theme(legend.position = "none")
pt32 <- ggplot(mtgxmmse, aes(x=last_mmse_test_score,
                             y=nctx_temporal_factor)) +
        geom_boxplot(outlier.color=NA) +
        geom_jitter() +
        scale_x_reverse() +
        theme_bw() +
        theme(axis.text.y = element_blank(),
              axis.title.y = element_blank())
plot_grid(pt11, pt12, pt13, pt21, pt22, pt23, pt31, pt32, NULL, ncol=3, align='v')
dev.off()

# sp_pt <- ggpairs(meta_clinical[, c("diagnosis", "brain_stem_sn", "nctx_temporal_fix", "motor_updrs_score", "last_mmse_test_score")]) + theme_bw()
# print(sp_pt)
# sp_pt <- ggpairs(meta_clinical[, c("diagnosis", "brain_stem_sn_factor", "nctx_temporal_fix_factor", "motor_updrs_score", "last_mmse_test_score")]) + theme_bw()
# print(sp_pt)
# dev.off()

# ---------------------------
# old code below - remove it (9/24/2024)







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
