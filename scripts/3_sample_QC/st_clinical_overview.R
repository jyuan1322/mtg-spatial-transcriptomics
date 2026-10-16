# Clinical overview panels in Fig. 1 and per-sample metadata.

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

# input parameters
comb_seurat_file <- "/mnt/data/projects/ASAP/ST/clustering_pipeline/layer_smoothing_V4_seed101_11302023_replication_test/comb_seurat_post_smoothed_label_s5_11302023.rds"

# output parameters
clinical_overview_grid_pdf <- "results/3_sample_QC/st_clinical_overview_grid.pdf"  # Fig. 1b
clinical_features_pdf <- "results/3_sample_QC/clinical_features_plotgrid_V4.pdf"   # Fig. 1d
clinical_mosaic_pdf <- "st_clinical_3d_grid.pdf"                                    # exploratory, not in results
clinical_meta_csv <- "/path/to/data/st_clinical_meta_data_V4_96samples.csv"        # input to st_clinical_overview_3dplot.py


# Comparison of clinical features
comb_seurat <- readRDS(comb_seurat_file)

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
pdf(clinical_overview_grid_pdf)
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
pdf(clinical_mosaic_pdf)
struct <- structable(~ nctx_temporal_factor_disp + brain_stem_sn_factor + diagnosis, data = temp)
mosaic(struct, data = temp, pop=F)
labeling_cells(text = as.table(struct), margin = 0)(as.table(struct))
# vcd::mosaic(diagnosis ~ nctx_temporal_factor + brain_stem_sn_factor, data = meta_clinical_sort)
dev.off()


# pdf("clinical_features_V4_12142023.pdf")
random_seed <- 101
set.seed(random_seed)
pdf(clinical_features_pdf, width=12.5, height=10)

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
