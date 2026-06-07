library(Seurat)
library(testit)
library(stringr)
library(dplyr)
library(tibble)
library(tidyr)

library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(corrplot)

# comb_seurat <- readRDS("/mnt/data0/projects/ASAP/ST/clustering_pipeline/layer_smoothing_V4_seed101_11302023_replication_test/comb_seurat_post_smoothed_label_s5_11302023.rds")
# comb_seurat@meta.data$spot_id <- rownames(comb_seurat@meta.data)
# comb_seurat@meta.data$nctx_temporal_fix <- comb_seurat@meta.data$nctx_temporal
# comb_seurat@meta.data$nctx_temporal_fix[comb_seurat@meta.data$nctx_temporal_fix > 3] <- 3
# all_clinical <- read.csv("/mnt/data0/projects/ASAP/ST/data/all_samples_clinical_data_11062022.csv")



# meta1 <- comb_seurat@meta.data[, c("orig.ident", "sample_name", "batch", "PMI_num", "RIN_num", "sex", "expired_age_num", "diagnosis", "nctx_temporal")]
# meta1 <- meta1[match(unique(meta1$sample_name), meta1$sample_name),]


# meta2 <- merge(meta1, all_clinical, by.x="sample_name", by.y="CaseID")

# nctx_temporal
#  0  1  2  3  4
# 58 22  6  8  1

#       nctx_temporal
# gender  0  1  2  3  4
#      1 34 16  5  6  1
#      2 24  6  1  2  0

#              Braak.score
# nctx_temporal  0  I II III IV
#             0  1  8 11  19 19
#             1  0  1  3   9  9
#             2  0  0  1   2  3
#             3  0  0  3   4  1
#             4  0  0  0   1  0

# [14] "gender"
# [15] "expired_age"
# [23] "AD"
# [27] "Unified.LB.Stage"
# [57] "Braak.score"
# [58] "NIA.R.for.AD"
# [60] "neurological_dx"
# [61] "neurological_dx_years"
# [62] "neurological_dx_age"
# [63] "dementia_yn"
# [64] "dementia_years"
# [65] "demential_age"


#
# ---
#

# boxplot for major and minor cell type proportion vs mmse score
comb_seurat <- readRDS("/mnt/data0/projects/ASAP/ST/clustering_pipeline/layer_smoothing_V4_seed101_11302023_replication_test/comb_seurat_post_smoothed_label_s5_11302023.rds")
comb_seurat@meta.data$spot_id <- rownames(comb_seurat@meta.data)
comb_seurat@meta.data$nctx_temporal_fix <- comb_seurat@meta.data$nctx_temporal
comb_seurat@meta.data$nctx_temporal_fix[comb_seurat@meta.data$nctx_temporal_fix > 3] <- 3

# cell_prop_dir <- "/mnt/data0/projects/ASAP/ST/st_sc_deconvolution/results/cell2location_map"
cell_prop_dir <- "/mnt/data0/projects/ASAP/ST/st_sc_deconvolution/results_V4_02202024/cell2location_map"
# cell_type_strs <- c("Astro", "Endo", "GABA", "GLU", "Micro", "OPC", "Oligo", "TEMRA")
# cell_type_strs <- paste(cell_type_strs, collapse="|")
# easier to ignore the first N columns
num_non_cell_type_cols <- 35
last_non_cell_type_col <- "X_scvi_labels"

comb_meta_data <- comb_seurat@meta.data
comb_meta_cell2loc <- c()
all_cell_type_cols <- c() # these can be numerous and not known beforehand
for(cellprop_file in list.files(cell_prop_dir, pattern="st_cellprop_table.*\\.csv")) {
    print(cellprop_file)
    cell2loc_sample <- gsub("\\.csv", "", gsub("st_cellprop_table_", "", cellprop_file))
    cellprops <- read.csv(file.path(cell_prop_dir, cellprop_file))
    colnames(cellprops)[colnames(cellprops) == "X"] <- "spot_id"
    
    # these may be slightly different if some cell types are missing
    # cell_type_cols <- colnames(cellprops)[grepl(cell_type_strs, colnames(cellprops))]
    assert(colnames(cellprops)[num_non_cell_type_cols] == last_non_cell_type_col)
    cell_type_cols <- colnames(cellprops)[(num_non_cell_type_cols+1):ncol(cellprops)]
    cellprops <- cellprops[, c("spot_id", cell_type_cols)]
    all_cell_type_cols <- sort(unique(c(all_cell_type_cols, cell_type_cols)))
    
    # comb_sub <- subset(comb_seurat, subset = sample_name == cell2loc_sample)
    # comb_sub@meta.data$Row.names <- rownames(comb_sub@meta.data)
    comb_meta_sub <- comb_meta_data[comb_meta_data$sample_name == cell2loc_sample, ]
    comb_meta_sub <- merge(comb_meta_sub, cellprops, by="spot_id", all.x=TRUE, sort=FALSE)
    
    comb_meta_cell2loc <- dplyr::bind_rows(comb_meta_cell2loc, comb_meta_sub)
    
    # temp <- merge(comb_sub@meta.data, cellprops, by="new_spot_ids")
    # temp <- merge(comb_sub@meta.data, cellprops, by="spot_ids")
    # rownames(temp) <- temp$Row.names
    # temp <- temp[ colnames(comb_sub) ,]
    # comb_sub@meta.data <- temp
}
rownames(comb_meta_cell2loc) <- comb_meta_cell2loc$spot_id
comb_meta_data <- merge(comb_meta_data, comb_meta_cell2loc[, c("spot_id", all_cell_type_cols)], by="spot_id", all.x=TRUE, sort=FALSE)
rownames(comb_meta_data) <- comb_meta_data$spot_id
# re-sort the meta rows to match data columns
comb_meta_data <- comb_meta_data[colnames(comb_seurat), ]
assert(all(colnames(comb_seurat) == rownames(comb_meta_data)))
assert(length(unique(rownames(comb_meta_data))) == nrow(comb_meta_data))

normalize_counts_per_spot <- TRUE
# normalize_counts_per_spot <- FALSE

# NOTE: for samples with no sc deconvolution data, the entries will be kept, but because of the divide by 0 the rows will have all NaN values
if(normalize_counts_per_spot) {
    temp <- comb_meta_data[, all_cell_type_cols]
    temp[is.na(temp)] <- 0
    temp <- sweep(temp, 1, rowSums(temp), `/`)
    comb_meta_data[, all_cell_type_cols] <- temp[, all_cell_type_cols]
}
comb_seurat@meta.data <- comb_meta_data

# if needed, grab additional variables from clinical meta table
# all_clinical <- read.csv("/mnt/data0/projects/ASAP/ST/data/all_samples_clinical_data_11062022.csv")
#
# meta2 <- merge(comb_seurat@meta.data, all_clinical[, c("CaseID", "Braak.score", "Unified.LB.Stage")], by.x="sample_name", by.y="CaseID", all.x=TRUE)
# rownames(meta2) <- meta2$spot_id
# assert(all(rownames(comb_seurat@meta.data) == rownames(meta2)))
#
# comb_seurat@meta.data <- meta2

# 6/6/2024: just run cont+ilbd vs case
meta2 <- comb_seurat@meta.data

# traits to regress: nctx_temporal_fix, last_mmse_test_score, Braak.score
# remove samples with no sc, and pseudobulk
meta2 <- meta2[rowSums(is.na(meta2[, all_cell_type_cols])) != length(all_cell_type_cols), ]
temp_meta2 <- meta2[, all_cell_type_cols]
temp_meta2[is.na(temp_meta2)] <- 0
meta2[, all_cell_type_cols] <- temp_meta2
# create major clusters
meta2$Pericytes_major <- rowSums(meta2[, c("Pericytes.1", "Pericytes.2")])
gaba_cell_types <- all_cell_type_cols[grepl("^ADARB2|^GAD1|LAMP5\\.CA1|^LAMP5\\.LCP2|LAMP5\\.NMBR|^PAX6|^PVALB|^SST|^VIP", all_cell_type_cols)]
meta2$GABA_major <- rowSums(meta2[, gaba_cell_types])
glu_cell_types <- all_cell_type_cols[grepl("^FEZF2|LAMP5\\.GREB1L|LAMP5\\.PARD3|^LINC00507|^RORB|^SLC17A7|^THEMIS", all_cell_type_cols)]
meta2$GLU_major <- rowSums(meta2[, glu_cell_types])

# pseudobulk: average proportion per sample per layer
# meta3 <- meta2[, c("sample_name", "batch", "sex", "expired_age_num", "PMI_num", "RIN_num", "smoothed_label_s5", "nctx_temporal_fix", "last_mmse_test_score", "Braak.score", all_cell_type_cols, "Pericytes_major", "GABA_major", "GLU_major")]
# meta3 <- meta3 %>% group_by(sample_name, batch, sex, expired_age_num, PMI_num, RIN_num, smoothed_label_s5, nctx_temporal_fix, last_mmse_test_score, Braak.score) %>% summarize_all(mean) %>% as.data.frame()
# meta3$Braak_num <- as.numeric(factor(meta3$Braak.score, levels=c("0", "I", "II", "III", "IV"))) - 1

# st_layers <- sort(unique(meta3$smoothed_label_s5))
# lnres_all <- c()
# for(quant_trait in c("nctx_temporal_fix", "last_mmse_test_score", "Braak_num")) {
#     for(st_layer in st_layers) {
#         lnres <- c()
#         meta3_sub <- meta3[meta3$smoothed_label_s5 == st_layer, ]
#         meta3_sub$expired_age_num <- scale(meta3_sub$expired_age_num)
#         meta3_sub$PMI_num <- scale(meta3_sub$PMI_num)
#         meta3_sub$RIN_num <- scale(meta3_sub$RIN_num)
#         meta3_sub[, quant_trait] <- scale(meta3_sub[, quant_trait])
#         for(ct_to_test in c(all_cell_type_cols, "Pericytes_major", "GABA_major", "GLU_major")) {
#             meta3_sub$ct <- meta3_sub[, ct_to_test]
#             meta3_sub$ct <- scale(meta3_sub$ct)
#             mdl <- lm(data = meta3_sub, as.formula(sprintf("ct ~ batch + sex + expired_age_num + PMI_num + RIN_num + %s", quant_trait)))
#             res <- summary(mdl)$coefficients[quant_trait, ]
#             lnres <- rbind(lnres, data.frame("cell_type" = ct_to_test,
#                                              "estimate" = res[["Estimate"]],
#                                              "st_err" = res[["Std. Error"]],
#                                              "pval" = res[["Pr(>|t|)"]]))
#         }
#         lnres$padj <- p.adjust(lnres$pval, method="BH")
#         lnres$st_layer <- st_layer
#         lnres$quant_trait <- quant_trait
#         lnres_all <- rbind(lnres_all, lnres)
#     }
# }
# write.csv(lnres_all, "lnres_all_st_sc_cellprops_quant_traits_scaled_05302024.csv", row.names=FALSE)

meta2$diagnosis_contilbd_case <- as.character(meta2$diagnosis)
meta2[meta2$diagnosis == "Control", "diagnosis_contilbd_case"] <- "contilbd"
meta2[meta2$diagnosis == "ILBD", "diagnosis_contilbd_case"] <- "contilbd"
meta2[meta2$diagnosis == "Case", "diagnosis_contilbd_case"] <- "case"
meta2$diagnosis_contilbd_case <- factor(meta2$diagnosis_contilbd_case, levels=c("contilbd", "case"))
meta3 <- meta2[, c("sample_name", "batch", "sex", "expired_age_num", "PMI_num", "RIN_num", "smoothed_label_s5", "nctx_temporal_fix", "diagnosis", "diagnosis_contilbd_case", all_cell_type_cols, "Pericytes_major", "GABA_major", "GLU_major")]
meta3 <- meta3 %>% group_by(sample_name, batch, sex, expired_age_num, PMI_num, RIN_num, smoothed_label_s5, nctx_temporal_fix, diagnosis, diagnosis_contilbd_case) %>% summarize_all(mean) %>% as.data.frame()

st_layers <- sort(unique(meta3$smoothed_label_s5))
lnres_all <- c()
lnres_wilcox <- c()
for(binary_trait in c("diagnosis_contilbd_case")) {
    for(st_layer in st_layers) {
        lnres <- c()
        lnresw <- c()
        meta3_sub <- meta3[meta3$smoothed_label_s5 == st_layer, ]
        meta3_sub$expired_age_num <- scale(meta3_sub$expired_age_num)
        meta3_sub$PMI_num <- scale(meta3_sub$PMI_num)
        meta3_sub$RIN_num <- scale(meta3_sub$RIN_num)
        # meta3_sub[, quant_trait] <- scale(meta3_sub[, quant_trait])
        assert(is.factor(meta3_sub[, binary_trait]))
        for(ct_to_test in c(all_cell_type_cols, "Pericytes_major", "GABA_major", "GLU_major")) {
            meta3_sub$ct <- meta3_sub[, ct_to_test]
            meta3_sub$ct <- scale(meta3_sub$ct)
            mdl <- lm(data = meta3_sub, as.formula(sprintf("ct ~ batch + sex + expired_age_num + PMI_num + RIN_num + %s", binary_trait)))
            highest_factor_level <-  levels(meta3_sub[, binary_trait])[max(as.numeric(meta3_sub[, binary_trait]))] # e.g. diagnosis_contilbd_case: "case"
            res <- summary(mdl)$coefficients[paste0(binary_trait, highest_factor_level), ]
            lnres <- rbind(lnres, data.frame("cell_type" = ct_to_test,
                                             "estimate" = res[["Estimate"]],
                                             "st_err" = res[["Std. Error"]],
                                             "pval" = res[["Pr(>|t|)"]]))
            # wilcoxon rank-sum test
            resw <- wilcox.test(as.formula(sprintf("ct ~ %s", binary_trait)), data = meta3_sub, exact = FALSE)
            lnresw <- rbind(lnresw, data.frame("cell_type" = ct_to_test,
                                             "statistic" = resw[["statistic"]],
                                             "pval" = resw[["p.value"]]))
            
        }
        lnres$padj <- p.adjust(lnres$pval, method="BH")
        lnres$st_layer <- st_layer
        lnres$binary_trait <- binary_trait
        lnres_all <- rbind(lnres_all, lnres)
        
        lnresw$padj <- p.adjust(lnresw$pval, method="BH")
        lnresw$st_layer <- st_layer
        lnresw$binary_trait <- binary_trait
        lnres_wilcox <- rbind(lnres_wilcox, lnresw)
    }
}
write.csv(lnres_all, sprintf("lnres_all_st_sc_cellprops_diag_contilbd_case_norm%s_06072024.csv", normalize_counts_per_spot), row.names=FALSE)
write.csv(lnres_wilcox, sprintf("lnres_wilcox_st_sc_cellprops_diag_contilbd_case_norm%s_06072024.csv", normalize_counts_per_spot), row.names=FALSE)




plot_cellprop_vs_clinical <- function(meta3, st_layer, clinical_feat, quant_feature, ct_to_test) {
    meta3_sub <- meta3[meta3$smoothed_label_s5 == st_layer, ]
    meta3_sub$ct <- meta3_sub[, ct_to_test]
    meta3_sub$clinical_feat <- meta3_sub[, clinical_feat]
    if(quant_feature) { # scatterplot
        st_pt <- ggplot(meta3_sub, aes(x=clinical_feat, y=ct)) + geom_point() + geom_smooth(method = "lm") + theme_bw() + ggtitle(paste(st_layer, clinical_feat, ct_to_test))
    } else { # boxplot
        meta3_sub <- meta3_sub[!is.na(meta3_sub$clinical_feat), ]
        meta3_sub$clinical_feat <- factor(meta3_sub$clinical_feat, levels=sort(unique(meta3_sub$clinical_feat)))
        st_pt <- ggplot(meta3_sub, aes(x=clinical_feat, y=ct)) + geom_boxplot(outlier.color=NA) + geom_jitter() + theme_bw() + ggtitle(paste(st_layer, clinical_feat, ct_to_test))
    }
    print(st_pt)
}

pdf("lnres_all_st_sc_cellprops_quant_plots_05302024.pdf")
plot_cellprop_vs_clinical(meta3, st_layer="Layer 1", clinical_feat="nctx_temporal_fix", quant_feature=FALSE, ct_to_test="RORB.LAMA2")
plot_cellprop_vs_clinical(meta3, st_layer="Layer 4", clinical_feat="nctx_temporal_fix", quant_feature=FALSE, ct_to_test="THEMIS.C1QL3")
plot_cellprop_vs_clinical(meta3, st_layer="Layer 5", clinical_feat="nctx_temporal_fix", quant_feature=FALSE, ct_to_test="THEMIS.C1QL3")
plot_cellprop_vs_clinical(meta3, st_layer="Layer 2", clinical_feat="last_mmse_test_score", quant_feature=FALSE, ct_to_test="VIP.SERPINF1")
plot_cellprop_vs_clinical(meta3, st_layer="Layer 3", clinical_feat="last_mmse_test_score", quant_feature=FALSE, ct_to_test="VIP.SERPINF1")
plot_cellprop_vs_clinical(meta3, st_layer="Layer 4", clinical_feat="last_mmse_test_score", quant_feature=FALSE, ct_to_test="VIP.SERPINF1")
plot_cellprop_vs_clinical(meta3, st_layer="Layer 5", clinical_feat="last_mmse_test_score", quant_feature=FALSE, ct_to_test="LAMP5.GREB1L")
plot_cellprop_vs_clinical(meta3, st_layer="Layer 6", clinical_feat="last_mmse_test_score", quant_feature=FALSE, ct_to_test="LAMP5.GREB1L")
plot_cellprop_vs_clinical(meta3, st_layer="Layer 5", clinical_feat="last_mmse_test_score", quant_feature=FALSE, ct_to_test="LAMP5.PARD3")
plot_cellprop_vs_clinical(meta3, st_layer="Layer 6", clinical_feat="last_mmse_test_score", quant_feature=FALSE, ct_to_test="LAMP5.PARD3")
plot_cellprop_vs_clinical(meta3, st_layer="Layer 5", clinical_feat="last_mmse_test_score", quant_feature=FALSE, ct_to_test="Pericytes.2")
plot_cellprop_vs_clinical(meta3, st_layer="Layer 6", clinical_feat="last_mmse_test_score", quant_feature=FALSE, ct_to_test="Pericytes.2")
plot_cellprop_vs_clinical(meta3, st_layer="Layer 5", clinical_feat="last_mmse_test_score", quant_feature=FALSE, ct_to_test="SLC17A7.NEFM")
plot_cellprop_vs_clinical(meta3, st_layer="Layer 6", clinical_feat="last_mmse_test_score", quant_feature=FALSE, ct_to_test="SLC17A7.NEFM")
plot_cellprop_vs_clinical(meta3, st_layer="Layer 5", clinical_feat="last_mmse_test_score", quant_feature=FALSE, ct_to_test="VIP.SERPINF1")
plot_cellprop_vs_clinical(meta3, st_layer="Layer 6", clinical_feat="last_mmse_test_score", quant_feature=FALSE, ct_to_test="VIP.SERPINF1")
dev.off()


pdf("lnres_all_st_sc_cellprops_quant_plots_major_05302024.pdf")
for(ct_to_test in c("Pericytes_major", "GLU_major", "GABA_major")) {
    for(clinical_feat in c("nctx_temporal_fix", "last_mmse_test_score")) {
        for(st_layer in st_layers) {
            plot_cellprop_vs_clinical(meta3, st_layer=st_layer, clinical_feat=clinical_feat, quant_feature=FALSE, ct_to_test=ct_to_test)
        }
    }
}
dev.off()










# bin by contour, and plot differences in average gene fraction
# contilbd
# case
# case - contilbd
# see ~/Downloads/spatial_contours_KLD_contour_test.R
random_seed <- 101
contour_score_table <- list(
"Layer 1" = 1, "Layer 2" = 2, "Layer 3" = 3, "Layer 4" = 4,
"Layer 5" = 5, "Layer 6" = 6, "WM" = 7)
# contour_score_table <- list(
# "Layer 1" = 7, "Layer 2" = 6, "Layer 3" = 5, "Layer 4" = 4,
# "Layer 5" = 3, "Layer 6" = 2, "WM" = 1)
# scores of 0 are omitted from the kernel weighted average
comb_seurat@meta.data[, "contour_score"] <- 0
for(st_layer in names(contour_score_table)) {
    comb_seurat@meta.data[ (!is.na(comb_seurat@meta.data$smoothed_label_s5) & comb_seurat@meta.data$smoothed_label_s5 == st_layer), "contour_score" ] <- contour_score_table[[st_layer]]
}


cluster_labels <- comb_seurat@meta.data[, c("orig.ident", "contour_score")]
# rbf_gamma <- 0.1
# rbf_gamma <- 0.1 # gamma: low value means more long-range influence
nnbrs <- 18
rbf_gamma <- 0.05
# pdf(sprintf("layer_contour_map_gamma%s_02122024.pdf", rbf_gamma))
for (img in names(comb_seurat@images)) {
    print(paste("Processing contour score: ", img))
    set.seed(random_seed)

    meta_row <- rownames(comb_seurat@images[[img]]@coordinates)[1]
    sample_name_filt <- comb_seurat@meta.data[meta_row, "orig.ident"]
    diagnosis_filt <- comb_seurat@meta.data[meta_row, "diagnosis"]
    
    coord_table <- comb_seurat@images[[img]]@coordinates
    coord_table <- merge(coord_table, cluster_labels, by=0)
    rownames(coord_table) <- coord_table$Row.names

    # hexagonal
    coord_table$euclid_x <- (coord_table$row %% 2) / 2 + (coord_table$col - (coord_table$col %% 2)) / 2
    coord_table$euclid_y <- sqrt(3) * ( (coord_table$row %% 2) / 2 + (coord_table$row - (coord_table$row %% 2)) / 2)
    coord_dists <- dist(coord_table[, c("euclid_x", "euclid_y")])
    coord_dists <- as.matrix(coord_dists, nrow=nrow(coord_table))
    adjA <- exp( - rbf_gamma * (coord_dists**2))
    diag(adjA) <- 0
    # normalize rows to sum to 1
    # set adjacency of na values to 0
    zero_score_rows <- which(coord_table$contour_score == 0)
    adjA[, zero_score_rows] <- 0
    adjA <- sweep(adjA, 1, rowSums(adjA), FUN="/")
    
    coord_table$contour_score2 <- adjA %*% coord_table$contour_score
    comb_seurat@meta.data[rownames(coord_table), "contour_score2"] <- coord_table$contour_score2
}
comb_seurat@meta.data$contour_bin <- cut(comb_seurat@meta.data$contour_score2, breaks=35)

meta2 <- comb_seurat@meta.data

# remove samples with no sc, and pseudobulk
meta2 <- meta2[rowSums(is.na(meta2[, all_cell_type_cols])) != length(all_cell_type_cols), ]
temp_meta2 <- meta2[, all_cell_type_cols]
temp_meta2[is.na(temp_meta2)] <- 0
meta2[, all_cell_type_cols] <- temp_meta2
# create major clusters
meta2$Pericytes_major <- rowSums(meta2[, c("Pericytes.1", "Pericytes.2")])
gaba_cell_types <- all_cell_type_cols[grepl("^ADARB2|^GAD1|LAMP5\\.CA1|^LAMP5\\.LCP2|LAMP5\\.NMBR|^PAX6|^PVALB|^SST|^VIP", all_cell_type_cols)]
meta2$GABA_major <- rowSums(meta2[, gaba_cell_types])
glu_cell_types <- all_cell_type_cols[grepl("^FEZF2|LAMP5\\.GREB1L|LAMP5\\.PARD3|^LINC00507|^RORB|^SLC17A7|^THEMIS", all_cell_type_cols)]
meta2$GLU_major <- rowSums(meta2[, glu_cell_types])

meta2$diagnosis_contilbd_case <- as.character(meta2$diagnosis)
meta2[meta2$diagnosis == "Control", "diagnosis_contilbd_case"] <- "contilbd"
meta2[meta2$diagnosis == "ILBD", "diagnosis_contilbd_case"] <- "contilbd"
meta2[meta2$diagnosis == "Case", "diagnosis_contilbd_case"] <- "case"
meta2$diagnosis_contilbd_case <- factor(meta2$diagnosis_contilbd_case, levels=c("contilbd", "case"))
meta3 <- meta2[, c("sample_name", "diagnosis_contilbd_case", "contour_bin", all_cell_type_cols, "Pericytes_major", "GABA_major", "GLU_major")]
# # meta3 <- meta3 %>% group_by(sample_name, diagnosis_contilbd_case, contour_bin) %>% summarize_all(mean) %>% as.data.frame()
# meta3_contilbd <- meta3 %>% filter(diagnosis_contilbd_case == "contilbd") %>%
#                         select (-c(sample_name, diagnosis_contilbd_case)) %>%
#                         group_by(contour_bin) %>% summarize_all(mean) %>% as.data.frame()
# meta3_case <- meta3 %>% filter(diagnosis_contilbd_case == "case") %>%
#                         select (-c(sample_name, diagnosis_contilbd_case)) %>%
#                         group_by(contour_bin) %>% summarize_all(mean) %>% as.data.frame()
# group by sample first, take mean, then take mean again


meta3_contilbd <- meta3 %>% filter(diagnosis_contilbd_case == "contilbd") %>%
                    select (-c(diagnosis_contilbd_case)) %>%
                    group_by(across(all_of(c("sample_name", "contour_bin")))) %>%
                    summarize_all(mean) %>%
                    ungroup() %>%
                    select(-c(sample_name)) %>%
                    group_by(contour_bin) %>%
                    summarize_all(mean) %>% as.data.frame()
meta3_case <- meta3 %>% filter(diagnosis_contilbd_case == "case") %>%
                    select (-c(diagnosis_contilbd_case)) %>%
                    group_by(across(all_of(c("sample_name", "contour_bin")))) %>%
                    summarize_all(mean) %>%
                    ungroup() %>%
                    select(-c(sample_name)) %>%
                    group_by(contour_bin) %>%
                    summarize_all(mean) %>% as.data.frame()
assert(all(rownames(meta3_contilbd) == rownames(meta3_case)))
assert(all(colnames(meta3_contilbd) == colnames(meta3_case)))
assert(all(meta3_contilbd$contour_bin == meta3_case$contour_bin))
meta3_diff <- meta3_case
meta3_diff[, 2:ncol(meta3_diff)] <- meta3_case[, 2:ncol(meta3_case)] - meta3_contilbd[, 2:ncol(meta3_contilbd)]

meta3_contilbd <- meta3_contilbd %>%
  pivot_longer(!contour_bin, names_to = "cell_type", values_to = "mean_pct") %>% as.data.frame()
meta3_case <- meta3_case %>%
  pivot_longer(!contour_bin, names_to = "cell_type", values_to = "mean_pct") %>% as.data.frame()
meta3_diff <- meta3_diff %>%
  pivot_longer(!contour_bin, names_to = "cell_type", values_to = "mean_pct") %>% as.data.frame()

major_cell_types <- c("Astrocytes", "Endothelial.Cells", "Fibroblast.Like.Cells", "GABA_major", "GLU_major", "Microglia", "Oligodendrocytes", "OPCs", "Pericytes_major", "T.Cells")

pdf(sprintf("layer_contour_geomcol_major_gamma%s_contilbd_norm%s_06072024.pdf", rbf_gamma, normalize_counts_per_spot), width=10, height=5)
num_levels <- length(levels(meta3_contilbd$contour_bin))
pt <- ggplot() + geom_col(data = meta3_contilbd[meta3_contilbd$cell_type %in% major_cell_types, ], aes(x = contour_bin, y = mean_pct, fill = cell_type), position = "stack") + theme_bw() +
    scale_x_discrete(breaks = levels(meta3_contilbd$contour_bin)[c(1, seq(from=5, to=num_levels, by=5))])
print(pt)
dev.off()
pdf(sprintf("layer_contour_geomcol_major_gamma%s_case_norm%s_06072024.pdf", rbf_gamma, normalize_counts_per_spot), width=10, height=5)
num_levels <- length(levels(meta3_case$contour_bin))
pt <- ggplot() + geom_col(data = meta3_case[meta3_case$cell_type %in% major_cell_types, ], aes(x = contour_bin, y = mean_pct, fill = cell_type), position = "stack") + theme_bw() +
    scale_x_discrete(breaks = levels(meta3_case$contour_bin)[c(1, seq(from=5, to=num_levels, by=5))])
print(pt)
dev.off()
pdf(sprintf("layer_contour_geomcol_major_gamma%s_ci_c_diff_norm%s_06072024.pdf", rbf_gamma, normalize_counts_per_spot), width=10, height=5)
num_levels <- length(levels(meta3_diff$contour_bin))
pt <- ggplot() + geom_col(data = meta3_diff[meta3_diff$cell_type %in% major_cell_types, ], aes(x = contour_bin, y = mean_pct, fill = cell_type), position = "stack") + theme_bw() +
    scale_x_discrete(breaks = levels(meta3_diff$contour_bin)[c(1, seq(from=5, to=num_levels, by=5))])
print(pt)
dev.off()




meta3 <- meta2[, c("sample_name", "batch", "sex", "expired_age_num", "PMI_num", "RIN_num", "contour_bin", "nctx_temporal_fix", "diagnosis", "diagnosis_contilbd_case", all_cell_type_cols, "Pericytes_major", "GABA_major", "GLU_major")]
meta3 <- meta3 %>% group_by(sample_name, batch, sex, expired_age_num, PMI_num, RIN_num, contour_bin, nctx_temporal_fix, diagnosis, diagnosis_contilbd_case) %>% summarize_all(mean) %>% as.data.frame()
st_contour_bins <- sort(unique(meta3$contour_bin))
lnres_all <- c()
lnres_wilcox <- c()
for(binary_trait in c("diagnosis_contilbd_case")) {
    for(st_bin in st_contour_bins) {
        lnres <- c()
        lnresw <- c()
        meta3_sub <- meta3[meta3$contour_bin == st_bin, ]
        meta3_sub$expired_age_num <- scale(meta3_sub$expired_age_num)
        meta3_sub$PMI_num <- scale(meta3_sub$PMI_num)
        meta3_sub$RIN_num <- scale(meta3_sub$RIN_num)
        # meta3_sub[, quant_trait] <- scale(meta3_sub[, quant_trait])
        assert(is.factor(meta3_sub[, binary_trait]))
        for(ct_to_test in c(all_cell_type_cols, "Pericytes_major", "GABA_major", "GLU_major")) {
            meta3_sub$ct <- meta3_sub[, ct_to_test]
            # rarely, there will be only zero values for a cell type in a particular bin
            if(sum(meta3_sub$ct) == 0) {
                next
            }
            meta3_sub$ct <- scale(meta3_sub$ct)
            mdl <- lm(data = meta3_sub, as.formula(sprintf("ct ~ batch + sex + expired_age_num + PMI_num + RIN_num + %s", binary_trait)))
            highest_factor_level <-  levels(meta3_sub[, binary_trait])[max(as.numeric(meta3_sub[, binary_trait]))] # e.g. diagnosis_contilbd_case: "case"
            res <- summary(mdl)$coefficients[paste0(binary_trait, highest_factor_level), ]
            lnres <- rbind(lnres, data.frame("cell_type" = ct_to_test,
                                             "estimate" = res[["Estimate"]],
                                             "st_err" = res[["Std. Error"]],
                                             "pval" = res[["Pr(>|t|)"]]))
            # wilcoxon rank-sum test
            resw <- wilcox.test(as.formula(sprintf("ct ~ %s", binary_trait)), data = meta3_sub, exact = FALSE)
            lnresw <- rbind(lnresw, data.frame("cell_type" = ct_to_test,
                                             "statistic" = resw[["statistic"]],
                                             "pval" = resw[["p.value"]]))
            
        }
        lnres$padj <- p.adjust(lnres$pval, method="BH")
        lnres$st_contour_bin <- st_bin
        lnres$binary_trait <- binary_trait
        lnres_all <- rbind(lnres_all, lnres)
        
        lnresw$padj <- p.adjust(lnresw$pval, method="BH")
        lnresw$st_contour_bin <- st_bin
        lnresw$binary_trait <- binary_trait
        lnres_wilcox <- rbind(lnres_wilcox, lnresw)
    }
}
write.csv(lnres_all, sprintf("lnres_all_st_sc_cellprops_contourbin_diag_contilbd_case_norm%s_06072024.csv", normalize_counts_per_spot), row.names=FALSE)
write.csv(lnres_wilcox, sprintf("lnres_wilcox_st_sc_cellprops_contourbin_diag_contilbd_case_norm%s_06072024.csv", normalize_counts_per_spot), row.names=FALSE)
