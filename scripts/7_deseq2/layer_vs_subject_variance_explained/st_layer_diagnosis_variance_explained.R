library(Seurat)
library(stringr)
library(testit)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggbeeswarm)

comb_seurat <- readRDS("/mnt/data0/projects/ASAP/ST/clustering_pipeline/layer_smoothing_V4_seed101_11302023_replication_test/comb_seurat_post_smoothed_label_s5_11302023.rds")
comb_seurat@meta.data$spot_id <- rownames(comb_seurat@meta.data)
comb_seurat@meta.data$nctx_temporal_fix <- comb_seurat@meta.data$nctx_temporal
comb_seurat@meta.data$nctx_temporal_fix[comb_seurat@meta.data$nctx_temporal_fix > 3] <- 3

agg_seurat <- AggregateExpression(comb_seurat, slot="counts", group.by = c("sample_name", "smoothed_label_s5"), return.seurat=FALSE)
agg_seurat <- agg_seurat$Spatial
# # keep only genes which are nonzero in at least 20 spots in at least 30 samples
# num_nonzero <- rowSums(agg_seurat > 20)
# vargenes <- names(num_nonzero[num_nonzero > 30])
# Use the roughly 2000 spatial-filtered variable features
# (don't perform correlation over a huge number of variables
vargenes <- VariableFeatures(comb_seurat)
# lognormalize the count matrix
agg_seurat_ln <- sweep(agg_seurat,2,colSums(agg_seurat),`/`)
agg_seurat_ln <- log1p(agg_seurat_ln * 10000)
agg_seurat_filt <- agg_seurat_ln[vargenes, ]
agg_corr <- cor(agg_seurat_filt, method="spearman") # correlate samples in columns
agg_corr <- as.data.frame(agg_corr)
agg_corr$layer_sampleA <- rownames(agg_corr)

agg_corr <- agg_corr %>% as.data.frame() %>% pivot_longer(!layer_sampleA, names_to="layer_sampleB", values_to="spearman_corr") %>% as.data.frame()
agg_corr[, c("sampleA", "layerA")] <- str_split_fixed(agg_corr$layer_sampleA, '_', n=2)
agg_corr[, c("sampleB", "layerB")] <- str_split_fixed(agg_corr$layer_sampleB, '_', n=2)

meta_data <- comb_seurat@meta.data
meta_data <- meta_data[match(unique(meta_data$sample_name), meta_data$sample_name),]
meta_data <- meta_data[, c("sample_name", "diagnosis")]

# append one column without merge
# https://stackoverflow.com/a/17580007
agg_corr$diagnosisA <- meta_data$diagnosis[match(agg_corr$sampleA, meta_data$sample_name)]
agg_corr$diagnosisB <- meta_data$diagnosis[match(agg_corr$sampleB, meta_data$sample_name)]
# agg_corr <- merge(agg_corr, meta_data, by.x="sampleA", by.y="sample_name", all.x=TRUE)
# colnames(agg_corr)[colnames(agg_corr) == "diagnosis"] <- "diagnosisA"
# get dataframe

# between layers, between indivs
agg_corr[(agg_corr$sampleA != agg_corr$sampleB) & (agg_corr$layerA != agg_corr$layerB), "layer_x_indiv"] <- "btwn_layer_btwn_indiv"
# between layers, within indivs
agg_corr[(agg_corr$sampleA == agg_corr$sampleB) & (agg_corr$layerA != agg_corr$layerB), "layer_x_indiv"] <- "btwn_layer_wthn_indiv"
# within layers, between indivs
agg_corr[(agg_corr$sampleA != agg_corr$sampleB) & (agg_corr$layerA == agg_corr$layerB), "layer_x_indiv"] <- "wthn_layer_btwn_indiv"

# between layers, between diagnoses
agg_corr[(agg_corr$diagnosisA != agg_corr$diagnosisB) & (agg_corr$layerA != agg_corr$layerB), "layer_x_diagn"] <- "btwn_layer_btwn_diagn"
# between layers, within diagnoses
agg_corr[(agg_corr$diagnosisA == agg_corr$diagnosisB) & (agg_corr$layerA != agg_corr$layerB), "layer_x_diagn"] <- "btwn_layer_wthn_diagn"
# within layers, between diagnoses
agg_corr[(agg_corr$diagnosisA != agg_corr$diagnosisB) & (agg_corr$layerA == agg_corr$layerB), "layer_x_diagn"] <- "wthn_layer_btwn_diagn"

# add jitter to outliers
outliers_indiv <-
  agg_corr %>% group_by(layer_x_indiv) %>%
  filter(spearman_corr > quantile(spearman_corr, 0.75) + 1.5 * IQR(spearman_corr) |
         spearman_corr < quantile(spearman_corr, 0.25) - 1.5 * IQR(spearman_corr))
outliers_diagn <-
  agg_corr %>% group_by(layer_x_diagn) %>%
  filter(spearman_corr > quantile(spearman_corr, 0.75) + 1.5 * IQR(spearman_corr) |
         spearman_corr < quantile(spearman_corr, 0.25) - 1.5 * IQR(spearman_corr))


# create combined version with four categories:
# btwn_layer_btwn_indiv
# btwn_layer_wthn_diag
# btwn_layer_wthn_indiv
# wthn_layer_btwn_indiv
agg_corr_temp <- agg_corr[!is.na(agg_corr$layer_x_diagn) & agg_corr$layer_x_diagn=="btwn_layer_wthn_diagn", ]
agg_corr_temp$layer_x_indiv <- "btwn_layer_wthn_diagn"
agg_corr_comb <- rbind(agg_corr, agg_corr_temp)
agg_corr_comb$layer_x_diagn <- NULL
colnames(agg_corr_comb)[colnames(agg_corr_comb) == "layer_x_indiv"] <- "layer_x_grp"

outliers_comb <-
  agg_corr_comb %>% group_by(layer_x_grp) %>%
  filter(spearman_corr > quantile(spearman_corr, 0.75) + 1.5 * IQR(spearman_corr) |
         spearman_corr < quantile(spearman_corr, 0.25) - 1.5 * IQR(spearman_corr))

set.seed(101)
pdf("st_btwn_layer_btwn_indiv_06102024.pdf")
# ggplot(agg_corr, aes(x=layer_x_indiv, y=spearman_corr)) + geom_violin() + geom_boxplot() + theme_bw()
ggplot(agg_corr %>% filter(!is.na(layer_x_indiv)), aes(x=layer_x_indiv, y=spearman_corr)) + geom_violin() + geom_boxplot(width=0.1, outlier.shape=NA) + geom_beeswarm(data=outliers_indiv, corral="wrap", corral.width=0.5, cex=0.1) + theme_bw()

ggplot(agg_corr %>% filter(!is.na(layer_x_diagn)), aes(x=layer_x_diagn, y=spearman_corr)) + geom_violin() + geom_boxplot(width=0.1, outlier.shape=NA) + geom_beeswarm(data=outliers_diagn, corral="wrap", corral.width=0.5, cex=0.1) + theme_bw()

ggplot(agg_corr_comb %>% filter(!is.na(layer_x_grp)), aes(x=layer_x_grp, y=spearman_corr)) + geom_violin() + geom_boxplot(width=0.1, outlier.shape=NA) + geom_beeswarm(data=outliers_comb, corral="wrap", corral.width=0.5, cex=0.1) + theme_bw()
dev.off()

# within each sample, within each layer

# calculate correlations across cell type proportions instead of gene expression
# between layers between diagnoses
# between layers within diagnoses
# within layers between diagnoses



cell_prop_dir <- "/mnt/data0/projects/ASAP/ST/st_sc_deconvolution/results_V4_02202024/cell2location_map"
# cell_prop_dir <- "/mnt/data0/projects/ASAP/ST/st_sc_deconvolution/results_V4_cellsubtypes_02202024/cell2location_map"
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

comb_meta_cell2loc <- comb_meta_cell2loc[, c("sample_name", "diagnosis", "smoothed_label_s5", all_cell_type_cols)]
comb_meta_cell2loc[is.na(comb_meta_cell2loc)] <- 0
comb_meta_cell2loc <- comb_meta_cell2loc %>% group_by(sample_name, diagnosis, smoothed_label_s5) %>% summarise(across(everything(), mean)) %>% as.data.frame()

comb_meta_mat <- comb_meta_cell2loc
rownames(comb_meta_mat) <- paste(comb_meta_mat$sample_name, comb_meta_mat$diagnosis, comb_meta_mat$smoothed_label_s5, sep="_")
comb_meta_mat <- comb_meta_mat[, all_cell_type_cols]
c2loc_corr <- cor(t(comb_meta_mat), method="spearman") %>% as.data.frame()
c2loc_corr$layer_sampleA <- rownames(c2loc_corr)
c2loc_corr <- c2loc_corr %>% pivot_longer(!layer_sampleA, names_to="layer_sampleB", values_to="spearman_corr") %>% as.data.frame()

c2loc_corr[, c("sampleA", "diagnosisA", "layerA")] <- str_split_fixed(c2loc_corr$layer_sampleA, '_', n=3)
c2loc_corr[, c("sampleB", "diagnosisB", "layerB")] <- str_split_fixed(c2loc_corr$layer_sampleB, '_', n=3)

# between layers, between indivs
c2loc_corr[(c2loc_corr$sampleA != c2loc_corr$sampleB) & (c2loc_corr$layerA != c2loc_corr$layerB), "layer_x_indiv"] <- "btwn_layer_btwn_indiv"
# between layers, within indivs
c2loc_corr[(c2loc_corr$sampleA == c2loc_corr$sampleB) & (c2loc_corr$layerA != c2loc_corr$layerB), "layer_x_indiv"] <- "btwn_layer_wthn_indiv"
# within layers, between indivs
c2loc_corr[(c2loc_corr$sampleA != c2loc_corr$sampleB) & (c2loc_corr$layerA == c2loc_corr$layerB), "layer_x_indiv"] <- "wthn_layer_btwn_indiv"

# between layers, between diagnoses
c2loc_corr[(c2loc_corr$diagnosisA != c2loc_corr$diagnosisB) & (c2loc_corr$layerA != c2loc_corr$layerB), "layer_x_diagn"] <- "btwn_layer_btwn_diagn"
# between layers, within diagnoses
c2loc_corr[(c2loc_corr$diagnosisA == c2loc_corr$diagnosisB) & (c2loc_corr$layerA != c2loc_corr$layerB), "layer_x_diagn"] <- "btwn_layer_wthn_diagn"
# within layers, between diagnoses
c2loc_corr[(c2loc_corr$diagnosisA != c2loc_corr$diagnosisB) & (c2loc_corr$layerA == c2loc_corr$layerB), "layer_x_diagn"] <- "wthn_layer_btwn_diagn"

# add jitter to outliers
outliers_indiv <-
  c2loc_corr %>% group_by(layer_x_indiv) %>%
  filter(spearman_corr > quantile(spearman_corr, 0.75) + 1.5 * IQR(spearman_corr) |
         spearman_corr < quantile(spearman_corr, 0.25) - 1.5 * IQR(spearman_corr))
outliers_diagn <-
  c2loc_corr %>% group_by(layer_x_diagn) %>%
  filter(spearman_corr > quantile(spearman_corr, 0.75) + 1.5 * IQR(spearman_corr) |
         spearman_corr < quantile(spearman_corr, 0.25) - 1.5 * IQR(spearman_corr))

pdf("st_cellprops_btwn_layer_btwn_indiv_06102024.pdf")
ggplot(c2loc_corr %>% filter(!is.na(layer_x_indiv)), aes(x=layer_x_indiv, y=spearman_corr)) + geom_violin() + geom_boxplot(width=0.1, outlier.shape=NA) + geom_beeswarm(data=outliers_indiv, corral="wrap", corral.width=0.5, cex=0.1) + theme_bw()

ggplot(c2loc_corr %>% filter(!is.na(layer_x_diagn)), aes(x=layer_x_diagn, y=spearman_corr)) + geom_violin() + geom_boxplot(width=0.1, outlier.shape=NA) + geom_beeswarm(data=outliers_diagn, corral="wrap", corral.width=0.5, cex=0.1) + theme_bw()
dev.off()
