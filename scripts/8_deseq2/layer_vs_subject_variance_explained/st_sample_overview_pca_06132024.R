library(Seurat)
library(testit)
library(stringr)
library(dplyr)
library(tibble)

library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(corrplot)
library(umap)

# comb_sc_all <- readRDS("/data/neurogen/Jacob/PD5D_Full_Integration_Files/FullIntegrationOSR_MTG_DietSeuratFinal.rds")

# panda server
# comb_seurat <- readRDS("/mnt/data0/projects/ASAP/ST/data/comb_seurat_smoothed_final_01102023.rds")
comb_seurat <- readRDS("/mnt/data0/projects/ASAP/ST/clustering_pipeline/layer_smoothing_V4_seed101_11302023_replication_test/comb_seurat_post_smoothed_label_s5_11302023.rds")
# remove the "_rep1" from some samples
# comb_seurat@meta.data$sample_name <- gsub("_rep[0-9]*", "", comb_seurat@meta.data$orig.ident)
comb_seurat@meta.data$spot_id <- rownames(comb_seurat@meta.data)



# aggregate expression and lognormalize it
comb_agg <- AggregateExpression(comb_seurat, group.by="sample_name", return.seurat=TRUE)
comb_agg <- FindVariableFeatures(comb_agg)
comb_agg <- ScaleData(comb_agg)
comb_agg <- RunPCA(comb_agg, npcs=30)
comb_agg <- RunUMAP(comb_agg, dims=1:30)
comb_agg <- FindNeighbors(comb_agg, reduction="pca", dims=1:30)
comb_agg <- FindClusters(comb_agg, resolution=0.5)



meta_fields <- comb_seurat@meta.data
meta_fields <- meta_fields[match(unique(meta_fields$sample_name), meta_fields$sample_name),]
meta_fields <- meta_fields[, c("sample_name", "batch", "PMI_num", "RIN_num", "sex", "expired_age_num", "diagnosis", "last_mmse_test_score", "motor_updrs_score", "nctx_temporal")]
comb_agg@meta.data$sample_name <- rownames(comb_agg@meta.data)
meta_fields <- merge(comb_agg@meta.data, meta_fields, by="sample_name", sort=FALSE)
rownames(meta_fields) <- meta_fields$sample_name
# comb_agg@meta.data <- meta_fields

umap_embed <- comb_agg@reductions$umap@cell.embeddings %>% as.data.frame()
umap_embed$sample_name <- rownames(umap_embed)
meta_fields <- merge(meta_fields, umap_embed, by="sample_name", sort=FALSE)
pca_embed <- comb_agg@reductions$pca@cell.embeddings %>% as.data.frame()
pca_embed$sample_name <- rownames(pca_embed)
meta_fields <- merge(meta_fields, pca_embed[, c("sample_name", "PC_1", "PC_2")], by="sample_name", sort=FALSE)
meta_fields$nCount_Spatial <- meta_fields$nCount_Spatial / 10000 # is this being multiplied by the scale.factor in AggregateExpression()?

pdf("st_subject_level_umap_06132024.pdf")
ggplot(meta_fields, aes(x = PC_1, y = PC_2, colour = nCount_Spatial)) +
  geom_point() + theme_bw()
ggplot(meta_fields, aes(x = PC_1, y = PC_2, colour = batch)) +
  geom_point() + theme_bw()
ggplot(meta_fields, aes(x = PC_1, y = PC_2, colour = PMI_num)) +
  geom_point() + theme_bw()
ggplot(meta_fields, aes(x = PC_1, y = PC_2, colour = RIN_num)) +
  geom_point() + theme_bw()
ggplot(meta_fields, aes(x = PC_1, y = PC_2, colour = sex)) +
  geom_point() + theme_bw()
ggplot(meta_fields, aes(x = PC_1, y = PC_2, colour = expired_age_num)) +
  geom_point() + theme_bw()
ggplot(meta_fields, aes(x = PC_1, y = PC_2, colour = diagnosis)) +
  geom_point() + theme_bw()
ggplot(meta_fields, aes(x = PC_1, y = PC_2, colour = last_mmse_test_score)) +
  geom_point() + theme_bw()
ggplot(meta_fields, aes(x = PC_1, y = PC_2, colour = motor_updrs_score)) +
  geom_point() + theme_bw()
ggplot(meta_fields, aes(x = PC_1, y = PC_2, colour = nctx_temporal)) +
  geom_point() + theme_bw()
dev.off()

# instead of variablefeatures, use the set of DE genes across all clinical features
# (run this locally, then export to server)
source("/Users/jy1008/Downloads/pseudobulk_DE_DESeq2_GSEA_12292023_utils.R")
layers <- c("Layer 1", "Layer 2", "Layer 3", "Layer 4", "Layer 5", "Layer 6", "WM")
effects_path <- "/Users/jy1008/Downloads/pseudobulk_DE_DESeq2_V4_12052023"
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
write.csv(combined_res, "st_de_combined_res_06132024.csv")


combined_res <- read.csv("~/st_de_combined_res_06132024.csv")
combined_res <- combined_res[combined_res$padj < 0.1, ]
de_gene_list <- unique(combined_res$gene)

VariableFeatures(comb_agg) <- de_gene_list
comb_agg <- ScaleData(comb_agg)
comb_agg <- RunPCA(comb_agg, npcs=30)
comb_agg <- RunUMAP(comb_agg, dims=1:30)
comb_agg <- FindNeighbors(comb_agg, reduction="pca", dims=1:30)
comb_agg <- FindClusters(comb_agg, resolution=0.5)


meta_fields <- comb_seurat@meta.data
meta_fields <- meta_fields[match(unique(meta_fields$sample_name), meta_fields$sample_name),]
meta_fields <- meta_fields[, c("sample_name", "batch", "PMI_num", "RIN_num", "sex", "expired_age_num", "diagnosis", "last_mmse_test_score", "motor_updrs_score", "nctx_temporal")]
comb_agg@meta.data$sample_name <- rownames(comb_agg@meta.data)
meta_fields <- merge(comb_agg@meta.data, meta_fields, by="sample_name", sort=FALSE)
rownames(meta_fields) <- meta_fields$sample_name
# comb_agg@meta.data <- meta_fields

umap_embed <- comb_agg@reductions$umap@cell.embeddings %>% as.data.frame()
umap_embed$sample_name <- rownames(umap_embed)
meta_fields <- merge(meta_fields, umap_embed, by="sample_name", sort=FALSE)
pca_embed <- comb_agg@reductions$pca@cell.embeddings %>% as.data.frame()
pca_embed$sample_name <- rownames(pca_embed)
meta_fields <- merge(meta_fields, pca_embed[, c("sample_name", "PC_1", "PC_2")], by="sample_name", sort=FALSE)
meta_fields$nCount_Spatial <- meta_fields$nCount_Spatial / 10000 # is this being multiplied by the scale.factor in AggregateExpression()?

pdf("st_subject_level_umap_de_06132024.pdf")
ggplot(meta_fields, aes(x = PC_1, y = PC_2, colour = nCount_Spatial)) +
  geom_point() + theme_bw()
ggplot(meta_fields, aes(x = PC_1, y = PC_2, colour = batch)) +
  geom_point() + theme_bw()
ggplot(meta_fields, aes(x = PC_1, y = PC_2, colour = PMI_num)) +
  geom_point() + theme_bw()
ggplot(meta_fields, aes(x = PC_1, y = PC_2, colour = RIN_num)) +
  geom_point() + theme_bw()
ggplot(meta_fields, aes(x = PC_1, y = PC_2, colour = sex)) +
  geom_point() + theme_bw()
ggplot(meta_fields, aes(x = PC_1, y = PC_2, colour = expired_age_num)) +
  geom_point() + theme_bw()
ggplot(meta_fields, aes(x = PC_1, y = PC_2, colour = diagnosis)) +
  geom_point() + theme_bw()
ggplot(meta_fields, aes(x = PC_1, y = PC_2, colour = last_mmse_test_score)) +
  geom_point() + theme_bw()
ggplot(meta_fields, aes(x = PC_1, y = PC_2, colour = motor_updrs_score)) +
  geom_point() + theme_bw()
ggplot(meta_fields, aes(x = PC_1, y = PC_2, colour = nctx_temporal)) +
  geom_point() + theme_bw()
dev.off()

