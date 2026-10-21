library(Seurat)


comb_sc_all <- readRDS("/data/neurogen/Jacob/PD5D_Full_Integration_Files/FullIntegrationOSR_MTG_DietSeuratFinal.rds")

out_dir <- "/data/neurogen/ASAP/ST/data/sc_for_deconvolution"

sample_list <- unique(comb_sc_all@meta.data$sample_id)
for(sample in sample_list) {
    comb_sub <- subset(comb_sc_all, subset = sample_id == sample)
    
    write.table(comb_sub@meta.data,
                file.path(out_dir, sprintf("sc_count_%s_meta_02062023.csv", sample)),
                sep = ',', row.names = T, col.names = NA, quote = F)
    write.table(as.matrix(GetAssayData(object = comb_sub, slot = "counts")),
                file.path(out_dir, sprintf("sc_count_%s_02062023.csv", sample)),
                sep = ',', row.names = T, col.names = NA, quote = F)
}


# Convert spatial data from Seurat to AnnData
# https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html
library(SeuratDisk)
library(Seurat)
library(dplyr)

comb_seurat <- readRDS("/data/neurogen/jy1008/scRNAseq_integration/layer_smoothing_06032022/comb_seurat_smoothed_final_01102023.rds")
SaveH5Seurat(comb_seurat, filename = "comb_seurat_smoothed_final_01102023.h5Seurat")
Convert("comb_seurat_smoothed_final_01102023.h5Seurat", dest = "comb_seurat_smoothed_final_01102023.h5ad")

comb_seurat <- readRDS("/data/neurogen/jy1008/scRNAseq_integration/temp_sandbox_02052024/layer_smoothing_V4_seed101_11302023_replication_test/comb_seurat_post_smoothed_label_s5_11302023.rds")
comb_seurat_diet <- DietSeurat(comb_seurat)
comb_seurat_diet@meta.data <- comb_seurat_diet@meta.data %>% mutate(across(everything(), as.character))
SaveH5Seurat(comb_seurat_diet, filename = "comb_seurat_post_smoothed_label_s5_diet_11302023.h5Seurat")
Convert("comb_seurat_post_smoothed_label_s5_diet_11302023.h5Seurat", dest = "comb_seurat_post_smoothed_label_s5_diet_11302023.h5ad")
