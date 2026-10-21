library(Seurat)
library(DESeq2)
library(stringr)
library(testit)

# layer-specific DE
# -- controls
# -- all samples
# -- smoothed or unsmoothed labels
# -- one layer vs all others

# layer-wise clinical DE
# -- nctx_temporal (principal phenotype)
# -- brain_stem_sn (early PD)
# -- last MMSE score (mental)
# -- motor UPDRS (motor)
# -- control/ILBD/case


# [1] "orig.ident"           "nCount_Spatial"       "nFeature_Spatial"
# [4] "sample_name"          "batch"                "PMI"
# [7] "RIN"                  "sex"                  "expired_age"
# [10] "diagnosis"            "last_mmse_test_score" "motor_updrs_score"
# [13] "nctx_temporal"        "brain_stem_sn"        "selected_spot"
# [16] "PMI_num"              "RIN_num"              "expired_age_num"
# [19] "Spatial_snn_res.0.3"  "seurat_clusters"      "Spatial_snn_res.0.35"
# [22] "Spatial_snn_res.0.4"  "Spatial_snn_res.0.45" "Spatial_snn_res.0.5"
# [25] "Spatial_snn_res.0.55" "Spatial_snn_res.0.6"  "layer_label_v1"
# [28] "layer_label_v2"       "smoothed_label_s1"    "s2_subclust"
# [31] "s2_subclust_score"    "s2_nbr_score"         "smoothed_label_s2"
# [34] "smoothed_label_s3"    "smoothed_label_s4"    "smoothed_label_s5"

comb_seurat <- readRDS("/mnt/data0/projects/ASAP/ST/clustering_pipeline/layer_smoothing_V4_seed101_11302023_replication_test/comb_seurat_post_smoothed_label_s5_11302023.rds")

random_seed <- 101

# Map diagnosis to pseudotime: Control = 0, ILBD = 1, Case = 2
comb_seurat@meta.data$diagnosis.int <- NA
comb_seurat@meta.data$diagnosis.int[comb_seurat@meta.data$diagnosis == "Control"] <- 0
comb_seurat@meta.data$diagnosis.int[comb_seurat@meta.data$diagnosis == "ILBD"] <- 1
comb_seurat@meta.data$diagnosis.int[comb_seurat@meta.data$diagnosis == "Case"] <- 2
comb_seurat@meta.data$diagnosis.int <- as.numeric(comb_seurat@meta.data$diagnosis.int)

# map diagnosis to ordered factor
comb_seurat@meta.data$diagnosis <- factor(comb_seurat@meta.data$diagnosis, ordered=TRUE, levels=c("Control", "ILBD", "Case"))

# in nctx_temporal, there is only 1 sample with a score of 4: set this to 3 so as not to overly bias the regression
comb_seurat@meta.data$nctx_temporal_fix <- comb_seurat@meta.data$nctx_temporal
comb_seurat@meta.data$nctx_temporal_fix[comb_seurat@meta.data$nctx_temporal_fix > 3] <- 3


run_layerwise <- function(input_seurat, layer_annot, min_layer_spot_count=20, controls_only=FALSE, random_seed=101) {
    
    if(controls_only) {
        comb_seurat <- subset(input_seurat, subset = diagnosis == "Control")
    } else {
        comb_seurat <- input_seurat
    }
    layers <- sort(unique(comb_seurat@meta.data[, layer_annot]))
    
    # loop over layers: create seurat object for each
    for(layer in layers) {
        layer_label <- sprintf("is_%s", gsub(' ', '', layer))
        comb_seurat@meta.data[, layer_label] <- comb_seurat@meta.data[, layer_annot] == layer
        Idents(object = comb_seurat) <- layer_label
        
        # get countData and filter out lowly expressed genes
        countData <- AggregateExpression(comb_seurat, assays="Spatial", slot="counts", group.by=c(layer_label, "sample_name"))
        countData <- countData$Spatial
        keep_genes <- rowSums(countData >= 5) >= 3
        countData <- countData[keep_genes, ]
        
        # filter out samples with fewer than N spots belonging to the given layer
        comb_seurat@meta.data$layer_label_sample_name <- paste(comb_seurat@meta.data[, layer_label], comb_seurat@meta.data$sample_name, sep='_')
        keep_samples <- names( which( table(comb_seurat@meta.data$layer_label_sample_name) > min_layer_spot_count ) )
        countData <- countData[, keep_samples]

        count_out_csv <- sprintf("agg_counts_filtered_layer_vs_layer_%s_%s_contsonly%s_V4_12052023.csv", layer_annot, gsub(' ', '', layer), controls_only)
        if(!file.exists(count_out_csv)) {
            write.csv(countData, file = count_out_csv)
        }
        
        sample_meta <- cbind(colnames(countData), str_split_fixed(colnames(countData), '_', n=2))
        sample_meta <- as.data.frame(sample_meta)
        rownames(sample_meta) <- sample_meta[, 1]
        sample_meta <- sample_meta[, c(2,3)]
        colnames(sample_meta) <- c(layer_label, "sample_name")
        sample_meta[, layer_label] <- relevel(as.factor(sample_meta[, layer_label]), ref="FALSE")
        

        # sample_meta <- comb_seurat@meta.data[, c("sample_name", "batch", "sex", "expired_age_num", "PMI_num", "RIN_num", layer_label)]
        # sample_meta <- sample_meta[!duplicated(sample_meta$sample_name), ]
        # sample_meta <- sample_meta[sample_meta$sample_name %in% colnames(countData), ]
        # filter out samples whose feature label is NA
        # sample_meta <- sample_meta[!is.na(sample_meta[, layer_label]), ]
        # rownames(sample_meta) <- sample_meta$sample_name
        # countData_filt <- countData[, as.character(sample_meta$sample_name)]
        # print(paste("Samples passing filtering", ncol(countData_filt)))
        
        assert(all(rownames(sample_meta) == colnames(countData)))
        
        sample_meta$sample_name <- as.factor(sample_meta$sample_name)
        
        assert(is.factor(sample_meta$sample_name))
        assert(is.factor(sample_meta[, layer_label]))
        
        tryCatch(
            expr = {
                set.seed(random_seed)
                dds <- DESeqDataSetFromMatrix(countData=countData, colData=sample_meta, design= model.matrix(~ sample_meta$sample_name + sample_meta[, layer_label]), tidy=FALSE)
                saveRDS(dds, sprintf("pseudobulk_layer_deseq2_dds_%s_contsonly%s_12052023.rds", layer_label, controls_only))

                dds_lrt <- DESeq(dds, test="LRT", reduced = model.matrix(~ sample_meta$sample_name))
                saveRDS(dds_lrt, sprintf("pseudobulk_layer_deseq2_dds_lrt_%s_contsonly%s_12052023.rds", layer_label, controls_only))

                results_table <- results(dds_lrt)

                write.csv(results_table, sprintf("pseudobulk_layer_deseq2_%s_contsonly%s_12052023.csv", layer_label, controls_only))

            },
            error = function(e){
                # (Optional)
                # Do this if an error is caught...
                print(paste("Running", layer_label, "controls_only:", controls_only))
                message(e)
            }
        )

        
    }
}

# decide which layer annotation to use: smoothed or unsmoothed
layer_annot <- "smoothed_label_s5"
layers <- sort(unique(comb_seurat@meta.data[, layer_annot]))
# features_quantitative <- c("nctx_temporal_fix", "brain_stem_sn", "last_mmse_test_score", "motor_updrs_score", "diagnosis", "diagnosis.int")
min_layer_spot_count <- 20

run_layerwise(comb_seurat, layer_annot, min_layer_spot_count=20, controls_only=FALSE, random_seed=random_seed)
run_layerwise(comb_seurat, layer_annot, min_layer_spot_count=20, controls_only=TRUE, random_seed=random_seed)
