library(SingleCellExperiment)
library(Seurat)
library(testit)

library(data.table)
library(ggplot2)
library(gridExtra)

library(reshape2)
library(tidyr)
library(dplyr)
library(e1071)

maynard_dir <- "/mnt/data/externalData/ST/Maynard2021/raw_files_to_Visium"
sample_list <- c("151507", "151508", "151509", "151510", "151669", "151670", "151671", "151672", "151673", "151674", "151675", "151676")

maynard_sce <- readRDS("/mnt/data/externalData/ST/Maynard2021/Maynard_All_Samples_sce.rds")
head(colData(maynard_sce))
maynard_meta <- colData(maynard_sce)

seurat_objs <- c()
for(sample_id in sample_list) {
    samp_dir <- file.path(maynard_dir, sample_id)
    seurat_obj <- Load10X_Spatial(
                    data.dir = samp_dir,
                    filename = "filtered_feature_bc_matrix.h5",
                    assay = "Spatial",
                    slice = paste0('M', sample_id),
                    filter.matrix = TRUE,
                    to.upper = FALSE)
    seurat_objs <- c(seurat_objs, seurat_obj)
}
maynard_comb <- merge(seurat_objs[[1]], y = seurat_objs[2:length(seurat_objs)], add.cell.ids = sample_list)
maynard_comb@meta.data$Row.names <- rownames(maynard_comb@meta.data)
maynard_comb@meta.data$Row.order <- 1:nrow(maynard_comb@meta.data)
df_merged <- merge(maynard_comb@meta.data, maynard_meta, by.x="Row.names", by.y="key", all.x=TRUE)
df_merged <- df_merged[order(df_merged$Row.order), ]
rownames(df_merged) <- df_merged$Row.names
df_merged$Row.names <- NULL
df_merged$Row.order <- NULL
assert(all(rownames(df_merged) == colnames(maynard_comb)))
maynard_comb@meta.data <- as.data.frame(df_merged)

pdf("temp_test_Maynard2021_spatial_01132025.pdf")
for (img in names(maynard_comb@images)) {
    sp_pt <- SpatialPlot(maynard_comb, group.by="layer_guess", image.alpha=0.3, images=img, stroke=NA)
    print(sp_pt)
}
dev.off()



# plot some marker genes
maynard_comb <- NormalizeData(maynard_comb,
                              normalization.method = "LogNormalize",
                              scale.factor = 10000)
marker_gene <- "LAMP5"
pdf(sprintf("Maynard2021_marker_%s_01152025.pdf", marker_gene))
for(img in names(maynard_comb@images)) {
    print(img)
    sp_pt <- SpatialFeaturePlot(maynard_comb, features=marker_gene, image.alpha=0.3, images=img, stroke=NA, slot="data")
    print(sp_pt)
}
dev.off()




# from main clustering pipeline
# spatial_pipeline_filter_varfeat_duplicates_morecovars_Maynard_01132025.R

library(Seurat)
library(reticulate)
use_condaenv("r_leidenalg", required=T)
library(hdf5r)
library(SeuratWrappers)
library(harmony)
library(data.table)
library(glmpca)
library(ggplot2)

# remove duplicate columns from a sparse matrix
# source: https://stackoverflow.com/a/51457395
duplicated.dgCMatrix <- function (dgCMat, MARGIN, include.all.zero.vectors = TRUE) {
    MARGIN <- as.integer(MARGIN)
    J <- rep(1:ncol(dgCMat), diff(dgCMat@p))
    I <- dgCMat@i + 1
    x <- dgCMat@x
    if (MARGIN == 1L) {
        ## check duplicated rows
        names(x) <- J
        if (include.all.zero.vectors) {
            RowLst <- split(x, factor(I, levels = 1:nrow(dgCMat)))
        } else {
            RowLst <- split(x, I)  ## will do `factor(I)` internally in `split`
        }
        result <- duplicated.default(RowLst)
    } else if (MARGIN == 2L) {
        ## check duplicated columns
        names(x) <- I
        if (include.all.zero.vectors) {
            ColLst <- split(x, factor(J, levels = 1:ncol(dgCMat)))
        } else {
            ColLst <- split(x, J)  ## will do `factor(J)` internally in `split`
        }
        result <- duplicated.default(ColLst)
    } else {
        warning("invalid MARGIN; return NULL")
        result <- NULL
    }
    result
}


# remove variablefeatures if they appear too rarely in the spatial samples
# pass_genes <- filter_variablefeatures(combined_controls, VariableFeatures(combined_controls), nzero_spot_thresh=20, frac_sample_thresh=0.3, bin_size=50)
filter_variablefeatures <- function(comb_seurat, gene_list, nzero_spot_thresh=20,
                               frac_sample_thresh=0.3, bin_size=50) {
    num_samples <- length(unique(comb_seurat@meta.data$orig.ident))
    sample_thresh <- floor(num_samples * frac_sample_thresh)
    print(paste("Genes must be nonzero in >", nzero_spot_thresh,
                " spots in >", sample_thresh, " samples", sep=""))
    sample_ids <- comb_seurat@meta.data$orig.ident

    pass_genes <- c()
    ii <- 0
    
    for (gene_chunk in split(gene_list, ceiling(seq_along(gene_list)/bin_size))) {
        gene_matrix <- comb_seurat@assays$Spatial@counts[gene_chunk, ]
        gene_matrix_gt0 <- gene_matrix > 0
        for (row in 1:nrow(gene_matrix_gt0)) {
            gene <- rownames(gene_matrix_gt0)[row]
            idxs <- gene_matrix_gt0[row, ] # slow step
            nzero_spots <- sample_ids[idxs]
            nzero_spots <- table(nzero_spots)
            pass_samples <- names(nzero_spots)[nzero_spots > nzero_spot_thresh]
            ii <- ii+1
            if (length(pass_samples) > sample_thresh) {
                pass_genes <- c(pass_genes, gene)
                print(paste(length(pass_genes), "/", ii, "Passed:", gene, sep=" "))
            } else {
                print(paste(length(pass_genes), "/", ii, "Failed:", gene, sep=" "))
            }
        }
    }
    return(pass_genes)
}


cluster_pipeline <- function(comb_seurat, job_name,
                        harmony_vars = c("orig.ident", "batch", "sex", "expired_age"),
                        harmony_thetas = c(0.4, 0.4, 0.4, 0.4),
                        resolutions_list = c(0.5, 0.8),
                        selected_de_genes = NULL,
                        num_pcs = 30,
                        random_seed = 101,
                        remove_mito_ribo = FALSE) {
    dir.create(job_name)

    checkpoints_file <- file.path(job_name, "checkpoints.txt")
    write(paste(Sys.time(), job_name),file=checkpoints_file,append=FALSE)
    write(paste(Sys.time(), "Harmony covars:"),file=checkpoints_file,append=TRUE)
    write(harmony_vars,file=checkpoints_file,append=TRUE)
    write(harmony_thetas,file=checkpoints_file,append=TRUE)
    write(paste(Sys.time(), "random seed:", random_seed),file=checkpoints_file,append=TRUE)
    write(paste(Sys.time(), "cluster res:"),file=checkpoints_file,append=TRUE)
    write(resolutions_list,file=checkpoints_file,append=TRUE)
    
    comb_seurat <- NormalizeData(comb_seurat,
                       normalization.method = "LogNormalize",
                       scale.factor = 10000)
    comb_seurat <- FindVariableFeatures(comb_seurat,
                       selection.method = "vst",
                       nfeatures = 4000)
    comb_seurat <- ScaleData(comb_seurat,
                       features = rownames(comb_seurat))

    gene_list <- VariableFeatures(comb_seurat)
    pass_genes <- filter_variablefeatures(comb_seurat, gene_list,
                                nzero_spot_thresh=20,
                                frac_sample_thresh=0.3, bin_size=50)
    # almost no intersection here?
    # if (remove_mito_ribo) {
    #     mito_genes <- rownames(comb_seurat)[grep("^MT-",rownames(comb_seurat))]
    #     ribo_genes <- rownames(comb_seurat)[grep("^RP[SL]",rownames(comb_seurat))]
    # }
    VariableFeatures(comb_seurat) <- pass_genes

    # New step: only grab genes which are also DE genes from reference data
    if(!is.null(selected_de_genes)) {
        print(paste("Num variable genes:", length(VariableFeatures(comb_seurat))))
        VariableFeatures(comb_seurat) <- intersect(VariableFeatures(comb_seurat), unique(selected_de_genes$gene))
        print(paste("Num variable genes:", length(VariableFeatures(comb_seurat))))
    }

    # filter out spots that are all-zero in variable genes, or duplicates
    varfeat_counts <- comb_seurat@assays$Spatial@counts[VariableFeatures(comb_seurat), ]
    dup_spots <- which(colSums(varfeat_counts) == 0 | duplicated.dgCMatrix(varfeat_counts, 2))
    total_zero_count <- length(which(colSums(varfeat_counts) == 0))
    print(sprintf("Variable features passing filtering: %s", length(pass_genes)))
    write(sprintf("Variable features passing filtering: %s", length(pass_genes)),file=checkpoints_file,append=TRUE)
    print(sprintf("Total spots all zero in filtered variable features: %s", total_zero_count))
    write(sprintf("Total spots all zero in filtered variable features: %s", total_zero_count),file=checkpoints_file,append=TRUE)
    print(sprintf("Total dup spots removed: %s", length(dup_spots)))
    write(sprintf("Total dup spots removed: %s", length(dup_spots)),file=checkpoints_file,append=TRUE)

    # NOTE: if dup_spots is empty, this step fails
    if(length(dup_spots) > 0) {
        comb_seurat <- comb_seurat[, -dup_spots]
    }

    # To save memory, remove all genes which are not in VariableFeatures
    # comb_seurat <- subset(comb_seurat, features = VariableFeatures(comb_seurat))
    # gc(verbose = TRUE)

    # comb_seurat <- RunPCA(comb_seurat, pc.genes = comb_seurat@var.genes, npcs = 30)
    # eventually, you should switch over, but for now, confirm PCA still works
    set.seed(random_seed)
    print("Running GLMPCA")
    # comb_seurat <- RunGLMPCA(comb_seurat,
    #                    features = pass_genes, L = 30,
    #                    minibatch = "stochastic", batch_size = 10000)

    # select number of PCs
    # https://stats.stackexchange.com/q/254592
    # https://github.com/willtownes/glmpca/issues/32#issuecomment-859748250
    # eigs <- comb_seurat@reductions$glmpca@stdev^2
    
    comb_seurat <- RunPCA(comb_seurat, npcs=num_pcs)
    eigs <- comb_seurat@reductions$pca@stdev^2
    prop_variance <- eigs / sum(eigs)
    PC <- 1:length(prop_variance)

    # pdf(file.path(job_name, "glmpca_out.pdf"))
    pdf(file.path(job_name, "pca_out.pdf"))
    ve_pt <- ggplot(data.frame(PC, prop_variance), aes(x=PC, y=prop_variance)) + geom_line()
    print(ve_pt)
    # dm_pt <- DimPlot(object = comb_seurat, reduction = "glmpca",
    #              pt.size = .1, group.by = "orig.ident") + NoLegend()
    dm_pt <- DimPlot(object = comb_seurat, reduction = "pca",
                 pt.size = .1, group.by = "orig.ident") + NoLegend()
    print(dm_pt)
    # vn_pt <- VlnPlot(object = comb_seurat, features = "GLMPC_1",
    #              group.by = "orig.ident", pt.size = .1) + NoLegend()
    vn_pt <- VlnPlot(object = comb_seurat, features = "PC_1",
                 group.by = "orig.ident", pt.size = .1) + NoLegend()
    print(vn_pt)
    dev.off()

    # combined <- RunHarmony(combined, group.by.vars = c("orig.ident"), reduction = "pca", plot_convergence = TRUE, theta = c(1))
    set.seed(random_seed)
    print("Running Harmony")
    pdf(file.path(job_name, "harmony_out.pdf"))
    # comb_seurat <- RunHarmony(comb_seurat, reduction.use = "glmpca",
    comb_seurat <- RunHarmony(comb_seurat, reduction.use = "pca",
                       plot_convergence = TRUE,
                       max_iter = 20,
                       assay.use = "Spatial",
                       group.by.vars = harmony_vars,
                       theta = harmony_thetas,
                       lambda = rep(1, length(harmony_thetas)))
    dm_pt <- DimPlot(object = comb_seurat, reduction = "harmony",
                 pt.size = .1, group.by = "orig.ident") + NoLegend()
    print(dm_pt)
    vn_pt <- VlnPlot(object = comb_seurat, features = "harmony_1",
                 group.by = "orig.ident", pt.size = .1) + NoLegend()
    print(vn_pt)
    dev.off()

    # The last value for resolution is what is stored in seurat_clusters
    # Use Spatial_snn_res.X for clustering
    comb_seurat <- comb_seurat %>%
        RunUMAP(reduction = "harmony", dims = 1:num_pcs) %>%
        FindNeighbors(reduction = "harmony", dims = 1:num_pcs) %>%
        identity()
    for (res in resolutions_list) {
        set.seed(random_seed)
        print(paste("Finding Clusters with res=", res, sep=''))
        column_name <- paste("Spatial_snn_res.", res, sep='')
        # comb_seurat <- comb_seurat %>%
        #     RunUMAP(reduction = "harmony", dims = 1:30) %>%
        #     FindNeighbors(reduction = "harmony", dims = 1:30) %>%
        #     FindClusters(resolution = res, algorithm=4, method="igraph") %>%
        #     identity()

        # comb_seurat <- comb_seurat %>%
        #     FindClusters(resolution = res, algorithm=4, method="igraph") %>%
        #     identity()
        # Use Louvain for Maynard
        comb_seurat <- comb_seurat %>%
            FindClusters(resolution = res, algorithm=1, method="igraph") %>%
            identity()

        # pdf(file.path(job_name,
        #               paste("leiden_out_res", res, ".pdf", sep='')))
        pdf(file.path(job_name,
                      paste("louvain_out_res", res, ".pdf", sep='')))
        dm_pt <- DimPlot(comb_seurat, reduction = "umap",
                     label = TRUE, pt.size = .1,
                     group.by = column_name) + ggtitle(column_name)
        print(dm_pt)
        dev.off()

        pdf(file.path(job_name,
                      paste("spatial_mapped_res", res, ".pdf", sep='')))
        for (img in names(comb_seurat@images)) {
            sp_pt <- SpatialPlot(comb_seurat, group.by= column_name, image.alpha=0.3, images=img)
            print(sp_pt)
        }
        dev.off()
    }

    # Save the Seurat object
    # Regenerate data and scale.data upon reloading
    comb_diet <- DietSeurat(comb_seurat, counts=TRUE, data=TRUE,
                     scale.data=FALSE, dimreducs = c('glmpca', 'harmony', 'umap'))
    saveRDS(comb_diet, file.path(job_name, "comb_seurat_diet_out.rds"))
    return(comb_seurat)
}


# get list of DE genes and use them as variable genes
de_gene_path <- "/mnt/data/projects/ASAP/ST/pseudobulk_DE_DESeq2_V4_12052023_backup/layer_vs_layer"
st_layers <- c(paste0("Layer", 1:6), "WM")
de_genes <- c()
for(st_layer in st_layers) {
    st_file <- file.path(de_gene_path, sprintf("pseudobulk_layer_deseq2_is_%s_contsonlyTRUE_12052023.csv", st_layer))
    de_genes_layer <- read.csv(st_file)
    de_genes_layer$layer <- st_layer
    colnames(de_genes_layer)[colnames(de_genes_layer) == 'X'] <- "gene"
    de_genes <- rbind(de_genes, de_genes_layer)
}
de_genes[is.na(de_genes$padj), "padj"] <- 1
# get the N most-significant marker genes for each layer
num_de_sig <- 100
# num_de_sig <- 200
de_genes_filt <- de_genes %>%
  group_by(layer) %>%
  slice_min(padj, n = num_de_sig) %>%
  ungroup()




setwd("/home/jy743/maynard_smoothing_test_01142025/clustering_pipeline_run")
comb_seurat <- maynard_comb

# merged_seurat <- "/mnt/data0/projects/ASAP/ST/data/MTG_visium/comb_seurat_V4_postqc_10272023.rds"
# combined <- readRDS(merged_seurat)
comb_seurat@meta.data$orig.ident <- comb_seurat@meta.data$subject

resolutions_list <- c(0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6)
# harmony_vars <- c("orig.ident", "sex", "expired_age", "PMI", "RIN")
# harmony_thetas <- c(0.4, 0.4, 0.4, 0.4, 0.4)

# all samples, with variablefeatures intersecting marker genes
harmony_vars <- c("sample_name")
# harmony_vars <- c("orig.ident")
harmony_thetas <- c(2.0)
for(random_seed in c(101)) {
    # combined_looprun <- combined
    combined_looprun <- comb_seurat
    job_name <- sprintf("all_data_1covar2p0_V4_nfsig%s_seed%s_01152025", num_de_sig, random_seed)
    cluster_pipeline(combined_looprun, job_name,
                     harmony_vars = harmony_vars,
                     harmony_thetas = harmony_thetas,
                     resolutions_list = resolutions_list,
                     selected_de_genes = de_genes_filt,
                     random_seed = random_seed)
}

# load and subcluster
# 0 and 1 -> layers 2, 3, and 4
comb_seurat <- readRDS("/home/jy743/maynard_smoothing_test_01142025/clustering_pipeline_run/all_data_1covar2p0_V4_nfsig100_seed101_01152025/comb_seurat_diet_out.rds")

comb_sub_L234 <- subset(comb_seurat, Spatial_snn_res.0.4 %in% c(0,1))

num_de_sig <- 100
de_genes_filt_234 <- de_genes %>%
  group_by(layer) %>%
  filter(layer %in% c("Layer2", "Layer3", "Layer4")) %>%
  slice_min(padj, n = num_de_sig) %>%
  ungroup()

# all samples, with variablefeatures intersecting marker genes
harmony_vars <- c("sample_name")
harmony_thetas <- c(2.0)
random_seed <- 101
job_name <- sprintf("all_data_1covar2p0_V4_subL234_nfsig%s_seed%s_01152025", num_de_sig, random_seed)
comb_sub_L234 <- cluster_pipeline(comb_sub_L234, job_name,
                 harmony_vars = harmony_vars,
                 harmony_thetas = harmony_thetas,
                 resolutions_list = resolutions_list,
                 selected_de_genes = de_genes_filt_234,
                 num_pcs = 10,
                 random_seed = random_seed)
comb_sub_L234 <- readRDS("/home/jy743/maynard_smoothing_test_01142025/clustering_pipeline_run/all_data_1covar2p0_V4_subL234_nfsig100_seed101_01152025/comb_seurat_diet_out.rds")

comb_seurat@meta.data$Row.names <- rownames(comb_seurat@meta.data)
comb_seurat@meta.data$row_index <- 1:nrow(comb_seurat@meta.data)
comb_sub_L234@meta.data$Row.names <- rownames(comb_sub_L234@meta.data)
comb_sub_L234@meta.data$subclust <- as.character(comb_sub_L234@meta.data$Spatial_snn_res.0.4)
# combine clusters 0 and 1 into a single cluster
comb_sub_L234@meta.data[comb_sub_L234@meta.data$subclust == "0", "subclust"] <- "1"
temp_meta <- merge(comb_seurat@meta.data, comb_sub_L234@meta.data[, c("Row.names", "subclust")], by = "Row.names", all.x=TRUE)
temp_meta <- temp_meta[order(temp_meta$row_index), ]
temp_meta$row_index <- NULL
rownames(temp_meta) <- temp_meta$Row.names
temp_meta$Row.names <- NULL
temp_meta$maj_clust <- as.character(temp_meta$Spatial_snn_res.0.4)
temp_meta$combined_clust <- temp_meta$maj_clust
temp_meta[!is.na(temp_meta$subclust), "combined_clust"] <- paste0("sub_", temp_meta[!is.na(temp_meta$subclust), "subclust"])
assert(all(colnames(comb_seurat) == rownames(temp_meta)))
comb_seurat@meta.data <- temp_meta
comb_seurat@meta.data$combined_clust <- factor(comb_seurat@meta.data$combined_clust)

setwd("/home/jy743/maynard_smoothing_test_01142025/clustering_pipeline_run/all_data_1covar2p0_V4_subL234_nfsig100_seed101_01152025")
pdf("temp_test_Maynard2021_spatial_subclust_comb_01152025.pdf")
for (img in names(comb_seurat@images)) {
    sp_pt <- SpatialPlot(comb_seurat, group.by="combined_clust", image.alpha=0.3, images=img, stroke=NA)
    print(sp_pt)
}
dev.off()
saveRDS(comb_seurat, "comb_seurat_subclust_01152025.rds")



# from gauss_kernel_convolution_cluster_smoothing_V4_maynard_01132025.R

group.colors <- c("Layer 1"="#f43545", "Layer 2"="#ff8901", "Layer 3"="#fad717",
                  "Layer 4"="#00ba71", "Layer 5"="#00c2de", "Layer 6"="#ba55d3",
                  "WM"="#652175")
group.colors.unsmooth <- c("Layer 1"="#f43545", "Layer 2"="#ff8901", "Layer 3"="#fad717",
                  "Layer 4"="#00ba71", "Layer 5"="#00c2de", "Layer 6"="#ba55d3",
                  "WM_1"="#3b2175", "WM_2"="#75215b", "WM_3"="#652175", "Noise"="#c8c8c8")

random_seed <- 101

svm_smooth_labels <- function(comb_seurat, cluster_label_name, svm_gamma = 3, svm_cost = 50, class_weights = NULL, image_name = NULL, output_label_name="smoothed_label_out", random_seed = 101) {
    
    # cluster_label_name needs to have levels
    assert(is.factor(comb_seurat@meta.data[, cluster_label_name]))
    
    cluster_labels <- comb_seurat@meta.data[, c("orig.ident", cluster_label_name)]
    comb_seurat@meta.data$smoothed_label <- NA

    if(is.null(image_name)) {
        img_list <- names(comb_seurat@images)
    } else {
        img_list <- c(image_name)
    }
    for (img in img_list) {
        # print(img)
        set.seed(random_seed)

        coord_table <- comb_seurat@images[[img]]@coordinates
        coord_table <- merge(coord_table, cluster_labels, by=0)
        rownames(coord_table) <- coord_table$Row.names

        # hexagonal
        coord_table$euclid_x <- (coord_table$row %% 2) / 2 + (coord_table$col - (coord_table$col %% 2)) / 2
        coord_table$euclid_y <- sqrt(3) * ( (coord_table$row %% 2) / 2 + (coord_table$row - (coord_table$row %% 2)) / 2)

        coord_table$smoothed_label <- as.factor(coord_table[, cluster_label_name])
        if(is.null(svm_gamma) | is.null(svm_cost)) {
            model <- svm(x = coord_table[,c("euclid_x","euclid_y")],
                         y = coord_table$smoothed_label)
        } else {
            model <- svm(x = coord_table[,c("euclid_x","euclid_y")],
                         y = coord_table$smoothed_label,
                         gamma = svm_gamma,
                         cost = svm_cost,
                         class.weights = class_weights)
                         # for L4
                         # class.weights = c("Layer 2" = 1, "Layer 3" = 1, "Layer 4" = 2,
                         #                   "Layer 5" = 1, "Layer 6" = 1))
                         # class.weights = c("TRUE" = svm_true_weight, "FALSE" = 1))
        }
        coord_table$svm_pred <- predict(model, x=coord_table[,c("row","col")])

        # comb_seurat@meta.data[rownames(coord_table), "smoothed_label"] <- coord_table$smoothed_label
        comb_seurat@meta.data[rownames(coord_table), output_label_name] <- coord_table$svm_pred
    }
    comb_seurat@meta.data$smoothed_label <- NULL
    label_integers <- levels(comb_seurat@meta.data[, cluster_label_name])
    comb_seurat@meta.data[, output_label_name] <- label_integers[ comb_seurat@meta.data[, output_label_name] ]
    return(comb_seurat)
}


score_layer_fit <- function(comb_seurat, cluster_label_name, image_name=NULL, verbose=FALSE) {
    allowed_neighbors <- list(
        "Layer 1" = c("Layer 1", "Layer 2"),
        "Layer 2" = c("Layer 1", "Layer 2", "Layer 3"),
        "Layer 3" = c("Layer 2", "Layer 3", "Layer 4"),
        "Layer 4" = c("Layer 3", "Layer 4", "Layer 5"),
        "Layer 5" =  c("Layer 4", "Layer 5", "Layer 6"),
        "Layer 6" = c("Layer 5", "Layer 6", "WM"),
        "WM" = c("Layer 6", "WM"))

    cluster_labels <- comb_seurat@meta.data[, c("orig.ident", cluster_label_name)]

    layer_scores <- list()
    if(is.null(image_name)) {
        img_list <- names(comb_seurat@images)
    } else {
        img_list <- c(image_name)
    }
    for (img in img_list) {
        coord_table <- comb_seurat@images[[img]]@coordinates
        coord_table <- merge(coord_table, cluster_labels, by=0)
        rownames(coord_table) <- coord_table$Row.names
        
        # hexagonal
        coord_table$euclid_x <- (coord_table$row %% 2) / 2 + (coord_table$col - (coord_table$col %% 2)) / 2
        coord_table$euclid_y <- sqrt(3) * ( (coord_table$row %% 2) / 2 + (coord_table$row - (coord_table$row %% 2)) / 2)

        coord_dists <- dist(coord_table[, c("euclid_x", "euclid_y")])
        coord_dists <- as.matrix(coord_dists, nrow=nrow(coord_table))
        # round to nearest 0.1
        coord_dists <- round(coord_dists * 10) / 10
        good_adjs <- 0
        total_adjs <- 0
        for(spot_id in rownames(coord_table)) {
            spot_label <- coord_table[spot_id, cluster_label_name]
            nbr_ids <- colnames(coord_dists)[ which( abs(coord_dists[spot_id, ] - 1) < 0.01 ) ]
            nbr_labels <- coord_table[nbr_ids, cluster_label_name]
            good_adjs <- good_adjs + sum(nbr_labels %in% allowed_neighbors[[spot_label]])
            total_adjs <- total_adjs + length(nbr_labels)
            if(verbose) {
                print('-----')
                print(paste(spot_id, spot_label))
                print(nbr_ids)
                print(nbr_labels)
                print(sum(nbr_labels %in% allowed_neighbors[[spot_label]]))
            }
        }
        layer_scores[[img]] <- (good_adjs / total_adjs)
    }
    return(layer_scores)
}

score_layer_fraction <- function(comb_seurat, orig_col, new_col, image_name, layer_set, penalty_thresh=1.4) {
    img_rows <- rownames(comb_seurat@images[[image_name]]@coordinates)
    meta_sub <- comb_seurat@meta.data[img_rows, c(orig_col, new_col)]
    total_score <- 0
    for(layer in layer_set) {
        new_count <- sum(meta_sub[, new_col] %in% c(layer))
        old_count <- sum(meta_sub[, orig_col] %in% c(layer))
        if(old_count > 0) {
            total_score <- total_score + max(0, abs(new_count/old_count - 1) - (penalty_thresh-1))
        }
    }
    return( total_score )
}


comb_seurat <- readRDS("~/maynard_smoothing_test_01142025/clustering_pipeline_run/all_data_1covar2p0_V4_subL234_nfsig100_seed101_01152025/comb_seurat_subclust_01152025.rds")
svm_gamma <- 3
svm_cost <- 50
cluster_label_name <- "combined_clust"

layer_smoothing_path <- "/home/jy743/maynard_smoothing_test_01142025/layer_smoothing_maynard/V3"
comb_seurat@meta.data$orig.ident <- comb_seurat@meta.data$sample_name
# comb_seurat@meta.data[, cluster_label_name] <- comb_seurat@meta.data$Spatial_snn_res.0.4
smoothed_path <- file.path(layer_smoothing_path)
dir.create(smoothed_path, recursive = TRUE)
setwd(smoothed_path)




# smooth clusters
col_name <- sprintf("smoothed_s1_svmg%s_svmc%s", svm_gamma, svm_cost)
comb_seurat <- svm_smooth_labels(comb_seurat,
                            cluster_label_name,
                            svm_gamma = svm_gamma,
                            svm_cost = svm_cost,
                            class_weights = NULL,
                            image_name = NULL,
                            output_label_name = col_name,
                            random_seed = random_seed)
comb_seurat@meta.data[, col_name] <- factor(comb_seurat@meta.data[, col_name])


# pairs <- list(
#         c("Layer1", "3"),
#         c("Layer2", "sub_2"),
#         c("Layer3", "sub_1"),
#         c("Layer4", "sub_1"),
#         c("Layer5", "2"),
#         c("Layer6", "4"),
#         c("WM", "5")
# )
comb_seurat@meta.data$s2_subclust_score <- 0
comb_seurat@meta.data[comb_seurat@meta.data$smoothed_s1_svmg3_svmc50 == "sub_2", "s2_subclust_score"] <- 2
comb_seurat@meta.data[comb_seurat@meta.data$smoothed_s1_svmg3_svmc50 == "sub_1", "s2_subclust_score"] <- 3
comb_seurat@meta.data[comb_seurat@meta.data$smoothed_s1_svmg3_svmc50 == "2", "s2_subclust_score"] <- 5
comb_seurat@meta.data[comb_seurat@meta.data$smoothed_s1_svmg3_svmc50 == "4", "s2_subclust_score"] <- 6
# cluster_labels <- comb_seurat@meta.data[, c("orig.ident", "s2_subclust_score")]
rbf_gamma <- 0.01 # gamma: low value means more long-range influence

run_contour_score <- function(comb_seurat, score_column="s2_subclust_score", rbf_gamma=0.01, random_seed=101) {
    comb_seurat@meta.data[, "s2_nbr_score"] <- NA
    for (img in names(comb_seurat@images)) {
        print(paste("Processing stage 2 subcluster score: ", img))
        set.seed(random_seed)

        # meta_row <- rownames(comb_seurat@images[[img]]@coordinates)[1]
        # sample_name_filt <- comb_seurat@meta.data[meta_row, "orig.ident"]
        # diagnosis_filt <- comb_seurat@meta.data[meta_row, "diagnosis"]

        cluster_labels <- comb_seurat@meta.data[, c("orig.ident", score_column)]
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
        zero_s2_score_rows <- which(coord_table[, score_column] == 0)
        adjA[, zero_s2_score_rows] <- 0
        adjA <- sweep(adjA, 1, rowSums(adjA), FUN="/")

        coord_table$s2_nbr_score <- adjA %*% coord_table$s2_subclust_score
        comb_seurat@meta.data[rownames(coord_table), "s2_nbr_score"] <- coord_table$s2_nbr_score
    }
    return(comb_seurat)
}
comb_seurat <- run_contour_score(comb_seurat, score_column="s2_subclust_score", rbf_gamma=rbf_gamma)

pdf(sprintf("Maynard2021_spatial_smoothing_s2_contour_rbfg%s_01172025.pdf", rbf_gamma))
for (img in names(comb_seurat@images)) {
    sp_pt <- SpatialFeaturePlot(comb_seurat, features="s2_nbr_score", image.alpha=0.3, images=img, stroke=NA)
    print(sp_pt)
}
dev.off()

# Set L1 to NA if s2_nbr_score greater than 4
#         c("Layer1", "3"),
#         c("WM", "5")
L1_gt_thresh <- 4
comb_seurat@meta.data$smoothed_label_s2 <- comb_seurat@meta.data$smoothed_s1_svmg3_svmc50
comb_seurat@meta.data[comb_seurat@meta.data$smoothed_s1_svmg3_svmc50 %in% c("3", "5") & comb_seurat@meta.data$s2_nbr_score < L1_gt_thresh, "smoothed_label_s2"] <- "3"
comb_seurat@meta.data[comb_seurat@meta.data$smoothed_s1_svmg3_svmc50 %in% c("3", "5") & comb_seurat@meta.data$s2_nbr_score >= L1_gt_thresh, "smoothed_label_s2"] <- "5"

saveRDS(comb_seurat, "comb_seurat_maynard_post_smoothing_01152025.rds")

# pdf(sprintf("Maynard2021_spatial_smoothing_svmg%s_svmc%s_01152025.pdf", svm_gamma, svm_cost))
# for (img in names(comb_seurat@images)) {
#     sp_pt <- SpatialPlot(comb_seurat, group.by=col_name, image.alpha=0.3, images=img, stroke=NA)
#     print(sp_pt)
# }
# dev.off()

pdf("Maynard2021_raw_cluster_res0.4_final_01172025.pdf")
for (img in names(comb_seurat@images)) {
    sp_pt <- SpatialPlot(comb_seurat, group.by="combined_clust", image.alpha=0.3, images=img, stroke=NA)
    print(sp_pt)
}
dev.off()

pdf(sprintf("Maynard2021_spatial_smoothing_s2_rbfg%s_final_01172025.pdf", rbf_gamma))
for (img in names(comb_seurat@images)) {
    sp_pt <- SpatialPlot(comb_seurat, group.by="smoothed_label_s2", image.alpha=0.3, images=img, stroke=NA)
    print(sp_pt)
}
dev.off()


# pdf(sprintf("Maynard2021_spatial_smoothing_mats_svmg%s_svmc%s_01152025.pdf", svm_gamma, svm_cost))
mat1 <- table(comb_seurat@meta.data[, c("layer_guess", cluster_label_name)])
# pheatmap(mat1, color = colorRampPalette(c("blue", "white", "red"))(50),
#     show_rownames = TRUE,
#     show_colnames = TRUE,
#     display_numbers = mat1)
# mat2 <- table(comb_seurat@meta.data[, c("layer_guess", col_name)])
mat2 <- table(comb_seurat@meta.data[, c("layer_guess", "smoothed_label_s2")])
# pheatmap(mat2, color = colorRampPalette(c("blue", "white", "red"))(50),
#     show_rownames = TRUE,
#     show_colnames = TRUE,
#     diplay_numbers = mat2)
# mat3 <- mat2 - mat1
# pheatmap(mat3, color = colorRampPalette(c("blue", "white", "red"))(50),
#     show_rownames = TRUE,
#     show_colnames = TRUE,
#     diplay_numbers = mat3)
# dev.off()
write.csv(mat1, sprintf("mat_raw_%s_%s.csv", svm_gamma, svm_cost))
write.csv(mat2, sprintf("mat_sth_%s_%s.csv", svm_gamma, svm_cost))
# write.csv(mat3, sprintf("mat_dif_%s_%s_%s.csv", subj_name, svm_gamma, svm_cost))

# plot before and after
# Br5292 Br 5595 Br8100
pairs <- list(
        c("Layer1", "3"),
        c("Layer2", "sub_2"),
        c("Layer3", "sub_1"),
        c("Layer4", "sub_1"),
        c("Layer5", "2"),
        c("Layer6", "4"),
        c("WM", "5")
)
pair_strings <- sapply(pairs, function(x) paste(x, collapse = "|"))

# sample layer_guess smoothed num_match_before num_match_after
vals <- comb_seurat@meta.data %>%
    select(sample_name, layer_guess, combined_clust, smoothed_label_s2) %>%
    pivot_longer(
        cols = c(combined_clust, smoothed_label_s2),
        names_to = "column_name",
        values_to = "value") %>%
    group_by(sample_name, layer_guess, column_name, value) %>%
    summarize(count = n(), .groups = 'drop') %>%
# get the count per sample (group also by column_name so you don't double-count)
    group_by(sample_name, column_name, layer_guess) %>%
    mutate(sample_sum = sum(count)) %>%
    mutate(count_pct = count / sample_sum) %>%
    ungroup() %>%
# keep only the matching values
    mutate(layer_match = paste(layer_guess, value, sep = "|")) %>%
    filter(layer_match %in% pair_strings) %>%
    as.data.frame()

# confirm against table view
vals_check <- vals %>% group_by(layer_guess, column_name, value) %>%
    summarise( grp_count=sum(count)) %>% as.data.frame()


# V2
# https://stackoverflow.com/a/67996373
ttest_results <- vals %>%
group_by(layer_match) %>%
t_test(count_pct ~ column_name, paired = TRUE) %>%
add_significance()
pdf("mat_box_plots_paired_V2_01152025.pdf")
set.seed(101)
vals$layer_guess <- factor(vals$layer_guess, levels=c(paste0("Layer", 1:6), "WM"))
vals$sample_name <- factor(vals$sample_name)

ggplot(vals, aes(x = column_name, y = count_pct)) +
geom_boxplot(aes(fill = column_name), alpha = .2) +
geom_line(aes(group = sample_name)) +
geom_point(size = 2) +
facet_wrap(~ layer_match) +
geom_text(
  data = ttest_results,  # Use the t-test results
  aes(x = 1.5, y = max(vals$count_pct) * 1.1, label = paste(p, p.signif)),  # Position the annotation
  size = 4, color = "black"
) +
theme_minimal()
dev.off()






# ARI plot
# sample_name comparison ARI
# https://www.datanovia.com/en/blog/how-to-add-p-values-onto-basic-ggplots/
library(aricode)
library(rstatix)
library(ggpubr)
comparisons <- c("combined_clust", "smoothed_label_s2")
ari_df <- c()
for(sid in unique(comb_seurat@meta.data$sample_name)) {
    meta_sub <- comb_seurat@meta.data[comb_seurat@meta.data$sample_name == sid, ]
    meta_sub <- meta_sub[!is.na(meta_sub$layer_guess), ]
    for(comp in comparisons) {
        ari_val <- ARI(meta_sub$layer_guess, meta_sub[, comp])
        ari_df <- rbind(ari_df, data.frame(sample_name = sid,
                                           comparison = comp,
                                           ari = ari_val))
    }
}
ari_df$comparison <- factor(ari_df$comparison, levels=comparisons)

stat.test <- ari_df %>%
  t_test(ari ~ comparison, paired = TRUE) %>%
  add_significance()

pdf("ari_plot_paired_01172025.pdf")
# Box plots with p-values
bxp <-  ggpaired(ari_df, x = "comparison", y = "ari", id="sample_name", fill = "#E7B800",
                 line.color = "gray", line.size = 0.4)
stat.test <- stat.test %>% add_xy_position(x = "comparison")
bxp + stat_pvalue_manual(stat.test, label = "p.signif")

# Show the p-value combined with the significance level
bxp +
  stat_pvalue_manual(stat.test, label = "{p}{p.signif}") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10)))
dev.off()

write.csv(ari_df, "layer_smoothing_maynard_ari_df_01152025.csv")
saveRDS(comb_seurat, "comb_seurat_maynard_smoothed_validation_01152025.rds")
