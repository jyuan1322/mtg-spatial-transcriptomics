# see ~/Downloads/gauss_kernel_convolution_cluster_smoothing_V4_11102023.R

library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(stringr)
library(testit)

library(philentropy) # KL divergence

# for contour gradient plot
# library(pracma)
# library(akima)
# library(ggquiver)

random_seed <- 101

# on panda server
comb_seurat <- readRDS("/mnt/data0/projects/ASAP/ST/clustering_pipeline/layer_smoothing_V4_seed101_11302023_replication_test/comb_seurat_post_smoothed_label_s5_11302023.rds")

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
gene_ave_epr <- AggregateExpression(comb_seurat, group.by = c("sample_name", "contour_bin"))
gene_df <- t(gene_ave_epr$Spatial)

gene_aves <- rowMeans(comb_seurat@assays$Spatial@data)
gene_aves <- gene_aves[gene_aves > 0]
thresh <- 0.75
gene_aves <- gene_aves[gene_aves > quantile(gene_aves, thresh)]
gene_list <- names(gene_aves)

gene_df <- as.data.frame(gene_df[, gene_list])
gene_df$sample_name <- str_split_fixed(rownames(gene_df), '_', 2)[, 1]
gene_df$contour_bin <- str_split_fixed(rownames(gene_df), '_', 2)[, 2]

gene_kl_scores <- c()
i <- 0
for(gene in gene_list) {
    i <- i + 1
    print(paste(i, gene))

    gene_df$gene <- gene_df[, gene]
    temp <- gene_df %>% group_by(sample_name, contour_bin) %>%
                summarise(gene = sum(gene)) %>% ungroup %>%
                complete(sample_name, contour_bin, fill = list(gene = 0))
    # gene must be greater than 0 in at least 3 contour bins in at least 25% of samples
    num_samples <- length(unique(temp$sample_name))
    sample_thresh <- round(num_samples * 0.25) # num samples passing thresh
    indiv_thresh <- 3 # num nonzero bins needed for a sample to pass thresh
    temp_count <- temp[temp$gene > 0, ] %>% count(sample_name)
    num_samples_passing_thresh <- nrow(temp_count[temp_count$n >= indiv_thresh, ])
    if( num_samples_passing_thresh < sample_thresh) {
        print(paste(gene, "Failed expr thresh:", num_samples_passing_thresh, "<", sample_thresh))
        next
    }
    temp <- merge( temp, unique(meta_data[, c("sample_name", "diagnosis")]), by="sample_name" )
    # create normalized histograms over contour bins, one per sample
    temp <- temp %>% group_by(sample_name) %>% mutate(gene = gene/sum(gene))
    # calculate average per contour bin and diagnosis
    temp <- temp %>% group_by(contour_bin, diagnosis) %>%
                summarize(gene = mean(gene, na.rm=TRUE))
    temp <- temp %>% pivot_wider(names_from = contour_bin, values_from = gene) %>% as.data.frame()
    rownames(temp) <- temp$diagnosis
    temp$diagnosis <- NULL
    
    tryCatch(
        expr = {
            kl_cont_ilbd <- KL(as.matrix(temp[c("Control", "ILBD"), ]))
            kl_cont_case <- KL(as.matrix(temp[c("Control", "Case"), ]))
            gene_kl_scores <- rbind(gene_kl_scores, c(gene, kl_cont_ilbd, kl_cont_case))
        },
        error = function(e){
            print(paste(gene, "KL error"))
        }
    )
}
colnames(gene_kl_scores) <- c("gene", "kl_cont_ilbd", "kl_cont_case")
# Error in if (any(colSums(x) > 1.00001)) stop("Please make sure that all vectors sum up to 1.0 ...",  :
# missing value where TRUE/FALSE needed




# UMAP plot with contour coloring, split by diagnosis
# center of each layer per diagnosis is marked
# NOTE: generate contour_score2 first
meta_data <- comb_seurat@meta.data
umap_coords <- comb_seurat@reductions$umap@cell.embeddings
assert( all(rownames(meta_data) == rownames(umap_coords)) )
meta_data <- cbind(meta_data, umap_coords)
special_points <- c()
special_points <- rbind(special_points, c("global mean", NA, as.vector(colMeans(meta_data[, c("UMAP_1", "UMAP_2")]))))
diagnosis_list <- unique(meta_data$diagnosis)
layer_list <- unique(meta_data$smoothed_label_s5)
for(diag in diagnosis_list) {
    meta_data_sub <- meta_data[meta_data$diagnosis == diag, ]
    for(st_layer in layer_list) {
        meta_data_sub <- meta_data[meta_data$diagnosis == diag & meta_data$smoothed_label_s5 == st_layer, ]
        special_points <- rbind(special_points, c(diag, st_layer, as.vector(colMeans(meta_data_sub[, c("UMAP_1", "UMAP_2")]))))
    }
}
colnames(special_points) <- c("diagnosis", "st_layer", "UMAP_1", "UMAP_2")
special_points <- as.data.frame(special_points)
special_points$UMAP_1 <- as.numeric(special_points$UMAP_1)
special_points$UMAP_2 <- as.numeric(special_points$UMAP_2)

layer_colors <- c("Layer 1" = "#f43545", "Layer 2"="#ff8901", "Layer 3"="#fad717", "Layer 4"="#00ba71", "Layer 5"="#00c2de", "Layer 6"="#ba55d3", "WM"="#652175")

pdf("contour_umap_plot_02222024.pdf", width=10, height=10)
for(diag in diagnosis_list) {
    meta_data_sub <- meta_data[meta_data$diagnosis == diag, ]
    pt <- ggplot(meta_data_sub, aes(UMAP_1, UMAP_2)) +
            geom_point(aes(color = contour_score2), size = 0.1) +
            scale_color_gradientn(colors=c("#edf8b1", "#7fcdbb", "#2c7fb8")) +
            geom_point(data=special_points, aes(x=UMAP_1,y=UMAP_2, shape=factor(diagnosis))) +
            scale_fill_manual(values = layer_colors) +
            theme_bw() + ggtitle(diag)
    print(pt)
}
dev.off()






pdf(sprintf("layer_contour_scores_gamma%s_02152024.pdf", rbf_gamma))
for(img in names(comb_seurat@images)) {
    print(img)
    meta_row <- rownames(comb_seurat@images[[img]]@coordinates)[1]
    sample_name_filt <- comb_seurat@meta.data[meta_row, "orig.ident"]
    diagnosis_filt <- comb_seurat@meta.data[meta_row, "diagnosis"]
    sp_pt <- SpatialFeaturePlot(comb_seurat, features="contour_score2", image.alpha=0.3, images=img, stroke=NA) +
        scale_fill_gradientn(colors=c("#f43545", "#ff8901", "#fad717", "#00ba71", "#00c2de", "#ba55d3", "#652175")) +
        ggtitle(paste(sample_name_filt, diagnosis_filt))
    print(sp_pt)
    sp_pt <- SpatialFeaturePlot(comb_seurat, features="contour_score2", image.alpha=0.3, images=img, stroke=NA) +
        scale_fill_gradientn(colors=c("#edf8b1", "#7fcdbb", "#2c7fb8")) +
        ggtitle(paste(sample_name_filt, diagnosis_filt))
    print(sp_pt)
}
dev.off()

# TODO: controls only?
meta_data <- comb_seurat@meta.data
meta_data$spot_id <- rownames(meta_data)

# merge meta data with cell convolution information
cell_type_prefixes <- c("^Astro", "^Endoth", "^GABA", "^GLU", "^Immun", "^Micro", "^OPC", "^Oligo", "^Peri","^Unknown")
sc_deconv_dir <- "/mnt/data0/projects/ASAP/ST/st_sc_deconvolution/results_V2_08032023/cell2location_map"
sc_deconv_files <- list.files(path=sc_deconv_dir, pattern = "*.csv")
sc_dfs <- c()
for(sc_file in sc_deconv_files) {
    print(sc_file)
    sc_df <- read.csv(file.path(sc_deconv_dir, sc_file))
    
    cols <- colnames(sc_df)
    cell_types <- c()
    for(ct_prefix in cell_type_prefixes) {
        cell_types <- c(cell_types, cols[grepl(ct_prefix, cols)])
    }
    sc_df <- sc_df[, c("spot_id", cell_types)]
    sc_dfs <- bind_rows(sc_dfs, sc_df)
}
meta_data <- merge(meta_data, sc_dfs, by="spot_id")
meta_data$contour_bin <- cut(meta_data$contour_score2, breaks=100)


cols <- colnames(meta_data)
cell_types <- c()
for(ct_prefix in cell_type_prefixes) {
    cell_types <- c(cell_types, cols[grepl(ct_prefix, cols)])
}
meta_data2 <- meta_data[, c("contour_bin", cell_types)]
meta_data2 <- meta_data2 %>%
  pivot_longer(!contour_bin, names_to = "cell_type", values_to = "pseudoabundance") %>% as.data.frame()

# meta_data2 <- meta_data2[!is.na(meta_data2$pseudoabundance), ]

meta_data2 <- meta_data2 %>%
  group_by(contour_bin, cell_type) %>%
  summarise(sum_pseudoabundance = sum(pseudoabundance, na.rm = TRUE)) %>% as.data.frame()

# major clusters only
meta_data2$major_cell_type <- meta_data2$cell_type
meta_data2[startsWith(meta_data2$cell_type, "GABA"), "major_cell_type"] <- "GABA_Neurons"
meta_data2[startsWith(meta_data2$cell_type, "GLU"), "major_cell_type"] <- "GLU_Neurons"


# only GABA and GLU
meta_data_gaba <- meta_data2[meta_data2$major_cell_type == "GABA_Neurons", ]
meta_data_glu <- meta_data2[meta_data2$major_cell_type == "GLU_Neurons", ]


# pdf(sprintf("layer_contour_geomcol_gamma%s_02022024.pdf", rbf_gamma))
# pt <- ggplot() + geom_col(data = meta_data2, aes(x = contour_bin, y = pseudoabundance, fill = cell_type), position = "fill")
# print(pt)
# dev.off()

pdf(sprintf("layer_contour_geomcol_major_gamma%s_02152024.pdf", rbf_gamma))
pt <- ggplot() + geom_col(data = meta_data2, aes(x = contour_bin, y = sum_pseudoabundance, fill = major_cell_type), position = "fill") +
    scale_x_discrete(breaks = levels(meta_data2$contour_bin)[c(1, seq(from=5, to=100, by=5))])
print(pt)
dev.off()
pdf(sprintf("layer_contour_geomcol_gaba_gamma%s_02152024.pdf", rbf_gamma))
pt <- ggplot() + geom_col(data = meta_data_gaba, aes(x = contour_bin, y = sum_pseudoabundance, fill = cell_type), position = "fill") +
    scale_x_discrete(breaks = levels(meta_data2$contour_bin)[c(1, seq(from=5, to=100, by=5))])
print(pt)
dev.off()
pdf(sprintf("layer_contour_geomcol_glu_gamma%s_02152024.pdf", rbf_gamma))
pt <- ggplot() + geom_col(data = meta_data_glu, aes(x = contour_bin, y = sum_pseudoabundance, fill = cell_type), position = "fill") +
    scale_x_discrete(breaks = levels(meta_data2$contour_bin)[c(1, seq(from=5, to=100, by=5))])
print(pt)

temp <- meta_data_glu %>%
    group_by(cell_type) %>%
    mutate(pct_pseudoabundance = sum_pseudoabundance/sum(sum_pseudoabundance)) %>%
    as.data.frame()
temp$cell_type <- factor(temp$cell_type, levels = paste0("GLU_Neurons_", 1:19))
# https://sashamaps.net/docs/resources/20-colors/
pt <- ggplot(data = temp, aes(x = contour_bin, y = pct_pseudoabundance, fill=cell_type)) +
    geom_bar(stat="identity", position="stack", alpha=1.0) +
    scale_fill_manual("Cell Types", values =
        c(
        "GLU_Neurons_1" = "#e6194B",
        "GLU_Neurons_2" = "#3cb44b",
        "GLU_Neurons_3" = "#ffe119",
        "GLU_Neurons_4" = "#4363d8",
        "GLU_Neurons_5" = "#f58231",
        "GLU_Neurons_6" = "#911eb4",
        "GLU_Neurons_7" = "#42d4f4",
        "GLU_Neurons_8" = "#f032e6",
        "GLU_Neurons_9" = "#bfef45",
        "GLU_Neurons_10" = "#fabed4",
        "GLU_Neurons_11" = "#469990",
        "GLU_Neurons_12" = "#dcbeff",
        "GLU_Neurons_13" = "#9A6324",
        "GLU_Neurons_14" = "#fffac8",
        "GLU_Neurons_15" = "#800000",
        "GLU_Neurons_16" = "#aaffc3",
        "GLU_Neurons_17" = "#808000",
        "GLU_Neurons_18" = "#ffd8b1",
        "GLU_Neurons_19" = "#000075"))
print(pt)
dev.off()


# histogram of layers by the same bins
temp <- table(meta_data[, c("smoothed_label_s5", "contour_bin")]) %>% as.data.frame()
pdf(sprintf("layer_contour_geomcol_layers_gamma%s_02152024.pdf", rbf_gamma))
pt <- ggplot() + geom_col(data = temp, aes(x = contour_bin, y = Freq, fill = smoothed_label_s5), position = "fill") +
    scale_x_discrete(breaks = levels(meta_data2$contour_bin)[c(1, seq(from=5, to=100, by=5))])
print(pt)
pt <- ggplot() + geom_col(data = temp, aes(x = contour_bin, y = Freq, fill = smoothed_label_s5), position = "dodge") +
    scale_x_discrete(breaks = levels(meta_data2$contour_bin)[c(1, seq(from=5, to=100, by=5))])
print(pt)
dev.off()

