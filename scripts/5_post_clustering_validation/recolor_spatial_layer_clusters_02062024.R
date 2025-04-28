# see "~/Downloads/recolor_spatial_layer_clusters_02062024.R"

library(Seurat)
library(dplyr)
library(ggplot2)
library(testit)

# eristwo
comb_seurat <- readRDS("/data/neurogen/jy1008/scRNAseq_integration/temp_sandbox_02052024/comb_seurat_post_smoothed_label_s5_11302023.rds")

# generate umap plot
group.colors.unsmoothed <- c("Layer 1" = "#f43545", "Layer 2"="#ff8901", "Layer 3"="#fad717", "Layer 4"="#00ba71", "Layer 5"="#00c2de", "Layer 6"="#ba55d3", "WM_1"="#3b2175", "WM_2"="#75215b", "WM_3"="#652175", "Noise"="#c8c8c8")

pdf("umap_major_clusters_V4_02062024.pdf")
DimPlot(comb_seurat, reduction="umap", label=TRUE, pt.size=.1, group.by="layer_label_v1", cols=group.colors.unsmoothed)
DimPlot(comb_seurat, reduction="umap", label=TRUE, pt.size=.1, group.by="layer_label_v1", cols=group.colors.unsmoothed) + NoLegend()
dev.off()

pdf("umap_v2_clusters_newL4_V4_01082025.pdf")
DimPlot(comb_seurat, reduction="umap", label=TRUE, group.by="layer_label_v2", cols=group.colors.unsmoothed)
DimPlot(comb_seurat, reduction="umap", label=TRUE, group.by="layer_label_v2", cols=group.colors.unsmoothed) + NoLegend()
DimPlot(comb_seurat, reduction="umap", label=TRUE, group.by="smoothed_label_s5", cols=group.colors)
DimPlot(comb_seurat, reduction="umap", label=TRUE, group.by="smoothed_label_s5", cols=group.colors) + NoLegend()
dev.off()

# Assignment of Layer 5 subclusters to new Layer 4 and Layer 5.
# see on eristwo: /data/neurogen/jy1008/scRNAseq_integration/spatial_clustering_validation_06302022/validate_new_layer4_12132022.R
# # extract Layer 4 from Layer 5; this is stored in layer_label_v2
# layer <- "Layer 5"
# subcluster_res <- 0.2
# cluster_col <- paste0("Spatial_snn_res.", subcluster_res)
# comb_sub <- readRDS(sprintf("/mnt/data0/projects/ASAP/ST/clustering_pipeline/post_clustering_validation_V4_seed101_10172023/seurat_subclust_5covar0p4_V4_%s_seed101_10272023/comb_seurat_diet_out.rds", gsub(' ','',layer)))
# Idents(comb_sub) <- comb_sub@meta.data[, cluster_col]
# comb_sub@meta.data[, "layer_label_v2"] <- as.character(comb_sub@meta.data[, cluster_col])
# comb_sub@meta.data[comb_sub@meta.data$layer_label_v2 == 3, "layer_label_v2"] <- "Layer 4"
# comb_sub@meta.data[comb_sub@meta.data$layer_label_v2 == 1, "layer_label_v2"] <- "Layer 5"
# comb_sub@meta.data[comb_sub@meta.data$layer_label_v2 == 2, "layer_label_v2"] <- "Layer 5"
# comb_sub@meta.data[comb_sub@meta.data$layer_label_v2 == 4, "layer_label_v2"] <- "Layer 5"
# comb_sub@meta.data[comb_sub@meta.data$layer_label_v2 == 5, "layer_label_v2"] <- "Layer 5"



comb_L5 <- readRDS("/data/neurogen/jy1008/scRNAseq_integration/temp_sandbox_02052024/seurat_subclust_5covar0p4_V4_Layer5_seed101_10272023/comb_seurat_diet_out.rds")
comb_L5@meta.data$L45_clusts <- as.character(comb_L5@meta.data$Spatial_snn_res.0.2)
comb_L5@meta.data[comb_L5@meta.data$L45_clusts == 3, "L45_clusts"] <- "Layer 4"
comb_L5@meta.data[comb_L5@meta.data$L45_clusts == 1, "L45_clusts"] <- "Layer 5"
comb_L5@meta.data[comb_L5@meta.data$L45_clusts == 2, "L45_clusts"] <- "Layer 5"
comb_L5@meta.data[comb_L5@meta.data$L45_clusts == 4, "L45_clusts"] <- "Layer 5"
comb_L5@meta.data[comb_L5@meta.data$L45_clusts == 5, "L45_clusts"] <- "Layer 5"

pdf("umap_minor_clusters_V4_02062024.pdf")
DimPlot(comb_L5, reduction="umap", label=TRUE, pt.size=.1, group.by="L45_clusts", cols=group.colors.unsmoothed)
DimPlot(comb_L5, reduction="umap", label=FALSE, pt.size=.1, group.by="L45_clusts", cols=group.colors.unsmoothed) + NoLegend()
dev.off()










comb_seurat <- readRDS("/Users/jy1008/Downloads/all_controls_5covars0p4_seed101_05082022/comb_seurat_diet_out.rds")

# 5 covars 0.4 theta, resolution 0.4
comb_seurat@meta.data["display_label"] <- comb_seurat@meta.data$Spatial_snn_res.0.4
comb_seurat@meta.data$display_label <- recode_factor(comb_seurat@meta.data$display_label, "6"="Layer 1", "4"="Layer 2/3", "1"="Layer 4", "2"="Layer 5", "5"="Layer 6", "7"="WM_1", "3"="WM_2")
# basic rainbow colors
# group.colors <- c("Layer 1" = "#ff0000", "Layer 2/3"="#ff7f00", "Layer 4"="#ffff00", "Layer 5"="#00f00", "Layer 6"="#0000ff", "WM_1"="#4b0082", "WM_2"="#9400d3")
# a nicer color set
group.colors <- c("Layer 1" = "#f43545", "Layer 2/3"="#ff8901", "Layer 4"="#fad717", "Layer 5"="#00ba71", "Layer 6"="#00c2de", "WM_1"="#00418d", "WM_2"="#ba55d3")

# test a sample
SpatialPlot(comb_seurat, group.by="display_label", image.alpha=0.3, images="slice1.28") + scale_fill_manual(values=group.colors)
SpatialFeaturePlot(comb_seurat, features="SNCA", images="slice1.22", slot="counts")
SpatialFeaturePlot(comb_seurat, features="SNCB", images="slice1.22", slot="counts")
SpatialFeaturePlot(comb_seurat, features="SNCG", images="slice1.22", slot="counts")

# two color continuous scale
# SpatialFeaturePlot(comb_seurat, features="SNCA", images="slice1.22", slot="counts") + theme(legend.title = element_blank()) + ggplot2::scale_fill_continuous(limits = c(0.0,1.0),, breaks = c(0.0, 0.5, 1.0))
# SpatialFeaturePlot(comb_seurat, features="SNCA", images="slice1.28", slot="counts") + theme(legend.title = element_blank())+ ggplot2::scale_fill_gradient(low="#ffffff", high="#00008b", limits=c(0,35))

pdf("recolored_layers_all_controls_5covars0p4_seed101_res0p4_05112022.pdf")
for (img in names(comb_seurat@images)) {
    sp_pt <- SpatialPlot(comb_seurat, group.by="display_label", image.alpha=0.3, images=img) + scale_fill_manual(values=group.colors) + theme(legend.title = element_blank())
    print(sp_pt)
}
dev.off()

pdf("all_controls_snca_05112022.pdf")
for (img in names(comb_seurat@images)) {
    sp_pt <- SpatialFeaturePlot(comb_seurat, features="SNCA", images=img, slot="counts")
    print(sp_pt)
}
dev.off()

pdf("all_controls_sncb_05112022.pdf")
for (img in names(comb_seurat@images)) {
    sp_pt <- SpatialFeaturePlot(comb_seurat, features="SNCB", images=img, slot="counts")
    print(sp_pt)
}
dev.off()

pdf("all_controls_sncg_05112022.pdf")
for (img in names(comb_seurat@images)) {
    sp_pt <- SpatialFeaturePlot(comb_seurat, features="SNCG", images=img, slot="counts")
    print(sp_pt)
}
dev.off()




# 5 covars 0.2 theta, resolution 0.4
comb_seurat <- readRDS("/Users/jy1008/Downloads/all_controls_5covars0p2_seed101_05072022/comb_seurat_diet_out.rds")

comb_seurat@meta.data["display_label"] <- comb_seurat@meta.data$Spatial_snn_res.0.4
comb_seurat@meta.data$display_label <- recode_factor(comb_seurat@meta.data$display_label, "7"="Layer 1", "5"="Layer 2/3", "1"="Layer 4", "3"="Layer 5", "2"="Layer 6", "6"="WM_1", "4"="WM_2", "8"="Noise")
# basic rainbow colors
# group.colors <- c("Layer 1" = "#ff0000", "Layer 2/3"="#ff7f00", "Layer 4"="#ffff00", "Layer 5"="#00f00", "Layer 6"="#0000ff", "WM_1"="#4b0082", "WM_2"="#9400d3")
# a nicer color set
group.colors <- c("Layer 1" = "#f43545", "Layer 2/3"="#ff8901", "Layer 4"="#fad717", "Layer 5"="#00ba71", "Layer 6"="#00c2de", "WM_1"="#00418d", "WM_2"="#ba55d3", "Noise"="#D3D3D3")

pdf("recolored_layers_all_controls_5covars0p2_seed101_res0p4_05112022.pdf")
for (img in names(comb_seurat@images)) {
    sp_pt <- SpatialPlot(comb_seurat, group.by="display_label", image.alpha=0.3, images=img) + scale_fill_manual(values=group.colors) + theme(legend.title = element_blank())
    print(sp_pt)
}
dev.off()



# 5 covars 0.2 theta, resolution 0.45
comb_seurat <- readRDS("/Users/jy1008/Downloads/all_controls_5covars0p2_seed101_05072022/comb_seurat_diet_out.rds")

# 5 covars 0.2 theta, resolution 0.4
comb_seurat@meta.data["display_label"] <- comb_seurat@meta.data$Spatial_snn_res.0.45
comb_seurat@meta.data$display_label <- recode_factor(comb_seurat@meta.data$display_label, "7"="Layer 1", "4"="Layer 2/3", "3"="Layer 4", "2"="Layer 5", "1"="Layer 6", "5"="WM_1", "6"="WM_2", "8"="Noise")
# basic rainbow colors
# group.colors <- c("Layer 1" = "#ff0000", "Layer 2/3"="#ff7f00", "Layer 4"="#ffff00", "Layer 5"="#00f00", "Layer 6"="#0000ff", "WM_1"="#4b0082", "WM_2"="#9400d3")
# a nicer color set
group.colors <- c("Layer 1" = "#f43545", "Layer 2/3"="#ff8901", "Layer 4"="#fad717", "Layer 5"="#00ba71", "Layer 6"="#00c2de", "WM_1"="#00418d", "WM_2"="#ba55d3", "Noise"="#D3D3D3")

pdf("recolored_layers_all_controls_5covars0p2_seed101_res0p45_05112022.pdf")
for (img in names(comb_seurat@images)) {
    sp_pt <- SpatialPlot(comb_seurat, group.by="display_label", image.alpha=0.3, images=img) + scale_fill_manual(values=group.colors) + theme(legend.title = element_blank())
    print(sp_pt)
}
dev.off()








# 6/22/2022 layers of all samples
comb_seurat_all <- readRDS("/Users/jy1008/Downloads/all_data_5covars0p4_seed101_05232022/comb_seurat_diet_no_reducs.rds") # local directory
comb_seurat_all@meta.data["num_label"] <- as.character(comb_seurat_all@meta.data$Spatial_snn_res.0.45)
comb_seurat_all@meta.data["layer_label"] <- 'Noise'
comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 8] <- 'Layer 1'
comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 4] <- 'Layer 2/3'
comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 2] <- 'Layer 4'
comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 1] <- 'Layer 5'
comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 3] <- 'Layer 6'
comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 6] <- 'WM_1'
comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 5] <- 'WM_2'
comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 7] <- 'WM_3'
comb_seurat_all@meta.data["num_label"] <- NULL
comb_seurat_all@meta.data$layer_label <- as.factor(comb_seurat_all@meta.data$layer_label)
Idents(comb_seurat_all) <- comb_seurat_all@meta.data$layer_label
group.colors <- c("Layer 1" = "#f43545", "Layer 2/3"="#ff8901", "Layer 4"="#fad717", "Layer 5"="#00ba71", "Layer 6"="#00c2de", "WM_1"="#00418d", "WM_2"="#ba55d3", "WM_3"="#652175", "Noise"="#c8c8c8")
#e03694
# test a sample
pdf("/Users/jy1008/Downloads/all_data_5covars0p4_seed101_05232022/recolored_layers_all_data_5covars0p4_seed101_res0p45_05112022.pdf")
for (img in names(comb_seurat_all@images)) {
    sp_pt <- SpatialPlot(comb_seurat_all, group.by="layer_label", image.alpha=0.3, images=img) + scale_fill_manual(values=group.colors) + theme(legend.title = element_blank())
    print(sp_pt)
}
dev.off()






# 7/7/2022 layers of all samples (fixed)
### 0.35
comb_seurat_all <- readRDS("/Users/jy1008/Downloads/all_data_5covars0p4_seed101_06252022/comb_seurat_diet_out_no_reducs.rds")
comb_seurat_all@meta.data["num_label"] <- as.character(comb_seurat_all@meta.data$Spatial_snn_res.0.35)
comb_seurat_all@meta.data["layer_label"] <- 'Noise'
comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 6] <- 'Layer 1'
comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 4] <- 'Layer 2/3'
comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 3] <- 'Layer 4'
comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 1] <- 'Layer 5'
comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 2] <- 'Layer 6'
comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 7] <- 'WM_1'
comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 5] <- 'WM_2'
comb_seurat_all@meta.data["num_label"] <- NULL
comb_seurat_all@meta.data$layer_label <- as.factor(comb_seurat_all@meta.data$layer_label)
Idents(comb_seurat_all) <- comb_seurat_all@meta.data$layer_label
group.colors <- c("Layer 1" = "#f43545", "Layer 2/3"="#ff8901", "Layer 4"="#fad717", "Layer 5"="#00ba71", "Layer 6"="#00c2de", "WM_1"="#00418d", "WM_2"="#ba55d3", "Noise"="#c8c8c8")
pdf("/Users/jy1008/Downloads/all_data_5covars0p4_seed101_06252022/recolored_layers_all_data_5covars0p4_seed101_res0p35_07072022.pdf")
for (img in names(comb_seurat_all@images)) {
    sp_pt <- SpatialPlot(comb_seurat_all, group.by="layer_label", image.alpha=0.3, images=img) + scale_fill_manual(values=group.colors) + theme(legend.title = element_blank())
    print(sp_pt)
}
dev.off()



### 0.4
comb_seurat_all <- readRDS("/Users/jy1008/Downloads/all_data_5covars0p4_seed101_06252022/comb_seurat_diet_out_no_reducs.rds")
comb_seurat_all@meta.data["num_label"] <- as.character(comb_seurat_all@meta.data$Spatial_snn_res.0.4)
comb_seurat_all@meta.data["layer_label"] <- 'Noise'
comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 7] <- 'Layer 1'
comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 3] <- 'Layer 2/3'
comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 2] <- 'Layer 4'
comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 1] <- 'Layer 5'
comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 4] <- 'Layer 6'
comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 6] <- 'WM_1'
comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 8] <- 'WM_2'
comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 5] <- 'WM_3'
comb_seurat_all@meta.data["num_label"] <- NULL
comb_seurat_all@meta.data$layer_label <- as.factor(comb_seurat_all@meta.data$layer_label)
Idents(comb_seurat_all) <- comb_seurat_all@meta.data$layer_label
group.colors <- c("Layer 1" = "#f43545", "Layer 2/3"="#ff8901", "Layer 4"="#fad717", "Layer 5"="#00ba71", "Layer 6"="#00c2de", "WM_1"="#00418d", "WM_2"="#ba55d3", "WM_3"="#652175", "Noise"="#c8c8c8")
pdf("/Users/jy1008/Downloads/all_data_5covars0p4_seed101_06252022/recolored_layers_all_data_5covars0p4_seed101_res0p4_07072022.pdf")
for (img in names(comb_seurat_all@images)) {
    sp_pt <- SpatialPlot(comb_seurat_all, group.by="layer_label", image.alpha=0.3, images=img) + scale_fill_manual(values=group.colors) + theme(legend.title = element_blank())
    print(sp_pt)
}
dev.off()




### 0.45
comb_seurat_all <- readRDS("/Users/jy1008/Downloads/all_data_5covars0p4_seed101_06252022/comb_seurat_diet_out_no_reducs.rds")
comb_seurat_all@meta.data["num_label"] <- as.character(comb_seurat_all@meta.data$Spatial_snn_res.0.45)
comb_seurat_all@meta.data["layer_label"] <- 'Noise'
# comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 9] <- 'Layer 1A'
# NOTE: cluster 9 looks like Layer 1, but it is a small cluster and widely dispersed spatially
comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 7] <- 'Layer 1'
comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 3] <- 'Layer 2'
comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 2] <- 'Layer 3'
comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 1] <- 'Layer 5'
comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 4] <- 'Layer 6'
comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 8] <- 'WM_1'
comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 6] <- 'WM_2'
comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 5] <- 'WM_3'
comb_seurat_all@meta.data["num_label"] <- NULL
comb_seurat_all@meta.data$layer_label <- as.factor(comb_seurat_all@meta.data$layer_label)
Idents(comb_seurat_all) <- comb_seurat_all@meta.data$layer_label
group.colors <- c("Layer 1" = "#f43545", "Layer 2/3"="#ff8901", "Layer 4"="#fad717", "Layer 5"="#00ba71", "Layer 6"="#00c2de", "WM_1"="#00418d", "WM_2"="#ba55d3", "WM_3"="#652175", "Noise"="#c8c8c8")
pdf("/Users/jy1008/Downloads/all_data_5covars0p4_seed101_06252022/recolored_layers_all_data_5covars0p4_seed101_res0p45_07072022.pdf")
for (img in names(comb_seurat_all@images)) {
    sp_pt <- SpatialPlot(comb_seurat_all, group.by="layer_label", image.alpha=0.3, images=img) + scale_fill_manual(values=group.colors) + theme(legend.title = element_blank())
    print(sp_pt)
}
dev.off()

# The newest layers: Layer 4 is a sub-cluster of Layer 5 and should be shown in a separate plot
comb_seurat_all@meta.data$layer_label2 <- as.character(comb_seurat_all@meta.data$layer_label)
comb_seurat_all@meta.data["layer_label2"][comb_seurat_all@meta.data["layer_label"] == 'WM_1'] <- 'WM'
comb_seurat_all@meta.data["layer_label2"][comb_seurat_all@meta.data["layer_label"] == 'WM_2'] <- 'WM'
comb_seurat_all@meta.data["layer_label2"][comb_seurat_all@meta.data["layer_label"] == 'WM_3'] <- 'WM'
comb_seurat_all@meta.data$layer_label2 <- as.factor(comb_seurat_all@meta.data$layer_label2)
Idents(comb_seurat_all) <- comb_seurat_all@meta.data$layer_label2
group.colors <- c("Layer 1" = "#f43545", "Layer 2"="#ff8901", "Layer 3"="#fad717", "Layer 4"="#00ba71", "Layer 5"="#00c2de", "Layer 6"="#ba55d3", "WM"="#652175", "Noise"="#c8c8c8")

# comb_seurat_all <- readRDS("/Users/jy1008/Downloads/all_data_5covars0p4_seed101_06252022/comb_seurat_diet_out.rds")
# pdf("/Users/jy1008/Downloads/all_data_5covars0p4_seed101_06252022/recolored_umap_all_data_5covars0p4_seed101_res0p45_07082022.pdf")
pdf("recolored_umap_all_data_5covars0p4_seed101_res0p45_L4inL5_04062023.pdf")
DimPlot(comb_seurat_all, reduction="umap", label=FALSE, pt.size=.1, group.by="layer_label2", cols=group.colors) + labs(title="") + theme(legend.position = "none")
DimPlot(comb_seurat_all, reduction="umap", label=TRUE, pt.size=.1, group.by="layer_label2", cols=group.colors)
dev.off()

# Clustering layer 5 to extract layer 4
comb_layer5 <- readRDS("/data/neurogen/jy1008/scRNAseq_integration/sublayers_test_11102022/all_data_5covars0p4_Layer5_11142022/comb_seurat_diet_out.rds")
comb_layer5@meta.data["layer_label"] <- 'Noise'
comb_layer5@meta.data["layer_label"][comb_layer5@meta.data["Spatial_snn_res.0.2"] == 1] <- 'Layer 5'
comb_layer5@meta.data["layer_label"][comb_layer5@meta.data["Spatial_snn_res.0.2"] == 2] <- 'Layer 5'
comb_layer5@meta.data["layer_label"][comb_layer5@meta.data["Spatial_snn_res.0.2"] == 3] <- 'Layer 4'
comb_layer5@meta.data["layer_label"][comb_layer5@meta.data["Spatial_snn_res.0.2"] == 4] <- 'Layer 5'
comb_layer5@meta.data["layer_label"][comb_layer5@meta.data["Spatial_snn_res.0.2"] == 5] <- 'Layer 5'
comb_layer5@meta.data$layer_label <- as.factor(comb_layer5@meta.data$layer_label)
Idents(comb_layer5) <- comb_layer5@meta.data$layer_label
pdf("recolored_umap_5covars0p4_seed101_res0p45_L5only_L4sublayer_04062023.pdf")
DimPlot(comb_layer5, reduction="umap", label=FALSE, pt.size=.1, group.by="layer_label", cols=group.colors) + labs(title="") + theme(legend.position = "none")
DimPlot(comb_layer5, reduction="umap", label=TRUE, pt.size=.1, group.by="layer_label", cols=group.colors)
dev.off()







### seed 201 res 0.5
# comb_seurat_all <- readRDS("/Users/jy1008/Downloads/all_data_5covars0p4_seed201_06252022/comb_seurat_diet_out_no_reducs.rds")
# comb_seurat_all@meta.data["num_label"] <- as.character(comb_seurat_all@meta.data$Spatial_snn_res.0.5)
# comb_seurat_all@meta.data["layer_label"] <- 'Noise'
# comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 9] <- 'Layer 1A'
# comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 7] <- 'Layer 1B'
# comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 10] <- 'Layer 1/3'
# comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 2] <- 'Layer 2/3A'
# comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 14] <- 'Layer 2/3B'
# comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 1] <- 'Layer 4'
# comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 4] <- 'Layer 5A'
# comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 5] <- 'Layer 5B'
# comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 6] <- 'Layer 5C'
# comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 11] <- 'Layer 5D'
# comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 8] <- 'Layer 6'
# comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 12] <- 'WM_A'
# comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 13] <- 'WM_B'
# comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 15] <- 'WM_C'
# comb_seurat_all@meta.data["layer_label"][comb_seurat_all@meta.data["num_label"] == 3] <- 'WM_D'
# comb_seurat_all@meta.data["num_label"] <- NULL
# comb_seurat_all@meta.data$layer_label <- as.factor(comb_seurat_all@meta.data$layer_label)
# Idents(comb_seurat_all) <- comb_seurat_all@meta.data$layer_label
# # https://colorswall.com/palette/102
# group.colors <- c(
#     "Layer 1A" = "#ff6666",
#     "Layer 1B" = "#ff0000",
#     "Layer 1/3" = "#bf360c",
#     "Layer 2/3A" = "#ffc266",
#     "Layer 2/3B" = "#ff9900",
#     "Layer 4" = "#ffff00",
#     "Layer 5A" = "#99cc99",
#     "Layer 5B" = "#66b366",
#     "Layer 5C" = "#339933",
#     "Layer 5D" = "#008000",
#     "Layer 6" = "#0000ff",
#     "WM_A" = "#825f87",
#     "WM_B" = "#704876",
#     "WM_C" = "#5e3165",
#     "WM_D" = "#4c1a54")
# pdf("/Users/jy1008/Downloads/all_data_5covars0p4_seed201_06252022/recolored_layers_all_data_5covars0p4_seed201_res0p5_07082022.pdf")
# pdf("/Users/jy1008/Downloads/all_data_5covars0p4_seed201_06252022/recolored_layers_all_data_5covars0p4_seed201_res0p5_09122022.pdf")
# for (img in names(comb_seurat_all@images)) {
#     sp_pt <- SpatialPlot(comb_seurat_all, group.by="layer_label", image.alpha=0.3, images=img) + scale_fill_manual(values=group.colors) + theme(legend.title = element_blank())
#     print(sp_pt)
# }
# dev.off()

# comb_seurat_all <- readRDS("/Users/jy1008/Downloads/all_data_5covars0p4_seed201_06252022/comb_seurat_diet_out.rds")
# pdf("/Users/jy1008/Downloads/all_data_5covars0p4_seed201_06252022/recolored_umap_all_data_5covars0p4_seed201_res0p5_07082022.pdf")
# DimPlot(comb_seurat_all, reduction="umap", label=FALSE, pt.size=.1, group.by="layer_label", cols=group.colors) + labs(title="") + theme(legend.position = "none")
# dev.off()



# from "~/Downloads/gauss_kernel_convolution_cluster_smoothing_gmsvm2_V1_01032023.R"

# recolor the layers
comb_seurat@meta.data$final_smoothed_label_plusL1 <- as.factor(comb_seurat@meta.data$final_smoothed_label_plusL1)
Idents(comb_seurat) <- comb_seurat@meta.data$final_smoothed_label_plusL1
group.colors <- c("Layer 1" = "#f43545", "Layer 2"="#ff8901", "Layer 3"="#fad717", "Layer 4"="#00ba71", "Layer 5"="#00c2de", "Layer 6"="#ba55d3", "WM"="#652175")
# pdf("recolored_layers_comb_seurat_smoothed_final_01102023.pdf")
pdf("recolored_layers_comb_seurat_smoothed_final_nospotborder_01302023.pdf")
for (img in names(comb_seurat@images)) {
    sp_pt <- SpatialPlot(comb_seurat, group.by="final_smoothed_label_plusL1", image.alpha=0.3, images=img, stroke=NA) + scale_fill_manual(values=group.colors) + theme(legend.title = element_blank())
    print(sp_pt)
}
dev.off()

# the original layers with a color scheme matching the smoothed output
Idents(comb_seurat) <- comb_seurat@meta.data$new_label
group.colors2 <- c("Layer 1" = "#f43545", "Layer 2"="#ff8901", "Layer 3"="#fad717", "Layer 4"="#00ba71", "Layer 5"="#00c2de", "Layer 6"="#ba55d3", "WM_1"="#3b2175", "WM_2"="#75215b", "WM_3"="#652175", "Noise"="#c8c8c8")
pdf("recolored_layers_plusnewL4_unsmoothed_comb_seurat_final_01302023.pdf")
for (img in names(comb_seurat@images)) {
    sp_pt <- SpatialPlot(comb_seurat, group.by="new_label", image.alpha=0.3, images=img, stroke=NA) + scale_fill_manual(values=group.colors2) + theme(legend.title = element_blank())
    print(sp_pt)
}
dev.off()
Idents(comb_seurat) <- comb_seurat@meta.data$final_smoothed_label_plusL1
