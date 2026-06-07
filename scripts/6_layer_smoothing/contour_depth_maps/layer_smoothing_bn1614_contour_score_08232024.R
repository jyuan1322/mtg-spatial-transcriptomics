library(Seurat)
library(ggplot2)
library(cowplot)

comb_seurat <- readRDS("/media/Mercury2LV/panda/projects/ASAP/ST/clustering_pipeline/layer_smoothing_V4_seed101_11302023_replication_test/comb_seurat_post_smoothed_label_s5_11302023.rds")
bn1614_meta <- comb_seurat@meta.data[comb_seurat@meta.data$sample_name == "BN1614", ]
sample_name_filt <- "BN1614"
diagnosis_filt <- "Case"

setwd("~")
bn1614_meta$layer_label_v3 <- as.character(bn1614_meta$layer_label_v2)
bn1614_meta[startsWith(bn1614_meta$layer_label_v3, "WM"), "layer_label_v3"] <- "WM"

group.colors <- c("Layer 1"="#f43545", "Layer 2"="#ff8901", "Layer 3"="#fad717",
                  "Layer 4"="#00ba71", "Layer 5"="#00c2de", "Layer 6"="#ba55d3",
                  "WM"="#652175")

pdf("layer_smoothing_bn1614_contour_score_08232024.pdf", width=5,height=3)
# https://stackoverflow.com/a/22181949
# normalize each histogram group to sum to 1
ggplot(bn1614_meta, aes(x=s2_nbr_score, fill=layer_label_v3)) + geom_histogram(aes(y=0.1*after_stat(density)), alpha=0.5, position="identity", binwidth=0.1) + ggtitle(paste(sample_name_filt, diagnosis_filt)) + theme_bw()
ggplot(bn1614_meta, aes(x=s2_nbr_score, fill=layer_label_v3)) + geom_histogram(aes(y=0.1*after_stat(density)), alpha=1.0, position="stack", binwidth=0.1) + ggtitle(paste(sample_name_filt, diagnosis_filt)) + theme_bw()

df_temp1 <- bn1614_meta[bn1614_meta$layer_label_v3 %in% c("Layer 1", "WM"), ]
df_temp2 <- bn1614_meta[!(bn1614_meta$layer_label_v3 %in% c("Layer 1", "WM")), ]
hp1 <- ggplot(df_temp1, aes(x=s2_nbr_score, fill=layer_label_v3)) +
    geom_histogram(aes(y=0.1*after_stat(density)), alpha=0.5,
                   position="identity", breaks = seq(2.4, 6, 0.1)) +
    scale_fill_manual(values=group.colors) +
    scale_x_continuous(breaks = c(2.4, 3, 4, 5, 6)) +
    ggtitle(paste(sample_name_filt, diagnosis_filt)) +
    labs(x=NULL, y=NULL) + theme_bw()
hp2 <- ggplot(df_temp2, aes(x=s2_nbr_score, fill=layer_label_v3)) +
    geom_histogram(aes(y=0.1*after_stat(density)), alpha=0.5,
                   position="identity", breaks = seq(2.4, 6, 0.1)) +
    scale_fill_manual(values=group.colors) +
    scale_x_continuous(breaks = c(2.4, 3, 4, 5, 6)) +
    labs(x=NULL, y=NULL) + theme_bw()
plot_grid(hp1, hp2, ncol = 1, align = "v")
plot_grid(hp1+theme(legend.position="none"),
          hp2+theme(legend.position="none"), ncol=1, align="v")
dev.off()
