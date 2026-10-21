library(ggplot2)
library(corrplot)
library(RColorBrewer)

# NOTE: this file should be named 01102023?
# layer_corrs <- read.csv("~/Downloads/maynard_validation_pearson_svm_smoothed_01102022.csv")

generate_corr_heatmap <- function(corr_file) {
    layer_corrs <- read.csv(corr_file)
    rownames(layer_corrs) <- layer_corrs$X
    layer_corrs$X <- NULL

    # column_order <- c(paste0("Layer.", 1:6), "WM")
    # row_order <- c(paste0("t_stat_Layer", 1:6), "t_stat_WM")
    cols <- colnames(layer_corrs)
    column_order <- c(sort(cols[startsWith(cols, "Layer")]),
                      sort(cols[startsWith(cols, "WM")]),
                      sort(cols[startsWith(cols, "Noise")]))

    row_order <- sort(rownames(layer_corrs))

    layer_corrs <- layer_corrs[row_order, column_order]

    pt <- corrplot(t(as.matrix(layer_corrs)), method="color", is.corr=TRUE, col=COL2("PiYG", 200))
    return(pt)
}

layer_corrs_unsmoothed <- "/Users/jy1008/Downloads/layer_smoothing_V4_seed101_11302023/maynard_validation_pearson_V4_layer_label_v2_all_11302023.csv"
pdf("layer_validation_maynard_V4_unsmoothed_seed101_02012024.pdf")
pt <- generate_corr_heatmap(layer_corrs_unsmoothed)
print(pt)
dev.off()

layer_corrs_smoothed <- "/Users/jy1008/Downloads/layer_smoothing_V4_seed101_11302023/maynard_validation_pearson_V4_smoothed_label_s5_all_11302023.csv"
pdf("layer_validation_maynard_V4_smoothed_seed101_02012024.pdf")
pt <- generate_corr_heatmap(layer_corrs_smoothed)
print(pt)
dev.off()
