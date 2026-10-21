# Post-Fold-Removal Visualization
#
# Generates a three-panel PDF comparing spatial gene count distributions
# before and after fold removal, alongside the high-resolution tissue image.
# Intended to be run after app.R has produced a cleaned Seurat RDS.
#
# Usage:
#   Rscript post_fold_removal_visualization.R \
#     -f <path_to_rds> \
#     -i <path_to_image.png> \
#     -t <sample_title> \
#     -o <output.pdf>
#
# Author: Jie Yuan

library(Seurat)
library(ggplot2)
library(optparse)
library(png)
library(grid)

option_list <- list(
  make_option(c("-f", "--rdsfile"), type = "character", default = NULL,
              help = "Input Seurat RDS with fold spots marked in 'selected_spot' column (required)",
              metavar = "character"),
  make_option(c("-i", "--imgfile"), type = "character", default = NULL,
              help = "Input high-resolution slide image (.png, required)",
              metavar = "character"),
  make_option(c("-t", "--title"), type = "character", default = "test",
              help = "Sample name used as plot title [default: %default]",
              metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = "fold_vis_out.pdf",
              help = "Output PDF filename [default: %default]",
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt        <- parse_args(opt_parser)

if (is.null(opt$rdsfile)) {
  print_help(opt_parser)
  stop("--rdsfile is required.", call. = FALSE)
}
if (is.null(opt$imgfile)) {
  print_help(opt_parser)
  stop("--imgfile is required.", call. = FALSE)
}


set.seed(101)

message("Loading Seurat object from: ", opt$rdsfile)
data <- readRDS(opt$rdsfile)

# Panel 1: nFeature_Spatial for all spots (pre-removal baseline)
p1 <- SpatialFeaturePlot(data, features = "nFeature_Spatial") +
        theme(legend.position = "right") +
        labs(title = opt$title)

# Panel 2: nFeature_Spatial with fold spots excluded (selected_spot == TRUE marks folds)
p2 <- SpatialFeaturePlot(
        subset(data, subset = selected_spot == FALSE),
        features = "nFeature_Spatial"
      ) +
        theme(legend.position = "right")

# Panel 3: High-resolution tissue image for visual reference
img <- readPNG(opt$imgfile)
p3  <- rasterGrob(img, interpolate = TRUE)

# Combine panels side-by-side and write to PDF
message("Saving visualization to: ", opt$output)
ggsave(
  filename = opt$output,
  plot = p1 + p2 + p3,
  width = 24,
  height = 8,
  units = "in"
)
message("Done.")
