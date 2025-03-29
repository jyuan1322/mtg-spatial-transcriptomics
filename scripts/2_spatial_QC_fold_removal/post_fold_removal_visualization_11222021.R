# This script is for visualizing the results of removing spots which overlap a fold in the tissue using the Shiny app.
# Jie Yuan
# 11/22/2021

library(Seurat)
library(ggplot2)
library(optparse)
library(png)
library(grid)

option_list = list(
  make_option(c("-f", "--rdsfile"), type="character", default=NULL,
              help="input Seurat object with fold spots marked", metavar="character"),
  make_option(c("-i", "--imgfile"), type="character", default=NULL,
              help="input high-definition slide image", metavar="character"),
  make_option(c("-t", "--title"), type="character", default="test",
              help="sample name [default= %default]", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="fold_vis_out.pdf",
              help="output pdf file name [default= %default]", metavar="character")
);
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$rdsfile)){
  print_help(opt_parser)
  stop("input rds is required", call.=FALSE)
}
if (is.null(opt$imgfile)){
  print_help(opt_parser)
  stop("input img is required", call.=FALSE)
}

set.seed(101)

data <- readRDS(opt$rdsfile)
p1 <- SpatialFeaturePlot(data, features="nFeature_Spatial") + theme(legend.position="right") + labs(title=opt$title)
p2 <- SpatialFeaturePlot(subset(data, subset=selected_spot == FALSE), features="nFeature_Spatial") + theme(legend.position="right")

img <- readPNG(opt$imgfile)
p3 <- rasterGrob(img, interpolate=TRUE)

ggsave(filename=opt$output, plot=p1 + p2 + p3, width=24, height=8, units="in")
