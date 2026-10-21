# see previous version: sc_seurat_merge_cell2location_props_04102024.R

library(Seurat)
library(testit)
library(stringr)
library(dplyr)
library(tidyr)
library(tibble)

library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(pheatmap)
library(corrplot)
library(umap)

random_seed <- 101

# panda server
# comb_seurat <- readRDS("/mnt/data0/projects/ASAP/ST/clustering_pipeline/layer_smoothing_V4_seed101_11302023_replication_test/comb_seurat_post_smoothed_label_s5_11302023.rds")
# tiger server
comb_seurat <- readRDS("/mnt/data/projects/ASAP/ST/clustering_pipeline/layer_smoothing_V4_seed101_11302023_replication_test/comb_seurat_post_smoothed_label_s5_11302023.rds")

# remove the "_rep1" from some samples
# comb_seurat@meta.data$sample_name <- gsub("_rep[0-9]*", "", comb_seurat@meta.data$orig.ident)
comb_seurat@meta.data$spot_id <- rownames(comb_seurat@meta.data)

# panda server
# cell_prop_dir <- "/mnt/data0/projects/ASAP/ST/st_sc_deconvolution/results_V4_02202024/cell2location_map"
# tiger server
cell_prop_dir <- "/mnt/data/projects/ASAP/ST/st_sc_deconvolution/results_V4_02202024/cell2location_map"

# cell_type_strs <- c("Astro", "Endo", "GABA", "GLU", "Micro", "OPC", "Oligo", "TEMRA")
# cell_type_strs <- paste(cell_type_strs, collapse="|")
# easier to ignore the first N columns
num_non_cell_type_cols <- 35
last_non_cell_type_col <- "X_scvi_labels"

comb_meta_data <- comb_seurat@meta.data
comb_meta_cell2loc <- c()
all_cell_type_cols <- c() # these can be numerous and not known beforehand
for(cellprop_file in list.files(cell_prop_dir, pattern="st_cellprop_table.*\\.csv")) {
    print(cellprop_file)
    cell2loc_sample <- gsub("\\.csv", "", gsub("st_cellprop_table_", "", cellprop_file))
    cellprops <- read.csv(file.path(cell_prop_dir, cellprop_file))
    colnames(cellprops)[colnames(cellprops) == "X"] <- "spot_id"
    
    # these may be slightly different if some cell types are missing
    # cell_type_cols <- colnames(cellprops)[grepl(cell_type_strs, colnames(cellprops))]
    assert(colnames(cellprops)[num_non_cell_type_cols] == last_non_cell_type_col)
    cell_type_cols <- colnames(cellprops)[(num_non_cell_type_cols+1):ncol(cellprops)]
    cellprops <- cellprops[, c("spot_id", cell_type_cols)]
    all_cell_type_cols <- sort(unique(c(all_cell_type_cols, cell_type_cols)))
    
    # comb_sub <- subset(comb_seurat, subset = sample_name == cell2loc_sample)
    # comb_sub@meta.data$Row.names <- rownames(comb_sub@meta.data)
    comb_meta_sub <- comb_meta_data[comb_meta_data$sample_name == cell2loc_sample, ]
    comb_meta_sub <- merge(comb_meta_sub, cellprops, by="spot_id", all.x=TRUE, sort=FALSE)
    
    comb_meta_cell2loc <- dplyr::bind_rows(comb_meta_cell2loc, comb_meta_sub)
    
    # temp <- merge(comb_sub@meta.data, cellprops, by="new_spot_ids")
    # temp <- merge(comb_sub@meta.data, cellprops, by="spot_ids")
    # rownames(temp) <- temp$Row.names
    # temp <- temp[ colnames(comb_sub) ,]
    # comb_sub@meta.data <- temp
}
rownames(comb_meta_cell2loc) <- comb_meta_cell2loc$spot_id
comb_meta_data <- merge(comb_meta_data, comb_meta_cell2loc[, c("spot_id", all_cell_type_cols)], by="spot_id", all.x=TRUE, sort=FALSE)
rownames(comb_meta_data) <- comb_meta_data$spot_id
# re-sort the meta rows to match data columns
comb_meta_data <- comb_meta_data[colnames(comb_seurat), ]
assert(all(colnames(comb_seurat) == rownames(comb_meta_data)))
assert(length(unique(rownames(comb_meta_data))) == nrow(comb_meta_data))

comb_seurat@meta.data <- comb_meta_data

# # print layers to validate merged data tables are correct
# group.colors <- c("Layer 1" = "#f43545", "Layer 2"="#ff8901", "Layer 3"="#fad717", "Layer 4"="#00ba71", "Layer 5"="#00c2de", "Layer 6"="#ba55d3", "WM"="#652175")
# pdf("cell2loc_layer_check_V4_03052024.pdf")
# for(img in names(comb_seurat@images)) {
#     # skip if this image has 0 coordinates (sample is not included)
#     if(nrow(comb_seurat@images[[img]]@coordinates) == 0) {
#         next
#     }
#
#     meta_row <- rownames(comb_seurat@images[[img]]@coordinates)[1]
#     sample_name_filt <- comb_seurat@meta.data[meta_row, "orig.ident"]
#     diagnosis_filt <- comb_seurat@meta.data[meta_row, "diagnosis"]
#
#     sp_pt <- SpatialPlot(comb_seurat, group.by="smoothed_label_s5", image.alpha=0.3, images=img, stroke=NA) + scale_fill_manual(values=group.colors) + theme(legend.title = element_blank())
#     print(sp_pt)
# }
# dev.off()

# in cell type columns, fill in NA values with 0
for(cell_type_col in all_cell_type_cols) {
    comb_seurat@meta.data[is.na(comb_seurat@meta.data[, cell_type_col]), cell_type_col] <- 0
}
# calculate per-spot cell2loc pseudoabundance sum
comb_seurat@meta.data$cell2loc_sum <- rowSums(comb_seurat@meta.data[, all_cell_type_cols])


# grab major cell types
major_cell_types <- c()
sc_meta_file_path <- "/mnt/data/projects/ASAP/ST/data/sc_for_deconvolution/sc_for_deconvolution_V4_02202024"
for(fl in list.files(path=sc_meta_file_path, pattern="*_Metadata.csv", full.names=TRUE, recursive=FALSE)) {
    sc_meta_data <- read.csv(fl)
    print(fl)
    major_cell_types <- rbind(major_cell_types, sc_meta_data[, c("MajorCellTypes", "CellSubtypes", "Complex_Assignment")])
}
major_cell_types$MajorCellTypes <- gsub("_", ".", major_cell_types$MajorCellTypes)
# replace parentheses with ..
major_cell_types$Complex_Assignment <- gsub("\\(|\\)", "..", major_cell_types$Complex_Assignment)
# replace space and all other punctuation with .
major_cell_types$Complex_Assignment <- gsub("[^[:alnum:] ]", ".", gsub(" ", ".", major_cell_types$Complex_Assignment))
# write.csv(table(major_cell_types[, c("CellSubtypes", "Complex_Assignment")]), "cell_subtypes_to_complex_assignment_03132024.csv")

# filter out very small clusters (already removed in all_cell_type_cols)
major_cell_types <- major_cell_types[major_cell_types$Complex_Assignment %in% all_cell_type_cols, ]
major_cell_types_unq <- unique(major_cell_types)
major_cell_type_cols <- sort(unique(major_cell_types_unq$MajorCellTypes))
for(ct in major_cell_type_cols) {
    ct_colnames <- unique(major_cell_types_unq[major_cell_types_unq$MajorCellTypes == ct, "Complex_Assignment"])
    # ct_colnames <- gsub("[^[:alnum:] ]", ".", gsub(" ", ".", ct_colnames))
    ct_colnames <- intersect(ct_colnames, all_cell_type_cols)
    if(length(ct_colnames) > 1) {
        comb_seurat@meta.data[, ct] <- rowSums(comb_seurat@meta.data[, ct_colnames])
    } else {
        assert(gsub("_", ".", ct) %in% colnames(comb_seurat@meta.data))
    }
}
# major_cell_type_cols contains only major cell types
# all_cell_type_cols contains all cell types but non-redundant, e.g. "GLU.Neurons" is not counted again


// # create grouped complex assignments based on shared first gene
// complex_group_cols <- sort(unique(str_split_fixed(all_cell_type_cols, pattern="\\.", n=2)[,1]))
// complex_group_cols[complex_group_cols == "T"] <- "T.Cells"
// for(cg in complex_group_cols) {
//     cg_colnames <- all_cell_type_cols[startsWith(all_cell_type_cols, cg)]
//     if(length(cg_colnames) > 1) {
//         comb_seurat@meta.data[, cg] <- rowSums(comb_seurat@meta.data[, cg_colnames])
//     } else {
//         # keep the original column name, but be sure to update the list of complex group columns
//         complex_group_cols[complex_group_cols == cg] <- cg_colnames[1]
//     }
// }



# Cell type pseudoabundances mapped to spatial coordinates
custom_colors=c("#ffffff", "#ff0000", "#e60000", "#cd0000")
single_cell_type <- "RORB.ESR1"
pdf(sprintf("st_sc_cell2location_mapped_normed_1clr_%s_04102024.pdf", single_cell_type))
for(img in names(comb_seurat@images)) {
    # skip if this image has 0 coordinates (sample is not included)
    if(nrow(comb_seurat@images[[img]]@coordinates) == 0) {
        next
    }
    
    meta_row <- rownames(comb_seurat@images[[img]]@coordinates)[1]
    sample_name_filt <- comb_seurat@meta.data[meta_row, "sample_name"]
    diagnosis_filt <- comb_seurat@meta.data[meta_row, "diagnosis"]

    
    # don't print cell types which are empty
    # nonzero_cell_types <- comb_seurat@meta.data %>%
    #                         filter(sample_name == as.character(sample_name_filt)) %>%
    #                         dplyr::select(all_cell_type_cols) %>%
    #                         colSums()
    # nonzero_cell_types <- names(nonzero_cell_types[nonzero_cell_types > 0])
    # if(length(nonzero_cell_types) == 0) {
    #     print(paste(sample_name_filt, "Skipped: no cell2loc data"))
    #     next
    # }
    
    # pdf(sprintf("st_sc_cell2location_mapped_%s_03052024.pdf", sample_name_filt))
    # for(cell_type_col in nonzero_cell_types) {
    #     print(paste(sample_name_filt, cell_type_col))
    #     sp_pt <- SpatialFeaturePlot(comb_seurat, features=cell_type_col, images=img)
    #     print(sp_pt)
    # }
    # dev.off()
    
    # print another version normalized by sum
    # pdf(sprintf("st_sc_cell2location_mapped_normed_%s_03052024.pdf", sample_name_filt))
    # for(cell_type_col in nonzero_cell_types) {
    #     print(paste(sample_name_filt, cell_type_col, "normed"))
    #     normed_name <- paste0(cell_type_col, "_normed")
    #     comb_seurat@meta.data[, normed_name] <- comb_seurat@meta.data[, cell_type_col] / comb_seurat@meta.data$cell2loc_sum
    #     sp_pt <- SpatialFeaturePlot(comb_seurat, features=normed_name, images=img)
    #     print(sp_pt)
    #     comb_seurat@meta.data[, normed_name] <- NULL
    # }
    # dev.off()
    
    # pdf(sprintf("st_sc_cell2location_mapped_normed_1clr_%s_04102024.pdf", sample_name_filt))
    # for(cell_type_col in nonzero_cell_types) {
    cell_type_col <- single_cell_type
    print(paste(sample_name_filt, cell_type_col, "normed"))
    normed_name <- paste0(cell_type_col, "_normed")
    comb_seurat@meta.data[, normed_name] <- comb_seurat@meta.data[, cell_type_col] / comb_seurat@meta.data$cell2loc_sum
    max_val <- ceiling(max(comb_seurat@meta.data[, normed_name]))
    sp_pt <- SpatialFeaturePlot(comb_seurat, features=normed_name, images=img) +
             ggplot2::scale_fill_gradientn(colours=custom_colors,
                limits=c(0,max_val)) + ggtitle(paste(sample_name_filt, diagnosis_filt))
    print(sp_pt)
    comb_seurat@meta.data[, normed_name] <- NULL
    # }
    # dev.off()
    
    # print image only if not empty
    # for( img_name in Images(comb_sub) ) {
    #     if( nrow(comb_sub@images[[img_name]]@coordinates) > 0 ) {
    #         pdf(sprintf("st_sc_cell2location_mapped_%s_02252023.pdf", cell2loc_sample))
    #         for (cell_type_col in cell_type_cols) {
    #             sp_pt <- SpatialFeaturePlot(comb_sub, features=cell_type_col, images=img_name)
    #             print(sp_pt)
    #         }
    #         dev.off()
    #     }
    # }
}
dev.off()


gaba_cell_types <- all_cell_type_cols[grepl("^ADARB2|^GAD1|LAMP5\\.CA1|^LAMP5\\.LCP2|LAMP5\\.NMBR|^PAX6|^PVALB|^SST|^VIP", all_cell_type_cols)]
glu_cell_types <- all_cell_type_cols[grepl("^FEZF2|LAMP5\\.GREB1L|LAMP5\\.PARD3|^LINC00507|^RORB|^SLC17A7|^THEMIS", all_cell_type_cols)]

# https://stackoverflow.com/a/60533772
entropy <- function(p) rowSums(-(p * log(p)), na.rm = TRUE)

cell_types_mat <- major_cell_types
calculate_celltype_entropy <- function(comb_seurat, all_cell_type_cols, major_cell_type_cols, cell_types_mat) {
    
    all_cols <- union(all_cell_type_cols, major_cell_type_cols)
    cell_props <- comb_seurat@meta.data[, c("sample_name", "diagnosis", "smoothed_label_s5", "cell2loc_sum", all_cols)]
    # remove samples which have 0 cell2loc info
    temp <- cell_props %>%
        group_by(sample_name) %>%
        summarise(cell2loc_sample_sum = sum(cell2loc_sum)) %>%
        as.data.frame()
    include_samples <- as.character(temp[temp$cell2loc_sample_sum > 0, "sample_name"])
    rm(temp)
    print(paste("Num samples", length(unique(cell_props$sample_name))))
    cell_props <- cell_props[cell_props$sample_name %in% include_samples, ]
    print(paste("Num samples after filtering", length(unique(cell_props$sample_name))))

    # include sc cell count info
    cell_types_mat$spot_id <- rownames(cell_types_mat)
    cell_types_mat$sample_id <- str_split_fixed(cell_types_mat$spot_id, "_", n=3)[, 2]
    cell_types_mat_lnsum <- cell_types_mat %>%
                            filter(sample_id %in% include_samples) %>%
                            filter(Complex_Assignment %in% all_cols) %>%
                            group_by(sample_id, Complex_Assignment) %>%
                            summarise(ct_count = log10(n()+1)) %>%
                            as.data.frame()
    # cell_types_mat_ave <- cell_types_mat_lnsum %>%
    #                         group_by(Complex_Assignment) %>%
    #                         summarise(ct_mean = mean(ct_count), ct_sd = sd(ct_count)) %>%
    #                         as.data.frame()
    cell_types_mat_maj <- cell_types_mat %>%
                            filter(sample_id %in% include_samples) %>%
                            filter(MajorCellTypes %in% all_cols) %>%
                            group_by(sample_id, MajorCellTypes) %>%
                            summarise(ct_count = log10(n()+1)) %>%
                            as.data.frame()
    colnames(cell_types_mat_lnsum)[colnames(cell_types_mat_lnsum) == "Complex_Assignment"] <- "cell_type"
    colnames(cell_types_mat_maj)[colnames(cell_types_mat_maj) == "MajorCellTypes"] <- "cell_type"
    cell_types_mat_comb <- rbind(cell_types_mat_lnsum, cell_types_mat_maj)
    cell_types_mat_comb <- unique(cell_types_mat_comb) # remove duplicate major cell types (e.g. Astrocytes)
    
    # normalize by spot
    # in each spot (row), divide pseudocount by spot count (row sum)
    # in spot sets under each sample and layer, calculate the mean
    # take the mean of the means across samples for each layer
    cell_props_spotnorm <- cell_props
    cell_props_spotnorm[, all_cell_type_cols] <- sweep(
            cell_props_spotnorm[, all_cell_type_cols], 1,
            rowSums(cell_props_spotnorm[, all_cell_type_cols]), `/`)
    # For GLU.Neurons, GABA.Neurons, and Pericytes, get sum over sub-types.
    # Don't double-count these clusters when normalizing.
    cell_props_spotnorm_maj <- cell_props
    cell_props_spotnorm_maj[, major_cell_type_cols] <- sweep(
            cell_props_spotnorm_maj[, major_cell_type_cols], 1,
            rowSums(cell_props_spotnorm_maj[, major_cell_type_cols]), `/`)
    add_maj_clusters <- setdiff(major_cell_type_cols, all_cell_type_cols)
    cell_props_spotnorm[, add_maj_clusters] <- cell_props_spotnorm_maj[, add_maj_clusters]

    cell_props_spotnorm_psamp_player <- cell_props_spotnorm %>%
            group_by(sample_name, smoothed_label_s5, diagnosis) %>%
            summarise_at(all_cols, mean, na.rm = FALSE) %>%
            as.data.frame()

    # mean over all indivs
    cell_props_spotnorm_all <- cell_props_spotnorm_psamp_player %>%
            group_by(smoothed_label_s5) %>%
            summarize_at(all_cols, mean, na.rm = FALSE) %>%
            as.data.frame()
    rownames(cell_props_spotnorm_all) <- cell_props_spotnorm_all$smoothed_label_s5
    cell_props_spotnorm_all$smoothed_label_s5 <- NULL
    column_order <- t(cell_props_spotnorm_all) %*% c(-100, -10, -1, 0, 1, 10, 100) %>% as.data.frame()
    column_order$group <- rownames(column_order)
    column_order <- column_order[order(column_order$V1), "group"]
    cell_props_spotnorm_all <- cell_props_spotnorm_all[, column_order]

    cell_props_spotnorm_all_rowsumnormed <- sweep(
            cell_props_spotnorm_all, 2,
            colSums(cell_props_spotnorm_all), `/`)

    # for each cell type, get the max mean percentage composition, and the layer that max occurs in
    cpmat <- as.data.frame(t(cell_props_spotnorm_all))
    temp <- cpmat %>%
             rowwise() %>%
             mutate(max_meanpct=max(c_across(colnames(cpmat))),
                    max_layer=colnames(cpmat)[
                                  which.max(c_across(colnames(cpmat)))
                              ]) %>%
             as.data.frame()
    rownames(temp) <- rownames(cpmat)
    cpmat <- temp
    rm(temp)

    # add major cell type group
    cpmat$major_cell_type <- NA
    cpmat[rownames(cpmat) %in% glu_cell_types, "major_cell_type"] <- "GLU"
    cpmat[rownames(cpmat) %in% gaba_cell_types, "major_cell_type"] <- "GABA"
    cpmat[is.na(cpmat$major_cell_type), "major_cell_type"] <- rownames(cpmat[is.na(cpmat$major_cell_type), ])

    # entropy calculation
    cpmat_normed <- t(cell_props_spotnorm_all_rowsumnormed)
    # cpmat has 70 rows (cell types) and 7 columns (layers)
    # each row (cell type) sums to 1
    entropy_vals <- entropy(cpmat_normed)
    entropy_vals <- data.frame('cell_type' = names(entropy_vals),
                               'entropy' = entropy_vals)
    cpmat$cell_type <- rownames(cpmat)
    cpmat <- merge(cpmat, entropy_vals, by="cell_type")
    cpmat$vis_cell_type <- cpmat$major_cell_type
    cpmat[cpmat$major_cell_type %in% c("GLU", "GABA"), "vis_cell_type"] <-
        sprintf("[%s] %s",
            cpmat[cpmat$major_cell_type %in% c("GLU", "GABA"), "major_cell_type"],
            cpmat[cpmat$major_cell_type %in% c("GLU", "GABA"), "cell_type"])
    
    # sort cpmat according to layer, and then increasing entropy
    # flip the entropy scale so lower entropy is on the right
    cpmat <- cpmat %>% arrange(major_cell_type, max_layer, -entropy)
    # fix the order
    cpmat$vis_cell_type <- as.character(cpmat$vis_cell_type)
    y_label_order <- c(
            cpmat$vis_cell_type[startsWith(cpmat$vis_cell_type, "[GABA]")],
            cpmat$vis_cell_type[startsWith(cpmat$vis_cell_type, "[GLU]")],
            "Pericytes.2", "Pericytes.1", "OPCs", "Pericytes", "Endothelial.Cells", "Astrocytes", "Microglia",
            "GABA.Neurons", "GLU.Neurons", "T.Cells", "Fibroblast.Like.Cells", "Oligodendrocytes")
    cpmat$vis_cell_type <- factor(cpmat$vis_cell_type, levels=y_label_order)
    write.csv(cpmat, "cell2loc_glu_gaba_entropy_calcs_08272024.csv")
    sp_pt <- ggplot(cpmat, aes(x=entropy, y=vis_cell_type,
                                color=max_layer, size=max_meanpct)) +
        geom_point() +
        scale_size(range = c(1,6), limits = c(0, 0.6)) +
        scale_x_reverse() +
        theme_bw()
    
    
    # actual cell count from single-cell data, mean across samples
    # rownames(cell_types_mat_comb) <- cell_types_mat_comb$cell_type
    # cell_types_mat_comb <- cell_types_mat_comb[cpmat$cell_type, ]
    cell_types_mat_comb$vis_cell_type <- cell_types_mat_comb$cell_type
    cell_types_mat_comb[cell_types_mat_comb$vis_cell_type %in% glu_cell_types, "vis_cell_type"] <-
        sprintf("[GLU] %s", cell_types_mat_comb[cell_types_mat_comb$vis_cell_type %in% glu_cell_types, "vis_cell_type"])
    cell_types_mat_comb[cell_types_mat_comb$vis_cell_type %in% gaba_cell_types, "vis_cell_type"] <-
        sprintf("[GABA] %s", cell_types_mat_comb[cell_types_mat_comb$vis_cell_type %in% gaba_cell_types, "vis_cell_type"])
    sc_bp_pt <- ggplot(cell_types_mat_comb,
                       aes(x=ct_count, y=factor(vis_cell_type, levels=y_label_order))) +
                geom_boxplot(outlier.size=1) +
                annotation_logticks(sides="tb") +
                theme_bw()
    # insert zoomed image for GABA
    # cell_types_mat_temp <- cell_types_mat_comb
    # cell_types_mat_temp[!(cell_types_mat_temp$cell_type %in% gaba_cell_types), "ct_count"] <- 0
    # sc_bp_pt_gaba <- ggplot(cell_types_mat_temp, aes(x=ct_count, y=factor(vis_cell_type, levels=y_label_order))) + geom_boxplot(outlier.size=1) + theme_bw()
    # cell_types_mat_temp <- cell_types_mat_comb
    # cell_types_mat_temp[!(cell_types_mat_temp$cell_type %in% glu_cell_types), "ct_count"] <- 0
    # sc_bp_pt_glu <- ggplot(cell_types_mat_temp, aes(x=ct_count, y=factor(vis_cell_type, levels=y_label_order))) + geom_boxplot(outlier.size=1) + theme_bw()
    # rm(cell_types_mat_temp)

    # boxplot of sum (per sample) pseudoabundance
    # spatial_sum_persamp <- cell_props %>%
    #                         select(-c(diagnosis, smoothed_label_s5, cell2loc_sum)) %>%
    #                         pivot_longer(!sample_name, names_to="cell_type", values_to="proportion") %>%
    #                         group_by(sample_name, cell_type) %>%
    #                         summarise(sum_psct = log10(sum(proportion)+1)) %>%
    #                         as.data.frame()
    # spatial_sum_persamp$vis_cell_type <- spatial_sum_persamp$cell_type
    # spatial_sum_persamp[spatial_sum_persamp$vis_cell_type %in% glu_cell_types, "vis_cell_type"] <-
    #     sprintf("[GLU] %s", spatial_sum_persamp[spatial_sum_persamp$vis_cell_type %in% glu_cell_types, "vis_cell_type"])
    #     spatial_sum_persamp[spatial_sum_persamp$vis_cell_type %in% gaba_cell_types, "vis_cell_type"] <-
    #     sprintf("[GABA] %s", spatial_sum_persamp[spatial_sum_persamp$vis_cell_type %in% gaba_cell_types, "vis_cell_type"])
    # st_bp_pt <- ggplot(spatial_sum_persamp,
    #                     aes(x=sum_psct, y=factor(vis_cell_type, levels=y_label_order))) +
    #             geom_boxplot(outlier.size=1) + theme_bw()
    
    pdf("cell2loc_glu_gaba_entropy_calcs_09022024.pdf", width=18, height=10)
    # all_pt <- cowplot::plot_grid(sc_bp_pt, st_bp_pt, sp_pt,
    #           nrow = 1, rel_widths = c(0.6, 0.6, 1), align = 'h', axis = 'bt')
    # all_pt <- cowplot::plot_grid(sc_bp_pt, sc_bp_pt_glu, sc_bp_pt_gaba, sp_pt,
    #           nrow = 1, rel_widths = c(0.6, 0.6, 0.6, 1.0), align = 'h', axis = 'bt')
    all_pt <- cowplot::plot_grid(sc_bp_pt, sp_pt,
              nrow = 1, rel_widths = c(0.6, 1), align = 'h', axis = 'bt')
    print(all_pt)
    dev.off()

    # subpanel for glu/gaba comparison
    pdf("cell2loc_glu_vs_gaba_subpanel_09092024.pdf", width=2.5, height=2.5)
    # "GLU" is the mean entropy of subclusters, "Glu.Neurons" is the combined entropy
    cpmat_temp <- cpmat
    cpmat_temp$major_cell_type <- str_replace_all(cpmat_temp$major_cell_type,
                                        c("Pericytes.1" = "Peri",
                                          "Pericytes.2" = "Peri"))
    cpmat_temp$major_cell_type <- str_replace_all(cpmat_temp$major_cell_type,
                                        c("Pericytes" = "Pericyte.Neurons")) # do this separately so it doesn't replace the above
    subpanel_disp_order <- c("GLU", "GLU.Neurons", "GABA", "GABA.Neurons", "Peri", "Pericyte.Neurons", "Oligodendrocytes", "Fibroblast.Like.Cells", "T.Cells", "Microglia", "Astrocytes", "Endothelial.Cells", "OPCs")
    set.seed(random_seed)
    en_pt <- ggplot(cpmat_temp, aes(x=entropy,
                                    y=factor(major_cell_type, levels=rev(subpanel_disp_order)))) +
                geom_boxplot(outlier.color=NA) +
                geom_jitter() +
                scale_x_reverse() +
                theme_bw()
    print(en_pt)
    dev.off()
    
    # get lowest N GLU and GABA entropy scores (most layer specificity)
    heatmap_disp_order <- list(
                "Major" = c("Oligodendrocytes", "Fibroblast.Like.Cells", "T.Cells", "GLU.Neurons", "GABA.Neurons", "Microglia", "Astrocytes", "Endothelial.Cells", "Pericytes", "OPCs"),
                "GLU" = c("LINC00507.FREM3.LINC01090", "LAMP5.GREB1L", "LAMP5.PARD3", "RORB.ADARB2", "LINC00507.RORB.SV2C", "RORB.TTC12", "RORB.ESR1", "FEZF2.IL26", "RORB.GRIN3A", "RORB.SEMA3E", "SLC17A7.IL15..FEZF2.SCUBE1", "FEZF2.ABO", "THEMIS.HS3ST5"),
                "GABA" = c("LAMP5.NMBR", "PVALB.CALM1", "PVALB.LGR5", "PVALB.SCUBE3...Chandelier..", "PVALB.WFDC2.ROBO2", "SST.HPGD", "PVALB.SULF1", "SST.ZFPM2.AS1", "PVALB.MEPE", "PVALB.WFDC2.CTXND1", "SST.KLHDC8")
                )
    pdf("cell2loc_glu_gaba_heatmaps_09092024.pdf", width=10, height=10)
    for(maj_ct in c("Major", "GLU", "GABA")) {
        # topN <- 10
        # subset table, and specify order
        # topN_mat <- cpmat %>% filter(major_cell_type == maj_ct) %>%
        #                 slice(which(rank(entropy) %in% 1:topN))
        # topN_mat <- topN_mat %>% arrange(max_layer, -entropy)
        assert(nrow(cpmat) == length(unique(cpmat$cell_type)))
        cpmat_temp <- cpmat
        rownames(cpmat_temp) <- cpmat_temp$cell_type
        topN_mat <- cpmat_temp[heatmap_disp_order[[maj_ct]], ]
        rownames(topN_mat) <- topN_mat$vis_cell_type
        topN_mat_sub <- topN_mat[, c(paste("Layer", 1:6), "WM")]
        # normalize each row to sum to 1
        topN_mat_norm <- sweep(topN_mat_sub,1,rowSums(topN_mat_sub),`/`)
        # shared legend scale
        breaks = seq(0.0, 0.5, 0.01)
        pheatmap(topN_mat_norm,
                    color=colorRampPalette(rev(brewer.pal(11, "Spectral")))(length(breaks)),
                    breaks = breaks,
                    cluster_rows = FALSE,
                    cluster_cols = FALSE,
                    display_numbers = round(topN_mat_norm, 4))
    }
    dev.off()
}
calculate_entropy(comb_seurat, all_cell_type_cols=all_cell_type_cols, major_cell_type_cols=major_cell_type_cols)

# add boxplot of total amount of GLU and GABA expression


# generate proportion charts for cell types across layers
# NOTE: NA values have already been replaced with 0 above
make_grid_spotnorm <- function(comb_seurat, cols, outname, width=10, height=NA) {
    cell_props <- comb_seurat@meta.data[, c("sample_name", "diagnosis", "smoothed_label_s5", "cell2loc_sum", cols)]
    # remove samples which have 0 cell2loc info
    temp <- cell_props %>%
        group_by(sample_name) %>%
        summarise(cell2loc_sample_sum = sum(cell2loc_sum)) %>%
        as.data.frame()
    include_samples <- as.character(temp[temp$cell2loc_sample_sum > 0, "sample_name"])
    print(paste("Num samples", length(unique(cell_props$sample_name))))
    cell_props <- cell_props[cell_props$sample_name %in% include_samples, ]
    print(paste("Num samples after filtering", length(unique(cell_props$sample_name))))

    # normalize by spot
    # in each spot (row), divide pseudocount by spot count (row sum)
    # in spot sets under each sample and layer, calculate the mean
    # take the mean of the means across samples for each layer
    cell_props_spotnorm <- cell_props
    cell_props_spotnorm[, 5:ncol(cell_props_spotnorm)] <- sweep(
            cell_props_spotnorm[, 5:ncol(cell_props_spotnorm)], 1,
            rowSums(cell_props_spotnorm[, 5:ncol(cell_props_spotnorm)]), `/`)
    cell_props_spotnorm <- cell_props_spotnorm %>%
            group_by(sample_name, smoothed_label_s5, diagnosis) %>%
            summarise_at(cols, mean, na.rm = FALSE) %>%
            as.data.frame()

    # mean over all indivs
    cell_props_spotnorm_all <- cell_props_spotnorm %>%
            group_by(smoothed_label_s5) %>%
            summarize_at(cols, mean, na.rm = FALSE) %>%
            as.data.frame()
    rownames(cell_props_spotnorm_all) <- cell_props_spotnorm_all$smoothed_label_s5
    cell_props_spotnorm_all$smoothed_label_s5 <- NULL
    column_order <- t(cell_props_spotnorm_all) %*% c(-100, -10, -1, 0, 1, 10, 100) %>% as.data.frame()
    column_order$group <- rownames(column_order)
    column_order <- column_order[order(column_order$V1), "group"]
    cell_props_spotnorm_all <- cell_props_spotnorm_all[, column_order]

    # # rearrange
    # distance.cell<- dist(t(cell_props_spotnorm_all), method = 'euclidean') # the default
    # cluster.cell <- hclust(distance.cell, method = 'complete') # the default
    # # order_score <- rowSums(pathway_table_wide * c(10^6, 10^5, 10^4, 10^3, 10^2, 10, 1))
    # column_order <- t(cell_props_spotnorm_all) %*% c(-100, -10, -1, 0, 1, 10, 100) %>% as.data.frame()
    # column_order <- column_order %>%
    #                     rownames_to_column() %>% deframe()
    # reordered.cluster.cell <- reorder(as.dendrogram(cluster.cell), column_order)
    # cell_display_order <- colnames(cell_props_spotnorm_all)[order.dendrogram(reordered.cluster.cell)]
    # cell_props_spotnorm_all <- cell_props_spotnorm_all[, cell_display_order]
    
    # reorder row dendrogram in pheatmap
    # https://www.rdocumentation.org/packages/pheatmap/versions/1.0.12/topics/pheatmap
    callback = function(hc, mat){
        # sv = svd(t(mat))$v[,1]
        # sv = mat %*% 10^c(1:ncol(mat))
        mat = apply(mat, 2, function(x) x==max(x, na.rm = TRUE))
        sv = mat %*% 10^c(ncol(mat):1)
        dend = reorder(as.dendrogram(hc), wts = sv)
        as.hclust(dend)
    }

    # largest_num_temp <- ceiling(max(cell_props_spotnorm_all) * 10) / 10
    if(is.na(height)) {
        height <- max(10, nrow(t(cell_props_spotnorm_all)) / 5)
    }
    pdf(sprintf("st_sc_corrplot_spotnorm_V4_%s.pdf", outname), width=10, height=20)
    pheatmap(t(cell_props_spotnorm_all),
                color=colorRampPalette(rev(brewer.pal(11, "Spectral")))(101),
                cluster_distance_rows = "euclidean",
                cluster_cols = FALSE,
                display_numbers = TRUE,
                clustering_callback = callback)
    # corrplot(t(as.matrix(cell_props_spotnorm_all)), tl.col="black", method="color", is.corr=FALSE, col.lim = c(0,largest_num_temp), tl.srt=45, cl.ratio = 0.25,
    #     col=colorRampPalette(rev(brewer.pal(11, "Spectral")))(101))
    # plot(reordered.cluster.cell, horiz=TRUE)
    cell_props_spotnorm_all_rowsumnormed <- sweep(
            cell_props_spotnorm_all, 2,
            colSums(cell_props_spotnorm_all), `/`)
    
    pheatmap(t(cell_props_spotnorm_all_rowsumnormed),
                color=colorRampPalette(rev(brewer.pal(11, "Spectral")))(101),
                cluster_distance_rows = "euclidean",
                cluster_cols = FALSE,
                display_numbers = TRUE,
                clustering_callback = callback)
    
    # keep only cols (cell types) whose max value is greater than a threshold
    ct_filter <- 0.3
    max_cts <- apply(cell_props_spotnorm_all_rowsumnormed, 2, function(x) max(x, na.rm = TRUE))
    cell_props_spotnorm_all_rowsumnormed_filt <- cell_props_spotnorm_all_rowsumnormed[, max_cts > ct_filter]
    pheatmap(t(cell_props_spotnorm_all_rowsumnormed_filt),
                color=colorRampPalette(rev(brewer.pal(11, "Spectral")))(101),
                cluster_distance_rows = "euclidean",
                cluster_cols = FALSE,
                display_numbers = TRUE,
                clustering_callback = callback)
    dev.off()
}
make_grid_spotnorm(comb_seurat, cols=all_cell_type_cols, outname="complex_all_03142024")
make_grid_spotnorm(comb_seurat, cols=major_cell_type_cols, outname="major_only_03142024")
make_grid_spotnorm(comb_seurat, cols=complex_group_cols, outname="complex_group_03142024")



# normalize by cell type
# for each sample and layer, get pseudocount sum
# for each sample, divide pseudocount sum in each layer by sum across layers
# across all samples, take the mean for each layer; only include samples in the mean which have nonzero expression
make_grid_cellnorm <- function(comb_seurat, cols, outname, width=10, height=NA, normalize_per_sample=FALSE) {
    cell_props <- comb_seurat@meta.data[, c("sample_name", "diagnosis", "smoothed_label_s5", "cell2loc_sum", cols)]
    # remove samples which have 0 cell2loc info
    temp <- cell_props %>%
        group_by(sample_name) %>%
        summarise(cell2loc_sample_sum = sum(cell2loc_sum)) %>%
        as.data.frame()
    include_samples <- as.character(temp[temp$cell2loc_sample_sum > 0, "sample_name"])
    print(paste("Num samples", length(unique(cell_props$sample_name))))
    cell_props <- cell_props[cell_props$sample_name %in% include_samples, ]
    print(paste("Num samples after filtering", length(unique(cell_props$sample_name))))
    
    if(normalize_per_sample) {
        cell_props_cellnorm <- cell_props
        cell_props_cellnorm <- cell_props_cellnorm %>%
            group_by(sample_name, smoothed_label_s5, diagnosis) %>%
            summarise_at(cols, sum, na.rm = FALSE) %>%
            as.data.frame()
        cell_props_cellnorm <- cell_props_cellnorm %>%
            group_by(sample_name, diagnosis) %>%
            mutate_at(cols, function(x) x/sum(x)) %>%
            mutate(across(everything(), ~replace(.x, is.nan(.x), 0))) %>%
            as.data.frame()

        # mean over all individuals
        cell_props_cellnorm_all <- cell_props_cellnorm %>%
            group_by(smoothed_label_s5) %>%
            summarise_at(cols, mean, na.rm = FALSE) %>%
            as.data.frame()
        rownames(cell_props_cellnorm_all) <- cell_props_cellnorm_all$smoothed_label_s5
        cell_props_cellnorm_all$smoothed_label_s5 <- NULL
    } else { # get pseudocount across layers, divide by column sum
        cell_props_cellnorm <- cell_props
        cell_props_cellnorm <- cell_props_cellnorm %>%
            group_by(smoothed_label_s5) %>%
            summarise_at(cols, sum, na.rm = FALSE) %>%
            as.data.frame()
        rownames(cell_props_cellnorm) <- cell_props_cellnorm$smoothed_label_s5
        cell_props_cellnorm$smoothed_label_s5 <- NULL
        cell_props_cellnorm_all <- sweep(
                cell_props_cellnorm, 2,
                colSums(cell_props_cellnorm), `/`)
    }
    
    if(is.na(height)) {
        height <- max(10, nrow(t(cell_props_cellnorm_all)) / 5)
    }
    largest_num_temp <- ceiling(max(cell_props_cellnorm_all) * 10) / 10
    pdf(sprintf("st_sc_corrplot_cellnorm_V4_%s.pdf", outname), width=10, height=height)
    pheatmap(t(cell_props_cellnorm_all),
                color=colorRampPalette(rev(brewer.pal(11, "Spectral")))(101),
                cluster_distance_rows = "canberra",
                cluster_cols = FALSE,
                display_numbers = TRUE)
    # corrplot(t(as.matrix(cell_props_cellnorm_all)), tl.col="black", method="color", is.corr=FALSE, col.lim = c(0,largest_num_temp), tl.srt=45, cl.ratio = 0.25,
    #     col=colorRampPalette(c("white","#31a354"))(101))
    # plot(reordered.cluster.cell, horiz=TRUE)
    dev.off()
}
make_grid_cellnorm(comb_seurat, cols=all_cell_type_cols, outname="complex_all_normpersampFALSE_03142024", width=10, height=NA, normalize_per_sample=FALSE)
make_grid_cellnorm(comb_seurat, cols=major_cell_type_cols, outname="major_only_normpersampFALSE_03142024", width=10, height=NA, normalize_per_sample=FALSE)
make_grid_cellnorm(comb_seurat, cols=complex_group_cols, outname="complex_group_normpersampFALSE_03142024", width=10, height=NA, normalize_per_sample=FALSE)



# -----
# Cell type composition vs PD progression
# -----

# TODO:
# - change in cell type proportion in each layer across diagnosis and across LB score

# in nctx_temporal, there is only 1 sample with a score of 4: set this to 3 so as not to overly bias the regression
comb_seurat@meta.data$nctx_temporal_fix <- comb_seurat@meta.data$nctx_temporal
comb_seurat@meta.data$nctx_temporal_fix[comb_seurat@meta.data$nctx_temporal_fix > 3] <- 3

cell_props <- comb_seurat@meta.data[, c("sample_name", "diagnosis", "brain_stem_sn", "nctx_temporal_fix", "smoothed_label_s5", "cell2loc_sum", "batch", "sex", "expired_age_num", "PMI_num", "RIN_num", all_cell_type_cols)]
# remove samples which have 0 cell2loc info
temp <- cell_props %>%
    group_by(sample_name) %>%
    summarise(cell2loc_sample_sum = sum(cell2loc_sum)) %>%
    as.data.frame()
include_samples <- as.character(temp[temp$cell2loc_sample_sum > 0, "sample_name"])
print(paste("Num samples", length(unique(cell_props$sample_name))))
cell_props <- cell_props[cell_props$sample_name %in% include_samples, ]
print(paste("Num samples after filtering", length(unique(cell_props$sample_name))))

# normalize by spot
# in each spot (row), divide pseudocount by spot count (row sum)
# in spot sets under each sample and layer, calculate the mean
# take the mean of the means across samples for each layer
cell_props_spotnorm <- cell_props
num_nonnumeric_cols <- 12
assert(colnames(cell_props_spotnorm)[num_nonnumeric_cols-1] == "RIN_num")
cell_props_spotnorm[, num_nonnumeric_cols:ncol(cell_props_spotnorm)] <- sweep(
        cell_props_spotnorm[, num_nonnumeric_cols:ncol(cell_props_spotnorm)], 1,
        rowSums(cell_props_spotnorm[, num_nonnumeric_cols:ncol(cell_props_spotnorm)]), `/`)
cell_props_spotnorm <- cell_props_spotnorm %>%
        group_by(sample_name, smoothed_label_s5, diagnosis, brain_stem_sn, nctx_temporal_fix, batch, sex, expired_age_num, PMI_num, RIN_num) %>%
        summarise_at(all_cell_type_cols, mean, na.rm = FALSE) %>%
        as.data.frame()

# outname_base <- "03142024"
outname_base <- "04102024"
layers <- sort(unique(cell_props_spotnorm$smoothed_label_s5))
for(ct in all_cell_type_cols) {
    
    for(clinical_feat_discrt in c("diagnosis", "brain_stem_sn", "nctx_temporal_fix")) {
        outname <- paste(ct, clinical_feat_discrt, outname_base, sep = "_")
        pdf(sprintf("st_sc_scatterplot_spotnorm_V4_%s.pdf", outname))
        for(layer in layers) {
            cell_props_sub <- cell_props_spotnorm[cell_props_spotnorm$smoothed_label_s5 == layer, ]
            cell_props_sub[, clinical_feat_discrt] <- factor(cell_props_sub[, clinical_feat_discrt])
            cell_props_sub <- cell_props_sub[!is.na(cell_props_sub[, clinical_feat_discrt]), ]
            cell_props_sub$diagnosis <- factor(cell_props_sub$diagnosis, levels = c("Control", "ILBD", "Case"))
            pt <- ggplot(cell_props_sub, aes_string(x = clinical_feat_discrt, y = ct, fill = clinical_feat_discrt)) + geom_boxplot() + geom_jitter() + ggtitle(paste(ct, layer, clinical_feat_discrt)) + theme_bw()
            print(pt)
        }
        dev.off()
    }
}






# regression
# control ILBD
# lnres_cont_ilbd <- c()
# for(ct in all_cell_type_cols) {
#     for(layer in layers) {
#         cell_props_sub <- cell_props_spotnorm[cell_props_spotnorm$smoothed_label_s5 == layer, ]
#         cell_props_sub <- cell_props_sub[cell_props_sub$diagnosis %in% c("Control", "ILBD"), ]
#         cell_props_sub$diagnosis <- factor(cell_props_sub$diagnosis, levels=c("Control", "ILBD"))
#         cell_props_sub$batch <- factor(cell_props_sub$batch)
#         cell_props_sub$sex <- factor(cell_props_sub$sex)
#         cell_props_sub$ct_to_test <- cell_props_sub[, ct]
#
#         mdl <- lm(data = cell_props_sub,
#                 ct_to_test ~ batch + sex + expired_age_num + PMI_num + RIN_num + diagnosis)
#         res <- summary(mdl)$coefficients["diagnosisILBD", ]
#         lnres_cont_ilbd <- rbind(lnres_cont_ilbd, data.frame(
#                                         "clinical_feature" = "Control-ILBD",
#                                         "cell_type" = ct,
#                                         "layer" = layer,
#                                         "estimate" = res[["Estimate"]],
#                                         "st_err" = res[["Std. Error"]],
#                                         "pval" = res[["Pr(>|t|)"]]))
#     }
# }
# lnres_cont_ilbd$padj <- p.adjust(lnres_cont_ilbd$pval, method="BH")


# control ilbd case: 0 1 2
lnres_diagnosis <- c()
for(ct in all_cell_type_cols) {
    for(layer in layers) {
        cell_props_sub <- cell_props_spotnorm[cell_props_spotnorm$smoothed_label_s5 == layer, ]
        cell_props_sub$diagnosis_int <- NA
        cell_props_sub$diagnosis_int[cell_props_sub$diagnosis == "Control"] <- 0
        cell_props_sub$diagnosis_int[cell_props_sub$diagnosis == "ILBD"] <- 1
        cell_props_sub$diagnosis_int[cell_props_sub$diagnosis == "Case"] <- 2
        assert(is.numeric(cell_props_sub$diagnosis_int))
        cell_props_sub$batch <- factor(cell_props_sub$batch)
        cell_props_sub$sex <- factor(cell_props_sub$sex)
        cell_props_sub$ct_to_test <- cell_props_sub[, ct]

        mdl <- lm(data = cell_props_sub,
                ct_to_test ~ batch + sex + expired_age_num + PMI_num + RIN_num + diagnosis_int)
        res <- summary(mdl)$coefficients["diagnosis_int", ]
        lnres_diagnosis <- rbind(lnres_diagnosis, data.frame(
                                        "clinical_feature" = "diagnosis_int",
                                        "cell_type" = ct,
                                        "layer" = layer,
                                        "estimate" = res[["Estimate"]],
                                        "st_err" = res[["Std. Error"]],
                                        "pval" = res[["Pr(>|t|)"]]))
    }
}
lnres_diagnosis$padj <- p.adjust(lnres_diagnosis$pval, method="BH")


# linear nctx_temporal_fix
lnres_nctx <- c()
for(ct in all_cell_type_cols) {
    for(layer in layers) {
        cell_props_sub <- cell_props_spotnorm[cell_props_spotnorm$smoothed_label_s5 == layer, ]
        assert(is.numeric(cell_props_sub$nctx_temporal_fix))
        cell_props_sub$batch <- factor(cell_props_sub$batch)
        cell_props_sub$sex <- factor(cell_props_sub$sex)
        cell_props_sub$ct_to_test <- cell_props_sub[, ct]

        mdl <- lm(data = cell_props_sub,
                ct_to_test ~ batch + sex + expired_age_num + PMI_num + RIN_num + nctx_temporal_fix)
        res <- summary(mdl)$coefficients["nctx_temporal_fix", ]
        lnres_nctx <- rbind(lnres_nctx, data.frame(
                                        "clinical_feature" = "nctx_temporal_fix",
                                        "cell_type" = ct,
                                        "layer" = layer,
                                        "estimate" = res[["Estimate"]],
                                        "st_err" = res[["Std. Error"]],
                                        "pval" = res[["Pr(>|t|)"]]))
    }
}
lnres_nctx$padj <- p.adjust(lnres_nctx$pval, method="BH")

# linear brain_stem_sn


pdf("st_sc_ct_pca_04112024.pdf")
for(layer in layers) {
    cell_props_sub <- cell_props_spotnorm[cell_props_spotnorm$smoothed_label_s5 == layer, ]
    num_pc_components <- 10
    ct_pca <- prcomp(cell_props_sub[, all_cell_type_cols], scale=TRUE)
    ct_pca <- ct_pca$x
    cell_props_sub <- cbind(cell_props_sub, ct_pca[, 1:10])
    
    ct_umap <- umap(ct_pca[, 1:10])
    ct_umap <- ct_umap$layout
    colnames(ct_umap) <- c("UMAP1", "UMAP2")
    cell_props_sub <- cbind(cell_props_sub, ct_umap)
    
    pt <- ggplot(cell_props_sub, aes(PC1, PC2)) + geom_point(aes(colour = factor(diagnosis))) + ggtitle(layer)
    print(pt)
    pt <- ggplot(cell_props_sub, aes(UMAP1, UMAP2)) + geom_point(aes(colour = factor(diagnosis))) + ggtitle(layer)
    print(pt)
}
dev.off()

gaba_cell_types <- all_cell_type_cols[grepl("^ADARB2|^GAD1|LAMP5\\.CA1|^LAMP5\\.LCP2|LAMP5\\.NMBR|^PAX6|^PVALB|^SST|^VIP", all_cell_type_cols)]
glu_cell_types <- all_cell_type_cols[grepl("^FEZF2|LAMP5\\.GREB1L|LAMP5\\.PARD3|^LINC00507|^RORB|^SLC17A7|^THEMIS", all_cell_type_cols)]
assert(length(intersect(glu_cell_types, gaba_cell_types)) == 0)
# ct_pca <- prcomp(cell_props[, gaba_cell_types], scale=TRUE)
# ct_pca <- prcomp(cell_props[, gaba_cell_types], scale=FALSE)
ct_pca <- prcomp(cell_props_spotnorm[, gaba_cell_types], scale=FALSE)
ct_pca <- ct_pca$x
# ct_umap <- umap(ct_pca[, 1:30])
# ct_umap <- ct_umap$layout
# colnames(ct_umap) <- c("UMAP1", "UMAP2")
# cell_props_temp <- cell_props
cell_props_temp <- cell_props_spotnorm
cell_props_temp <- cbind(cell_props_temp, ct_pca[, 1:10])
pdf("st_sc_ct_pca_gaba_only_grid_04162024.pdf")
pairs.panels(ct_pca[, 1:5],
             gap=0,
             bg = cell_props_temp$diagnosis,
             pch=21)
dev.off()

pdf("st_sc_ct_pca_gaba_only_04162024.pdf")
pt <- ggplot(cell_props_temp, aes(PC1, PC2)) + geom_point(aes(colour = factor(diagnosis)), size=1) + ggtitle("GABA only") + theme_bw()
print(pt)
# for(sid in unique(cell_props_temp$sample_name)) {
#     print(sid)
#     cell_props_temp2 <- cell_props_temp[cell_props_temp$sample_name == sid, ]
#     pt <- ggplot(cell_props_temp2, aes(PC1, PC2)) + geom_point(aes(colour = factor(diagnosis)), size=1) + ggtitle(sid) + theme_bw()
#     print(pt)
# }
dev.off()

# -----
# Cell type composition vs PD progression (using raw pseudoabundance from cell2location)
# -----
# in nctx_temporal, there is only 1 sample with a score of 4: set this to 3 so as not to overly bias the regression
comb_seurat@meta.data$nctx_temporal_fix <- comb_seurat@meta.data$nctx_temporal
comb_seurat@meta.data$nctx_temporal_fix[comb_seurat@meta.data$nctx_temporal_fix > 3] <- 3

cell_props <- comb_seurat@meta.data[, c("sample_name", "diagnosis", "brain_stem_sn", "nctx_temporal_fix", "smoothed_label_s5", "cell2loc_sum", "batch", "sex", "expired_age_num", "PMI_num", "RIN_num", all_cell_type_cols)]
# remove samples which have 0 cell2loc info
temp <- cell_props %>%
    group_by(sample_name) %>%
    summarise(cell2loc_sample_sum = sum(cell2loc_sum)) %>%
    as.data.frame()
include_samples <- as.character(temp[temp$cell2loc_sample_sum > 0, "sample_name"])
print(paste("Num samples", length(unique(cell_props$sample_name))))
cell_props <- cell_props[cell_props$sample_name %in% include_samples, ]
print(paste("Num samples after filtering", length(unique(cell_props$sample_name))))

cell_props_psabd <- cell_props %>%
        group_by(sample_name, smoothed_label_s5, diagnosis, brain_stem_sn, nctx_temporal_fix, batch, sex, expired_age_num, PMI_num, RIN_num) %>%
        summarise_at(all_cell_type_cols, mean, na.rm = FALSE) %>%
        as.data.frame()

# outname_base <- "03142024"
# outname_base <- "04102024"
outname_base <- "04152024"
layers <- sort(unique(cell_props_psabd$smoothed_label_s5))
for(ct in all_cell_type_cols) {
    
    for(clinical_feat_discrt in c("diagnosis", "brain_stem_sn", "nctx_temporal_fix")) {
        outname <- paste(ct, clinical_feat_discrt, outname_base, sep = "_")
        pdf(sprintf("st_sc_scatterplot_psabd_V4_%s.pdf", outname))
        for(layer in layers) {
            cell_props_sub <- cell_props_psabd[cell_props_psabd$smoothed_label_s5 == layer, ]
            cell_props_sub[, clinical_feat_discrt] <- factor(cell_props_sub[, clinical_feat_discrt])
            cell_props_sub <- cell_props_sub[!is.na(cell_props_sub[, clinical_feat_discrt]), ]
            cell_props_sub$diagnosis <- factor(cell_props_sub$diagnosis, levels = c("Control", "ILBD", "Case"))
            pt <- ggplot(cell_props_sub, aes_string(x = clinical_feat_discrt, y = ct, fill = clinical_feat_discrt)) + geom_boxplot() + geom_jitter() + ggtitle(paste(ct, layer, clinical_feat_discrt)) + theme_bw()
            print(pt)
        }
        dev.off()
    }
}
# linear nctx_temporal_fix
lnres_nctx <- c()
for(ct in all_cell_type_cols) {
    for(layer in layers) {
        cell_props_sub <- cell_props_psabd[cell_props_psabd$smoothed_label_s5 == layer, ]
        assert(is.numeric(cell_props_sub$nctx_temporal_fix))
        cell_props_sub$batch <- factor(cell_props_sub$batch)
        cell_props_sub$sex <- factor(cell_props_sub$sex)
        cell_props_sub$ct_to_test <- cell_props_sub[, ct]

        mdl <- lm(data = cell_props_sub,
                ct_to_test ~ batch + sex + expired_age_num + PMI_num + RIN_num + nctx_temporal_fix)
        res <- summary(mdl)$coefficients["nctx_temporal_fix", ]
        lnres_nctx <- rbind(lnres_nctx, data.frame(
                                        "clinical_feature" = "nctx_temporal_fix",
                                        "cell_type" = ct,
                                        "layer" = layer,
                                        "estimate" = res[["Estimate"]],
                                        "st_err" = res[["Std. Error"]],
                                        "pval" = res[["Pr(>|t|)"]]))
    }
}
lnres_nctx$padj <- p.adjust(lnres_nctx$pval, method="BH")























# for each cell type, in each layer, plot the mean spot-normalized cell count for all individuals according to clinical feature on x-axis

pt <- ggplot(all_pathways, aes(x=nctx_temporal, y=leading_edge_ave_mean, group=interaction(layer, pathway),  color=layer, linetype=pathway)) + geom_line(position=position_dodge(0.1), size=1) + geom_point(position=position_dodge(0.1)) +
    geom_errorbar(aes(ymin=leading_edge_ave_mean-leading_edge_ave_sd, ymax=leading_edge_ave_mean+leading_edge_ave_sd), width=.4, size=1,
        position=position_dodge(0.1)) + xlab("MTG Lewy body load") + ylab("GeomMean log2norm expression of leading edges") + theme_bw()




# out_dir <- "/data/neurogen/ASAP/ST/data/sc_for_deconvolution"

# sample_list <- unique(comb_sc_all@meta.data$sample_id)
# for(sample in sample_list) {
#     comb_sub <- subset(comb_sc_all, subset = sample_id == sample)
#
#     write.table(comb_sub@meta.data,
#                 file.path(out_dir, sprintf("sc_count_%s_meta_02062023.csv", sample)),
#                 sep = ',', row.names = T, col.names = NA, quote = F)
#     write.table(as.matrix(GetAssayData(object = comb_sub, slot = "counts")),
#                 file.path(out_dir, sprintf("sc_count_%s_02062023.csv", sample)),
#                 sep = ',', row.names = T, col.names = NA, quote = F)
# }


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
