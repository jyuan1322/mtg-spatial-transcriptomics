library(Seurat)
library(data.table)
library(testit)
library(ggplot2)
library(gridExtra)

library(reshape2)
library(tidyr)
library(e1071)

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

comb_seurat <- readRDS("/mnt/data0/projects/ASAP/ST/clustering_pipeline/post_clustering_validation_V4_seed101_10172023/comb_seurat_diet_out_plusL4_V4_seed101_11092023.rds")


pdf("layer_smoothing_s0_all_layer_label_v2_11302023.pdf")
for(img in names(comb_seurat@images)) {
    print(img)
    meta_row <- rownames(comb_seurat@images[[img]]@coordinates)[1]
    sample_name_filt <- comb_seurat@meta.data[meta_row, "orig.ident"]
    diagnosis_filt <- comb_seurat@meta.data[meta_row, "diagnosis"]
    
    sp_pt <- SpatialPlot(comb_seurat, group.by="layer_label_v2", image.alpha=0.3, images=img, stroke=NA) + scale_fill_manual(values=group.colors.unsmooth) + ggtitle(paste(sample_name_filt, diagnosis_filt))
    print(sp_pt)
}
dev.off()


# Stage 1: Smoothing L2-6
# default gamma is 1/(data dim); here dim is 2 (x and y coordinate)
# large gamma -> low bias, high variance
# large cost -> low bias, high variance
# order matters here, since if the best score is achieved, the loop breaks
# we prefer lower L4 weights, low bias, high variance models, which are closer to the original data
# svm_gamma_range <- c(5, 4, 3, 2)
svm_gamma_range <- c(3)
# svm_cost_range <- c(100, 75, 50, 25)
svm_cost_range <- c(50)
L4_weight_range <- c(1.0, 1.25, 1.5, 1.75, 2.0)
cluster_label_name <- "layer_label_v2"

comb_seurat@meta.data[, "smoothed_label_s1"] <- NULL
comb_seurat_L4 <- subset(comb_seurat, subset = layer_label_v2 %in% c("Layer 2", "Layer 3", "Layer 4", "Layer 5", "Layer 6"))
comb_seurat@meta.data[, "smoothed_label_s1"] <- NA
scan_params_all <- c()
pdf("layer_smoothing_s1_L2-6_all_11302023.pdf")
for( img in names(comb_seurat_L4@images) ) {
    meta_row <- rownames(comb_seurat_L4@images[[img]]@coordinates)[1]
    sample_name_filt <- comb_seurat_L4@meta.data[meta_row, "sample_name"]
    diagnosis_filt <- comb_seurat_L4@meta.data[meta_row, "diagnosis"]
    
    scan_params <- data.frame()
    best_score_achieved <- FALSE
    for(svm_gamma in svm_gamma_range) {
        for(svm_cost in svm_cost_range) {
            for(L4_weight in L4_weight_range) {
                class_weights <- c("Layer 2" = 1, "Layer 3" = 1, "Layer 4" = L4_weight,
                                   "Layer 5" = 1, "Layer 6" = 1)
                comb_seurat_L4_temp <- svm_smooth_labels(comb_seurat_L4,
                                            cluster_label_name,
                                            svm_gamma = svm_gamma,
                                            svm_cost = svm_cost,
                                            class_weights = class_weights,
                                            image_name = img,
                                            output_label_name = "smoothed_label_s1",
                                            random_seed = random_seed)
                layer_scores <- score_layer_fit(comb_seurat_L4_temp,
                                               cluster_label_name = "smoothed_label_s1",
                                               image_name=img, verbose=FALSE)
                layer_score <- layer_scores[[img]]
                # l4_penalty <- score_Layer4_fraction(comb_seurat_L4_temp,
                #                                orig_col = "layer_label_v2",
                #                                new_col = "smoothed_label_s1",
                #                                image_name = img,
                #                                penalty_thresh = 1.4)
                l4_penalty <- score_layer_fraction(comb_seurat_L4_temp,
                                                 orig_col = "layer_label_v2",
                                                 new_col = "smoothed_label_s1",
                                                 image_name = img,
                                                 layer_set = c("Layer 4"),
                                                 penalty_thresh = 1.4)
                
                scan_params <- rbind(scan_params, data.frame(sample_name = sample_name_filt,
                                                             svm_gamma = svm_gamma,
                                                             svm_cost = svm_cost,
                                                             L4_weight = L4_weight,
                                                             layer_score = layer_score,
                                                             l4_penalty = l4_penalty))
                print(paste(sample_name_filt, svm_gamma, svm_cost, L4_weight, layer_score, l4_penalty))
                # if optimal score is reached, break the loop
                if(layer_score == 1 && l4_penalty == 0) {
                    best_score_achieved <- TRUE
                    break
                }
            }
            if(best_score_achieved) {break}
        }
        if(best_score_achieved) {break}
    }
    # get best-performing params
    scan_params <- scan_params[order(scan_params[, "layer_score"] - 10*scan_params[, "l4_penalty"], decreasing=TRUE), ]
    scan_params_all <- rbind(scan_params_all, scan_params)
    best_svm_gamma <- scan_params[1, "svm_gamma"]
    best_svm_cost <- scan_params[1, "svm_cost"]
    best_L4_weight <- scan_params[1, "L4_weight"]
    class_weights <- c("Layer 2" = 1, "Layer 3" = 1, "Layer 4" = best_L4_weight,
                       "Layer 5" = 1, "Layer 6" = 1)
    comb_seurat_L4_temp <- svm_smooth_labels(comb_seurat_L4,
                                cluster_label_name,
                                svm_gamma = best_svm_gamma,
                                svm_cost = best_svm_cost,
                                class_weights = class_weights,
                                image_name = img,
                                output_label_name = "smoothed_label_s1",
                                random_seed = random_seed)

    
    sp_pt <- SpatialPlot(comb_seurat_L4_temp, group.by="smoothed_label_s1", image.alpha=0.3, images=img, stroke=NA) + scale_fill_manual(values=group.colors) + ggtitle(paste(sample_name_filt, diagnosis_filt))
    print(sp_pt)
    

    # # add best smoothed annotation to combined seurat object
    # smoothed_label_name <- "smoothed_label_v1"
    # if( !(smoothed_label_name %in% colnames(comb_seurat@meta.data)) ) {
    #     comb_seurat@meta.data[, smoothed_label_name] <- NA
    # }
    # csL4_rownames <- rownames(comb_sub_L4@meta.data)
    # comb_seurat@meta.data[csL4_rownames, smoothed_label_name] <- comb_sub_L4@meta.data[csL4_rownames, "smoothed_label"]
    img_rows <- rownames(comb_seurat_L4_temp@images[[img]]@coordinates)
    comb_seurat@meta.data[img_rows, "smoothed_label_s1"] <- comb_seurat_L4_temp@meta.data[img_rows, "smoothed_label_s1"]
}
dev.off()
write.csv(scan_params_all, "layer_smoothing_s1_L2-6_all_scan_params_11302023.csv")



comb_seurat@meta.data[, "s2_subclust"] <- comb_seurat@meta.data[, "smoothed_label_s1"]

# 11/30/2023 Don't sub-cluster WM and L1
comb_seurat@meta.data[comb_seurat@meta.data$layer_label_v2 == "Layer 1", "s2_subclust"] <- "Layer 1"
comb_seurat@meta.data[comb_seurat@meta.data$layer_label_v2 == "WM_2", "s2_subclust"] <- "WM_2"
comb_seurat@meta.data[comb_seurat@meta.data$layer_label_v2 == "WM_3", "s2_subclust"] <- "WM_3"

pdf("layer_smoothing_s2_prefilter_11302023.pdf")
for(img in names(comb_seurat@images)) {
    meta_row <- rownames(comb_seurat@images[[img]]@coordinates)[1]
    sample_name_filt <- comb_seurat@meta.data[meta_row, "orig.ident"]
    diagnosis_filt <- comb_seurat@meta.data[meta_row, "diagnosis"]
    sp_pt <- SpatialPlot(comb_seurat, group.by="s2_subclust", image.alpha=0.3, images=img, stroke=NA) + scale_fill_manual(values=group.colors.unsmooth) + ggtitle(paste(sample_name_filt, diagnosis_filt))
    print(sp_pt)
}
dev.off()

# WM_2 and L1 subclusters (res=0.1) beside cluster 1 are all subject to re-labeling

# for each spot, calculate score based on gaussian kernel
# L1_1: 1, 2, 3, 4, 5, 6, WM_3: 7

s2_score_table <- list("Layer 2" = 2, "Layer 3" = 3, "Layer 4" = 4, "Layer 5" = 5, "Layer 6" = 6)
# scores of 0 are omitted from the kernel weighted average
comb_seurat@meta.data[, "s2_subclust_score"] <- 0
for(s2_layer in names(s2_score_table)) {
    comb_seurat@meta.data[ (!is.na(comb_seurat@meta.data$s2_subclust) & comb_seurat@meta.data$s2_subclust == s2_layer), "s2_subclust_score" ] <- s2_score_table[[s2_layer]]
}

cluster_labels <- comb_seurat@meta.data[, c("orig.ident", "s2_subclust_score")]
rbf_gamma <- 0.01 # gamma: low value means more long-range influence
for (img in names(comb_seurat@images)) {
    print(paste("Processing stage 2 subcluster score: ", img))
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
    zero_s2_score_rows <- which(coord_table$s2_subclust_score == 0)
    adjA[, zero_s2_score_rows] <- 0
    adjA <- sweep(adjA, 1, rowSums(adjA), FUN="/")

    coord_table$s2_nbr_score <- adjA %*% coord_table$s2_subclust_score
    comb_seurat@meta.data[rownames(coord_table), "s2_nbr_score"] <- coord_table$s2_nbr_score
}
pdf(sprintf("layer_smoothing_s2_nbr_score_all_gamma%s_11302023.pdf", rbf_gamma))
for(img in names(comb_seurat@images)) {
    print(img)
    meta_row <- rownames(comb_seurat@images[[img]]@coordinates)[1]
    sample_name_filt <- comb_seurat@meta.data[meta_row, "orig.ident"]
    diagnosis_filt <- comb_seurat@meta.data[meta_row, "diagnosis"]
    sp_pt <- SpatialFeaturePlot(comb_seurat, features="s2_nbr_score", image.alpha=0.3, images=img, stroke=NA) + ggtitle(paste(sample_name_filt, diagnosis_filt))
    print(sp_pt)
}
dev.off()

# barplots of smoothed score split by layer
pdf(sprintf("layer_smoothing_s2_nbr_score_perlayer_hists_gamma%s_11302023.pdf", rbf_gamma))
for(img in names(comb_seurat@images)) {
    img_rows <- rownames(comb_seurat@images[[img]]@coordinates)
    sample_name_filt <- comb_seurat@meta.data[img_rows[1], "orig.ident"]
    diagnosis_filt <- comb_seurat@meta.data[img_rows[1], "diagnosis"]
    
    meta_sub <- comb_seurat@meta.data[img_rows, ]
    sp_pt <- ggplot(meta_sub, aes(x=s2_nbr_score, fill=layer_label_v2)) + geom_histogram(alpha=0.2, position="identity", binwidth=0.1) + ggtitle(paste(sample_name_filt, diagnosis_filt))
    print(sp_pt)
}
dev.off()

# filter WM_2 and WM_3
pdf(sprintf("layer_smoothing_s2_nbr_score_WM_hists_gamma%s_11302023.pdf", rbf_gamma))
for(img in names(comb_seurat@images)) {
    img_rows <- rownames(comb_seurat@images[[img]]@coordinates)
    sample_name_filt <- comb_seurat@meta.data[img_rows[1], "orig.ident"]
    diagnosis_filt <- comb_seurat@meta.data[img_rows[1], "diagnosis"]
    
    meta_sub <- comb_seurat@meta.data[img_rows, ]
    # meta_sub_wm2 <- meta_sub[meta_sub$s2_subclust == "WM_2", ]
    # wm2_pt <- ggplot(meta_sub_wm2, aes(x=s2_nbr_score)) + geom_histogram() + ggtitle(paste(sample_name_filt, diagnosis_filt, "WM_2"))
    # meta_sub_wm3 <- meta_sub[meta_sub$s2_subclust == "WM_3", ]
    # wm3_pt <- ggplot(meta_sub_wm3, aes(x=s2_nbr_score)) + geom_histogram() + ggtitle(paste(sample_name_filt, diagnosis_filt, "WM_3"))
    
    # sp_pt <- grid.arrange(wm2_pt, wm3_pt, ncol=2)
    
    meta_sub_wm <- meta_sub[meta_sub$s2_subclust %in% c("WM_2", "WM_3"), ]
    sp_pt <- ggplot(meta_sub_wm, aes(x=s2_nbr_score)) + geom_histogram() + ggtitle(paste(sample_name_filt, diagnosis_filt)) + theme_bw()
    print(sp_pt)
}
dev.off()


# filter L1_1, L1_2, L1_3
pdf(sprintf("layer_smoothing_s2_nbr_score_L1_hists_gamma%s_11302023.pdf", rbf_gamma))
for(img in names(comb_seurat@images)) {
    img_rows <- rownames(comb_seurat@images[[img]]@coordinates)
    sample_name_filt <- comb_seurat@meta.data[img_rows[1], "orig.ident"]
    diagnosis_filt <- comb_seurat@meta.data[img_rows[1], "diagnosis"]
    
    meta_sub <- comb_seurat@meta.data[img_rows, ]
    # meta_sub_L1_1 <- meta_sub[meta_sub$s2_subclust == "L1_1", ]
    # L1_1_pt <- ggplot(meta_sub_L1_1, aes(x=s2_nbr_score)) + geom_histogram() + ggtitle(paste(sample_name_filt, diagnosis_filt, "L1_1"))
    # meta_sub_L1_2 <- meta_sub[meta_sub$s2_subclust == "L1_2", ]
    # L1_2_pt <- ggplot(meta_sub_L1_2, aes(x=s2_nbr_score)) + geom_histogram() + ggtitle(paste(sample_name_filt, diagnosis_filt, "L1_2"))
    # meta_sub_L1_3 <- meta_sub[meta_sub$s2_subclust == "L1_3", ]
    # L1_3_pt <- ggplot(meta_sub_L1_3, aes(x=s2_nbr_score)) + geom_histogram() + ggtitle(paste(sample_name_filt, diagnosis_filt, "L1_3"))
    
    # sp_pt <- grid.arrange(L1_1_pt, L1_2_pt, L1_3_pt, ncol=3)
    
    meta_sub_L1 <- meta_sub[meta_sub$s2_subclust %in% c("Layer 1"), ]
    sp_pt <- ggplot(meta_sub_L1, aes(x=s2_nbr_score)) + geom_histogram() + ggtitle(paste(sample_name_filt, diagnosis_filt)) + theme_bw()
    print(sp_pt)
}
dev.off()

comb_seurat@meta.data$smoothed_label_s2 <- comb_seurat@meta.data$s2_subclust
# Set WM to NA if s2_nbr_score less than 3
WM_lt_thresh <- 3
# WM_3 tends to confound with L1, whereas WM_2 tends to confound with GM
# set WM_3 spots failing threshold to L1, but set WM_2 spots failing threshold to NA (to be imputed later)
comb_seurat@meta.data[comb_seurat@meta.data$s2_subclust %in% c("WM_2", "WM_3"), "smoothed_label_s2"] <- "WM"

comb_seurat@meta.data[comb_seurat@meta.data$s2_subclust %in% c("WM_2", "WM_3") & comb_seurat@meta.data$s2_nbr_score < WM_lt_thresh, "smoothed_label_s2"] <- NA

# Set L1 to NA if s2_nbr_score greater than 4
L1_gt_thresh <- 4
comb_seurat@meta.data[comb_seurat@meta.data$s2_subclust %in% c("Layer 1"), "smoothed_label_s2"] <- "Layer 1"
comb_seurat@meta.data[comb_seurat@meta.data$s2_subclust %in% c("Layer 1") & comb_seurat@meta.data$s2_nbr_score > L1_gt_thresh, "smoothed_label_s2"] <- NA

pdf(sprintf("layer_smoothing_s2_postfilter_label_all_gamma%s_L1%s_WM%s_11302023.pdf", rbf_gamma, L1_gt_thresh, WM_lt_thresh))
for(img in names(comb_seurat@images)) {
    print(img)
    meta_row <- rownames(comb_seurat@images[[img]]@coordinates)[1]
    sample_name_filt <- comb_seurat@meta.data[meta_row, "orig.ident"]
    diagnosis_filt <- comb_seurat@meta.data[meta_row, "diagnosis"]
    sp_pt <- SpatialPlot(comb_seurat, group.by="smoothed_label_s2", image.alpha=0.3, images=img, stroke=NA) + scale_fill_manual(values=group.colors) + ggtitle(paste(sample_name_filt, diagnosis_filt))
    print(sp_pt)
}
dev.off()



# stage 3 smoothing: L2-6 and WM
svm_gamma_range <- c(1)
svm_cost_range <- c(50)
WM_weight_range <- c(1.0, 0.75, 0.5, 0.25)
cluster_label_name <- "smoothed_label_s2"

comb_seurat@meta.data[, "smoothed_label_s3"] <- NULL
comb_seurat_WM <- subset(comb_seurat, subset = smoothed_label_s2 %in% c("Layer 2", "Layer 3", "Layer 4", "Layer 5", "Layer 6", "WM"))
comb_seurat@meta.data[, "smoothed_label_s3"] <- NA

comb_seurat_WM@meta.data$smoothed_label_s2 <- as.factor(comb_seurat_WM@meta.data$smoothed_label_s2)

scan_params_all <- c()
pdf("layer_smoothing_s3_L2-6WM_all_11302023.pdf")
for( img in names(comb_seurat_WM@images) ) {
    meta_row <- rownames(comb_seurat_WM@images[[img]]@coordinates)[1]
    sample_name_filt <- comb_seurat_WM@meta.data[meta_row, "sample_name"]
    diagnosis_filt <- comb_seurat_WM@meta.data[meta_row, "diagnosis"]
    
    print(paste(sample_name_filt, diagnosis_filt))
    
    scan_params <- data.frame()
    best_score_achieved <- FALSE
    for(svm_gamma in svm_gamma_range) {
        for(svm_cost in svm_cost_range) {
            for(WM_weight in WM_weight_range) {
                class_weights <- c("Layer 2" = 1, "Layer 3" = 1, "Layer 4" = 1,
                                   "Layer 5" = 1, "Layer 6" = 1, "WM" = WM_weight)
                comb_seurat_WM_temp <- svm_smooth_labels(comb_seurat_WM,
                                            cluster_label_name,
                                            svm_gamma = svm_gamma,
                                            svm_cost = svm_cost,
                                            class_weights = class_weights,
                                            image_name = img,
                                            output_label_name = "smoothed_label_s3",
                                            random_seed = random_seed)
                layer_scores <- score_layer_fit(comb_seurat_WM_temp,
                                               cluster_label_name = "smoothed_label_s3",
                                               image_name=img, verbose=FALSE)
                layer_score <- layer_scores[[img]]
                # WM_penalty <- score_WM_fraction(comb_seurat_WM_temp,
                #                                 orig_col = "smoothed_label_s2",
                #                                 new_col = "smoothed_label_s3",
                #                                 image_name = img,
                #                                 penalty_thresh = 1.4)
                WM_penalty <- score_layer_fraction(comb_seurat_WM_temp,
                                                orig_col = "smoothed_label_s2",
                                                new_col = "smoothed_label_s3",
                                                image_name = img,
                                                layer_set = c("Layer 6", "WM"),
                                                penalty_thresh = 1.4)

                scan_params <- rbind(scan_params, data.frame(sample_name = sample_name_filt,
                                                             svm_gamma = svm_gamma,
                                                             svm_cost = svm_cost,
                                                             WM_weight = WM_weight,
                                                             layer_score = layer_score,
                                                             WM_penalty = WM_penalty))
                print(paste(sample_name_filt, svm_gamma, svm_cost, WM_weight, layer_score, WM_penalty))
                # if optimal score is reached, break the loop
                if(layer_score == 1 && WM_penalty == 0) {
                    best_score_achieved <- TRUE
                    break
                }
            }
            if(best_score_achieved) {break}
        }
        if(best_score_achieved) {break}
    }
    # get best-performing params
    scan_params <- scan_params[order(scan_params[, "layer_score"] - 10*scan_params[, "WM_penalty"], decreasing=TRUE), ]
    scan_params_all <- rbind(scan_params_all, scan_params)
    best_svm_gamma <- scan_params[1, "svm_gamma"]
    best_svm_cost <- scan_params[1, "svm_cost"]
    best_WM_weight <- scan_params[1, "WM_weight"]
    class_weights <- c("Layer 2" = 1, "Layer 3" = 1, "Layer 4" = 1,
                       "Layer 5" = 1, "Layer 6" = 1, "WM" = best_WM_weight)
    comb_seurat_WM_temp <- svm_smooth_labels(comb_seurat_WM,
                                cluster_label_name,
                                svm_gamma = best_svm_gamma,
                                svm_cost = best_svm_cost,
                                class_weights = class_weights,
                                image_name = img,
                                output_label_name = "smoothed_label_s3",
                                random_seed = random_seed)

    
    sp_pt <- SpatialPlot(comb_seurat_WM_temp, group.by="smoothed_label_s3", image.alpha=0.3, images=img, stroke=NA) + scale_fill_manual(values=group.colors) + ggtitle(paste(sample_name_filt, diagnosis_filt))
    print(sp_pt)
    

    # # add best smoothed annotation to combined seurat object
    # smoothed_label_name <- "smoothed_label_v1"
    # if( !(smoothed_label_name %in% colnames(comb_seurat@meta.data)) ) {
    #     comb_seurat@meta.data[, smoothed_label_name] <- NA
    # }
    # csL4_rownames <- rownames(comb_sub_L4@meta.data)
    # comb_seurat@meta.data[csL4_rownames, smoothed_label_name] <- comb_sub_L4@meta.data[csL4_rownames, "smoothed_label"]
    img_rows <- rownames(comb_seurat_WM_temp@images[[img]]@coordinates)
    comb_seurat@meta.data[img_rows, "smoothed_label_s3"] <- comb_seurat_WM_temp@meta.data[img_rows, "smoothed_label_s3"]
}
dev.off()
write.csv(scan_params_all, "layer_smoothing_s3_L2-6WM_all_scan_params_11302023.csv")








# stage 4 smoothing: L1 and all others
svm_gamma_range <- c(10)
svm_cost_range <- c(50)
L1_weight_range <- c(2.0, 1.75, 1.5, 1.25, 1.0)
cluster_label_name <- "smoothed_label_s3"

comb_seurat@meta.data[, "smoothed_label_s4"] <- NA
# transfer Layer 1 label to latest smoothed label set
comb_seurat@meta.data$smoothed_label_s3[comb_seurat@meta.data$smoothed_label_s2 == "Layer 1"] <- "Layer 1"
comb_seurat@meta.data$smoothed_label_s3 <- as.factor(comb_seurat@meta.data$smoothed_label_s3)

pdf("layer_smoothing_s3_L2-6WM_all_plusL1temp_11302023.pdf")
for( img in names(comb_seurat@images) ) {
    meta_row <- rownames(comb_seurat@images[[img]]@coordinates)[1]
    sample_name_filt <- comb_seurat@meta.data[meta_row, "sample_name"]
    diagnosis_filt <- comb_seurat@meta.data[meta_row, "diagnosis"]

    sp_pt <- SpatialPlot(comb_seurat, group.by="smoothed_label_s3", image.alpha=0.3, images=img, stroke=NA) + scale_fill_manual(values=group.colors) + ggtitle(paste(sample_name_filt, diagnosis_filt))
    print(sp_pt)
}
dev.off()

# ignore NA for now
comb_seurat_L1 <- subset(comb_seurat, subset = smoothed_label_s3 %in% c("Layer 1", "Layer 2", "Layer 3", "Layer 4", "Layer 5", "Layer 6", "WM"))

scan_params_all <- c()
pdf("layer_smoothing_s4_L1-6WM_all_11302023.pdf")
for( img in names(comb_seurat_L1@images) ) {
    meta_row <- rownames(comb_seurat_L1@images[[img]]@coordinates)[1]
    sample_name_filt <- comb_seurat_L1@meta.data[meta_row, "sample_name"]
    diagnosis_filt <- comb_seurat_L1@meta.data[meta_row, "diagnosis"]
    
    print(paste(sample_name_filt, diagnosis_filt))
    
    scan_params <- data.frame()
    best_score_achieved <- FALSE
    for(svm_gamma in svm_gamma_range) {
        for(svm_cost in svm_cost_range) {
            for(L1_weight in L1_weight_range) {
                class_weights <- c("Layer 1" = L1_weight, "Layer 2" = 1, "Layer 3" = 1,
                                   "Layer 4" = 1, "Layer 5" = 1, "Layer 6" = 1, "WM" = 1)
                comb_seurat_L1_temp <- svm_smooth_labels(comb_seurat_L1,
                                            cluster_label_name,
                                            svm_gamma = svm_gamma,
                                            svm_cost = svm_cost,
                                            class_weights = class_weights,
                                            image_name = img,
                                            output_label_name = "smoothed_label_s4",
                                            random_seed = random_seed)
                layer_scores <- score_layer_fit(comb_seurat_L1_temp,
                                               cluster_label_name = "smoothed_label_s4",
                                               image_name=img, verbose=FALSE)
                layer_score <- layer_scores[[img]]
                L1_penalty <- score_layer_fraction(comb_seurat_L1_temp,
                                                orig_col = "smoothed_label_s3",
                                                new_col = "smoothed_label_s4",
                                                image_name = img,
                                                layer_set = c("Layer 1", "Layer 2"),
                                                penalty_thresh = 1.4)

                scan_params <- rbind(scan_params, data.frame(sample_name = sample_name_filt,
                                                             svm_gamma = svm_gamma,
                                                             svm_cost = svm_cost,
                                                             L1_weight = L1_weight,
                                                             layer_score = layer_score,
                                                             L1_penalty = L1_penalty))
                print(paste(sample_name_filt, svm_gamma, svm_cost, L1_weight, layer_score, L1_penalty))
                # if optimal score is reached, break the loop
                if(layer_score == 1 && L1_penalty == 0) {
                    best_score_achieved <- TRUE
                    break
                }
            }
            if(best_score_achieved) {break}
        }
        if(best_score_achieved) {break}
    }
    # get best-performing params
    # we reduced the layer proportion penalty here, becaues in many tissues there is no true Layer 1 but a lot of single mislabeled spots - now it only serves as a tie breaker for multiple runs with the same layer score
    scan_params <- scan_params[order(scan_params[, "layer_score"] - 0.001*scan_params[, "L1_penalty"], decreasing=TRUE), ]
    scan_params_all <- rbind(scan_params_all, scan_params)
    best_svm_gamma <- scan_params[1, "svm_gamma"]
    best_svm_cost <- scan_params[1, "svm_cost"]
    best_L1_weight <- scan_params[1, "L1_weight"]
    class_weights <- c("Layer 1" = best_L1_weight, "Layer 2" = 1, "Layer 3" = 1,
                       "Layer 4" = 1, "Layer 5" = 1, "Layer 6" = 1, "WM" = 1)
    comb_seurat_L1_temp <- svm_smooth_labels(comb_seurat_L1,
                                cluster_label_name,
                                svm_gamma = best_svm_gamma,
                                svm_cost = best_svm_cost,
                                class_weights = class_weights,
                                image_name = img,
                                output_label_name = "smoothed_label_s4",
                                random_seed = random_seed)

    
    sp_pt <- SpatialPlot(comb_seurat_L1_temp, group.by="smoothed_label_s4", image.alpha=0.3, images=img, stroke=NA) + scale_fill_manual(values=group.colors) + ggtitle(paste(sample_name_filt, diagnosis_filt))
    print(sp_pt)

    img_rows <- rownames(comb_seurat_L1_temp@images[[img]]@coordinates)
    comb_seurat@meta.data[img_rows, "smoothed_label_s4"] <- comb_seurat_L1_temp@meta.data[img_rows, "smoothed_label_s4"]
}
dev.off()
write.csv(scan_params_all, "layer_smoothing_s4_L1-6WM_all_scan_params_11302023.csv")

pdf("layer_smoothing_s4_L1-6WM_all_withNA_11302023.pdf")
for( img in names(comb_seurat@images) ) {
    meta_row <- rownames(comb_seurat@images[[img]]@coordinates)[1]
    sample_name_filt <- comb_seurat@meta.data[meta_row, "sample_name"]
    diagnosis_filt <- comb_seurat@meta.data[meta_row, "diagnosis"]

    sp_pt <- SpatialPlot(comb_seurat, group.by="smoothed_label_s4", image.alpha=0.3, images=img, stroke=NA) + scale_fill_manual(values=group.colors) + ggtitle(paste(sample_name_filt, diagnosis_filt))
    print(sp_pt)
}
dev.off()


# impute NA values using nearest neighbors
nnbrs <- 18
rbf_gamma <- 0.5 # default 1/(num data dimensions)
cluster_labels <- comb_seurat@meta.data[, c("orig.ident", "smoothed_label_s4")]
for (img in names(comb_seurat@images)) {
    print(paste("(Stage 5) Imputing NA spots: ", img))
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
    # set adjacency of na values to 0
    # grab the N nearest neighbors, zero all others, and normalize those
    na_rows <- which(is.na(coord_table$smoothed_label_s4))
    adjA[, na_rows] <- 0
    for(i in 1:nrow(adjA)) {
        adjA[i, order(adjA[i, ], decreasing=TRUE)[(nnbrs+1):ncol(adjA)]] <- 0
    }
    adjA <- sweep(adjA, 1, rowSums(adjA), FUN="/")
    for(layer in sort(unique(coord_table$smoothed_label_s4))) {
            coord_table[, paste0(gsub(' ', '', layer), "_nbr_score")] <- adjA %*% (coord_table$smoothed_label_s4 %in% c(layer))
    }
    nbr_scores_table <- coord_table[, endsWith(colnames(coord_table), "_nbr_score")]
    # NOTE: Taking random assignments in cases of ties so as not to systematically bias one layer over another. Make sure to set a seed.
    coord_table$imputed_NN_layer <- colnames(nbr_scores_table)[max.col(nbr_scores_table,ties.method="random")]
    coord_table$imputed_NN_layer <- gsub("_nbr_score", "", coord_table$imputed_NN_layer)
    coord_table$imputed_NN_layer <- gsub("Layer", "Layer ", coord_table$imputed_NN_layer)
    coord_table$smoothed_label_s5 <- coord_table$smoothed_label_s4
    coord_table[is.na(coord_table$smoothed_label_s5), "smoothed_label_s5"] <- coord_table[is.na(coord_table$smoothed_label_s5), "imputed_NN_layer"]
    comb_seurat@meta.data[rownames(coord_table), "smoothed_label_s5"] <- coord_table$smoothed_label_s5
}
pdf("layer_smoothing_s5_L1-6WM_all_NAimputed_11302023.pdf")
for(img in names(comb_seurat@images)) {
    print(img)
    meta_row <- rownames(comb_seurat@images[[img]]@coordinates)[1]
    sample_name_filt <- comb_seurat@meta.data[meta_row, "orig.ident"]
    diagnosis_filt <- comb_seurat@meta.data[meta_row, "diagnosis"]
    
    sp_pt <- SpatialPlot(comb_seurat, group.by="smoothed_label_s5", image.alpha=0.3, images=img, stroke=NA) + scale_fill_manual(values=group.colors) + ggtitle(paste(sample_name_filt, diagnosis_filt))
    print(sp_pt)
}
dev.off()

saveRDS(comb_seurat, "comb_seurat_post_smoothed_label_s5_11302023.rds")

# ------
