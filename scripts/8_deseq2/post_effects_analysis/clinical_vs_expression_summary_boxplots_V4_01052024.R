library(stringr)
library(ggplot2)
library(dplyr)
library(matrixStats)

clinical_data <- read.csv("/Users/jy1008/Downloads/all_samples_meta_data_V4_01042024.csv")
clinical_data$X <- NULL
# in nctx_temporal, there is only 1 sample with a score of 4: set this to 3 so as not to overly bias the regression
clinical_data$nctx_temporal_fix <- clinical_data$nctx_temporal
clinical_data$nctx_temporal_fix[clinical_data$nctx_temporal_fix > 3] <- 3

agg_counts <- list()
agg_counts_norm <- list() # normalized by colsum
agg_count_dir <- "/Users/jy1008/Downloads/pseudobulk_DE_DESeq2_V4_12052023/agg_counts"
layers <- c(paste("Layer", 1:6), "WM")
for(layer in layers) {
    layer_nospace <- gsub(' ', '', layer)
    agg_count_csv <- file.path(agg_count_dir, sprintf("agg_counts_filtered_smoothed_label_s5_%s_V4_12052023.csv", layer_nospace))
    agg_count_df <- read.csv(agg_count_csv)
    colnames(agg_count_df)[colnames(agg_count_df) == "X"] <- "gene"
    agg_counts[[layer]] <- agg_count_df

    count_scaling_factor <- 10^4
    agg_count_dfnorm <- agg_count_df
    temp_geneonly <- agg_count_df[, 2:ncol(agg_count_df)]
    agg_count_dfnorm[, 2:ncol(agg_count_dfnorm)] <- sweep(temp_geneonly, 2, colSums(temp_geneonly),`/`) * count_scaling_factor
    agg_counts_norm[[layer]] <- agg_count_dfnorm
}







# For one pathway and one count matrix, return average expression for leading edges of the pathway
# Plotted on the y-axis should be the normalized counts prior to log transformation
get_leading_edge_ave <- function(agg_counts_norm, leading_edges, pathway) {
    layers <- sort(unique(leading_edges$layer))
    scores_per_layer <- c()
    for(layer in layers) {
        genes <- leading_edges[leading_edges$pathway == pathway & leading_edges$layer == layer, "leading_edge"]
        if(length(genes) > 0) {
            norm_counts <- agg_counts_norm[[layer]]
            norm_counts <- norm_counts[norm_counts$gene %in% genes, ]
            
            scores <- data.frame(mean_expr = colMeans(norm_counts[, 2:ncol(norm_counts)]),
                                 layer = layer)
            scores$sample_name <- rownames(scores)
            scores_per_layer <- rbind(scores_per_layer, scores)
        }
    }
    return(scores_per_layer)
}

plot_leading_edge_expr <- function(clinical_feature, pathway_db, pathway) {
    # use leading edge matrix to generate mean expression scores for pathways
    # leading_edges <- readRDS("gsea_leading_edges_nctx_temporal_04172023.rds")
    gsea_dir <- "/Users/jy1008/Downloads/pseudobulk_DE_GSEA_V4_12292023"
    gsea_dir2 <- sprintf("gsea_clustering_%s_%s", pathway_db, clinical_feature)
    leading_edge_file <- sprintf("gsea_leading_edges_%s_%s.rds", pathway_db, clinical_feature)
    leading_edges <- readRDS(file.path(gsea_dir, gsea_dir2, leading_edge_file))
    leading_edges$layer <- gsub('_', ' ', leading_edges$layer)

    scores_per_layer <- get_leading_edge_ave(agg_counts_norm, leading_edges=leading_edges, pathway=pathway)
    scores_per_layer <- merge(scores_per_layer, clinical_data[, c("sample_name", clinical_feature)], by="sample_name")
    scores_per_layer <- scores_per_layer[!is.na(scores_per_layer[, clinical_feature]), ]

    scores_per_layer_ave <- scores_per_layer %>% group_by_at(c(clinical_feature, "layer")) %>% dplyr::summarize(across(c(mean_expr), list(mean = mean, sd = sd, se=~sd(.)/sqrt(n())), na.rm = TRUE)) %>% as.data.frame
    scores_per_layer_ave[, "clinical_feature_var"] <- scores_per_layer_ave[, clinical_feature] # for ggplot

    scores_per_layer[, "clinical_feature_var"] <- scores_per_layer[, clinical_feature] # for ggplot

    outname <- paste(clinical_feature, pathway, "01122024", sep='_')

    set.seed(101)
    pdf(sprintf("gsea_pathway_pseudoscores_%s.pdf", outname), width=4, height=4)
    # standard deviation
    pt <- ggplot(scores_per_layer_ave, aes(x=clinical_feature_var, y=mean_expr_mean, group=layer,  color=layer)) + geom_line(position=position_dodge(0.1), size=1) + geom_point(position=position_dodge(0.1)) +
        geom_errorbar(aes(ymin=mean_expr_mean-mean_expr_sd, ymax=mean_expr_mean+mean_expr_sd), width=.4, size=1,
            position=position_dodge(0.1)) + xlab(clinical_feature) + ylab("Mean scaled expression of leading edges") + theme_bw()
    pt <- pt + geom_jitter(data=scores_per_layer, aes(x=clinical_feature_var, y=mean_expr, group=layer, color=layer), width=0.1, height=0, size=0.5)
    print(pt)
    dev.off()
}

clinical_feature <- "nctx_temporal_fix"
pathway_db <- "kegg"
# pathway <- "KEGG_OXIDATIVE_PHOSPHORYLATION"
# pathway <- "KEGG_RIBOSOME"
# pathway <- "KEGG_MAPK_SIGNALING_PATHWAY"
pathway <- "KEGG_LYSOSOME"

plot_leading_edge_expr(clinical_feature, pathway_db, pathway)












all_pathways <- c()
indiv_scores <- c()
for(pathway in pathways) {
    comb_pt_data_mult <- c()
    for(layer in layers_to_test) {
        genes <- leading_edges[leading_edges$pathway == pathway & leading_edges$layer == layer, "leading_edge"]
        
        layer_cols <- colnames(st_countnorm)[grepl(gsub('_', '.', layer), colnames(st_countnorm))]

        st_countnorm_sub <- st_countnorm[genes, layer_cols]
        colnames(st_countnorm_sub) <- str_split_fixed(colnames(st_countnorm_sub), '_', n=2)[, 1]

        
        sample_cols <- colnames(st_countnorm_sub)
        sample_meta <- cbind(clinical_data[sample_cols, c(clinical_features, "Unified.LB.Stage.int")],
                             st_input_files[sample_cols, c("batch", "PMI", "RIN", "sex", "expired_age", "diagnosis", "diagnosis.int", "diagnosis.cont.ilbd", "diagnosis.cont.case")])
        
        # there is a single sample with nctx_temporal = 4: fold this into 3
        sample_meta$nctx_temporal_orig <- sample_meta$nctx_temporal
        sample_meta$nctx_temporal[sample_meta$nctx_temporal > 3] <- 3

        
        st_countnorm_sub <- t(st_countnorm_sub)
        st_countnorm_sub <- as.data.frame(st_countnorm_sub)
        st_countnorm_sub$CaseID <- rownames(st_countnorm_sub)
        # arithmetic mean: the data is counts scaled by total count
        st_countnorm_sub$leading_edge_ave <- rowMeans( st_countnorm_sub[, 1:(ncol(st_countnorm_sub)-1)] )
        
        st_countnorm_sub <- st_countnorm_sub[, c("CaseID", "leading_edge_ave")]
        comb_plot_data <- merge(sample_meta, st_countnorm_sub, by="CaseID")
        comb_plot_data$nctx_temporal <- as.factor(comb_plot_data$nctx_temporal)
        comb_plot_data$layer <- layer

        comb_pt_data_mult <- rbind(comb_pt_data_mult, comb_plot_data)
    }

    comb_pt_data_mult <- comb_pt_data_mult[!is.na(comb_pt_data_mult[, clinical_feature]), ]

    comb_pt_data_mult <- comb_pt_data_mult[, c("leading_edge_ave", clinical_feature, "layer")]
    comb_pt_data_mult_ave <- comb_pt_data_mult %>% group_by_at(c(clinical_feature, "layer")) %>% dplyr::summarize(across(c(leading_edge_ave), list(mean = mean, sd = sd, se=~sd(.)/sqrt(n())), na.rm = TRUE)) %>% as.data.frame
    
    comb_pt_data_mult_ave$leading_edge_ave_sd[is.na(comb_pt_data_mult_ave$leading_edge_ave_sd)] <- 0
    comb_pt_data_mult_ave$leading_edge_ave_se[is.na(comb_pt_data_mult_ave$leading_edge_ave_se)] <- 0
    
    comb_pt_data_mult_ave$pathway <- pathway
    all_pathways <- rbind(all_pathways, comb_pt_data_mult_ave)
    
    comb_pt_data_mult$pathway <- pathway
    indiv_scores <- rbind(indiv_scores, comb_pt_data_mult)
}




set.seed(101)
pdf(sprintf("gsea_pathway_pseudoscores_%s_ypoints_05162023.pdf", outname), width=4, height=4)
# standard deviation
pt <- ggplot(all_pathways, aes(x=nctx_temporal, y=leading_edge_ave_mean, group=interaction(layer, pathway),  color=layer, linetype=pathway)) + geom_line(position=position_dodge(0.1), size=1) + geom_point(position=position_dodge(0.1)) +
    geom_errorbar(aes(ymin=leading_edge_ave_mean-leading_edge_ave_sd, ymax=leading_edge_ave_mean+leading_edge_ave_sd), width=.4, size=1,
        position=position_dodge(0.1)) + xlab("MTG Lewy body load") + ylab("Mean scaled expression of leading edges") + theme_bw() + theme(legend.position = "none")
pt <- pt + geom_jitter(data=indiv_scores, aes(x=nctx_temporal, y=leading_edge_ave, group=interaction(layer, pathway), color=layer), width=0.1, height=0, size=0.5)
print(pt)
dev.off()
# pdf(sprintf("gsea_pathway_pseudoscores_%s_npoints_05162023.pdf", outname), width=8, height=4)
# # standard deviation
# pt <- ggplot(all_pathways, aes(x=nctx_temporal, y=leading_edge_ave_mean, group=interaction(layer, pathway),  color=layer, linetype=pathway)) + geom_line(position=position_dodge(0.1)) + geom_point(position=position_dodge(0.1)) +
#     geom_errorbar(aes(ymin=leading_edge_ave_mean-leading_edge_ave_sd, ymax=leading_edge_ave_mean+leading_edge_ave_sd), width=.4,
#         position=position_dodge(0.1)) + xlab("MTG Lewy body load") + ylab("GeomMean log2norm expression of leading edges") + theme_bw()
# print(pt)
# dev.off()
# pdf(sprintf("gsea_pathway_pseudoscores_%s_se_05162023.pdf", outname), width=8, height=4)
# # standard error
# pt <- ggplot(all_pathways, aes(x=nctx_temporal, y=leading_edge_ave_mean, group=interaction(layer, pathway),  color=layer, linetype=pathway)) + geom_line(position=position_dodge(0.1)) + geom_point(position=position_dodge(0.1)) +
#     geom_errorbar(aes(ymin=leading_edge_ave_mean-leading_edge_ave_se, ymax=leading_edge_ave_mean+leading_edge_ave_se), width=.4,
#         position=position_dodge(0.1)) + xlab("MTG Lewy body load") + ylab("GeomMean log2norm expression of leading edges") + theme_bw()
# print(pt)
# dev.off()



# scatterplot of score vs motor UPDRS score and last MMSE score
# leading_edges <- readRDS("/Users/jy1008/Downloads/gsea_leiden_analysis_kegg_motor_updrs_score_04232023/gsea_leading_edges_kegg_motor_updrs_score_04172023.rds")

clinical_feature <- "motor_updrs_score"
# clinical_feature <- "last_mmse_test_score"
layers_to_test <- c("Layer_5", "Layer_6")
# pathways <- c("KEGG_OXIDATIVE_PHOSPHORYLATION", "KEGG_PARKINSONS_DISEASE")
pathways <- c("KEGG_OXIDATIVE_PHOSPHORYLATION")
outname <- "kegg_oxphos"
outname <- paste(clinical_feature, outname, sep='_')

# all_pathways <- c()
# for(pathway in pathways) {
#     pt <- NULL
#     comb_pt_data_mult <- c()
#     for(layer in layers_to_test) {
#         genes <- leading_edges[leading_edges$pathway == pathway & leading_edges$layer == layer, "leading_edge"]
#         if(length(genes) == 0) {
#             print(paste("no leading edge genes:", clinical_feature, pathway, layer))
#             next
#         }
#         layer_cols <- colnames(st_lognorm)[grepl(gsub('_', '.', layer), colnames(st_lognorm))]
#
#         st_lognorm_sub <- st_lognorm[genes, layer_cols]
#         colnames(st_lognorm_sub) <- str_split_fixed(colnames(st_lognorm_sub), '_', n=2)[, 1]
#
#
#         sample_cols <- colnames(st_lognorm_sub)
#         sample_meta <- cbind(clinical_data[sample_cols, c(clinical_features, "Unified.LB.Stage.int")],
#                              st_input_files[sample_cols, c("batch", "PMI", "RIN", "sex", "expired_age", "diagnosis", "diagnosis.int", "diagnosis.cont.ilbd", "diagnosis.cont.case")])
#
#         # there is a single sample with nctx_temporal = 4: fold this into 3
#         sample_meta$nctx_temporal_orig <- sample_meta$nctx_temporal
#         sample_meta$nctx_temporal[sample_meta$nctx_temporal > 3] <- 3
#
#
#         st_lognorm_sub <- t(st_lognorm_sub)
#         st_lognorm_sub <- as.data.frame(st_lognorm_sub)
#         st_lognorm_sub$CaseID <- rownames(st_lognorm_sub)
#         # st_lognorm_sub$leading_edge_ave <- rowMeans(st_lognorm_sub[, 1:(ncol(st_lognorm_sub)-1)])
#         # geometric
#         st_lognorm_sub$leading_edge_ave <- exp(rowMeans(st_lognorm_sub[, 1:(ncol(st_lognorm_sub)-1)]))
#         # this is the std per subject (for scatterplot), not per subject group (boxplot)
#         # this should be switched to geometric, but also it's unused
#         st_lognorm_sub$leading_edge_std <- rowSds(as.matrix(st_lognorm_sub[, 1:(ncol(st_lognorm_sub)-2)]))
#
#         st_lognorm_sub <- st_lognorm_sub[, c("CaseID", "leading_edge_ave", "leading_edge_std")]
#         comb_plot_data <- merge(sample_meta, st_lognorm_sub, by="CaseID")
#         comb_plot_data$nctx_temporal <- as.factor(comb_plot_data$nctx_temporal)
#         comb_plot_data$layer <- layer
#
#         # linear regression of average lognorm expression against clinical feature
#         comb_plot_data_temp <- comb_plot_data
#         # comb_plot_data_temp$nctx_temporal <- as.numeric(comb_plot_data_temp$nctx_temporal)
#         comb_plot_data_temp[, clinical_feature] <- as.numeric(comb_plot_data_temp[, clinical_feature])
#         comb_plot_data_temp[, clinical_feature] <- scale(comb_plot_data_temp[, clinical_feature],  center=TRUE, scale=TRUE)
#         comb_plot_data_temp$expired_age <- scale(comb_plot_data_temp$expired_age, center=TRUE, scale=TRUE)
#         comb_plot_data_temp$PMI <- scale(comb_plot_data_temp$PMI, center=TRUE, scale=TRUE)
#         comb_plot_data_temp$RIN <- scale(comb_plot_data_temp$RIN, center=TRUE, scale=TRUE)
#         comb_plot_data_temp$leading_edge_ave <- scale(comb_plot_data_temp$leading_edge_ave, center=TRUE, scale=TRUE)
#
#         lmd <- lm(as.formula(paste0("leading_edge_ave ~ batch + sex + expired_age + PMI + RIN + ", clinical_feature)), data=comb_plot_data_temp)
#         print(paste("pathway:", pathway, "layer:", layer, "feature:", clinical_feature))
#         print(paste("lm estimate:", summary(lmd)$coefficients[clinical_feature, "Estimate"]))
#         print(paste("lm Std. Error:", summary(lmd)$coefficients[clinical_feature, "Std. Error"]))
#         print(paste("lm Pr(>|t|):", summary(lmd)$coefficients[clinical_feature, "Pr(>|t|)"]))
#
#         comb_pt_data_mult <- rbind(comb_pt_data_mult, comb_plot_data)
#     }
#     comb_pt_data_mult <- comb_pt_data_mult[!is.na(comb_pt_data_mult[, clinical_feature]), ]
#     comb_pt_data_mult <- comb_pt_data_mult[, c("leading_edge_ave", "leading_edge_std", clinical_feature, "layer", "nctx_temporal")]
#     comb_pt_data_mult$pathway <- pathway
#     all_pathways <- rbind(all_pathways, comb_pt_data_mult)
# }


all_pathways <- c()
for(pathway in pathways) {
    pt <- NULL
    comb_pt_data_mult <- c()
    for(layer in layers_to_test) {
        genes <- leading_edges[leading_edges$pathway == pathway & leading_edges$layer == layer, "leading_edge"]
        if(length(genes) == 0) {
            print(paste("no leading edge genes:", clinical_feature, pathway, layer))
            next
        }
        layer_cols <- colnames(st_countnorm)[grepl(gsub('_', '.', layer), colnames(st_countnorm))]

        st_countnorm_sub <- st_countnorm[genes, layer_cols]
        colnames(st_countnorm_sub) <- str_split_fixed(colnames(st_countnorm_sub), '_', n=2)[, 1]

        
        sample_cols <- colnames(st_countnorm_sub)
        sample_meta <- cbind(clinical_data[sample_cols, c(clinical_features, "Unified.LB.Stage.int")],
                             st_input_files[sample_cols, c("batch", "PMI", "RIN", "sex", "expired_age", "diagnosis", "diagnosis.int", "diagnosis.cont.ilbd", "diagnosis.cont.case")])
        
        # there is a single sample with nctx_temporal = 4: fold this into 3
        sample_meta$nctx_temporal_orig <- sample_meta$nctx_temporal
        sample_meta$nctx_temporal[sample_meta$nctx_temporal > 3] <- 3
        
        
        st_countnorm_sub <- t(st_countnorm_sub)
        st_countnorm_sub <- as.data.frame(st_countnorm_sub)
        st_countnorm_sub$CaseID <- rownames(st_countnorm_sub)
        # geometric
        # st_countnorm_sub$leading_edge_ave <- exp(rowMeans(st_countnorm_sub[, 1:(ncol(st_countnorm_sub)-1)]))
        # arithmetic
        st_countnorm_sub$leading_edge_ave <- rowMeans(st_countnorm_sub[, 1:(ncol(st_countnorm_sub)-1)])
        # this is the std per subject (for scatterplot), not per subject group (boxplot)
        st_countnorm_sub$leading_edge_std <- rowSds(as.matrix(st_countnorm_sub[, 1:(ncol(st_countnorm_sub)-2)]))
        
        st_countnorm_sub <- st_countnorm_sub[, c("CaseID", "leading_edge_ave", "leading_edge_std")]
        comb_plot_data <- merge(sample_meta, st_countnorm_sub, by="CaseID")
        comb_plot_data$nctx_temporal <- as.factor(comb_plot_data$nctx_temporal)
        comb_plot_data$layer <- layer

        # linear regression of average lognorm expression against clinical feature
        comb_plot_data_temp <- comb_plot_data
        # comb_plot_data_temp$nctx_temporal <- as.numeric(comb_plot_data_temp$nctx_temporal)
        comb_plot_data_temp[, clinical_feature] <- as.numeric(comb_plot_data_temp[, clinical_feature])
        comb_plot_data_temp[, clinical_feature] <- scale(comb_plot_data_temp[, clinical_feature],  center=TRUE, scale=TRUE)
        comb_plot_data_temp$expired_age <- scale(comb_plot_data_temp$expired_age, center=TRUE, scale=TRUE)
        comb_plot_data_temp$PMI <- scale(comb_plot_data_temp$PMI, center=TRUE, scale=TRUE)
        comb_plot_data_temp$RIN <- scale(comb_plot_data_temp$RIN, center=TRUE, scale=TRUE)
        comb_plot_data_temp$leading_edge_ave <- scale(comb_plot_data_temp$leading_edge_ave, center=TRUE, scale=TRUE)
        
        lmd <- lm(as.formula(paste0("leading_edge_ave ~ batch + sex + expired_age + PMI + RIN + ", clinical_feature)), data=comb_plot_data_temp)
        print(paste("pathway:", pathway, "layer:", layer, "feature:", clinical_feature))
        print(paste("lm estimate:", summary(lmd)$coefficients[clinical_feature, "Estimate"]))
        print(paste("lm Std. Error:", summary(lmd)$coefficients[clinical_feature, "Std. Error"]))
        print(paste("lm Pr(>|t|):", summary(lmd)$coefficients[clinical_feature, "Pr(>|t|)"]))
        
        comb_pt_data_mult <- rbind(comb_pt_data_mult, comb_plot_data)
    }
    comb_pt_data_mult <- comb_pt_data_mult[!is.na(comb_pt_data_mult[, clinical_feature]), ]
    comb_pt_data_mult <- comb_pt_data_mult[, c("leading_edge_ave", "leading_edge_std", clinical_feature, "layer", "nctx_temporal")]
    comb_pt_data_mult$pathway <- pathway
    all_pathways <- rbind(all_pathways, comb_pt_data_mult)
}




# pdf(sprintf("gsea_pathway_pseudoscores_%s_05122023.pdf", outname))
# for(layer in unique(all_pathways$layer)) {
#     for(pathway in unique(all_pathways$pathway)) {
#         all_pathways_sub <- all_pathways[all_pathways$layer == layer & all_pathways$pathway == pathway, ]
#         print(paste(layer, pathway))
#         dt_pt <- ggplot(all_pathways_sub, aes_string(x=clinical_feature, y="leading_edge_ave")) + geom_point(aes(col=nctx_temporal)) +
#             geom_smooth(method = "lm", color="black") +
#             scale_color_manual(values = c("#984ea3", "#2b83ba", "#abdda4", "#fdae61", "#d7191c")) +
#             # geom_errorbar(aes(ymin=leading_edge_ave - leading_edge_std, ymax=leading_edge_ave + leading_edge_std), width=.4, position=position_dodge(0.1)) +
#             xlab("Motor UPDRS score") + ylab("Mean log2norm expression of leading edges") + theme_bw() + ggtitle(paste(layer, pathway))
#         print(dt_pt)
#     }
# }
# dev.off()

# layers 5 and 6 combined
# https://stackoverflow.com/a/57989074
pdf(sprintf("gsea_pathway_pseudoscores_L5L6combined_%s_05162023.pdf", outname), width=6, height=4)
layer <- "Layer_6"
all_pathways_sub <- all_pathways[all_pathways$layer == layer & all_pathways$pathway == pathway, ]
dt_pt <- ggplot(all_pathways_sub, aes_string(x=clinical_feature, y="leading_edge_ave")) +
                    geom_point(aes(col=nctx_temporal), stroke=0) +
                    geom_smooth(method = "lm", color="black", fill=NA) +
                    scale_color_manual(values = c("#984ea3", "#2b83ba", "#abdda4", "#fdae61", "#d7191c")) +
                    # geom_errorbar(aes(ymin=leading_edge_ave - leading_edge_std, ymax=leading_edge_ave + leading_edge_std), width=.4, position=position_dodge(0.1)) +
                    xlab(clinical_feature) + ylab("Mean scaled expression of leading edges") + theme_bw() + ggtitle(paste(layer, pathway))
ggg <- ggplot_build(dt_pt)
dat_temp <- ggg$data[[2]]
dt_pt <- dt_pt + geom_line(aes(x, ymin), data=dat_temp, linetype="dashed") +
  geom_line(aes(x, ymax), data=dat_temp, linetype="dashed")

layer <- "Layer_5"
all_pathways_sub <- all_pathways[all_pathways$layer == layer & all_pathways$pathway == pathway, ]
dt_pt <- dt_pt + geom_smooth(data = all_pathways_sub, method = "lm", color="gray", fill=NA)
ggg <- ggplot_build(dt_pt)
dat_temp <- ggg$data[[5]]
dt_pt <- dt_pt + geom_line(aes(x, ymin), data=dat_temp, linetype="dashed", color="gray") +
  geom_line(aes(x, ymax), data=dat_temp, linetype="dashed", color="gray")
print(dt_pt)
dev.off()
