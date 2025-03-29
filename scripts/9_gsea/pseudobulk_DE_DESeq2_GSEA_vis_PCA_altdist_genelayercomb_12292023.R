# see "pseudobulk_DE_DESeq2_post_analysis_subcats_02142023.R

# https://lashlock.github.io/compbio/R_presentation.html
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#differential-expression-analysis

library(Seurat)
# library(DESeq2)
library(stringr)
library(testit)
library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(cluster)

library(ape)
library(phytools)
library(magrittr)
library(ggtreeExtra)
library(ggtree)
library(treeio)
library(tidytree)
library(ggstar)
library(ggnewscale)
source("/Users/jy1008/Downloads/pseudobulk_DE_DESeq2_GSEA_12292023_utils.R")
source("/Users/jy1008/Downloads/gheatmap_custom.R")

prep_gsea_df <- function(combined_gt, pathway_db, num_sig_tests_cutoff = 1, nes_cutoff = 0, goslim_cutoff=0.7) {
    # create dotplot across layers
    # ignore NA p-values: these are produced by fgsea when estimation is poor
    combined_gt <- combined_gt[!is.na(combined_gt$pval), ]
    
    # keep only genes which are significant in more than one category?
    # num_sig_tests_cutoff <- 1
    combined_gt_filt <- combined_gt[combined_gt$padj < 0.05, ]
    pathway_sig_counts <- table(combined_gt_filt$pathway)
    sig_pathways <- names(pathway_sig_counts)[pathway_sig_counts >= num_sig_tests_cutoff]
    combined_gt <- combined_gt[combined_gt$pathway %in% sig_pathways, ]
    
    # keep only genes passing log fold change filter
    # nes_cutoff <- 0
    combined_gt_filt <- combined_gt[abs(combined_gt$NES) > nes_cutoff, ]
    sig_pathways <- unique(combined_gt_filt$pathway)
    combined_gt <- combined_gt[combined_gt$pathway %in% sig_pathways, ]

    # -log10(pval=0) = Inf, which doesn't show up in dotplot
    # set Inf values to equal the highest finite value
    combined_gt$nlog10_pval <- -log10(combined_gt$pval)
    combined_gt$nlog10_pval[is.infinite(combined_gt$nlog10_pval)] <- max(combined_gt$nlog10_pval[is.finite(combined_gt$nlog10_pval)])
    
    
    # run GO slim on filtered pathways
    # https://rdrr.io/bioc/clusterProfiler/src/R/simplify.R
    if(grepl("GO", pathway_db, fixed=TRUE)) {
        ont <- gsub("^GO-", "", gsub("-slim$", "", pathway_db))
        semData <- GOSemSim::godata(ont = ont)
        pathway_id_list <- sort(unique(combined_gt$pathway_id))
        print(sprintf("calculating GO slim on %s with %s pathways", ont, length(pathway_id_list)))
        sim <- GOSemSim::mgoSim(pathway_id_list, pathway_id_list, semData = semData,
                                measure = "Wang", combine = NULL)
        
        ## to satisfy codetools for calling gather
        go1 <- go2 <- similarity <- NULL

        sim.df <- as.data.frame(sim)
        sim.df$go1 <- row.names(sim.df)
        sim.df <- gather(sim.df, go2, similarity, -go1)

        sim.df <- sim.df[!is.na(sim.df$similarity),]
        
        ## feature 'by' is attached to 'go1'
        combined_gt_min_padjs <- combined_gt[, c("pathway_id", "padj")] %>%
                                    group_by(pathway_id) %>%
                                    slice(which.min(padj)) %>%
                                    as.data.frame()
        sim.df <- merge(sim.df, combined_gt_min_padjs, by.x="go1", by.y="pathway_id")
        sim.df$go2 <- as.character(sim.df$go2)
        
        ID <- pathway_id_list
        GO_to_remove <- character()
        cutoff <- goslim_cutoff
        # cutoff <- 0.7
        # # cutoff <- 0.6
        for (i in seq_along(ID)) {
            ii <- which(sim.df$go2 == ID[i] & sim.df$similarity > cutoff)
            ## if length(ii) == 1, then go1 == go2
            if (length(ii) < 2)
                next

            sim_subset <- sim.df[ii,]
            # Prevent the following scenario from happening: neither pathway is removed despite similarity
            # go1        go2 similarity       padj
            # 65064   GO:0001539 GO:0001539      1.000 0.03600563
            # 2228763 GO:0060285 GO:0001539      0.877 0.03600563
            # > jj <- which(sim_subset[, "padj"] == min(sim_subset[, "padj"]))
            # > jj
            # [1] 1 2
            jj <- which(sim_subset[, "padj"] == min(sim_subset[, "padj"]))
            if(length(sim_subset$go1[-jj]) == 0) {
                print("---")
                print(sim_subset)
                print(sim_subset$go1[-jj])
                print("---")
            }
            ## sim.df <- sim.df[-ii[-jj]]
            GO_to_remove <- c(GO_to_remove, sim_subset$go1[-jj]) %>% unique
        }
        GO_to_keep <- setdiff(pathway_id_list, GO_to_remove)
        combined_gt <- combined_gt[combined_gt$pathway_id %in% GO_to_keep, ]
    }
    return(combined_gt)
}


#####################################################
# Create dotplot with leading edge-based clustering #
#####################################################


pathway_hclust <- function(combined_gt) {

    # generate distance matrix across pathways
    tot_distances <- 0
    pathway_order <- sort(unique(combined_gt$pathway))
    for(comp in sort(unique(combined_gt$comparison))) {
        leading_edges <- data.frame(matrix(ncol = 3, nrow = 0))
        combined_gt_sub <- combined_gt[combined_gt$comparison == comp, ]
        for(i in 1:nrow(combined_gt_sub)) {
            gsea_row <- combined_gt_sub[i, ]
            leading_edge_set <- strsplit(gsea_row$leadingEdge, ', *')[[1]]
            # print(paste(i, "/", nrow(combined_gt_sub),
            #             gsea_row$pathway, length(leading_edge_set)))
            for(gene in leading_edge_set) {
                new_row <- data.frame(
                                pathway = gsea_row$pathway,
                                leading_edge = gene,
                                NES_sign = sign(gsea_row$NES))
                leading_edges <- rbind(leading_edges, new_row)
            }
        }
        leading_edge_grid <- as.data.frame(pivot_wider(leading_edges, names_from = leading_edge, values_from = NES_sign, values_fill = 0))
        for(missing_pathway in setdiff(pathway_order, leading_edge_grid$pathway)) {
            print(paste("adding missing pathway", missing_pathway))
            leading_edge_grid <- rbind(leading_edge_grid,
                                    c(missing_pathway, rep( 0, (ncol(leading_edge_grid)-1) )))
            # Note that the 0s are added as characters: reconvert all columns to numeric afterward
        }
        rownames(leading_edge_grid) <- leading_edge_grid$pathway
        leading_edge_grid$pathway <- NULL
        # convert all columns to numeric
        # https://stackoverflow.com/a/63365119
        leading_edge_grid[] <- sapply(leading_edge_grid, as.numeric)
        
        # confirm the pathway order is the same across all comparisons
        leading_edge_grid <- leading_edge_grid[sort(rownames(leading_edge_grid)), ]
        assert(all(rownames(leading_edge_grid) == pathway_order))

        # remove columns with zero variance prior to PCA
        leading_edge_grid <- leading_edge_grid[ , which(apply(leading_edge_grid, 2, var) != 0)]
        
        # apply PCA before calculating distance
        num_pc_components <- 30
        leading_edge_pca <- prcomp(leading_edge_grid, scale=TRUE)
        leading_edge_pca <- leading_edge_pca$x
        if(ncol(leading_edge_pca) > num_pc_components) {
            leading_edge_pca <- leading_edge_pca[, 1:num_pc_components]
        }
        temp_dist <- dist(leading_edge_pca, method="canberra")
        # temp_dist <- dist(leading_edge_grid, method="canberra")
        # temp_dist <- dist(leading_edge_pca, method="euclidean")
        # temp_dist <- dist(leading_edge_grid, method="euclidean")
        tot_distances <- tot_distances + temp_dist
        # tot_distances <- tot_distances + temp_dist**2
    }
    # tot_distances <- sqrt(tot_distances)
    
    # append layer-wise NES score for the current layer
    # gsea_subset_temp <- combined_gt[, c("pathway", "comparison", "NES")]
    # gsea_subset_temp <- gsea_subset_temp %>% pivot_wider(names_from = comparison, values_from = NES, values_fill = 0) %>% as.data.frame()
    gsea_subset_temp <- combined_gt[, c("pathway", "comparison", "NES", "padj")]
    gsea_subset_temp$NES_sig <- sign(gsea_subset_temp$NES) * (gsea_subset_temp$padj < 0.05)
    gsea_subset_temp <- gsea_subset_temp[, c("pathway", "comparison", "NES_sig")]
    gsea_subset_temp <- gsea_subset_temp %>% pivot_wider(names_from = comparison, values_from = NES_sig, values_fill = 0) %>% as.data.frame()
    rownames(gsea_subset_temp) <- gsea_subset_temp$pathway
    gsea_subset_temp$pathway <- NULL
    gsea_subset_temp <- gsea_subset_temp[pathway_order, ]
    assert(all(rownames(gsea_subset_temp) == pathway_order))
    nes_dist <- dist(gsea_subset_temp, method="canberra")
    nes_dist_scaling_factor <- 5
    
    tot_distances <- tot_distances + nes_dist_scaling_factor * nes_dist

    # create dist object
    attributes(tot_distances) <- list(Size = length(pathway_order),
                                  Labels = pathway_order,
                                  Diag = FALSE,
                                  Upper = FALSE,
                                  method = "sum canberra L1",
                                  call = 'dist(leading_edge_grid, method="canberra")',
                                  # method = "euclidean",
                                  # call = 'dist(leading_edge_grid, method="euclidean")',
                                  class = "dist")
    # replace NAs in tot_distances with Inf
    tot_dist_temp <- tot_distances
    tot_dist_temp[is.na(tot_dist_temp)] <- max(tot_distances[!is.na(tot_distances)]) + 100
    cluster.row <- hclust(tot_dist_temp, method = 'complete') # the default
    return(list(cluster.row, tot_dist_temp))
}

get_best_cutree <- function(cluster.row, tot_dist_temp, desired_num_k = NA, min_k = 5) {
    # use silhouette score to combine very similar labels
    # tot_dist_temp
    # cluster.row
    sil_widths <- c()
    for(clust_k in min_k:length(cluster.row$labels)) {
        cutree_labels <- cutree(cluster.row, k=clust_k)
        ss <- silhouette(cutree_labels, tot_dist_temp)
        if(length(ss)==1 && is.na(ss)) { # length is to avoid Warning()
            mean_width <- 0
        } else {
            mean_width <- mean(ss[, 3])
        }
        sil_widths <- rbind(sil_widths, data.frame("clust_k" = clust_k, "mean_sil_width" = mean_width))
    }
    # score by cluster
    # https://stackoverflow.com/a/72166839
    # ss[ , ] %>% as.data.frame() %>% group_by(cluster) %>% mutate(mean_clust_sw = mean(sil_width), clust_count = n()) %>% arrange(clust_count) %>% as.data.frame()
    max_clust_ks <- localMaxima(sil_widths[order(sil_widths$clust_k), "mean_sil_width"])
    # max_clust_ks <- sil_widths[max_clust_ks, "mean_sil_width"]
    if(is.na(desired_num_k)) {
        best_clust_k <- sil_widths[sil_widths$mean_sil_width == max(sil_widths$mean_sil_width), "clust_k"]
        best_clust_k <- max(best_clust_k) # in case of ties
    } else {
        # desired_num_k <- 100
        best_clust_k <- max_clust_ks[which( abs(max_clust_ks - desired_num_k) == min(abs(max_clust_ks - desired_num_k)) )]
        # in the case of a tie, get the smaller k
        best_clust_k <- min(best_clust_k)
    }
    best_cutree_labels <- cutree(cluster.row, k=best_clust_k)
    return(list(sil_widths, best_cutree_labels))
}






run_gsea_visualization <- function(pathway_db, clinical_feat, gsea_dir, out_dir, outname, desired_num_k_range=c(20, 50, 100), goslim_cutoff=0.7) {
    # Main function to generate all visualizations

    if(!dir.exists(out_dir)) { dir.create(out_dir) }
    if(!dir.exists(file.path(out_dir, pathway_db))) {dir.create(file.path(out_dir, pathway_db))}
    
    clinical_out_dir <- file.path(out_dir, pathway_db, clinical_feat)
    if(!dir.exists(clinical_out_dir)) {dir.create(clinical_out_dir)}
    setwd(clinical_out_dir)

    gt_file <- sprintf("de_genes_DESeq2_pseudobulk_%s_%s_[0-9]+\\.fGSEA.xls$", clinical_feat, pathway_db)
    file_match <- list.files( path=file.path( gsea_dir, pathway_db), pattern=gt_file )
    if(length(file_match) != 1) {
        print(paste("Failed to read file", gt_file))
        next()
    }

    gt_file <- file_match[1]
    print(paste("Reading file", gt_file))


    # nes_filter_list <- c(0, 0.5, 1.0, 1.5, 2.0)
    # for(nes_filter in nes_filter_list) {
    nes_filter <- 0
    combined_gt <- read.table(file.path(gsea_dir, pathway_db, gt_file), sep='\t', header=TRUE)
    combined_gt <- prep_gsea_df(combined_gt, pathway_db, nes_cutoff = nes_filter, goslim_cutoff=goslim_cutoff)
    print(paste("GO Slim cutoff:", goslim_cutoff))
    print(dim(combined_gt))
    
    # # separate pathways of +/- NES, averaging across all padj sig tests
    # comb_gt_signs <- combined_gt %>% filter(padj < 0.05) %>%
    #                     group_by(pathway) %>%
    #                     mutate(mean_NES = mean(NES),
    #                            mean_NES_sign = sign(mean(NES))) %>%
    #                     as.data.frame()
    # pos_pathways <- unique(comb_gt_signs[comb_gt_signs$mean_NES_sign >= 0, "pathway"])
    # neg_pathways <- unique(comb_gt_signs[comb_gt_signs$mean_NES_sign < 0, "pathway"])
    # assert(length(pos_pathways) + length(neg_pathways) == length(unique(combined_gt$pathway)))

    # combined_gt_pos <- combined_gt[combined_gt$pathway %in% pos_pathways, ]
    # pathway_clust_out <- pathway_hclust(combined_gt_pos)
    # cluster.row <- pathway_clust_out[[1]]
    # tot_dist_temp <- pathway_clust_out[[2]]
    # ggtree_all <- ggtree::ggtree(cluster.row) +
    #                      ggtree::geom_tiplab()
    # pathway_display_order_all <- rev(ggtree::get_taxa_name(ggtree_all))
    # layer_order <- sort(unique(combined_gt_pos$comparison))
    
    # # use silhouette score to combine very similar labels
    # # desired_num_k_range <- c(20, 50, 100, NA)
    # # for(desired_num_k in desired_num_k_range) {
    # temp_res <- get_best_cutree(cluster.row, tot_dist_temp, desired_num_k = NA)
    # sil_widths <- temp_res[[1]]
    # best_cutree_labels <- temp_res[[2]]
    # combined_gt_pos[, paste0("best_cutree_labels_", desired_num_k)] <- best_cutree_labels[combined_gt_pos$pathway]
    # # }
    # # pdf(sprintf("leading_edge_dendro_test_silw_%s_%s_%s.pdf", pathway_db, clinical_feat, outname))
    # sp_pt <- ggplot(sil_widths, aes(x=clust_k, y=mean_sil_width)) + geom_point(size=1, shape=19) + theme_bw() + xlab("cutree number of clusters") + ylab("mean silhouette width")
    # print(sp_pt)
    # # dev.off()
    
    
    # generate distance matrix across pathways
    pathway_clust_out <- pathway_hclust(combined_gt)
    cluster.row <- pathway_clust_out[[1]]
    tot_dist_temp <- pathway_clust_out[[2]]
    # ggtree_all <- ggtree::ggtree(cluster.row, branch.length = "none") +
    ggtree_all <- ggtree::ggtree(cluster.row) +
                         ggtree::geom_tiplab()
    # NOTE: the pathway order is changed by ggtree
    # NOTE: also that the ggtree order is reversed, requiring rev()
    # pathway_display_order_all <- attr(tot_dist_temp, "Labels")[order.dendrogram(as.dendrogram(cluster.row))]
    pathway_display_order_all <- rev(ggtree::get_taxa_name(ggtree_all))
    layer_order <- sort(unique(combined_gt$comparison))


    # use silhouette score to combine very similar labels
    # desired_num_k_range <- c(20, 50, 100, NA)
    desired_num_k_range <- c(desired_num_k_range, NA)
    for(desired_num_k in desired_num_k_range) {
        temp_res <- get_best_cutree(cluster.row, tot_dist_temp, desired_num_k = desired_num_k)
        sil_widths <- temp_res[[1]]
        best_cutree_labels <- temp_res[[2]]
        combined_gt[, paste0("best_cutree_labels_", desired_num_k)] <- best_cutree_labels[combined_gt$pathway]
    }
    print(sil_widths)

    pdf(sprintf("leading_edge_dendro_test_silw_%s_%s_%s.pdf", pathway_db, clinical_feat, outname))
    sp_pt <- ggplot(sil_widths, aes(x=clust_k, y=mean_sil_width)) + geom_point(size=1, shape=19) + theme_bw() + xlab("cutree number of clusters") + ylab("mean silhouette width")
    print(sp_pt)
    dev.off()

    # all pathways, no summarization
    width <- 10
    height <- max(10, length(pathway_display_order_all) / 5)

    print(paste0("number of rows in dotplot: ", length(pathway_display_order_all)))
    pdf(sprintf("leading_edge_dendro_test_%s_%s_%s.pdf", pathway_db, clinical_feat, outname), width=width, height=height)
    dt_pt_all <- ggplot(combined_gt, aes(x = factor(comparison, levels = layer_order), y = factor(pathway, levels = pathway_display_order_all), color = NES, size = nlog10_pval)) +
    geom_point(shape = 21,  stroke = 1, aes(colour = padj < 0.05, fill = NES)) +
    scale_fill_gradient2(low = "#0571b0", mid = "#f7f7f7", high = "#ca0020", midpoint = 0, name = 'NES') +
    scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = rgb(0,0,0,0)), name = "p_adj < 0.05") +
    cowplot::theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          panel.grid = element_line(color = "#dddddd", size = 0.75, linetype = 2)) +
    scale_y_discrete(label = function(x) stringr::str_trunc(x, 60))
    print(dt_pt_all)
    # dendrogram
    # plot(reordered.cluster.row, horiz=TRUE)
    plot(as.dendrogram(cluster.row), horiz=TRUE)
    dev.off()

    # pdf(sprintf("leading_edge_dendro_test_siltiles_%s_%s_%s.pdf", pathway_db, clinical_feat, outname), width=width*3, height=height)
    # # for(cluster_label_to_tile in c("best_cutree_labels", "best_desired_cutree_labels")) {
    # for(desired_num_k in desired_num_k_range) {
    #     combined_gt$clust_plus_size <- combined_gt[, paste0("best_cutree_labels_", desired_num_k)]
    #     combined_gt$clust_plus_size <- factor(combined_gt$clust_plus_size, levels = str_sort(unique(combined_gt$clust_plus_size), numeric=TRUE))
    #     ld_pt <- ggplot(combined_gt, aes(x=1, y=factor(pathway, levels = pathway_display_order_all), fill=clust_plus_size)) +
    #                     geom_tile() + scale_y_discrete(label = function(x) stringr::str_trunc(x, 60)) +
    #                     ggtitle(paste("num pathway groups:", desired_num_k))
    #     print(ld_pt)
    # }
    # combined_gt$clust_plus_size <- NULL
    # dev.off()


    for(desired_num_k in desired_num_k_range) {
        # for cowplot, pick at most 50 clusters
        # desired_num_k <- 50
        # if(length(unique(combined_gt$best_cutree_labels_NA)) < desired_num_k) {
        #     desired_num_k <- "NA"
        # }
        combined_gt$clust_plus_size <- combined_gt[, paste0("best_cutree_labels_", desired_num_k)]
        combined_gt$clust_plus_size <- factor(combined_gt$clust_plus_size, levels = str_sort(unique(combined_gt$clust_plus_size), numeric=TRUE))
        ld_pt <- ggplot(combined_gt, aes(x=1, y=factor(pathway, levels = pathway_display_order_all), fill=clust_plus_size)) +
                        geom_tile() + theme(legend.position = "none") +
                        scale_y_discrete(label = function(x) stringr::str_trunc(x, 5))
        actual_num_k <- length(unique(combined_gt$clust_plus_size))

        # cowplot dotplot with dendrogram
        # https://wilkelab.org/cowplot/articles/aligning_plots.html
        width <- 10
        height <- max(10, length(pathway_display_order_all) / 5)
        pdf(sprintf("dotplot_combined_all_pathways_%s_%s_%sdk-%sak_%s.pdf", pathway_db, clinical_feat, desired_num_k, actual_num_k, outname), width=width, height=height)
        all_pt <- cowplot::plot_grid(ggtree_all, ld_pt, dt_pt_all,
                  nrow = 1, rel_widths = c(0.5, 0.2, 2), align = 'h', axis = 'bt')
        print(all_pt)
        dev.off()
        
        # combine samples falling under the same dendro cluster into a single row
        # use a single pathway to name the cluster
        # combined_gt$clust_plus_size <- combined_gt[, paste0("best_cutree_labels_", desired_num_k)]
        # combined_gt$clust_plus_size <- factor(combined_gt$clust_plus_size, levels = str_sort(unique(combined_gt$clust_plus_size), numeric=TRUE))
        
        combined_gt_temp <- combined_gt %>%
                            group_by(clust_plus_size) %>%
                            mutate(clust_pathway_count = n_distinct(pathway),
                                   clust_pathway_nlog10_sum = sum(nlog10_pval)) %>%
                            as.data.frame()
        # filter by smallest p-value
        comb_gt_summ <- combined_gt_temp %>%
                        group_by(pathway) %>%
                        mutate(nlog10_pval_sum = sum(nlog10_pval)) %>%
                        group_by(clust_plus_size) %>%
                        # filter(pval == min(pval)) %>%
                        filter(nlog10_pval_sum == max(nlog10_pval_sum)) %>%
                        as.data.frame()
        keep_pathways <- comb_gt_summ$pathway
        comb_gt_summ <- combined_gt_temp[combined_gt_temp$pathway %in% keep_pathways, ]

        # redo distances
        # generate distance matrix across pathways
        pathway_clust_out <- pathway_hclust(comb_gt_summ)
        cluster.row.pt2 <- pathway_clust_out[[1]]
        tot_dist_temp <- pathway_clust_out[[2]]
        # ggtree_pt <- ggtree::ggtree(cluster.row.pt2, branch.length = "none") +
        ggtree_pt <- ggtree::ggtree(cluster.row.pt2) +
                             ggtree::geom_tiplab()
        # NOTE: the pathway order is changed by ggtree
        # NOTE: also that the ggtree order is reversed, requiring rev()
        # pathway_display_order_sub <- attr(tot_dist_temp, "Labels")[order.dendrogram(as.dendrogram(cluster.row.pt2))]
        pathway_display_order_sub <- rev(ggtree::get_taxa_name(ggtree_pt))
        layer_order <- sort(unique(comb_gt_summ$comparison))

        dt_pt <- ggplot(comb_gt_summ, aes(x = factor(comparison, levels = layer_order), y = factor(pathway, levels = pathway_display_order_sub), color = NES, size = nlog10_pval)) +
        geom_point(shape = 21,  stroke = 1, aes(colour = padj < 0.05, fill = NES)) +
        scale_fill_gradient2(low = "#0571b0", mid = "#f7f7f7", high = "#ca0020", midpoint = 0, name = 'NES') +
        scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = rgb(0,0,0,0)), name = "p_adj < 0.05") +
        cowplot::theme_cowplot() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              panel.grid = element_line(color = "#dddddd", size = 0.75, linetype = 2)) +
        scale_y_discrete(label = function(x) stringr::str_trunc(x, 60))

        # depict size of each pathway cluster
        ld_pt <- ggplot(comb_gt_summ, aes(x=1, y=factor(pathway, levels = pathway_display_order_sub), size=clust_pathway_count, label=clust_pathway_count)) +
                        geom_point(shape = 21, stroke = 1) +
                        geom_text(hjust=-1, vjust=0, size=3) +
                        cowplot::theme_cowplot() +
                        theme(legend.title=element_blank()) +
                        scale_y_discrete(label = function(x) stringr::str_trunc(x, 5))


        # clustered pathways
        width <- 15
        height <- max(10, length(pathway_display_order_sub) / 5)
        pdf(sprintf("dotplot_combined_mega_pathways_%s_%s_%sdk-%sak_%s.pdf", pathway_db, clinical_feat, desired_num_k, actual_num_k, outname), width=width, height=height)
        # all_pt <- cowplot::plot_grid(ggtree::ggtree(as.dendrogram(cluster.row.pt2)), ld_pt, dt_pt,
        all_pt <- cowplot::plot_grid(ggtree_pt, ld_pt, dt_pt,
                  nrow = 1, rel_widths = c(0.5, 0.5, 2), align = 'h', axis = 'bt')
        print(all_pt)
        dev.off()
        
        
        
        # mega pathways filtered and reclustered by NES (zeroed if not padj < 0.05)
        comp_list <- sort(unique(comb_gt_summ$comparison))
        comb_gt_summ_wide <- comb_gt_summ %>% mutate(NESsig = as.numeric(NES * (padj<0.05))) %>%
                             pivot_wider(names_from = comparison,
                                         values_from = NESsig, values_fill = 0) %>%
                             select(all_of(c("pathway", comp_list))) %>%
                             group_by(pathway) %>%
                             summarise(across(all_of(comp_list), sum, .names = "sum_{col}")) %>%
                             as.data.frame()
        rownames(comb_gt_summ_wide) <- comb_gt_summ_wide$pathway
        comb_gt_summ_wide$pathway <- NULL
        temp_dist <- dist(comb_gt_summ_wide, method="canberra")
        temp_clust <- hclust(temp_dist, method = 'complete') # the default
        temp_ggtree <- ggtree::ggtree(temp_clust) +
                             ggtree::geom_tiplab()
        pathway_display_order_temp <- rev(ggtree::get_taxa_name(temp_ggtree))
        layer_order <- sort(unique(combined_gt$comparison))

        temp_dt_pt <- ggplot(comb_gt_summ, aes(x = factor(comparison, levels = layer_order), y = factor(pathway, levels = pathway_display_order_temp), color = NES, size = nlog10_pval)) +
        geom_point(shape = 21,  stroke = 1, aes(colour = padj < 0.05, fill = NES)) +
        scale_fill_gradient2(low = "#0571b0", mid = "#f7f7f7", high = "#ca0020", midpoint = 0, name = 'NES') +
        scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = rgb(0,0,0,0)), name = "p_adj < 0.05") +
        cowplot::theme_cowplot() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              panel.grid = element_line(color = "#dddddd", size = 0.75, linetype = 2)) +
        scale_y_discrete(label = function(x) stringr::str_trunc(x, 60))

        # depict size of each pathway cluster
        comb_gt_summ$clust_pathway_count_gt1 <- comb_gt_summ$clust_pathway_count
        comb_gt_summ$clust_pathway_count_gt1[comb_gt_summ$clust_pathway_count_gt1 == 1] <- NA
        temp_ld_pt <- ggplot(comb_gt_summ, aes(x=1, y=factor(pathway, levels = pathway_display_order_temp), size=clust_pathway_count_gt1, label=clust_pathway_count_gt1)) +
                        geom_point(shape = 21, stroke = 1) +
                        geom_text(hjust=-1, vjust=0, size=3) +
                        cowplot::theme_cowplot() +
                        theme(legend.title=element_blank()) +
                        scale_y_discrete(label = function(x) stringr::str_trunc(x, 5))


        # clustered pathways
        width <- 15
        height <- max(10, length(pathway_display_order_temp) / 5)
        pdf(sprintf("dotplot_combined_mega_resorted_pathways_%s_%s_%sdk-%sak_%s.pdf", pathway_db, clinical_feat, desired_num_k, actual_num_k, outname), width=width, height=height)
        # all_pt <- cowplot::plot_grid(ggtree::ggtree(as.dendrogram(cluster.row.pt2)), ld_pt, dt_pt,
        all_pt <- cowplot::plot_grid(temp_ggtree, temp_ld_pt, temp_dt_pt,
                  nrow = 1, rel_widths = c(0.5, 0.5, 2), align = 'h', axis = 'bt')
        print(all_pt)
        dev.off()
        
    }
    return(list(combined_gt, cluster.row, pathway_display_order_all))
}


# you may have to run this multiple times to adjust the dimensions of the rings
draw_circle_plot2 <- function(combined_gt, cluster.row, pathway_display_order_all, label_column = "best_cutree_labels_NA") {
    # circular plot
    # ggtree_all
    # pathway_display_order_all
    # save cluster.row into a unique object
    # colors: scales::hue_pal()(300)
    temp_circ_meta <- combined_gt
    temp_circ_meta$pathway <- factor(temp_circ_meta$pathway, levels = pathway_display_order_all)
    temp_circ_meta <- temp_circ_meta[order(temp_circ_meta$pathway), ]
    temp_circ_meta <- temp_circ_meta %>% group_by(pathway) %>% filter(row_number()==1) %>% as.data.frame()

    # colors
    color_list <- scales::hue_pal()(length(unique(temp_circ_meta[, label_column])))
    temp_circ_meta$label_column_color <- color_list[temp_circ_meta[, label_column]]

    cluster.row.ape <- ape::as.phylo(cluster.row)
    temp_circ_meta$nodeid <- nodeid(cluster.row.ape, temp_circ_meta$pathway)

    ggtree_circ <- ggtree::ggtree(cluster.row.ape, layout="fan", open.angle=15, size=0.1)
    for(clust_col in unique(temp_circ_meta$label_column_color)) {
        nodeid_set <- temp_circ_meta[
                    temp_circ_meta$label_column_color==clust_col, "nodeid"]
        rootnode <- phytools::findMRCA(cluster.row.ape, tips=nodeid_set)
        if(is.null(rootnode)) {
            nodedf <- data.frame(node=nodeid_set)
        } else {
            nodedf <- data.frame(node=rootnode)
        }
        ggtree_circ <- ggtree_circ + geom_hilight(data=nodedf, mapping=aes(node=node),
        alpha=0.3, fill=clust_col, size=0.05)
    }

    # cluster.row.ape$tip.label
    combined_gt_temp <- combined_gt
    combined_gt_temp <- combined_gt_temp[, c("pathway", "comparison", "NES")]
    combined_gt_temp <- pivot_wider(combined_gt_temp, names_from = comparison, values_from = NES, values_fill=0) %>% as.data.frame()
    rownames(combined_gt_temp) <- combined_gt_temp$pathway
    combined_gt_temp$pathway <- NULL
    combined_gt_temp <- combined_gt_temp[cluster.row.ape$tip.label, ]
    colnames(combined_gt_temp) <- str_split_fixed(
                                    colnames(combined_gt_temp), " * ", n=3)[, 3]
    # combined_gt_temp$nodeid <- nodeid(cluster.row.ape, combined_gt_temp$pathway)
    ggtree_circ2 <- ggtree_circ + ggtree::geom_tiplab(offset=48, size=5) + new_scale_fill()
    ggtree_circ3 <- gheatmap(ggtree_circ2, combined_gt_temp, offset=-2, width=.3,
             colnames_angle=90, colnames_offset_y = 0.25,
             colnames_position="top", font.size=2) +
        scale_fill_gradient2(low = "#0571b0", mid = "#f7f7f7", high = "#ca0020", midpoint = 0, name = 'NES')
    
    combined_gt_tempb <- combined_gt
    combined_gt_tempb$LE_size <- lengths(str_split(combined_gt_tempb$leadingEdge, ", "))
    combined_gt_tempb <- combined_gt_tempb[, c("pathway", "size", "LE_size")]
    combined_gt_tempb$LE_frac <- combined_gt_tempb$LE_size / combined_gt_tempb$size
    combined_gt_tempb <- combined_gt_tempb %>% group_by(pathway) %>%
                            summarize(mean_size = mean(size, na.rm=TRUE),
                                      mean_LE_frac = mean(LE_frac, na.rm=TRUE),
                                      mean_LE_size = mean(LE_size, na.rm=TRUE)) %>%
                            as.data.frame()
    rownames(combined_gt_tempb) <- combined_gt_tempb$pathway
    combined_gt_tempb <- combined_gt_tempb[cluster.row.ape$tip.label, ]
    ggtree_circ3b <- ggtree_circ3 + geom_fruit(data=combined_gt_tempb, geom=geom_bar,
                                                mapping=aes(y=pathway, x=mean_LE_size),
                                                orientation="y", stat="identity",
                                                offset=0.31, pwidth=0.1)
    # black border for padj < 0.05
    combined_gt_temp2 <- combined_gt
    combined_gt_temp2$padj_sig <- combined_gt_temp2$padj < 0.05
    combined_gt_temp2 <- combined_gt_temp2[, c("pathway", "comparison", "padj_sig")]
    combined_gt_temp2 <- pivot_wider(combined_gt_temp2, names_from = comparison, values_from = padj_sig, values_fill=0) %>% as.data.frame()
    rownames(combined_gt_temp2) <- combined_gt_temp2$pathway
    combined_gt_temp2$pathway <- NULL
    combined_gt_temp2 <- combined_gt_temp2[cluster.row.ape$tip.label, ]
    colnames(combined_gt_temp2) <- str_split_fixed(
                                    colnames(combined_gt_temp2), " * ", n=3)[, 3]
    ggtree_circ4 <- gheatmap(ggtree_circ2, combined_gt_temp2, offset=-2, width=.3,
    colnames_angle=90, colnames_offset_y = 0.25,
    colnames_position="top", font.size=2)
    ggtree_circ4b <- ggtree_circ4 + geom_fruit(data=combined_gt_tempb, geom=geom_bar,
                                                mapping=aes(y=pathway, x=mean_LE_size),
                                                orientation="y", stat="identity",
                                                offset=0.31, pwidth=0.1)
    # ggtree_circ4 <- ggtree_circ3 + ggtree::geom_tiplab(offset=30, size=2)
    # pdf("circular_ggtree_test_04042024.pdf", width=15, height=15)
    desired_num_k <- gsub("^best_cutree_labels_", "", label_column)
    actual_num_k <- length(unique(combined_gt[, label_column]))
    pdf(sprintf("circular_ggtree_%s_%s_%sdk_%sak_%s.pdf", pathway_db, clinical_feat, desired_num_k, actual_num_k, outname), width=15, height=15)
    print(ggtree_circ3b)
    print(ggtree_circ3b + theme(legend.position = "none"))
    print(ggtree_circ4b)
    print(ggtree_circ4b + theme(legend.position = "none"))
    dev.off()
}

draw_circle_plot <- function(combined_gt, cluster.row, pathway_display_order_all, label_column = "best_cutree_labels_NA") {
    # circular plot
    # ggtree_all
    # pathway_display_order_all
    # save cluster.row into a unique object
    # colors: scales::hue_pal()(300)
    temp_circ_meta <- combined_gt
    temp_circ_meta$pathway <- factor(temp_circ_meta$pathway, levels = pathway_display_order_all)
    temp_circ_meta <- temp_circ_meta[order(temp_circ_meta$pathway), ]
    temp_circ_meta <- temp_circ_meta %>% group_by(pathway) %>% filter(row_number()==1) %>% as.data.frame()

    # colors
    color_list <- scales::hue_pal()(length(unique(temp_circ_meta[, label_column])))
    temp_circ_meta$label_column_color <- color_list[temp_circ_meta[, label_column]]

    cluster.row.ape <- ape::as.phylo(cluster.row)
    temp_circ_meta$nodeid <- nodeid(cluster.row.ape, temp_circ_meta$pathway)

    ggtree_circ <- ggtree::ggtree(cluster.row.ape, layout="fan", open.angle=15, size=0.1)
    for(clust_col in unique(temp_circ_meta$label_column_color)) {
        nodeid_set <- temp_circ_meta[
                    temp_circ_meta$label_column_color==clust_col, "nodeid"]
        rootnode <- phytools::findMRCA(cluster.row.ape, tips=nodeid_set)
        if(is.null(rootnode)) {
            nodedf <- data.frame(node=nodeid_set)
        } else {
            nodedf <- data.frame(node=rootnode)
        }
        ggtree_circ <- ggtree_circ + geom_hilight(data=nodedf, mapping=aes(node=node),
        alpha=0.3, fill=clust_col, size=0.05)
    }

    # cluster.row.ape$tip.label
    combined_gt_temp <- combined_gt
    combined_gt_temp <- combined_gt_temp[, c("pathway", "comparison", "NES")]
    combined_gt_temp <- pivot_wider(combined_gt_temp, names_from = comparison, values_from = NES, values_fill=0) %>% as.data.frame()
    rownames(combined_gt_temp) <- combined_gt_temp$pathway
    combined_gt_temp$pathway <- NULL
    combined_gt_temp <- combined_gt_temp[cluster.row.ape$tip.label, ]
    colnames(combined_gt_temp) <- str_split_fixed(
                                    colnames(combined_gt_temp), " * ", n=3)[, 3]
    # combined_gt_temp$nodeid <- nodeid(cluster.row.ape, combined_gt_temp$pathway)
    ggtree_circ2 <- ggtree_circ + ggtree::geom_tiplab(offset=36, size=2) + new_scale_fill()
    ggtree_circ3 <- gheatmap(ggtree_circ2, combined_gt_temp, offset=-2, width=.3,
             colnames_angle=90, colnames_offset_y = 0.25,
             colnames_position="top", font.size=2) +
        scale_fill_gradient2(low = "#0571b0", mid = "#f7f7f7", high = "#ca0020", midpoint = 0, name = 'NES')
    # black border for padj < 0.05
    combined_gt_temp2 <- combined_gt
    combined_gt_temp2$padj_sig <- combined_gt_temp2$padj < 0.05
    combined_gt_temp2 <- combined_gt_temp2[, c("pathway", "comparison", "padj_sig")]
    combined_gt_temp2 <- pivot_wider(combined_gt_temp2, names_from = comparison, values_from = padj_sig, values_fill=0) %>% as.data.frame()
    rownames(combined_gt_temp2) <- combined_gt_temp2$pathway
    combined_gt_temp2$pathway <- NULL
    combined_gt_temp2 <- combined_gt_temp2[cluster.row.ape$tip.label, ]
    colnames(combined_gt_temp2) <- str_split_fixed(
                                    colnames(combined_gt_temp2), " * ", n=3)[, 3]
    ggtree_circ4 <- gheatmap(ggtree_circ2, combined_gt_temp2, offset=-2, width=.3,
    colnames_angle=90, colnames_offset_y = 0.25,
    colnames_position="top", font.size=2)
    # ggtree_circ4 <- ggtree_circ3 + ggtree::geom_tiplab(offset=30, size=2)
    # pdf("circular_ggtree_test_04042024.pdf", width=15, height=15)
    pdf(sprintf("circular_ggtree_%s_%s_%sdk_%sak_%s.pdf", pathway_db, clinical_feat, desired_num_k, actual_num_k, outname), width=15, height=15)
    print(ggtree_circ3)
    print(ggtree_circ3 + theme(legend.position = "none"))
    print(ggtree_circ4)
    print(ggtree_circ4 + theme(legend.position = "none"))
    dev.off()
}



gsea_dir <- "/Users/jy1008/Downloads/pseudobulk_DE_GSEA_V4_01102024/"
# clinical_features <- c("layer_vs", "nctx_temporal_fix", "brain_stem_sn", "diagnosis_cont_ilbd", "diagnosis_cont_case", "last_mmse_test_score", "motor_updrs_score")
clinical_features <- c("nctx_temporal_fix")

# pathway_dbs <- c("c2.cp.kegg.v7.2", "GO-BP-slim", "GO-MF-slim", "GO-CC-slim")

# pathway_dbs <- c("c2.cp.kegg.v7.2", "GO-BP", "GO-MF", "GO-CC")
pathway_dbs <- c("GO-BP")

# only apply group label to clusters with more than N members
# pathway_grp_thresh <- 3
# outname <- "03222024"
# out_dir <- "/Users/jy1008/Downloads/pseudobulk_DE_GSEA_V4_01102024/combined_pathway_dotplots"
# out_dir <- "/Users/jy1008/Downloads/pseudobulk_DE_GSEA_V4_01102024/combined_pathway_dotplots_V3"

# cutoff = 0.7
# outname <- "03282024"
# out_dir <- "/Users/jy1008/Downloads/pseudobulk_DE_GSEA_V4_01102024/combined_pathway_dotplots_V6"

# cutoff = 0.6 test
# outname <- "cutoff0.6_04232024"
out_dir <- "/Users/jy1008/Downloads/pseudobulk_DE_GSEA_V4_01102024/combined_pathway_dotplots_V8_circular"

pathway_db <- "GO-BP"
clinical_feat <- "nctx_temporal_fix"
goslim_cutoff <- 0.6
outname <- sprintf("cutoff%s_04232024", goslim_cutoff)
res <- run_gsea_visualization(pathway_db, clinical_feat, gsea_dir, out_dir, outname, desired_num_k_range=c(20, 50, 100), goslim_cutoff=goslim_cutoff)
combined_gt <- res[[1]]
cluster.row <- res[[2]]
pathway_display_order_all <- res[[3]]
draw_circle_plot2(combined_gt, cluster.row, pathway_display_order_all, label_column = "best_cutree_labels_NA")



pathway_db <- "GO-MF"
clinical_feat <- "nctx_temporal_fix"
goslim_cutoff <- 0.7
outname <- sprintf("cutoff%s_04232024", goslim_cutoff)
res <- run_gsea_visualization(pathway_db, clinical_feat, gsea_dir, out_dir, outname, desired_num_k_range=c(25, 75, 87), goslim_cutoff=goslim_cutoff)
combined_gt <- res[[1]]
cluster.row <- res[[2]]
pathway_display_order_all <- res[[3]]
draw_circle_plot2(combined_gt, cluster.row, pathway_display_order_all, label_column = "best_cutree_labels_NA")
draw_circle_plot2(combined_gt, cluster.row, pathway_display_order_all, label_column = "best_cutree_labels_75")


pathway_db <- "GO-CC"
clinical_feat <- "nctx_temporal_fix"
goslim_cutoff <- 0.7
outname <- sprintf("cutoff%s_04232024", goslim_cutoff)
res <- run_gsea_visualization(pathway_db, clinical_feat, gsea_dir, out_dir, outname, desired_num_k_range=c(23, 53, 74), goslim_cutoff=goslim_cutoff)
combined_gt <- res[[1]]
cluster.row <- res[[2]]
pathway_display_order_all <- res[[3]]
draw_circle_plot2(combined_gt, cluster.row, pathway_display_order_all, label_column = "best_cutree_labels_NA")
draw_circle_plot2(combined_gt, cluster.row, pathway_display_order_all, label_column = "best_cutree_labels_53")





pathway_db <- "GO-BP"
clinical_feat <- "diagnosis_cont_ilbd"
goslim_cutoff <- 0.6
outname <- sprintf("cutoff%s_04232024", goslim_cutoff)
res <- run_gsea_visualization(pathway_db, clinical_feat, gsea_dir, out_dir, outname, desired_num_k_range=c(25, 44), goslim_cutoff=goslim_cutoff)
combined_gt <- res[[1]]
cluster.row <- res[[2]]
pathway_display_order_all <- res[[3]]
draw_circle_plot2(combined_gt, cluster.row, pathway_display_order_all, label_column = "best_cutree_labels_25")
draw_circle_plot2(combined_gt, cluster.row, pathway_display_order_all, label_column = "best_cutree_labels_NA")

pathway_db <- "GO-BP"
clinical_feat <- "diagnosis_cont_case"
goslim_cutoff <- 0.6
outname <- sprintf("cutoff%s_04232024", goslim_cutoff)
res <- run_gsea_visualization(pathway_db, clinical_feat, gsea_dir, out_dir, outname, desired_num_k_range=c(24, 36, 50), goslim_cutoff=goslim_cutoff)
combined_gt <- res[[1]]
cluster.row <- res[[2]]
pathway_display_order_all <- res[[3]]
draw_circle_plot2(combined_gt, cluster.row, pathway_display_order_all, label_column = "best_cutree_labels_24")
draw_circle_plot2(combined_gt, cluster.row, pathway_display_order_all, label_column = "best_cutree_labels_NA")

pathway_db <- "GO-BP"
clinical_feat <- "layer_vs"
goslim_cutoff <- 0.6
outname <- sprintf("cutoff%s_04232024", goslim_cutoff)
res <- run_gsea_visualization(pathway_db, clinical_feat, gsea_dir, out_dir, outname, desired_num_k_range=c(24, 36, 50), goslim_cutoff=goslim_cutoff)

pathway_db <- "GO-BP"
clinical_feat <- "last_mmse_test_score"
goslim_cutoff <- 0.6
outname <- sprintf("cutoff%s_04232024", goslim_cutoff)
res <- run_gsea_visualization(pathway_db, clinical_feat, gsea_dir, out_dir, outname, desired_num_k_range=c(39, 73), goslim_cutoff=goslim_cutoff)
combined_gt <- res[[1]]
cluster.row <- res[[2]]
pathway_display_order_all <- res[[3]]
draw_circle_plot2(combined_gt, cluster.row, pathway_display_order_all, label_column = "best_cutree_labels_39")
draw_circle_plot2(combined_gt, cluster.row, pathway_display_order_all, label_column = "best_cutree_labels_73")

pathway_db <- "GO-BP"
clinical_feat <- "motor_updrs_score"
goslim_cutoff <- 0.6
outname <- sprintf("cutoff%s_04232024", goslim_cutoff)
res <- run_gsea_visualization(pathway_db, clinical_feat, gsea_dir, out_dir, outname, desired_num_k_range=c(24, 36, 50), goslim_cutoff=goslim_cutoff)
