# see script "pseudobulk_DE_DESeq2_post_analysis_gseaclusterprofiler_subcats_02142023"

# GSEA with cluster profiler
# https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html
# https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/

library(clusterProfiler)
library(enrichplot)
library(ggplot2)

organism = "org.Hs.eg.db"
# BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

# DOSE: an R/Bioconductor package for Disease Ontology Semantic and Enrichment analysis
# Yu et al., 2015
require(DOSE)
library(ggnewscale)

source("/Users/jy1008/Downloads/pseudobulk_DE_DESeq2_GSEA_12292023_utils.R")

reformat_gsea_df <- function(gse_obj) {
    df <- gse_obj@result
    setType <- paste0("GO:", gse_obj@setType)
    new_df <- data.frame(
                gene_set = setType,
                pathway_id = df$ID,
                pathway = df$Description,
                pval = df$pvalue,
                padj = df$p.adjust,
                log2err = NA,
                ES = df$enrichmentScore,
                NES = df$NES,
                size = df$setSize,
                leadingEdge = NA
                    )
    for(i in 1:nrow(df)) {
        new_df[i, "leadingEdge"] <- paste(strsplit(df$core_enrichment[i], '/')[[1]], collapse=', ')
    }
    return(new_df)
}

run_clusterprofiler <- function(de_files_list, clinical_feat, ont, outname, slim=TRUE) {
    # input_files <- list("comparison name" = "file_name")
    
    gse_results_all_layers <- c()

    for(st_layer in names(de_files_list)) {
        print(sprintf("%s * %s", clinical_feat, st_layer))

        results <- read.csv(de_files_list[[st_layer]])
        results$gene <- results$X
        results$X <- NULL
        results$comparison <- sprintf("%s * %s", clinical_feat, gsub(' ', '_', st_layer))
        
        # switch to sorting by pvalue, not padj, as many padjs are NA
        # results[is.na(results$padj), "padj"] <- 1
        results$nlog10_pval <- -log10(results$pval)
        max_val <- max(results[!is.infinite(results$nlog10_pval), ]$nlog10_pval)
        inf_subst <- round(max_val + 50, -2)
        if ( any(is.infinite(results$nlog10_pval)) ) {
            results[is.infinite(results$nlog10_pval), ]$nlog10_pval <- inf_subst
        }
        
        # switch back to using just the sign of the fold change
        ### A novel significance score for gene selection and ranking
        ### Xiao et al., 2014
        results$pi_score <- results$nlog10_pval * sign(results$log2FoldChange)
        
        # https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/
        # we want the log2 fold change
        original_gene_list <- results$pi_score
        # name the vector
        names(original_gene_list) <- results$gene
        # omit any NA values
        gene_list <- na.omit(original_gene_list)
        # gene_list <- gene_list[gene_list != 0]
        # sort the list in decreasing order (required for clusterProfiler)
        gene_list = sort(gene_list, decreasing = TRUE)

        
        # ont ="ALL",
        # minGSSize = 3, maxGSSize = 800
        gse <- gseGO(geneList=gene_list,
                     ont = ont,
                     keyType = "SYMBOL",
                     minGSSize = 10,
                     maxGSSize = 500,
                     # pvalueCutoff = 0.05,
                     pvalueCutoff = 1.0,
                     verbose = TRUE,
                     OrgDb = get(organism),
                     pAdjustMethod = "BH",
                     eps = 0)
        
        # if no sig pathways, skip to next test
        if(dim(gse@result)[1] == 0) {
            print("No significant pathways")
            next
        }
        
        if(slim) {
            # NOTE: this doesn't work with gseGO(ont = "ALL", ...)
            gse_slim <- simplify(gse, cutoff=0.7, by="p.adjust", select_fun=min)
            gse_slim_new_df <- reformat_gsea_df(gse_slim)
            gse_slim_new_df$comparison <- sprintf("%s * %s", clinical_feat, gsub(' ', '_', st_layer))
            gse_results_all_layers <- rbind(gse_results_all_layers, gse_slim_new_df)
            gse <- gse_slim # for plotting
        } else {
            gse_new_df <- reformat_gsea_df(gse)
            gse_new_df$comparison <- sprintf("%s * %s", clinical_feat, gsub(' ', '_', st_layer))
            gse_results_all_layers <- rbind(gse_results_all_layers, gse_new_df)
        }

        # showCategory is required
        pdf(sprintf("clusterProfiler_gseGO_dotplot_%s_%s_ont%s_slim%s_%s.pdf", clinical_feat, gsub(' ', '_', st_layer), ont, slim, outname), width=10, height=60)
        pt <- dotplot(gse, showCategory=(dim(gse@result)[1]+1), split=".sign", color="pvalue") + facet_grid(.~.sign) + theme(axis.text.y = element_text(size = 5))
        print(pt)
        dev.off()
        
        
        x2 <- pairwise_termsim(gse)
        emap_sizes <- c(5, 20, 50, 100)
        for(emap_size in emap_sizes) {
            pdf(sprintf("clusterProfiler_gseGO_emap_%s_%s_size%s_ont%s_slim%s_%s.pdf", clinical_feat, gsub(' ', '_', st_layer), emap_size, ont, slim, outname), width=50, height=50)
            pt <- emapplot(x2, showCategory=emap_size, color="pvalue", cex_label_category=1)
            print(pt)
            dev.off()
        }
        
        cnet_sizes <- c(3, 5, 10, 20)
        for(cnet_size in cnet_sizes) {
            pdf(sprintf("clusterProfiler_gseGO_cnet_%s_%s_size%s_ont%s_slim%s_%s.pdf", clinical_feat, gsub(' ', '_', st_layer), cnet_size, ont, slim, outname), width=50, height=50)
            pt <- cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = cnet_size, cex_label_category=1)
            print(pt)
            dev.off()
        }
    }
    return(gse_results_all_layers)
}



random_seed <- 101
layers <- c("Layer 1", "Layer 2", "Layer 3", "Layer 4", "Layer 5", "Layer 6", "WM")
effects_path <- "/Users/jy1008/Downloads/pseudobulk_DE_DESeq2_V4_12052023"

use_go_slim <- FALSE
outname <- "01162024"
out_go_path <- "/Users/jy1008/Downloads/pseudobulk_DE_GSEA_V4_01102024"
clinical_features <- c("layer_vs",
                       "diagnosis_cont_ilbd",
                       "diagnosis_cont_case",
                       "last_mmse_test_score",
                       "motor_updrs_score",
                       "nctx_temporal_fix",
                       "brain_stem_sn")

# Biological Process, Molecular Function, Cellular Component
ont_list <- c("BP", "MF", "CC")

for(clinical_feat in clinical_features) {
    for(ont in ont_list) {
        
        if(clinical_feat == "layer_vs") {
            de_files_list <- gsea_layer_vs_layer_filelist(effects_path, layers, contsonly=FALSE)
        } else {
            de_files_list <- gsea_clinical_per_layer_filelist(effects_path, clinical_feat, layers)
        }

        if(use_go_slim) {
            gsea_dir <- file.path(out_go_path, sprintf("GO-%s-slim", ont))
        } else {
            gsea_dir <- file.path(out_go_path, sprintf("GO-%s", ont))
        }
        if(!file.exists(gsea_dir)) {
            dir.create(gsea_dir)
        }
        setwd(gsea_dir)
        set.seed(random_seed)
        gsea_results <- run_clusterprofiler(de_files_list=de_files_list,
                                            clinical_feat=clinical_feat,
                                            ont=ont, outname=outname, slim=use_go_slim)
        is_slim <- if (use_go_slim) "-slim" else ""
        write.table(gsea_results, sprintf("de_genes_DESeq2_pseudobulk_%s_GO-%s%s_%s.fGSEA.xls", clinical_feat, ont, is_slim, outname), sep='\t', row.names=FALSE)
    }
}
