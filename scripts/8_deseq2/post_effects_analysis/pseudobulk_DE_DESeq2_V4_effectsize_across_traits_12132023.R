library(stringr)
library(dplyr)
library(tidyr)
library(testit)
library(ggplot2)


results_list <- list()
de_dir <- "/mnt/data0/projects/ASAP/ST/clustering_pipeline/pseudobulk_DE_DESeq2_V4_12052023"
layers <- c(paste("Layer", 1:6), "WM")
for(layer in layers) {
    results_list[[layer]] <- list()
}
# quantitative features
for(trait in c("nctx_temporal_fix", "brain_stem_sn", "last_mmse_test_score", "motor_updrs_score")) {
    for(layer in layers) {
        layer_nospace <- gsub(' ', '', layer)
        filename <- file.path(de_dir, trait, sprintf("pseudobulk_clinical_deseq2_%s_%s_12052023.csv", trait, layer_nospace))
        if(file.exists(filename)) {
            results <- read.csv(filename)
            colnames(results)[colnames(results) == "X"] <- "gene"
            results_list[[layer]][[trait]] <- results
        } else {
            print(paste("File doesn't exist:", filename))
        }
    }
}
# ordinal cont/ILBD/cases
trait <- "diagnosis_ordinal"
for(layer in layers) {
    layer_nospace <- gsub(' ', '', layer)
    filename <- file.path(de_dir, "diagnosis_ordinal_12122023", sprintf("pseudobulk_clinical_deseq2_diagnosis_%s_12122023.csv", layer_nospace))
    if(file.exists(filename)) {
        results <- read.csv(filename)
        colnames(results)[colnames(results) == "X"] <- "gene"
        results_list[[layer]][[trait]] <- results
    } else {
        print(paste("File doesn't exist:", filename))
    }
}
# cont/ILBD and cont/case
for(trait in c("diagnosis_cont_ilbd", "diagnosis_cont_case")) {
    for(layer in layers) {
        layer_nospace <- gsub(' ', '', layer)
        filename <- file.path(de_dir, "diagnosis_contilbd_contcase", sprintf("pseudobulk_clinical_deseq2_%s_%s_12132023.csv", trait, layer_nospace))
        if(file.exists(filename)) {
            results <- read.csv(filename)
            colnames(results)[colnames(results) == "X"] <- "gene"
            results_list[[layer]][[trait]] <- results
        } else {
            print(paste("File doesn't exist:", filename))
        }
    }
}
# single-cell DE

plot_effects_2d <- function(results_list, traitA, traitB) {
    layers <- c(paste("Layer", 1:6), "WM")
    pdf(sprintf("pseudobulk_clinical_deseq2_%s_vs_%s_12212023.pdf", traitA, traitB))
    for(layer in layers) {
        if( !all(c(traitA, traitB) %in% names(results_list[[layer]])) ) {
            print(paste("Not present in results list:", traitA, traitB, layer))
            next
        }
        resultsA <- results_list[[layer]][[traitA]]
        resultsB <- results_list[[layer]][[traitB]]
        res_comb <- merge(resultsA, resultsB, by="gene", all=TRUE, suffixes = c(paste0('.', traitA), paste0('.', traitB)))
        res_comb[ is.na( res_comb[, paste0("log2FoldChange.", traitA)] ), paste0("log2FoldChange.", traitA) ] <- 0
        res_comb[ is.na( res_comb[, paste0("log2FoldChange.", traitB)] ), paste0("log2FoldChange.", traitB) ] <- 0
        res_comb[ is.na( res_comb[, paste0("pvalue.", traitA)] ), paste0("pvalue.", traitA) ] <- 1
        res_comb[ is.na( res_comb[, paste0("pvalue.", traitB)] ), paste0("pvalue.", traitB) ] <- 1
        res_comb[ is.na( res_comb[, paste0("padj.", traitA)] ), paste0("padj.", traitA) ] <- 1
        res_comb[ is.na( res_comb[, paste0("padj.", traitB)] ), paste0("padj.", traitB) ] <- 1
        
        res_comb[, "padj_sig"] <- "Neither"
        res_comb[res_comb[, paste0("padj.", traitA)] < 0.05 & res_comb[, paste0("padj.", traitB)] >= 0.05, "padj_sig"] <- traitA
        res_comb[res_comb[, paste0("padj.", traitA)] >= 0.05 & res_comb[, paste0("padj.", traitB)] < 0.05, "padj_sig"] <- traitB
        res_comb[res_comb[, paste0("padj.", traitA)] < 0.05 & res_comb[, paste0("padj.", traitB)] < 0.05, "padj_sig"] <- "Both"
        
        # label only if significant in one or more
        res_comb$gene_label <- res_comb$gene
        res_comb[!(res_comb$padj_sig == "Both"), "gene_label"] <- ""
        
        sp_pt <- ggplot(res_comb, aes_string(x=paste0("log2FoldChange.", traitA), y=paste0("log2FoldChange.", traitB))) + geom_point(data=res_comb, aes(color=padj_sig)) + geom_text(aes(label = gene_label), check_overlap=FALSE) + ggtitle(layer) + theme_bw()
        print(sp_pt)
    }
    dev.off()
}
plot_effects_2d(results_list, "diagnosis_cont_ilbd", "diagnosis_cont_case")
plot_effects_2d(results_list, "nctx_temporal_fix", "brain_stem_sn")
plot_effects_2d(results_list, "nctx_temporal_fix", "last_mmse_test_score")
plot_effects_2d(results_list, "nctx_temporal_fix", "motor_updrs_score")
plot_effects_2d(results_list, "brain_stem_sn", "last_mmse_test_score")
plot_effects_2d(results_list, "brain_stem_sn", "motor_updrs_score")
plot_effects_2d(results_list, "last_mmse_test_score", "motor_updrs_score")

plot_effects_2d_ordinal <- function(results_list) {
    traitA <- "diagnosis_cont_ilbd"
    traitB <- "diagnosis_cont_case"
    trait_ord <- "diagnosis_ordinal" # use this trait to filter significant genes
    
    layers <- c(paste("Layer", 1:6), "WM")
    pdf(sprintf("pseudobulk_clinical_deseq2_%s_vs_%s_ordsigonly_12212023.pdf", traitA, traitB))
    for(layer in layers) {
        if( !all(c(traitA, traitB) %in% names(results_list[[layer]])) ) {
            print(paste("Not present in results list:", traitA, traitB, layer))
            next
        }
        resultsA <- results_list[[layer]][[traitA]]
        resultsB <- results_list[[layer]][[traitB]]
        results_ord <- results_list[[layer]][[trait_ord]]
        res_comb <- merge(resultsA, resultsB, by="gene", all=TRUE, suffixes = c(paste0('.', traitA), paste0('.', traitB)))
        res_comb[ is.na( res_comb[, paste0("log2FoldChange.", traitA)] ), paste0("log2FoldChange.", traitA) ] <- 0
        res_comb[ is.na( res_comb[, paste0("log2FoldChange.", traitB)] ), paste0("log2FoldChange.", traitB) ] <- 0
        res_comb[ is.na( res_comb[, paste0("pvalue.", traitA)] ), paste0("pvalue.", traitA) ] <- 1
        res_comb[ is.na( res_comb[, paste0("pvalue.", traitB)] ), paste0("pvalue.", traitB) ] <- 1
        res_comb[ is.na( res_comb[, paste0("padj.", traitA)] ), paste0("padj.", traitA) ] <- 1
        res_comb[ is.na( res_comb[, paste0("padj.", traitB)] ), paste0("padj.", traitB) ] <- 1
        
        res_comb[, "padj_sig"] <- "Neither"
        res_comb[res_comb[, paste0("padj.", traitA)] < 0.05 & res_comb[, paste0("padj.", traitB)] >= 0.05, "padj_sig"] <- traitA
        res_comb[res_comb[, paste0("padj.", traitA)] >= 0.05 & res_comb[, paste0("padj.", traitB)] < 0.05, "padj_sig"] <- traitB
        res_comb[res_comb[, paste0("padj.", traitA)] < 0.05 & res_comb[, paste0("padj.", traitB)] < 0.05, "padj_sig"] <- "Both"
        
        results_ord[is.na(results_ord$padj), "padj"] <- 1
        sig_gene_list <- results_ord[results_ord[, "padj"] < 0.05, "gene"]
        print(layer)
        print(sig_gene_list)
        res_comb <- res_comb[res_comb$gene %in% sig_gene_list, ]
        print(dim(res_comb))
        
        sp_pt <- ggplot(res_comb, aes_string(x=paste0("log2FoldChange.", traitA), y=paste0("log2FoldChange.", traitB))) + geom_point(data=res_comb, aes(color=padj_sig)) + theme_bw()
        print(sp_pt)
    }
    dev.off()
}
plot_effects_2d_ordinal(results_list)
