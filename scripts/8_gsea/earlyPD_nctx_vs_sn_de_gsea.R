# see ~/Downloads/process_layer_output_top5_heatmap_02012024.R
# see ~/Downloads/process_layer_output_MAplots_02012024.R

library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)
# library(matrixStats)
# library(corrplot)
library(testit)
# library(ggdendro)
# library(DESeq2)
library(ggrepel)


# read each layer DE file in MTG and SN and plot in scatterplot
pdf("earlyPD_nctx_vs_sn_de_genes_10032024.pdf")
st_layers <- c(paste0("Layer",1:6), "WM")
for(st_layer in st_layers) {
    mtg_csv_path <- file.path("/Users/jy1008/Downloads/pseudobulk_DE_DESeq2_V4_12052023", "nctx_temporal_fix")
    sn_csv_path <- file.path("/Users/jy1008/Downloads/pseudobulk_DE_DESeq2_V4_12052023", "brain_stem_sn")

    mtg_de <- read.csv(file.path(mtg_csv_path,
                        sprintf("pseudobulk_clinical_deseq2_%s_%s_12052023.csv",
                        "nctx_temporal_fix", st_layer)), row.names = 1)
    mtg_de[is.na(mtg_de$padj), "padj"] <- 1
    mtg_de$gene <- rownames(mtg_de)
    sn_de <- read.csv(file.path(sn_csv_path,
                        sprintf("pseudobulk_clinical_deseq2_%s_%s_12052023.csv",
                        "brain_stem_sn", st_layer)), row.names = 1)
    sn_de[is.na(sn_de$padj), "padj"] <- 1
    sn_de$gene <- rownames(sn_de)

    de_merged <- merge(mtg_de, sn_de, by="gene", suffixes = c(".mtg",".sn"))
    # keep only genes which are singificant in at least one comparison
    de_merged$is_sig <- NA
    de_merged$is_sig[de_merged$padj.mtg < 0.05] <- "mtg_only"
    de_merged$is_sig[de_merged$padj.sn < 0.05] <- "sn_only"
    de_merged$is_sig[(de_merged$padj.mtg < 0.05) &
                     (de_merged$padj.sn < 0.05)] <- "mtg_and_sn"

    de_filt <- de_merged[!is.na(de_merged$is_sig), ]
    pt <- ggplot(de_filt, aes(x=log2FoldChange.mtg, y=log2FoldChange.sn, label=gene, color=is_sig)) +
            geom_point() +
            geom_text(hjust=0, vjust=0) +
            coord_fixed() +
            geom_abline(slope=1, intercept=0) +
            theme_bw() + labs(title=st_layer)
    print(pt)
    pt <- ggplot(de_filt, aes(x=log2FoldChange.mtg, y=log2FoldChange.sn, label=gene, color=is_sig)) +
            geom_point() +
            # geom_text(hjust=0, vjust=0) +
            coord_fixed() +
            geom_abline(slope=1, intercept=0) +
            theme_bw() + labs(title=st_layer)
    print(pt)
}
dev.off()




# read each layer DE file in MTG and control vs ILBD and plot in scatterplot
pdf("earlyPD_nctx_vs_contilbd_de_genes_10032024.pdf")
st_layers <- c(paste0("Layer",1:6), "WM")
for(st_layer in st_layers) {
    mtg_csv_path <- file.path("/Users/jy1008/Downloads/pseudobulk_DE_DESeq2_V4_12052023", "nctx_temporal_fix")
    contilbd_csv_path <- file.path("/Users/jy1008/Downloads/pseudobulk_DE_DESeq2_V4_12052023", "diagnosis_cont_ilbd")

    mtg_de <- read.csv(file.path(mtg_csv_path,
                        sprintf("pseudobulk_clinical_deseq2_%s_%s_12052023.csv",
                        "nctx_temporal_fix", st_layer)), row.names = 1)
    mtg_de[is.na(mtg_de$padj), "padj"] <- 1
    mtg_de$gene <- rownames(mtg_de)
    contilbd_de <- read.csv(file.path(contilbd_csv_path,
                        sprintf("pseudobulk_clinical_deseq2_%s_%s_12132023.csv",
                        "diagnosis_cont_ilbd", st_layer)), row.names = 1)
    contilbd_de[is.na(sn_de$padj), "padj"] <- 1
    contilbd_de$gene <- rownames(contilbd_de)

    de_merged <- merge(mtg_de, contilbd_de, by="gene", suffixes = c(".mtg",".contilbd"))
    # keep only genes which are singificant in at least one comparison
    de_merged$is_sig <- NA
    de_merged$is_sig[de_merged$padj.mtg < 0.05] <- "mtg_only"
    de_merged$is_sig[de_merged$padj.contilbd < 0.05] <- "ILBD/cont only"
    de_merged$is_sig[(de_merged$padj.mtg < 0.05) &
                     (de_merged$padj.contilbd < 0.05)] <- "mtg_and_ILBD/cont"

    de_filt <- de_merged[!is.na(de_merged$is_sig), ]
    pt <- ggplot(de_filt, aes(x=log2FoldChange.mtg, y=log2FoldChange.contilbd, label=gene, color=is_sig)) +
            geom_point() +
            geom_text(hjust=0, vjust=0) +
            coord_fixed() +
            geom_abline(slope=1, intercept=0) +
            theme_bw() + labs(title=st_layer)
    print(pt)
    pt <- ggplot(de_filt, aes(x=log2FoldChange.mtg, y=log2FoldChange.contilbd, label=gene, color=is_sig)) +
            geom_point() +
            # geom_text(hjust=0, vjust=0) +
            coord_fixed() +
            geom_abline(slope=1, intercept=0) +
            theme_bw() + labs(title=st_layer)
    print(pt)
}
dev.off()



# read each layer DE file in MTG and controlILBDvscase and plot in scatterplot
pdf("earlyPD_nctx_vs_contilbdvscase_de_genes_10032024.pdf")
st_layers <- c(paste0("Layer",1:6), "WM")
for(st_layer in st_layers) {
    mtg_csv_path <- file.path("/Users/jy1008/Downloads/pseudobulk_DE_DESeq2_V4_12052023", "nctx_temporal_fix")
    contilbd_csv_path <- file.path("/Users/jy1008/Downloads/pseudobulk_DE_DESeq2_V4_12052023", "diagnosis_contilbd_case")

    mtg_de <- read.csv(file.path(mtg_csv_path,
                        sprintf("pseudobulk_clinical_deseq2_%s_%s_12052023.csv",
                        "nctx_temporal_fix", st_layer)), row.names = 1)
    mtg_de[is.na(mtg_de$padj), "padj"] <- 1
    mtg_de$gene <- rownames(mtg_de)
    contilbd_de <- read.csv(file.path(contilbd_csv_path,
                        sprintf("pseudobulk_clinical_deseq2_%s_%s_06032024.csv",
                        "diagnosis_contilbd_case", st_layer)), row.names = 1)
    contilbd_de[is.na(sn_de$padj), "padj"] <- 1
    contilbd_de$gene <- rownames(contilbd_de)

    de_merged <- merge(mtg_de, contilbd_de, by="gene", suffixes = c(".mtg",".contilbd"))
    # keep only genes which are singificant in at least one comparison
    de_merged$is_sig <- NA
    de_merged$is_sig[de_merged$padj.mtg < 0.05] <- "mtg_only"
    de_merged$is_sig[de_merged$padj.contilbd < 0.05] <- "contILBD/case only"
    de_merged$is_sig[(de_merged$padj.mtg < 0.05) &
                     (de_merged$padj.contilbd < 0.05)] <- "mtg_and_contILBD/case"

    de_filt <- de_merged[!is.na(de_merged$is_sig), ]
    pt <- ggplot(de_filt, aes(x=log2FoldChange.mtg, y=log2FoldChange.contilbd, label=gene, color=is_sig)) +
            geom_point() +
            geom_text(hjust=0, vjust=0) +
            coord_fixed() +
            geom_abline(slope=1, intercept=0) +
            theme_bw() + labs(title=st_layer)
    print(pt)
    pt <- ggplot(de_filt, aes(x=log2FoldChange.mtg, y=log2FoldChange.contilbd, label=gene, color=is_sig)) +
            geom_point() +
            # geom_text(hjust=0, vjust=0) +
            coord_fixed() +
            geom_abline(slope=1, intercept=0) +
            theme_bw() + labs(title=st_layer)
    print(pt)
}
dev.off()


# read each layer DE file in controlvsILBD and controlILBDvscase and plot in scatterplot
pdf("earlyPD_contvsILBD_vs_contilbdvscase_de_genes_10032024.pdf")
st_layers <- c(paste0("Layer",1:6), "WM")
for(st_layer in st_layers) {
    mtg_csv_path <- file.path("/Users/jy1008/Downloads/pseudobulk_DE_DESeq2_V4_12052023", "diagnosis_cont_ilbd")
    contilbd_csv_path <- file.path("/Users/jy1008/Downloads/pseudobulk_DE_DESeq2_V4_12052023", "diagnosis_contilbd_case")

    mtg_de <- read.csv(file.path(mtg_csv_path,
                        sprintf("pseudobulk_clinical_deseq2_%s_%s_12132023.csv",
                        "diagnosis_cont_ilbd", st_layer)), row.names = 1)
    mtg_de[is.na(mtg_de$padj), "padj"] <- 1
    mtg_de$gene <- rownames(mtg_de)
    contilbd_de <- read.csv(file.path(contilbd_csv_path,
                        sprintf("pseudobulk_clinical_deseq2_%s_%s_06032024.csv",
                        "diagnosis_contilbd_case", st_layer)), row.names = 1)
    contilbd_de[is.na(sn_de$padj), "padj"] <- 1
    contilbd_de$gene <- rownames(contilbd_de)

    de_merged <- merge(mtg_de, contilbd_de, by="gene", suffixes = c(".mtg",".contilbd"))
    # keep only genes which are singificant in at least one comparison
    de_merged$is_sig <- NA
    de_merged$is_sig[de_merged$padj.mtg < 0.05] <- "ILBD/cont only"
    de_merged$is_sig[de_merged$padj.contilbd < 0.05] <- "case/contILBD only"
    de_merged$is_sig[(de_merged$padj.mtg < 0.05) &
                     (de_merged$padj.contilbd < 0.05)] <- "mtg_and_case/contILBD"

    de_filt <- de_merged[!is.na(de_merged$is_sig), ]
    pt <- ggplot(de_filt, aes(x=log2FoldChange.mtg, y=log2FoldChange.contilbd, label=gene, color=is_sig)) +
            geom_point() +
            geom_text(hjust=0, vjust=0) +
            coord_fixed() +
            geom_abline(slope=1, intercept=0) +
            theme_bw() + labs(title=st_layer)
    print(pt)
    pt <- ggplot(de_filt, aes(x=log2FoldChange.mtg, y=log2FoldChange.contilbd, label=gene, color=is_sig)) +
            geom_point() +
            # geom_text(hjust=0, vjust=0) +
            coord_fixed() +
            geom_abline(slope=1, intercept=0) +
            theme_bw() + labs(title=st_layer)
    print(pt)
}
dev.off()




# read each layer DE file in cont/ilbd redef and cont/case redef and plot in scatterplot
pdf("earlyPD_contilbd_redef_vs_contcase_redef_de_genes_10102024.pdf")
st_layers <- c(paste0("Layer",1:6), "WM")
for(st_layer in st_layers) {
    contilbd_csv_path <- file.path("/Users/jy1008/Downloads/pseudobulk_DE_DESeq2_V4_12052023", "diagnosis_cont_ilbd_redef")
    contcase_csv_path <- file.path("/Users/jy1008/Downloads/pseudobulk_DE_DESeq2_V4_12052023", "diagnosis_cont_case_redef")

    contilbd_de <- read.csv(file.path(contilbd_csv_path,
                        sprintf("pseudobulk_clinical_deseq2_%s_%s_10102024.csv",
                        "diagnosis_cont_ilbd_redef", st_layer)), row.names = 1)
    contilbd_de[is.na(contilbd_de$padj), "padj"] <- 1
    contilbd_de$gene <- rownames(contilbd_de)
    contcase_de <- read.csv(file.path(contcase_csv_path,
                        sprintf("pseudobulk_clinical_deseq2_%s_%s_10102024.csv",
                        "diagnosis_cont_case_redef", st_layer)), row.names = 1)
    contcase_de[is.na(contcase_de$padj), "padj"] <- 1
    contcase_de$gene <- rownames(contcase_de)

    de_merged <- merge(contilbd_de, contcase_de, by="gene", suffixes = c(".contilbd",".contcase"))
    # keep only genes which are singificant in at least one comparison
    de_merged$is_sig <- NA
    de_merged$is_sig[de_merged$padj.contilbd < 0.05] <- "contilbd_only"
    de_merged$is_sig[de_merged$padj.contcase < 0.05] <- "contcase_only"
    de_merged$is_sig[(de_merged$padj.contilbd < 0.05) &
                     (de_merged$padj.contcase < 0.05)] <- "contilbd_and_contcase"

    de_filt <- de_merged[!is.na(de_merged$is_sig), ]
    pt <- ggplot(de_filt, aes(x=log2FoldChange.contilbd, y=log2FoldChange.contcase, label=gene, color=is_sig)) +
            geom_point() +
            geom_text(hjust=0, vjust=0) +
            coord_fixed() +
            geom_abline(slope=1, intercept=0) +
            theme_bw() + labs(title=st_layer)
    print(pt)
    pt <- ggplot(de_filt, aes(x=log2FoldChange.contilbd, y=log2FoldChange.contcase, label=gene, color=is_sig)) +
            geom_point() +
            # geom_text(hjust=0, vjust=0) +
            coord_fixed() +
            geom_abline(slope=1, intercept=0) +
            theme_bw() + labs(title=st_layer)
    print(pt)
}
dev.off()


# read each layer DE file in cont/ilbd redef and cont/case redef and plot in scatterplot
pdf("earlyPD_contilbdvscase_vs_contcaseredef_de_genes_10102024.pdf")
st_layers <- c(paste0("Layer",1:6), "WM")
for(st_layer in st_layers) {
    contilbd_csv_path <- file.path("/Users/jy1008/Downloads/pseudobulk_DE_DESeq2_V4_12052023", "diagnosis_contilbd_case")
    contcase_csv_path <- file.path("/Users/jy1008/Downloads/pseudobulk_DE_DESeq2_V4_12052023", "diagnosis_cont_case_redef")

    contilbd_de <- read.csv(file.path(contilbd_csv_path,
                        sprintf("pseudobulk_clinical_deseq2_%s_%s_06032024.csv",
                        "diagnosis_contilbd_case", st_layer)), row.names = 1)
    contilbd_de[is.na(contilbd_de$padj), "padj"] <- 1
    contilbd_de$gene <- rownames(contilbd_de)
    contcase_de <- read.csv(file.path(contcase_csv_path,
                        sprintf("pseudobulk_clinical_deseq2_%s_%s_10102024.csv",
                        "diagnosis_cont_case_redef", st_layer)), row.names = 1)
    contcase_de[is.na(contcase_de$padj), "padj"] <- 1
    contcase_de$gene <- rownames(contcase_de)

    de_merged <- merge(contilbd_de, contcase_de, by="gene", suffixes = c(".contilbd",".contcase"))
    # keep only genes which are singificant in at least one comparison
    de_merged$is_sig <- NA
    de_merged$is_sig[de_merged$padj.contilbd < 0.05] <- "ctilbd_vs_case_only"
    de_merged$is_sig[de_merged$padj.contcase < 0.05] <- "contcase_redef_only"
    de_merged$is_sig[(de_merged$padj.contilbd < 0.05) &
                     (de_merged$padj.contcase < 0.05)] <- "both"

    de_filt <- de_merged[!is.na(de_merged$is_sig), ]
    pt <- ggplot(de_filt, aes(x=log2FoldChange.contilbd, y=log2FoldChange.contcase, label=gene, color=is_sig)) +
            geom_point() +
            geom_text(hjust=0, vjust=0) +
            coord_fixed() +
            geom_abline(slope=1, intercept=0) +
            theme_bw() + labs(title=st_layer)
    print(pt)
    pt <- ggplot(de_filt, aes(x=log2FoldChange.contilbd, y=log2FoldChange.contcase, label=gene, color=is_sig)) +
            geom_point() +
            # geom_text(hjust=0, vjust=0) +
            coord_fixed() +
            geom_abline(slope=1, intercept=0) +
            theme_bw() + labs(title=st_layer)
    print(pt)
}
dev.off()
























# /Users/jy1008/Downloads/pseudobulk_DE_GSEA_V4_01102024/GO-BP/de_genes_DESeq2_pseudobulk_nctx_temporal_fix_GO-BP_01162024.fGSEA.xls


# read each layer DE file in MTG and control vs ILBD and plot in scatterplot
pdf("earlyPD_nctx_vs_contilbd_gsea_10072024.pdf")
st_layers_gsea <- c(paste0("Layer_",1:6), "WM")

csv_path <- "/Users/jy1008/Downloads/pseudobulk_DE_GSEA_V4_01102024/GO-BP/"
# mtg_csv_path <- file.path("/Users/jy1008/Downloads/pseudobulk_DE_DESeq2_V4_12052023", "nctx_temporal_fix")
# contilbd_csv_path <- file.path("/Users/jy1008/Downloads/pseudobulk_DE_DESeq2_V4_12052023", "diagnosis_cont_ilbd")

mtg_gsea <- read.table(file.path(csv_path,
                    sprintf("de_genes_DESeq2_pseudobulk_%s_GO-BP_01162024.fGSEA.xls", "nctx_temporal_fix")), header=TRUE)
contilbd_gsea <- read.table(file.path(csv_path,
                    sprintf("de_genes_DESeq2_pseudobulk_%s_GO-BP_01162024.fGSEA.xls", "diagnosis_cont_ilbd")), header=TRUE)

for(st_layer in st_layers_gsea) {
    mtg_sub <- mtg_gsea[grepl(st_layer, mtg_gsea$comparison, fixed=TRUE), ]
    contilbd_sub <- contilbd_gsea[grepl(st_layer, contilbd_gsea$comparison, fixed=TRUE), ]

    de_merged <- merge(mtg_sub, contilbd_sub, by="pathway", suffixes = c(".mtg",".contilbd"))
    # keep only genes which are singificant in at least one comparison
    de_merged$is_sig <- NA
    de_merged$is_sig[de_merged$padj.mtg < 0.05] <- "mtg_only"
    de_merged$is_sig[de_merged$padj.contilbd < 0.05] <- "ILBD/cont only"
    de_merged$is_sig[(de_merged$padj.mtg < 0.05) &
                     (de_merged$padj.contilbd < 0.05)] <- "mtg_and_ILBD/cont"

    de_filt <- de_merged[!is.na(de_merged$is_sig), ]
    de_filt$pathway_label <- de_filt$pathway
    de_filt[de_filt$is_sig != "mtg_and_ILBD/cont", "pathway_label"] <- ""
    pt <- ggplot(de_filt, aes(x=NES.mtg, y=NES.contilbd, label=pathway_label, color=is_sig)) +
            geom_point() +
            # geom_text(hjust=0, vjust=0) +
            geom_text_repel(max.overlaps=Inf) +
            coord_fixed() +
            geom_abline(slope=1, intercept=0) +
            theme_bw() + labs(title=st_layer)
    print(pt)
    pt <- ggplot(de_filt, aes(x=NES.mtg, y=NES.contilbd, label=pathway, color=is_sig)) +
            geom_point() +
            # geom_text(hjust=0, vjust=0) +
            coord_fixed() +
            geom_abline(slope=1, intercept=0) +
            theme_bw() + labs(title=st_layer)
    print(pt)
}
dev.off()






pdf("earlyPD_nctx_vs_contilbdvscase_gsea_10072024.pdf")
st_layers_gsea <- c(paste0("Layer_",1:6), "WM")

csv_path <- "/Users/jy1008/Downloads/pseudobulk_DE_GSEA_V4_01102024/GO-BP/"
# mtg_csv_path <- file.path("/Users/jy1008/Downloads/pseudobulk_DE_DESeq2_V4_12052023", "nctx_temporal_fix")
# contilbd_csv_path <- file.path("/Users/jy1008/Downloads/pseudobulk_DE_DESeq2_V4_12052023", "diagnosis_cont_ilbd")

mtg_gsea <- read.table(file.path(csv_path,
                    sprintf("de_genes_DESeq2_pseudobulk_%s_GO-BP_01162024.fGSEA.xls", "nctx_temporal_fix")), header=TRUE)
contilbd_gsea <- read.table(file.path(csv_path,
                    sprintf("de_genes_DESeq2_pseudobulk_%s_GO-BP_06032024.fGSEA.xls", "diagnosis_contilbd_case")), header=TRUE)

for(st_layer in st_layers_gsea) {
    mtg_sub <- mtg_gsea[grepl(st_layer, mtg_gsea$comparison, fixed=TRUE), ]
    contilbd_sub <- contilbd_gsea[grepl(st_layer, contilbd_gsea$comparison, fixed=TRUE), ]

    de_merged <- merge(mtg_sub, contilbd_sub, by="pathway", suffixes = c(".mtg",".contilbd"))
    # keep only genes which are singificant in at least one comparison
    de_merged$is_sig <- NA
    de_merged$is_sig[de_merged$padj.mtg < 0.05] <- "mtg_only"
    de_merged$is_sig[de_merged$padj.contilbd < 0.05] <- "case/contILBD only"
    de_merged$is_sig[(de_merged$padj.mtg < 0.05) &
                     (de_merged$padj.contilbd < 0.05)] <- "mtg_and_case/contILBD"

    de_filt <- de_merged[!is.na(de_merged$is_sig), ]
    pt <- ggplot(de_filt, aes(x=NES.mtg, y=NES.contilbd, label=pathway, color=is_sig)) +
            geom_point() +
            geom_text(hjust=0, vjust=0) +
            coord_fixed() +
            geom_abline(slope=1, intercept=0) +
            theme_bw() + labs(title=st_layer)
    print(pt)
    pt <- ggplot(de_filt, aes(x=NES.mtg, y=NES.contilbd, label=pathway, color=is_sig)) +
            geom_point() +
            # geom_text(hjust=0, vjust=0) +
            coord_fixed() +
            geom_abline(slope=1, intercept=0) +
            theme_bw() + labs(title=st_layer)
    print(pt)
}
dev.off()





pdf("earlyPD_contilbd_redef_vs_contcase_redef_gsea_10072024.pdf")
st_layers_gsea <- c(paste0("Layer_",1:6), "WM")

csv_path <- "/Users/jy1008/Downloads/pseudobulk_DE_GSEA_V4_01102024/GO-BP/"
# mtg_csv_path <- file.path("/Users/jy1008/Downloads/pseudobulk_DE_DESeq2_V4_12052023", "nctx_temporal_fix")
# contilbd_csv_path <- file.path("/Users/jy1008/Downloads/pseudobulk_DE_DESeq2_V4_12052023", "diagnosis_cont_ilbd")

mtg_gsea <- read.table(file.path(csv_path,
                    sprintf("de_genes_DESeq2_pseudobulk_%s_GO-BP_10102024.fGSEA.xls", "diagnosis_cont_ilbd_redef")), header=TRUE)
contilbd_gsea <- read.table(file.path(csv_path,
                    sprintf("de_genes_DESeq2_pseudobulk_%s_GO-BP_10102024.fGSEA.xls", "diagnosis_cont_case_redef")), header=TRUE)

for(st_layer in st_layers_gsea) {
    mtg_sub <- mtg_gsea[grepl(st_layer, mtg_gsea$comparison, fixed=TRUE), ]
    contilbd_sub <- contilbd_gsea[grepl(st_layer, contilbd_gsea$comparison, fixed=TRUE), ]

    de_merged <- merge(mtg_sub, contilbd_sub, by="pathway", suffixes = c(".mtg",".contilbd"))
    # keep only genes which are singificant in at least one comparison
    de_merged$is_sig <- NA
    de_merged$is_sig[de_merged$padj.mtg < 0.05] <- "ilbd/cont"
    de_merged$is_sig[de_merged$padj.contilbd < 0.05] <- "case/cont"
    de_merged$is_sig[(de_merged$padj.mtg < 0.05) &
                     (de_merged$padj.contilbd < 0.05)] <- "both"

    de_filt <- de_merged[!is.na(de_merged$is_sig), ]
    pt <- ggplot(de_filt, aes(x=NES.mtg, y=NES.contilbd, label=pathway, color=is_sig)) +
            geom_point() +
            geom_text(hjust=0, vjust=0) +
            coord_fixed() +
            geom_abline(slope=1, intercept=0) +
            theme_bw() + labs(title=st_layer)
    print(pt)
    pt <- ggplot(de_filt, aes(x=NES.mtg, y=NES.contilbd, label=pathway, color=is_sig)) +
            geom_point() +
            # geom_text(hjust=0, vjust=0) +
            coord_fixed() +
            geom_abline(slope=1, intercept=0) +
            theme_bw() + labs(title=st_layer)
    print(pt)
    de_merged[is.na(de_merged$is_sig), "is_sig"] <- "neither"
    pt <- ggplot(de_merged, aes(x=NES.mtg, y=NES.contilbd, label=pathway, color=is_sig)) +
            geom_point() +
            # geom_text(hjust=0, vjust=0) +
            coord_fixed() +
            geom_abline(slope=1, intercept=0) +
            theme_bw() + labs(title=st_layer)
    print(pt)
}
dev.off()






temp_file <- "/Users/jy1008/Downloads/pseudobulk_DE_GSEA_V4_12292023/kegg/de_genes_DESeq2_pseudobulk_nctx_temporal_fix_c2.cp.kegg.v7.2.symbols.gmt12292023.fGSEA.xls"
temp_gsea <- read.table(temp_file, header=TRUE, sep='\t')





# dotplot with four conditions: ctrl/ilbd redef, ctrl/case redef, (ctrl+ilbd)/case, nctx_temporal_fix
library(ggtree)
library(cluster)
source("/Users/jy1008/Downloads/pseudobulk_DE_DESeq2_GSEA_12292023_utils.R")

comb_path <- "/Users/jy1008/Downloads/pseudobulk_DE_GSEA_V4_01102024/GO-BP"
comb_files <- list(
"nctx_temporal_fix" = "de_genes_DESeq2_pseudobulk_nctx_temporal_fix_GO-BP_01162024.fGSEA.xls",
"diag_contilbd_case" = "de_genes_DESeq2_pseudobulk_diagnosis_contilbd_case_GO-BP_06032024.fGSEA.xls",
"diag_cont_case_redef" = "de_genes_DESeq2_pseudobulk_diagnosis_cont_case_redef_GO-BP_10152024.fGSEA.xls",
"diag_cont_ilbd_redef" = "de_genes_DESeq2_pseudobulk_diagnosis_cont_ilbd_redef_GO-BP_10152024.fGSEA.xls"
)
combined_gt <- c()
for(clin_trait in names(comb_files)) {
    df_path <- file.path(comb_path, comb_files[[clin_trait]])
    df <- read.table(df_path, header=TRUE, sep='\t')
    # df <- df %>% tidyr::separate(col = comparison, into = c("trait", "layer"), sep = " \\* ")
    combined_gt <- rbind(combined_gt, df)
}
pathway_db <- "GO-BP"
nes_filter <- 0
goslim_cutoff <- 0.7
combined_gt <- prep_gsea_df(combined_gt, pathway_db, nes_cutoff = nes_filter, goslim_cutoff=goslim_cutoff)

# grab subset of pathways which are significant in early PD
early_sig_pathways <- unique(combined_gt[combined_gt$padj < 0.05 &
                                         grepl("diagnosis_cont_ilbd_redef \\* ", combined_gt$comparison),
                                        "pathway"])
combined_gt_orig <- combined_gt
combined_gt <- combined_gt[combined_gt$pathway %in% early_sig_pathways, ]

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
# layer_order <- sort(unique(combined_gt$comparison))
st_layers <- c(paste0("Layer_", 1:6), "WM")
layer_order <- c(paste0("diagnosis_cont_ilbd_redef * ", st_layers),
                 paste0("diagnosis_cont_case_redef * ", st_layers),
                 paste0("diagnosis_contilbd_case * ", st_layers),
                 paste0("nctx_temporal_fix * ", st_layers))

# use silhouette score to combine very similar labels
# desired_num_k_range <- c(50, 100, 150, NA)
desired_num_k_range <- c(25, 50, 65, NA)
for(desired_num_k in desired_num_k_range) {
    temp_res <- get_best_cutree(cluster.row, tot_dist_temp, desired_num_k = desired_num_k)
    sil_widths <- temp_res[[1]]
    best_cutree_labels <- temp_res[[2]]
    combined_gt[, paste0("best_cutree_labels_", desired_num_k)] <- best_cutree_labels[combined_gt$pathway]
}
print(sil_widths)

clinical_feat <- "earlyPD_redef"
outname <- "earlyonly_10152024"
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
    width <- 15
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
    # layer_order <- sort(unique(comb_gt_summ$comparison))

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
    # layer_order <- sort(unique(combined_gt$comparison))

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





# from: pseudobulk_DE_DESeq2_GSEA_vis_PCA_altdist_genelayercomb_12292023.R
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




# panel B: plot leading edge genes based on their DE summary  stats for select pathways
# pathway_of_interest <- "cotranslational protein targeting to membrane"
# pathway_of_interest <- "mitochondrial translation"
pathway_of_interest <- "modulation of chemical synaptic transmission"
# pathway_of_interest <- "mitochondrial respiratory chain complex assembly"
# pathway_of_interest <- "covalent chromatin modification"
leading_edges_sub1 <- combined_gt[combined_gt$pathway == pathway_of_interest & grepl("diagnosis_cont_ilbd_redef", combined_gt$comparison), ]
leading_edges_sub2 <- combined_gt[combined_gt$pathway == pathway_of_interest & grepl("diagnosis_cont_case_redef", combined_gt$comparison), ]
le_gene_list <- unique( c(unlist(strsplit(leading_edges_sub1$leadingEdge, ', *')),
                          unlist(strsplit(leading_edges_sub2$leadingEdge, ', *'))) )
# read each layer DE file in MTG and SN and plot in scatterplot

st_layers <- c(paste0("Layer",1:6), "WM")
de_all <- c()
for(st_layer in st_layers) {
    contilbd_csv_path <- file.path("/Users/jy1008/Downloads/pseudobulk_DE_DESeq2_V4_12052023", "diagnosis_cont_ilbd_redef")
    contcase_csv_path <- file.path("/Users/jy1008/Downloads/pseudobulk_DE_DESeq2_V4_12052023", "diagnosis_cont_case_redef")

    contilbd_de <- read.csv(file.path(contilbd_csv_path,
                        sprintf("pseudobulk_clinical_deseq2_%s_%s_10152024.csv",
                        "diagnosis_cont_ilbd_redef", st_layer)), row.names = 1)
    contilbd_de[is.na(contilbd_de$padj), "padj"] <- 1
    contilbd_de$gene <- rownames(contilbd_de)
    contcase_de <- read.csv(file.path(contcase_csv_path,
                        sprintf("pseudobulk_clinical_deseq2_%s_%s_10152024.csv",
                        "diagnosis_cont_case_redef", st_layer)), row.names = 1)
    contcase_de[is.na(contcase_de$padj), "padj"] <- 1
    contcase_de$gene <- rownames(contcase_de)

    de_merged <- merge(contilbd_de, contcase_de, by="gene", suffixes = c(".contilbd",".contcase"))
    de_merged$st_layer <- st_layer
    de_all <- rbind(de_all, de_merged)
}
de_filt <- de_all[de_all$gene %in% le_gene_list, ]
de_filt$is_sig <- "neither"
de_filt[de_filt$pvalue.contilbd < 0.05, "is_sig"] <- "contilbd"
de_filt[de_filt$pvalue.contcase < 0.05, "is_sig"] <- "contcase"
de_filt[de_filt$pvalue.contilbd < 0.05 & de_filt$pvalue.contcase < 0.05, "is_sig"] <- "both"
de_filt$label_repel <- de_filt$gene
de_filt[de_filt$is_sig == "neither", "label_repel"] <- ""

pdf(sprintf("earlyPD_contilbd_vs_contcase_redef_de_genes_%s_10152024.pdf", gsub(" ", "_", pathway_of_interest)), width=7, height=5)
for(st_layer in st_layers) {
    de_sub <- de_filt[de_filt$st_layer == st_layer, ]
    pt <- ggplot(de_sub, aes(x=log2FoldChange.contilbd, y=log2FoldChange.contcase, label=label_repel, color=is_sig)) +
            geom_point() +
            geom_text_repel(max.overlaps=Inf) +
            # coord_fixed() +
            geom_abline(slope=1, intercept=0) +
            theme_bw(base_size=20) + labs(title=paste(pathway_of_interest, st_layer))
    print(pt)
}
dev.off()


# scatterplot for all genes, with select few highlighted
de_all_disp <- de_all
de_all_disp$is_sig <- "neither"
de_all_disp[de_all_disp$pvalue.contilbd < 0.05, "is_sig"] <- "contilbd"
de_all_disp[de_all_disp$pvalue.contcase < 0.05, "is_sig"] <- "contcase"
de_all_disp[de_all_disp$pvalue.contilbd < 0.05 & de_all_disp$pvalue.contcase < 0.05, "is_sig"] <- "both"
de_all_disp$basemean_ave <- log((de_all_disp$baseMean.contilbd + de_all_disp$baseMean.contcase)/2 + 1)
disp_gene_list <- c("VGF", "DAP3", "GSTM1", "SFN", "HBB", "FOXJ1", "ADAMTS4", "HSPA6")
de_all_disp$label_repel <- ""
de_all_disp[de_all_disp$gene %in% disp_gene_list, "label_repel"] <- de_all_disp[de_all_disp$gene %in% disp_gene_list, "gene"]
de_all_disp$is_sig_disp <- "bg"
de_all_disp[de_all_disp$gene %in% disp_gene_list, "is_sig_disp"] <- de_all_disp[de_all_disp$gene %in% disp_gene_list, "is_sig"]

# just show the listed genes
de_all_disp <- de_all_disp[de_all_disp$is_sig_disp != "bg", ]

# de_all_disp[abs(de_all_disp$log2FoldChange.contcase) > 2, "label_repel"] <- de_all_disp[abs(de_all_disp$log2FoldChange.contcase) > 2, "gene"]
# de_all_disp[abs(de_all_disp$log2FoldChange.contilbd) > 1, "label_repel"] <- de_all_disp[abs(de_all_disp$log2FoldChange.contilbd) > 1, "gene"]
pdf("earlyPD_contilbd_vs_contcase_redef_de_genes_all_10152024.pdf", width=7, height=5)
for(st_layer in st_layers) {
    de_sub <- de_all_disp[de_all_disp$st_layer == st_layer, ]
    pt <- ggplot(de_sub, aes(x=log2FoldChange.contilbd, y=log2FoldChange.contcase, label=label_repel, color=is_sig_disp)) +
            geom_point() +
            geom_text() +
            # coord_fixed() +
            geom_abline(slope=1, intercept=0) +
            theme_bw(base_size=20) + labs(title=paste("all", st_layer))
    print(pt)
}
dev.off()









library(Seurat)
library(testit)
library(dplyr)
library(ggplot2)

comb_seurat <- readRDS("/mnt/data/projects/ASAP/ST/clustering_pipeline/layer_smoothing_V4_seed101_11302023_replication_test/comb_seurat_post_smoothed_label_s5_11302023.rds")

comb_seurat@meta.data$nctx_temporal_fix <- comb_seurat@meta.data$nctx_temporal
comb_seurat@meta.data$nctx_temporal_fix[comb_seurat@meta.data$nctx_temporal_fix > 3] <- 3
assert(all(colnames(comb_seurat@assays$Spatial@data) == rownames(comb_seurat@meta.data)))


subject_meta <- comb_seurat@meta.data[match(unique(comb_seurat@meta.data$sample_name), comb_seurat@meta.data$sample_name),]

# from pseudobulk_DE_DESeq2_V4_ilbd_pd_relabel_10102024.R
# use a limited set of equal-size cont/ilbd and cont/case comparisons
# controls: keep all
control_redef <- subject_meta[!is.na(subject_meta$nctx_temporal) &
             subject_meta$nctx_temporal == 0 &
             subject_meta$brain_stem_sn == 0 &
             subject_meta$diagnosis == "Control",
             c("sample_name", "diagnosis", "nctx_temporal", "brain_stem_sn")]
# "true" prodromal: nonzero SN, zero MTG
ilbd_redef <- subject_meta[!is.na(subject_meta$nctx_temporal) &
             subject_meta$diagnosis == "ILBD" &
             subject_meta$nctx_temporal == 0 &
             subject_meta$brain_stem_sn > 0 &
             subject_meta$brain_stem_sn < 4,
             c("sample_name", "diagnosis", "nctx_temporal", "brain_stem_sn")]
# "true" late-stage case:
case_redef <- subject_meta[!is.na(subject_meta$nctx_temporal) &
             subject_meta$diagnosis == "Case" &
             subject_meta$nctx_temporal > 1 &
             # subject_meta$nctx_temporal < 4 &
             subject_meta$brain_stem_sn > 3,
             c("sample_name", "diagnosis", "nctx_temporal", "brain_stem_sn")]
control_redef_samps <- as.character(control_redef$sample_name)
ilbd_redef_samps <- as.character(ilbd_redef$sample_name)
case_redef_samps <- as.character(case_redef$sample_name)
comb_seurat@meta.data$diagnosis_redef <- NA
comb_seurat@meta.data[comb_seurat@meta.data$sample_name %in% control_redef_samps, "diagnosis_redef"] <- "control_redef"
comb_seurat@meta.data[comb_seurat@meta.data$sample_name %in% ilbd_redef_samps, "diagnosis_redef"] <- "ilbd_redef"
comb_seurat@meta.data[comb_seurat@meta.data$sample_name %in% case_redef_samps, "diagnosis_redef"] <- "case_redef"

plot_gene_diagnosis_all_layers <- function(gene_list, comb_seurat, outname) {
    meta_data <- comb_seurat@meta.data[, c("sample_name", "nCount_Spatial", "diagnosis", "diagnosis_redef", "smoothed_label_s5", "nctx_temporal_fix")]
    meta_data$diagnosis <- factor(meta_data$diagnosis, levels=c("Control", "ILBD", "Case"))
    meta_data$nctx_temporal_fix <- factor(meta_data$nctx_temporal_fix)
    meta_data$diagnosis_redef <- factor(meta_data$diagnosis_redef, levels=c("control_redef", "ilbd_redef", "case_redef"))
    
    if(length(gene_list) == 1) {
        meta_data <- cbind(meta_data, comb_seurat@assays$Spatial@counts[gene_list[1], ])
        colnames(meta_data)[ncol(meta_data)] <- "gene"
    } else {
        meta_data <- cbind(meta_data, t(as.matrix(comb_seurat@assays$Spatial@counts[gene_list, ])))
        meta_data$gene <- rowSums(meta_data[, gene_list])
        meta_data <- meta_data[, c("sample_name", "nCount_Spatial", "diagnosis", "diagnosis_redef", "smoothed_label_s5", "nctx_temporal_fix", "gene")]
    }
    layers <- sort(unique(meta_data$smoothed_label_s5))
    meta_data <- meta_data[!is.na(meta_data$nctx_temporal_fix), ]

    # get aggregated sum counts
    meta_data_agg <- meta_data %>%
                    group_by(sample_name, diagnosis_redef, smoothed_label_s5, nctx_temporal_fix) %>%
                    summarise(gene_sum = sum(gene), count_sum = sum(nCount_Spatial)) %>%
                    as.data.frame()
    scale_factor <- 10^-4
    lognorm_const <- 1.0
    meta_data_agg$lognorm_gene_sum <- log(meta_data_agg$gene_sum / meta_data_agg$count_sum / scale_factor + lognorm_const)

    group.colors <- c("Layer 1" = "#f43545", "Layer 2"="#ff8901", "Layer 3"="#fad717", "Layer 4"="#00ba71", "Layer 5"="#00c2de", "Layer 6"="#ba55d3", "WM"="#652175")
    pdf(sprintf("pseudobulk_V4_nctx_temporal_fix_vs_%s.pdf", outname))
    meta_sub_agg <- meta_data_agg %>% group_by(smoothed_label_s5, diagnosis_redef) %>%
                        summarise(diag_mean_sum = mean(lognorm_gene_sum),
                                  diag_sd_sum = sd(lognorm_gene_sum)) %>%
                        as.data.frame()
    sp_pt <- ggplot(meta_sub_agg, aes(x=diagnosis_redef, y=diag_mean_sum, color=smoothed_label_s5, group=smoothed_label_s5)) +
    geom_line(position=position_dodge(0.3)) +
    geom_point(position=position_dodge(0.3)) +
    geom_errorbar(aes(ymin=diag_mean_sum-diag_sd_sum,
                      ymax=diag_mean_sum+diag_sd_sum),
                  position=position_dodge(0.3),
                  width=.5, linewidth=1) +
    scale_color_manual(values=group.colors) + theme_bw()
    print(sp_pt)

    dev.off()
}
plot_gene_diagnosis_all_layers(c("GSTM1"), comb_seurat, outname="earlyPD_test_GSTM1_10152024")
plot_gene_diagnosis_all_layers(c("SERPINA3"), comb_seurat, outname="earlyPD_test_SERPINA3_10152024")
plot_gene_diagnosis_all_layers(c("NPAS4"), comb_seurat, outname="earlyPD_test_NPAS4_10152024")
plot_gene_diagnosis_all_layers(c("VGF"), comb_seurat, outname="earlyPD_test_VGF_10152024")
plot_gene_diagnosis_all_layers(c("DAP3"), comb_seurat, outname="earlyPD_test_DAP3_10152024")


plot_gene_diagnosis_all_layers_bars <- function(gene_list, comb_seurat, outname) {
    meta_data <- comb_seurat@meta.data[, c("sample_name", "nCount_Spatial", "diagnosis", "diagnosis_redef", "smoothed_label_s5", "nctx_temporal_fix")]
    meta_data$diagnosis <- factor(meta_data$diagnosis, levels=c("Control", "ILBD", "Case"))
    meta_data$nctx_temporal_fix <- factor(meta_data$nctx_temporal_fix)
    meta_data$diagnosis_redef <- factor(meta_data$diagnosis_redef, levels=c("control_redef", "ilbd_redef", "case_redef"))
    
    if(length(gene_list) == 1) {
        meta_data <- cbind(meta_data, comb_seurat@assays$Spatial@counts[gene_list[1], ])
        colnames(meta_data)[ncol(meta_data)] <- "gene"
    } else {
        meta_data <- cbind(meta_data, t(as.matrix(comb_seurat@assays$Spatial@counts[gene_list, ])))
        meta_data$gene <- rowSums(meta_data[, gene_list])
        meta_data <- meta_data[, c("sample_name", "nCount_Spatial", "diagnosis", "diagnosis_redef", "smoothed_label_s5", "nctx_temporal_fix", "gene")]
    }
    layers <- sort(unique(meta_data$smoothed_label_s5))
    meta_data <- meta_data[!is.na(meta_data$nctx_temporal_fix), ]

    # get aggregated sum counts
    meta_data_agg <- meta_data %>%
                    group_by(sample_name, diagnosis_redef, smoothed_label_s5, nctx_temporal_fix) %>%
                    summarise(gene_sum = sum(gene), count_sum = sum(nCount_Spatial)) %>%
                    as.data.frame()
    scale_factor <- 10^-4
    lognorm_const <- 1.0
    meta_data_agg$lognorm_gene_sum <- log(meta_data_agg$gene_sum / meta_data_agg$count_sum / scale_factor + lognorm_const)

    group.colors <- c("Layer 1" = "#f43545", "Layer 2"="#ff8901", "Layer 3"="#fad717", "Layer 4"="#00ba71", "Layer 5"="#00c2de", "Layer 6"="#ba55d3", "WM"="#652175")
    pdf(sprintf("pseudobulk_V4_nctx_temporal_fix_bar_vs_%s.pdf", outname), width=7, height=5)
    
    sp_pt <- ggplot(meta_data_agg[!is.na(meta_data_agg$diagnosis_redef), ], aes(x=diagnosis_redef, y=lognorm_gene_sum, fill=smoothed_label_s5)) +
             geom_boxplot() +
             # geom_boxplot(outlier.color=NA) +
             # geom_jitter() +
             scale_fill_manual(values=group.colors) +
             theme_bw(base_size=20) +
             labs(title=paste(gene_list, collapse=" "))
    print(sp_pt)

    dev.off()
}
plot_gene_diagnosis_all_layers_bars(c("GSTM1"), comb_seurat, outname="earlyPD_test_GSTM1_10152024")
plot_gene_diagnosis_all_layers_bars(c("SERPINA3"), comb_seurat, outname="earlyPD_test_SERPINA3_10152024")
plot_gene_diagnosis_all_layers_bars(c("NPAS4"), comb_seurat, outname="earlyPD_test_NPAS4_10152024")
plot_gene_diagnosis_all_layers_bars(c("VGF"), comb_seurat, outname="earlyPD_test_VGF_10152024")
plot_gene_diagnosis_all_layers_bars(c("DAP3"), comb_seurat, outname="earlyPD_test_DAP3_10152024")
plot_gene_diagnosis_all_layers_bars(c("EGFR"), comb_seurat, outname="earlyPD_test_EGFR_10152024")
plot_gene_diagnosis_all_layers_bars(c("STAR"), comb_seurat, outname="earlyPD_test_STAR_10152024")
plot_gene_diagnosis_all_layers_bars(c("SORCS3"), comb_seurat, outname="earlyPD_test_SORCS3_10152024")
plot_gene_diagnosis_all_layers_bars(c("RPUSD3"), comb_seurat, outname="earlyPD_test_RPUSD3_10152024")
plot_gene_diagnosis_all_layers_bars(c("CHCHD1"), comb_seurat, outname="earlyPD_test_CHCHD1_10152024")
plot_gene_diagnosis_all_layers_bars(c("MRPS17"), comb_seurat, outname="earlyPD_test_MRPS17_10152024")
plot_gene_diagnosis_all_layers_bars(c("MTRNR2L13"), comb_seurat, outname="earlyPD_test_MTRNR2L13_10152024")
plot_gene_diagnosis_all_layers_bars(c("AD000090.1"), comb_seurat, outname="earlyPD_test_AD000090.1_10152024")
plot_gene_diagnosis_all_layers_bars(c("ADAMTS4"), comb_seurat, outname="earlyPD_test_ADAMTS4_10152024")
plot_gene_diagnosis_all_layers_bars(c("CCL3L1"), comb_seurat, outname="earlyPD_test_CCL3L1_10152024")
plot_gene_diagnosis_all_layers_bars(c("SFN"), comb_seurat, outname="earlyPD_test_SFN_10152024")
plot_gene_diagnosis_all_layers_bars(c("APOL4"), comb_seurat, outname="earlyPD_test_APOL4_10152024")
plot_gene_diagnosis_all_layers_bars(c("HSPA6"), comb_seurat, outname="earlyPD_test_HSPA6_10152024")
plot_gene_diagnosis_all_layers_bars(c("SLCO4A1"), comb_seurat, outname="earlyPD_test_SLCO4A1_10152024")
plot_gene_diagnosis_all_layers_bars(c("GJC3"), comb_seurat, outname="earlyPD_test_GJC3_10152024")
plot_gene_diagnosis_all_layers_bars(c("CXCL10"), comb_seurat, outname="earlyPD_test_CXCL10_10152024")
plot_gene_diagnosis_all_layers_bars(c("HBB"), comb_seurat, outname="earlyPD_test_HBB_10152024")
plot_gene_diagnosis_all_layers_bars(c("CD36"), comb_seurat, outname="earlyPD_test_CD36_10152024")
plot_gene_diagnosis_all_layers_bars(c("KLK7"), comb_seurat, outname="earlyPD_test_KLK7_10152024")
plot_gene_diagnosis_all_layers_bars(c("KANSL1"), comb_seurat, outname="earlyPD_test_KANSL1_10152024")
plot_gene_diagnosis_all_layers_bars(c("SNCG"), comb_seurat, outname="earlyPD_test_SNCG_10152024")
plot_gene_diagnosis_all_layers_bars(c("SNCA"), comb_seurat, outname="earlyPD_test_SNCA_10152024")
plot_gene_diagnosis_all_layers_bars(c("ECM2"), comb_seurat, outname="earlyPD_test_ECM2_10152024")
plot_gene_diagnosis_all_layers_bars(c("CHI3L2"), comb_seurat, outname="earlyPD_test_CHI3L2_10152024")
plot_gene_diagnosis_all_layers_bars(c("TGFB2"), comb_seurat, outname="earlyPD_test_TGFB2_10152024")
plot_gene_diagnosis_all_layers_bars(c("IGKC"), comb_seurat, outname="earlyPD_test_IGKC_10152024")
plot_gene_diagnosis_all_layers_bars(c("FOXJ1"), comb_seurat, outname="earlyPD_test_FOXJ1_10152024")

# 12/7/2025
# Group on the x-axis first by layer, then by diagnosis
plot_gene_layers_then_diagnosis_bars <- function(gene_list, comb_seurat, outname) {
    meta_data <- comb_seurat@meta.data[, c("sample_name", "nCount_Spatial", "diagnosis", "diagnosis_redef", "smoothed_label_s5", "nctx_temporal_fix")]
    meta_data$diagnosis <- factor(meta_data$diagnosis, levels=c("Control", "ILBD", "Case"))
    meta_data$nctx_temporal_fix <- factor(meta_data$nctx_temporal_fix)
    meta_data$diagnosis_redef <- factor(meta_data$diagnosis_redef, levels=c("control_redef", "ilbd_redef", "case_redef"))
    
    if(length(gene_list) == 1) {
        meta_data <- cbind(meta_data, comb_seurat@assays$Spatial@counts[gene_list[1], ])
        colnames(meta_data)[ncol(meta_data)] <- "gene"
    } else {
        meta_data <- cbind(meta_data, t(as.matrix(comb_seurat@assays$Spatial@counts[gene_list, ])))
        meta_data$gene <- rowSums(meta_data[, gene_list])
        meta_data <- meta_data[, c("sample_name", "nCount_Spatial", "diagnosis", "diagnosis_redef", "smoothed_label_s5", "nctx_temporal_fix", "gene")]
    }
    layers <- sort(unique(meta_data$smoothed_label_s5))
    meta_data <- meta_data[!is.na(meta_data$nctx_temporal_fix), ]

    # get aggregated sum counts
    meta_data_agg <- meta_data %>%
                    group_by(sample_name, diagnosis_redef, smoothed_label_s5, nctx_temporal_fix) %>%
                    summarise(gene_sum = sum(gene), count_sum = sum(nCount_Spatial)) %>%
                    as.data.frame()
    scale_factor <- 10^-4
    lognorm_const <- 1.0
    meta_data_agg$lognorm_gene_sum <- log(meta_data_agg$gene_sum / meta_data_agg$count_sum / scale_factor + lognorm_const)

    group.colors <- c("Layer 1" = "#f43545", "Layer 2"="#ff8901", "Layer 3"="#fad717", "Layer 4"="#00ba71", "Layer 5"="#00c2de", "Layer 6"="#ba55d3", "WM"="#652175")
    pdf(sprintf("pseudobulk_V4_nctx_temporal_fix_bar2_vs_%s.pdf", outname), width=7, height=5)
    
    sp_pt <- ggplot(
        meta_data_agg[!is.na(meta_data_agg$diagnosis_redef), ],
        aes(x=smoothed_label_s5, y=lognorm_gene_sum, fill=diagnosis_redef)) +
        geom_boxplot() +
        # geom_boxplot(outlier.color=NA) +
        # geom_jitter() +
        # scale_fill_manual(values=group.colors) +
        scale_fill_manual(values = c(
            "control_redef" = "#1f77b4",
            "ilbd_redef" = "#ff7f0e",
            "case_redef" = "#2ca02c"
        )) +
        theme_bw(base_size=20) +
        labs(title=paste(gene_list, collapse=" "))
    print(sp_pt)

    dev.off()
}
plot_gene_layers_then_diagnosis_bars(c("VGF"), comb_seurat, outname="earlyPD_VGF_12072025")
plot_gene_layers_then_diagnosis_bars(c("DAP3"), comb_seurat, outname="earlyPD_DAP3_12072025")
plot_gene_layers_then_diagnosis_bars(c("GSTM1"), comb_seurat, outname="earlyPD_GSTM1_12072025")
plot_gene_layers_then_diagnosis_bars(c("SFN"), comb_seurat, outname="earlyPD_SFN_12072025")
plot_gene_layers_then_diagnosis_bars(c("HBB"), comb_seurat, outname="earlyPD_HBB_12072025")
plot_gene_layers_then_diagnosis_bars(c("FOXJ1"), comb_seurat, outname="earlyPD_FOXJ1_12072025")
plot_gene_layers_then_diagnosis_bars(c("ADAMTS4"), comb_seurat, outname="earlyPD_ADAMTS4_12072025")
plot_gene_layers_then_diagnosis_bars(c("HSPA6"), comb_seurat, outname="earlyPD_HSPA6_12072025")

plot_gene_layers_then_diagnosis_bars(c("ADRA2A"), comb_seurat, outname="earlyPD_ADRA2A_12072025")
plot_gene_layers_then_diagnosis_bars(c("GRIN2D"), comb_seurat, outname="earlyPD_GRIN2D_12072025")
plot_gene_layers_then_diagnosis_bars(c("GRIN2B"), comb_seurat, outname="earlyPD_GRIN2B_12072025")
plot_gene_layers_then_diagnosis_bars(c("GRM7"), comb_seurat, outname="earlyPD_GRM7_12072025")
plot_gene_layers_then_diagnosis_bars(c("GRM2"), comb_seurat, outname="earlyPD_GRM2_12072025")
plot_gene_layers_then_diagnosis_bars(c("SNCA"), comb_seurat, outname="earlyPD_SNCA_12072025")
plot_gene_layers_then_diagnosis_bars(c("SNCG"), comb_seurat, outname="earlyPD_SNCG_12072025")

# 2/1/2026
# for a particular pathway
# get a list of leading edge genes from the redef comparisons
# cross reference with gene DE
# identify a gene involved in the pathway which is downregulated in ILBD and upregulated in PD
pathway_of_interest <- "covalent chromatin modification"
gsea1 <- read.delim("~/de_genes_DESeq2_pseudobulk_diagnosis_cont_ilbd_redef_GO-BP_10102024.fGSEA.xls", header = TRUE, stringsAsFactors = FALSE)
gsea2 <- read.delim("~/de_genes_DESeq2_pseudobulk_diagnosis_cont_case_redef_GO-BP_10102024.fGSEA.xls", header = TRUE, stringsAsFactors = FALSE)

dfs <- list(gsea1, gsea2)

gene_counts_list <- lapply(dfs, function(df) {
  df_sub <- df[df$pathway == pathway_of_interest & df$padj < 0.05, ]

  genes <- trimws(unlist(strsplit(df_sub$leadingEdge, ",\\s*")))
  genes <- genes[genes != "" & !is.na(genes)]

  as.data.frame(table(genes))
})

combined_gene_counts <- Reduce(function(x, y) {
  merge(x, y, by = "genes", all = TRUE)
}, gene_counts_list)

combined_gene_counts[is.na(combined_gene_counts)] <- 0
combined_gene_counts$total <- rowSums(combined_gene_counts[-1])
# > combined_gene_counts
#        genes Freq.x Freq.y total
# 1       AASS      1      0     1
# 2     ACTL6A      2      3     5
# 3     ACTL6B      1      0     1
# 4      AKAP8      1      0     1
# 5     AKAP8L      1      6     7

# read DE genes

# folder containing your CSVs

# function to read each file and add a 'layer' column
read_with_layer <- function(file) {
  df <- read.csv(file, stringsAsFactors = FALSE)
  
  # extract Layer info from filename
  layer <- sub(".*_redef_(Layer[0-9]+|WM)_.*\\.csv$", "\\1", basename(file))
  
  df$layer <- layer
  if (colnames(df)[1] == "X") colnames(df)[1] <- "gene"
  df
}

# read and combine all files
files <- list.files(
  path = "~/pseudobulk_clinical_deseq2_diagnosis_cont_ilbd_redef_10102024",
  pattern = "pseudobulk_clinical_deseq2.*\\.csv$",
  full.names = TRUE
)
de_ilbd <- do.call(rbind, lapply(files, read_with_layer))

files <- list.files(
  path = "~/pseudobulk_clinical_deseq2_diagnosis_cont_case_redef_10102024",
  pattern = "pseudobulk_clinical_deseq2.*\\.csv$",
  full.names = TRUE
)
de_case <- do.call(rbind, lapply(files, read_with_layer))

# add source column
de_ilbd$group <- "ilbd"
de_case$group <- "case"

# combine
combined_de <- rbind(de_ilbd, de_case)
# filter for leading edge genes
combined_de <- combined_de[combined_de$gene %in% as.character(combined_gene_counts$gene), ]
# combined_de <- combined_de[!is.na(combined_de$padj) & combined_de$padj < 0.05, ]

library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)

wide_de <- combined_de %>%
  select(gene, layer, group, log2FoldChange) %>%
  pivot_wider(
    names_from = group,
    values_from = log2FoldChange
  )
wide_de <- wide_de %>%
  mutate(diff = abs(ilbd - case))
top_genes <- wide_de %>%
  group_by(gene) %>%
  summarise(max_diff = max(diff, na.rm = TRUE)) %>%
  arrange(desc(max_diff))
top_n_genes <- top_genes %>% slice_head(n = 40)
plot_df <- combined_de %>% filter(gene %in% top_n_genes$gene)

plot_df <- plot_df %>%
  mutate(
    layer = factor(layer, levels = c(paste0("Layer", 1:6), "WM")),
    gene = factor(gene, levels = top_n_genes$gene)  # preserve order
  )

p <- ggplot(plot_df, aes(x = layer, y = gene, fill = log2FoldChange)) +
  geom_tile(color = "white") +
  facet_wrap(~ group) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal()
ggsave("temp_gene_layer_heatmap_chromatin_02012026.pdf", p, width = 10, height = 8)


genes <- c("MIER2", "KDM6B", "JADE3", "USP49", "KMT2D", "TRERF1", "BAZ1A", "BCOR", "JARID2", "GATA2", "SIRT7", "N6AMT1", "VEGFA", "DNMT3B", "BEND3", "PADI4", "KANSL1L", "JDP2", "FLCN", "JMJD6", "MECOM", "EYA1", "EYA2", "UHRF1", "CCNA2")
for (gene in genes) {
  outname <- paste0("earlyPD_", gene, "_02012026")
  plot_gene_layers_then_diagnosis_bars(c(gene), comb_seurat, outname = outname)
}
# JPD2, JMJD6, KANSL1L, KMT2D, MIER2, SIRT7*


pathway_of_interest <- "positive regulation of immune response"
# ...
ggsave("temp_gene_layer_heatmap_immune_02012026.pdf", p, width = 10, height = 8)
genes <- c("MICB", "C2")
for (gene in genes) {
  outname <- paste0("earlyPD_", gene, "_02012026")
  plot_gene_layers_then_diagnosis_bars(c(gene), comb_seurat, outname = outname)
}
