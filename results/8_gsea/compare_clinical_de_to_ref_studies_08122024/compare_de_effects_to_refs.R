library(ggplot2)
library(ggrepel)


# https://stackoverflow.com/a/7549819
lm_eqn <- function(df, x, y){
    m <- lm(as.formula(paste(y, '~', x)), df);
    eq <- sprintf("y == %s + %s * x, r^2 = %s \n SP corr = %s",
              format(unname(coef(m)[1]), digits = 2),
              format(unname(coef(m)[2]), digits = 2),
              format(summary(m)$r.squared, digits = 3),
              format(cor(df[, x], df[, y], method="spearman"), digits = 3))
    # as.character(as.expression(eq));
    return(eq)
}

count_point_quadrants <- function(x, y) {
    temp_xy <- data.frame(x=x, y=y)
    q1 <- nrow(temp_xy[temp_xy$x > 0 & temp_xy$y > 0, ])
    q2 <- nrow(temp_xy[temp_xy$x < 0 & temp_xy$y > 0, ])
    q3 <- nrow(temp_xy[temp_xy$x < 0 & temp_xy$y < 0, ])
    q4 <- nrow(temp_xy[temp_xy$x > 0 & temp_xy$y < 0, ])
    out_string <- sprintf("Q1:%s, Q2:%s, Q3:%s, Q4:%s", q1, q2, q3, q4)
    return(out_string)
}

# ref df should have 4 columns: gene, log2FoldChange, pvalue, padj
reference_scatter <- function(clin_feat, ref_df, ref_df_name) {
    
    # clin_feat <- "nctx_temporal_fix"
    de_dir <- sprintf("/Users/jy1008/Downloads/pseudobulk_DE_DESeq2_V4_12052023/%s", clin_feat)

    sp_layers <- c(paste0("Layer", 1:6), "WM")

    de_df_list <- list()
    for(sp_layer in sp_layers) {
        de_file <- sprintf("pseudobulk_clinical_deseq2_%s_%s_12052023.csv", clin_feat, sp_layer)
        if(file.exists(file.path(de_dir, de_file))) {
            de_df <- read.csv(file.path(de_dir, de_file))
            colnames(de_df)[1] <- "gene"
            de_df[is.na(de_df$padj), "padj"] <- 1
            if(clin_feat == "last_mmse_test_score") {
                de_df$log2FoldChange <- -1 * de_df$log2FoldChange
            }
            
            de_df_list[[sp_layer]] <- de_df
        } else {
            print(paste("skipping", de_file))
        }
    }
    
    pdf(sprintf("compare_de_to_ref_40478_2020_932_%s_%s_08112024.pdf", clin_feat, ref_df_name))
    for(sp_layer in names(de_df_list)) {
        de_df <- de_df_list[[sp_layer]]
        de_df_filt <- de_df[de_df$pvalue < 0.1, c("gene", "log2FoldChange", "baseMean", "padj")]
        ref_df_filt <- ref_df[ref_df$pvalue < 0.1, c("gene", "log2FoldChange")]
        temp_df <- merge(de_df_filt, ref_df_filt, by="gene", suffixes=c(".ours", ".ref"))
        temp_df$nlog10_padj <- -log10(temp_df$padj)

        sp_pt_label <- paste(
            lm_eqn(temp_df, "log2FoldChange.ours", "log2FoldChange.ref"), '\n',
            count_point_quadrants(temp_df$log2FoldChange.ours, temp_df$log2FoldChange.ref))
        
        sp_pt <- ggplot(temp_df, aes(x=log2FoldChange.ours, y=log2FoldChange.ref, label=gene, size=baseMean, colour=nlog10_padj)) +
        geom_point() +
        scale_colour_viridis_c() +
        geom_smooth(method = "lm", formula = y ~ x) +
        geom_text_repel(hjust=0, vjust=0) +
        ggtitle(paste(sp_layer, clin_feat)) + theme_bw()
        
        sp_pt <- sp_pt + annotate(geom='text', label=sp_pt_label, x=-Inf, y=Inf, hjust=0, vjust=1)
        print(sp_pt)
    }
    dev.off()
}

# read file from reference
ref_file <- "/Users/jy1008/Downloads/40478_2020_932_MOESM4_ESM.csv"
ref_df_name <- "40478_Nido_2020"
ref_df <- read.csv(ref_file)
colnames(ref_df)[colnames(ref_df) == "GeneSymbol"] <- "gene"
reference_scatter(clin_feat="nctx_temporal_fix", ref_df=ref_df, ref_df_name=ref_df_name)
reference_scatter(clin_feat="last_mmse_test_score", ref_df=ref_df, ref_df_name=ref_df_name)



ref_file <- "/Users/jy1008/Downloads/41467_2024_47867_MOESM10_ESM.csv"
ref_df_name <- "41467_Flones_2024_CIPDm1"
ref_df <- read.csv(ref_file)
ref_df <- ref_df[, c("gene_name", "log2FoldChange_CI.PDmodel1", "pvalue_CI.PDmodel1", "padj_CI.PDmodel1")]
colnames(ref_df) <- c("gene", "log2FoldChange", "pvalue", "padj")
reference_scatter(clin_feat="nctx_temporal_fix", ref_df=ref_df, ref_df_name=ref_df_name)
reference_scatter(clin_feat="last_mmse_test_score", ref_df=ref_df, ref_df_name=ref_df_name)


ref_file <- "/Users/jy1008/Downloads/41467_2024_47867_MOESM10_ESM.csv"
ref_df_name <- "41467_Flones_2024_CIPDm2"
ref_df <- read.csv(ref_file)
ref_df <- ref_df[, c("gene_name", "log2FoldChange_CI.PDmodel2", "pvalue_CI.PDmodel2", "padj_CI.PDmodel2")]
colnames(ref_df) <- c("gene", "log2FoldChange", "pvalue", "padj")
reference_scatter(clin_feat="nctx_temporal_fix", ref_df=ref_df, ref_df_name=ref_df_name)
reference_scatter(clin_feat="last_mmse_test_score", ref_df=ref_df, ref_df_name=ref_df_name)
