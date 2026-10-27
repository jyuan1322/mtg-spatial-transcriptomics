# function: in a vector of numbers, find local maxima
# https://stackoverflow.com/a/6836924
localMaxima <- function(x) {
    # Use -Inf instead if x is numeric (non-integer)
    y <- diff(c(-.Machine$integer.max, x)) > 0L
    rle(y)$lengths
    y <- cumsum(rle(y)$lengths)
    y <- y[seq.int(1L, length(y), 2L)]
    if (x[[1]] == x[[2]]) {
    y <- y[-1]
    }
    y
}

gsea_clinical_per_layer_filelist <- function(effects_path, clinical_feat, layers) {

    de_files_list <- list()
    for(st_layer in layers) {
        clinical_feat_file <- sprintf("pseudobulk_clinical_deseq2_%s_%s_[0-9]+\\.csv$", clinical_feat, gsub(' ', '', st_layer))
        
        file_match <- list.files( path=file.path( effects_path, clinical_feat), pattern=clinical_feat_file )
        if(length(file_match) == 1) {
            clinical_feat_file <- file_match[1]
            print(paste("Reading file", clinical_feat_file))
            de_files_list[st_layer] <- file.path(effects_path, clinical_feat, clinical_feat_file)
        } else {
            print(paste("Failed to read file", clinical_feat_file))
            next()
        }
    }
    return(de_files_list)
}
# layer-wise DE requires different file handling
gsea_layer_vs_layer_filelist <- function(effects_path, layers, contsonly=FALSE) {
    
    de_files_list <- list()
    for(st_layer in layers) {
        layer_file <- sprintf("pseudobulk_layer_deseq2_is_%s_contsonly%s_[0-9]+\\.csv$", gsub(' ', '', st_layer), contsonly)
        
        file_match <- list.files( path=file.path( effects_path, "layer_vs_layer"), pattern=layer_file )
        if(length(file_match) == 1) {
            layer_file <- file_match[1]
            print(paste("Reading file", layer_file))
            de_files_list[st_layer] <- file.path(effects_path, "layer_vs_layer", layer_file)
        } else {
            print(paste("Failed to read file", layer_file))
            next()
        }
    }
    return(de_files_list)
}
