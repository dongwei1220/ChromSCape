## Authors : Pacôme Prompsy Title : Wrappers & functions to preprocess & reduce
## dimensionality of sc Matrix Description : Funtions to load, filter, normalize
## & reduce sc Epigenetic Matrices prior to analysis

# Wrapper to create the single cell experiment from sparse datamatrix and annot,
# remove Zero count Features then Cells and calculate QC Metrics
# (scater::calculateQCMetrics)
#' Title
#'
#' @param datamatrix 
#' @param annot 
#' @param remove_zero_cells 
#' @param remove_zero_features 
#'
#' @return
#' @export
#'
#' @examples
create_scExp <- function(datamatrix, annot, remove_zero_cells = TRUE, remove_zero_features = TRUE, remove_non_canonical = TRUE, remove_chr_M = TRUE) {
    
    stopifnot(is.data.frame(annot), remove_zero_cells %in% c(T, F), remove_zero_features %in% c(T, F))
    
    if (ncol(datamatrix) != nrow(annot)) 
        stop("ChromSCape::create_scExp - datamatrix and annot should contain the same number of cells")
    if (length(match(c("cell_id", "sample_id"), colnames(annot))) < 2) 
        stop("ChromSCape::create_scExp - annot should contain cell_id & sample_id as column names")
  
    scExp <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = datamatrix), colData = annot)
    
    if (has_genomic_coordinates(scExp)) {
        if (remove_non_canonical) {
            # Removing non-canonical chromosomes
            splitID <- sapply(rownames(scExp), function(x) strsplit(as.character(x),
                                                                    split = "_|:|-", fixed = F))
            normal_chr <- which(sapply(splitID, length) <= 3)  # weird chromosomes contain _/:/- in the name
            nrow_init = nrow(scExp)
            scExp <- scExp[normal_chr, ]
            if (length(normal_chr) < nrow_init && verbose) 
                cat("ChromSCape::create_scExp -",nrow(scExp) - length(normal_chr), "non canonical regions were removed.")
        }
        if (remove_chr_M) {
            # Remove chrM from mat if it is inside
            chrM_regions = grep("chrM", rownames(scExp))
            if (length(chrM_regions) > 0) {
                scExp <- scExp[-chrM_regions, ]
                if (verbose) 
                  cat("ChromSCape::create_scExp -",length(chrM_regions), "chromosome M regions were removed.")
            }
        }
    }
    dim_b = dim(scExp)
    if (remove_zero_features) 
        scExp <- scExp[(Matrix::rowSums(counts(scExp) > 0) > 0), ]  # remove windows that do not have any read in any cells
    if (remove_zero_cells) 
        scExp <- scExp[, (Matrix::colSums(counts(scExp) > 0) > 0)]  # remove cells that do not have any read in any cells
    
    if (dim(scExp)[2] != dim_b[2]) 
        cat("ChromSCape::create_scExp -",dim_b[2] - dim(scExp)[2], "cells with 0 signals were removed.")
    if (dim(scExp)[1] != dim_b[1]) 
        cat("ChromSCape::create_scExp -",dim_b[1] - dim(scExp)[1], "features with 0 signals were removed.")
    
    
    if (has_genomic_coordinates(scExp)) {
        rows = rownames(scExp)
        rowRanges(scExp) = get_genomic_coordinates(scExp)
        rownames(scExp) = rows
    }
    
    scExp <- scater::calculateQCMetrics(scExp)
    
    return(scExp)
}

# Function to filter out cells & features from sparse matrix based on total
# count per cell, number of cells 'ON' (count >= 2) in features and top covered
# cells that might be doublets
#' Title
#'
#' @param scExp 
#' @param min_cov_cell 
#' @param quant_removal 
#' @param percentMin 
#' @param bin_min_count 
#'
#' @return
#' @export
#'
#' @examples
filter_scExp <- function(scExp, min_cov_cell = 1600, quant_removal = 95, percentMin = 1, bin_min_count = 2, verbose = T) {
    
    stopifnot(is(scExp, "SingleCellExperiment"), is.numeric(min_cov_cell), is.numeric(quant_removal), is.numeric(percentMin), is.numeric(bin_min_count), 
        verbose %in% c(F, T))
    
    if (is.null(scExp)) 
        warning("ChromSCape::filter_scExp - Please specify a SingleCellExperiment")
    
    cellCounts = Matrix::colSums(counts(scExp))
    
    thresh <- quantile(cellCounts, probs = seq(0, 1, 0.01))
    
    sel1000 = (cellCounts > 1000 & cellCounts <= thresh[quant_removal + 1])
    sel <- (cellCounts > min_cov_cell & cellCounts <= thresh[quant_removal + 1])
    
    if (verbose) 
        cat("ChromSCape::filter_scExp -",length(which(sel)), "cells pass the threshold of", min_cov_cell, "minimum reads and are lower than the ", quant_removal, 
            "th centile of library size ~=", round(thresh[quant_removal + 1]), "reads.\n")
    
    # Create the binary matrix of cells >= 1000 counts if count >= bin_min_count, count = 1 else
    bina_counts <- counts(scExp)[,sel1000]
    sel_below_2 <- (bina_counts < bin_min_count)
    bina_counts[sel_below_2] <- 0
    sel_above_1 = (bina_counts >= bin_min_count)
    bina_counts[sel_above_1] <- 1
    
    nCells_in_feature = Matrix::rowSums(bina_counts)
    fixedFeature <- names(which(nCells_in_feature > ((percentMin/100) * (ncol(bina_counts)))))  # Feature selection
    if (verbose) 
        cat("ChromSCape::filter_scExp -",length(fixedFeature), "features pass the threshold of", percentMin, " % of total cells 'ON', representing a minimum of", 
            round((percentMin/100) * (ncol(bina_counts))), "cells.\n")
    
    scExp <- scExp[, sel]
    scExp <- scExp[fixedFeature, ]
    scExp <- scater::calculateQCMetrics(scExp)
    return(scExp)
}

has_genomic_coordinates <- function(scExp) {
    stopifnot(is(scExp, "SingleCellExperiment"), !is.null(rownames(scExp)))
    ID <- rownames(scExp)[1:min(10, length(rownames(scExp)))]
    chr <- unlist(lapply(strsplit(ID, split = "_|-|:"), FUN = function(x) unlist(x)[1]))
    
    if (length(grep("chr|(^\\d+$|^X$|^Y$)", chr[1:min(10, length(chr))], ignore.case = T)) >= min(10, length(chr))) 
        return(T) else return(F)
}

get_genomic_coordinates <- function(scExp) {
    
    stopifnot(is(scExp, "SingleCellExperiment"))
    if (!has_genomic_coordinates(scExp)) 
        stop("Feature names are not genomic coordinates")
    
    ID <- rownames(scExp)
    chr <- unlist(lapply(strsplit(ID, split = "_|-|:"), FUN = function(x) unlist(x)[1]))
    start <- unlist(lapply(strsplit(ID, split = "_|-|:"), FUN = function(x) unlist(x)[2]))
    end <- unlist(lapply(strsplit(ID, split = "_|-|:"), FUN = function(x) unlist(x)[3]))
    feature <- GRanges(data.frame(ID = ID, chr = chr, start = start, end = end))
    return(feature)
}

exclude_features_scExp <- function(scExp, features_to_exclude, by = c("region", "feature_name"), verbose = T) {
    
    stopifnot(is(scExp, "SingleCellExperiment"), is.data.frame(features_to_exclude), is.character(by[1]))
    if (!by[1] %in% c("region", "feature_name")) 
        stop("ChromSCape::exclude_features_scExp - by must be either 'region' or 'feature_name'")
    
    if (by[1] == "region") {
        if (!has_genomic_coordinates(scExp)) 
            stop("ChromSCape::exclude_features_scExp - Feature names are not genomic coordinates")
        regions <- rowRanges(scExp)
        colnames(features_to_exclude)[1:3] = c("chr", "start", "stop")
        excl_gr <- makeGRangesFromDataFrame(features_to_exclude, ignore.strand = TRUE, seqnames.field = c("chr"), start.field = c("start"), 
            end.field = c("stop"))
        ovrlps <- as.data.frame(findOverlaps(regions, excl_gr))[, 1]
        if (length(unique(ovrlps) > 0)) 
            scExp <- scExp[-unique(ovrlps), ]
        if (verbose) 
            cat("ChromSCape::exclude_features_scExp - Removed", length(unique(ovrlps)), "regions from the analysis.\n")
    }
    if (by[1] == "feature_name") {
        if (has_genomic_coordinates(scExp)) 
            warning("ChromSCape::exclude_features_scExp - Excluding by feature name while object feature names are genomic coordinates !")
        features = rownames(scExp)
        features_to_exclude = as.character(features_to_exclude[, 1])
        ovrlps = intersect(features, features_to_exclude)
        if (length(unique(ovrlps) > 0)) 
            scExp <- scExp[-which(rownames(scExp) %in% ovrlps), ]
        if (verbose) 
            cat("ChromSCape::exclude_features_scExp - Removed", length(unique(ovrlps)), " features from the analysis.\n")
    }
    return(scExp)
}

preprocess_TPM <- function(scExp) {
    
    size = BiocGenerics::width(SummarizedExperiment::rowRanges(scExp))
    normcounts(scExp) = SingleCellExperiment::counts(scExp)/size
    normcounts(scExp) = 10^6 * Matrix::t(Matrix::t(normcounts(scExp))/Matrix::colSums(normcounts(scExp)))
    return(scExp)
}

preprocess_RPKM <- function(scExp) {
    
    normcounts(scExp) = 10^9 * Matrix::t(Matrix::t(SingleCellExperiment::counts(scExp))/Matrix::colSums(SingleCellExperiment::counts(scExp)))
    size = BiocGenerics::width(SummarizedExperiment::rowRanges(scExp))
    normcounts(scExp) = normcounts(scExp)/size
    
    return(scExp)
}

preprocess_CPM <- function(scExp) {
    normcounts(scExp) = 10^6 * Matrix::t(Matrix::t(SingleCellExperiment::counts(scExp))/Matrix::colSums(SingleCellExperiment::counts(scExp)))
    return(scExp)
}

preprocess_feature_size_only <- function(scExp) {
    size = BiocGenerics::width(SummarizedExperiment::rowRanges(scExp))
    normcounts(scExp) = SingleCellExperiment::counts(scExp)/size
    return(scExp)
}

normalize_scExp <- function(scExp, type = c("RPKM", "CPM", "TPM", "feature_size_only")) {
    stopifnot(type[1] %in% c("RPKM", "CPM", "TPM", "feature_size_only"), is(scExp, "SingleCellExperiment"))
    if (!has_genomic_coordinates(scExp)) {
        warning("ChromSCape::normalize_scExp - Switching to CPM normalization as features are not genomic coordinates.")
        type = "CPM"
    }
    switch(type[1], RPKM = return(preprocess_RPKM(scExp)), TPM = return(preprocess_TPM(scExp)), feature_size_only = return(preprocess_feature_size_only(scExp)), 
        CPM = return(preprocess_CPM(scExp)))
}

feature_annotation_scExp <- function(scExp, reference_annotation) {
    stopifnot(is(scExp, "SingleCellExperiment"))
    
    if (is.null(rowRanges(scExp))) 
        stop("ChromSCape::feature_annotation_scExp - The object doesn't have ranges of coordinates as rowData")
    if (is.data.frame(reference_annotation)) 
        reference_annotation = makeGRangesFromDataFrame(reference_annotation, keep.extra.columns = T)
    
    feature_ranges = rowRanges(scExp)
    hits = GenomicRanges::distanceToNearest(feature_ranges, reference_annotation, ignore.strand = T, select = "all")
    
    annotFeat = data.frame(chr = as.character(seqnames(feature_ranges[queryHits(hits)])), start = as.character(start(feature_ranges[queryHits(hits)])), 
        end = as.character(end(feature_ranges[queryHits(hits)])), Gene = as.character(reference_annotation@elementMetadata$Gene)[subjectHits(hits)], 
        distance = hits@elementMetadata$distance) %>% mutate(ID = paste(chr, start, end, sep = "_")) %>% select(ID, chr, start, 
        end, Gene, distance)
    
    annotFeat <- annotFeat %>% group_by(ID) %>% summarise_all(funs(paste(unique(.), collapse = ", "))) %>% as.data.frame()
    rowData(scExp) = annotFeat
    return(scExp)
}

choose_perplexity <- function(dataset) {
    stopifnot(!is.null(dataset), !is.null(dim(dataset)))
    perplexity = 30
    if (nrow(dataset) <= 200) {
        perplexity = 20
    }
    if (nrow(dataset) <= 250) {
        perplexity = 25
    }
    if (nrow(dataset) <= 150) {
        perplexity = 15
    }
    if (nrow(dataset) <= 100) {
        perplexity = 10
    }
    if (nrow(dataset) <= 50) {
        perplexity = 5
    }
    perplexity
}

reduce_dims_scExp <- function(scExp, dimension_reductions = c("PCA", "TSNE"), n = 50, verbose = T) {
    stopifnot(is(scExp, "SingleCellExperiment"), is.numeric(n), dimension_reductions[1] %in% c("PCA", "TSNE"))
    
    if (!"normcounts" %in% names(assays(scExp))) {
        warning("ChromSCape::reduce_dims_scExp - The raw counts are not normalized, running dimensionality reduction on raw counts.")
        mat = counts(scExp)
    } else {
        mat = normcounts(scExp)
    }
    
    pca <- stats::prcomp(Matrix::t(mat), center = T, scale. = F)
    pca = pca$x[, 1:n]
    
    # Reduce the perplexity if the number of samples is too low to avoid perplexity error
    if ("TSNE" %in% dimension_reductions) 
        tsne <- Rtsne(pca, dims = 2, pca = FALSE, theta = 0, perplexity = choose_perplexity(pca), verbose = verbose, max_iter = 1000)
    tsne = as.data.frame(tsne$Y)
    
    # save PCA & T-SNE in scExp object
    if ("TSNE" %in% dimension_reductions) 
        reducedDims(scExp) <- list(PCA = as.data.frame(pca), TSNE = tsne) else reducedDims(scExp) <- list(PCA = pca)
    return(scExp)
}

num_cell_after_QC_filt_scExp <- function(scExp, annot){
  
  stopifnot(is(scExp, "SingleCellExperiment"), !is.null(datamatrix))
  
  table <- as.data.frame(table(annot$sample_id))
  table_filtered <- as.data.frame(table(colData(scExp)$sample_id))
  
  colnames(table) = c("Sample","#Cells Before Filtering")
  rownames(table) = NULL 
  colnames(table_filtered) = c("Sample","#Cells After Filtering")
  rownames(table_filtered) = NULL 
  
  table_both = left_join(table,table_filtered, by=c("Sample"))
  table_both[,1] = as.character(table_both[,1])
  table_both = table_both %>% 
    bind_rows(., tibble(Sample="",`#Cells Before Filtering`=sum(table_both[,2]),`#Cells After Filtering`=sum(table_both[,3]) ) )
  table_both %>% kable(escape=F, align="c") %>% kable_styling(c("striped", "condensed"), full_width = T) %>% group_rows("Total cell count", dim(table_both)[1], dim(table_both)[1])
}