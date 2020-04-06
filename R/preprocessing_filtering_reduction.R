## Authors : Pacôme Prompsy
## Title : Wrappers & functions to preprocess & reduce dimensionality of sc Matrix
## Description : Funtions to load, filter, normalize & reduce sc Epigenetic
## Matrices prior to analysis

# Wrapper to create the single cell experiment from sparse datamatrix and annot, remove Zero 
# count Features then Cells and calculate QC Metrics (scater::calculateQCMetrics)
#' Title
#'
#' @param datamatrix 
#' @param annot 
#' @param removeZeroCells 
#' @param removeZeroFeatures 
#'
#' @return
#' @export
#'
#' @examples
create_scExp <- function(datamatrix, annot, removeZeroCells = TRUE, removeZeroFeatures=TRUE){
  
  stopifnot(is.data.frame(annot),removeZeroCells %in% c(T,F), removeZeroFeatures %in% c(T,F))
  
  if(ncol(datamatrix) != nrow(annot))  stop('datamatrix and annot should contain the same number of cells')
  if(length(match(c("cell_id","sample_id"),colnames(annot))) < 2)  stop('annot should contain cell_id & sample_id as column names')
  
  scExp <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = datamatrix), colData = annot)
  
  dim_b = dim(scExp)
  if(removeZeroFeatures) scExp <- scExp[(Matrix::rowSums(counts(scExp)>0)>0),] # remove windows that do not have any read in any cells
  if(removeZeroCells) scExp <- scExp[,(Matrix::colSums(counts(scExp)>0)>0) ] # remove cells that do not have any read in any cells
  
  if( dim(scExp)[2] != dim_b[2] ) cat( dim_b[2] - dim(scExp)[2], "cells with 0 signals were removed.")
  if( dim(scExp)[1] != dim_b[1] ) cat( dim_b[1] - dim(scExp)[1], "features with 0 signals were removed.")
  
  scExp <- scater::calculateQCMetrics(scExp)
  
  if(has_genomic_coordinates(scExp)) {
    rows = rownames(scExp) 
    rowRanges(scExp) = get_genomic_coordinates(scExp)
    rownames(scExp) = rows
  }
  
  return(scExp)
}

# Function to filter out cells & features from sparse matrix based on total count per cell,
# number of cells 'ON' (count >= 2) in features and top covered cells that might be doublets
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
filter_scExp <- function(scExp, min_cov_cell = 1600, quant_removal = 95,
                         percentMin = 1, bin_min_count = 2, verbose = T){
  
  stopifnot(is(scExp, "SingleCellExperiment"), is.numeric(min_cov_cell),
            is.numeric(quant_removal), is.numeric(percentMin), is.numeric(bin_min_count), verbose %in% c(F,T))
  
  if(is.null(scExp)) warn("Please specify a SingleCellExperiment")
  
  cellCounts = Matrix::colSums(counts(scExp))
                               
  thresh <- quantile(cellCounts, probs=seq(0, 1, 0.01))
  
  sel1000 =  ( cellCounts > 1000 & cellCounts <= thresh[quant_removal+1] ) 
  sel <- ( cellCounts > min_cov_cell & cellCounts <= thresh[quant_removal+1] )
  if(verbose) cat(length(which(sel)),"cells pass the threshold of", min_cov_cell,
                  "minimum reads and are lower than the ", quant_removal, "th centile of library size ~=",round(thresh[quant_removal+1]),"reads.\n")
    
  # Create the binary matrix of cells >= 1000 counts
  # if count >= bin_min_count, count = 1 else
  bina_counts <- counts(scExp)[sel1000,]
  sel_below_2  <- (bina_counts < bin_min_count)
  bina_counts[sel_below_2] <- 0
  sel_above_1 = (bina_counts >= bin_min_count)
  bina_counts[sel_above_1] <- 1
  
  nCells_in_feature = Matrix::rowSums(bina_counts) 
  fixedFeature <- ( nCells_in_feature > ((percentMin/100)*(nrow(bina_counts))) ) # Feature selection
  if(verbose) cat(length(which(fixedFeature)),"features pass the threshold of", percentMin,
                  " % of total cells 'ON', representing a minimum of", round((percentMin/100)*(nrow(bina_counts))), "cells.\n")
  
  scExp <- scExp[,sel]
  scExp <- scExp[fixedFeature,]
  scExp <- scater::calculateQCMetrics(scExp) 
  return(scExp)
}

has_genomic_coordinates <- function(scExp){
  stopifnot(is(scExp, "SingleCellExperiment"), !is.null(rownames(scExp)))
  ID <- rownames(scExp)[1:min(10,length(rownames(scExp)))]
  chr <- unlist(lapply(strsplit(ID, split= "_|-|:"), FUN=function(x) unlist(x)[1]))

  if(length(grep("chr|(^\\d+$|^X$|^Y$)",chr[1:min(10,length(chr))],
                 ignore.case = T)) >= min(10,length(chr))) return(T) else return(F)
}

get_genomic_coordinates <- function(scExp){
  
  stopifnot(is(scExp, "SingleCellExperiment"))
  if(!has_genomic_coordinates(scExp)) stop("Feature names are not genomic coordinates")
  
  ID <- rownames(scExp)
  chr <- unlist(lapply(strsplit(ID, split= "_|-|:"), FUN=function(x) unlist(x)[1]))
  start <- unlist(lapply(strsplit(ID, split= "_|-|:"), FUN=function(x) unlist(x)[2]))
  end <- unlist(lapply(strsplit(ID, split= "_|-|:"), FUN=function(x) unlist(x)[3]))
  feature <- GRanges(data.frame(ID=ID, chr=chr, start=start, end=end))
  return(feature)
}

exclude_features <- function(scExp, features_to_exclude, by = c("region","feature_name"), verbose =T) {
  
  stopifnot(is(scExp,"SingleCellExperiment"), is.data.frame(features_to_exclude),
            is.character(by[1]))
  if(!by[1] %in% c("region","feature_name")) stop("by must be either 'region' or 'feature_name'")
  
  if(by[1] == "region"){
    if(!has_genomic_coordinates(scExp)) stop("Feature names are not genomic coordinates")
    regions <- rowRanges(scExp)
    colnames(features_to_exclude)[1:3] = c("chr","start","stop")
    excl_gr <- makeGRangesFromDataFrame(features_to_exclude, ignore.strand=TRUE,
                  seqnames.field=c("chr"), start.field=c("start"), end.field=c("stop"))
    ovrlps <- as.data.frame(findOverlaps(regions, excl_gr))[, 1]
    if(length(unique(ovrlps)>0)) scExp <- scExp[-unique(ovrlps), ]
    if(verbose) cat("Removed",length(unique(ovrlps)),"regions from the analysis.\n")
  }
  if(by[1] == "feature_name"){
    if(has_genomic_coordinates(scExp)) warning("Excluding by feature name while object feature names are genomic coordinates !") 
    features = rownames(scExp)
    features_to_exclude = as.character(features_to_exclude[,1])
    ovrlps = intersect(features,features_to_exclude)
    if(length(unique(ovrlps)>0)) scExp <- scExp[-which(rownames(scExp) %in% ovrlps), ]
    if(verbose) cat("Removed",length(unique(ovrlps))," features from the analysis.\n")
  }
  return(scExp)
}


preprocess_TPM <- function(scExp){
  
  size = BiocGenerics::width(SummarizedExperiment::rowRanges(scExp))
  assay(scExp) = assay(scExp) / size
  assay(scExp) = 10^6 * t( t(assay(scExp)) / Matrix::colSums(assay(scExp)) ) 
  
  return(scExp)
}

preprocess_RPKM <- function(scExp){
  
  assay(scExp) = 10^9 * t( t(assay(scExp)) / Matrix::colSums(assay(scExp)) ) 
  size = BiocGenerics::width(SummarizedExperiment::rowRanges(scExp))
  assay(scExp) = assay(scExp) / size
  
  return(scExp)
}

preprocess_CPM <- function(scExp){
  assay(scExp) = 10^6 * t( t(assay(scExp)) / Matrix::colSums(assay(scExp))) 
  return(scExp)
}

preprocess_feature_size_only <- function(scExp){
  size = BiocGenerics::width(SummarizedExperiment::rowRanges(scExp))
  assay(scExp) = assay(scExp) / size
  return(scExp)
}

normalize_scExp <- function(scExp, type = c("RPKM","CPM","TPM","feature_size_only")){
  stopifnot(type[1] %in% c("RPKM","CPM","TPM","feature_size_only"),
            is(scExp,"SingleCellExperiment"))
  if(!has_genomic_coordinates(scExp)) {
    warning("Switching to CPM normalization as features are not genomic coordinates.")
    type = "CPM"
  }
  switch (type,
          "TPM" = return(preprocess_TPM(scExp)),
          "RPKM" = return(preprocess_RPKM(scExp)),
          "feature_size_only" = return(preprocess_feature_size_only(scExp)),
          "CPM" = return(preprocess_CPM(scExp))
          )
}

feature_annotation_scExp <- function(scExp, reference_annotation){
  stopifnot(is(scExp,"SingleCellExperiment"))
  
  if(is.null(rowRanges(scExp)) ) stop("The object doesn't have ranges of coordinates as rowData") 
  if(is.data.frame(reference_annotation)) reference_annotation = makeGRangesFromDataFrame(reference_annotation,keep.extra.columns = T) 
  
  feature_ranges = rowRanges(scExp)
  hits = GenomicRanges::distanceToNearest(feature_ranges,reference_annotation,ignore.strand=T,select="all")
  
  annotFeat = data.frame(chr = as.character(seqnames(feature_ranges[queryHits(hits)])),
                         start = as.character(start(feature_ranges[queryHits(hits)])),
                         end = as.character(end(feature_ranges[queryHits(hits)])),
                         Gene =  as.character(reference_annotation@elementMetadata$Gene)[subjectHits(hits)],
                         distance = hits@elementMetadata$distance
  ) %>% mutate(ID = paste(chr,start,end,sep="_")) %>% select(ID,chr,start,end,Gene,distance)
  
  annotFeat <- annotFeat %>% group_by(ID) %>% summarise_all(funs(paste(unique(.), collapse = ', '))) %>% as.data.frame()
  IDs = annotFeat$ID
  rowData(scExp) = makeGRangesFromDataFrame(annotFeat,keep.extra.columns = T)
  return(scExp)
}
  
