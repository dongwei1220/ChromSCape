context("Testing preprocessing, filtering & reduction functions")



# Functions for testing purposes
create_scDataset_raw <- function(cells=300,features=600,
                                 featureType = c("window","peak","gene"),
                                sparse=T, nsamp=4, ref = "hg38",batch_id = rep(1,nsamp)) {
  
  stopifnot(featureType %in% c("window","peak","gene"), ref %in% c("mm10","hg38"),
            nsamp >= 1, cells >= nsamp, features >=1, length(batch_id) == nsamp)
  
  stopifnot()
  
  set.seed(147)
  
  # Create cell names
  cell_counts = sapply(split( 1:cells , sample(nsamp, cells , repl = TRUE) ), length)
  cell_names = sample = batches = list()
  for(i in 1:length(cell_counts)) {
    cell_names[[i]] = paste0("sample_",i,"_c", 1:cell_counts[i])
    sample[[i]] = rep(paste0("sample_",i),cell_counts[i])
    batches[[i]] = rep(batch_id[i],cell_counts[i])
  }
  cell_names = as.character(unlist(cell_names))
  sample = as.character(unlist(sample))
  batches = as.numeric(unlist(batches))
  
  # Create feature names
  eval(parse(text = paste0("chr <- ChromSCape::",ref,".chromosomes"))) #load species chromsizes
  chr = GRanges(chr)
  
  if(featureType[1] == "window") {
    chr_ranges = unlist(tileGenome(setNames(width(chr),seqnames(chr)),
                                      ntile = features))[1:features] # ~constant window size
    features_names = paste(as.data.frame(chr_ranges)$seqnames,
                           as.data.frame(chr_ranges)$start,
                           as.data.frame(chr_ranges)$end, sep="_")
  }
  if(featureType[1] == "peak") {
    size_peaks = c(1000,2500,7999,10000,150000,10^6) #Different size of peaks
    peaks = sapply(split( 1:features , sample(length(size_peaks), features , repl = TRUE) ), length)
    chr_ranges_list = GRangesList()
    for(i in 1:length(peaks)){
      chr_ranges = unlist(tileGenome(setNames(width(chr),seqnames(chr)),
                                     tilewidth = size_peaks[i], cut.last.tile.in.chrom = F))
      chr_ranges_list[[i]] = chr_ranges[sample(1:length(chr_ranges),size = peaks[i]),]
    }
    chr_ranges = GenomicRanges::sort.GenomicRanges(unlist(chr_ranges_list))[1:features]
    
    features_names = paste(as.data.frame(chr_ranges)$seqnames,
                           as.data.frame(chr_ranges)$start,
                           as.data.frame(chr_ranges)$end, sep="_")
  }
  if(featureType[1] == "gene") {
    eval(parse(text = paste0("chr <- ChromSCape::",ref,".GeneTSS"))) #load species chromsizes
    features_names = as.character(sample(chr$gene,features,replace = F))
  }
  vec = rpois(cells*features,0.5) #Add count to values > 0, iteratively
  for(i in 1:10) vec[vec >= i] = vec[vec >= i]  +  i^2*rpois(length(vec[vec >= i]),0.5)
  mat = matrix(vec, nrow = cells, ncol = features, 
               dimnames = list(cell_names, features_names))
  annot = data.frame(cell_id = cell_names,
                     sample_id = sample,
                     batch_id = batches,
                     total_counts = Matrix::rowSums(mat))
  return(list("mat" =  mat, "annot" = annot))
}

out = create_scDataset_raw(featureType = "window")
mat = out$mat
annot = out$annot

#### create_scExp
# Wrapper to create the single cell experiment from sparse datamatrix and annot, remove Zero 
# count Features then Cells and calculate QC Metrics (scater::calculateQCMetrics)

test_that("Wrong input - basic", {
  expect_error(create_scExp(c(),list()))
  expect_error(create_scExp(-1, ))
  expect_error(create_scExp(NA, NA))
})

test_that("Wrong input - advanced", {
  expect_error(create_scExp(annot, mat))
  expect_error(create_scExp(mat, annot, removeZeroCells = 3))
  expect_error(create_scExp(mat, annot, removeZeroFeatures = 3))
  expect_error(create_scExp(mat, annot[1:5,]))
})

test_that("Some cells are empty", {
  mat. = mat
  mat.[sample(1:nrow(mat.),3),] = 0
  annot. = annot
  annot.$total_counts = Matrix::rowSums(mat.)
  expect_output(create_scExp(mat.,annot.), "cells with 0 signals were removed." )
})

test_that("Some features are empty", {
  mat. = mat
  mat.[,sample(1:ncol(mat.),3)] = 0
  annot. = annot
  annot.$total_counts = Matrix::rowSums(mat.)
  expect_output(create_scExp(mat.,annot.), "features with 0 signals were removed." )
})

#### filter_scExp
# Function to filter out cells & features from sparse matrix based on total count per cell,
# number of cells 'ON' (count >= 2) in features and top covered cells that might be doublets

scExp = create_scExp(mat,annot)

test_that("Wrong input - basic", {
  expect_error(filter_scExp(c(),list()))
  expect_error(filter_scExp(-1, ))
  expect_error(filter_scExp(NA, NA))
})

test_that("No cell filter doesn't change number cells", {
  
  expect_equal(nrow(filter_scExp(scExp,
                            percentMin = 0,
                            quant_removal = 100,
                            min_cov_cell = 0)), nrow(scExp) )
})

test_that("No feature filter doesn't change number features", {
  
  expect_equal(ncol(filter_scExp(scExp,
                                 percentMin = 0,
                                 quant_removal = 100,
                                 min_cov_cell = 0)), ncol(scExp) )
})

test_that("Max cell filters remove all cells", {
  
  expect_equal(nrow(filter_scExp(scExp,
                    min_cov_cell = max(Matrix::rowSums(counts(scExp))) )),0 )
})

test_that("Max feature filters remove all features", {
  expect_equal(ncol(filter_scExp(scExp,percentMin = 101)),0 )
})


test_that("Verbose is on /off", {
  expect_output(filter_scExp(scExp,verbose=T), "cells pass the threshold of" )
  expect_output(filter_scExp(scExp,verbose=T), "features pass the threshold of" )
  expect_invisible({scExp. =filter_scExp(scExp,verbose=F)} )
})

test_that("Some cells are empty", {
  mat. = mat
  mat.[sample(1:nrow(mat.),3),] = 0
  annot. = annot
  annot.$total_counts = Matrix::rowSums(mat.)
  scExp. = create_scExp(mat., annot.)
  expect_type(filter_scExp(scExp.),typeof(scExp))
})

test_that("Some features are empty", {
  mat. = mat
  mat.[,sample(1:ncol(mat.),3)] = 0
  annot. = annot
  annot.$total_counts = Matrix::rowSums(mat.)
  scExp. = create_scExp(mat., annot.)
  expect_type(filter_scExp(scExp.),typeof(scExp))
})



