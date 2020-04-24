# This script is for testing reproducibility between
# the branch "package" & the branch "master" of ChromSCape
# Step 1: filtering & reduction

context("Testing reproducibility of preprocessing & filt with master")

library(testthat)
library(ChromSCape)
create_scDataset_raw <- function(cells=300,features=600,
                                 featureType = c("window","peak","gene"),
                                 sparse=T, nsamp=4, ref = "hg38",batch_id = rep(1,nsamp)) {
  
  stopifnot(featureType %in% c("window","peak","gene"), ref %in% c("mm10","hg38"),
            nsamp >= 1, cells >= nsamp, features >=1, length(batch_id) == nsamp)
  
  stopifnot()
  
  set.seed(47)
  
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
    chr_ranges = unlist(GenomicRanges::tileGenome(setNames(width(chr),GenomicRanges::seqnames(chr)),
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
      chr_ranges = unlist(GenomicRanges::tileGenome(setNames(GenomicRanges::width(chr),GenomicRanges::seqnames(chr)),
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
  for(i in 1:10) vec[vec >= i] = vec[vec >= i]  +  i^2*stats::rpois(length(vec[vec >= i]),0.5)
  mat = matrix(vec, nrow = features, ncol = cells, 
               dimnames = list( features_names,cell_names))
  annot = data.frame(cell_id = cell_names,
                     sample_id = sample,
                     batch_id = batches,
                     total_counts = Matrix::colSums(mat))
  if(sparse) return(list("mat" =  as(mat,"dgCMatrix"), "annot" = annot)) else return(list("mat" =  mat, "annot" = annot))
}

out = create_scDataset_raw(featureType = "window", sparse = T)
datamatrix = out$mat
annot_raw = out$annot

master = new.env()

test_that("Re checking that seed is the same", {
  load("tests/test_scChIP/scChIP_raw.RData", master)
  expect_equal(master$datamatrix, as.matrix(datamatrix) )
  expect_equal(master$annot_raw, annot_raw )
})

scExp = create_scExp(datamatrix,annot_raw)

test_that("Step 1 : creating scExp", {
  load("tests/test_scChIP/filter_red_01.RData", master)
  expect_equal(SummarizedExperiment::assay(master$umi), as.matrix(SummarizedExperiment::assay(scExp)) )
  expect_equal(colData(master$umi), colData(scExp))
})

scExp = filter_scExp(scExp)

test_that("Step 2 : filtering ", {
  load("tests/test_scChIP/filter_red_02.RData", master)
  expect_equal(master$SelMatCov, as.matrix(SummarizedExperiment::assay(scExp)))
})

regions_to_exclude = read.table("../../Data/Annotation/bed/MM468_5FU3_all_5FU5_initial_CNV_from_ChIP_input.bed")
scExp = exclude_features_scExp(scExp,features_to_exclude = regions_to_exclude, by ="region")

test_that("Step 3 : remove specific features ", {
  load("tests/test_scChIP/filter_red_03.RData", master)
  expect_equal(master$SelMatCov, as.matrix(SummarizedExperiment::assay(scExp)))
})

scExp = normalize_scExp(scExp, "CPM")
test_that("Step 4 : Normalize ", {
  load("tests/test_scChIP/filter_red_04.RData", master)
  expect_equivalent(master$norm_mat, as.matrix(SingleCellExperiment::normcounts(scExp)) )
})

scExp = feature_annotation_scExp(scExp,ref="hg38")

test_that("Step 5 : feature annotation ", {
  load("tests/test_scChIP/Simulated_window_300_600_not_sparse_seed47_1600_1_95_uncorrected_annotFeat.RData",master)
  to_test = as.numeric(as.data.frame(SummarizedExperiment::rowData(scExp))$distance)
  original = as.numeric(GenomicRanges::makeGRangesFromDataFrame(master$annotFeat,
                                                 keep.extra.columns = T)$distance)
  original[original > 0] = original[original > 0] -2 # Difference between bedtools & GenomicRanges of 2 (start = 1 vs start = 0 ?)
  expect_equal(to_test,original[match(rownames(scExp),master$annotFeat$ID)]) #re order !
  
  expect_equal(master$annotFeat$Gene[match(rownames(scExp),master$annotFeat$ID)],rowData(scExp)$Gene)
  
})

# test_that("Step 6 : batch correction ", {
#   scExp = ChromSCape::filter_scExp(scExp)
#   load("tests/test_scChIP/filter_red_02.RData", master)
#   expect_equal(master$SelMatCov, assay(scExp))
# })

scExp = reduce_dims_scExp(scExp)

test_that("Step 7 : dimensionality reduction", {
  
  load("tests/test_scChIP/Simulated_window_300_600_not_sparse_seed47_1600_1_95_uncorrected.RData", master)
  #expect_equal(as.data.frame(master$pca), reducedDim(scExp,"PCA")) # at a sign
  expect_equal(cor(master$pca[,1], reducedDim(scExp,"PCA")[,1]), -1)
  expect_equal(cor(master$pca[,2], reducedDim(scExp,"PCA")[,2]), 1)
  expect_equal(cor(master$pca[,3], reducedDim(scExp,"PCA")[,3]), 1)

})

# BECAUSE OF PCA signs are random, all the downstream steps are diverging. Therefore we put
# PCA from master into reducedDim(scExp,"PCA") to keep signs to test all the
# downstream steps
SingleCellExperiment::reducedDim(scExp, "PCA") = as.data.frame(master$pca)
scExp = correlation_and_hierarchical_clust_scExp(scExp)

test_that("Step 8 : correlation & hiearchical clust", {
  
  load("tests/test_scChIP/Simulated_window_300_600_not_sparse_seed47_1600_1_95_uncorrected_99_1_cor_filtered.RData", master)
  expect_equal(master$hc_cor$merge,scExp@metadata$hc_cor$merge) # correlation not same
  expect_equal(master$hc_cor$height,scExp@metadata$hc_cor$height)
  expect_equal(master$hc_cor$order,scExp@metadata$hc_cor$order)
  expect_equal(cor(t(master$pca)),reducedDim(scExp,"Cor") ) # Depending on the sign of the PCA, pearson correlation matrix will give different correlation results
  expect_equal(colnames(cor(t(master$pca))),colnames(reducedDim(scExp,"Cor")) )
  expect_equal(rownames(cor(t(master$pca))),rownames(reducedDim(scExp,"Cor")) )

})

scExp_cf = filter_correlated_cell_scExp(scExp,random_iter = 50)

test_that("Step 9 : correlation filterin clust", {
  
  load("tests/test_scChIP/Simulated_window_300_600_not_sparse_seed47_1600_1_95_uncorrected_99_1_cor_filtered.RData", master)
  load("tests/test_scChIP/z.RData", master)
  expect_equal(nrow(master$annot_sel),ncol(scExp_cf))
  expect_equal(master$annot_sel$cell_id,colData(scExp_cf)$cell_id)
  expect_equal(table(master$annot_sel$sample_id),table(colData(scExp_cf)$sample_id))
})


scExp_cf = consensus_clustering_scExp(scExp_cf,prefix = "",reps = 1000, seed = 3.14)

test_that("Step 10 : correlation consensus hierarchical clustering", {
  
  load("tests/test_scChIP/Simulated_window_300_600_not_sparse_seed47_1600_1_95_uncorrected_99_1_consclust.RData", master)
  
  # consclust list
  expect_equal(master$consclust[[2]]$consensusMatrix, scExp_cf@metadata$consclust[[2]]$consensusMatrix)
  expect_equal(master$consclust[[5]]$consensusMatrix, scExp_cf@metadata$consclust[[5]]$consensusMatrix)
  expect_equal(master$consclust[[8]]$consensusClass, scExp_cf@metadata$consclust[[8]]$consensusClass)
  
  #icl list
  expect_equal(master$icl$clusterConsensus, scExp_cf@metadata$icl$clusterConsensus)
  expect_equal(master$icl$itemConsensus, scExp_cf@metadata$icl$itemConsensus)
})

scExp_cf = choose_cluster_scExp(scExp_cf, nclust = 2)

test_that("Step 11 : correlation consensus hierarchical clustering - choose cluster", {
  
  load("tests/test_scChIP/Simulated_window_300_600_not_sparse_seed47_1600_1_95_uncorrected_99_1_affectation_k2.RData", master)
  
  # cell affectation to clusters
  affectation. = as.data.frame(colData(scExp_cf))[,c("cell_id","sample_id","chromatin_group")]
  colnames(affectation.)[3] = "ChromatinGroup"
  all_equal(master$affectation,  affectation. )
  all_equal(table(master$affectation[,c("sample_id","ChromatinGroup")]),  table(affectation.[,c("sample_id","ChromatinGroup")]) )
  
  # tsne
  load("tests/test_scChIP/Simulated_window_300_600_not_sparse_seed47_1600_1_95_uncorrected_99_1_tsne_filtered.RData", master)
  all_equal(master$tsne_filtered$Y, reducedDim(scExp_cf,"TSNE") )
  
})

scExp_cf = differential_analysis_scExp(scExp_cf, qval.th = 0.4, cdiff.th = 0.3)

test_that("Step 12 : differential analysis between clusters", {
  
  load("tests/test_scChIP/Simulated_window_300_600_not_sparse_seed47_1600_1_95_uncorrected_99_1_2_0.4_0.3_one_vs_rest.RData", master)
  
  my.res_original = as.data.frame(apply(master$my.res_save, MARGIN = 2, as.character), stringsAsFactors=F)
  my.res_new = as.data.frame(apply(scExp_cf@metadata$diff$res, MARGIN = 2, as.character), stringsAsFactors=F)
  all_equal(my.res_original, my.res_new)
  
  all_equal(master$summary_save, scExp_cf@metadata$diff$summary)
  all_equal(my.res_original, my.res_new)
  
})

scExp_cf = gene_set_enrichment_analysis_scExp(scExp_cf,ref = "hg38", qval.th = 0.4, cdiff.th = 0.3,
                                              use_peaks = F )

test_that("Step 12 : GSEA of genes associated to differential loci:", {
  
  load("tests/test_scChIP/Simulated_window_300_600_not_sparse_seed47_1600_1_95_uncorrected_99_1_2_0.4_0.3_one_vs_rest_GSEA.RData", master)
  
  annotFeat_long_original = as.data.frame(cSplit(master$annotFeat[match(rownames(scExp_cf), 
                                                                        master$annotFeat$ID),], splitCols="Gene", sep=", ", direction="long"))
  annotFeat_long_new = as.data.frame(cSplit(as.data.frame(rowData(scExp_cf)), splitCols="Gene", sep=", ", direction="long"))
  
  all_equal(annotFeat_long_original$chr,annotFeat_long_new$chr)
  all_equal(annotFeat_long_original$start,annotFeat_long_new$start)
  all_equal(annotFeat_long_original$end,annotFeat_long_new$end)
  all_equal(annotFeat_long_original$Gene,annotFeat_long_new$Gene)
  
  annotFeat_long_original$distance[as.numeric(annotFeat_long_original$distance)>0] = as.character(as.numeric(
    annotFeat_long_original$distance[as.numeric(annotFeat_long_original$distance)>0]) -2)
  all_equal(annotFeat_long_original$distance,annotFeat_long_new$distance) # diff of -2 for >0 regions (#annotFeat)
  
  # GSEA
  all_equal(master$Both,scExp_cf@metadata$enr$Both)
  all_equal(master$Overexpressed[[1]],scExp_cf@metadata$enr$Overexpressed[[1]])
  all_equal(master$Overexpressed[[2]],scExp_cf@metadata$enr$Overexpressed[[2]])
  all_equal(master$Underexpressed[[1]],scExp_cf@metadata$enr$Underexpressed[[1]])
  all_equal(master$Underexpressed[[2]],scExp_cf@metadata$enr$Underexpressed[[2]])

})
