# This script is for testing reproducibility between
# the branch "package" & the branch "master" of ChromSCape
# Step 1: filtering & reduction

context("Testing reproducibility of preprocessing & filt with master")

out = create_scDataset_raw(featureType = "window", sparse = F)
datamatrix = out$mat
annot_raw = out$annot

master = new.env()

test_that("Re checking that seed is the same", {
  load("tests/test_scChIP/scChIP_raw.RData", master)
  expect_equal(master$datamatrix, datamatrix )
  expect_equal(master$annot_raw, annot_raw )
})

scExp = create_scExp(datamatrix,annot_raw)

test_that("Step 1 : creating scExp", {
  load("tests/test_scChIP/filter_red_01.RData", master)
  expect_equal(assay(master$umi), assay(scExp))
  expect_equal(colData(master$umi), colData(scExp))
})

scExp = filter_scExp(scExp)

test_that("Step 2 : filtering ", {
  load("tests/test_scChIP/filter_red_02.RData", master)
  expect_equal(master$SelMatCov, assay(scExp))
})

regions_to_exclude = read.table("../../Data/Annotation/bed/MM468_5FU3_all_5FU5_initial_CNV_from_ChIP_input.bed")
scExp = exclude_features_scExp(scExp,features_to_exclude = regions_to_exclude, by ="region")

test_that("Step 3 : remove specific features ", {
  load("tests/test_scChIP/filter_red_03.RData", master)
  expect_equal(master$SelMatCov, assay(scExp))
})

scExp = normalize_scExp(scExp, "CPM")
test_that("Step 4 : Normalize ", {
  load("tests/test_scChIP/filter_red_04.RData", master)
  expect_equivalent(master$norm_mat, normcounts(scExp))
})

# reference_annotation = read.table("annotation/hg38/Gencode_TSS_pc_lincRNA_antisense.bed",col.names = c("chr","start","end","Gene"))
# scExp = feature_annotation_scExp(scExp,reference_annotation)
# 
# test_that("Step 5 : feature annotation ", {
#   load("tests/test_scChIP/Simulated_window_300_600_not_sparse_seed47_1600_1_95_uncorrected_annotFeat.RData",master)
#   to_test = as.numeric(as.data.frame(rowData(scExp))$distance)
#   original = as.numeric(makeGRangesFromDataFrame(master$annotFeat,
#                                                  keep.extra.columns = T)$distance)
#   original[original > 0] = original[original > 0] -2
#   expect_equal(to_test,original)
# })

# test_that("Step 6 : batch correction ", {
#   scExp = ChromSCape::filter_scExp(scExp)
#   load("tests/test_scChIP/filter_red_02.RData", master)
#   expect_equal(master$SelMatCov, assay(scExp))
# })

scExp = reduce_dims_scExp(scExp)

test_that("Step 7 : dimensionality reduction", {
 
  load("tests/test_scChIP/Simulated_window_300_600_not_sparse_seed47_1600_1_95_uncorrected.RData", master)
  
  expect_equal(as.data.frame(master$pca),reducedDim(scExp,"PCA"))
})

scExp = correlation_and_hierarchical_clust_scExp(scExp)

test_that("Step 8 : correlation & hiearchical clust", {
  
  load("tests/test_scChIP/Simulated_window_300_600_not_sparse_seed47_1600_1_95_uncorrected_99_1_cor_filtered.RData", master)
  expect_equal(master$hc_cor$merge,scExp@metadata$hc_cor$merge)
  expect_equal(master$hc_cor$height,scExp@metadata$hc_cor$height)
  expect_equal(master$hc_cor$order,scExp@metadata$hc_cor$order)
  expect_equal(master$hc_cor$labels,scExp@metadata$hc_cor$labels)
  expect_equal(master$hc_cor$dist.method,scExp@metadata$hc_cor$dist.method)
  expect_equal(master$hc_cor$call$method,scExp@metadata$hc_cor$call$method)
  expect_equal(cor(t(master$pca)),reducedDim(scExp,"Cor") )
  expect_equal(colnames(cor(t(master$pca))),colnames(reducedDim(scExp,"Cor")) )
  expect_equal(rownames(cor(t(master$pca))),rownames(reducedDim(scExp,"Cor")) )
})

scExp. = filter_correlated_cell_scExp(scExp,random_iter = 500)

test_that("Step 9 : correlation filterin clust", {
  
  load("tests/test_scChIP/Simulated_window_300_600_not_sparse_seed47_1600_1_95_uncorrected_99_1_cor_filtered.RData", master)
  load("tests/test_scChIP/z.RData", master)
  expect_equal(nrow(master$annot_sel),ncol(scExp.))
  expect_equal(master$annot_sel$cell_id,colData(scExp.)$cell_id)
  table(master$annot_sel$sample_id)
  table(colData(scExp.)$sample_id)
  
})

