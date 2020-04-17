# This script is for testing reproducibility between
# the branch "package" & the branch "master" of ChromSCape
# Step 1: filtering & reduction

context("Testing reproducibility of preprocessing & filt with master")

master = new.env()
load("../../Data/ChromSCape_Data/Reproducibility_new/datasets/HBCx_95_mm10_001/scChIP_raw.RData",master)
datamatrix = master$datamatrix
annot_raw = master$annot_raw

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

scExp = normalize_scExp(scExp, "RPKM")
test_that("Step 4 : Normalize ", {
  load("tests/test_scChIP/filter_red_04.RData", master)
  expect_equal(master$norm_mat, assay(scExp))
})

reference_annotation = read.table("annotation/hg38/Gencode_TSS_pc_lincRNA_antisense.bed",col.names = c("chr","start","end","Gene"))
scExp = feature_annotation_scExp(scExp,reference_annotation)

test_that("Step 5 : feature annotation ", {
  load("tests/test_scChIP/Simulated_window_300_600_not_sparse_seed47_1600_1_95_uncorrected_annotFeat.RData",master)
  
  to_test = as.numeric(as.data.frame(rowData(scExp)$X)$distance)
  original = as.numeric(makeGRangesFromDataFrame(master$annotFeat,
                                                 keep.extra.columns = T)$distance)
  original[original > 0] = original[original > 0] -2
  expect_equal(to_test,original)
})

test_that("Step 6 : batch correction ", {
  scExp = ChromSCape::filter_scExp(scExp)
  load("tests/test_scChIP/filter_red_02.RData", master)
  expect_equal(master$SelMatCov, assay(scExp))
})

test_that("Step 7 : dimensionality reduction", {
  scExp = ChromSCape::filter_scExp(scExp)
  load("tests/test_scChIP/filter_red_02.RData", master)
  expect_equal(master$SelMatCov, assay(scExp))
})