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
scExp = exclude_features(scExp,features_to_exclude = regions_to_exclude, by ="region")

test_that("Step 3 : remove specific features ", {
  load("tests/test_scChIP/filter_red_03.RData", master)
  expect_equal(master$SelMatCov, assay(scExp))
})

scExp = normalize_scExp(scExp, "RPKM")
test_that("Step 4 : Normalize ", {
  load("tests/test_scChIP/filter_red_04.RData", master)
  expect_equal(master$norm_mat, assay(scExp))
})

test_that("Step 5 : feature annotation ", {
  scExp = ChromSCape::filter_scExp(scExp)
  load("tests/test_scChIP/filter_red_02.RData", master)
  expect_equal(master$SelMatCov, assay(scExp))
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