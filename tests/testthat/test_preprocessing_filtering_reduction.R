context("Testing preprocessing, filtering & reduction functions")



# Functions for testing purposes
create_scDatamatrix <- function(cells=300,features=600,
                                featureType = c("window","peak","gene"),
                                sparse=T, nsamp=4, ref = "hg38") {
  
  set.seed(147)
  cell_counts = rep(cells, nsamp)
  cell_names = list(paste0("sample_",nsamp,"_c", 1:cells))
  total = cells
  if(nsamp > 1){
    for(i in 1:(nsamp-1)){
      cell_counts[i] = round(runif(1)*total) + 1
      total = total - cell_counts[i]
      cell_names[[i]] = paste0("sample_",i,"_c", 1:cell_counts[i])
    }
    cell_counts[nsamp] = total
    cell_names[[nsamp]] = paste0("sample_",nsamp,"_c", 1:cell_counts[nsamp])
  }
  cell_names = as.character(unlist(cell_names))
  
  if(featureType[1] == "window") invisible(paste0(ref,".chromosomes")) features_names = paste()
  if(featureType[1] == "peak") features_names = paste()
  if(featureType[1] == "gene") features_names = paste() 
  
  return(matrix(rpois(cells*features,0.5),
          nrow = cells, ncol = features, 
          dimnames = list(cell_names, feature_names)))
}
create_annot <- function(datamatrix)

create_scExp
filter_scExp

test_that("non matrix given as input", {
  expect_error(create_scExp(-1, ))
  expect_error(create_scExp(list(a = c(1,1,1), )))
  expect_error(create_scExp(NA), )
})

test_that("vectors", {
  expect_equal(increment(c(0,1)), c(1,2))
})

test_that("empty vector", {
  expect_equal(increment(c()), c())
})

test_that("test NA", {
  expect_true(is.na(increment(NA)))
})
# Function to filter out cells & features from sparse matrix based on total count per cell,
# number of cells 'ON' (count >= 2) in features and top covered cells that might be doublets
 
