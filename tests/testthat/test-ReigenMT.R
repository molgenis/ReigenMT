# load the h5 files
test_that("h5 loading works", {
  # location of the file
  h5_loc <- 'resources/interaction_qtl.h5'
  # read the file
  h5_loaded <- limix_h5_to_sumstats_format(h5_loc)
  # check if the number of columns and rows match
  expect_equal((ncol(h5_loaded) == 6 & nrow(h5_loaded) == 25), T)
})

# test eigenvalue decomposition works
test_that("eigenvector decomposition works", {
  # Generate 100 numbers that are either 0, 1, or 2, representing genotypes
  set.seed(1)
  number_variant1 <- sample(0:2, 100, replace = TRUE)
  set.seed(2)
  number_variant2 <- sample(0:2, 100, replace = TRUE)
  set.seed(3)
  number_variant3 <- sample(0:2, 100, replace = TRUE)
  set.seed(4)
  number_variant4 <- sample(0:2, 100, replace = TRUE)
  set.seed(5)
  number_variant5 <- sample(0:2, 100, replace = TRUE)
  # create datatable
  genotypes <- data.table::data.table(var1 = number_variant1, 
                                      var2 = number_variant2, 
                                      var3 = number_variant3, 
                                      var4 = number_variant4, 
                                      var5 = number_variant5, 
                                      var6 = number_variant1, 
                                      var7 = number_variant2, 
                                      var8 = number_variant2)
  # transpose, because that is what the function expects
  genotypes <- data.table::data.table(t(genotypes))
  genotype_eigenvalues <- genotypes_to_eigenvalues(genotypes)
  # the number of eigenvalues should be the number of variants we put in
  expect_equal((length(genotype_eigenvalues) == 8), T)
})

# test eigenvalue thresholding works
test_that("eigenvector thresholding works", {
  # Generate 100 numbers that are either 0, 1, or 2, representing genotypes
  set.seed(1)
  number_variant1 <- sample(0:2, 100, replace = TRUE)
  set.seed(2)
  number_variant2 <- sample(0:2, 100, replace = TRUE)
  set.seed(3)
  number_variant3 <- sample(0:2, 100, replace = TRUE)
  set.seed(4)
  number_variant4 <- sample(0:2, 100, replace = TRUE)
  set.seed(5)
  number_variant5 <- sample(0:2, 100, replace = TRUE)
  # create datatable
  genotypes <- data.table::data.table(var1 = number_variant1, 
                                      var2 = number_variant2, 
                                      var3 = number_variant3, 
                                      var4 = number_variant4, 
                                      var5 = number_variant5, 
                                      # the last three genotypes are the same as three previous ones
                                      var6 = number_variant1, 
                                      var7 = number_variant2, 
                                      var8 = number_variant2)
  # transpose, because that is what the function expects
  genotypes <- data.table::data.table(t(genotypes))
  genotype_eigenvalues <- genotypes_to_eigenvalues(genotypes)
  # the number of eigenvalues should be the number of variants we put in
  neigen <- find_number_of_eigen(genotype_eigenvalues, 8)
  # check if the n eigen is five
  expect_equal((neigen == 5), T)
})

