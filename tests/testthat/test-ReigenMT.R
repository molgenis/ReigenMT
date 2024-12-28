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

# get positions of variants in MatrixEQTL format
test_that("variant position extractions MatrixEQTL format works", {
  # numeric genotypes
  numeric_gens_loc <- 'resources/genotypes.txt.gz'
  numeric_gens <- data.table::fread(numeric_gens_loc, header = T)
  # summary stats
  #summary_stats_loc <- 'resources/cis.eqtls.txt'
  #summary_stats <- data.table::fread(summary_stats_loc, header = T)
  # variant positions
  var_positions_loc <- 'resources/gen.positions.txt.gz'
  var_positions <- data.table::fread(var_positions_loc, header = T)[1:100, ]
  # make into hashtable
  var_positions_hm <- get_position_hashtab(var_positions, 'snp', 'pos')
  # check if read the last one correctly one
  expect_equal((var_positions_hm[['snp_19_273384']] == 273384), T)
})

# test eigenMT 
test_that("eigenMT works on MatrixEQTL format", {
  # numeric genotypes
  numeric_gens_loc <- 'resources/genotypes.txt.gz'
  numeric_gens <- data.table::fread(numeric_gens_loc, header = T)
  # summary stats
  summary_stats_loc <- 'resources/cis.eqtls.txt.gz'
  summary_stats <- data.table::fread(summary_stats_loc, header = T)
  # get the unique genes
  summary_stats_genes <- unique(summary_stats[['gene']])
  # set a seed
  set.seed(1337)
  # select ten random genes
  summary_stats_genes_random100 <- summary_stats_genes[sample.int(length(summary_stats_genes), 100)]
  # just do a couple of genes
  summary_stats <- summary_stats[summary_stats[['gene']] %in% summary_stats_genes_random100, ]
  # variant positions
  var_positions_loc <- 'resources/gen.positions.txt.gz'
  var_positions <- data.table::fread(var_positions_loc, header = T)
  # and those couple of variants
  var_positions <- var_positions[var_positions[['snp']] %in% summary_stats[['SNP']], ]
  # make into hashtable
  var_positions_hm <- get_position_hashtab(var_positions, 'snp', 'pos')
  # get corrected values
  eigenmt_matrixeqtl <- eigenmt(summary_stats = summary_stats, genotypes = numeric_gens, genotype_to_position = var_positions_hm, var_explained_threshold = 0.975, window_size = 100)
  # read the emperical p values as well
  emperical_stats_loc <- 'resources/empiricalPvalues.txt.gz'
  emperical_stats <- data.table::fread(emperical_stats_loc, header = T, sep = '\t')
  # get the top value per feature
  eigenmt_matrixeqtl <- eigenmt_matrixeqtl[order(eigenmt_matrixeqtl[['p-value']]), ]
  eigenmt_matrixeqtl <- eigenmt_matrixeqtl[!duplicated(eigenmt_matrixeqtl[['gene']]), ]
  # join eigenMT and emperical p values
  both <- merge(eigenmt_matrixeqtl[, c('gene', 'feature_bf_eigen')], emperical_stats[, c('GENE', 'EMPIRICAL_PVALUE')], by.x = 'gene', by.y = 'GENE')
  # do -log10 of the p-value
  both[['minlog10_eigen']] <- -log10(both[['feature_bf_eigen']])
  both[['minlog10_perm']] <- -log10(both[['EMPIRICAL_PVALUE']])
  # correlate the values
  pval_cor <- cor(both$minlog10_eigen, both$minlog10_perm, method = 'spearman')
  # this should be more than .95
  expect_equal((pval_cor > .95), T)
})
