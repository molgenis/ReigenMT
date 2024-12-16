get_position_hashtab <- function(feature_to_position, var_column, pos_column) {
  # make the hash table
  feature_to_position_hash <- r2r::hashmap()
  # data.table
  if (is.data.table(feature_to_position) | is.data.frame(feature_to_position)) {
    for (i in 1 : nrow(feature_to_position)) {
      # get feature
      feature <- feature_to_position[i, ..var_column][[var_column]][1]
      pos <- feature_to_position[i, ..pos_column][[pos_column]][1]
      feature_to_position_hash[[feature]] <- pos
    }
  }
  # if a list convert from that
  else if (is.list(feature_to_position)) {
    # convert
    for (var_name in names(feature_to_position)) {
      feature_to_position_hash[[var_name]] <- feature_to_position[[var_name]]
    }
  }
  # other formats
  else if (F == T){
    
  }
  return(feature_to_position_hash)
}

subset_genotypes <- function(genotypes, variants, id_column='ID') {
  # init variable
  genotypes_subset <- NULL
  # if already a datatable
  if (is.data.table(genotypes)) {
    # just subset
    genotypes_subset <- genotypes[genotypes[[id_column]] %in% variants, ]
    # and remove the ID column
    genotypes_subset[, c(id_column):=NULL]
  }
  # if matrix of data.frame, convert first
  else if(is.matrix(genotypes) | is.data.frame(genotypes)) {
    genotypes_subset <- data.table(genotypes)
    genotypes_subset <- genotypes[genotypes[[id_column]] %in% variants, ]
    # and remove the ID column
    genotypes_subset[, c(id_column):=NULL]
  }
  # if plink
  else if (F == T) {
    
  }
  # otherwise assume we can convert to data.table
  else {
    warning(paste('not recognizing genotype format, assuming matrix-compatible'))
    genotypes_subset <- data.table(genotypes)
    genotypes_subset <- genotypes[genotypes[[id_column]] %in% variants, ]
    # and remove the ID column
    genotypes_subset[, c(id_column):=NULL]
  }
  return(genotypes_subset)
}


genotypes_to_eigenvalues <- function(genotypes, verbose = F) {
  genotypes_t <- NULL
  # Transpose the genotype data if in a numeric format
  if (is.matrix(genotypes) | is.data.table(genotypes) | is.data.frame(genotypes)) {
    genotypes_t <- data.table(t(genotypes))
  }
  # or if in plink format
  else if(F == T) {
    
  }
  # we'll assume is a matrix-compatible format
  else {
    warning(paste('not recognizing genotype format, assuming matrix-compatible'))
    genotypes_t <- data.table(t(genotypes))
  }
  # make into double
  genotypes_t[, (names(genotypes_t)) := lapply(.SD, function(x){as.double(x)})]
  # get the means of each column
  var_means <- colMeans(genotypes_t, na.rm = T)
  # make var means into list
  var_means_list <- as.list(var_means)
  names(var_means_list) <- colnames(genotypes_t)
  # replace all the NAs
  for (col in names(var_means_list)) {
    setnafill(genotypes_t, type=c("const","locf","nocb"), fill=var_means_list[[col]], cols=col)
  }
  # Perform the fit using the corpcor package
  suppressMessages(fitted <- corpcor::cov.shrink(genotypes_t))
  # Extract the alpha (shrinkage intensity)
  alpha <- attributes(fitted)$lambda
  # get the correlation matrix
  cor_shrink <- cov2cor(fitted)
  # get the sum of values that are perfect
  number_of_perfects <- sum(cor_shrink == 1)
  # that means perfect LD, which is just one effect
  if (number_of_perfects == (ncol(cor_shrink) * nrow(cor_shrink))) {
    return(NA)
  }
  # Extract the shrinkage covariance matrix
  shrunk_cov <- matrix(fitted, ncol = ncol(genotypes_t), nrow = ncol(genotypes_t))
  # Get the variances of the variants (diagonal of shrinkage cov mat)
  variances <- diag(shrunk_cov)
  # Calculate inverse square root of the variances
  inv_sqrt_variances <- 1 / sqrt(variances)
  # Create a diagonal matrix of the inverse square root of the variances
  shrunk_precision <- diag(inv_sqrt_variances)
  # Calculate the shrunk correlation matrix
  shrunk_cor <- shrunk_precision %*% shrunk_cov %*% shrunk_precision
  # Compute the eigenvalues of the regularized (shrinkage) covariance matrix
  #eigenvalues <- eigen(shrunk_cor, symmetric = T)$values
  eigenvalues <- RSpectra::eigs_sym(as(shrunk_cor, "dgCMatrix"), k = ncol(shrunk_cor)-1, retvec = F)$values # k is the number of observations-1 so that we get almost all eigenvectors without it defaulting to eigs
  # order them in opposite
  eigenvalues <- sort(eigenvalues, decreasing = T)
  # Set negative eigenvalues to zero
  eigenvalues[eigenvalues < 0] <- 0
  return(eigenvalues)
}



find_number_of_eigen <- function(eigenvalues, n_variants, var_explained_threshold=.99, verbose=F) {
  # get the total variance
  total_variance_explained <- n_variants * var_explained_threshold
  # keep a sum of eigenvalues
  eigenvalues_sum <- 0
  # and keep track of how many eigenvalues we need
  eigenvalue_i <- 0
  while(eigenvalues_sum < total_variance_explained) {
    # update number first
    eigenvalue_i <- eigenvalue_i + 1
    # then add to the sum
    eigenvalues_sum <- eigenvalues_sum + eigenvalues[eigenvalue_i]
  }
  return(eigenvalue_i)
}

eigenmt <- function(summary_stats, genotypes, genotype_to_position, variant_column_summary_stats='SNP', feature_column_summary_stats='gene', pvalue_column='p-value', genotype_position_column='pos', genotype_var_position_column='snp', var_explained_threshold=.99, window_size=200) {
  # if not a data.table, make it one
  if (!is.data.table(summary_stats)) {
    summary_stats <- data.table(summary_stats)
  }
  # check if we already have a hashmap of positions
  genotype_to_pos_hm <- NULL
  if (is.hashtab(genotype_to_position)) {
    genotype_to_pos_hm <- genotype_to_position
  }
  # otherwise make genotype to position into hashmap
  else{
    genotype_to_pos_hm <- get_position_hashtab(genotype_to_position, genotype_var_position_column, genotype_position_column)
  }
  # get the unique features
  unique_features <- unique(summary_stats[[feature_column_summary_stats]])
  # remove any NA we might have
  unique_features <- unique_features[!is.na(unique_features)]
  message(paste('found', length(unique_features)), 'unique features')
  # do a pblapply
  corrected_per_feature <- pbapply::pblapply(unique_features, function(x) {
    # get the entries with this feature
    sumstats_feature <- summary_stats[!is.na(summary_stats[[feature_column_summary_stats]]) & summary_stats[[feature_column_summary_stats]] == x, ]
    # extract variants for this feature
    variants_feature <- sumstats_feature[[variant_column_summary_stats]]
    # get the position
    variant_positions <- as.vector(unlist(genotype_to_position[variants_feature]))
    # order the sumstats by this position
    sumstats_feature <- sumstats_feature[order(variant_positions), ]
    # extract variants for this feature
    variants_feature <- sumstats_feature[[variant_column_summary_stats]]
    # we'll check multiple windows, so we need to take into account the number of effects per window
    n_effects_windows <- 0
    # initialize the start and stop
    start <- 1
    end <- start + window_size - 1
    nvariants <- length(variants_feature)
    while(start < nvariants) {
      # if we are at the last variant set, we cannot do the full window
      if (end > nvariants) {
        end <- nvariants
      }
      # get the variants
      variants_window <- variants_feature[start : end]
      # now subset the genotype data to those variants
      variant_genotypes <- subset_genotypes(genotypes, variants_window)
      # get the eigenvalues for those variants
      variant_eigenvalues <- genotypes_to_eigenvalues(variant_genotypes)
      # check for NA, which is perfect LD
      n_tests <- NULL
      if (length(variant_eigenvalues) == 1) {
        n_tests <- 1
      }
      else {
        # get the number of eigenvalues we need
        n_tests <- NULL
        if (start ==1) {
          n_tests <- find_number_of_eigen(eigenvalues = variant_eigenvalues, n_variants = length(variants_window), var_explained_threshold = var_explained_threshold, verbose = T)
        }
        else {
          n_tests <- find_number_of_eigen(eigenvalues = variant_eigenvalues, n_variants = length(variants_window), var_explained_threshold = var_explained_threshold, verbose = F)
        }
        
      }
      # add that to the number of tests
      n_effects_windows <- n_effects_windows + n_tests
      # increase the start and stop
      start <- start + window_size
      end <- start + window_size - 1
    }
    # add this number of tests to the summary stats
    sumstats_feature[, ('n_tests_feature') := rep(n_effects_windows, times = nrow(sumstats_feature))]
  })
  # combine across features
  corrected_all <- do.call('rbind', corrected_per_feature)
  return(corrected_all)
}

