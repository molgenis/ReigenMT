#' Convert HDF5 to Summary Statistics Format
#'
#' This function reads an HDF5 file and converts its contents into a data.table
#' with an additional column indicating the feature name. This assumes the format of chunked LIMIX-QTL outputs
#'
#' @param h5_loc A string specifying the location of the HDF5 file.
#' @return A data.table containing the merged data from all features in the HDF5 file,
#' with an additional column named 'feature' indicating the feature name.
#' @import data.table
#' @importFrom rhdf5 h5ls h5read
#' @export
#' @examples
#' \dontrun{
#' h5_file <- "path/to/your/file.h5"
#' result <- limix_h5_to_sumstats_format(h5_file)
#' print(result)
#' }
limix_h5_to_sumstats_format <- function(h5_loc) {
  # get the contents
  h5_contents <- data.frame(rhdf5::h5ls(h5_loc))
  # put in a list
  all_features <- list()
  # check each feature
  for (feature in h5_contents[['name']]) {
    # get the data
    feature_dt <- data.table::data.table(rhdf5::h5read(h5_loc, feature))
    # add the feature as a column
    feature_dt[, ('feature') := rep(feature, times = nrow(feature_dt))]
    # put in the list
    all_features[[feature]] <- feature_dt
  }
  # merge all of them
  all_features_dt <- do.call('rbind', all_features)
  return(all_features_dt)
}


#' Create a Position Hash Table
#'
#' This function creates a hash table mapping features to their positions from various input formats. This table is required in some functions, to link variants to genomic positions.
#'
#' @param feature_to_position A data.table, data.frame, list, or other recognized format containing feature and position information.
#' @param var_column The name of the column containing feature names.
#' @param pos_column The name of the column containing position values.
#' @return A hash table mapping each feature to its corresponding position.
#' @import data.table
#' @importFrom r2r hashmap
#' @export
#' @examples
#' \dontrun{
#' feature_data <- data.table(feature = c("feature1", "feature2"), position = c(100, 200))
#' hash_table <- get_position_hashtab(feature_data, "feature", "position")
#' print(hash_table)
#' }
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
  # other formats
  else if ((is.list(feature_to_position)) & (class(feature_to_position$feature_to_position) == 'SnpMatrix')) {
    # get this information from the map
    for (i in 1 : nrow(feature_to_position$map)) {
      feature_to_position_hash[[feature_to_position$map$snp.name[i]]] <- feature_to_position$map$position[i]
    }
  }
  # if a conventional list convert from that
  else if (is.list(feature_to_position)) {
    # convert
    for (var_name in names(feature_to_position)) {
      feature_to_position_hash[[var_name]] <- feature_to_position[[var_name]]
    }
  }
  else {
    warning('could not interpret as a known filetype, will try a list-like approach')
    # convert
    for (var_name in names(feature_to_position)) {
      feature_to_position_hash[[var_name]] <- feature_to_position[[var_name]]
    }
  }
  return(feature_to_position_hash)
}


#' Subset Genotypes by Variants
#'
#' This function subsets genotype data based on a list of variants and optionally removes the ID column.
#'
#' @param genotypes A data.table, data.frame, matrix, or list containing genotype data.
#' @param variants A vector of variant IDs to subset the genotypes.
#' @param id_column A string specifying the name of the column containing variant IDs. Default is 'ID'.
#' @param strip_id A boolean whether or not to remove the identifier from the table. Default is True.
#' @return A subset of the genotype data containing only the specified variants, with the ID column removed.
#' @import data.table
#' @importFrom methods is
#' @export
#' @examples
#' \dontrun{
#' genotypes <- data.table(ID = c("var1", "var2", "var3"), value = c(1, 2, 3))
#' variants <- c("var1", "var3")
#' subset <- subset_genotypes(genotypes, variants, "ID")
#' print(subset)
#' }
subset_genotypes <- function(genotypes, variants, id_column='ID', strip_id=T) {
  # init variable
  genotypes_subset <- NULL
  # if already a datatable
  if (is.data.table(genotypes)) {
    # just subset
    genotypes_subset <- genotypes[genotypes[[id_column]] %in% variants, ]
    # and remove the ID column
    if (strip_id) {
      genotypes_subset[, c(id_column):=NULL]
    }
  }
  # if matrix of data.frame, convert first
  else if(is.matrix(genotypes) | is.data.frame(genotypes)) {
    genotypes_subset <- data.table(genotypes)
    genotypes_subset <- genotypes[genotypes[[id_column]] %in% variants, ]
    # and remove the ID column
    if (strip_id) {
      genotypes_subset[, c(id_column):=NULL]
    }
  }
  # if plink
  else if ((is.list(genotypes)) & (class(genotypes$genotypes) == 'SnpMatrix')) {
    # get the indices that we want to keep
    keep_indices <- genotypes$map$snp.name %in% variants
    genotypes_subset <- list(
      'fam' = genotypes$fam,
      'map' = genotypes$map[keep_indices, ],
      'genotypes' = genotypes$genotypes[, keep_indices]
    )
  }
  # otherwise assume we can convert to data.table
  else {
    warning(paste('not recognizing genotype format, assuming matrix-compatible'))
    genotypes_subset <- data.table(genotypes)
    genotypes_subset <- genotypes[genotypes[[id_column]] %in% variants, ]
    # and remove the ID column
    if (strip_id) {
      genotypes_subset[, c(id_column):=NULL]
    }
  }
  return(genotypes_subset)
}


#' Compute Eigenvalues from Genotype Data
#'
#' This function computes the eigenvalues of the regularized (shrinkage) covariance matrix from genotype data.
#'
#' @param genotypes A data.table, data.frame, matrix, or list containing genotype data.
#' @param verbose A logical value indicating whether to print detailed messages. Default is FALSE.
#' @return A vector of eigenvalues of the shrinkage covariance matrix.
#' @import data.table
#' @importFrom corpcor cov.shrink
#' @importFrom spam eigen.spam
#' @importFrom stats cov2cor
#' @export
#' @examples
#' \dontrun{
#' genotypes <- data.table(ID = c("var1", "var2", "var3"), value = c(1, 2, 3))
#' eigenvalues <- genotypes_to_eigenvalues(genotypes)
#' print(eigenvalues)
#' }
genotypes_to_eigenvalues <- function(genotypes, verbose = F) {
  genotypes_t <- NULL
  # Transpose the genotype data if in a numeric format
  if (is.matrix(genotypes) | is.data.table(genotypes) | is.data.frame(genotypes)) {
    genotypes_t <- data.table::data.table(t(genotypes))
  }
  # or if in plink format
  else if(is.list(genotypes) & class(genotypes$genotypes) == 'SnpMatrix') {
    genotypes_t <- data.table::data.table(as(genotypes$genotypes, "numeric"))
  }
  # we'll assume is a matrix-compatible format
  else {
    warning(paste('not recognizing genotype format, assuming matrix-compatible'))
    genotypes_t <- data.table::data.table(t(genotypes))
  }
  # Redirect output to null to suppress print statements
  sink(tempfile())
  # Ensure to reset sink after function execution
  on.exit(sink())
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
  cor_shrink <- stats::cov2cor(fitted)
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
  eigenvalues <- NULL
  # eigen_sym needs at least 3x3
  # if (nrow(shrunk_cor) > 100) {
  #   eigenvalues <- RSpectra::eigs_sym(as(shrunk_cor, "dgCMatrix"), k = ncol(shrunk_cor)-1, retvec = F)$values # k is the number of observations-1 so that we get almost all eigenvectors without it defaulting to eigs
  # }
  # else {
  #   eigenvalues <- eigen(shrunk_cor, symmetric = T)$values
    eigenvalues <- spam::eigen.spam(shrunk_cor, nev = ncol(shrunk_cor), symmetric = T)$values
  # }
  # order them in opposite
  eigenvalues <- sort(eigenvalues, decreasing = T)
  # Set negative eigenvalues to zero
  eigenvalues[eigenvalues < 0] <- 0
  return(eigenvalues)
}


#' Find Number of Eigenvalues to Explain Variance
#'
#' This function determines the number of eigenvalues needed to explain a specified threshold of variance.
#'
#' @param eigenvalues A numeric vector of eigenvalues.
#' @param n_variants An integer specifying the number of variants.
#' @param var_explained_threshold A numeric value indicating the threshold of variance to be explained. Default is 0.975.
#' @param verbose A logical value indicating whether to print detailed messages. Default is FALSE.
#' @return An integer representing the number of eigenvalues needed to explain the specified variance threshold.
#' @export
#' @examples
#' \dontrun{
#' eigenvalues <- c(2.5, 1.5, 1.0, 0.5)
#' n_variants <- 4
#' threshold <- 0.975
#' num_eigen <- find_number_of_eigen(eigenvalues, n_variants, threshold)
#' print(num_eigen)
#' }
find_number_of_eigen <- function(eigenvalues, n_variants, var_explained_threshold=.975, verbose=F) {
  # get the total variance
  total_variance_explained <- n_variants * var_explained_threshold
  # keep a sum of eigenvalues
  eigenvalues_sum <- 0
  # and keep track of how many eigenvalues we need
  eigenvalue_i <- 0
  # now keep adding eigenvalues until we either get to the variance threshold, or run out of eigenvalues
  while((eigenvalue_i != length(eigenvalues)) & (eigenvalues_sum < total_variance_explained)) {
    # update number first
    eigenvalue_i <- eigenvalue_i + 1
    # then add to the sum
    eigenvalues_sum <- eigenvalues_sum + eigenvalues[eigenvalue_i]
  }
  # if we didnt explain all variance, all effects were independent, and we need one more
  if (eigenvalues_sum < total_variance_explained) {
    eigenvalue_i <- eigenvalue_i + 1
  }
  return(eigenvalue_i)
}


#' EigenMT: Multiple Testing Correction for Summary Statistics
#'
#' This function performs multiple testing correction on summary statistics using eigenvalues of genotype data.
#'
#' @param summary_stats A data.table or data.frame containing summary statistics.
#' @param genotypes A data.table, data.frame, matrix, or list containing genotype data.
#' @param genotype_to_position A data.table, data.frame, list, or hash table mapping genotypes to positions.
#' @param variant_column_summary_stats A string specifying the column name for variants in the summary statistics. Default is 'SNP'.
#' @param feature_column_summary_stats A string specifying the column name for features in the summary statistics. Default is 'gene'.
#' @param pvalue_column A string specifying the column name for p-values in the summary statistics. Default is 'p-value'.
#' @param genotype_position_column A string specifying the column name for positions in the genotype data. Default is 'pos'.
#' @param genotype_var_position_column A string specifying the column name for variant positions in the genotype data. Default is 'snp'.
#' @param var_explained_threshold A numeric value indicating the threshold of variance to be explained. Default is 0.975.
#' @param window_size An integer specifying the window size for checking multiple variants. Default is 100.
#' @return A data.table with the corrected summary statistics, including the number of tests per feature.
#' @import data.table
#' @importFrom pbapply pblapply
#' @export
#' @examples
#' \dontrun{
#' summary_stats <- data.table(SNP = c("rs1", "rs2"), gene = c("gene1", "gene2"), `p-value` = c(0.01, 0.05))
#' genotypes <- data.table(ID = c("var1", "var2"), value = c(1, 2))
#' genotype_to_position <- data.table(snp = c("rs1", "rs2"), pos = c(100, 200))
#' corrected_stats <- eigenmt(summary_stats, genotypes, genotype_to_position)
#' print(corrected_stats)
#' }
eigenmt <- function(summary_stats, genotypes, genotype_to_position, variant_column_summary_stats='SNP', feature_column_summary_stats='gene', pvalue_column='p-value', genotype_position_column='pos', genotype_var_position_column='snp', var_explained_threshold=.975, window_size=100) {
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
  else if(!is.null(genotype_to_position) & !is.hashtab(genotype_to_position)){
    genotype_to_pos_hm <- get_position_hashtab(genotype_to_position, genotype_var_position_column, genotype_position_column)
  }
  # if no positions supplied, but we have plink format, we can interpret it ourselves
  else if(is.null(genotype_to_position) & is.list(genotypes) & class(genotypes$genotypes) == 'SnpMatrix') {
    genotype_to_pos_hm <- get_position_hashtab(genotypes)
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
    # get the number of variants
    nvariants <- length(variants_feature)
    # we'll check multiple windows, so we need to take into account the number of effects per window
    n_effects_windows <- 0
    # but only if we have multiple variants
    if (nvariants > 1) {
      # initialize the start and stop
      start <- 1
      end <- start + window_size - 1

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
          n_tests <- find_number_of_eigen(eigenvalues = variant_eigenvalues, n_variants = length(variants_window), var_explained_threshold = var_explained_threshold, verbose = F)
        }
        # add that to the number of tests
        n_effects_windows <- n_effects_windows + n_tests
        # increase the start and stop
        start <- start + window_size
        end <- start + window_size - 1
      }
    }
    # otherwise it is just one effect
    else {
      n_effects_windows <- 1
    }
    # add this number of tests to the summary stats
    sumstats_feature[, ('n_tests_feature') := rep(n_effects_windows, times = nrow(sumstats_feature))]
    return(sumstats_feature)
  })
  # combine across features
  corrected_all <- do.call('rbind', corrected_per_feature)
  return(corrected_all)
}
