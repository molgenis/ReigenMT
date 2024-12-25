# ReigenMT
R implementation of eigenMT algorithm for multiple testing correction


## about EigenMT

EigenMT multiple testing correction is based on quantifying the number of independent for each feature by performing eigenvalue decomposition on variants tested for that feature. The number of eigenvalues required to explain (almost) all variance across these variants is considered the number of independent tests that were performed (sort of a proxy for LD structure) for the feature. This number of independent tests can then be used to perform bonferroni correction, instead of canonical bonferroni correction where the alpha is determined on the basis of all tests. This method of multiple testing correction is less strict than canonical bonferroni correction and has been shown to be comparable to permutation-based multiple testing correction, while being computationally much less expensive. Some caution is advised though when dealing with data that has large outliers, as in these cases bonferroni can even result in lower p-values than permutations, and as such the bonferroni-based eigenMT implementation will likely do the same.

The original paper can be found here: https://doi.org/10.1016/j.ajhg.2015.11.021, and all credits with regards to the original implementation go to the original authors.


## about ReigenMT

ReigenMT is an R implementation of the original eigenMT algorithm, based on the code here: https://github.com/molgenis/eigenMT, which was an extension of the original source code here: https://github.com/joed3/eigenMT

ReigenMT currently supports plink format genotypes, as well as MatrixEQTL format genotypes. Annotation files for genomic positions of genetic variants are also required in the case of MatrixEQTL (MatrixEQTL format), for plink genotypes, these positions are inferred from the genotype files themselves. ReigenMT (like the original implementation) is currently set up to be used exclusively in the context of cis-eQTL mapping, or more specifically, the genetic variants tested per feature are assumed to be on the same chromosome. Support for variants on multiple chromosomes for the same feature, may be implemented at a later date.


## installation

ReigenMT can be installed via devtools::install_github in the following manner:

```r
devtools::install_github("https://github.com/molgenis/ReigenMT");
```

## usage

For MatrixEQTL data, the package can be used as follows:

```r
# numeric genotypes
numeric_gens_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/interaction_eqtl/scripts/eigenMT_testdata/genotypes.txt'
numeric_gens <- data.table::fread(numeric_gens_loc, header = T)
# summary stats
summary_stats_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/interaction_eqtl/scripts/eigenMT_testdata/cis.eqtls.txt'
summary_stats <- data.table::fread(summary_stats_loc, header = T)
# variant positions
var_positions_loc <- '/groups/umcg-franke-scrna/tmp04/projects/multiome/ongoing/qtl/interaction_eqtl/scripts/eigenMT_testdata/gen.positions.txt'
var_positions <- data.table::fread(var_positions_loc, header = T)
# make into hashtable
var_positions_hm <- get_position_hashtab(var_positions, 'snp', 'pos')
# get corrected values
eigenmt_matrixeqtl <- eigenmt(summary_stats = summary_stats, genotypes = numeric_gens, genotype_to_position = var_positions_hm, var_explained_threshold = 0.975, window_size = 200)
```

For example LIMIX-QTL h5 inputs and plink genotypes, the package can be used as follows:

```r
# genotypes
genotypes_loc <- '/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_multiome/genotype_input/EUR_imputed_hg38_varFiltered_chr7'

# read the chromsome annotations
genotypes <- read.plink(
  bed = paste(genotypes_loc, '.bed', sep = ''),
  bim = paste(genotypes_loc, '.bim', sep = ''),
  fam = paste(genotypes_loc, '.fam', sep = ''),
)
# get the variant locations
var_locations_plink <- get_position_hashtab(genotypes)
# get the summary stats
summary_stats_h5_loc <- '/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_multiome/output/L1/UT/Mono/qtl/qtl_results_7_67421326_78719150.h5'
summary_stats_h5 <- limix_h5_to_sumstats_format(summary_stats_h5_loc)
# get the corrected values
eigenmt_limix <- eigenmt(summary_stats = summary_stats_h5, genotypes = genotypes, genotype_to_position = var_locations_plink, variant_column_summary_stats = 'snp_id', feature_column_summary_stats = 'feature', var_explained_threshold = 0.975, pvalue_column = 'p_value')
```

The get_position_hashtab is optional, if the 'genotype_to_position' parameter is not supplied, this hastab will be calculated automatically. This step however might however be lengthy, depending on the number of variants, and as such it can be better to calculate this beforehand, to prevent having to create it multiple times.

The 'eigenmt' function will return the original summary stats, with the number of tests per feature ('n_tests_feature'), the p-value corrected for the number of independent tests of a feature ('feature_bf_eigen') and the p-value corrected for the number of independent tests across all features ('total_bf_eigen')
