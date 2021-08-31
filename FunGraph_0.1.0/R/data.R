#' genotypic data
#'
#' A dataset containing genotypic data of 106 Populus euphratica samples
#'
#'@format A data frame with 8305 rows and 108 variables, contain row names as SNP ID:
#'@docType data
#'@name geno
#'@usage data(geno)
"geno"

#' phenotypic data
#'
#' A dataset containing phenotypic data of 106 Populus euphratica samples
#' under salt-free and salt-stress environment measured 14 times
#'
#'@format A data frame with 106 rows and 28 variables, contain row names as sample ID:
#'@docType data
#'@name pheno
#'@usage data(pheno)
"pheno"

#' LR(log-likelihood-ratio) data
#'
#' A dataset containing calculated LR value from geno and pheno in this dataset
#'
#' @format A vector length = 8305
"LR"

#' generic effect data
#'
#' A dataset containing calculated generic effect from geno and pheno in this dataset
#'
#' @format A data frame of 28 columns and 8305 SNPs
"generic_effect"
