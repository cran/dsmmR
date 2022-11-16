#' @title lambda genome
#' @description Contains the complete genome of the Escherichia phage Lambda.
#' @name lambda
#' @docType data
#' @usage data("lambda", package = "dsmmR")
#' data(lambda, package = "dsmmR") # equivalent.
#' # The following requires the package to be loaded,
#' # e.g. through `library(dsmmR)`.
#' data("lambda")
#' data(lambda)
#' @format A vector object of type \code{"character"} and length of 48502.
#'         It has class of \code{"Rdata"}.
#' @keywords datasets
#' @seealso \code{\link[utils:data]{data}}
#' @references
#' Sanger, F., Coulson, A. R., Hong, G. F., Hill, D. F., & Petersen, G. B.
#' (1982). Nucleotide sequence of bacteriophage \eqn{\lambda} DNA.
#' Journal of molecular biology, 162(4), 729-773.
#'
#' @examples
#' data("lambda", package = "dsmmR")
#' class(lambda)
#' sequence <- c(lambda) # Convert to "character" class
#' str(sequence)
NULL

