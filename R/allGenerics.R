
##############################################################
#' Obtain concordance coefficients.
#'
#' \code{getPsi} returns the concordance coefficient(s) from a matrix or a result obtained by the \code{\link{concordance.test}} function.
#'
#' @param x A numeric matrix or an object \code{\link{ConcordanceTest-class}}
#' @param y A numeric matrix (optional)
#' @param ... Not used
#'
#' @return A numeric vector with coefficient(s)
#'
#' @family  concordance functions
#'
#' @export
#' @docType methods
#' @rdname getPsi-methods
#'
#' @examples
#'
#' matRandom <- matrix(rnorm(30),10,3)
#' testResult <- concordance.test(matRandom)
#' getPsi(testResult)
#' getPsi(matRandom)
#' @export
setGeneric("getPsi", function (x, y,...) {standardGeneric("getPsi")});

