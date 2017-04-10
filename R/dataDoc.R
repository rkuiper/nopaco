#' Hypothetical data of outcomes for two risk models
#'
#' This data is generated as explained in the nopaco vignette. 
#' It represents the outcomes of the risk models (model A and model B). 
#' Both models were applied to gene expression profiles 100 subjects, each 
#' run in duplo. 
#' 
#'
#' @docType data
#'
#' @usage data(scores)
#'
#' @format A list with two elements named 'modelA' and 'modelB' both containing a 
#' dataframe with outcome scores for 100 subjects in the rows each having two 
#' replicate measurements in the columns.
#'
#' @keywords datasets
#'
#'
#' @source \code{vignette("nopaco", package = "nopaco")}
#'
#' @examples
#' data(scores)
#' str(scores)
#' plot(scores[['modelA']])
#' plot(scores[['modelB']])
"scores"

