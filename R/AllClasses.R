################################################################################
## Models

#' Class to hold data relating to multiple models fitted to SNP data
#'
#' @title The snpmod class
#' @section Slots: 
#'  \describe{
#'    \item{\code{snps}:}{data.frame containing marginal probabilities of inclusion for the SNPs}
#'    \item{\code{models}:}{data.frame containing summaries for each model}
#'    \item{\code{model.snps}:}{list containing the SNPs for each model. May be removed.}
#'  }
#'
#' @name snpmod
#' @rdname snpmod-class
#' @aliases snpmod-class
#' @export
#' @author Chris Wallace
setClass("snpmod",
         representation(snps="data.frame",
                        models="data.frame",
                        model.snps="list"),
         validity=function(object) {
           if(nrow(object@models)!=length(object@model.snps))
             stop("Model summary should contain same number of models as model.snps decodes")
         })

setClass("snppicker",
         representation(groups="list",
                        plotsdata="list"),
         validity=function(object) {
           if(length(object@groups)!=length(object@plotsdata))
             stop("groups and plotsdata should be lists of equal length")
         })

setClass("groups",
         representation(tags="character",groups="list"),
         validity=function(object) {
           if(length(object@tags)!=length(object@groups)) {
             stop("groups must be named by their tag")
           }
         })
setClass("tags",
         representation(tags="character",snps="character"),
         validity=function(object) {
           if(length(object@tags)!=length(object@snps))
             stop("tags must be in tags vector, tagging themselves")
         })

