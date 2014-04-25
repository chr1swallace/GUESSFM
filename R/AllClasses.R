################################################################################
## Models

#' Class to hold data relating to multiple models fitted to SNP data
#'
#' @section Slots: 
#' @slot snps data.frame containing marginal probabilities of inclusion for the SNPs
#' @slot models data.frame containing summaries for each model
#' @slot model.snps list containing the SNPs for each model. May be removed.
#'
#' @export
setClass("snpmod",
         slots=c(snps="data.frame",
                        models="data.frame",
                        model.snps="list"),
         validity=function(object) {
           if(nrow(object@models)!=length(object@model.snps))
             stop("Model summary should contain same number of models as model.snps decodes")
         })
#' groups class
#'
#' @section Slots:
#' @slot tags character vector giving tag SNPs.  Each tag indexes one group of SNPs
#' @slot .Data list of character vectors giving the SNP membership of each group
setClass("groups",
         slots=c(tags="character"),
         contains="list",
         validity=function(object) {
           if(length(object@tags)!=length(object@.Data)) {
             stop("groups must be named by their tag")
           }
         })
#' tags class
#'
#' @section Slots:
#' @slot tags character vector giving tag SNPs, one per SNP in \code{snps}, repeated as necessary
#' @slot .Data character vector giving SNPs included in this tags object
setClass("tags",
         slots=c(tags="character"),
         contains="character",
         validity=function(object) {
           if(length(object@tags)!=length(object@.Data))
             stop("tags must be in tags vector, tagging themselves")
         })

#' Class to hold results of snp.picker algorithm
#'
#' \code{snp.picker} groups SNPs according to LD and model inclusion
#' and outputs objects of class snppicker.
#'
#' @section Slots:
#' @slot groups list of data.frames describing SNPs in each group
#' @slot plotsdata list of additional data relating to the snp.picker
#' process that allows a summary of that process to be plotted via
#' \code{plot}.
setClass("snppicker",
         slots=c(plotsdata="list",groups="list"),
         validity=function(object) {
           if(length(object@.Data)!=length(object@plotsdata))
             stop("groups and plotsdata should be lists of equal length")
         })
