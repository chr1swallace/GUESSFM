################################################################################
## Models

#' Class to hold data relating to multiple models fitted to SNP data
#'
#' @slot snps data.frame containing marginal probabilities of inclusion for the SNPs
#' @slot models data.frame containing summaries for each model
#' @slot model.snps list containing the SNPs for each model. May be removed.
#'
#' @export
#' @examples
#' new("snpmod")
setClass("snpmod",
         slots=c(snps="data.frame",
                        models="data.frame",
                        model.snps="list"),
         validity=function(object) {
           if(nrow(object@models)!=length(object@model.snps))
             stop("Model summary should contain same number of models as model.snps decodes")
         })
#' Group focused class for holding information about sets of SNPs
#' defined by their mutual LD
#'
#' \code{groups} and \code{tags} are two structures for holding the
#' same information, depending on whether your focus is on the sets of
#' SNPs or their index members.  It is easy to convert one to another
#' and perhaps, in future, one class may be deprecated.
#'
#' @slot tags character vector giving tag SNPs.  Each tag indexes one group of SNPs
#' @slot .Data list of character vectors giving the SNP membership of each group
#' 
#' @export
#' @examples
#' new("groups")
setClass("groups",
         slots=c(tags="character"),
         contains="list",
         validity=function(object) {
           if(length(object@tags)!=length(object@.Data)) {
             stop("groups must be named by their tag")
           }
         })
#' Tags focused class for holding information about sets of SNPs
#' defined by their mutual LD
#'
#' \code{groups} and \code{tags} are two structures for holding the
#' same information, depending on whether your focus is on the sets of
#' SNPs or their index members.  It is easy to convert one to another
#' and perhaps, in future, one class may be deprecated.
#'
#' @slot tags character vector giving tag SNPs, one per SNP in
#' in \code{.Data}, repeated as necessary
#' @slot .Data character vector giving SNPs included in this tags object
#' @rdname groups-class
#' @examples
#' new("tags")
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
#' @slot groups list of data.frames describing SNPs in each group
#' @slot plotsdata list of additional data relating to the snp.picker
#' process that allows a summary of that process to be plotted via
#' \code{plot}.
#' new("snppicker")
setClass("snppicker",
         slots=c(plotsdata="list",groups="list"),
         validity=function(object) {
           if(length(object@groups)!=length(object@plotsdata))
             stop("groups and plotsdata should be lists of equal length")
         })

#' Class to hold results of pp.nsnp
#'
#' \code{pp.nsnp} summarises prior and posterior support for the
#' number of SNPs required to model a trait or several traits
#'
#' @slot .Data contains a named list of numeric vectors giving the support for each number
#' @slot plot contains a ggplot object
#' @slot traits names of traits, corresponding to items in .Data
#'
#' @export
#'
#' @examples
#' new("ppnsnp")
setClass("ppnsnp",
         slots=c(.Data="list",plot="ANY",traits="character"),
         validity=function(object) {
           if(length(object@.Data)!=length(object@traits))
             stop("traits should be same length as .Data")
           for(i in seq_along(object@.Data)) {
             if(!is(object[[i]],"array"))
               stop("ppnsnp should contain a list of arrays")
           }})
