setGeneric("length")
setGeneric("plot")
setGeneric("summary")
#' Create a union of groups, snppicker or tags objects
#'
#' First, objects are converted to class groups.  Then any groups
#' which overlap are identified, and replaced by their union.  Groups
#' which do not overlap are unchanged.  The combined set of groups is
#' returned.
#' @docType methods
#' @rdname union-methods
#' @author Chris Wallace
#' @export
#' @param x object of class \code{groups}, \code{snppicker} or \code{tags}
#' @param y object of same class as x
#' @return object of class groups
setGeneric("union", function(x,y) standardGeneric("union"))
#' Check whether a snp is in a snppicker, groups or tags object
#'
#' @docType methods
#' @rdname snpin-methods
#' @author Chris Wallace
#' @export
#' @param x character vector of SNPs to be checked
#' @param y object of class snppicker, groups or tags
#' @return logical matrix of nrow equal to length(x) and ncol equal to number of groups in y
setGeneric("snpin", function(x,y) standardGeneric("snpin"))
#setGeneric("snpdrop", function(x,y) standardGeneric("drop"))

##' Convert from old definitions of groups, tags classes to new
##' DON'T USE THIS FUNCTION UNLESS YOU HAVE OBJECTS STORED FROM A PREVIOUS PRE-DEVELOPMENT VERSION!
##' 
##' @param object 
##' @return new S4 structure
setGeneric("convert",function(object) standardGeneric("convert"))
##' Accessor function to show list of snps in an object of class groups
##'
##' @title snps
##' @param object 
##' @return a list of character vectors
##' @author Chris Wallace
setGeneric("snps",function(object) standardGeneric("snps"))
##' Accessor function to show tags from an object of class groups
##'
##' @title tags
##' @param object 
##' @return character vector of tag SNPs
##' @author Chris Wallace
setGeneric("tags",function(object) standardGeneric("tags"))
