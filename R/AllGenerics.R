## generics from S3
setGeneric("length")
setGeneric("plot")
setGeneric("summary")

#' Create a union of groups, snppicker or tags objects
#'
#' First, objects are converted to class groups.  Then any groups
#' which overlap are identified, and replaced by their union.  Groups
#' which do not overlap are unchanged.  The combined set of groups is
#' returned.
#' @rdname union
#' @author Chris Wallace
#' @export
#' @param x object of class \code{groups}, \code{snppicker} or \code{tags}
#' @param y object of same class as x
#' @return object of class groups
setGeneric("union", function(x,y) standardGeneric("union"))

#' Check whether a snp is in a snppicker, groups or tags object
#'
#' @rdname snpin
#' @author Chris Wallace
#' @export
#' @param x character vector of SNPs to be checked
#' @param y object of class snppicker, groups or tags
#' @return logical matrix of nrow equal to length(x) and ncol equal to number of groups in y
setGeneric("snpin", function(x,y) standardGeneric("snpin"))
#setGeneric("snpdrop", function(x,y) standardGeneric("drop"))

##' Convert from old definitions of groups, tags classes to new
##' 
##' DON'T USE THIS FUNCTION UNLESS YOU HAVE OBJECTS STORED FROM A
##' PREVIOUS PRE-DEVELOPMENT VERSION!
##' 
#' @rdname conversion
##' @param object GUESSFM object from pre-development version
##' @return new S4 structure
setGeneric("convert",function(object) standardGeneric("convert"))

##' Accessor functions
##'
##' \code{snps} shows list of snps in an object of class groups and
##' returns a list of character vectors
##'
##' @title Accessors for groups objects
##' @param object object from which items should be extracted
##' @return a list of character vectors
##' @author Chris Wallace
##' @rdname accessors
##' @export
setGeneric("snps",function(object) standardGeneric("snps"))

##' @details \code{tags} shows tags from an object of class groups and
##' returns a character vector of tag SNPs
##' @rdname accessors
##' @export
setGeneric("tags",function(object) standardGeneric("tags"))

##' @details \code{tagsof} shows tags for a named character vector of SNPs
##' @rdname groups-subset
##' @export
##' @param object tags object
##' @param i character vector of SNPs
##' @return data.frame of tags and their tagged SNPs
setGeneric("tagsof", function(object,i) standardGeneric("tagsof"))

##' @details \code{taggedby} shows SNPs tagged by a named character vector of tag SNPs
##' @rdname groups-subset
##' @export
##' @param object tags object
##' @param i character vector of tag SNPs
##' @return data.frame of tags and their tagged SNPs
setGeneric("taggedby", function(object,i) standardGeneric("taggedby"))
