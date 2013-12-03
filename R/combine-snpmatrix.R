##' calculate intersection and setdiffs for two vectors
##'
##' Originally from randomFunctions
##' @title textvenn
##' @param A a vector
##' @param B a vector of the same type as A
##' @param quiet 
##' @return a list with three elements: A is setdiff(A,B).  int is intersect(A,B).  B is setdiff(B,A)
##' @author Chris Wallace
##' @export
##' @examples
##' A <- 1:3
##' B <- 2:4
##' textvenn(A,B)
##' tv <- textvenn(A,B,quiet=TRUE)
##' tv
textvenn <- function(A,B,quiet=FALSE) {
  if(!is(B,class(A)))
    stop("A and B should be of same class")
  AnotB <- setdiff(A,B)
  AandB <- intersect(B,A)
  BnotA <- setdiff(B,A)
  if(!quiet) {
    cat("set A:\n")
    cat(AnotB,"\n")
    cat("intersection:\n")
    cat(AandB,"\n")
    cat("set B:\n")
    cat(BnotA,"\n")
  }
  invisible(list(A=AnotB,int=AandB,B=BnotA))
}
##' Bind two SnpMatrices with overlapping SNPs but distinct samples.
##'
##' BEWARE: it is assumed that alleles are aligned, and NO CHECK is
##' made.  Combining SnpMatrices with alleles aligned to opposite
##' strands will produce a meaningless result.
##'
##' @title snpmatrix.combine
##' @param X a SnpMatrix
##' @param Y a SnpMatrix with some snps overlapping X, but no samples overlapping
##' @return a new SnpMatrix formed from X and Y.  SNPs found in only one dataset will have missing genotypes in the other
##' @author Chris Wallace
##' @export
snpmatrix.combine <- function(X,Y) {
  if(!is(X,"SnpMatrix") || !is(Y,"SnpMatrix"))
    stop("X and Y must be of class SnpMatrix")
  if(length(intersect(rownames(X),rownames(Y))))
    stop("X and Y should have distinct samples")
  sv <- textvenn(colnames(X),colnames(Y),quiet=TRUE)
  if(!length(sv$int))
    warning("no overlapping SNPs found.")
  int <- rbind(X[,sv$int],Y[,sv$int])
  A <- rbind(X[,sv$A],
           new("SnpMatrix",
               matrix(as.raw(0),nrow=nrow(Y),ncol=length(sv$A),
                      dimnames=list(rownames(Y),sv$A))))
  B <- rbind(new("SnpMatrix",
                 matrix(as.raw(0),nrow=nrow(X),ncol=length(sv$B),
                        dimnames=list(rownames(X),sv$B))),
             Y[,sv$B])
  gdata <- new("SnpMatrix",cbind(int,A,B))

  return(gdata)
}
