##' Uses snpStats' snp.imputation function
##'
##' @title Fill in missing values in SnpMatrix object
##' @param X SnpMatrix object
##' @param bp optional vector of base pair positions, length == number of SNPs
##' @param strata optional vector, length == number of samples, as.factor(strata) defines distinct strata of samples, eg by geography
##' @param numeric default FALSE
##' @param ... arguments passed on to \code{impute.snps}
##' @return matrix, either numeric (if numeric=TRUE) or SnpMatrix (if numeric=FALSE)
##' @author Chris Wallace
##' @export
impute.missing <- function(X,bp=1:ncol(X),strata=NULL, numeric=FALSE, ...) {
  N <- as(X,"numeric")
  if(!is.null(strata)) {
    strata <- as.factor(strata)
    if(length(levels(strata))>10)
      stop("too many levels in strata\n")
    for(i in levels(strata)) {
      cat("\nstrata",i,"\n")
      wh <- which(strata==i)
      N[wh,] <- impute.missing(X[wh,,drop=FALSE],bp, numeric=TRUE, ...)
    }
  } else {
    csumm <- col.summary(X)
    use <- csumm[,"Certain.calls"]==1
    X2 <- X[,use]
    bp <- bp[use]
    imp <- (csumm[,"Call.rate"]<1)[use]
    cat(sum(imp),"to impute\n")
    for(i in which(imp)) {
      cat(i,".")
      rule <- snp.imputation(X2[,-i,drop=FALSE],X2[,i,drop=FALSE],bp[-i],bp[i])
      if(is.null(rule@.Data[[1]]))
        next
      imp <- impute.snps(rules=rule,snps=X2[,rule@.Data[[1]]$snps,drop=FALSE], ...)
      wh.na <- which(is.na(N[,i]))
      N[wh.na, colnames(X2)[i]] <- imp[wh.na]
    }
    cat("\n")
  }
  if(numeric) {
    return(N)
  } else {
    return(new("SnpMatrix",
               data=round(N)+1, nrow=nrow(N), ncol=ncol(N), dimnames=dimnames(N)))
  }
}
