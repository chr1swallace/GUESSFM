## https://gist.github.com/mrdwab/4205477
## LinearizeNestedList:
##
## https://sites.google.com/site/akhilsbehl/geekspace/
##         articles/r/linearize_nested_lists_in_r
##
## Akhil S Bhel
## 
## Implements a recursive algorithm to linearize nested lists upto any
## arbitrary level of nesting (limited by R's allowance for recursion-depth).
## By linearization, it is meant to bring all list branches emanating from
## any nth-nested trunk upto the top-level trunk s.t. the return value is a
## simple non-nested list having all branches emanating from this top-level
## branch.
##
## Since dataframes are essentially lists a boolean option is provided to
## switch on/off the linearization of dataframes. This has been found
## desirable in the author's experience.
##
## Also, one'd typically want to preserve names in the lists in a way as to
## clearly denote the association of any list element to it's nth-level
## history. As such we provide a clean and simple method of preserving names
## information of list elements. The names at any level of nesting are
## appended to the names of all preceding trunks using the `NameSep` option
## string as the seperator. The default `/` has been chosen to mimic the unix
## tradition of filesystem hierarchies. The default behavior works with
## existing names at any n-th level trunk, if found; otherwise, coerces simple
## numeric names corresponding to the position of a list element on the
## nth-trunk. Note, however, that this naming pattern does not ensure unique
## names for all elements in the resulting list. If the nested lists had
## non-unique names in a trunk the same would be reflected in the final list.
## Also, note that the function does not at all handle cases where `some`
## names are missing and some are not.
##
## Clearly, preserving the n-level hierarchy of branches in the element names
## may lead to names that are too long. Often, only the depth of a list
## element may only be important. To deal with this possibility a boolean
## option called `ForceNames` has been provided. ForceNames shall drop all
## original names in the lists and coerce simple numeric names which simply
## indicate the position of an element at the nth-level trunk as well as all
## preceding trunk numbers.
##
## Returns:
## LinearList: Named list.
##
LinearizeNestedList <- function(NList, LinearizeDataFrames=FALSE,
                                NameSep="/", ForceNames=FALSE) {
    ## Sanity checks:
    ##
    stopifnot(is.character(NameSep), length(NameSep) == 1)
    stopifnot(is.logical(LinearizeDataFrames), length(LinearizeDataFrames) == 1)
    stopifnot(is.logical(ForceNames), length(ForceNames) == 1)
    if (! is.list(NList)) return(NList)
    ##
    ## If no names on the top-level list coerce names. Recursion shall handle
    ## naming at all levels.
    ##
    if (is.null(names(NList)) | ForceNames == TRUE)
        names(NList) <- as.character(1:length(NList))
    ##
    ## If simply a dataframe deal promptly.
    ##
    if (is.data.frame(NList) & LinearizeDataFrames == FALSE)
        return(NList)
    if (is.data.frame(NList) & LinearizeDataFrames == TRUE)
        return(as.list(NList))
    ##
    ## Book-keeping code to employ a while loop.
    ##
    A <- 1
    B <- length(NList)
    ##
    ## We use a while loop to deal with the fact that the length of the nested
    ## list grows dynamically in the process of linearization.
    ##
    while (A <= B) {
        Element <- NList[[A]]
        EName <- names(NList)[A]
        if (is.list(Element)) {
            ##
            ## Before and After to keep track of the status of the top-level trunk
            ## below and above the current element.
            ##
            if (A == 1) {
                Before <- NULL
            } else {
                Before <- NList[1:(A - 1)]
            }
            if (A == B) {
                After <- NULL
            } else {
                After <- NList[(A + 1):B]
            }
            ##
            ## Treat dataframes specially.
            ##
            if (is.data.frame(Element)) {
                if (LinearizeDataFrames == TRUE) {
                    ##
                    ## `Jump` takes care of how much the list shall grow in this step.
                    ##
                    Jump <- length(Element)
                    NList[[A]] <- NULL
                    ##
                    ## Generate or coerce names as need be.
                    ##
                    if (is.null(names(Element)) | ForceNames == TRUE)
                        names(Element) <- as.character(1:length(Element))
                    ##
                    ## Just throw back as list since dataframes have no nesting.
                    ##
                    Element <- as.list(Element)
                    ##
                    ## Update names
                    ##
                    names(Element) <- paste(EName, names(Element), sep=NameSep)
                    ##
                    ## Plug the branch back into the top-level trunk.
                    ##
                    NList <- c(Before, Element, After)
                }
                Jump <- 1
            } else {
                NList[[A]] <- NULL
                ##
                ## Go recursive! :)
                ##
                if (is.null(names(Element)) | ForceNames == TRUE)
                    names(Element) <- as.character(1:length(Element))
                Element <- LinearizeNestedList(Element, LinearizeDataFrames,
                                               NameSep, ForceNames)
                names(Element) <- paste(EName, names(Element), sep=NameSep)
                Jump <- length(Element)
                NList <- c(Before, Element, After)
            }
        } else {
            Jump <- 1
        }
        ##
        ## Update book-keeping variables.
        ##
        A <- A + Jump
        B <- length(NList)
    }
    return(NList)
}

sparse.cor <- function(x){
    n <- nrow(x)
    cMeans <- colMeans(x)
    cSums <- colSums(x)
    ## Calculate the population covariance matrix.
    ## There's no need to divide by (n-1) as the std. dev is also calculated the same way.
    ## The code is optimized to minize use of memory and expensive operations
    covmat <- tcrossprod(cMeans, (-2*cSums+n*cMeans))
    crossp <- as.matrix(crossprod(x))
    covmat <- covmat+crossp
    sdvec <- sqrt(diag(covmat)) # standard deviations of columns
    covmat/crossprod(t(sdvec)) # correlation matrix
}

##' @importFrom Matrix sparseMatrix
makex <- function(obj) {
    snps <- sort(rownames(obj@snps))
    nums <- seq_along(snps)
    names(nums) <- snps
    npermod <- sapply(obj@model.snps,length)
    x <- sparseMatrix(i=rep(1:nrow(obj@models),times=npermod),
                      j=nums[unlist(obj@model.snps)],
                      x=rep(1,length(unlist(obj@model.snps))))
                                        #x=rep(obj@models$PP,times=npermod)^0.5)
    dimnames(x)[[2]] <- snps
    w <- obj@models$PP
    if(snps[[1]]=="1")
        x <- x[,-1]
    return(list(x=x,w=w))
}
makemppi <- function(x) {
    wx <- x$w*x$x
    colSums(wx)        
}
maker <- function(x) {
    sx <- sampx(x,n=1000000)
    r=sparse.cor(sx) #sqrt(w) * x) # sx
    diag(r) <-0 
    r[is.na(r)] <- 0
    r
}
maked <- function(r,r2) {
    ## rd <- dist(scale(cbind(r,r2)))
    ## rd <- dist(scale(r2))
    ## rd <- dist((sign(r)*(r2)))
    if(any(is.na(r2))) 
        r2[ which(is.na(r2)) ] <- 0
    if(any(r==0))
        r[ which(r==0) ] <- -0.001
    rd <- as.dist( ( sign(r) * (r2) + 1 ) / 2)
}
sampx <- function(d,n=1000) {
    wh <- sample(1:nrow(d$x),n,replace=TRUE,prob=d$w)
    d$x[wh,]
}
##' @title calculate max or min of subset of a matrix
##' @param r2 matrix, with row/column names which are members of L1, L2
##' @param L1 first list of vectors of snps
##' @param L2 second list vectors of snps
##' @param fun max or min, but could be another function that returns a scalar. unquoted.
##' @return smaller matrix indexed by elements of L1 (rows) and L2 (columns)
##' @export
##' @author Chris Wallace
calc.maxmin <- function(r2,L1,L2=L1,fun=max) {
    maxr2 <- matrix(0,length(L1),length(L2))
    for(i in seq_along(L1)) {
        for(j in seq_along(L2)) {
            maxr2[i,j] <- fun(r2[ L1[[i]], L2[[j]] ])
        }
    }
    maxr2
}

##' @title Group SNPs
##' @param SM2 snpmod object
##' @param snp.data SnpMatrix for r2 calculation
##' @param min.mppi trim snp groups with total MPPI < min.mppi in all
##'     diseases
##' @param r2.minmerge merge groups with minimum between-group r2 >
##'     r2.minmerge
##' @return list with three components.
##'
##' First is a data.frame with each row giving summary statistics for
##'     each group.
##'
##' Second is a groups object, each elements ordered according to the rows of the summary
##' 
##' Third is the r2 matrix calculated.
##' @export
##' @author Chris Wallace
group.multi <- function(SM2,snp.data,min.mppi=0.01,r2.minmerge=0.5) {
    stopifnot(is.list(SM2))
    bs <- best.snps(SM2,pp.thr=0)
    bs <- do.call("rbind",bs)
    snps <- unique(bs[bs$Marg_Prob_Incl>0.001,]$var) %>% setdiff(., "1")
    snp.data <- snp.data[,snps]
    r2 <- ld(snp.data,snp.data,stat="R.squared",symmetric = TRUE)
    X <- lapply(SM2,makex)
    MPPI <- lapply(X,makemppi) %>% do.call("cbind",.)
    R <- lapply(X, function(x) maker(x)[snps,snps])
    rmax <- rmin <- R[[1]]
    if(length(R)>1)
        for(i in 2:length(R))
            rmin <- pmin(rmin,R[[i]])
    rmax <- R[[1]]
    if(length(R)>1)
        for(i in 2:length(R))
            rmax <- pmax(rmax,R[[i]])
    ## plot(rmin,rmax)
    r <- ifelse(abs(rmin) > abs(rmax), rmin, rmax)
    rd <- maked(r,r2)
    h <- hclust(rd,method="complete")
    d <- as.dendrogram(h)
    r.tol=quantile(r,0.9)

    ## utility functions in local environment
    mem.sum <- function(members) { (colSums(MPPI[members,,drop=FALSE])) }
    mem.marg <- function(members) { sapply(X, function(x) {
        sum( pmin( apply(x$x[,members,drop=FALSE],1,sum), 1 ) * x$w ) })
    }
    mem.ab <- function(members) {list(a=mem.sum(members), b=mem.marg(members))}
    obj.ab <- function(object) {members <- labels(object); mem.ab(members)}
    mem.maxr.minr2 <- function(members) {
        r.sub <- r[members,members,drop=FALSE];
        r2.sub <- r2[members,members,drop=FALSE];
        mx <- max(r.sub[lower.tri(r.sub)],na.rm=TRUE) # max r
        mn <- min(r2.sub[upper.tri(r2.sub)],na.rm=TRUE) # min r2
        c(mx,mn)
    }
    cutter <- function (object,mppi.max=1.01,max.size=50,marg.sum.ratio=1.1,max.r=0,min.r2=0.5) {
        if(is.leaf(object))                 # Already at leaf, return
            return(labels(object))
        members <- labels(object)
        ab <- mem.ab(members)
        if(max(ab[[1]]) < min.mppi)             # no support for this sub, don't split further
            return(labels(object))
        if(max(ab[[1]]) > mppi.max)                # sum(mppi) must be <=1
            return(list(cutter(object[[1]]), cutter(object[[2]])))
        mxmn <- mem.maxr.minr2(members)
        if(mxmn[1] > r.tol || mxmn[2] < min.r2)                # within sub r2 too low
            return(list(cutter(object[[1]]), cutter(object[[2]])))
        if(min(c(mem.sum(labels(object[[1]])),mem.sum(labels(object[[2]])))<min.mppi)) # sub-sub is ignorable
            return(list(cutter(object[[1]]), cutter(object[[2]])))
        if(max(ab[[1]])<=mppi.max & all(ab[[1]]<ab[[2]]*marg.sum.ratio)) # no reason to continue
            return(labels(object))
        return(list(cutter(object[[1]]), cutter(object[[2]])))
    }
    
    mem.summ <- function(members) {
        n <- length(members)
        ab <- mem.ab(members)
        mppi.min <- apply(MPPI[members,,drop=FALSE],2,min)
        mppi.max <- apply(MPPI[members,,drop=FALSE],2,max)
        r2.sub <- r2[members,members]
        r2.summ <- summary(r2.sub[ upper.tri(r2.sub) ])
        r.sub <- r[members,members]
        r.summ <- summary(r.sub[ upper.tri(r.sub) ])
        c(n=n,sum.mppi=ab[[1]],r2=r2.summ["Min."],r2=r2.summ["Max."],
          r=r.summ["Min."],r=r.summ["Max."],
          mppi.min=mppi.min,mppi.max=mppi.max)
    }

    ## apply these functions
    ret <- cutter(d)
    if(!is.list(ret))
        ret <- list(ret)
    ret <- LinearizeNestedList(ret)
    ret.mppi <- sapply(ret, mem.sum) %>% t()
    use <- apply(ret.mppi,1,max) > 0.001
    df <- sapply(ret,mem.summ)
    df <- t(df)
    union.summary <- df[use,,drop=FALSE]
    union.content <- ret[use]

    ## step 1 - keep groups with sum.mppi > 0.01
    use <- apply(union.summary[,grep("sum.mppi",colnames(union.summary)),drop=FALSE],1,max) > min.mppi
    G1 <- union.summary[use,,drop=FALSE]
    G2 <- union.content[use]

    rownames(G1) <- NULL

    merger <- function(G1,G2) {
        maxr2 <- calc.maxmin(r2,G2,fun=max)
        minr2 <- calc.maxmin(r2,G2,fun=min)
        maxr <- calc.maxmin(r,G2,fun=max)
        diag(maxr2) <- 0
  
        tomerge <- maxr2>r2.minmerge  & maxr < r.tol #& minr2>0.4

        if(any(tomerge,na.rm=TRUE)) {
            wh <- which(tomerge,arr.ind=TRUE)
            wh <- wh[wh[,1] < wh[,2],,drop=FALSE ]
            wh <- cbind(wh,maxr2[wh])
            wh <- wh[order(wh[,3],decreasing=TRUE),,drop=FALSE]
            for(k in 1:nrow(wh)) {
                a <- wh[k,1]
                b <- wh[k,2]
                sumcols <- grep("sum.mppi",colnames(G1))
                if(any(colSums(G1[c(a,b),sumcols,drop=FALSE]) > 1.01))
                    next
                G2[[a]] <- c(G2[[a]], G2[[b]])
                G2[[b]] <- NULL
                for(nm in c(1,sumcols)) 
                    G1[a,nm] <- sum(G1[c(a,b),nm])
                for(nm in setdiff(1:ncol(G1), c(1,sumcols)))
                    G1[a,nm] <- max(G1[c(a,b),nm])
                G1 <- G1[-b,,drop=FALSE]
                return(merger(G1,G2))
            }
        }
        return(list(G1,G2))
    } 
    G <- merger(G1,G2)
    tmp <- G[[2]]
    names(tmp) <- sapply(tmp,"[[",1)
    newgroups <- new("groups",
                     tmp,
                     tags=names(tmp))
    return(list(summary=G[[1]],groups=newgroups,r2=r2))
}
