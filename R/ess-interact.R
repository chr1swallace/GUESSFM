##' Create an object of class ESS given a base file name, and assuming the filenames created by run.bvs()
##'
##' A sensible base file name is the file ending in
##' _sweeps_features.txt.  If you used run.bvs() to run GUESS, then
##' this function should find all the correct files.  If you didn't,
##' then this function probably won't work, and you should use
##' as.ESS.object from the R2GUESS package directly.
##' @title ess.read
##' @param f base file name, a character string
##' @return object of class ESS
##' @export
##' @author Chris Wallace
ess.read <- function(f) {
    message("Reading from base file ",f)
    DIR <- paste0(dirname(f),"/")
    n <- gsub(".*/out_|_sweeps_features.txt","",f)
    suff <- function(str) {
        paste(str,n,sep="_")
    }
    as.ESS.object(dataY=suff("Y"),dataX=suff("X"),file.par="par.xml",command=FALSE,
                  path.input=DIR,path.output=DIR,path.par=DIR,
                  root.file.output=sprintf("out_%s_sweeps",n))
}

