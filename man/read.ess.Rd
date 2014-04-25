\name{read.ess}
\alias{read.ess}
\title{ess.read}
\usage{
read.ess(f, ...)
}
\arguments{
  \item{f}{base file name, a character string}
}
\value{
object of class ESS
}
\description{
Create an object of class ESS given a base file name
}
\details{
A sensible base file name is the file ending in
_sweeps_features.txt.  If you used run.bvs() to run GUESS,
then this function should find all the correct files.  If
you didn't, then this function probably won't work, and you
should use as.ESS.object from the R2GUESS package directly.
}
\author{
Chris Wallace
}

