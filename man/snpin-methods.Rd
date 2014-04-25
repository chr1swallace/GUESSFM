\docType{methods}
\name{snpin}
\alias{snpin}
\alias{snpin,character,groups-method}
\alias{snpin,character,snppicker-method}
\alias{snpin,character,tags-method}
\title{Check whether a snp is in a snppicker, groups or tags object}
\usage{
snpin(x, y)

\S4method{snpin}{character,snppicker}(x, y)

\S4method{snpin}{character,tags}(x, y)

\S4method{snpin}{character,groups}(x, y)
}
\arguments{
  \item{x}{character vector of SNPs to be checked}

  \item{y}{object of class snppicker, groups or tags}
}
\value{
logical vector of length equal to length(x)
}
\description{
Check whether a snp is in a snppicker, groups or tags
object
}
\author{
Chris Wallace
}

