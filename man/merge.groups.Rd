\name{merge.groups}
\alias{merge.groups}
\title{merge groups}
\usage{
\method{merge}{groups}(groups, tags)
}
\arguments{
  \item{groups}{object of class groups}

  \item{tags}{index tags indicating the groups which should
  be merged}
}
\value{
object of class groups
}
\description{
Merge groups in an object of class groups
}
\details{
Groups identified by the index tag SNPs given by
\code{tags} will be merged, and the number of groups in the
object reduced.
}

