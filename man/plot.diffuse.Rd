\name{plot.diffuse}
\alias{plot.diffuse}
\title{Diffusion plot}
\usage{
  plot.diffuse(results, maxrank = 1000, thin = 500)
}
\arguments{
  \item{results}{an object of class snpmod, or a list of
  such objects}

  \item{maxrank}{truncate the x axis at this value}
}
\value{
  no return value
}
\description{
  A summary plot to indicate how diffuse the posterior is.
}
\details{
  This is measured by how rapidly the cumulative posterior
  probability increases with number of models.  For one
  dataset, this may be informative in an absolute sense,
  but comparing different datasets, by supplying a list of
  results, can be useful too.
}
\author{
  Chris Wallace
}

