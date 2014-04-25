\name{abf.calc}
\alias{abf.calc}
\title{Calculate Approximate Bayes Factors}
\usage{
abf.calc(y, x, models, family = "binomial", q = NULL,
  method = c("glm.fit", "glm"), R2 = NULL, snp.data = NULL,
  return.R2 = FALSE, verbose = FALSE, parallel.file = NULL)
}
\arguments{
  \item{y}{response vector}

  \item{x}{explanatory variables, a SnpMatrix object,
  data.frame or numeric matrix}

  \item{models}{vector of models to consider}

  \item{family}{family argument to pass to glm.  Currently
  only "binomial" and "gaussian" are implemented for the
  \code{glm.fit} method.}

  \item{backend}{currently ignored}

  \item{method}{use glm.fit or glm to fit the models.
  Default is glm.fit which should be faster, but is less
  forgiving about the odd missing value etc, so if your
  code is giving potentially glm related errors, try
  running with \code{method="glm"} in the first instance.}

  \item{R2}{matrix giving pairwise r-squared measures of LD
  from which tag SNPs will be calculated when not directly
  available}

  \item{snp.data}{if R2 is missing, it is calculated from
  this SnpMatrix object}

  \item{return.R2}{if true, return the calculated R2
  matrix.  Useful if you are analysing several strata of a
  population and you wish to avoid repeating the
  calculation.}

  \item{verbose}{print lots of progress messages if TRUE.
  Default is FALSE.}
}
\value{
a data.frame containing model name, ABF, and an indicator
of whether the ABF was calculating directly or via a tag
SNP
}
\description{
Uses the BIC approximation to calculate approximate Bayes
factors for specified models
}
\author{
Chris Wallace
}

