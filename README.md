# GUESSFM

R package for fine mapping genetic associations in imputed GWAS data, detailed in http://biorxiv.org/content/early/2015/02/12/015164.

## Installation


Depends on GUESS available at http://www.bgx.org.uk/software/guess.html and described in the paper http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003657

To install, you may first need some R package dependencies, the packages VGAM, reshape, ggplot2, grid, ggbio, snpStats, parallel, R2GUESS.  Eg, if you don't have ggbio, from inside R, do

```R
install.packages("ggbio") 
```

and repeat for each package you don't have.  Then, to install GUESSFM direct from github, do

```R
install.packages("devtools") # if you don't already have the package
library(devtools)
install_github("chr1swallace/GUESSFM")
```

If you want the vignettes to build as well, you could do 
```R
install_github("chr1swallace/GUESSFM", build_vignettes=TRUE)
```

[Introduction](http://rawgit.com/chr1swallace/GUESSFM/master/inst/doc/introduction.html)

## Documentation

To load the GUESSFM library and see the vignette, do

```R
library(GUESSFM)
vignette(package="GUESSFM")
```

## Issues

If you run into any problems, please do raise an issue on the github
by clicking "Issues" on the top right of this page, and then the green
"New Issue" button.  Raising issues is better than email, because I
will have a reminder of what I need to fix every time I see the page,
and your issue may help others.  If you still prefer, you could email
me at chris.wallace at cimr.cam.ac.uk.


