# GUESSFM

R package for fine mapping genetic associations in imputed GWAS data, detailed in 
> Wallace et al. (2015) Dissection of a Complex Disease Susceptibility Region Using a Bayesian Stochastic Search Approach to Fine Mapping. PLoS Genet 11(6): e1005272. doi: 10.1371/journal.pgen.1005272

which is available (open access) at http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005272



## Installation

GUESSFM is purely an R package and so is platform independent.  It depends on the software GUESS, which is available at http://www.bgx.org.uk/software/guess.html and should be installed according to the instructions there.  GUESS is described in the paper http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003657

To install, you may first need some R package dependencies, the packages VGAM, reshape, ggplot2, grid, ggbio, snpStats, parallel, R2GUESS.  E.g., if you don't have ggpplot, from inside R, do

```R
install.packages("ggplot") 
```

and repeat for each package you don't have.  

Some packages (e.g. ggbio, snpStats) are from Bioconductor.  For these, you need to do

```R
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("ggbio")
biocLite("snpStats")
```

Then, to install GUESSFM direct from github, do

```R
install.packages("devtools") # if you don't already have the package
library(devtools)
install_github("chr1swallace/GUESSFM")
```

If you want the vignettes to build as well, this *should* work, but doesn't (even though they build locally):
```R
install_github("chr1swallace/GUESSFM", build_vignettes=TRUE)
```

Instead, see the vignettes direct from github:

[Introduction](http://rawgit.com/chr1swallace/GUESSFM/master/inst/doc/introduction.html) | 
[Tags/groups](http://rawgit.com/chr1swallace/GUESSFM/master/inst/doc/groups.html) | 
[Plotting](http://rawgit.com/chr1swallace/GUESSFM/master/inst/doc/plotting.html)

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
me at cew54 at cam.ac.uk.


