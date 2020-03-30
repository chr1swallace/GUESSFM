# GUESSFM

R package for fine mapping genetic associations in imputed GWAS data, detailed in 
> Wallace et al. (2015) Dissection of a Complex Disease Susceptibility Region Using a Bayesian Stochastic Search Approach to Fine Mapping. PLoS Genet 11(6): e1005272. doi: 10.1371/journal.pgen.1005272

which is available (open access) at http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005272



## Installatiou

### Short version

GUESSFM is purely an R package and so is platform independent.  It depends on the software GUESS, which is available at http://www.bgx.org.uk/software/guess.html.
GUESS is described in the paper http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003657

This is most easily installed via the R package R2GUESS, which has been removed from CRAN, but can still be installed via

```{sh}
wget https://cran.r-project.org/src/contrib/Archive/R2GUESS/R2GUESS_2.0.tar.gz
R CMD install R2GUESS_2.0.tar.gz 
```

Then GUESSFM can be installed from R by

```R
install.packages("devtools") # if you don't already have the package
library(devtools)
install_github("chr1swallace/GUESSFM")
```

### Longer version (if this fails)

<!-- GUESSFM is purely an R package and so is platform independent.  It depends on the software GUESS, which is available at http://www.bgx.org.uk/software/guess.html. -->

<!-- This can be more easily installed via the R package R2GUESS, but this has been removed from CRAN, and so this branch requires GUESS to be installed as described at the link above.  -->
 <!-- and should be installed according to the instructions there.   -->


To install GUESSFM you may first need some R package dependencies, the packages VGAM, reshape, ggplot2, grid, ggbio, snpStats, parallel.  E.g., if you don't have ggpplot, from inside R, do

```R
install.packages("ggplot2") 
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


