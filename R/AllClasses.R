################################################################################
## Models

setClass("snpmod",
         representation(snps="data.frame",
                        models="data.frame",
                        model.snps="list"),
         validity=function(object) {
           if(nrow(object@models)!=length(object@model.snps))
             stop("Model summary should contain same number of models as model.snps decodes")
         })

setClass("snppicker",
         representation(groups="list",
                        plotsdata="list"),
         validity=function(object) {
           if(length(object@groups)!=length(object@plotsdata))
             stop("groups and plotsdata should be lists of equal length")
         })

setClass("groups",
         representation(tags="character",groups="list"),
         validity=function(object) {
           if(length(object@tags)!=length(object@groups)) {
             stop("groups must be named by their tag")
           }
         })
setClass("tags",
         representation(tags="character",snps="character"),
         validity=function(object) {
           if(length(object@tags)!=length(object@snps))
             stop("tags must be in tags vector, tagging themselves")
         })

