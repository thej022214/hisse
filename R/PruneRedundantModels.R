PruneRedundantModels <- function(..., precision=1e-5) {
	models <- list(...)
    ## Check if we have a single list with all the models or the arguments are a list.
    if( !inherits(models[[1]], what = c("geohisse.fit", "hisse.fit", "misse.fit", "muhisse.fit")) ){
        ## This is a list of models. Need to extract.
        models <- models[[1]]
    }
    ## Check if elements of the list are HiSSE, GeoHiSSE, or MiSSE models.
    mod.class.geohisse <- sapply(models, function(x) inherits(x, what = "geohisse.fit"))
    mod.class.hisse <- sapply(models, function(x) inherits(x, what = "hisse.fit"))
    mod.class.misse <- sapply(models, function(x) inherits(x, what = "misse.fit"))
    mod.class.muhisse <- sapply(models, function(x) inherits(x, what = "muhisse.fit"))

    if( all(mod.class.geohisse) & all(mod.class.hisse) & all(mod.class.misse) & all(mod.class.muhisse)){
        ## Strange! Break.
        stop( "list of models need to be only HiSSE, MuHiSSE, GeoHiSSE, or MiSSE fits." )
    }
    if( !all(mod.class.geohisse) & !all(mod.class.hisse) & !all(mod.class.misse) & !all(mod.class.muhisse)){
        ## Strange! Break.
        stop( "list of models need to be only HiSSE, MuHiSSE, GeoHiSSE, or MiSSE fits." )
    }
	mod.nparameters <- simplify2array(lapply(lapply(models, "[[", "starting.vals"),length))
	models <- models[order(mod.nparameters, decreasing=FALSE)]
	mod.loglik <- simplify2array(lapply(models, "[[", "loglik"))
	models_to_delete <- c()
	isTrueAllEqual <- function(...) {
		return(isTRUE(all.equal(...)))
	}
	if(length(models)>1) {
		for(i in 2:(length(models))) {
			if(any(sapply(mod.loglik[1:(i-1)], isTrueAllEqual, mod.loglik[i], tolerance=precision))) {
				models_to_delete <- c(models_to_delete, i)
			}
		}
	}
	if(length(models_to_delete)>0) {
		models <- models[-models_to_delete]
	}
    return( models )
}
