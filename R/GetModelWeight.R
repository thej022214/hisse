## Compute model weights for GeoHiSSE, HiSSE, or MiSSE models given a list of models.

# BCO: DEPRECATED -- duplicates functionality in GetAICWeights, while not allowing for the possibility of AICc. Do not use. Kept here for historical reference for now

GetModelWeight <- function(...){
    models <- list(...)
    ## Check if we have a single list with all the models or the arguments are a list.
    if( !inherits(models[[1]], what = c("geohisse.fit", "hisse.fit", "misse.fit")) ){
        ## This is a list of models. Need to extract.
        models <- models[[1]]
    }
    ## Check if elements of the list are HiSSE, GeoHiSSE, or MiSSE models.
    mod.class.geohisse <- sapply(models, function(x) inherits(x, what = "geohisse.fit"))
    mod.class.hisse <- sapply(models, function(x) inherits(x, what = "hisse.fit"))
    mod.class.misse <- sapply(models, function(x) inherits(x, what = "misse.fit"))

    if( all(mod.class.geohisse) & all(mod.class.hisse) & all(mod.class.misse) ){
        ## Strange! Break.
        stop( "list of models need to be only HiSSE, or GeoHiSSE, or MiSSE fits." )
    }
    if( !all(mod.class.geohisse) & !all(mod.class.hisse & !all(mod.class.misse)) ){
        ## Strange! Break.
        stop( "list of models need to be only HiSSE, or GeoHiSSE, or MiSSE fits." )
    }
    mod.names <- names( models )
    mod.AIC <- sapply(models, function(x) x$AIC )
    delta <- mod.AIC - min( mod.AIC )
    aicw <- exp( -0.5 * delta) / sum( exp( -0.5 * delta) )
    names( aicw ) <- mod.names
    return( aicw )
}
