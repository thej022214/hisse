## Compute model weights for GeoHiSSE and HiSSE models given a list of models.

GetModelWeight <- function(...){
    models <- list(...)
    ## Check if we have a single list with all the models or the arguments are a list.
    if( !inherits(models[[1]], what = c("geohisse.fit", "hisse.fit")) ){
        ## This is a list of models. Need to extract.
        models <- models[[1]]
    }
    ## Check if elements of the list are HiSSE or GeoHiSSE models.
    mod.class.geohisse <- sapply(models, function(x) inherits(x, what = "geohisse.fit"))
    mod.class.hisse <- sapply(models, function(x) inherits(x, what = "hisse.fit"))
    if( all(mod.class.geohisse) & all(mod.class.hisse) ){
        ## Strange! Break.
        stop( "list of models need to be only HiSSE OR GeoHiSSE fits." )
    }
    if( !all(mod.class.geohisse) & !all(mod.class.hisse) ){
        ## Strange! Break.
        stop( "list of models need to be only HiSSE OR GeoHiSSE fits." )
    }    
    mod.names <- names( models )
    mod.AIC <- sapply(models, function(x) x$AIC )
    delta <- mod.AIC - min( mod.AIC )
    aicw <- exp( -0.5 * delta) / sum( exp( -0.5 * delta) )
    names( aicw ) <- mod.names
    return( aicw )
}

