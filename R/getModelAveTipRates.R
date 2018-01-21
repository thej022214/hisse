GetModelAveTipRates <- function(x, AIC.weights=NULL){
    hisse.results <- x

    if( !inherits(hisse.results, what = c("list", "hisse.states", "hisse.geosse.states")) ) stop("x needs to be a list of model reconstructions or a single model reconstruction object of class 'hisse.states' or 'hisse.geosse.states'.")

    ## If hisse.results is a list of model reconstructions, then test if they have $aic. Return error message otherwise.
    ## There is no need for the $aic element if AIC.weigths argument is provided.
    ## AIC.weights is a substitute for AIC. See help page.
    if(class(hisse.results) == "list"){
        if( is.null( AIC.weights ) ){
            empty.aic <- sapply(hisse.results, function(x) !is.null(x$aic) )
            if( sum( empty.aic ) != length(hisse.results) ) stop("All elements of the list need to have a '$aic' element.")
            model.class <- sapply(hisse.results, function(x) !inherits(x, what = c("hisse.states", "hisse.geosse.states")) )
            if( as.logical( sum( model.class ) ) ) stop("x needs to be a list of model reconstruction with class 'hisse.states' or 'hisse.geosse.states' ")
        } else{
            if( !length(hisse.results) == length(AIC.weights) ){
                stop( " AIC.weights needs to be NULL or a numeric vector with length equal to the number of models in 'x'. " )
            }
        }
    }

    if( inherits(hisse.results, what = c("hisse.states", "hisse.geosse.states")) ){ # we have to make a list so we can run this generally
        if( !is.null( AIC.weights ) ){
            stop( "If using a single model, then 'AIC.weights' needs to be NULL. " )
        }
        
        if(is.null(hisse.results$aic)){
            ## If a user forgot to include the aic, then we add a random value in for them
            hisse.results$aic = 42
        }
        tmp.list <- list()
        tmp.list[[1]] <- hisse.results
        hisse.results <- tmp.list
    }
    
    rates.tips.turnover <- ConvertManyToRate(hisse.results, rate.param="turnover", which.element="tip.mat", AIC.weights=AIC.weights)
    rates.tips.net.div <- ConvertManyToRate(hisse.results, rate.param="net.div", which.element="tip.mat", AIC.weights=AIC.weights)
    rates.tips.speciation <- ConvertManyToRate(hisse.results, rate.param="speciation", which.element="tip.mat", AIC.weights=AIC.weights)
    rates.tips.extinct.fraction <- ConvertManyToRate(hisse.results, rate.param="extinction.fraction", which.element="tip.mat", AIC.weights=AIC.weights)
    rates.tips.extinction <- ConvertManyToRate(hisse.results, rate.param="extinction", which.element="tip.mat", AIC.weights=AIC.weights)
    
    ## Objects will always be of list class here.
    if(class(hisse.results[[1]])=="hisse.states"){
        states.tips <- ConvertManyToBinaryState(hisse.results, which.element="tip.mat", AIC.weights=AIC.weights)
    }
    if(class(hisse.results[[1]])=="hisse.geosse.states"){
        states.tips <- ConvertManyToMultiState(hisse.results, which.element="tip.mat", AIC.weights=AIC.weights)
    }

    final.df <- data.frame(taxon=hisse.results[[1]]$phy$tip.label, state=states.tips, turnover=rates.tips.turnover, net.div=rates.tips.net.div, speciation=rates.tips.speciation, extinct.frac=rates.tips.extinct.fraction, extinction=rates.tips.extinction)
    return(final.df)
}

GetModelAveNodeRates <- function(x, AIC.weights=NULL){
    hisse.results <- x

    if( !inherits(hisse.results, what = c("list", "hisse.states", "hisse.geosse.states")) ) stop("x needs to be a list of model reconstructions or a single model reconstruction object of class 'hisse.states' or 'hisse.geosse.states'.")

    ## If hisse.results is a list of model reconstructions, then test if they have $aic. Return error message otherwise.
    ## There is no need for the $aic element if AIC.weigths argument is provided.
    ## AIC.weights is a substitute for AIC. See help page.
    if(class(hisse.results) == "list"){
        if( is.null( AIC.weights ) ){
            empty.aic <- sapply(hisse.results, function(x) !is.null(x$aic) )
            if( sum( empty.aic ) != length(hisse.results) ) stop("All elements of the list need to have a '$aic' element.")
            model.class <- sapply(hisse.results, function(x) !inherits(x, what = c("hisse.states", "hisse.geosse.states")) )
            if( as.logical( sum( model.class ) ) ) stop("x needs to be a list of model reconstruction with class 'hisse.states' or 'hisse.geosse.states' ")
        } else{
            if( !length(hisse.results) == length(AIC.weights) ){
                stop( " AIC.weights needs to be NULL or a numeric vector with length equal to the number of models in 'x'. " )
            }
        }
    }

    if( inherits(hisse.results, what = c("hisse.states", "hisse.geosse.states")) ){ # we have to make a list so we can run this generally
        if( !is.null( AIC.weights ) ){
            stop( "If using a single model, then 'AIC.weights' needs to be NULL. " )
        }
        
        if(is.null(hisse.results$aic)){
            ## If a user forgot to include the aic, then we add a random value in for them
            hisse.results$aic = 42
        }
        tmp.list <- list()
        tmp.list[[1]] <- hisse.results
        hisse.results <- tmp.list
    }
    
    rates.tips.turnover <- ConvertManyToRate(hisse.results, rate.param="turnover", which.element="node.mat", AIC.weights=AIC.weights)
    rates.tips.net.div <- ConvertManyToRate(hisse.results, rate.param="net.div", which.element="node.mat", AIC.weights=AIC.weights)
    rates.tips.speciation <- ConvertManyToRate(hisse.results, rate.param="speciation", which.element="node.mat", AIC.weights=AIC.weights)
    rates.tips.extinct.fraction <- ConvertManyToRate(hisse.results, rate.param="extinction.fraction", which.element="node.mat", AIC.weights=AIC.weights)
    rates.tips.extinction <- ConvertManyToRate(hisse.results, rate.param="extinction", which.element="node.mat", AIC.weights=AIC.weights)

    ## Objects will always be of list class here.
    if(class(hisse.results[[1]])=="hisse.states"){
        states.internal <- ConvertManyToBinaryState(hisse.results, which.element="node.mat", AIC.weights=AIC.weights)
    }
    if(class(hisse.results[[1]])=="hisse.geosse.states"){
        states.internal <- ConvertManyToMultiState(hisse.results, which.element="node.mat", AIC.weights=AIC.weights)
    }
    
    final.df <- data.frame(id=hisse.results[[1]]$node.mat[,1], state=states.internal, turnover=rates.tips.turnover, net.div=rates.tips.net.div, speciation=rates.tips.speciation, extinct.frac=rates.tips.extinct.fraction, extinction=rates.tips.extinction)
    return(final.df)
}
