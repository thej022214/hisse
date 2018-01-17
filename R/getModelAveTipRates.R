GetModelAveTipRates <- function(x){
    hisse.results <- x

    if( !inherits(hisse.results, what = c("list", "hisse.states", "hisse.geosse.states")) ) stop("x needs to be a list of model reconstructions or a single model reconstruction object of class 'hisse.states' or 'hisse.geosse.states'.")

    ## If hisse.results is a list of model reconstructions, then test if they have $aic. Return error message otherwise.
    if(class(hisse.results) == "list"){
        empty.aic <- sapply(hisse.results, function(x) !is.null(x$aic) )
        if( as.logical( sum( empty.aic ) ) ) stop("All elements of the list need to have a '$aic' element.")
        model.class <- sapply(hisse.results, function(x) !inherits(x, what = c("hisse.states", "hisse.geosse.states")) )
        if( as.logical( sum( model.class ) ) ) stop("x needs to be a list of model reconstruction with class 'hisse.states' or 'hisse.geosse.states' ")
    }

    if( inherits(hisse.results, what = c("hisse.states", "hisse.geosse.states")) ){ # we have to make a list so we can run this generally
        if(is.null(hisse.results$aic)){
            ## If a user forgot to include the aic, then we add a random value in for them
            hisse.results$aic = 42
        }
        tmp.list <- list()
        tmp.list[[1]] <- hisse.results
        hisse.results <- tmp.list
    }
    
    rates.tips.turnover <- ConvertManyToRate(hisse.results, rate.param="turnover", "tip.mat")
    rates.tips.net.div <- ConvertManyToRate(hisse.results, rate.param="net.div", "tip.mat")
    rates.tips.speciation <- ConvertManyToRate(hisse.results, rate.param="speciation", "tip.mat")
    rates.tips.extinct.fraction <- ConvertManyToRate(hisse.results, rate.param="extinction.fraction", "tip.mat")
    rates.tips.extinction <- ConvertManyToRate(hisse.results, rate.param="extinction", "tip.mat")
    
    ## Objects will always be of list class here.
    if(class(hisse.results[[1]])=="hisse.states"){
        states.tips <- ConvertManyToBinaryState(hisse.results, "tip.mat")
    }
    if(class(hisse.results[[1]])=="hisse.geosse.states"){
        states.tips <- ConvertManyToMultiState(hisse.results, "tip.mat")
    }

    final.df <- data.frame(taxon=hisse.results[[1]]$phy$tip.label, state=states.tips, turnover=rates.tips.turnover, net.div=rates.tips.net.div, speciation=rates.tips.speciation, extinct.frac=rates.tips.extinct.fraction, extinction=rates.tips.extinction)
    return(final.df)
}

GetModelAveNodeRates <- function(x){
    hisse.results <- x

    if( !inherits(hisse.results, what = c("list", "hisse.states", "hisse.geosse.states")) ) stop("x needs to be a list of model reconstructions or a single model reconstruction object of class 'hisse.states' or 'hisse.geosse.states'.")

    ## If hisse.results is a list of model reconstructions, then test if they have $aic. Return error message otherwise.
    if(class(hisse.results) == "list"){
        empty.aic <- sapply(hisse.results, function(x) !is.null(x$aic) )
        if( as.logical( sum( empty.aic ) ) ) stop("All elements of the list need to have a '$aic' element.")
        model.class <- sapply(hisse.results, function(x) !inherits(x, what = c("hisse.states", "hisse.geosse.states")) )
        if( as.logical( sum( model.class ) ) ) stop("x needs to be a list of model reconstruction with class 'hisse.states' or 'hisse.geosse.states' ")
    }

    if( inherits(hisse.results, what = c("hisse.states", "hisse.geosse.states")) ){ # we have to make a list so we can run this generally
        if(is.null(hisse.results$aic)){
            ## If a user forgot to include the aic, then we add a random value in for them
            hisse.results$aic = 42
        }
        tmp.list <- list()
        tmp.list[[1]] <- hisse.results
        hisse.results <- tmp.list
    }
    
    rates.tips.turnover <- ConvertManyToRate(hisse.results, rate.param="turnover", "node.mat")
    rates.tips.net.div <- ConvertManyToRate(hisse.results, rate.param="net.div", "node.mat")
    rates.tips.speciation <- ConvertManyToRate(hisse.results, rate.param="speciation", "node.mat")
    rates.tips.extinct.fraction <- ConvertManyToRate(hisse.results, rate.param="extinction.fraction", "node.mat")
    rates.tips.extinction <- ConvertManyToRate(hisse.results, rate.param="extinction", "node.mat")

    ## Objects will always be of list class here.
    if(class(hisse.results[[1]])=="hisse.states"){
        states.internal <- ConvertManyToBinaryState(hisse.results, "node.mat")
    }
    if(class(hisse.results[[1]])=="hisse.geosse.states"){
        states.internal <- ConvertManyToMultiState(hisse.results, "node.mat")
    }
    
    final.df <- data.frame(id=hisse.results[[1]]$node.mat[,1], state=states.internal, turnover=rates.tips.turnover, net.div=rates.tips.net.div, speciation=rates.tips.speciation, extinct.frac=rates.tips.extinct.fraction, extinction=rates.tips.extinction)
    return(final.df)
}
