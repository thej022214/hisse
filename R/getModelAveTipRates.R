######################################################################################################################################
######################################################################################################################################
### Model Average Rates -- calculates a model-averaged rates for each observed state under a set of parameters
######################################################################################################################################
######################################################################################################################################

GetModelAveRates <- function(x, AIC.weights=NULL, type=c("tips", "nodes", "both"), bound.par.matrix=cbind(c(0,-1000,0,0,0),c(1000,1000,1000,3,1000))){
    hisse.results <- x

    if (!type %in% c("tips", "nodes", "both")) stop("Argument 'type' needs to be one of 'tips', 'nodes', or 'both'.")

    ## Check if the format of bound par matrix is correct.
    ## This needs to be a numeric matrix with two columns and 5 rows.
    ## First column are minimums and second colum are maximums.
    if( !inherits(bound.par.matrix, what = c("matrix","data.frame")) ){
        stop(" 'bound.par.matrix' needs to be a matrix or data.frame ")
    }
    if( !ncol(bound.par.matrix) == 2 | !nrow(bound.par.matrix) == 5 ){
        stop( " 'bound.par.matrix' needs to be a matrix with 5 rows and 2 columns. See help.")
    }
    ## Beware that the test below might misbehave if the lengths are different.
    ## But since we tested the dimensions, then this seems fine.
    if( !all( bound.par.matrix[,1] < bound.par.matrix[,2] ) ){
        stop("Min bounds larger than max! First column of 'bound.par.matrix' are the minimum bounds and second column are the maximum.")
    }
    
    ## Create flag to compute the model average for the both nodes and tips or just one of them.
    if( type == "both" ){
        to.mod.ave <- c("tips", "nodes")
    } else{
        to.mod.ave <- type
    }    

    if( !inherits(hisse.results, what = c("list", "hisse.states", "hisse.geosse.states", "muhisse.states", "misse.states")) ) stop("x needs to be a list of model reconstructions or a single model reconstruction object of class 'hisse.states', 'hisse.geosse.states', 'muhisse.states', or 'misse.states'.")

    ## If hisse.results is a list of model reconstructions, then test if they have $aic. Return error message otherwise.
    ## There is no need for the $aic element if AIC.weigths argument is provided.
    ## AIC.weights is a substitute for AIC. See help page.
    if(class(hisse.results) == "list"){
        if( is.null( AIC.weights ) ){
            empty.aic <- sapply(hisse.results, function(x) !is.null(x$aic) )
            if( sum( empty.aic ) != length(hisse.results) ) stop("All elements of the list need to have a '$aic' element.")
            model.class <- sapply(hisse.results, function(x) !inherits(x, what = c("hisse.states", "hisse.geosse.states", "muhisse.states")) )
            if( as.logical( sum( model.class ) ) ) stop("x needs to be a list of model reconstruction with class 'hisse.states' or 'hisse.geosse.states' ")
        } else{
            if( !length(hisse.results) == length(AIC.weights) ){
                stop( " AIC.weights needs to be NULL or a numeric vector with length equal to the number of models in 'x'. " )
            }
        }
    }

    if( inherits(hisse.results, what = c("hisse.states", "hisse.geosse.states", "muhisse.states", "misse.states")) ){ # we have to make a list so we can run this generally
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

    ## This will compute the AIC weights just once now. ;)
    if( is.null(AIC.weights) ){
        AIC.weights <- GetAICWeights(hisse.results)
    }
    
    ## ##################################################################################################################
    ## Compute the model average for the nodes:
    if( "tips" %in% to.mod.ave ){        
        rates.tips.turnover <- ConvertManyToRate_ModelAve(hisse.results, rate.param="turnover", which.element="tip.mat")
        rates.tips.net.div <- ConvertManyToRate_ModelAve(hisse.results, rate.param="net.div", which.element="tip.mat")
        rates.tips.speciation <- ConvertManyToRate_ModelAve(hisse.results, rate.param="speciation", which.element="tip.mat")
        rates.tips.extinct.fraction <- ConvertManyToRate_ModelAve(hisse.results, rate.param="extinction.fraction", which.element="tip.mat")
        rates.tips.extinction <- ConvertManyToRate_ModelAve(hisse.results, rate.param="extinction", which.element="tip.mat")

        ## Here 'averaged.rates' is a list with the vector of results.
        averaged.tip.rates <- CheckReconBounds(x=list(rates.tips.turnover, rates.tips.net.div, rates.tips.speciation, rates.tips.extinct.fraction, rates.tips.extinction), n.models = length(hisse.results), AIC.weights, bound.par.matrix)
        
        ## Objects will always be of list class here.
        if(class(hisse.results[[1]])=="hisse.states"){
            states.tips <- ConvertManyToBinaryState(hisse.results, which.element="tip.mat", AIC.weights=AIC.weights)
        }
     
        if(class(hisse.results[[1]])=="hisse.geosse.states"){
            states.tips <- ConvertManyToMultiState(hisse.results, which.element="tip.mat", AIC.weights=AIC.weights)
        }

        if(class(hisse.results[[1]])=="muhisse.states"){
            states.tips <- ConvertManyToMultiState(hisse.results, which.element="tip.mat", AIC.weights=AIC.weights)
        }

        if(class(hisse.results[[1]])=="misse.states"){
            states.tips <- rep(0, dim(hisse.results[[1]][["tip.mat"]])[1])
        }


        ## Note that the columns of 'averaged.tip.rates' need to be in the correct order here. They should be.
        final.df.tips <- data.frame(taxon=hisse.results[[1]]$phy$tip.label, state=states.tips, turnover=averaged.tip.rates[[1]], net.div=averaged.tip.rates[[2]], speciation=averaged.tip.rates[[3]], extinct.frac=averaged.tip.rates[[4]], extinction=averaged.tip.rates[[5]])
        ## Check if the function need to return only the reconstruction for the tips.
        if( type == "tips" ) return(final.df.tips)
    }
    
    ## ##################################################################################################################
    ## Compute the model average for the nodes:
    if( "nodes" %in% to.mod.ave ){
        rates.nodes.turnover <- ConvertManyToRate_ModelAve(hisse.results, rate.param="turnover", which.element="node.mat")
        rates.nodes.net.div <- ConvertManyToRate_ModelAve(hisse.results, rate.param="net.div", which.element="node.mat")
        rates.nodes.speciation <- ConvertManyToRate_ModelAve(hisse.results, rate.param="speciation", which.element="node.mat")
        rates.nodes.extinct.fraction <- ConvertManyToRate_ModelAve(hisse.results, rate.param="extinction.fraction", which.element="node.mat")
        rates.nodes.extinction <- ConvertManyToRate_ModelAve(hisse.results, rate.param="extinction", which.element="node.mat")

        ## Here 'averaged.node.rates' is a list with the vector of results.
        averaged.node.rates <- CheckReconBounds(x=list(rates.nodes.turnover, rates.nodes.net.div, rates.nodes.speciation, rates.nodes.extinct.fraction, rates.nodes.extinction), n.models = length(hisse.results), AIC.weights, bound.par.matrix)
        
        ## Objects will always be of list class here.
        if(class(hisse.results[[1]])=="hisse.states"){
            states.internal <- ConvertManyToBinaryState(hisse.results, which.element="node.mat", AIC.weights=AIC.weights)
        }

        if(class(hisse.results[[1]])=="hisse.geosse.states"){
            states.internal <- ConvertManyToMultiState(hisse.results, which.element="node.mat", AIC.weights=AIC.weights)
        }

        if(class(hisse.results[[1]])=="muhisse.states"){
            states.internal <- ConvertManyToMultiState(hisse.results, which.element="node.mat", AIC.weights=AIC.weights)
        }

        if(class(hisse.results[[1]])=="misse.states"){
            states.internal <- rep(0, dim(hisse.results[[1]][["node.mat"]])[1])
        }

        final.df.nodes <- data.frame(id=hisse.results[[1]]$node.mat[,1], state=states.internal, turnover=averaged.node.rates[[1]], net.div=averaged.node.rates[[2]], speciation=averaged.node.rates[[3]], extinct.frac=averaged.node.rates[[4]], extinction=averaged.node.rates[[5]])
        ## Check if the function need to return only the reconstruction for the nodes:
        if( type == "nodes" ) return(final.df.nodes)
    }

    ## If both node and tips need to be returned:
    if( type == "both" ){
        final.both <- list(tips = final.df.tips, nodes = final.df.nodes)
        return( final.both )
    }
}
