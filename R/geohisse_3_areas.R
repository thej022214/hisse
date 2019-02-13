
## ##############################################################################################################
## GeoHiSSE -- Expanded set of GeoSSE models for examining diversification in relation to geographic range evolution
## This function allows for 3 endemic ranges on the model contrasted with the 2 hidden areas of the original model.
## ##############################################################################################################

## Implementation remark:
## Here we are using functions to call different types of models. To reduce unnecessary opperations inside the likelihood evaluation this code make primary checks for the size of the model and restrict the ODE's only to the potential scenarios under a particular number of hidden states.
## This reduces the number of ODE numerical integrations per branch and can significatly improve time to fit the model.

GeoHiSSE_Plus <- function(phy, data, f=c(1,1,1,1,1,1), speciation=c(1,2,3,4,5,6), extirpation=c(1,2,3), hidden.areas=FALSE, trans.rate=NULL, assume.cladogenetic=TRUE, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, sann=FALSE, sann.its=10000, bounded.search=TRUE, max.tol=.Machine$double.eps^.50, mag.san.start=0.5, starting.vals=NULL, speciation.upper=1000, extirpation.upper=1000, trans.upper=100, ode.eps=0){

    ## This is the higher-level function. Here we will get the data, perform tests and pass to more specific functions.

    ## ########################################################
    ## Data testing and parameter format block
    ## ########################################################

    ## The series of checks below test for more common mismatches.
    ## The length of the speciation vector need to be a multiple of 6 and match the dimention of the trans.rate matrix.
    trans.dim <- dim(trans.rate)[1]
    if(hidden.areas == TRUE & any( length(speciation) < c(6,12) ) ){
        stop("You chose a hidden state but this is not reflected in the speciation vector.")
    }
    if(hidden.areas == TRUE & trans.dim < 6){
        stop("You chose a hidden state but this is not reflected in the transition matrix")
    }
    if(hidden.areas == FALSE & length(speciation) > 6){
        stop("You chose no hidden state but this is not reflected in the speciation vector.")
    }
    if(hidden.areas == FALSE & length(extirpation) > 3){
        stop("You chose no hidden state but this is not reflected in the extinction vector.")
    }
    if(hidden.areas == FALSE & trans.dim > 6){
        stop("You chose no hidden state but this is not reflected in the transition matrix")
    }
    ## Make a more detailed check if the vector lengths match:
    if( trans.dim != length( speciation ) ){
        stop("Number of speciation parameters and dimension of trans.rate do not match")
    }
    if( (trans.dim/2) != length( extirpation ) ){
        stop("Number of extirpation parameters and dimension of trans.rate do not match")
    }
    
    ## Temporary fix for the current BUG:
    if( !is.null(phy$node.label) ) phy$node.label <- NULL

    ## Check input vectors:
    if( !is.numeric(speciation) | !is.numeric(extirpation) ){
        stop("the speciation and extirpation vectors need to be numeric and integers.")
    }
    if( length(speciation) %% 6 != 0 ){
        stop("the speciation vector need to be a multiple of 6.")
    }
    if( length(extirpation) %% 3 != 0 ){
        stop("the extirpation vector need to be a multiple of 3.")
    }

    ## Also, need to make sure that the parameter indexes on the trans.rate matrix have some important things:
    rr.trans.rate <- range(trans.rate[trans.rate > 0 & !is.na(trans.rate)])
    if( !all( seq(from = rr.trans.rate[1], to = rr.trans.rate[2]) %in% unique( trans.rate[trans.rate > 0 & !is.na(trans.rate)] ) ) ){
        stop("The trans.rate matrix needs sequential parameter indexes. The matrix needs all integers between the min and max numbers.")
    }
    
    if(!is.null(root.p)) {
        ## The vector of the root.p need to be as long as the speciation vector.
        if( length( root.p ) < 6 ){
            stop("length of root.p vector need to be at least 6.")
        }
        if( length(speciation) %% 6 != 0 ){
            stop("the root.p vector need to be a multiple of 6.")
        }
        if( length( root.p ) > length( speciation ) ){
            stop("length of root.p vector cannot exceed length of speciation vector. Check parameters!")
        }
        if( length( root.p ) != length( speciation ) ){
            if( length( root.p ) == 6 ){
                warning("For hidden states, you need to specify the root.p for all hidden states. We have adjusted it so that there's equal probability among all hidden states.")
                nrates <- length( speciation ) / 6 ## Number of hidden states.
                root.p <- rep(root.p, times = nrates)
                root.p <- root.p / sum(root.p)
            } else{
                stop("User provided root probabilities need to be a vector of length 6 or 6 * # of hiddens rates.")
            }
        } else{
            ## All good:
            root.p <- root.p / sum(root.p)
        }
    }
    
    if(!root.type == "madfitz" & !root.type == "herr_als"){
        stop("Check that you specified a proper root.type option. Options are 'madfitz' or 'herr_als'. See help for more details.", call.=FALSE)
    }

    if(is.null(trans.rate)){
        ## We will use the same function to create the transition matrix, but with a separated option.
        stop("Rate matrix needed. See TransMatMakerGeoHiSSE() to create one.")
    }
    
    if(assume.cladogenetic == FALSE){
        stop("Sorry. 'assume.cladogenetic == FALSE' option is not implemented in this function." )
    }
    
    ## Return error message if the data is not in the correct format.
    if( !inherits(data, what = c("matrix","data.frame")) ){
        stop("'data' needs to be a matrix or data.frame with 3 columns. See help.")
    }
    if( !ncol( data ) == 3 ){
        stop("'data' needs to be a matrix or data.frame with 3 columns. See help.")
    }
    if( !all( apply(data, 2, is.numeric) ) ){
        stop("Wrong data input. Check help page.")
    }    
    if( any( data > 1 ) ){
        stop("Wrong data input. Check help page.")
    }
    check.area.presence <- rowSums(data)
    if( any( check.area.presence > 2 ) ){
        stop("Widespread areas can only be comprised of two endemic areas. Species occurring on all three areas are not allowed at the moment. Please contact the authors to ask about extensions.")
    }

    ## ########################################################
    ## Transform the data from matrix to vector format.
    ## ########################################################
    
    data.code <- apply(data, 1, function(x) paste0(x, collapse=""))
    
    ## The order of the states here, from 1 to 6, need to match the order used in the C code and the rest of
    ##     the pipeline. Note that this order is different from the order of the states at the transition matrix for
    ##     argument of the function.
    
    states.new <- sapply(data.code, function(x) switch(x, "100"=1, "010"=2, "001"=3, "110"=4, "011"=5, "101"=6) )
    
    ## data.new <- data.frame(unname(states.new), unname(states.new), row.names=names(states.new))
    ## data.new <- data.new[phy$tip.label,]
    
    ## ########################################################
    ## Block to filter models and call more specific functions.
    ## ########################################################

    ## The simpler model.
    if( hidden.areas == FALSE & assume.cladogenetic == TRUE ){
        ## Complete the pipeline and return the answer:
        ## Function will return a list with all model estimates.
        fit.out <- geohisse_3_one_rate(phy=phy, data=states.new, f=f, speciation=speciation, extirpation=extirpation
                                     , trans.rate=trans.rate, condition.on.survival=condition.on.survival
                                     , root.type=root.type, root.p=root.p, sann=sann, sann.its=sann.its
                                     , bounded.search=bounded.search, max.tol=max.tol
                                     , mag.san.start=mag.san.start, starting.vals=starting.vals
                                     , speciation.upper=speciation.upper, extirpation.upper=extirpation.upper
                                     , trans.upper=trans.upper, ode.eps=ode.eps)
        obj <- list(loglik = fit.out$loglik, AIC = fit.out$AIC, AICc = fit.out$AICc, solution= fit.out$solution
                  , index.par= fit.out$index.par, f=fit.out$f, hidden.areas=FALSE
                  , assume.cladogenetic=TRUE, condition.on.survival=condition.on.survival
                  , root.type=root.type, root.p=root.p, phy=phy, data=data, trans.matrix=trans.rate
                  , max.tol=max.tol, starting.vals=fit.out$starting.vals, upper.bounds=fit.out$upper.bounds
                  , lower.bounds=fit.out$lower.bounds, ode.eps=ode.eps, n.hidden.rates=1)
        class(obj) <- append(class(obj), "geohisse_plus.fit")
        return(obj)
    }
    
    ## Model with 2 hidden rates. The transition matrix will have 12 rows.
    ## At this point it does not matter if the model is a full or a null model. Parameters will all depend on
    ##    the transition matrix and on the speciation and extirpation vectors.
    if( hidden.areas == TRUE & assume.cladogenetic == TRUE & trans.dim == 12 ){
        fit.out <- geohisse_3_two_rate(phy=phy, data=states.new, f=f, speciation=speciation, extirpation=extirpation
                                     , trans.rate=trans.rate, condition.on.survival=condition.on.survival
                                     , root.type=root.type, root.p=root.p, sann=sann, sann.its=sann.its
                                     , bounded.search=bounded.search, max.tol=max.tol
                                     , mag.san.start=mag.san.start, starting.vals=starting.vals
                                     , speciation.upper=speciation.upper, extirpation.upper=extirpation.upper
                                     , trans.upper=trans.upper, ode.eps=ode.eps)
        obj <- list(loglik = fit.out$loglik, AIC = fit.out$AIC, AICc = fit.out$AICc, solution= fit.out$solution
                  , index.par= fit.out$index.par, f=fit.out$f, hidden.areas=FALSE
                  , assume.cladogenetic=TRUE, condition.on.survival=condition.on.survival
                  , root.type=root.type, root.p=root.p, phy=phy, data=data, trans.matrix=trans.rate
                  , max.tol=max.tol, starting.vals=fit.out$starting.vals, upper.bounds=fit.out$upper.bounds
                  , lower.bounds=fit.out$lower.bounds, ode.eps=ode.eps, n.hidden.rates=2)
        class(obj) <- append(class(obj), "geohisse_plus.fit")
        return(obj)
    }

    if( hidden.areas == TRUE & assume.cladogenetic == TRUE & trans.dim > 12 ){
        ## This will be the fit of the model with a general number of hidden states from 3 to 12.
        ## Need to check if the model is not trying to fit more hidden rates than possible here.
        if( trans.dim > 6*12 | length(speciation) > 6*12 | length(extirpation) > 3*12 ){
            stop("The maximum of hidden rates is 12!")
        }
        n.hidden.rates <- length(speciation) / 6
        fit.out <- geohisse_3_multi_rate(phy=phy, data=states.new, f=f, speciation=speciation, extirpation=extirpation
                                     , trans.rate=trans.rate, condition.on.survival=condition.on.survival
                                     , root.type=root.type, root.p=root.p, sann=sann, sann.its=sann.its
                                     , bounded.search=bounded.search, max.tol=max.tol
                                     , mag.san.start=mag.san.start, starting.vals=starting.vals
                                     , speciation.upper=speciation.upper, extirpation.upper=extirpation.upper
                                     , trans.upper=trans.upper, ode.eps=ode.eps)
        obj <- list(loglik = fit.out$loglik, AIC = fit.out$AIC, AICc = fit.out$AICc, solution= fit.out$solution
                  , index.par= fit.out$index.par, f=fit.out$f, hidden.areas=FALSE
                  , assume.cladogenetic=TRUE, condition.on.survival=condition.on.survival
                  , root.type=root.type, root.p=root.p, phy=phy, data=data, trans.matrix=trans.rate
                  , max.tol=max.tol, starting.vals=fit.out$starting.vals, upper.bounds=fit.out$upper.bounds
                  , lower.bounds=fit.out$lower.bounds, ode.eps=ode.eps, n.hidden.rates=n.hidden.rates)
        class(obj) <- append(class(obj), "geohisse_plus.fit")
        return(obj)
    }

}

## ##############################################################################################################
## Print function for the output class:
## ##############################################################################################################

print.geohisse_plus.fit <- function(x,...){
    ## Function to print a "geohisse.fit" object.
    ## Keep only the parameters estimated:
    keep.par <- x$index.par > 0
    par.list <- x$solution[ keep.par ]
    ntips <- Ntip( x$phy )
    output <- c(x$loglik, x$AIC, x$AICc, ntips, x$n.hidden.rates)
    names(output) <- c("-lnL", "AIC", "AICc", "n.taxa", "n.hidden.rates")
    cat("\n")
    cat("Fit \n")
    print(output)
    cat("\n")
    cat("Model parameters: \n")
    cat("\n")
    print(par.list)
    cat("\n")
}

## ##############################################################################################################
## Fitting function for no hidden states:
## ##############################################################################################################

geohisse_3_one_rate <- function(phy, data, f, speciation, extirpation, trans.rate, condition.on.survival, root.type, root.p, sann, sann.its, bounded.search, max.tol, mag.san.start, starting.vals, speciation.upper, extirpation.upper, trans.upper, ode.eps){

    ## Note that some parameters are not present at this function call.
    ## We are assuming no hidden states and that the model is cladogenetic.

    ## Diversification parameters for the model (ext'1 means local extinction):
    ## s1, s2, s3, ext1, ext2, ext3, (ext'1, ext'2, ext'3), s12, s13, s23
    ## The transition parameters are set using the transition matrix.

    ## This step works even if the 
    extirpation.tmp <- extirpation
    extirpation.tmp[which(extirpation.tmp > 0)] <- extirpation.tmp[which(extirpation.tmp > 0)] + max(speciation)
    pars.tmp <- c(speciation, extirpation.tmp)
    
    ## Here we need a vector with the transition indexes from the 'trans.rate' parameter in the same order as the vector we call in the C code.
    ## This matrix also informs the model about jumps between endemic regions and separated rates between extirpation and local extinction.
    ## We can use the values of this transition matrix to guess if the model is separating local extinction from extirpation. This information will be important to set appropriate starting values and bound values.
    trans.tmp <- c(trans.rate["(1)", "(2)"], trans.rate["(2)", "(1)"], trans.rate["(2)", "(3)"]
                 , trans.rate["(3)", "(2)"], trans.rate["(1)", "(3)"], trans.rate["(3)", "(1)"]
                 , trans.rate["(1)", "(12)"], trans.rate["(1)", "(13)"], trans.rate["(2)", "(12)"]
                 , trans.rate["(2)", "(23)"], trans.rate["(3)", "(13)"], trans.rate["(3)", "(23)"]
                 , trans.rate["(12)", "(1)"], trans.rate["(13)", "(1)"], trans.rate["(12)", "(2)"]
                 , trans.rate["(23)", "(2)"], trans.rate["(13)", "(3)"], trans.rate["(23)", "(3)"])

    ## In this case the vector will very often have 0 positions. Should just assume some 0s.
    trans.tmp[which(trans.tmp > 0)] <- trans.tmp[which(trans.tmp > 0)] + max(pars.tmp)
    pars.tmp <- c(pars.tmp, trans.tmp)
    ## Number of free parameters for the model:
    np <- max(pars.tmp)
    ## For this all to work we need to have all integers from 1 to np.
    if( !all( pars.tmp[pars.tmp > 0] %in% pars.tmp ) ){
        stop("Wrong internal set-up of parameters.")
    }

    ## Need to fix the parameter indexes for the extirpation and local extinction.
    ## Note that the rate of d12_1 is the same as ext2 in the original geosse model.
    ## So if extirpation is not separated from local extinction, then we need to make this constrain.
    ## Check if any extirpation is set to be estimated independently:
    ## The indexation below is not relative and will break if the order of the parameters change.

    ## Note that in some cases the widespread areas do not belong to the model.
    ## We can get this info by checking if the respective speciation exists. If it is set to 0, then the extirpation associated with that same widespread area also should be 0.
    extirpation.trans.par <- pars.tmp[22:27]
    if( any(extirpation.trans.par > 0) ){
        ## Set the ones to be estimated, but constrain the others.
        ## The transitions that are zero need to be set to the correspondent extinction parameter.
        ## Need to find a good way to scale this up when adding the hidden rates.
        if( pars.tmp[4] > 0 ){
            pars.tmp[ 24 ][ pars.tmp[ 24 ] == 0 ] <- pars.tmp[7]
            pars.tmp[ 22 ][ pars.tmp[ 22 ] == 0 ] <- pars.tmp[8]
        }
        if( pars.tmp[6] > 0 ){        
            pars.tmp[ 26 ][ pars.tmp[ 26 ] == 0 ] <- pars.tmp[7]
            pars.tmp[ 23 ][ pars.tmp[ 23 ] == 0 ] <- pars.tmp[9]
        }
        if( pars.tmp[5] > 0 ){        
            pars.tmp[ 27 ][ pars.tmp[ 27 ] == 0 ] <- pars.tmp[8]
            pars.tmp[ 25 ][ pars.tmp[ 25 ] == 0 ] <- pars.tmp[9]
        }
    } else{
        ## Constrain all the extirpation events to be the same as extinction. (13 to 18)
        if( pars.tmp[4] > 0 ){
            pars.tmp[ 24 ] <- pars.tmp[7]
            pars.tmp[ 22 ] <- pars.tmp[8]
        }
        if( pars.tmp[6] > 0 ){        
            pars.tmp[ 26 ] <- pars.tmp[7]
            pars.tmp[ 23 ] <- pars.tmp[9]
        }
        if( pars.tmp[5] > 0 ){        
            pars.tmp[ 27 ] <- pars.tmp[8]
            pars.tmp[ 25 ] <- pars.tmp[9]
        }
    }
   
    cat("Initializing...", "\n")

    ## Data here is a named vector with states from 1 to 6 in the same order as phy$tip.label.
    ## Any required transformation need to be perfomed prior to this function call.
    freqs <- table(data)
    ## This will be equal to 1 under default values.
    samp.freq.tree <- Ntip(phy) / sum(table(data) / f[as.numeric(names(freqs))])

    ## The starting point for the parameters is the same as the geosse model.
    ## Just need to distribute the speciatiation and the extirpation parameters.
    ## These rates are in normal space. Need to log!
    if( sum(extirpation) == 0 ){
        ## In the case of a Yule model.
        init.pars.tmp <- starting.point.geosse(phy, eps=0, samp.freq.tree=samp.freq.tree)
    }else{
        init.pars.tmp <- starting.point.geosse(phy, eps=mag.san.start, samp.freq.tree=samp.freq.tree)
    }
    ## Original code:
    ## init.pars <- unname( c( rep(init.pars.tmp[1], times = 6), rep(init.pars.tmp[4], times = 3)
    ##                      , rep(init.pars.tmp[6], times = 3) ) )
    ## Need to construct the parameter vector to be passed to the likelihood function.
    ## The map back to the total parameteter vector is simpler because the positions will match.
    pars.lik <- rep(NA, times = np)
    s_id <- unique(pars.tmp[1:6])[unique(pars.tmp[1:6]) > 0]
    ext_id <- unique(pars.tmp[7:9])[unique(pars.tmp[7:9]) > 0]
    d_id <- seq(from = max(s_id, ext_id)+1, to = length(pars.lik))
    pars.lik[s_id] <- log( init.pars.tmp[1] )
    pars.lik[ext_id] <- log( init.pars.tmp[4] )
    pars.lik[d_id] <- log( init.pars.tmp[6] )
    
    ## If the model separates extirpation and extinction, then we need to correct the tail of these vectors.
    ## Record the potential extirpation transition parameters:
    if( any(extirpation.trans.par > 0) ){
        ## The true extirpation parameters will be at the end of the likelihood vector:
        tail.id <- np - ( sum(extirpation.trans.par > 0) + 1 )
        pars.lik[tail.id:np] <- init.pars.tmp[4]
    }
    
    ## So here we have two important vectors:
    ## "pars.lik" used to estimate the parameters of the model
    ## indexes in pars.tmp are the key to translate "pars.lik" into a longer vector to be passed to the C code.
    
    if(bounded.search == TRUE){
        upper.pars <- rep(NA, times = length(pars.lik))
        upper.pars[s_id] <- log(speciation.upper)
        upper.pars[ext_id] <- log(extirpation.upper)
        upper.pars[d_id] <- log(trans.upper)

        if( any(extirpation.trans.par > 0) ){
            ## If this is true, then "tail.id" should exist in the working space.
            ## But better be safe than sorry!
            tail.id <- np - ( sum(extirpation.trans.par > 0) + 1 )
            upper.pars[tail.id:np] <- log(extirpation.upper)
        }
        
    } else{
        ## We set a really large bound to emulate no bounds. But be aware that the search cannot, of course, be done in an infinite parameter space.
        upper.pars <- rep(21, times = np)
    }

    ## Lower bound in log space:
    lower.pars <- rep(-20, times = np)

    ## ##########################
    ## Pre-likelihood steps (this are constants for the search):

    ## Data need to be a presence and absence matrix with columns ordered from state 1 to 6.
    states <- matrix(0, nrow = length(data), ncol = 6)
    for( i in 1:length(data) ){
        states[i,data[i]] <- 1
    }

    tot_time <- max(branching.times(phy))
    split.times <- sort(branching.times(phy), decreasing=TRUE)

    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
    phy <- reorder(phy, "pruningwise")
    anc <- unique(phy$edge[,1])
    TIPS <- 1:nb.tip

    ## The ODEs need to be integrated across the each of the branches of the tree.
    compD <- matrix(0, nrow=nb.tip + nb.node, ncol=6)
    compE <- matrix(0, nrow=nb.tip + nb.node, ncol=6)

    ## Initializes the tip sampling and sets internal nodes to be zero:
    ncols <- dim(compD)[2]
    for(i in 1:(nb.tip)){
        ## Probability to be in the observed state is 1. (With a weight given by the 'f' vector.)
        compD[i,] <- f * states[i,]
        ## Probability of extinction for the observed state is 0.
        compE[i,] <- 1 - f
    }

    ## ##########################
    
    if(sann == FALSE){
        if(bounded.search == TRUE){
            cat("Starting log-likelihood value:")
            start.lik <- opt_geohisse_3_one_rate(pars.lik, pars=pars.tmp, phy=phy
                                               , condition.on.survival=condition.on.survival
                                               , root.type=root.type, root.p=root.p, ode.eps=ode.eps
                                               , split.times=split.times, nb.tip=nb.tip, nb.node=nb.node
                                               , anc=anc, compD=compD, compE=compE, bad.likelihood=10000000)
            cat( -1 * start.lik )
            cat( "\n" )
            cat("Finished. Beginning bounded subplex routine...", "\n")
            opts <- list("algorithm" = "NLOPT_LN_SBPLX", "maxeval" = 100000, "ftol_rel" = max.tol)
            ## The evaluation function here need to be changed. We will use one for each of the models.
            ## Need to make this function: "opt_geohisse_3_one_rate"
            ## I think that the objects are fine and have all the needed information. Just need to adapt the downstream.
            out <- nloptr(x0=pars.lik, eval_f=opt_geohisse_3_one_rate, ub=upper.pars, lb=lower.pars
                        , opts=opts
                        , pars=pars.tmp, phy=phy, condition.on.survival=condition.on.survival
                        , root.type=root.type, root.p=root.p, ode.eps=ode.eps
                        , split.times=split.times, nb.tip=nb.tip, nb.node=nb.node
                        , anc=anc, compD=compD, compE=compE, bad.likelihood=10000000)
            ## out <- nloptr(x0=pars.lik, eval_f=opt_geohisse_3_one_rate, ub=upper.pars, lb=lower.pars, opts=opts
            ##             , pars=pars.tmp, phy=phy, data=data, f=f, condition.on.survival=condition.on.survival
            ##             , root.type=root.type, root.p=root.p, np=np, ode.eps=ode.eps)
            ## Here need to expand the pars.lik into the long format.
            solution <- expand.pars(lik = exp(out$solution), long = pars.tmp)
            loglik <- -out$objective
        }else{
            cat("Finished. Beginning subplex routine...", "\n")
            out <- subplex(pars.lik, fn=opt_geohisse_3_one_rate, control=list(reltol=max.tol, parscale=rep(0.1, np))
                         , pars=pars.tmp, phy=phy, condition.on.survival=condition.on.survival
                         , root.type=root.type, root.p=root.p, ode.eps=ode.eps
                         , split.times=split.times, nb.tip=nb.tips, nb.node=nb.node
                         , anc=anc, compD=compD, compE=compE, bad.likelihood=10000000)
            ## out <- subplex(pars.lik, fn=opt_geohisse_3_one_rate, control=list(reltol=max.tol, parscale=rep(0.1, np))
            ##              , pars=pars.tmp, phy=phy, data=data, f=f, hidden.states=hidden.areas
            ##              , condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p
            ##              , np=np, ode.eps=ode.eps)
            solution <- expand.pars(lik = exp(out$par), long = pars.tmp)
            loglik <- -out$value
        }
    }else{
        cat("Finished. Beginning simulated annealing...", "\n")
        out.sann <- GenSA(pars.lik, fn=opt_geohisse_3_one_rate, lower=lower.pars, upper=upper.pars
                        , control=list(max.call=sann.its)
                        , pars=pars.tmp, phy=phy, condition.on.survival=condition.on.survival
                        , root.type=root.type, root.p=root.p, ode.eps=ode.eps
                        , split.times=split.times, nb.tip=nb.tips, nb.node=nb.node
                        , anc=anc, compD=compD, compE=compE, bad.likelihood=10000000)
        ## out.sann <- GenSA(pars.lik, fn=opt_geohisse_3_one_rate, lower=lower.pars, upper=upper.pars
        ##                 , control=list(max.call=sann.its)
        ##                 , pars=pars.tmp, phy=phy, data=data, f=f
        ##                 , condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p
        ##                 , np=np, ode.eps=ode.eps)
        cat("Finished. Refining using subplex routine...", "\n")
        out <- nloptr(x0=out.sann$par, eval_f=opt_geohisse_3_one_rate, ub=upper.pars, lb=lower.pars, opts=opts
                    , pars=pars.tmp, phy=phy, condition.on.survival=condition.on.survival
                    , root.type=root.type, root.p=root.p, ode.eps=ode.eps
                    , split.times=split.times, nb.tip=nb.tips, nb.node=nb.node
                    , anc=anc, compD=compD, compE=compE, bad.likelihood=10000000)
        ## out <- nloptr(x0=out.sann$par, eval_f=opt_geohisse_3_one_rate, ub=upper.pars, lb=lower.pars, opts=opts
        ##             , pars=pars.tmp, phy=phy, data=data, f=f, condition.on.survival=condition.on.survival
        ##             , root.type=root.type, root.p=root.p, np=np, ode.eps=ode.eps)
        solution <- expand.pars(lik = exp(out$solution), long = pars.tmp)
        loglik <- -out$objective
    }

    ## A new cool thing and easy to do with the updated data input would be to return the vector with the name of the areas as informed by the user.
    names(solution) <- c("s1", "s2", "s3", "s12", "s23", "s13", "x1", "x2", "x3", "d1_2"
                       , "d2_1", "d2_3", "d3_2", "d1_3", "d3_1", "d1_12", "d1_13", "d2_12"
                       , "d2_23", "d3_13", "d3_23", "d12_1", "d13_1", "d12_2", "d23_2"
                       , "d13_3", "d23_3")
    
    cat("Finished. Summarizing results...", "\n")

    ## This is just a partial list. The complete list will be the output of the higher level function that is calling this function to run.
    ## obj <- list(loglik = loglik, AIC = -2*loglik+2*np, AICc = -2*loglik+(2*np*(Ntip(phy)/(Ntip(phy)-np-1)))
    ##           , solution=solution, index.par=pars, f=f, hidden.areas=hidden.areas
    ##           , assume.cladogenetic=assume.cladogenetic, condition.on.survival=condition.on.survival
    ##           , root.type=root.type, root.p=root.p, phy=phy, data=data, trans.matrix=trans.rate
    ##           , max.tol=max.tol, starting.vals=ip, upper.bounds=upper, lower.bounds=lower, ode.eps=ode.eps)
    ## class(obj) <- append(class(obj), "geohisse.fit")
    ## return(obj)
    fit.out <- list(loglik = loglik, AIC = -2*loglik+2*np, AICc = -2*loglik+(2*np*(Ntip(phy)/(Ntip(phy)-np-1)))
                  , solution=solution, index.par=pars.tmp, f=f, max.tol=max.tol, starting.vals=pars.lik
                  , upper.bounds=upper.pars, lower.bounds=lower.pars, ode.eps=ode.eps)
    return( fit.out )
}

## Run the fit of a GeoHiSSE model with two hidden rates.
## This function is generic and will work with both full and "null" models.
geohisse_3_two_rate <- function(phy, data, f, speciation, extirpation, trans.rate, condition.on.survival, root.type, root.p, sann, sann.its, bounded.search, max.tol, mag.san.start, starting.vals, speciation.upper, extirpation.upper, trans.upper, ode.eps){
    ## This function is called internally and should not need any default value.

    extirpation.tmp <- extirpation
    extirpation.tmp[which(extirpation.tmp > 0)] <- extirpation.tmp[which(extirpation.tmp > 0)] + max(speciation)
    pars.tmp <- c(speciation, extirpation.tmp)

    ## First increase the indexes of the trans.rate matrix:
    trans.rate[trans.rate > 0 & !is.na(trans.rate)] <- trans.rate[trans.rate > 0 & !is.na(trans.rate)] + max( pars.tmp )
    
    ## Before getting the parameters from trans.rate we need to check if local extirpation is 0.
    ## If it is, then set it to the respective extirpation parameter. Otherwise, leave it be.
    for( i in seq(from = 0, by = 6, length.out = nrow( trans.rate )/6) ){
        trans.rate[4+i,1+i] <- ifelse(trans.rate[4+i,1+i] == 0, yes = extirpation.tmp[2+(i/2)], no = trans.rate[4+i,1+i])
        trans.rate[4+i,2+i] <- ifelse(trans.rate[4+i,2+i] == 0, yes = extirpation.tmp[1+(i/2)], no = trans.rate[4+i,2+i])
        trans.rate[5+i,1+i] <- ifelse(trans.rate[5+i,1+i] == 0, yes = extirpation.tmp[3+(i/2)], no = trans.rate[5+i,1+i])
        trans.rate[5+i,3+i] <- ifelse(trans.rate[5+i,3+i] == 0, yes = extirpation.tmp[1+(i/2)], no = trans.rate[5+i,3+i])
        trans.rate[6+i,2+i] <- ifelse(trans.rate[6+i,2+i] == 0, yes = extirpation.tmp[3+(i/2)], no = trans.rate[6+i,2+i])
        trans.rate[6+i,3+i] <- ifelse(trans.rate[6+i,3+i] == 0, yes = extirpation.tmp[2+(i/2)], no = trans.rate[6+i,3+i])
    }    
    
    ## If we use the indexes instead of the names I can make a general call for the positions.
    ## First, let's make this work with the name calling way.
    trans.tmp <- c(trans.rate[ "(1A)", "(2A)"] , trans.rate["(2A)", "(1A)"] , trans.rate["(2A)", "(3A)"]
                 , trans.rate[ "(3A)", "(2A)"] , trans.rate["(1A)", "(3A)"] , trans.rate["(3A)", "(1A)"]
                 , trans.rate[ "(1A)","(12A)"], trans.rate[ "(1A)","(13A)"], trans.rate[ "(2A)","(12A)"]
                 , trans.rate[ "(2A)","(23A)"], trans.rate[ "(3A)","(13A)"], trans.rate[ "(3A)","(23A)"]
                 , trans.rate["(12A)", "(1A)"], trans.rate["(13A)", "(1A)"], trans.rate["(12A)", "(2A)"]
                 , trans.rate["(23A)", "(2A)"], trans.rate["(13A)", "(3A)"], trans.rate["(23A)", "(3A)"]
                 , trans.rate[ "(1B)", "(2B)"] , trans.rate["(2B)", "(1B)"] , trans.rate["(2B)", "(3B)"]
                 , trans.rate[ "(3B)", "(2B)"] , trans.rate["(1B)", "(3B)"] , trans.rate["(3B)", "(1B)"]
                 , trans.rate[ "(1B)","(12B)"], trans.rate[ "(1B)","(13B)"], trans.rate[ "(2B)","(12B)"]
                 , trans.rate[ "(2B)","(23B)"], trans.rate[ "(3B)","(13B)"], trans.rate[ "(3B)","(23B)"]
                 , trans.rate["(12B)", "(1B)"], trans.rate["(13B)", "(1B)"], trans.rate["(12B)", "(2B)"]
                 , trans.rate["(23B)", "(2B)"], trans.rate["(13B)", "(3B)"], trans.rate["(23B)", "(3B)"])

    ## Record the indexes for the extirpation parameters:
    sep_ext_tmp <- c(trans.rate["(12A)", "(1A)"], trans.rate["(13A)", "(1A)"], trans.rate["(12A)", "(2A)"]
                   , trans.rate["(23A)", "(2A)"], trans.rate["(13A)", "(3A)"], trans.rate["(23A)", "(3A)"]
                   , trans.rate["(12B)", "(1B)"], trans.rate["(13B)", "(1B)"], trans.rate["(12B)", "(2B)"]
                   , trans.rate["(23B)", "(2B)"], trans.rate["(13B)", "(3B)"], trans.rate["(23B)", "(3B)"])

    ## In this case we need more positions for the transitions between the hidden layers.
    hidden.tmp <- c(trans.rate[ "(1A)", "(1B)"] , trans.rate["(2A)", "(2B)"] , trans.rate["(3A)", "(3B)"]
                  , trans.rate[ "(1B)", "(1A)"] , trans.rate["(2B)", "(2A)"] , trans.rate["(3B)", "(3A)"]
                  , trans.rate[ "(12A)", "(12B)"] , trans.rate["(23A)", "(23B)"] , trans.rate["(13A)", "(13B)"]
                  , trans.rate[ "(12B)", "(12A)"] , trans.rate["(23B)", "(23A)"] , trans.rate["(13B)", "(13A)"])
    trans.tmp <- c(trans.tmp, hidden.tmp)

    pars.tmp <- c(pars.tmp, trans.tmp)
    ## Number of free parameters for the model:
    np <- max(pars.tmp)
    ## For this all to work we need to have all integers from 1 to np.
    if( !all( 1:np %in% pars.tmp ) ){
        stop("Wrong internal set-up of parameters.")
    }
   
    cat("Initializing...", "\n")

    ## Data here is a named vector with states from 1 to 6 in the same order as phy$tip.label.
    ## Any required transformation need to be perfomed prior to this function call.
    freqs <- table(data)
    ## This will be equal to 1 under default values.
    samp.freq.tree <- Ntip(phy) / sum(table(data) / f[as.numeric(names(freqs))])

    ## The starting point for the parameters is the same as the geosse model.
    ## Just need to distribute the speciatiation and the extirpation parameters.
    ## These rates are in normal space. Need to log!
    if( sum(extirpation) == 0 ){
        ## In the case of a Yule model.
        init.pars.tmp <- starting.point.geosse(phy, eps=0, samp.freq.tree=samp.freq.tree)
    }else{
        init.pars.tmp <- starting.point.geosse(phy, eps=mag.san.start, samp.freq.tree=samp.freq.tree)
    }

    pars.lik <- rep(NA, times = np)
    s_id <- unique(pars.tmp[1:ncol(trans.rate)])[unique(pars.tmp[1:ncol(trans.rate)]) > 0]
    seq_get_ext <- (ncol(trans.rate)+1):(ncol(trans.rate)*1.5) ## Get the id for extirpation.
    ext_id <- unique(pars.tmp[seq_get_ext])[unique(pars.tmp[seq_get_ext]) > 0]
    ## NOTE: Extirpations (e.g., d12_1) are set with the same rate as dispersions (e.g., d1_12).
    ##       This differs from normal GeoHiSSE. But it should be all good.
    d_id <- seq(from = max(s_id, ext_id)+1, to = length(pars.lik))
    pars.lik[s_id] <- log( init.pars.tmp[1] )
    pars.lik[ext_id] <- log( init.pars.tmp[4] )
    pars.lik[d_id] <- log( init.pars.tmp[6] )

    ## Find the id for the extirpation parameters if necessary.
    if( !all( sep_ext_tmp %in% pars.tmp[c(s_id,ext_id)] ) ){
        ## Then, some parameters are exclusive to extirpation.
        ## Find the parameters id and search on the vector.
        ext_sep_pars <- sep_ext_tmp[ !sep_ext_tmp %in% pars.tmp[c(s_id,ext_id)] ]
        set_ext <- unique( pars.tmp[ pars.tmp > 0 ] ) %in% ext_sep_pars ## These are the par values we need to change.
        pars.lik[set_ext] <- log( init.pars.tmp[4] )
    }
    
    if(bounded.search == TRUE){
        upper.pars <- rep(NA, times = length(pars.lik))
        upper.pars[s_id] <- log(speciation.upper)
        upper.pars[ext_id] <- log(extirpation.upper)
        upper.pars[d_id] <- log(trans.upper)

        if( !all( sep_ext_tmp %in% pars.tmp[c(s_id,ext_id)] ) ){
            ## Then, some parameters are exclusive to extirpation.
            ## Find the parameters id and search on the vector.
            ext_sep_pars <- sep_ext_tmp[ !sep_ext_tmp %in% pars.tmp[c(s_id,ext_id)] ]
            set_ext <- unique( pars.tmp[ pars.tmp > 0 ] ) %in% ext_sep_pars ## These are the par values we need to change.
            upper.pars[set_ext] <- log(extirpation.upper)
        }
        
        ## if( any(extirpation.trans.par > 0) ){
        ##     ## If this is true, then "tail.id" should exist in the working space.
        ##     ## But better be safe than sorry!
        ##     tail.id <- np - ( sum(extirpation.trans.par > 0) + 1 )
        ##     upper.pars[tail.id:np] <- log(extirpation.upper)
        ## }
        
    } else{
        ## We set a really large bound to emulate no bounds. But be aware that the search cannot, of course, be done in an infinite parameter space.
        upper.pars <- rep(21, times = np)
    }

    ## Lower bound in log space:
    lower.pars <- rep(-20, times = np)

    ## ##########################
    ## Pre-likelihood steps (this are constants for the search):

    ## Data need to be a presence and absence matrix with columns ordered from state 1 to 6.
    ## With hidden states, then we need a larger matrix. The observed state is present for all hidden layers.
    
    states <- matrix(0, nrow = length(data), ncol = 6)
    for( i in 1:length(data) ){
        states[i,data[i]] <- 1
    }
    ## We then just need to "duplicate" the matrix for each hidden rate.
    states.tmp <- states
    f.tmp <- f
    for( i in 1:((ncol(trans.rate)/6)-1) ){
        states <- cbind(states, states.tmp)
        f <- c(f, f.tmp)
    }

    tot_time <- max(branching.times(phy))
    split.times <- sort(branching.times(phy), decreasing=TRUE)

    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
    phy <- reorder(phy, "pruningwise")
    anc <- unique(phy$edge[,1])
    TIPS <- 1:nb.tip

    ## The ODEs need to be integrated across the each of the branches of the tree.
    compD <- matrix(0, nrow=nb.tip + nb.node, ncol=ncol(states) )
    compE <- matrix(0, nrow=nb.tip + nb.node, ncol=ncol(states) )

    ## Initializes the tip sampling and sets internal nodes to be zero:
    ncols <- dim(compD)[2]
    for(i in 1:(nb.tip)){
        ## Probability to be in the observed state is 1. (With a weight given by the 'f' vector.)
        compD[i,] <- f * states[i,]
        ## Probability of extinction for the observed state is 0.
        compE[i,] <- 1 - f
    }

    ## ##########################
    
    if(sann == FALSE){
        if(bounded.search == TRUE){
            cat("Starting log-likelihood value:")
            start.lik <- opt_geohisse_3_two_rate(pars.lik, pars=pars.tmp, phy=phy
                                               , condition.on.survival=condition.on.survival
                                               , root.type=root.type, root.p=root.p, ode.eps=ode.eps
                                               , split.times=split.times, nb.tip=nb.tip, nb.node=nb.node
                                               , anc=anc, compD=compD, compE=compE, bad.likelihood=10000000)
            cat( -1 * start.lik )
            cat( "\n" )
            cat("Finished. Beginning bounded subplex routine...", "\n")
            opts <- list("algorithm" = "NLOPT_LN_SBPLX", "maxeval" = 100000, "ftol_rel" = max.tol)
            ## The evaluation function here need to be changed. We will use one for each of the models.
            ## Need to make this function: "opt_geohisse_3_one_rate"
            ## I think that the objects are fine and have all the needed information. Just need to adapt the downstream.
            out <- nloptr(x0=pars.lik, eval_f=opt_geohisse_3_two_rate, ub=upper.pars, lb=lower.pars
                        , opts=opts
                        , pars=pars.tmp, phy=phy, condition.on.survival=condition.on.survival
                        , root.type=root.type, root.p=root.p, ode.eps=ode.eps
                        , split.times=split.times, nb.tip=nb.tip, nb.node=nb.node
                        , anc=anc, compD=compD, compE=compE, bad.likelihood=10000000)
            ## Here need to expand the pars.lik into the long format.
            solution <- expand.pars(lik = exp(out$solution), long = pars.tmp)
            loglik <- -out$objective
        }else{
            cat("Finished. Beginning subplex routine...", "\n")
            out <- subplex(pars.lik, fn=opt_geohisse_3_two_rate, control=list(reltol=max.tol, parscale=rep(0.1, np))
                         , pars=pars.tmp, phy=phy, condition.on.survival=condition.on.survival
                         , root.type=root.type, root.p=root.p, ode.eps=ode.eps
                         , split.times=split.times, nb.tip=nb.tips, nb.node=nb.node
                         , anc=anc, compD=compD, compE=compE, bad.likelihood=10000000)
            solution <- expand.pars(lik = exp(out$par), long = pars.tmp)
            loglik <- -out$value
        }
    }else{
        cat("Finished. Beginning simulated annealing...", "\n")
        out.sann <- GenSA(pars.lik, fn=opt_geohisse_3_two_rate, lower=lower.pars, upper=upper.pars
                        , control=list(max.call=sann.its)
                        , pars=pars.tmp, phy=phy, condition.on.survival=condition.on.survival
                        , root.type=root.type, root.p=root.p, ode.eps=ode.eps
                        , split.times=split.times, nb.tip=nb.tips, nb.node=nb.node
                        , anc=anc, compD=compD, compE=compE, bad.likelihood=10000000)
        cat("Finished. Refining using subplex routine...", "\n")
        out <- nloptr(x0=out.sann$par, eval_f=opt_geohisse_3_two_rate, ub=upper.pars, lb=lower.pars, opts=opts
                    , pars=pars.tmp, phy=phy, condition.on.survival=condition.on.survival
                    , root.type=root.type, root.p=root.p, ode.eps=ode.eps
                    , split.times=split.times, nb.tip=nb.tips, nb.node=nb.node
                    , anc=anc, compD=compD, compE=compE, bad.likelihood=10000000)
        ## out <- nloptr(x0=out.sann$par, eval_f=opt_geohisse_3_one_rate, ub=upper.pars, lb=lower.pars, opts=opts
        ##             , pars=pars.tmp, phy=phy, data=data, f=f, condition.on.survival=condition.on.survival
        ##             , root.type=root.type, root.p=root.p, np=np, ode.eps=ode.eps)
        solution <- expand.pars(lik = exp(out$solution), long = pars.tmp)
        loglik <- -out$objective
    }

    ## Need to make this vector general so it works with more hidden rates:
    hidden.rates <- ncol( trans.rate ) / 6
    
    add_letters <- function(x, n){
        return( paste0(x, rep(LETTERS[1:n], each = length(x)) ) )
    }

    names(solution) <- c( add_letters(x = c("s1", "s2", "s3", "s12", "s23", "s13"), n = hidden.rates)
                       , add_letters(x = c("x1", "x2", "x3"), n = hidden.rates)
                       , add_letters(x = c("d1_2", "d2_1", "d2_3", "d3_2", "d1_3", "d3_1", "d1_12", "d1_13"
                                         , "d2_12", "d2_23", "d3_13", "d3_23", "d12_1", "d13_1", "d12_2"
                                         , "d23_2", "d13_3", "d23_3"), n = hidden.rates)
                         ## The last elements are easier to just write down... :(
                       , c("d1A_1B", "d2A_2B", "d3A_3B", "d1B_1A", "d2B_2A", "d3B_3A", "d12A_12B", "d23A_23B", "d13A_13B"
                         , "d12B_12A", "d23B_23A", "d13B_13A") )
    
    cat("Finished. Summarizing results...", "\n")

    fit.out <- list(loglik = loglik, AIC = -2*loglik+2*np, AICc = -2*loglik+(2*np*(Ntip(phy)/(Ntip(phy)-np-1)))
                  , solution=solution, index.par=pars.tmp, f=f, max.tol=max.tol, starting.vals=pars.lik
                  , upper.bounds=upper.pars, lower.bounds=lower.pars, ode.eps=ode.eps)
    return( fit.out )
}

## #################################################################
## The general function that should fit any number of hidden rates to the model.
## #################################################################

geohisse_3_multi_rate <- function(phy, data, f, speciation, extirpation, trans.rate, condition.on.survival, root.type, root.p, sann, sann.its, bounded.search, max.tol, mag.san.start, starting.vals, speciation.upper, extirpation.upper, trans.upper, ode.eps){

    ## Another way to do this is to inflate all the parameter vector.
    ## We can add 0's to indicate the parameters that are not in use in this case.
    ## In this case we have vector with the full model. This follows the strategy used by Jeremy in GeoHiSSE.
    extirpation.full <- rep(0, times = 3*12)
    speciation.full <- rep(0, times = 6*12)
    extirpation.full[1:length(extirpation)] <- extirpation
    speciation.full[1:length(speciation)] <- speciation
    
    extirpation.tmp <- extirpation.full
    extirpation.tmp[which(extirpation.tmp > 0)] <- extirpation.tmp[which(extirpation.tmp > 0)] + max(speciation)
    pars.tmp <- c(speciation.full, extirpation.tmp)

    ## First increase the indexes of the trans.rate matrix:
    trans.rate.copy <- trans.rate
    dim.trans.rate <- dim( trans.rate )
    trans.rate <- matrix(data = 0, nrow = 6*12, ncol = 6*12)
    trans.rate[1:dim.trans.rate[1], 1:dim.trans.rate[2]] <- trans.rate.copy
    trans.rate[trans.rate > 0 & !is.na(trans.rate)] <- trans.rate[trans.rate > 0 & !is.na(trans.rate)] + max( pars.tmp )

    ## Before getting the parameters from trans.rate we need to check if local extirpation is 0.
    ## If it is, then set it to the respective extirpation parameter. Otherwise, leave it be.
    ## Note that, because trans.rate has the maximum size, this code will work with any number of hidden rates.
    for( i in seq(from = 0, by = 6, length.out = 12) ){
        trans.rate[4+i,1+i] <- ifelse(trans.rate[4+i,1+i] == 0, yes = extirpation.tmp[2+(i/2)], no = trans.rate[4+i,1+i])
        trans.rate[4+i,2+i] <- ifelse(trans.rate[4+i,2+i] == 0, yes = extirpation.tmp[1+(i/2)], no = trans.rate[4+i,2+i])
        trans.rate[5+i,1+i] <- ifelse(trans.rate[5+i,1+i] == 0, yes = extirpation.tmp[3+(i/2)], no = trans.rate[5+i,1+i])
        trans.rate[5+i,3+i] <- ifelse(trans.rate[5+i,3+i] == 0, yes = extirpation.tmp[1+(i/2)], no = trans.rate[5+i,3+i])
        trans.rate[6+i,2+i] <- ifelse(trans.rate[6+i,2+i] == 0, yes = extirpation.tmp[3+(i/2)], no = trans.rate[6+i,2+i])
        trans.rate[6+i,3+i] <- ifelse(trans.rate[6+i,3+i] == 0, yes = extirpation.tmp[2+(i/2)], no = trans.rate[6+i,3+i])
    }

    ## Here we are going to use the same system of multipliers above to get the indexes in the correct order.
    trans.tmp <- vector(mode="numeric")
    for( i in seq(from = 0, by = 6, length.out = 12) ){
        get.sub.mat <- c(trans.rate[1+i,2+i], trans.rate[2+i,1+i], trans.rate[2+i,3+i]
                     , trans.rate[3+i,2+i], trans.rate[1+i,3+i], trans.rate[3+i,1+i]
                     , trans.rate[1+i,4+i], trans.rate[1+i,5+i], trans.rate[2+i,4+i]
                     , trans.rate[2+i,6+i], trans.rate[3+i,5+i], trans.rate[3+i,6+i]
                     , trans.rate[4+i,1+i], trans.rate[5+i,1+i], trans.rate[4+i,2+i]
                     , trans.rate[6+i,2+i], trans.rate[5+i,3+i], trans.rate[6+i,3+i])
        trans.tmp <- c(trans.tmp, get.sub.mat)
    }

    ## Record the indexes for the extirpation parameters:
    sep_ext_tmp <- vector(mode="numeric")
    for( i in seq(from = 0, by = 6, length.out = 12) ){
        get.sub.mat <- c(trans.rate[4+i,1+i], trans.rate[5+i,1+i], trans.rate[4+i,2+i]
                       , trans.rate[6+i,2+i], trans.rate[5+i,3+i], trans.rate[6+i,3+i])
        sep_ext_tmp <- c(sep_ext_tmp, get.sub.mat)
    }

    pars.tmp <- c(pars.tmp, trans.tmp)
    
    ## Now we just have the rate of transition between the layers.
    ## Here we will assume a global rate of transition between the layers.
    n.hidden.rates <- length(speciation) / 6
    trans.hidden.par <- rep(0.0, times = 12)
    trans.hidden.par[1:n.hidden.rates] <- trans.rate[1,7] ## One of the hidden layer transitions.
    pars.tmp <- c(pars.tmp, trans.hidden.par) ## This is the last parameter on the vector.

    ## Number of free parameters for the model:
    np <- max(pars.tmp)
    ## For this all to work we need to have all integers from 1 to np.
    if( !all( 1:np %in% pars.tmp ) ){
        stop("Wrong internal set-up of parameters.")
    }
   
    cat("Initializing...", "\n")

    ## Data here is a named vector with states from 1 to 6 in the same order as phy$tip.label.
    ## Any required transformation need to be perfomed prior to this function call.
    freqs <- table(data)
    ## This will be equal to 1 under default values.
    samp.freq.tree <- Ntip(phy) / sum(table(data) / f[as.numeric(names(freqs))])

    ## The starting point for the parameters is the same as the geosse model.
    ## Just need to distribute the speciatiation and the extirpation parameters.
    ## These rates are in normal space. Need to log!
    if( sum(extirpation) == 0 ){
        ## In the case of a Yule model.
        init.pars.tmp <- starting.point.geosse(phy, eps=0, samp.freq.tree=samp.freq.tree)
    }else{
        init.pars.tmp <- starting.point.geosse(phy, eps=mag.san.start, samp.freq.tree=samp.freq.tree)
    }

    pars.lik <- rep(NA, times = np)
    s_id <- unique(pars.tmp[1:72])[unique(pars.tmp[1:72]) > 0]
    ext_id <- unique(pars.tmp[73:108])[unique(pars.tmp[73:108]) > 0]
    d_id <- seq(from = max(s_id, ext_id)+1, to = length(pars.lik))
    pars.lik[s_id] <- log( init.pars.tmp[1] )
    pars.lik[ext_id] <- log( init.pars.tmp[4] )
    pars.lik[d_id] <- log( init.pars.tmp[6] )
    
    if(bounded.search == TRUE){
        upper.pars <- rep(NA, times = length(pars.lik))
        upper.pars[s_id] <- log(speciation.upper)
        upper.pars[ext_id] <- log(extirpation.upper)
        upper.pars[d_id] <- log(trans.upper)

        if( !all( sep_ext_tmp %in% pars.tmp[73:108] ) ){
            ## Then, some parameters are exclusive to extirpation.
            ## Find the parameters id and search on the vector.
            ext_sep_pars <- sep_ext_tmp[ !sep_ext_tmp %in% pars.tmp[73:108] ]
            set_ext <- unique( pars.tmp[ pars.tmp > 0 ] ) %in% ext_sep_pars ## These are the par values we need to change.
            upper.pars[set_ext] <- log(extirpation.upper)
        }
        
    } else{
        ## We set a really large bound to emulate no bounds. But be aware that the search cannot, of course, be done in an infinite parameter space.
        upper.pars <- rep(21, times = np)
    }

    ## Lower bound in log space:
    lower.pars <- rep(-20, times = np)

    ## ##########################
    ## Pre-likelihood steps (this are constants for the search):

    ## Data need to be a presence and absence matrix with columns ordered from state 1 to 6.
    ## With hidden states, then we need a larger matrix. The observed state is present for all hidden layers.
    
    states <- matrix(0, nrow = length(data), ncol = 6)
    for( i in 1:length(data) ){
        states[i,data[i]] <- 1
    }
    ## We then just need to "duplicate" the matrix for each hidden rate.
    states.tmp <- states
    f.tmp <- f
    for( i in 1:((ncol(trans.rate)/6)-1) ){
        states <- cbind(states, states.tmp)
        f <- c(f, f.tmp)
    }

    tot_time <- max(branching.times(phy))
    split.times <- sort(branching.times(phy), decreasing=TRUE)

    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
    phy <- reorder(phy, "pruningwise")
    anc <- unique(phy$edge[,1])
    TIPS <- 1:nb.tip

    ## The ODEs need to be integrated across the each of the branches of the tree.
    compD <- matrix(0, nrow=nb.tip + nb.node, ncol=ncol(states) )
    compE <- matrix(0, nrow=nb.tip + nb.node, ncol=ncol(states) )

    ## Initializes the tip sampling and sets internal nodes to be zero:
    ncols <- dim(compD)[2]
    for(i in 1:(nb.tip)){
        ## Probability to be in the observed state is 1. (With a weight given by the 'f' vector.)
        compD[i,] <- f * states[i,]
        ## Probability of extinction for the observed state is 0.
        compE[i,] <- 1 - f
    }

    ## ##########################
    
    if(sann == FALSE){
        if(bounded.search == TRUE){
            cat("Starting log-likelihood value:")
            start.lik <- opt_geohisse_3_multi_rate(pars.lik, pars=pars.tmp, phy=phy
                                                 , condition.on.survival=condition.on.survival
                                                 , root.type=root.type, root.p=root.p, ode.eps=ode.eps
                                                 , split.times=split.times, nb.tip=nb.tip, nb.node=nb.node
                                                 , anc=anc, compD=compD, compE=compE, bad.likelihood=10000000)
            cat( -1 * start.lik )
            cat( "\n" )
            cat("Finished. Beginning bounded subplex routine...", "\n")
            opts <- list("algorithm" = "NLOPT_LN_SBPLX", "maxeval" = 100000, "ftol_rel" = max.tol)
            ## The evaluation function here need to be changed. We will use one for each of the models.
            ## Need to make this function: "opt_geohisse_3_one_rate"
            ## I think that the objects are fine and have all the needed information. Just need to adapt the downstream.
            out <- nloptr(x0=pars.lik, eval_f=opt_geohisse_3_multi_rate, ub=upper.pars, lb=lower.pars
                        , opts=opts
                        , pars=pars.tmp, phy=phy, condition.on.survival=condition.on.survival
                        , root.type=root.type, root.p=root.p, ode.eps=ode.eps
                        , split.times=split.times, nb.tip=nb.tip, nb.node=nb.node
                        , anc=anc, compD=compD, compE=compE, bad.likelihood=10000000)
            ## Here need to expand the pars.lik into the long format.
            solution <- expand.pars(lik = exp(out$solution), long = pars.tmp)
            loglik <- -out$objective
        }else{
            cat("Finished. Beginning subplex routine...", "\n")
            out <- subplex(pars.lik, fn=opt_geohisse_3_multi_rate, control=list(reltol=max.tol, parscale=rep(0.1, np))
                         , pars=pars.tmp, phy=phy, condition.on.survival=condition.on.survival
                         , root.type=root.type, root.p=root.p, ode.eps=ode.eps
                         , split.times=split.times, nb.tip=nb.tips, nb.node=nb.node
                         , anc=anc, compD=compD, compE=compE, bad.likelihood=10000000)
            solution <- expand.pars(lik = exp(out$par), long = pars.tmp)
            loglik <- -out$value
        }
    }else{
        cat("Finished. Beginning simulated annealing...", "\n")
        out.sann <- GenSA(pars.lik, fn=opt_geohisse_3_multi_rate, lower=lower.pars, upper=upper.pars
                        , control=list(max.call=sann.its)
                        , pars=pars.tmp, phy=phy, condition.on.survival=condition.on.survival
                        , root.type=root.type, root.p=root.p, ode.eps=ode.eps
                        , split.times=split.times, nb.tip=nb.tips, nb.node=nb.node
                        , anc=anc, compD=compD, compE=compE, bad.likelihood=10000000)
        cat("Finished. Refining using subplex routine...", "\n")
        out <- nloptr(x0=out.sann$par, eval_f=opt_geohisse_3_multi_rate, ub=upper.pars, lb=lower.pars, opts=opts
                    , pars=pars.tmp, phy=phy, condition.on.survival=condition.on.survival
                    , root.type=root.type, root.p=root.p, ode.eps=ode.eps
                    , split.times=split.times, nb.tip=nb.tips, nb.node=nb.node
                    , anc=anc, compD=compD, compE=compE, bad.likelihood=10000000)
        ## out <- nloptr(x0=out.sann$par, eval_f=opt_geohisse_3_one_rate, ub=upper.pars, lb=lower.pars, opts=opts
        ##             , pars=pars.tmp, phy=phy, data=data, f=f, condition.on.survival=condition.on.survival
        ##             , root.type=root.type, root.p=root.p, np=np, ode.eps=ode.eps)
        solution <- expand.pars(lik = exp(out$solution), long = pars.tmp)
        loglik <- -out$objective
    }

    ## Need to make this vector general so it works with more hidden rates:
    hidden.rates <- ncol( trans.rate ) / 6
    
    add_letters <- function(x, n){
        return( paste0(x, rep(LETTERS[1:n], each = length(x)) ) )
    }

    names(solution) <- c( add_letters(x = c("s1", "s2", "s3", "s12", "s23", "s13"), n = hidden.rates)
                       , add_letters(x = c("x1", "x2", "x3"), n = hidden.rates)
                       , add_letters(x = c("d1_2", "d2_1", "d2_3", "d3_2", "d1_3", "d3_1", "d1_12", "d1_13"
                                         , "d2_12", "d2_23", "d3_13", "d3_23", "d12_1", "d13_1", "d12_2"
                                         , "d23_2", "d13_3", "d23_3"), n = hidden.rates)
                         ## The last elements are transitions to each of the 12 hidden layers.
                       , paste0("dn_", LETTERS[1:12]) )
    
    cat("Finished. Summarizing results...", "\n")

    fit.out <- list(loglik = loglik, AIC = -2*loglik+2*np, AICc = -2*loglik+(2*np*(Ntip(phy)/(Ntip(phy)-np-1)))
                  , solution=solution, index.par=pars.tmp, f=f, max.tol=max.tol, starting.vals=pars.lik
                  , upper.bounds=upper.pars, lower.bounds=lower.pars, ode.eps=ode.eps)
    return( fit.out )
}

## #################################################################################### 
## Optimization function for one_rate model:
## ####################################################################################

opt_geohisse_3_one_rate <- function(p, pars, phy, bad.likelihood=10000000, condition.on.survival, root.type, root.p, ode.eps, split.times, nb.tip, nb.node, anc, compD, compE) {
    
    ## Generates the final vector with the appropriate parameter estimates in the right place:
    p.new <- exp(p)
    
    ## Expand the parameter vector into the long form for the C code.
    model.vec <- expand.pars(lik = p.new, long = pars)
    
    logcomp <- c() ## Keep track of the log of the compensation.
    
    ## Start the postorder traversal indexing lists by node number:
    
    ## This is the real core of the likelhood function.
    ## Need to make sure that we have all the objects in the correct format to run these.

    ## Define the run function that will be called at each iteration:
    runSilent <- function(yini, times, model.vec){
        ## This is silencing the integration warnings.
        options(warn = -1)
        on.exit(options(warn = 0))
        capture.output( res <- lsoda(yini, times, func = "geosse_3_areas_derivs"
                                   , model.vec, initfunc="initmod_geosse_3_areas", dllname = "hisse"
                                   , rtol=1e-8, atol=1e-8) )
        return( res )
    }

    ## A constant to be used by the C code.
    NUMELEMENTS <- 27
    
    for (i in seq(from = 1, length.out = nb.node)) {
        ## A vector of all the internal nodes:
        focal <- anc[i]
        desRows <- which(phy$edge[,1]==focal)
        desNodes <- phy$edge[desRows,2]
        ## Note: when the tree has been reordered branching.times are no longer valid. Fortunately, we extract this information in the initial cache setup. Also focal is the rootward node, whereas desNodes represent a vector of all descendant nodes:
        rootward.age <- as.numeric( split.times[which(names(split.times)==focal)] )

        ## Matrix objects to store the results.
        v <- matrix(nrow = length(desRows), ncol = 6)
        phi <- matrix(nrow = length(desRows), ncol = 6)
        ## A counter to go through the matrix.
        v_phy_count <- 1
        ## BEGIN FOR LOOP over the nodes of the tree.
        ## It would be simple and better to separate this function between the case with and without the hidden states. It will be simpler to debug and less if else clauses to check.
        for (desIndex in sequence(length(desRows))){
            focal.edge.length <- phy$edge.length[desRows[desIndex]]
            tipward.age <- as.numeric( rootward.age - focal.edge.length )
            ## Strange rounding errors. A tip age should be zero. This ensures that:
            if(tipward.age < .Machine$double.eps^0.5){
                tipward.age <- 0.0
            }
            node.D <- compD[desNodes[desIndex],]
            node.E <- compE[desNodes[desIndex],]
            ## Call to lsoda that utilizes C code. Requires a lot of inputs. Note that for now we hardcode the NUMELEMENTS arguments. The reason for this is because with lsoda we can only pass a vector of parameters.
            ## Need to make sure that 'node.E' and 'node.D' have the correct quantities.
            yini <- c(node.E, node.D)
            ## yini <- c(E_0=cache$node.E[1], E_1=cache$node.E[2], E_01=cache$node.E[3], D_N0=cache$node.D[1], D_N1=cache$node.D[2], D_N2=cache$node.D[3])
            ## This informs the heights of the tipward and rootward nodes of the branch, so that we can compute the length of the branch.
            times <- c(tipward.age, rootward.age)

            ## The matrix below has the probabilities for each of the states at the tipward node for each of the
            ##     descendants of the focal node. Note that the loop above visited all the branches imediatelly
            ##     descendants of the focal node.
            ## We need these to compute the probability of each of the states at the focal node.
            ## The probability of each of the states at the focal node is equal to the probability that having
            ##     a given state at the focal node will generate tipward probabilities equal to the ones observed
            ##     or estimated.
            prob.subtree.cal.full <- runSilent(yini=yini, times=times, model.vec=model.vec)

            ## NEED TO UPDATE THIS CHECK! NOT ADJUSTED FOR THE MODEL.
            ## It will not break like this, but it is not checking the model correctly.
            
            ## ###### THIS CHECKS TO ENSURE THAT THE INTEGRATION WAS SUCCESSFUL ###########
            if(attributes(prob.subtree.cal.full)$istate[1] < 0){
                return(bad.likelihood)
            }else{
                prob.subtree.cal <- prob.subtree.cal.full[-1,-1]
            }
            ## ############################################################################

            if( any( is.nan( prob.subtree.cal ) ) ){
                return(bad.likelihood)
            }
            ## This is default and cannot change, but if we get a negative probability, discard the results:
            ## Note that this is only testing the Ds part of the ODEs.
            if( any( prob.subtree.cal < 0 ) ){
                return(bad.likelihood)
            }
            ## This can be modified at the input, but if the sum of the D's at the end of a branch are less than some value, then discard the results. A little more stringent than diversitree, but with difficult problems, this stabilizes things immensely.
            if( sum(prob.subtree.cal[7:12]) < ode.eps ){
                return(bad.likelihood)
            }
            
            ## Designating phi here because of its relation to Morlon et al (2011) and using "e" would be confusing:
            ## Here phi is a vector with the results of the integration so far in the E's.
            ## Then 'v' is the D's part.
            phi[v_phy_count,] <- prob.subtree.cal[1:6]
            v[v_phy_count,] <- prob.subtree.cal[7:12]
            v_phy_count <- v_phy_count + 1
        }

        ## The integration above was made from the tipward node to the rootward node.
        ## We computed a probability that the rootward node will be in each of the states for each of the branches.
        ## Of course, the rootward node is a single entity. So we need to combine the probabilities.
        ## Each column of v gives the probability for that state at the rootward node.
        ## An speciation happened at that node to form the observed descendants. This is an observed speciation
        ##     event, so its probability (just like the observed states) is equal to 1.
        ## The events track the infinitesimal moment just after the observed speciation event (i.e., no transitions
        ##     are possible after speciation, because we are tracking a single event.
        
        ## In the case of the endemic states at the nodes it is simple. Only one thing can happen:
        compD[focal,1] <- v[1,1] * v[2,1] * model.vec[1]
        compD[focal,2] <- v[1,2] * v[2,2] * model.vec[2]
        compD[focal,3] <- v[1,3] * v[2,3] * model.vec[3]

        ## When the state at the node is widespread then there are more than one thing possible.
        ## Here we need to enumerate all possible speciation types given that the state at the node was widespread.
        ## state 12
        compD[focal,4] <- 0.5 * (v[1,4] * v[2,1] + v[1,1] * v[2,4]) * model.vec[1] + 0.5 * (v[1,4] * v[2,2] + v[1,2] * v[2,4]) * model.vec[2] + 0.5 * (v[1,1] * v[2,2] + v[1,2] * v[2,1]) * model.vec[4]
        ## state 23
        compD[focal,5] <- 0.5 * (v[1,5] * v[2,2] + v[1,2] * v[2,5]) * model.vec[2] + 0.5 * (v[1,5] * v[2,3] + v[1,3] * v[2,5]) * model.vec[3] + 0.5 * (v[1,2] * v[2,3] + v[1,3] * v[2,2]) * model.vec[5]
        ## state 13
        compD[focal,6] <- 0.5 * (v[1,6] * v[2,1] + v[1,1] * v[2,6]) * model.vec[1] + 0.5 * (v[1,6] * v[2,3] + v[1,3] * v[2,6]) * model.vec[3] + 0.5 * (v[1,1] * v[2,3] + v[1,3] * v[2,1]) * model.vec[6]
        
        ## This is the probability of extinction at the node associated with each of the states.
        compE[focal,] <- phi[1,]

        ## This bit of code is used to fix a state at a node of the tree.
        ## Note that we essentially set the probability of that state at that node to 1 and 0 to everything else.
        ## When implementing this, make sure to set to the known vector BEFORE computing the probabilities!
        ## if(!is.null(node)){
        ##     if(node == focal){
        ##         fixer <- rep(0, times = 6)
        ##         fixer[state] <- 1
        ##         compD[focal,] <- compD[focal,] * fixer
        ##     }
        ## }
        
        ## #########################
        ## Logcompensation bit for dealing with underflow issues.

        ## Having a real problem here.
        ## The 'tmp' value is something negative.
        ## It is underflow problem. Probabilities cannot be negative. The values for the probabilities are too small.
        tmp <- sum(compD[focal,])
        compD[focal,] <- compD[focal,] / tmp
        ## Keep track of elements for the log compensation:
        logcomp <- c(logcomp, log(tmp))
        
    } ## FINISH FOR LOOP over the nodes of the tree.
    ## Compute quantities for the root node.
    
    root.node <- nb.tip + 1L
    if ( is.na(sum(log(compD[root.node,]))) || is.na(log(sum(1 - compE[root.node,]))) ){
        return(bad.likelihood)
    }else{
        if(root.type == "madfitz" | root.type == "herr_als"){
            if(is.null(root.p)){
                root.p <- compD[root.node,] / sum(compD[root.node,])
                ## Strange here. Why to accept NA as a possible value at the root?
                root.p[ which(is.na(root.p)) ] <- 0
            }
        }
        if(condition.on.survival == TRUE){
            if(root.type == "madfitz"){
                ## Here the endemic speciation rates contribute to multiple widespread areas, but this is happening to the classe model too. Should not be the issue.
                lambda <- c(model.vec[1:3], sum(model.vec[c(1,2,4)]), sum(model.vec[c(2,3,5)])
                          , sum(model.vec[c(1,3,6)]) )
                compD[root.node,] <- compD[root.node,] / sum(root.p * lambda * (1 - compE[root.node,])^2)
                ## Corrects for possibility that you have 0/0:
                compD[root.node,which(is.na(compD[root.node,]))] <- 0
                loglik <- log(sum(compD[root.node,] * root.p)) + sum(logcomp)
            }else{
                lambda <- c(model.vec[1:3], sum(model.vec[c(1,2,4)]), sum(model.vec[c(2,3,5)])
                          , sum(model.vec[c(1,3,6)]) )
                compD[root.node,] <- (compD[root.node,] * root.p) / (lambda * (1 - compE[root.node,])^2)
                ## Corrects for possibility that you have 0/0:
                compD[root.node,which(is.na(compD[root.node,]))] <- 0
                loglik <- log(sum(compD[root.node,])) + sum(logcomp)
            }
        }
    }
    
    ## There is a option to return another quantity. Looks a debug feature.
    ## if(get.phi == TRUE){
    ##     obj <- NULL
    ##     obj$compD.root <- compD[root.node,]/sum(compD[root.node,])
    ##     obj$compE <- compE
    ##     obj$root.p <- root.p
    ##     return(obj)
    ## }else{
    ##     return(loglik)
    ## }

    ## IMPORTANT: Note that in the original hisse implementation this is NOT the output directly to nloptr. hisse outputs this lik to another function that does -1 * lik and then pass it to nloptr! Here we are droppping this extra layer and outputing directly to nloptr.
    return( -1 * loglik )
}

## #################################################################################### 
## Optimization function for two_rate model:
## ####################################################################################

opt_geohisse_3_two_rate <- function(p, pars, phy, bad.likelihood=10000000, condition.on.survival, root.type, root.p, ode.eps, split.times, nb.tip, nb.node, anc, compD, compE) {
    
    ## Generates the final vector with the appropriate parameter estimates in the right place:
    p.new <- exp(p)
    
    ## Expand the parameter vector into the long form for the C code.
    model.vec <- expand.pars(lik = p.new, long = pars)
    
    logcomp <- c() ## Keep track of the log of the compensation.
    
    ## Start the postorder traversal indexing lists by node number:
    
    ## This is the real core of the likelhood function.
    ## Need to make sure that we have all the objects in the correct format to run these.

    ## Define the run function that will be called at each iteration:
    runSilent <- function(yini, times, model.vec){
        ## This is silencing the integration warnings.
        options(warn = -1)
        on.exit(options(warn = 0))
        capture.output( res <- lsoda(yini, times, func = "geosse_3_areas_two_rates_derivs"
                                   , model.vec, initfunc="initmod_geosse_3_areas_two_rates", dllname = "hisse"
                                   , rtol=1e-8, atol=1e-8) )
        return( res )
    }

    ## A constant to be used by the C code.
    NUMELEMENTS <- 66
    
    for (i in seq(from = 1, length.out = nb.node)) {
        ## A vector of all the internal nodes:
        focal <- anc[i]
        desRows <- which(phy$edge[,1]==focal)
        desNodes <- phy$edge[desRows,2]
        ## Note: when the tree has been reordered branching.times are no longer valid. Fortunately, we extract this information in the initial cache setup. Also focal is the rootward node, whereas desNodes represent a vector of all descendant nodes:
        rootward.age <- as.numeric( split.times[which(names(split.times)==focal)] )

        ## Matrix objects to store the results.
        v <- matrix(nrow = length(desRows), ncol = ncol(compD) )
        phi <- matrix(nrow = length(desRows), ncol = ncol(compD) )
        ## A counter to go through the matrix.
        v_phy_count <- 1
        ## BEGIN FOR LOOP over the nodes of the tree.
        ## It would be simple and better to separate this function between the case with and without the hidden states. It will be simpler to debug and less if else clauses to check.
        for (desIndex in sequence(length(desRows))){
            focal.edge.length <- phy$edge.length[desRows[desIndex]]
            tipward.age <- as.numeric( rootward.age - focal.edge.length )
            ## Strange rounding errors. A tip age should be zero. This ensures that:
            if(tipward.age < .Machine$double.eps^0.5){
                tipward.age <- 0.0
            }
            node.D <- compD[desNodes[desIndex],]
            node.E <- compE[desNodes[desIndex],]
            ## Call to lsoda that utilizes C code. Requires a lot of inputs. Note that for now we hardcode the NUMELEMENTS arguments. The reason for this is because with lsoda we can only pass a vector of parameters.
            ## Need to make sure that 'node.E' and 'node.D' have the correct quantities.
            yini <- c(node.E, node.D)
            ## yini <- c(E_0=cache$node.E[1], E_1=cache$node.E[2], E_01=cache$node.E[3], D_N0=cache$node.D[1], D_N1=cache$node.D[2], D_N2=cache$node.D[3])
            ## This informs the heights of the tipward and rootward nodes of the branch, so that we can compute the length of the branch.
            times <- c(tipward.age, rootward.age)

            ## The matrix below has the probabilities for each of the states at the tipward node for each of the
            ##     descendants of the focal node. Note that the loop above visited all the branches imediatelly
            ##     descendants of the focal node.
            ## We need these to compute the probability of each of the states at the focal node.
            ## The probability of each of the states at the focal node is equal to the probability that having
            ##     a given state at the focal node will generate tipward probabilities equal to the ones observed
            ##     or estimated.
            prob.subtree.cal.full <- runSilent(yini=yini, times=times, model.vec=model.vec)

            ## NEED TO UPDATE THIS CHECK! NOT ADJUSTED FOR THE MODEL.
            ## It will not break like this, but it is not checking the model correctly.
            
            ## ###### THIS CHECKS TO ENSURE THAT THE INTEGRATION WAS SUCCESSFUL ###########
            if(attributes(prob.subtree.cal.full)$istate[1] < 0){
                return(bad.likelihood)
            }else{
                prob.subtree.cal <- prob.subtree.cal.full[-1,-1]
            }
            ## ############################################################################

            if( any( is.nan( prob.subtree.cal ) ) ){
                return(bad.likelihood)
            }
            ## This is default and cannot change, but if we get a negative probability, discard the results:
            ## Note that this is only testing the Ds part of the ODEs.
            if( any( prob.subtree.cal < 0 ) ){
                return(bad.likelihood)
            }
            ## This can be modified at the input, but if the sum of the D's at the end of a branch are less than some value, then discard the results. A little more stringent than diversitree, but with difficult problems, this stabilizes things immensely.
            if( sum(prob.subtree.cal[13:24]) < ode.eps ){
                return(bad.likelihood)
            }
            
            ## Designating phi here because of its relation to Morlon et al (2011) and using "e" would be confusing:
            ## Here phi is a vector with the results of the integration so far in the E's.
            ## Then 'v' is the D's part.
            phi[v_phy_count,] <- prob.subtree.cal[1:12]
            v[v_phy_count,] <- prob.subtree.cal[13:24]
            v_phy_count <- v_phy_count + 1
        }

        ## The integration above was made from the tipward node to the rootward node.
        ## We computed a probability that the rootward node will be in each of the states for each of the branches.
        ## Of course, the rootward node is a single entity. So we need to combine the probabilities.
        ## Each column of v gives the probability for that state at the rootward node.
        ## An speciation happened at that node to form the observed descendants. This is an observed speciation
        ##     event, so its probability (just like the observed states) is equal to 1.
        ## The events track the infinitesimal moment just after the observed speciation event (i.e., no transitions
        ##     are possible after speciation, because we are tracking a single event.
        
        ## In the case of the endemic states at the nodes it is simple. Only one thing can happen:
        compD[focal,1] <- v[1,1] * v[2,1] * model.vec[1]
        compD[focal,2] <- v[1,2] * v[2,2] * model.vec[2]
        compD[focal,3] <- v[1,3] * v[2,3] * model.vec[3]
        compD[focal,7] <- v[1,7] * v[2,7] * model.vec[7]
        compD[focal,8] <- v[1,8] * v[2,8] * model.vec[8]
        compD[focal,9] <- v[1,9] * v[2,9] * model.vec[9]

        ## When the state at the node is widespread then there are more than one thing possible.
        ## Here we need to enumerate all possible speciation types given that the state at the node was widespread.
        ## state 12A
        compD[focal,4] <- 0.5 * (v[1,4] * v[2,1] + v[1,1] * v[2,4]) * model.vec[1] + 0.5 * (v[1,4] * v[2,2] + v[1,2] * v[2,4]) * model.vec[2] + 0.5 * (v[1,1] * v[2,2] + v[1,2] * v[2,1]) * model.vec[4]
        ## state 23A
        compD[focal,5] <- 0.5 * (v[1,5] * v[2,2] + v[1,2] * v[2,5]) * model.vec[2] + 0.5 * (v[1,5] * v[2,3] + v[1,3] * v[2,5]) * model.vec[3] + 0.5 * (v[1,2] * v[2,3] + v[1,3] * v[2,2]) * model.vec[5]
        ## state 13A
        compD[focal,6] <- 0.5 * (v[1,6] * v[2,1] + v[1,1] * v[2,6]) * model.vec[1] + 0.5 * (v[1,6] * v[2,3] + v[1,3] * v[2,6]) * model.vec[3] + 0.5 * (v[1,1] * v[2,3] + v[1,3] * v[2,1]) * model.vec[6]
        ## state 12B
        compD[focal,10] <- 0.5 * (v[1,10] * v[2,7] + v[1,7] * v[2,10]) * model.vec[7] + 0.5 * (v[1,10] * v[2,8] + v[1,8] * v[2,10]) * model.vec[8] + 0.5 * (v[1,7] * v[2,8] + v[1,8] * v[2,7]) * model.vec[10]
        ## state 23B
        compD[focal,11] <- 0.5 * (v[1,11] * v[2,8] + v[1,8] * v[2,11]) * model.vec[8] + 0.5 * (v[1,11] * v[2,9] + v[1,9] * v[2,11]) * model.vec[9] + 0.5 * (v[1,8] * v[2,9] + v[1,9] * v[2,8]) * model.vec[11]
        ## state 13B
        compD[focal,12] <- 0.5 * (v[1,12] * v[2,7] + v[1,7] * v[2,12]) * model.vec[7] + 0.5 * (v[1,12] * v[2,9] + v[1,9] * v[2,12]) * model.vec[9] + 0.5 * (v[1,7] * v[2,9] + v[1,9] * v[2,7]) * model.vec[12]
        
        ## This is the probability of extinction at the node associated with each of the states.
        compE[focal,] <- phi[1,]
       
        ## #########################
        ## Logcompensation bit for dealing with underflow issues.

        ## Having a real problem here.
        ## The 'tmp' value is something negative.
        ## It is underflow problem. Probabilities cannot be negative. The values for the probabilities are too small.
        tmp <- sum(compD[focal,])
        compD[focal,] <- compD[focal,] / tmp
        ## Keep track of elements for the log compensation:
        logcomp <- c(logcomp, log(tmp))
        
    } ## FINISH FOR LOOP over the nodes of the tree.
    ## Compute quantities for the root node.
    
    root.node <- nb.tip + 1L
    if ( is.na(sum(log(compD[root.node,]))) || is.na(log(sum(1 - compE[root.node,]))) ){
        return(bad.likelihood)
    }else{
        if(root.type == "madfitz" | root.type == "herr_als"){
            if(is.null(root.p)){
                root.p <- compD[root.node,] / sum(compD[root.node,])
                ## Strange here. Why to accept NA as a possible value at the root?
                root.p[ which(is.na(root.p)) ] <- 0
            }
        }
        ## WARNING: THIS WILL CHANGE A BIT WITH THE HIDDEN STATES:
        ## TRY TO FOLLOW JEREMY'S CODE!
        if(condition.on.survival == TRUE){
            if(root.type == "madfitz"){
                ## Here the endemic speciation rates contribute to multiple widespread areas, but this is happening to the classe model too. Should not be the issue.
                ## For hidden rates this needs to be in the same order of the compD matrix.
                lambda <- c(model.vec[1:3], sum(model.vec[c(1,2,4)]), sum(model.vec[c(2,3,5)]), sum(model.vec[c(1,3,6)])
                          , model.vec[7:9], sum(model.vec[c(7,8,10)]), sum(model.vec[c(8,9,11)]), sum(model.vec[c(7,8,12)]) )
                compD[root.node,] <- compD[root.node,] / sum(root.p * lambda * (1 - compE[root.node,])^2)
                ## Corrects for possibility that you have 0/0:
                compD[root.node,which(is.na(compD[root.node,]))] <- 0
                loglik <- log(sum(compD[root.node,] * root.p)) + sum(logcomp)
            }else{
                lambda <- c(model.vec[1:3], sum(model.vec[c(1,2,4)]), sum(model.vec[c(2,3,5)]), sum(model.vec[c(1,3,6)])
                          , model.vec[7:9], sum(model.vec[c(7,8,10)]), sum(model.vec[c(8,9,11)]), sum(model.vec[c(7,8,12)]) )
                compD[root.node,] <- (compD[root.node,] * root.p) / (lambda * (1 - compE[root.node,])^2)
                ## Corrects for possibility that you have 0/0:
                compD[root.node,which(is.na(compD[root.node,]))] <- 0
                loglik <- log(sum(compD[root.node,])) + sum(logcomp)
            }
        }
    }

    return( -1 * loglik )
}

## #################################################################################### 
## Optimization function for multi_rate model:
## This will deal with hidden states from 2 up to 12.
## ####################################################################################

opt_geohisse_3_multi_rate <- function(p, pars, phy, bad.likelihood=10000000, condition.on.survival, root.type, root.p, ode.eps, split.times, nb.tip, nb.node, anc, compD, compE) {
    
    ## Generates the final vector with the appropriate parameter estimates in the right place:
    p.new <- exp(p)
    
    ## Expand the parameter vector into the long form for the C code.
    model.vec <- expand.pars(lik = p.new, long = pars)
    
    logcomp <- c() ## Keep track of the log of the compensation.
    
    ## Start the postorder traversal indexing lists by node number:
    
    ## This is the real core of the likelhood function.
    ## Need to make sure that we have all the objects in the correct format to run these.

    ## Define the run function that will be called at each iteration:
    runSilent <- function(yini, times, model.vec){
        ## This is silencing the integration warnings.
        options(warn = -1)
        on.exit(options(warn = 0))
        capture.output( res <- lsoda(yini, times, func = "geosse_3_areas_multi_rates_derivs"
                                   , model.vec, initfunc="initmod_geosse_3_areas_multi_rates", dllname = "hisse"
                                   , rtol=1e-8, atol=1e-8) )
        return( res )
    }

    ## A constant to be used by the C code.
    NUMELEMENTS <- 336
    
    for (i in seq(from = 1, length.out = nb.node)) {
        ## A vector of all the internal nodes:
        focal <- anc[i]
        desRows <- which(phy$edge[,1]==focal)
        desNodes <- phy$edge[desRows,2]
        ## Note: when the tree has been reordered branching.times are no longer valid. Fortunately, we extract this information in the initial cache setup. Also focal is the rootward node, whereas desNodes represent a vector of all descendant nodes:
        rootward.age <- as.numeric( split.times[which(names(split.times)==focal)] )

        ## Matrix objects to store the results.
        v <- matrix(nrow = length(desRows), ncol = ncol(compD) )
        phi <- matrix(nrow = length(desRows), ncol = ncol(compD) )
        ## A counter to go through the matrix.
        v_phy_count <- 1
        ## BEGIN FOR LOOP over the nodes of the tree.
        ## It would be simple and better to separate this function between the case with and without the hidden states. It will be simpler to debug and less if else clauses to check.
        for (desIndex in sequence(length(desRows))){
            focal.edge.length <- phy$edge.length[desRows[desIndex]]
            tipward.age <- as.numeric( rootward.age - focal.edge.length )
            ## Strange rounding errors. A tip age should be zero. This ensures that:
            if(tipward.age < .Machine$double.eps^0.5){
                tipward.age <- 0.0
            }
            node.D <- compD[desNodes[desIndex],]
            node.E <- compE[desNodes[desIndex],]
            ## Call to lsoda that utilizes C code. Requires a lot of inputs. Note that for now we hardcode the NUMELEMENTS arguments. The reason for this is because with lsoda we can only pass a vector of parameters.
            ## Need to make sure that 'node.E' and 'node.D' have the correct quantities.
            yini <- c(node.E, node.D)
            ## yini <- c(E_0=cache$node.E[1], E_1=cache$node.E[2], E_01=cache$node.E[3], D_N0=cache$node.D[1], D_N1=cache$node.D[2], D_N2=cache$node.D[3])
            ## This informs the heights of the tipward and rootward nodes of the branch, so that we can compute the length of the branch.
            times <- c(tipward.age, rootward.age)

            ## The matrix below has the probabilities for each of the states at the tipward node for each of the
            ##     descendants of the focal node. Note that the loop above visited all the branches imediatelly
            ##     descendants of the focal node.
            ## We need these to compute the probability of each of the states at the focal node.
            ## The probability of each of the states at the focal node is equal to the probability that having
            ##     a given state at the focal node will generate tipward probabilities equal to the ones observed
            ##     or estimated.
            prob.subtree.cal.full <- runSilent(yini=yini, times=times, model.vec=model.vec)

            ## NEED TO UPDATE THIS CHECK! NOT ADJUSTED FOR THE MODEL.
            ## It will not break like this, but it is not checking the model correctly.
            
            ## ###### THIS CHECKS TO ENSURE THAT THE INTEGRATION WAS SUCCESSFUL ###########
            if(attributes(prob.subtree.cal.full)$istate[1] < 0){
                return(bad.likelihood)
            }else{
                prob.subtree.cal <- prob.subtree.cal.full[-1,-1]
            }
            ## ############################################################################

            if( any( is.nan( prob.subtree.cal ) ) ){
                return(bad.likelihood)
            }
            ## This is default and cannot change, but if we get a negative probability, discard the results:
            ## Note that this is only testing the Ds part of the ODEs.
            if( any( prob.subtree.cal < 0 ) ){
                return(bad.likelihood)
            }
            ## This can be modified at the input, but if the sum of the D's at the end of a branch are less than some value, then discard the results. A little more stringent than diversitree, but with difficult problems, this stabilizes things immensely.
            if( sum(prob.subtree.cal[73:144]) < ode.eps ){
                return(bad.likelihood)
            }
            
            ## Designating phi here because of its relation to Morlon et al (2011) and using "e" would be confusing:
            ## Here phi is a vector with the results of the integration so far in the E's.
            ## Then 'v' is the D's part.
            phi[v_phy_count,] <- prob.subtree.cal[1:72]
            v[v_phy_count,] <- prob.subtree.cal[73:144]
            v_phy_count <- v_phy_count + 1
        }

        ## The integration above was made from the tipward node to the rootward node.
        ## We computed a probability that the rootward node will be in each of the states for each of the branches.
        ## Of course, the rootward node is a single entity. So we need to combine the probabilities.
        ## Each column of v gives the probability for that state at the rootward node.
        ## An speciation happened at that node to form the observed descendants. This is an observed speciation
        ##     event, so its probability (just like the observed states) is equal to 1.
        ## The events track the infinitesimal moment just after the observed speciation event (i.e., no transitions
        ##     are possible after speciation, because we are tracking a single event.

        ## Get the integration for each of the states in the model.
        for( i in 0:11 ){
            compD[focal,1+(i*6)] <- v[1,1+(i*6)] * v[2,1+(i*6)] * model.vec[1+(i*6)]
            compD[focal,2+(i*6)] <- v[1,2+(i*6)] * v[2,2+(i*6)] * model.vec[2+(i*6)]
            compD[focal,3+(i*6)] <- v[1,3+(i*6)] * v[2,3+(i*6)] * model.vec[3+(i*6)]

            compD[focal,4+(i*6)] <- 0.5 * (v[1,4+(i*6)] * v[2,1+(i*6)] + v[1,1+(i*6)] * v[2,4+(i*6)]) * model.vec[1+(i*6)]
            + 0.5 * (v[1,4+(i*6)] * v[2,2+(i*6)] + v[1,2+(i*6)] * v[2,4+(i*6)]) * model.vec[2+(i*6)]
            + 0.5 * (v[1,1+(i*6)] * v[2,2+(i*6)] + v[1,2+(i*6)] * v[2,1+(i*6)]) * model.vec[4+(i*6)]
            compD[focal,5+(i*6)] <- 0.5 * (v[1,5+(i*6)] * v[2,2+(i*6)] + v[1,2+(i*6)] * v[2,5+(i*6)]) * model.vec[2+(i*6)]
            + 0.5 * (v[1,5+(i*6)] * v[2,3+(i*6)] + v[1,3+(i*6)] * v[2,5+(i*6)]) * model.vec[3+(i*6)]
            + 0.5 * (v[1,2+(i*6)] * v[2,3+(i*6)] + v[1,3+(i*6)] * v[2,2+(i*6)]) * model.vec[5+(i*6)]
            compD[focal,6+(i*6)] <- 0.5 * (v[1,6+(i*6)] * v[2,1+(i*6)] + v[1,1+(i*6)] * v[2,6+(i*6)]) * model.vec[1+(i*6)]
            + 0.5 * (v[1,6+(i*6)] * v[2,3+(i*6)] + v[1,3+(i*6)] * v[2,6+(i*6)]) * model.vec[3+(i*6)]
            + 0.5 * (v[1,1+(i*6)] * v[2,3+(i*6)] + v[1,3+(i*6)] * v[2,1+(i*6)]) * model.vec[6+(i*6)]
        }
        ## This is the probability of extinction at the node associated with each of the states.
        compE[focal,] <- phi[1,]
       
        ## #########################
        ## Logcompensation bit for dealing with underflow issues.

        ## Having a real problem here.
        ## The 'tmp' value is something negative.
        ## It is underflow problem. Probabilities cannot be negative. The values for the probabilities are too small.
        tmp <- sum(compD[focal,])
        compD[focal,] <- compD[focal,] / tmp
        ## Keep track of elements for the log compensation:
        logcomp <- c(logcomp, log(tmp))
        
    } ## FINISH FOR LOOP over the nodes of the tree.
    ## Compute quantities for the root node.
    
    root.node <- nb.tip + 1L
    if ( is.na(sum(log(compD[root.node,]))) || is.na(log(sum(1 - compE[root.node,]))) ){
        return(bad.likelihood)
    }else{
        if(root.type == "madfitz" | root.type == "herr_als"){
            if(is.null(root.p)){
                root.p <- compD[root.node,] / sum(compD[root.node,])
                ## Strange here. Why to accept NA as a possible value at the root?
                root.p[ which(is.na(root.p)) ] <- 0
            }
        }
        if(condition.on.survival == TRUE){
            if(root.type == "madfitz"){
                ## Here the endemic speciation rates contribute to multiple widespread areas, but this is happening to the classe model too. Should not be the issue.
                ## For hidden rates this needs to be in the same order of the compD matrix.
                lambda <- vector(mode = "numeric")
                for( i in 0:11 ){
                    lambda <- c(lambda, model.vec[1:3+(i*6)], sum(model.vec[c(1,2,4)+(i*6)])
                              , sum(model.vec[c(2,3,5)+(i*6)]), sum(model.vec[c(1,3,6)+(i*6)]) )
                }
                
                compD[root.node,] <- compD[root.node,] / sum(root.p * lambda * (1 - compE[root.node,])^2)
                ## Corrects for possibility that you have 0/0:
                compD[root.node,which(is.na(compD[root.node,]))] <- 0
                loglik <- log(sum(compD[root.node,] * root.p)) + sum(logcomp)
            }else{
                lambda <- vector(mode = "numeric")
                for( i in 0:11 ){
                    lambda <- c(lambda, model.vec[1:3+(i*6)], sum(model.vec[c(1,2,4)+(i*6)])
                              , sum(model.vec[c(2,3,5)+(i*6)]), sum(model.vec[c(1,3,6)+(i*6)]) )
                }
                compD[root.node,] <- (compD[root.node,] * root.p) / (lambda * (1 - compE[root.node,])^2)
                ## Corrects for possibility that you have 0/0:
                compD[root.node,which(is.na(compD[root.node,]))] <- 0
                loglik <- log(sum(compD[root.node,])) + sum(logcomp)
            }
        }
    }

    return( -1 * loglik )
}

## #################################################################################### 
## Some helping functions:
## ####################################################################################

expand.pars <- function(lik, long){
    ## Function to expand the model likelihood pars to the long format with all the parameters of the model.
    out <- rep(0, times = length(long) )
    fit.id <- which(long > 0)
    for( i in 1:length(fit.id)){
        ## The 'long' vectors has integers that are indexes.
        out[ fit.id[i] ] <- lik[ long[ fit.id[i] ] ]
    }
    return( out )
}
