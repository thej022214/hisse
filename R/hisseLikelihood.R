## Exporting the likelihood function for the HiSSE model.
## The parameters of the model will be a single named vector.
## The names for the vector will be produced by a make likelihood function.
## The make function will return a likelihood function with the data already attached as well as the correct named vector of the parameters. The function will also return the suggested starting state for the search.

makeHiSSELikelihood <- function(phy, data, hidden.states = TRUE, null4 = FALSE, f=c(1,1), condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, ode.eps=0){

    ## Things to correct:
    ## The "turnover.beta.factor" are set to constant values. [ See DevOptimize and DevOptimizeNull ]
    ## The "DevOptimize" and "DevOptimizeNull" functions are different. Need to specify if the model is a null model. [ Need to test if I get the same likelihood just by constraining the parameters to be the same. ]
    ## Need to add a parameter to tell if the model should include the hidden state or not. [ The DevOptimize part of the model is different in the presence or absence of the hidden state. ]
    ## The size of the parameters vector and the names will likely need to be adjusted depending on these options.
    
    if(!is.null(root.p)) {
        root.p <- root.p / sum(root.p)	
        if(length(root.p) == 2){
            root.p <- rep(root.p, 2)
            root.p <- root.p / sum(root.p)	
            warning("For hidden states, you need to specify the root.p for all four hidden states. We have adjusted it so that there's equal chance for 0A as 0B, and for 1A as 1B")
        }
    }
    
    if( !root.type %in% c("madfitz", "herr_als") ){
        stop("Check that you specified a proper root.type option. Options are 'madfitz' or 'herr_als'. See help for more details.", call.=FALSE)
    }

    ## Erase the node labels in the phylogeny.
    phy$node.label <- NULL

    ## Correct format for the data.
    data.new <- data.frame(data[,2], data[,2], row.names=data[,1])
    data.new <- data.new[phy$tip.label,]

    ## This is used to scale starting values to account for sampling:
    if(length(f) == 2){
        ## The sampling of the whole tree:
        samp.freq.tree <- Ntip(phy) / sum(table(data.new[,1]) / f)
    }else{
        stop("Sampling frequency needs to be a numeric vector with length 2."
           , call.=FALSE)
    }
    
    ## Create the 'pars' vector for a full model:
    ## This will always have the same length and order of parameters:
    ## Only the first 20 entries will matter for the user.
    ## NOTE: Values for elements 1 to 56 need to be equal to 1.
    pars.full <- c(1:20, rep(1, times = 36))

    if( !null4 ){
        ## This option covers most of the models. From the BiSSE-like model to the full HiSSE model.
        ## This option can also be used to model a Null HiSSE with 2 states. [The NULL for the BiSSE model.]
    if( hidden.states ){
        ## The names for the 20 first parameters in the model.
        ## These are to be returned to the user.
        names.pars <- paste("turnover", c("0A","1A","0B","1B"), sep=".")
        names.pars <- c(names.pars, paste("eps", c("0A","1A","0B","1B"), sep=".") )
        names.pars <- c(names.pars, "q.1A-0A", "q.0B-0A", "q.1B-0A"
                      , "q.0A-1A", "q.0B-1A", "q.1B-1A"
                      , "q.0A-0B", "q.1A-0B", "q.1B-0B"
                      , "q.0A-1B", "q.1A-1B", "q.0B-1B")
        ## Make an empty numeric vector to return to the user.
        initial.pars <- setNames(object = rep(0.0, times = length(names.pars))
                               , nm = names.pars)

        ## Make the likelihood function for the model given the data:
        loglik <- function(p){
            ## Because this is the full model we don't need to adjust the formating for the p vector.
            model.vec <- c(p, rep(1, 36))
            cache <- ParametersToPass(phy=phy, data=data.new[,1], model.vec=model.vec, f=f
                                    , timeslice=NULL, hidden.states=TRUE)
            ## Set the value for constant not used in the model.
            cache$turnover.beta.factor0 <- 1
            cache$turnover.beta.factor1 <- 1
            cache$turnover.beta.factorA <- 1
            cache$turnover.beta.factorB <- 1
            cache$eps.beta.factor0 <- 1
            cache$eps.beta.factor1 <- 1
            cache$eps.beta.factorA <- 1
            cache$eps.beta.factorB <- 1

            logl <- DownPass(phy=phy, cache=cache, hidden.states=TRUE
                           , condition.on.survival=condition.on.survival
                           , root.type=root.type, root.p=root.p, ode.eps=ode.eps)
            ## This function is returning the true log likelihood for the model.
            return(logl)
        }
        
    } else{
        ## The names for the 20 first parameters in the model.
        ## These are to be returned to the user.
        names.pars <- paste("turnover", c("0","1"), sep=".")
        names.pars <- c(names.pars, paste("eps", c("0","1"), sep=".") )
        ## These are positions 1 and 4 of the q matrix vector part.
        names.pars <- c(names.pars, "q.10", "q.01")
        ## Make an empty numeric vector to return to the user.
        initial.pars <- setNames(object = rep(0.0, times = length(names.pars))
                               , nm = names.pars)
        
        ## Make the likelihood function for the model given the data:
        loglik <- function(p){
            ## In this case we don't have the hidden states, so the model is reduced.
            p.full <- c(p[1], p[2], 0, 0, p[3], p[4], 0, 0
                      , p[5], 0, 0, p[6], 0, 0, 0, 0, 0, 0, 0, 0)
            model.vec <- c(p.full, rep(1, 36))
            cache <- ParametersToPass(phy=phy, data=data.new[,1], model.vec=model.vec, f=f
                                    , timeslice=NULL, hidden.states=FALSE)
            ## Set the value for constant not used in the model.
            cache$turnover.beta.factor0 <- 1
            cache$turnover.beta.factor1 <- 1
            cache$turnover.beta.factorA <- 1
            cache$turnover.beta.factorB <- 1
            cache$eps.beta.factor0 <- 1
            cache$eps.beta.factor1 <- 1
            cache$eps.beta.factorA <- 1
            cache$eps.beta.factorB <- 1

            logl <- DownPass(phy=phy, cache=cache, hidden.states=FALSE
                           , condition.on.survival=condition.on.survival
                           , root.type=root.type, root.p=root.p, ode.eps=ode.eps)
            ## This function is returning the true log likelihood for the model.
            return(logl)
        }
        
    }
        return( list(log.lik = loglik, pars = initial.pars) )
    } else{
        ## This is special for the Null 4 model. This is the null model for the full HiSSE model.

        ## This model has a special format for the transition matrix. Need to inform the user what is the order for the parameters. For this I will return a guide matrix. This matrix is the same as used in the null.4 model and will serve only for reference. This need to be clear in the help page of the function.
        sub.mat1 <- sub.mat2 <- TransMatMaker(hidden.states=TRUE)
        sub.mat3 <- sub.mat4 <- matrix(NA, 4,4)
        diag(sub.mat3) <- diag(sub.mat4) <- 0.0
        trans.mat <-rbind(cbind(sub.mat1, sub.mat3),cbind(sub.mat4,sub.mat2))
        trans.mat[!is.na(trans.mat)] <- paste0("q.", 1:32)
        rownames(trans.mat) <- c("(0A)","(0B)","(0C)","(0D)","(1A)","(1B)","(1C)","(1D)")
        colnames(trans.mat) <- c("(0A)","(0B)","(0C)","(0D)","(1A)","(1B)","(1C)","(1D)")

        ## Create the vector with the name of the parameters:
        div.vec <- c("turnover.A", "turnover.B", "turnover.C", "turnover.D",
                     "eps.A", "eps.B", "eps.C", "eps.D")
        names.pars <- c(div.vec, paste0("q.", 1:32) )
        initial.pars <- setNames( rep(0.0, times = length(names.pars) ), names.pars )

        loglik <- function(p){
            ## Prepare the parameter vector
            ## Need to set the diversification parameters for 0 and 1.
            model.vec <- c(p[1], p[1], p[2], p[2], p[3], p[3], p[4], p[4],
                           p[5], p[5], p[6], p[6], p[7], p[7], p[8], p[8],
                           p[9:length(p)])
            cache <- ParametersToPassNull(phy=phy, data=data.new[,1], model.vec, f=f)
            logl <- DownPassNull(phy=phy, cache=cache, root.type=root.type
                               , condition.on.survival=condition.on.survival
                               , root.p=root.p, ode.eps=ode.eps)
            return(logl)
        }
        
        return( list(log.lik = loglik, pars = initial.pars, trans.mat.guide = trans.mat) )
    }
    
}
