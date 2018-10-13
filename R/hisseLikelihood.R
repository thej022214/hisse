## Exporting the likelihood function for the HiSSE model.
## The parameters of the model will be a single named vector.
## The names for the vector will be produced by a make likelihood function.
## The make function will return a likelihood function with the data already attached as well as the correct named vector of the parameters. The function will also return the suggested starting state for the search.

makeHiSSELikelihood <- function(phy, data, f=c(1,1), condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, ode.eps=0){
    if(!is.null(root.p)) {
        root.type="user"
        root.p <- root.p / sum(root.p)	
        if(length(root.p) == 2){
            root.p <- rep(root.p, 2)
            root.p <- root.p / sum(root.p)	
            warning("For hidden states, you need to specify the root.p for all four hidden states. We have adjusted it so that there's equal chance for 0A as 0B, and for 1A as 1B")
        }
    }
    
    if(!root.type == "madfitz" & !root.type == "equal" & !root.type == "user"){
        stop("Check that you specified a proper root.type option. Options are 'madfitz', 'equal', or 'user'.", call.=FALSE)
    }

    ## Erase the node labels in the phylogeny.
    phy$node.label <- NULL

    ## Create the 'pars' vector for a full model:
    ## This will always have the same length and order of parameters:
    ## Only the first 20 entries will matter for the user.
    pars.full <- c(1:20, rep(21, times = 36))

    ## The names for the 20 first parameters in the model.
    ## These are to be returned to the user.
    names.pars <- paste("turnover", c("0A","1A","0B","1B"), sep=".")
    names.pars <- c(names.pars, paste("eps", c("0A","1A","0B","1B"), sep=".") )
    names.pars <- c(names.pars, "q.1A-0A", "q.0B-0A", "q.1B-0A"
                  , "q.0A-1A", "q.0B-1A", "q.1B-1A"
                  , "q.0A-0B", "q.1A-0B", "q.1B-0B"
                  , "q.0A-1B", "q.1A-1B", "q.0B-1B")

    data.new <- data.frame(data[,2], data[,2], row.names=data[,1])
    data.new <- data.new[phy$tip.label,]

    ## This is used to scale starting values to account for sampling:
    if(length(f) == 2){
        samp.freq.tree <- Ntip(phy) / sum(table(data.new[,1]) / f)
    }else{
        stop("Sampling frequency needs to be a numeric vector with length 2.")
    }

    ## Prepare suggested starting point for the model.
    cat("Generating starting point suggestion...", "\n")
    init.pars <- starting.point.generator(phy, 2, samp.freq.tree, yule=FALSE)
    names(init.pars) <- NULL
    init.eps <- init.pars[3]/init.pars[1]
    if(init.eps == 0){ ## This can have unexpected behavior.
        init.eps <- 1e-6
    }
    def.set.pars <- c(rep(log(init.pars[1]+init.pars[3]), 4), rep(log(init.eps),4), rep(log(init.pars[5]), 12), rep(log(1), 36))
    np.sequence <- 1:20 ## The full model!
    ## ip will have the initial values for the important 20 parameters.
    ## But note that the likelihood function depends on more than that.
    ip <- sapply(np.sequence, function(x) def.set.pars[which(pars.full == np.sequence[x])[1]])
    ## Feed the initial parameters to the parameter vector:
    initial.pars <- setNames(object = ip, nm = names.pars)

    ## Make the likelihood function for the model given the data:
    loglik <- function(p){
        lik <- DevOptimize(p, pars=pars.full, phy=phy, data=data.new[,1], f=f, hidden.states=TRUE
                         , condition.on.survival=condition.on.survival, root.type=root.type
                         , root.p=root.p, timeslice=NULL, np=20, ode.eps=ode.eps)
        ## "DevOptimize" function returns -lik because of minimization.
        ## This returns the natural loglik value.
        return( -lik )
    }

    return( list(loglik = loglik, pars = initial.pars) )
}
