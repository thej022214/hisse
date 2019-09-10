######################################################################################################################################
######################################################################################################################################
### ANCESTRAL STATE RECONSTRUCTION -- calculates marginal probabilities
######################################################################################################################################
######################################################################################################################################

MarginRecon <- function(phy, data, f, pars, hidden.states=TRUE, four.state.null=FALSE, timeslice=NULL, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, aic=NULL, verbose=TRUE, n.cores=NULL){
    if(!is.null(root.p)) {
        root.type="user"
        root.p <- root.p / sum(root.p)
        if(hidden.states ==TRUE & length(root.p)==2){
            root.p <- rep(root.p, 2)
            root.p <- root.p / sum(root.p)
            warning("For hidden states, you need to specify the root.p for all four hidden states. We have adjusted it so that there's equal chance for 0A as 0B, and for 1A as 1B")
        }
    }
    
    if(is.null(n.cores)){
        n.cores=1
    }
    
    if(four.state.null == FALSE){
        phy$node.label = NULL
        data.new <- data.frame(data[,2], data[,2], row.names=data[,1])
        data.new <- data.new[phy$tip.label,]
        model.vec = pars
        #Prerequisites for running the downpass algorithm:
        cache = ParametersToPass(phy, data.new[,1], model.vec, f=f, timeslice=NULL, hidden.states=hidden.states)
        cache$turnover.beta.factor0 = 1 / dbeta(0.1, model.vec[21], model.vec[25])
        cache$turnover.beta.factor1 = 1 / dbeta(0.1, model.vec[22], model.vec[26])
        cache$turnover.beta.factorA = 1 / dbeta(0.1, model.vec[23], model.vec[27])
        cache$turnover.beta.factorB = 1 / dbeta(0.1, model.vec[24], model.vec[28])
        cache$eps.beta.factor0 = 1 / dbeta(0.1, model.vec[29], model.vec[33])
        cache$eps.beta.factor1 = 1 / dbeta(0.1, model.vec[30], model.vec[34])
        cache$eps.beta.factorA = 1 / dbeta(0.1, model.vec[31], model.vec[35])
        cache$eps.beta.factorB = 1 / dbeta(0.1, model.vec[32], model.vec[36])
        
        nb.tip <- length(phy$tip.label)
        nb.node <- phy$Nnode
        if(hidden.states == FALSE){
            nstates=2
        }else{
            nstates=4
        }
        nodes <- unique(phy$edge[,1])
        
        NodeEval <- function(node){
            focal <- node
            marginal.probs.tmp <- c()
            for (j in 1:nstates){
                marginal.probs.tmp <- c(marginal.probs.tmp, DownPass(phy, cache, hidden.states=hidden.states, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, node=focal, state=j))
            }
            best.probs = max(marginal.probs.tmp)
            marginal.probs.rescaled = marginal.probs.tmp - best.probs
            marginal.probs = exp(marginal.probs.rescaled) / sum(exp(marginal.probs.rescaled))
            return(c(node, marginal.probs))
        }
        node.marginals <- mclapply((nb.tip+1):(nb.tip+nb.node), NodeEval, mc.cores=n.cores)
        
        if(hidden.states==TRUE){
            TipEval <- function(tip){
                marginal.probs.tmp <- numeric(4)
                nstates = which(!cache$states[tip,] == 0)
                cache$states.keep = cache$states[tip,]
                for (j in nstates){
                    cache$states[tip,] = 0
                    cache$states[tip,j] = 1
                    marginal.probs.tmp[j] <- DownPass(phy, cache, hidden.states=hidden.states, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, node=NULL, state=j)
                }
                cache$states[tip,] = cache$states.keep
                best.probs = max(marginal.probs.tmp[nstates])
                marginal.probs.rescaled = marginal.probs.tmp[nstates] - best.probs
                marginal.probs <- numeric(4)
                marginal.probs[nstates] = exp(marginal.probs.rescaled) / sum(exp(marginal.probs.rescaled))
                return(c(tip, marginal.probs))
            }
            tip.marginals <- mclapply(1:nb.tip, TipEval, mc.cores=n.cores)
        }
        obj <- NULL
        if(hidden.states == TRUE){
            obj$node.mat <- matrix(unlist(node.marginals), ncol = 5, byrow = TRUE)
            obj$tip.mat = matrix(unlist(tip.marginals), ncol = 5, byrow = TRUE)
            raw.rates <- ParameterTransform(x=model.vec[1:4], y=model.vec[5:8])
            rates.mat <- matrix(0, 5, 4)
            rates.mat[1,] <- model.vec[1:4]
            rates.mat[2,] <- raw.rates[1:4] - raw.rates[5:8]
            rates.mat[3,] <- raw.rates[1:4]
            rates.mat[4,] <- model.vec[5:8]
            rates.mat[5,] <- raw.rates[5:8]
            rownames(rates.mat) <- c("turnover", "net.div", "speciation", "extinction.fraction", "extinction")
            colnames(rates.mat) <- c("0A", "1A", "0B", "1B")
            colnames(obj$node.mat) <- colnames(obj$tip.mat)  <- c("id", "0A", "1A", "0B", "1B")
        }else{
            obj$node.mat <- matrix(unlist(node.marginals), ncol = 3, byrow = TRUE)
            obj$tip.mat = cbind(1:Ntip(phy), cache$states)
            raw.rates <- ParameterTransform(x=model.vec[1:4], y=model.vec[5:8])
            rates.mat <- matrix(0, 5, 2)
            rates.mat[1,] <- model.vec[1:2]
            rates.mat[2,] <- raw.rates[1:2] - raw.rates[5:6]
            rates.mat[3,] <- raw.rates[1:2]
            rates.mat[4,] <- model.vec[5:6]
            rates.mat[5,] <- raw.rates[5:6]
            rownames(rates.mat) <- c("turnover", "net.div", "speciation", "extinction.fraction", "extinction")
            colnames(rates.mat) <- c("0", "1")
            colnames(obj$node.mat) <- colnames(obj$tip.mat) <- c("id", "0", "1")
        }
        obj$rates.mat = rates.mat
        phy$node.label = apply(obj$node.mat[,2:dim(obj$node.mat)[2]], 1, which.max)
        obj$phy = phy
    }else{
        phy$node.label = NULL
        data.new <- data.frame(data[,2], data[,2], row.names=data[,1])
        data.new <- data.new[phy$tip.label,]
        model.vec = pars
        #Prerequisites for running the downpass algorithm:
        cache = ParametersToPassNull(phy, data.new[,1], model.vec, f=f)
        
        nb.tip <- length(phy$tip.label)
        nb.node <- phy$Nnode
        nstates=8
        nodes <- unique(phy$edge[,1])
        
        NodeEval <- function(node){
            focal <- node
            marginal.probs.tmp <- c()
            for (j in 1:nstates){
                marginal.probs.tmp <- c(marginal.probs.tmp, DownPassNull(phy, cache, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, node=focal, state=j))
            }
            best.probs = max(marginal.probs.tmp)
            marginal.probs.rescaled = marginal.probs.tmp - best.probs
            marginal.probs = exp(marginal.probs.rescaled) / sum(exp(marginal.probs.rescaled))
            return(c(node, marginal.probs))
        }
        node.marginals <- mclapply((nb.tip+1):(nb.tip+nb.node), NodeEval, mc.cores=n.cores)
        
        TipEval <- function(tip){
            marginal.probs.tmp <- numeric(8)
            nstates = which(!cache$states[tip,] == 0)
            cache$states.keep = cache$states[tip,]
            for (j in nstates){
                cache$states[tip,] = 0
                cache$states[tip,j] = 1
                marginal.probs.tmp[j] <- DownPassNull(phy, cache, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, node=NULL, state=j)
            }
            cache$states[tip,] = cache$states.keep
            best.probs = max(marginal.probs.tmp[nstates])
            marginal.probs.rescaled = marginal.probs.tmp[nstates] - best.probs
            marginal.probs <- numeric(8)
            marginal.probs[nstates] = exp(marginal.probs.rescaled) / sum(exp(marginal.probs.rescaled))
            return(c(tip, marginal.probs))
        }
        tip.marginals <- mclapply(1:nb.tip, TipEval, mc.cores=n.cores)
        obj <- NULL
        obj$node.mat <- matrix(unlist(node.marginals), ncol = 9, byrow = TRUE)
        obj$tip.mat = matrix(unlist(tip.marginals), ncol = 9, byrow = TRUE)
        raw.rates <- ParameterTransform(x=model.vec[1:8], y=model.vec[9:16])
        rates.mat <- matrix(0, 5, 8)
        rates.mat[1,] <- model.vec[1:8]
        rates.mat[2,] <- raw.rates[1:8] - raw.rates[9:16]
        rates.mat[3,] <- raw.rates[1:8]
        rates.mat[4,] <- model.vec[9:16]
        rates.mat[5,] <- raw.rates[9:16]
        rownames(rates.mat) <- c("turnover", "net.div", "speciation", "extinction.fraction", "extinction")
        colnames(rates.mat) <- c("0A", "0B", "0C", "0D","1A", "1B", "1C", "1D")
        colnames(obj$node.mat) <- colnames(obj$tip.mat)  <- c("id", "0A", "0B", "0C", "0D","1A", "1B", "1C", "1D")
        obj$rates.mat = rates.mat
        phy$node.label = apply(obj$node.mat[,2:dim(obj$node.mat)[2]], 1, which.max)
        obj$phy = phy
    }
    
    if(!is.null(aic)){
        obj$aic = aic
    }
    class(obj) = "hisse.states"
    return(obj)
}



######################################################################################################################################
######################################################################################################################################
### ANCESTRAL STATE RECONSTRUCTION for GeoSSE -- calculates marginal probabilities for our set of GeoSSE models
######################################################################################################################################
######################################################################################################################################

MarginReconGeoSSE <- function(phy, data, f, pars, hidden.areas=TRUE, assume.cladogenetic=TRUE, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, aic=NULL, verbose=TRUE, n.cores=NULL){
    
    if( !is.null(phy$node.label) ) phy$node.label <- NULL
    
    if(!is.null(root.p)) {
        root.type="user"
        root.p <- root.p / sum(root.p)
        if(hidden.areas == TRUE & length(root.p) == 3){
            ## What if length(root.p) is 4 or 5? This code will easily break or allow wrong entries.
            root.p <- rep(root.p, 3)
            root.p <- root.p / sum(root.p)
            warning("For hidden states, you need to specify the root.p for all four hidden states. We have adjusted it so that there's equal chance among all hidden states.")
        }
    }
    
    data.new <- data.frame(data[,2], data[,2], row.names=data[,1])
    data.new <- data.new[phy$tip.label,]
    model.vec = pars
    
    #Prerequisites for running the downpass algorithm:
    if(assume.cladogenetic == TRUE){
        cache = ParametersToPassGeoHiSSE(phy, data.new[,1], model.vec, f=f, hidden.states=hidden.areas)
    }else{
        cache = ParametersToPassMuSSE(phy, data.new[,1], model.vec, f=f, hidden.states=hidden.areas)
    }
    
    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
    if(hidden.areas == FALSE){
        nstates = 3
    }else{
        nstates = 15
    }
    nodes <- unique(phy$edge[,1])
    
    if(is.null(n.cores)){
        n.cores=1
    }
    
    NodeEval <- function(node){
        focal <- node
        marginal.probs.tmp <- c()
        for (j in 1:nstates){
            if(assume.cladogenetic == TRUE){
                marginal.probs.tmp <- c(marginal.probs.tmp, DownPassGeoHisse(phy, cache, hidden.states=hidden.areas, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, node=focal, state=j))
            }else{
                marginal.probs.tmp <- c(marginal.probs.tmp, DownPassMusse(phy, cache, hidden.states=hidden.areas, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, node=focal, state=j))
            }
        }
        best.probs = max(marginal.probs.tmp)
        marginal.probs.rescaled = marginal.probs.tmp - best.probs
        marginal.probs = exp(marginal.probs.rescaled) / sum(exp(marginal.probs.rescaled))
        return(c(node, marginal.probs))
    }
    node.marginals <- mclapply((nb.tip+1):(nb.tip+nb.node), NodeEval, mc.cores=n.cores)
    
    if(hidden.areas==TRUE){
        TipEval <- function(tip){
            marginal.probs.tmp <- numeric(4)
            nstates = which(!cache$states[tip,] == 0)
            cache$states.keep = cache$states[tip,]
            for (j in nstates){
                cache$states[tip,] = 0
                cache$states[tip,j] = 1
                if(assume.cladogenetic == TRUE){
                    marginal.probs.tmp[j] <- DownPassGeoHisse(phy, cache, hidden.states=hidden.areas, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, node=NULL, state=j)
                }else{
                    marginal.probs.tmp[j] <- DownPassMusse(phy, cache, hidden.states=hidden.areas, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, node=NULL, state=j)
                }
            }
            cache$states[tip,] = cache$states.keep
            best.probs = max(marginal.probs.tmp[nstates])
            marginal.probs.rescaled = marginal.probs.tmp[nstates] - best.probs
            marginal.probs <- numeric(15)
            marginal.probs[nstates] = exp(marginal.probs.rescaled) / sum(exp(marginal.probs.rescaled))
            return(c(tip, marginal.probs))
        }
        tip.marginals <- mclapply(1:nb.tip, TipEval, mc.cores=n.cores)
    }
    
    obj <- NULL
    
    if(hidden.areas == TRUE){
        obj$node.mat <- matrix(unlist(node.marginals), ncol = 15+1, byrow = TRUE)
        obj$tip.mat = matrix(unlist(tip.marginals), ncol = 15+1, byrow = TRUE)
        rates.mat <- matrix(0, 2, 15)
        if(assume.cladogenetic == TRUE){
            rates.mat[1,] <- model.vec[c(1, 2, 3, 24, 25, 26, 47, 48, 49, 70, 71, 72, 93, 94, 95)]
            rates.mat[2,] <- c(model.vec[c(4, 5)], 0, model.vec[c(27, 28)], 0, model.vec[c(50, 51)], 0, model.vec[c(73, 74)], 0, model.vec[c(96, 97)], 0)
        }else{
            rates.mat[1,] <- model.vec[c(1, 2, 3, 25, 26, 27, 49, 50, 51, 73, 74, 75, 97, 98, 99)]
            rates.mat[2,] <- model.vec[c(4, 5, 6, 28, 29, 30, 52, 53, 54, 76, 77, 78, 100, 101, 102)]
        }
        rownames(rates.mat) <- c("speciation", "extinction")
        colnames(rates.mat) <- c("0A", "1A", "01A", "0B", "1B", "01B", "0C", "1C", "01C", "0D", "1D", "01D", "0E", "1E", "01E")
        colnames(obj$node.mat) <- colnames(obj$tip.mat)  <- c("id", "0A", "1A", "01A", "0B", "1B", "01B", "0C", "1C", "01C", "0D", "1D", "01D", "0E", "1E", "01E")
        rates.mat <- ParameterTransformGeoSSE(rates.mat, assume.cladogenetic=assume.cladogenetic)
    }else{
        obj$node.mat <- matrix(unlist(node.marginals), ncol = 3+1, byrow = TRUE)
        obj$tip.mat = cbind(1:Ntip(phy), cache$states)
        rates.mat <- matrix(0, 2, 3)
        rates.mat[1,] <- model.vec[c(1, 2, 3)]
        if(assume.cladogenetic == TRUE){
            rates.mat[2,] <- c(model.vec[c(4, 5)], 0)
        }else{
            rates.mat[2,] <- model.vec[c(4, 5, 6)]
        }
        rownames(rates.mat) <- c("speciation", "extinction")
        colnames(rates.mat) <- c("0", "1", "01")
        colnames(obj$node.mat) <- colnames(obj$tip.mat)  <- c("id", "0", "1", "01")
        rates.mat <- ParameterTransformGeoSSE(rates.mat, assume.cladogenetic=assume.cladogenetic)
    }
    obj$rates.mat = rates.mat
    phy$node.label = apply(obj$node.mat[,2:dim(obj$node.mat)[2]], 1, which.max)
    obj$phy = phy
    
    if(!is.null(aic)){
        obj$aic = aic
    }
    
    class(obj) = "hisse.geosse.states"
    return(obj)
}



######################################################################################################################################
######################################################################################################################################
### ANCESTRAL STATE RECONSTRUCTION for newer, faster HiSSE -- calculates marginal probabilities for our set of HiSSE models
######################################################################################################################################
######################################################################################################################################

MarginReconfHiSSE <- function(phy, data, f, pars, hidden.states=1, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, aic=NULL, get.tips.only=FALSE, verbose=TRUE, n.cores=NULL, dt.threads=1){
    
    if( !is.null(phy$node.label) ) phy$node.label <- NULL
    
    if(!is.null(root.p)) {
        root.type="user"
        root.p <- root.p / sum(root.p)
        if(hidden.states > 1 & length(root.p)==2){
            root.p <- rep(root.p, hidden.states)
            root.p <- root.p / sum(root.p)
            warning("For hidden states, you need to specify the root.p for all possible hidden states. We have adjusted it so that there's equal chance for each of the specified hidden states")
        }
    }
    
    setDTthreads(threads=dt.threads)
    
    model.vec = pars
    
    # Some new prerequisites #
    data.new <- data.frame(data[,2], data[,2], row.names=data[,1])
    data.new <- data.new[phy$tip.label,]
    gen <- FindGenerations(phy)
    dat.tab <- OrganizeDataHiSSE(data=data.new, phy=phy, f=f, hidden.states=TRUE)
    nb.tip <- Ntip(phy)
    nb.node <- phy$Nnode
    ##########################
    
    ### Ughy McUgherson. This is a must in order to pass CRAN checks: http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
    DesNode = NULL
    ##########################
    
    cache <- ParametersToPassfHiSSE(model.vec=model.vec, hidden.states=TRUE, nb.tip=nb.tip, nb.node=nb.node, bad.likelihood=exp(-500), ode.eps=0)
    
    nstates <- 8
    nstates.to.eval <- 2 * hidden.states
    nstates.not.eval <- 8 - nstates.to.eval
    nodes <- unique(phy$edge[,1])
    
    if(is.null(n.cores)){
        n.cores=1
    }
    
    NodeEval <- function(node){
        focal <- node
        marginal.probs.tmp <- c()
        for (j in 1:nstates.to.eval){
            marginal.probs.tmp <- c(marginal.probs.tmp, DownPassHiSSE(dat.tab=dat.tab, gen=gen, cache=cache, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, node=focal, state=j, fix.type="fix"))
        }
        marginal.probs.tmp <- c(marginal.probs.tmp, rep(log(cache$bad.likelihood)^13, nstates.not.eval))
        best.probs <- max(marginal.probs.tmp)
        marginal.probs.rescaled <- marginal.probs.tmp - best.probs
        marginal.probs <- exp(marginal.probs.rescaled) / sum(exp(marginal.probs.rescaled))
        return(c(node, marginal.probs))
    }
    
    if(get.tips.only == FALSE){
        cat(paste("Calculating marginal probabilities for ", length(nodes), " internal nodes...", sep=""), "\n")
        obj <- NULL
        node.marginals <- mclapply((nb.tip+1):(nb.tip+nb.node), NodeEval, mc.cores=n.cores)
        obj$node.mat <- matrix(unlist(node.marginals), ncol = 8+1, byrow = TRUE)
        colnames(obj$node.mat) <- c("id", "(0A)", "(1A)", "(0B)", "(1B)", "(0C)", "(1C)", "(0D)", "(1D)")
        phy$node.label <- apply(obj$node.mat[,2:dim(obj$node.mat)[2]], 1, which.max)
    }else{
        cat("Calculating marginal probabilities for internal nodes is turned off...", "\n")
        obj <- NULL
    }
    
    dat.tab <- OrganizeDataHiSSE(data=data.new, phy=phy, f=f, hidden.states=TRUE)
    TipEval <- function(tip){
        setkey(dat.tab, DesNode)
        marginal.probs.tmp <- numeric(8)
        nstates = which(!dat.tab[tip,7:14] == 0)
        cache$states.keep <- as.data.frame(dat.tab[tip,7:14])
        for (j in nstates[1:hidden.states]){
            cache$to.change <- cache$states.keep
            tmp.state <- 1 * c(cache$to.change[1,j])
            cache$to.change[1,] <- 0
            cache$to.change[1,j] <- tmp.state
            for (k in 1:dim(cache$to.change)[2]){
                dat.tab[tip, paste("compD", k, sep="_") := cache$to.change[,k]]
            }
            marginal.probs.tmp[j] <- DownPassHiSSE(dat.tab=dat.tab, gen=gen, cache=cache, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, node=NULL, state=j)
        }
        for (k in 1:dim(cache$to.change)[2]){
            dat.tab[tip, paste("compD", k, sep="_") := cache$states.keep[,k]]
        }
        best.probs <- max(marginal.probs.tmp[nstates[1:hidden.states]])
        marginal.probs.rescaled <- marginal.probs.tmp[nstates[1:hidden.states]] - best.probs
        marginal.probs <- numeric(8)
        marginal.probs[nstates[1:hidden.states]] <- exp(marginal.probs.rescaled) / sum(exp(marginal.probs.rescaled))
        return(c(tip, marginal.probs))
    }
    
    if(hidden.states>1){
        cat(paste("Finished. Calculating marginal probabilities for ", nb.tip, " tips...", sep=""), "\n")
        tip.marginals <- mclapply(1:nb.tip, TipEval, mc.cores=n.cores)
        obj$tip.mat <- matrix(unlist(tip.marginals), ncol = 8+1, byrow = TRUE)
    }else{
        obj$tip.mat <- matrix(0, ncol = 8+1, nrow = nb.tip)
        obj$tip.mat[,1] <- 1:nb.tip
        setkey(dat.tab, DesNode)
        obj$tip.mat[,2:3] <- matrix(unlist(dat.tab[1:nb.tip,7:8]), ncol = 2, byrow = FALSE)
    }
    
    cat("Done.","\n")
    
    colnames(obj$tip.mat)  <- c("id", "(0A)", "(1A)", "(0B)", "(1B)", "(0C)", "(1C)", "(0D)", "(1D)")
    rates.mat <- matrix(0, 2, 8)
    rates.mat[1,] <- model.vec[c(1:2, 13:14, 25:26, 37:38)]
    rates.mat[2,] <- model.vec[c(3:4, 15:16, 27:28, 39:40)]
    rownames(rates.mat) <- c("turnover", "extinction.fraction")
    colnames(rates.mat) <- c("(0A)", "(1A)", "(0B)", "(1B)", "(0C)", "(1C)", "(0D)", "(1D)")
    rates.mat <- ParameterTransformfHiSSE(rates.mat)
    obj$rates.mat = rates.mat
    obj$phy = phy
    
    if(!is.null(aic)){
        obj$aic = aic
    }
    
    class(obj) = "hisse.states"
    return(obj)
}



######################################################################################################################################
######################################################################################################################################
### ANCESTRAL STATE RECONSTRUCTION for MuHiSSE -- calculates marginal probabilities for our set of MuHiSSE models
######################################################################################################################################
######################################################################################################################################

MarginReconMuHiSSE <- function(phy, data, f, pars, hidden.states=2, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, aic=NULL, get.tips.only=FALSE, verbose=TRUE, n.cores=NULL, dt.threads=1){
    
    if( !is.null(phy$node.label) ) phy$node.label <- NULL
    
    if(!is.null(root.p)) {
        root.type="user"
        root.p <- root.p / sum(root.p)
        if(hidden.states > 1 & length(root.p)==4){
            root.p <- rep(root.p, hidden.states)
            root.p <- root.p / sum(root.p)
            warning("For hidden states, you need to specify the root.p for all possible hidden states. We have adjusted it so that there's equal chance for each of the specified hidden states")
        }
    }
    
    setDTthreads(threads=dt.threads)

    model.vec = pars
    
    # Some new prerequisites #
    data.new <- data.frame(data[,2], data[,3], row.names=data[,1])
    data.new <- data.new[phy$tip.label,]
    gen <- FindGenerations(phy)
    dat.tab <- OrganizeData(data=data.new, phy=phy, f=f, hidden.states=TRUE)
    nb.tip <- Ntip(phy)
    nb.node <- phy$Nnode
    ##########################
    
    ### Ughy McUgherson. This is a must in order to pass CRAN checks: http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
    DesNode = NULL
    ##########################
    
    cache <- ParametersToPassMuHiSSE(model.vec=model.vec, hidden.states=TRUE, nb.tip=nb.tip, nb.node=nb.node, bad.likelihood=exp(-500), ode.eps=0)
    
    nstates <- 32
    nstates.to.eval <- 4 * hidden.states
    nstates.not.eval <- 32 - nstates.to.eval
    nodes <- unique(phy$edge[,1])
    
    if(is.null(n.cores)){
        n.cores=1
    }
    
    NodeEval <- function(node){
        focal <- node
        marginal.probs.tmp <- c()
        for (j in 1:nstates.to.eval){
            marginal.probs.tmp <- c(marginal.probs.tmp, DownPassMuHisse(dat.tab=dat.tab, gen=gen, cache=cache, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, node=focal, state=j, fix.type="fix"))
        }
        marginal.probs.tmp <- c(marginal.probs.tmp, rep(log(cache$bad.likelihood)^13, nstates.not.eval))
        best.probs <- max(marginal.probs.tmp)
        marginal.probs.rescaled <- marginal.probs.tmp - best.probs
        marginal.probs <- exp(marginal.probs.rescaled) / sum(exp(marginal.probs.rescaled))
        return(c(node, marginal.probs))
    }
    
    if(get.tips.only == FALSE){
        cat(paste("Calculating marginal probabilities for ", length(nodes), " internal nodes...", sep=""), "\n")
        obj <- NULL
        node.marginals <- mclapply((nb.tip+1):(nb.tip+nb.node), NodeEval, mc.cores=n.cores)
        obj$node.mat <- matrix(unlist(node.marginals), ncol = 32+1, byrow = TRUE)
        colnames(obj$node.mat) <- c("id", "(00A)","(01A)","(10A)","(11A)", "(00B)","(01B)","(10B)","(11B)", "(00C)","(01C)","(10C)","(11C)", "(00D)","(01D)","(10D)","(11D)", "(00E)","(01E)","(10E)","(11E)", "(00F)","(01F)","(10F)","(11F)", "(00G)","(01G)","(10G)","(11G)", "(00H)","(01H)","(10H)","(11H)")
        phy$node.label <- apply(obj$node.mat[,2:dim(obj$node.mat)[2]], 1, which.max)
    }else{
        cat("Calculating marginal probabilities for internal nodes is turned off...", "\n")
        obj <- NULL
    }
    
    dat.tab <- OrganizeData(data=data.new, phy=phy, f=f, hidden.states=TRUE)
    TipEval <- function(tip){
        setkey(dat.tab, DesNode)
        marginal.probs.tmp <- numeric(32)
        nstates = which(!dat.tab[tip,7:38] == 0)
        cache$states.keep <- as.data.frame(dat.tab[tip,7:38])
        for (j in nstates[1:hidden.states]){
            cache$to.change <- cache$states.keep
            tmp.state <- 1 * c(cache$to.change[1,j])
            cache$to.change[1,] <- 0
            cache$to.change[1,j] <- tmp.state
            for (k in 1:dim(cache$to.change)[2]){
                dat.tab[tip, paste("compD", k, sep="_") := cache$to.change[,k]]
            }
            marginal.probs.tmp[j] <- DownPassMuHisse(dat.tab=dat.tab, gen=gen, cache=cache, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, node=NULL, state=j)
        }
        for (k in 1:dim(cache$to.change)[2]){
            dat.tab[tip, paste("compD", k, sep="_") := cache$states.keep[,k]]
        }
        best.probs <- max(marginal.probs.tmp[nstates[1:hidden.states]])
        marginal.probs.rescaled <- marginal.probs.tmp[nstates[1:hidden.states]] - best.probs
        marginal.probs <- numeric(32)
        marginal.probs[nstates[1:hidden.states]] <- exp(marginal.probs.rescaled) / sum(exp(marginal.probs.rescaled))
        return(c(tip, marginal.probs))
    }
    
    if(hidden.states>1){
        cat(paste("Finished. Calculating marginal probabilities for ", nb.tip, " tips...", sep=""), "\n")
        tip.marginals <- mclapply(1:nb.tip, TipEval, mc.cores=n.cores)
        obj$tip.mat <- matrix(unlist(tip.marginals), ncol = 32+1, byrow = TRUE)
    }else{
        obj$tip.mat <- matrix(0, ncol = 32+1, nrow = nb.tip)
        obj$tip.mat[,1] <- 1:nb.tip
        setkey(dat.tab, DesNode)
        obj$tip.mat[,2:5] <- matrix(unlist(dat.tab[1:nb.tip,7:10]), ncol = 4, byrow = FALSE)
    }
    
    cat("Done.","\n")
    
    colnames(obj$tip.mat)  <- c("id", "(00A)","(01A)","(10A)","(11A)", "(00B)","(01B)","(10B)","(11B)", "(00C)","(01C)","(10C)","(11C)", "(00D)","(01D)","(10D)","(11D)", "(00E)","(01E)","(10E)","(11E)", "(00F)","(01F)","(10F)","(11F)", "(00G)","(01G)","(10G)","(11G)", "(00H)","(01H)","(10H)","(11H)")
    rates.mat <- matrix(0, 2, 32)
    rates.mat[1,] <- model.vec[c(1:4, 49:52, 97:100, 145:148, 193:196, 241:244, 289:292, 337:340)]
    rates.mat[2,] <- model.vec[c(5:8, 53:56, 101:104, 149:152, 197:200, 245:248, 293:296, 341:344)]
    rownames(rates.mat) <- c("turnover", "extinction.fraction")
    colnames(rates.mat) <- c("(00A)","(01A)","(10A)","(11A)", "(00B)","(01B)","(10B)","(11B)", "(00C)","(01C)","(10C)","(11C)", "(00D)","(01D)","(10D)","(11D)", "(00E)","(01E)","(10E)","(11E)", "(00F)","(01F)","(10F)","(11F)", "(00G)","(01G)","(10G)","(11G)", "(00H)","(01H)","(10H)","(11H)")
    rates.mat <- ParameterTransformMuHiSSE(rates.mat)
    obj$rates.mat = rates.mat
    obj$phy = phy
    
    if(!is.null(aic)){
        obj$aic = aic
    }
    
    class(obj) = "muhisse.states"
    return(obj)
}



######################################################################################################################################
######################################################################################################################################
### ANCESTRAL STATE RECONSTRUCTION for fast GeoHiSSE -- calculates marginal probabilities for our set of GeoHiSSE models
######################################################################################################################################
######################################################################################################################################

MarginReconfGeoSSE <- function(phy, data, f, pars, hidden.areas=2, assume.cladogenetic=TRUE, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, aic=NULL, get.tips.only=FALSE, verbose=TRUE, n.cores=NULL, dt.threads=1){
    
    if( !is.null(phy$node.label) ) phy$node.label <- NULL
    
    if(!is.null(root.p)) {
        root.type="user"
        root.p <- root.p / sum(root.p)
        if(hidden.areas ==TRUE & length(root.p)==3){
            root.p <- rep(root.p, 3)
            root.p <- root.p / sum(root.p)
            warning("For hidden states, you need to specify the root.p for all possible hidden states. We have adjusted it so that there's equal chance for 0A as 0B, and for 1A as 1B")
        }
    }
    
    setDTthreads(threads=dt.threads)

    model.vec = pars
    
    # Some new prerequisites #
    data.new <- data.frame(data[,2], data[,2], row.names=data[,1])
    data.new <- data.new[phy$tip.label,]
    gen <- FindGenerations(phy)
    dat.tab <- OrganizeDataGeo(data=data.new[,1], phy=phy, f=f, hidden.states=TRUE)
    nb.tip <- Ntip(phy)
    nb.node <- phy$Nnode
    ### Ughy McUgherson. This is a must in order to pass CRAN checks: http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
    DesNode = NULL
    ##########################
    
    cache <- ParametersToPassGeoHiSSEfast(model.vec, hidden.states=TRUE, assume.cladogenetic=assume.cladogenetic, nb.tip=nb.tip, nb.node=nb.node, bad.likelihood=exp(-500), ode.eps=0)
    
    nstates = 30
    nstates.to.eval <- 3 * hidden.areas
    nstates.not.eval <- 30 - nstates.to.eval
    
    nodes <- unique(phy$edge[,1])
    
    if(is.null(n.cores)){
        n.cores=1
    }
    
    NodeEval <- function(node){
        focal <- node
        marginal.probs.tmp <- c()
        for (j in 1:nstates.to.eval){
            marginal.probs.tmp <- c(marginal.probs.tmp, DownPassGeoHissefast(dat.tab=dat.tab, gen=gen, cache=cache, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, node=focal, state=j))
        }
        marginal.probs.tmp <- c(marginal.probs.tmp, rep(log(cache$bad.likelihood)^13, nstates.not.eval))
        best.probs = max(marginal.probs.tmp)
        marginal.probs.rescaled = marginal.probs.tmp - best.probs
        marginal.probs = exp(marginal.probs.rescaled) / sum(exp(marginal.probs.rescaled))
        return(c(node, marginal.probs))
    }
    
    if(get.tips.only == FALSE){
        cat(paste("Calculating marginal probabilities for ", length(nodes), " internal nodes...", sep=""), "\n")
        obj <- NULL
        node.marginals <- mclapply((nb.tip+1):(nb.tip+nb.node), NodeEval, mc.cores=n.cores)
        obj$node.mat <- matrix(unlist(node.marginals), ncol = 30+1, byrow = TRUE)
        colnames(obj$node.mat) <- c("id", "(00A)", "(11A)", "(01A)", "(00B)", "(11B)", "(01B)", "(00C)", "(11C)", "(01C)", "(00D)", "(11D)", "(01D)", "(00E)", "(11E)", "(01E)", "(00F)", "(11F)", "(01F)", "(00G)", "(11G)", "(01G)", "(00H)", "(11H)", "(01H)", "(00I)", "(11I)", "(01I)", "(00J)", "(11J)", "(01J)")
        phy$node.label <- apply(obj$node.mat[,2:dim(obj$node.mat)[2]], 1, which.max)
    }else{
        cat("Calculating marginal probabilities for internal nodes is turned off...", "\n")
        obj <- NULL
    }
    
    dat.tab <- OrganizeDataGeo(data=data.new[,1], phy=phy, f=f, hidden.states=TRUE)
    TipEval <- function(tip){
        setkey(dat.tab, DesNode)
        marginal.probs.tmp <- numeric(30)
        nstates = which(!dat.tab[tip,7:36] == 0)
        cache$states.keep <- as.data.frame(dat.tab[tip,7:36])
        for (j in nstates[1:hidden.areas]){
            cache$to.change <- cache$states.keep
            tmp.state <- 1 * c(cache$to.change[1,j])
            cache$to.change[1,] <- 0
            cache$to.change[1,j] <- tmp.state
            for (k in 1:dim(cache$to.change)[2]){
                dat.tab[tip, paste("compD", k, sep="_") := cache$to.change[,k]]
            }
            cache$to.change[1,j] <- tmp.state
            marginal.probs.tmp[j] <- DownPassGeoHissefast(dat.tab=dat.tab, gen=gen, cache=cache, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, node=NULL, state=j)
        }
        for (k in 1:dim(cache$to.change)[2]){
            dat.tab[tip, paste("compD", k, sep="_") := cache$states.keep[,k]]
        }
        best.probs = max(marginal.probs.tmp[nstates[1:hidden.areas]])
        marginal.probs.rescaled = marginal.probs.tmp[nstates[1:hidden.areas]] - best.probs
        marginal.probs <- numeric(30)
        marginal.probs[nstates[1:hidden.areas]] = exp(marginal.probs.rescaled) / sum(exp(marginal.probs.rescaled))
        return(c(tip, marginal.probs))
    }
    
    if(hidden.areas>1){
        cat(paste("Finished. Calculating marginal probabilities for ", nb.tip, " tips...", sep=""), "\n")
        tip.marginals <- mclapply(1:nb.tip, TipEval, mc.cores=n.cores)
        obj$tip.mat = matrix(unlist(tip.marginals), ncol = 30+1, byrow = TRUE)
    }else{
        obj$tip.mat <- matrix(0, ncol = 30+1, nrow = nb.tip)
        obj$tip.mat[,1] <- 1:nb.tip
        setkey(dat.tab, DesNode)
        obj$tip.mat[,2:4] <- matrix(unlist(dat.tab[1:nb.tip,7:9]), ncol = 3, byrow = FALSE)
    }

    cat("Done.","\n")
    rates.mat <- matrix(0, 2, 30)
    rates.mat[1,] <- model.vec[c(1:3, 39:41, 77:79, 115:117, 153:155, 191:193, 229:231, 267:269, 305:307, 343:345)]
    rates.mat[2,] <- c(model.vec[c(4:5)],0, model.vec[c(42:43)],0,  model.vec[c(80:81)],0, model.vec[c(118:119)],0, model.vec[c(156:157)],0, model.vec[c(194:195)],0, model.vec[c(232:233)],0, model.vec[c(270:271)],0, model.vec[c(308:309)],0, model.vec[c(346:347)],0)
    rownames(rates.mat) <- c("turnover", "extinction.fraction")
    colnames(rates.mat) <- c("(00A)", "(11A)", "(01A)", "(00B)", "(11B)", "(01B)", "(00C)", "(11C)", "(01C)", "(00D)", "(11D)", "(01D)", "(00E)", "(11E)", "(01E)", "(00F)", "(11F)", "(01F)", "(00G)", "(11G)", "(01G)", "(00H)", "(11H)", "(01H)", "(00I)", "(11I)", "(01I)", "(00J)", "(11J)", "(01J)")
    
    colnames(obj$tip.mat) <- c("id", "(00A)", "(11A)", "(01A)", "(00B)", "(11B)", "(01B)", "(00C)", "(11C)", "(01C)", "(00D)", "(11D)", "(01D)", "(00E)", "(11E)", "(01E)", "(00F)", "(11F)", "(01F)", "(00G)", "(11G)", "(01G)", "(00H)", "(11H)", "(01H)", "(00I)", "(11I)", "(01I)", "(00J)", "(11J)", "(01J)")
    rates.mat <- ParameterTransformfGeoSSE(rates.mat)
    
    obj$rates.mat = rates.mat
    obj$phy = phy
    
    
    if(!is.null(aic)){
        obj$aic = aic
    }
    
    class(obj) = "hisse.geosse.states"
    return(obj)
}



######################################################################################################################################
######################################################################################################################################
### ANCESTRAL STATE RECONSTRUCTION for fast MiSSE -- calculates marginal probabilities for our set of MiSSE models
######################################################################################################################################
######################################################################################################################################


MarginReconMiSSE <- function(phy, f, pars, hidden.states=2, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, aic=NULL, get.tips.only=FALSE, verbose=TRUE, n.cores=NULL, dt.threads=1){
    
    if( !is.null(phy$node.label) ) phy$node.label <- NULL
    
    if(!is.null(root.p)) {
        root.type="user"
        root.p <- root.p / sum(root.p)
        if(hidden.states ==TRUE & length(root.p)==2){
            root.p <- rep(root.p, 2)
            root.p <- root.p / sum(root.p)
            warning("For hidden states, you need to specify the root.p for all possible hidden states. We have adjusted it so that there's equal chance for each of the specified hidden states")
        }
    }
    
    setDTthreads(threads=dt.threads)

    model.vec = pars
    
    # Some new prerequisites #
    gen <- FindGenerations(phy)
    dat.tab <- OrganizeDataMiSSE(phy=phy, f=f, hidden.states=hidden.states)
    nb.tip <- Ntip(phy)
    nb.node <- phy$Nnode
    ### Ughy McUgherson. This is a must in order to pass CRAN checks: http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
    DesNode = NULL
    ##########################
    
    cache <- ParametersToPassMiSSE(model.vec=model.vec, hidden.states=hidden.states, nb.tip=nb.tip, nb.node=nb.node, bad.likelihood=exp(-500), ode.eps=0)
    
    nstates = 26
    nstates.to.eval <- hidden.states
    nstates.not.eval <- 26 - nstates.to.eval
    nodes <- unique(phy$edge[,1])
    
    if(is.null(n.cores)){
        n.cores=1
    }
    
    NodeEval <- function(node){
        focal <- node
        marginal.probs.tmp <- c()
        for (j in 1:nstates.to.eval){
            marginal.probs.tmp <- c(marginal.probs.tmp, DownPassMisse(dat.tab=dat.tab, cache=cache, gen=gen, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, node=focal, state=j))
        }
        marginal.probs.tmp <- c(marginal.probs.tmp, rep(log(cache$bad.likelihood)^13, nstates.not.eval))
        best.probs = max(marginal.probs.tmp)
        marginal.probs.rescaled = marginal.probs.tmp - best.probs
        marginal.probs = exp(marginal.probs.rescaled) / sum(exp(marginal.probs.rescaled))
        return(c(node, marginal.probs))
    }
    
    if(get.tips.only == FALSE){
        cat(paste("Calculating marginal probabilities for ", length(nodes), " internal nodes...", sep=""), "\n")
        obj <- NULL
        node.marginals <- mclapply((nb.tip+1):(nb.tip+nb.node), NodeEval, mc.cores=n.cores)
        obj$node.mat <- matrix(unlist(node.marginals), ncol = 26+1, byrow = TRUE)
        colnames(obj$node.mat) <- c("id", "(0A)", "(0B)", "(0C)", "(0D)", "(0E)", "(0F)", "(0G)", "(0H)", "(0I)", "(0J)", "(0K)", "(0L)", "(0M)", "(0N)", "(0O)", "(0P)", "(0Q)", "(0R)", "(0S)", "(0T)", "(0U)", "(0V)", "(0W)", "(0X)", "(0Y)", "(0Z)")
        phy$node.label = apply(obj$node.mat[,2:dim(obj$node.mat)[2]], 1, which.max)
    }else{
        cat("Calculating marginal probabilities for internal nodes is turned off...", "\n")
        obj <- NULL
    }
    
    dat.tab <- OrganizeDataMiSSE(phy=phy, f=f, hidden.states=hidden.states)
    TipEval <- function(tip){
        setkey(dat.tab, DesNode)
        marginal.probs.tmp <- numeric(2)
        nstates = which(!dat.tab[tip,7:32] == 0)
        cache$states.keep <- as.data.frame(dat.tab[tip,7:32])
        for (j in nstates){
            cache$to.change <- cache$states.keep
            tmp.state <- 1 * c(cache$to.change[1,j])
            cache$to.change[1,] <- 0
            cache$to.change[1,j] <- tmp.state
            for (k in 1:dim(cache$to.change)[2]){
                dat.tab[tip, paste("compD", k, sep="_") := cache$to.change[,k]]
            }
            marginal.probs.tmp[j] <- DownPassMisse(dat.tab=dat.tab, cache=cache, gen=gen, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, node=NULL, state=j)
        }
        for (k in 1:dim(cache$to.change)[2]){
            dat.tab[tip, paste("compD", k, sep="_") := cache$states.keep[,k]]
        }
        best.probs = max(marginal.probs.tmp[nstates])
        marginal.probs.rescaled = marginal.probs.tmp[nstates] - best.probs
        marginal.probs <- numeric(26)
        marginal.probs[nstates] = exp(marginal.probs.rescaled) / sum(exp(marginal.probs.rescaled))
        return(c(tip, marginal.probs))
    }
    cat(paste("Finished. Calculating marginal probabilities for ", nb.tip, " tips...", sep=""), "\n")
    tip.marginals <- mclapply(1:nb.tip, TipEval, mc.cores=n.cores)
    obj$tip.mat = matrix(unlist(tip.marginals), ncol = 26+1, byrow = TRUE)
    colnames(obj$tip.mat)  <- c("id", "(0A)", "(0B)", "(0C)", "(0D)", "(0E)", "(0F)", "(0G)", "(0H)", "(0I)", "(0J)", "(0K)", "(0L)", "(0M)", "(0N)", "(0O)", "(0P)", "(0Q)", "(0R)", "(0S)", "(0T)", "(0U)", "(0V)", "(0W)", "(0X)", "(0Y)", "(0Z)")
    
    cat("Done.","\n")
    rates.mat <- matrix(0, 2, 26)
    index.vector <- 1:52
    rates.mat[1,] <- model.vec[index.vector %% 2 == 1][-27]
    rates.mat[2,] <- model.vec[index.vector %% 2 == 0]
    rownames(rates.mat) <- c("turnover", "extinction.fraction")
    colnames(rates.mat) <- c("(0A)", "(0B)", "(0C)", "(0D)", "(0E)", "(0F)", "(0G)", "(0H)", "(0I)", "(0J)", "(0K)", "(0L)", "(0M)", "(0N)", "(0O)", "(0P)", "(0Q)", "(0R)", "(0S)", "(0T)", "(0U)", "(0V)", "(0W)", "(0X)", "(0Y)", "(0Z)")
    rates.mat <- ParameterTransformMiSSE(rates.mat)
    obj$rates.mat = rates.mat
    obj$phy = phy
    
    if(!is.null(aic)){
        obj$aic = aic
    }
    
    class(obj) = "misse.states"
    return(obj)
}



######################################################################################################################################
######################################################################################################################################
### UTILITY FUNCTIONS
######################################################################################################################################
######################################################################################################################################

ParameterTransform <- function(x, y){
    speciation <- x / (1+y)
    extinction <- (x*y) / (1+y)
    return(c(speciation, extinction))
}


ParameterTransformGeoSSE <- function(x, assume.cladogenetic=TRUE){
    ## Also need to add the extirpation bit as well to the rate matrix. -- especially if we separate it. It is an event that "represents" extinction of a range. So should it count?
    if(assume.cladogenetic == TRUE){
        if(dim(x)[2] == 15){
            rates.mat <- matrix(0, 3, 15)
            rownames(rates.mat) <- c("turnover", "net.div", "extinction.fraction")
            colnames(rates.mat) <- c("0A", "1A", "01A", "0B", "1B", "01B", "0C", "1C", "01C", "0D", "1D", "01D", "0E", "1E", "01E")
            rates.mat[1,] <- x[1,] + x[2,]
            rates.mat[2,] <- x[1,] - x[2,]
            rates.mat[3,] <- x[2,] / x[1,]
            for(widespread.index in c(3,6,9,12,15)){
                rates.mat[1,widespread.index] <- sum(x[1,c(widespread.index-2,widespread.index-1,widespread.index)])
                rates.mat[2,widespread.index] <- sum(x[2,c(widespread.index-2,widespread.index-1,widespread.index)])
                rates.mat[3,widespread.index] <- 0
                #rates.mat[1,widespread.index] <- sum(x[1,c(widespread.index-2,widespread.index-2,widespread.index)]) + sum(x[2,c(widespread.index-2,widespread.index-2,widespread.index)])
                #rates.mat[2,widespread.index] <- sum(x[1,c(widespread.index-2,widespread.index-2,widespread.index)]) - sum(x[2,c(widespread.index-2,widespread.index-2,widespread.index)])
                #rates.mat[3,widespread.index] <- sum(x[2,c(widespread.index-2,widespread.index-2,widespread.index)]) / sum(x[1,c(widespread.index-2,widespread.index-2,widespread.index)])
            }
            rates.mat[3,is.na(rates.mat[3,])] = 0
        }else{
            rates.mat <- matrix(0, 3, 3)
            rownames(rates.mat) <- c("turnover", "net.div", "extinction.fraction")
            colnames(rates.mat) <- c("0", "1", "01")
            rates.mat[1,] <- x[1,] + x[2,]
            rates.mat[2,] <- x[1,]
            rates.mat[3,] <- x[2,] / x[1,]
            #rates.mat[1,3] <- sum(x[1,c(1,2,3)]) + sum(x[2,c(1,2,3)])
            #rates.mat[2,3] <- sum(x[1,c(1,2,3)]) - sum(x[2,c(1,2,3)])
            #rates.mat[3,3] <- sum(x[2,c(1,2,3)]) / sum(x[1,c(1,2,3)])
            rates.mat[3,is.na(rates.mat[3,])] = 0
            rates.mat[1,3] <- sum(x[1,c(1,2,3)])
            rates.mat[2,3] <- sum(x[1,c(1,2,3)])
            rates.mat[3,3] <- 0
        }
    }else{
        if(dim(x)[2] == 15){
            rates.mat <- matrix(0, 3, 15)
            rownames(rates.mat) <- c("turnover", "net.div", "extinction.fraction")
            colnames(rates.mat) <- c("0A", "1A", "01A", "0B", "1B", "01B", "0C", "1C", "01C", "0D", "1D", "01D", "0E", "1E", "01E")
            rates.mat[1,] <- x[1,] + x[2,]
            rates.mat[2,] <- x[1,] - x[2,]
            rates.mat[3,] <- x[2,] / x[1,]
            rates.mat[3,is.na(rates.mat[3,])] = 0
        }else{
            rates.mat <- matrix(0, 3, 3)
            rownames(rates.mat) <- c("turnover", "net.div", "extinction.fraction")
            colnames(rates.mat) <- c("0", "1", "01")
            rates.mat[1,] <- x[1,] + x[2,]
            rates.mat[2,] <- x[1,] - x[2,]
            rates.mat[3,] <- x[2,] / x[1,]
            rates.mat[3,is.na(rates.mat[3,])] = 0
        }
    }
    rates.mat <- rbind(x, rates.mat)
    return(rates.mat)
}


ParameterTransformfHiSSE <- function(x){
    if(dim(x)[2] == 8){
        rates.mat <- matrix(0, 3, 8)
        rownames(rates.mat) <- c("net.div", "speciation", "extinction")
        colnames(rates.mat) <- c("(0A)", "(1A)", "(0B)", "(1B)", "(0C)", "(1C)", "(0D)", "(1D)")
        rates.mat[1,] <- (x[1,] / (1 + x[2,])) - ((x[1,] * x[2,]) / (1 + x[2,]))
        rates.mat[2,] <- x[1,] / (1 + x[2,])
        rates.mat[3,] <- (x[1,] * x[2,]) / (1 + x[2,])
        rates.mat[3,is.na(rates.mat[3,])] = 0
    }else{
        rates.mat <- matrix(0, 3, 2)
        rownames(rates.mat) <- c("net.div", "speciation", "extinction")
        colnames(rates.mat) <- c("(0A)","(1A)")
        rates.mat[1,] <- (x[1,] / (1 + x[2,])) - ((x[1,] * x[2,]) / (1 + x[2,]))
        rates.mat[2,] <- x[1,] / (1 + x[2,])
        rates.mat[3,] <- (x[1,] * x[2,]) / (1 + x[2,])
        rates.mat[3,is.na(rates.mat[3,])] = 0
    }
    rates.mat <- rbind(x, rates.mat)
    return(rates.mat)
}


ParameterTransformMuHiSSE <- function(x){
    if(dim(x)[2] == 32){
        rates.mat <- matrix(0, 3, 32)
        rownames(rates.mat) <- c("net.div", "speciation", "extinction")
        colnames(rates.mat) <- c("(00A)","(01A)","(10A)","(11A)", "(00B)","(01B)","(10B)","(11B)", "(00C)","(01C)","(10C)","(11C)", "(00D)","(01D)","(10D)","(11D)", "(00E)","(01E)","(10E)","(11E)", "(00F)","(01F)","(10F)","(11F)", "(00G)","(01G)","(10G)","(11G)", "(00H)","(01H)","(10H)","(11H)")
        rates.mat[1,] <- (x[1,] / (1 + x[2,])) - ((x[1,] * x[2,]) / (1 + x[2,]))
        rates.mat[2,] <- x[1,] / (1 + x[2,])
        rates.mat[3,] <- (x[1,] * x[2,]) / (1 + x[2,])
        rates.mat[3,is.na(rates.mat[3,])] = 0
    }else{
        rates.mat <- matrix(0, 3, 4)
        rownames(rates.mat) <- c("net.div", "speciation", "extinction")
        colnames(rates.mat) <- c("(00A)","(01A)","(10A)","(11A)")
        rates.mat[1,] <- (x[1,] / (1 + x[2,])) - ((x[1,] * x[2,]) / (1 + x[2,]))
        rates.mat[2,] <- x[1,] / (1 + x[2,])
        rates.mat[3,] <- (x[1,] * x[2,]) / (1 + x[2,])
        rates.mat[3,is.na(rates.mat[3,])] = 0
    }
    rates.mat <- rbind(x, rates.mat)
    return(rates.mat)
}


ParameterTransformfGeoSSE <- function(x, assume.cladogenetic=TRUE){
    ## Also need to add the extirpation bit as well to the rate matrix. -- especially if we separate it. It is an event that "represents" extinction of a range. So should it count?
    if(assume.cladogenetic == TRUE){
        if(dim(x)[2] == 30){
            rates.mat <- matrix(0, 3, 30)
            rownames(rates.mat) <- c("speciation", "extinction", "net.div")
            colnames(rates.mat) <- c("00A", "11A", "01A", "00B", "11B", "01B", "00C", "11C", "01C", "00D", "11D", "01D", "00E", "11E", "01E", "00F", "11F", "01F", "00G", "11G", "01G", "00H", "11H", "01H", "00I", "11I", "01I", "00J", "11J", "01J")
            rates.mat[1,] <- x[1,] / (1 + x[2,])
            rates.mat[2,] <- (x[1,] * x[2,]) / (1 + x[2,])
            rates.mat[3,] <- rates.mat[1,] - rates.mat[2,]
            for(widespread.index in c(3,6,9,12,15,18,21,24,27,30)){
                rates.mat[1,widespread.index] <- x[1,widespread.index] - rates.mat[1,widespread.index-1] - rates.mat[1, widespread.index-2]
                rates.mat[2,widespread.index] <- 0
                rates.mat[3,widespread.index] <- sum(x[1,c(widespread.index-2,widespread.index-1,widespread.index)])
            }
            rates.mat[3,is.na(rates.mat[3,])] = 0
        }else{
            rates.mat <- matrix(0, 3, 3)
            rownames(rates.mat) <- c("speciation", "extinction", "net.div")
            colnames(rates.mat) <- c("(00A)", "(11A)", "(01A)")
            rates.mat[1,] <- x[1,] / (1 + x[2,])
            rates.mat[2,] <- (x[1,] * x[2,]) / (1 + x[2,])
            rates.mat[3,] <- rates.mat[1,] - rates.mat[2,]
            
            rates.mat[3,is.na(rates.mat[3,])] = 0
            rates.mat[1,3] <- x[1,3] - rates.mat[1,2] - rates.mat[1, 1]
            rates.mat[2,3] <- 0
            rates.mat[3,3] <- sum(x[1,c(1,2,3)])
        }
    }else{
        if(dim(x)[2] == 30){
            rates.mat <- matrix(0, 3, 30)
            rownames(rates.mat) <- c("speciation", "extinction", "net.div")
            colnames(rates.mat) <- c("00A", "11A", "01A", "00B", "11B", "01B", "00C", "11C", "01C", "00D", "11D", "01D", "00E", "11E", "01E", "00F", "11F", "01F", "00G", "11G", "01G", "00H", "11H", "01H", "00I", "11I", "01I", "00J", "11J", "01J")
            rates.mat[1,] <- x[1,] / (1 + x[2,])
            rates.mat[2,] <- (x[1,] * x[2,]) / (1 + x[2,])
            rates.mat[3,] <- x[1,] - x[2,]
            rates.mat[3,is.na(rates.mat[3,])] = 0
        }else{
            rates.mat <- matrix(0, 3, 3)
            rownames(rates.mat) <- c("speciation", "extinction", "net.div")
            colnames(rates.mat) <- c("(00A)", "(11A)", "(01A)")
            rates.mat[1,] <- x[1,] / (1 + x[2,])
            rates.mat[2,] <- (x[1,] * x[2,]) / (1 + x[2,])
            rates.mat[3,] <- x[1,] - x[2,]
            rates.mat[3,is.na(rates.mat[3,])] = 0
        }
    }
    rates.mat <- rbind(x, rates.mat)
    return(rates.mat)
}


ParameterTransformMiSSE <- function(x){
    rates.mat <- matrix(0, 3, 26)
    rownames(rates.mat) <- c("net.div", "speciation", "extinction")
    colnames(rates.mat) <- c("(0A)", "(0B)", "(0C)", "(0D)", "(0E)", "(0F)", "(0G)", "(0H)", "(0I)", "(0J)", "(0K)", "(0L)", "(0M)", "(0N)", "(0O)", "(0P)", "(0Q)", "(0R)", "(0S)", "(0T)", "(0U)", "(0V)", "(0W)", "(0X)", "(0Y)", "(0Z)")
    rates.mat[1,] <- (x[1,] / (1 + x[2,])) - ((x[1,] * x[2,]) / (1 + x[2,]))
    rates.mat[2,] <- x[1,] / (1 + x[2,])
    rates.mat[3,] <- (x[1,] * x[2,]) / (1 + x[2,])
    rates.mat[3,is.na(rates.mat[3,])] = 0
    rates.mat <- rbind(x, rates.mat)
    return(rates.mat)
}


ParameterTransform <- function(x, y){
    speciation <- x / (1+y)
    extinction <- (x*y) / (1+y)
    return(c(speciation, extinction))
}


print.hisse.states <- function(x,...){
    print(x$phy)
}


print.hisse.geosse.states <- function(x,...){
    print(x$phy)
}

print.muhisse.states <- function(x,...){
    print(x$phy)
}

print.misse.states <- function(x,...){
    print(x$phy)
}
