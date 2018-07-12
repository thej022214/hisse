## Functions associated with simulation of a Hidden areas GeoSSE model using the ClaSSE model. GeoSSE is a subset of a ClaSSE model.

SimulateGeoHiSSE <- function(pars, hidden.areas=1, x0="0A", max.taxa=Inf, max.time=Inf, add.jumps=FALSE, add.extinction=FALSE, include.extinct=FALSE, return.GeoHiSSE_pars=FALSE, override.safeties=FALSE){
    if(return.GeoHiSSE_pars){
        if( !add.jumps & !add.extinction ){
            row.names.mod.matrix <- c("s01", "s0", "s1", "x0", "x1", "d0", "d1")
        }
        if( add.jumps & !add.extinction ){
            row.names.mod.matrix <- c("s01", "s0", "s1", "x0", "x1", "d0", "d1", "jd0", "jd1")
        }
        if( !add.jumps & add.extinction ){
            row.names.mod.matrix <- c("s01", "s0", "s1", "x0", "x1", "d0", "d1", "x*0", "x*1")
        }
        if( add.jumps & add.extinction ){
            row.names.mod.matrix <- c("s01", "s0", "s1", "x0", "x1", "d0", "d1", "jd0", "jd1", "x*0", "x*1")
        }
        mod.matrix <- matrix(data=0, nrow=length(row.names.mod.matrix), ncol=hidden.areas+1)
        rownames( mod.matrix ) <- row.names.mod.matrix
        colnames( mod.matrix ) <- LETTERS[1:(hidden.areas+1)]
        qmat01 <- matrix(data=0, ncol=hidden.areas+1, nrow=hidden.areas+1)
        diag(qmat01) <- NA
        trans.names <- paste0("01", LETTERS[1:(hidden.areas+1)])
        rownames(qmat01) <- trans.names
        colnames(qmat01) <- trans.names
        qmat0 <- matrix(data=0, ncol=hidden.areas+1, nrow=hidden.areas+1)
        diag(qmat0) <- NA
        trans.names <- paste0("0", LETTERS[1:(hidden.areas+1)])
        rownames(qmat0) <- trans.names
        colnames(qmat0) <- trans.names
        qmat1 <- matrix(data=0, ncol=hidden.areas+1, nrow=hidden.areas+1)
        diag(qmat1) <- NA
        trans.names <- paste0("1", LETTERS[1:(hidden.areas+1)])
        rownames(qmat1) <- trans.names
        colnames(qmat1) <- trans.names
        par.list <- list( model.pars = mod.matrix, q.01 = qmat01, q.0 = qmat0, q.1 = qmat1 )
        class(par.list) <- append(class(par.list),"GeoHiSSE_pars")
        return( par.list )
    }

    ## Make a series of checks:
    if( is.infinite(max.taxa) & is.infinite(max.time)){
        if( override.safeties ) warning("WARNING: Simulation running without stopping criteria defined. \n")
        if( !override.safeties ) stop("No stopping criteria have been informed.")
    }
    ## The parameter list is too specific. Better to use the custom class.
    if( !inherits(pars, "GeoHiSSE_pars") ) stop("Argument 'pars' needs to be of class 'GeoHiSSE_pars'. Please check argument 'return.GeoHiSSE_pars'.")
    ## If there are no hidden areas in the model, there is a chance the vector of parameters is a vector instead of a matrix.
    if( !is.matrix(pars$model.pars) ) stop("The parameter 'pars$model.pars' need to be a matrix.")
    ## Check if jumps or extinction parameters were provided and, if positive, check if the correct option were selected.
    if( "jd0" %in% rownames( pars$model.pars ) & !add.jumps ) stop("Detected jump dispersal parameters. Please set 'add.jumps=TRUE'.")
    if( "jd1" %in% rownames( pars$model.pars ) & !add.jumps ) stop("Detected jump dispersal parameters. Please set 'add.jumps=TRUE'.")
    if( "x*0" %in% rownames( pars$model.pars ) & !add.extinction ) stop("Detected extra extinction parameters. Please set 'add.extinction=TRUE'.")
    if( "x*1" %in% rownames( pars$model.pars ) & !add.extinction ) stop("Detected extra extinction parameters. Please set 'add.extinction=TRUE'.")

    ## Generate the key for the translation of the parameters:
    model.pars.vec <- c( pars$model.pars )
    names(model.pars.vec) <- paste0(rownames( pars$model.pars ), rep(colnames(pars$model.pars), each=nrow(pars$model.pars)))

    ## Get the transtions for the hidden areas.
    ## Have a feeling this could be more elegant. But, well...
    tr.size <- ((hidden.areas+1)^2)-(hidden.areas+1)
    vec.01 <- vector("numeric", length=tr.size)
    vec.0 <- vector("numeric", length=tr.size)
    vec.1 <- vector("numeric", length=tr.size)
    nm.01 <- vector("character", length=tr.size)
    nm.0 <- vector("character", length=tr.size)
    nm.1 <- vector("character", length=tr.size)
    count <- 1 ## keep the count for the loop.
    areas.letters <- LETTERS[1:(hidden.areas+1)]
    for( i in 1:(hidden.areas+1) ){
        for( j in 1:(hidden.areas+1) ){
            if( i == j ) next
            vec.01[count] <- pars$q.01[i,j]
            nm.01[count] <- paste0("q01", areas.letters[i], areas.letters[j])
            vec.0[count] <- pars$q.0[i,j]
            nm.0[count] <- paste0("q0", areas.letters[i], areas.letters[j])
            vec.1[count] <- pars$q.1[i,j]
            nm.1[count] <- paste0("q1", areas.letters[i], areas.letters[j])
            count <- count + 1
        }
    }
    names(vec.01) <- nm.01
    names(vec.0) <- nm.0
    names(vec.1) <- nm.1

    ## Translate from GeoHiSSE par names to ClaSSE.
    ## Need to update 'TranslateParsMakerGeoHiSSE' for the simulations.
    parkey <- TranslateParsMakerGeoHiSSE(k=hidden.areas, add.extinction=add.extinction, add.jumps=add.jumps)
    geohisse.pars <- c(model.pars.vec, vec.01, vec.0, vec.1)
    classe.pars <- vector("numeric", length=nrow(parkey))
    names(classe.pars) <- parkey[,1]
    for( i in 1:length(geohisse.pars) ){
        classe.pars[parkey[,2] %in% names(geohisse.pars[i])] <- geohisse.pars[i]
    }

    ## Get a vector of parameters in the correct format for ClaSSE.
    full.classe.pars <- GetArgnamesClasse(k=(hidden.areas+1)*3)
    full.classe.pars <- setNames(rep(0, times=length(full.classe.pars)), full.classe.pars)
    mm <- match(names(classe.pars), table=names(full.classe.pars))
    full.classe.pars[mm] <- classe.pars

    ## Translate x0 to ClaSSE format.
    par.areas <- paste0(c("01", "0", "1"), rep(LETTERS[1:(hidden.areas+1)], each=3))
    if( !x0 %in% par.areas ) stop(paste0("x0 needs to be one of ", paste(par.areas, sep="", collapse=", "), " .", collapse=""))
    classe.x0 <- which( par.areas %in% x0 )

    ## Simulate the phylogenetic tree:
    ## Repeat until we get the sims done. Record how many times it failed.
    attempt <- 0
    sims <- NULL
    print( "Simulating the phylogeny..." )
    while( is.null(sims) ){
    sims <- diversitree::tree.classe(pars=full.classe.pars, max.taxa=max.taxa, max.t=max.time
                                   , include.extinct=include.extinct, x0=classe.x0)
    attempt <- attempt + 1
    if( is.null(sims) ) print( paste("Simulation attemp ", attempt," failed. Trying again...", collapse="") )
    }
    print( "Simulation finished!" )
    
    ## Translate the traits back to GeoHiSSE format.
    tip.state <- sapply(sims$tip.state, function(x) par.areas[x])
    sims$node.state <- sapply(sims$node.state, function(x) par.areas[x])

    ## Return a vector with the true ranges at the tips and a data matrix in the correct format for a call to the GeoHiSSE function.
    geohisse.par.areas <- rep(0:2, times=hidden.areas+1) ## Code in the order for the GeoHiSSE function.
    geohisse.tip.state <- sapply(sims$tip.state, function(x) geohisse.par.areas[x])
    geohisse.mat <- matrix(NA, ncol = 2, nrow = length(geohisse.tip.state) )
    geohisse.mat[,1] <- names( geohisse.tip.state )
    geohisse.mat[,2] <- as.numeric( geohisse.tip.state )
    colnames( geohisse.mat ) <- c("Species", "Range")
    
    return( list(phy=sims, data=geohisse.mat, hidden.areas=tip.state, sim.attempts=attempt, pars=pars, classe.pars=full.classe.pars) )
}

GetArgnamesClasse <- function(k){
    ## Function to generate the correct argument names for 'tree.classe'.
    ## This is a copy from package 'diversitree'. Function: diversitree:::default.argnames.classe
    ## Please find original function here: https://github.com/richfitz/diversitree/blob/master/R/model-classe.R
    fmt <- sprintf("%%0%dd", ceiling(log10(k + .5)))
    sstr <- sprintf(fmt, 1:k)
    lambda.names <- sprintf("lambda%s%s%s", rep(sstr, each=k*(k+1)/2),
                            rep(rep(sstr, times=seq(k,1,-1)), k), 
                            unlist(lapply(1:k, function(i) sstr[i:k])))
    mu.names <- sprintf("mu%s", sstr)
    q.names <- sprintf("q%s%s", rep(sstr, each=k-1), 
                       unlist(lapply(1:k, function(i) sstr[-i])))
    c(lambda.names, mu.names, q.names)
}

ParNamesKeyClaSSEtoGeoHiSSE <- function(k, add.extinction, add.jumps){
    ## Generate a matrix with translation between pars of GeoHiSSE to ClaSSE parameters.
    ## If k >= 3 then all the numbers need to be 0 padded. Need to modify this function.
    ## add.extinction and add.jumps allow for adding additional parameters to the models.
    if( !add.extinction & !add.jumps ){
        mod.par <- sapply(LETTERS[1:(k+1)], function(x) paste0(c("s0", "s0", "s1", "s1", "s01", "x0", "x1", "d0", "d1", "x0", "x1"), x) )
    }
    if( add.extinction & !add.jumps ){
        ## Here we just need to mark the extirpation and extinction in different ways so the parameters can receive different values.
        mod.par <- sapply(LETTERS[1:(k+1)], function(x) paste0(c("s0", "s0", "s1", "s1", "s01", "x0", "x1", "d0", "d1", "x*0", "x*1"), x) )
    }
    if( !add.extinction & add.jumps ){
        mod.par <- sapply(LETTERS[1:(k+1)], function(x) paste0(c("s0", "s0", "s1", "s1", "s01", "x0", "x1", "d0", "d1", "jd0", "jd1", "x0", "x1"), x) )
    }
    if( add.extinction & add.jumps ){
        mod.par <- sapply(LETTERS[1:(k+1)], function(x) paste0(c("s0", "s0", "s1", "s1", "s01", "x0", "x1", "d0", "d1", "jd0", "jd1", "x*0", "x*1"), x) )
    }
    
    if( k >= 3 ){
        lambda.classe.par <- sapply(0:k, function(x) paste0("lambda", c(paste0(sprintf("%02d", c(1,1,2)+(x*3)), collapse="")
                                                                      , paste0(sprintf("%02d", c(2,2,2)+(x*3)), collapse="")
                                                                      , paste0(sprintf("%02d", c(1,1,3)+(x*3)), collapse="")
                                                                      , paste0(sprintf("%02d", c(3,3,3)+(x*3)), collapse="")
                                                                      , paste0(sprintf("%02d", c(1,2,3)+(x*3)), collapse="")
                                                                        )
                                                            )
                                    )
        q.classe.par <- sapply(0:k, function(x) paste0("q", c(paste0(sprintf("%02d", c(1,3)+(x*3)), collapse="")
                                                            , paste0(sprintf("%02d", c(1,2)+(x*3)), collapse="")
                                                            , paste0(sprintf("%02d", c(2,1)+(x*3)), collapse="")
                                                            , paste0(sprintf("%02d", c(3,1)+(x*3)), collapse="")
                                                              )
                                                       )
                               )
        if( add.jumps ){
            ## This adds the parameters to allow for the jump between the endemic areas.
            q.classe.par.jumps <- sapply(0:k, function(x) paste0("q", c(paste0(sprintf("%02d", c(2,3)+(x*3)), collapse="")
                                                                , paste0(sprintf("%02d", c(3,2)+(x*3)), collapse="")
                                                                  )
                                                           )
                                   )
        }
        if( !add.jumps ){
            q.classe.par.jumps <- NULL
        }
        
        mu.classe.par <- sapply(0:k, function(x) paste0("mu", c(paste0(sprintf("%02d", c(2)+(x*3)), collapse="")
                                                              , paste0(sprintf("%02d", c(3)+(x*3)), collapse="")
                                                                )
                                                        )
                                )
    } else{
        lambda.classe.par <- sapply(0:k, function(x) paste0("lambda", c(112,222,113,333,123)+(111*(x*3))))
        q.classe.par <- sapply(0:k, function(x) paste0("q", c(13,12,21,31)+(11*(x*3))))
        if( add.jumps ){
            q.classe.par.jumps <- sapply(0:k, function(x) paste0("q", c(23,32)+(11*(x*3))))
        }
        if( !add.jumps ){
            q.classe.par.jumps <- NULL
        }
        mu.classe.par <- sapply(0:k, function(x) paste0("mu", c(2,3)+(x*3)))
    }
    classe.par <- rbind(lambda.classe.par, q.classe.par, q.classe.par.jumps, mu.classe.par)
    mod.par.vec <- as.vector(mod.par)
    classe.par.vec <- as.vector(classe.par)
    res <- cbind(classe.par.vec, mod.par.vec)
    colnames(res) <- c("ClaSSE", "GeoHiSSE")
    return(res)
}

qNamesKeyClaSSEtoGeoSSE <- function(k, area, init){
    ## k: The number of hidden states in the GeoHiSSE model.
    ## area: One of the three areas 01, 0 or 1.
    ## init: The initial number for the areas. This is the order that 01, 0 and 1 appears in the model. 01 = 1, 0 = 2, 1 = 3.
    ## Messing with this might change the order of the states in the simulation and you might be simulating stuff different from what you think.
    if( k >= 3 ){
        classe.id <- sprintf("%02d", seq(from=init, by=3, length.out=k+1))
        ## sprintf is formatting numbers to have fixed two digits.
    } else{
        classe.id <- seq(from=init, by=3, length.out=k+1)
    }
    geohisse.id <- LETTERS[1:(k+1)]

    q.hidden.classe <- vector(mode="character", length=((k+1)^2)-(k+1))
    q.hidden.geohisse <- vector(mode="character", length=((k+1)^2)-(k+1))
    ## Length is number of cells in a matrix without the diagonal
    count <- 1
    for( i in 1:(k+1)){
        for(j in 1:(k+1)){
            if( i == j ) next
            q.hidden.classe[count] <- paste0("q", classe.id[i], classe.id[j])
            q.hidden.geohisse[count] <- paste0("q", area, geohisse.id[i], geohisse.id[j])
            count <- count+1
        }
    }
    res <- cbind(q.hidden.classe, q.hidden.geohisse)
    colnames(res) <- c("ClaSSE", "GeoHiSSE")
    return(res)
}

TranslateParsMakerGeoHiSSE <- function(k, add.extinction=FALSE, add.jumps=FALSE){
    ## Returns a table with the names of the parameters for the ClaSSE model that corresponds to the
    ## GeoHiSSE model. All other parameters of the ClaSSE model are set to zero.
    ## k: number of hidden states in the model.
    ## add.extinction: if the extinction parameter should be separated from extirpation.
    ## add.jumps: if jumps between endemic areas should be added.

    divpars <- ParNamesKeyClaSSEtoGeoHiSSE(k, add.extinction=add.extinction, add.jumps=add.jumps)
    transpars.01 <- qNamesKeyClaSSEtoGeoSSE(k=k, area="01", init=1)
    transpars.0 <- qNamesKeyClaSSEtoGeoSSE(k=k, area="0", init=2)
    transpars.1 <- qNamesKeyClaSSEtoGeoSSE(k=k, area="1", init=3)

    parkey <- rbind( divpars, transpars.01, transpars.0, transpars.1 )
    colnames( parkey ) <- c("classe.pars", "geohisse.pars")
    return( parkey )
}
