## Functions associated with simulation of a Hidden areas GeoSSE model using the ClaSSE model. GeoSSE is a subset of a ClaSSE model.

SimulateHiGeoSSE <- function(pars, hidden.areas=1, x0="AB0", max.taxa=Inf, max.time=Inf, include.extinct=FALSE, return.HiGeoSSE_pars=FALSE, override.safeties=FALSE){
    if(return.HiGeoSSE_pars){
        mod.matrix <- matrix(data=0, nrow=7, ncol=hidden.areas+1)
        rownames( mod.matrix ) <- c("sAB", "sA", "sB", "xA", "xB", "dA", "dB")
        colnames( mod.matrix ) <- as.character( 0:hidden.areas )
        qmatAB <- matrix(data=0, ncol=hidden.areas+1, nrow=hidden.areas+1)
        diag(qmatAB) <- NA
        trans.names <- paste0("AB", 0:hidden.areas)
        rownames(qmatAB) <- colnames(qmatAB) <- trans.names
        qmatA <- matrix(data=0, ncol=hidden.areas+1, nrow=hidden.areas+1)
        diag(qmatA) <- NA
        trans.names <- paste0("A", 0:hidden.areas)
        rownames(qmatA) <- colnames(qmatA) <- trans.names
        qmatB <- matrix(data=0, ncol=hidden.areas+1, nrow=hidden.areas+1)
        diag(qmatB) <- NA
        trans.names <- paste0("B", 0:hidden.areas)
        rownames(qmatB) <- colnames(qmatB) <- trans.names
        par.list <- list( model.pars = mod.matrix, q.AB = qmatAB, q.A = qmatA, q.B = qmatB)
        class(par.list) <- append(class(par.list),"HiGeoSSE_pars")
        return( par.list )
    }

    ## Make a series of checks:
    if( is.infinite(max.taxa) & is.infinite(max.time)){
        if( override.safeties ) print("WARNING: Simulation running without stopping criteria defined. \n")
        if( !override.safeties ) stop("No stopping criteria have been informed.")
    }
    ## The parameter list is too specific. Better to use the custom class.
    if( !inherits(pars, "HiGeoSSE_pars") ) stop("Argument 'pars' needs to be of class 'HiGeoSSE_pars'. Please check argument 'return.HiGeoSSE_pars'.")

    ## Generate the key for the translation of the parameters:
    model.pars.vec <- c( pars$model.pars )
    names(model.pars.vec) <- paste0(rownames( pars$model.pars ), rep(colnames(pars$model.pars), each=7))

    ## Get the transtions for the hidden areas.
    ## Have a feeling this could be more elegant. But, well...
    tr.size <- ((hidden.areas+1)^2)-(hidden.areas+1)
    AB.vec <- A.vec <- B.vec <- vector("numeric", length=tr.size)
    AB.nm <- A.nm <- B.nm <- vector("character", length=tr.size)
    count <- 1 ## keep the count for the loop.
    for( i in 1:(hidden.areas+1) ){
        for( j in 1:(hidden.areas+1) ){
            if( i == j ) next
            AB.vec[count] <- pars$q.AB[i,j]
            AB.nm[count] <- paste0("qAB", i-1, j-1)
            A.vec[count] <- pars$q.A[i,j]
            A.nm[count] <- paste0("qA", i-1, j-1)
            B.vec[count] <- pars$q.B[i,j]
            B.nm[count] <- paste0("qB", i-1, j-1)
            count <- count + 1
        }
    }
    names(AB.vec) <- AB.nm
    names(A.vec) <- A.nm
    names(B.vec) <- B.nm

    ## Translate from HiGeoSSE par names to ClaSSE.
    parkey <- TranslateParsMakerHiGeoSSE(k=hidden.areas)
    higeosse.pars <- c(model.pars.vec, AB.vec, A.vec, B.vec)
    classe.pars <- vector("numeric", length=nrow(parkey))
    names(classe.pars) <- parkey[,1]
    for( i in 1:length(higeosse.pars) ){
        classe.pars[parkey[,2] %in% names(higeosse.pars[i])] <- higeosse.pars[i]
    }

    ## Get a vector of parameters in the correct format for ClaSSE.
    full.classe.pars <- GetArgnamesClasse(k=(hidden.areas+1)*3)
    full.classe.pars <- setNames(rep(0, times=length(full.classe.pars)), full.classe.pars)
    mm <- match(names(classe.pars), table=names(full.classe.pars))
    full.classe.pars[mm] <- classe.pars

    ## Translate x0 to ClaSSE format.
    par.areas <- paste0(c("AB", "A", "B"), rep(0:hidden.areas, each=3))
    if( !x0 %in% par.areas ) stop(paste0("x0 needs to be one of ", paste(par.areas, sep="", collapse=", "), " .", collapse=""))
    classe.x0 <- which( par.areas %in% x0 )

    ## Simulate the phylogenetic tree:
    sims <- diversitree::tree.classe(pars=full.classe.pars, max.taxa=max.taxa, max.t=max.t
                                   , include.extinct=include.extinct, x0=classe.x0)
    
    ## Need to elaborate the returning object. Now returns only the same output as 'tree.classe' function.
    return( sims )
    
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

ParNamesKeyClaSSEtoHiGeoSSE <- function(k){
    ## Generate a matrix with translation between pars of HiGeoSSE to ClaSSE parameters.
    ## If k >= 3 then all the numbers need to be 0 padded. Need to modify this function.
    mod.par <- sapply(0:k, function(x) paste0(c("sA", "sA", "sB", "sB", "sAB", "xA", "xB", "dA", "dB", "xA", "xB"), x) )
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
        mu.classe.par <- sapply(0:k, function(x) paste0("mu", c(paste0(sprintf("%02d", c(2)+(x*3)), collapse="")
                                                              , paste0(sprintf("%02d", c(3)+(x*3)), collapse="")
                                                                )
                                                        )
                                )
    } else{
        lambda.classe.par <- sapply(0:k, function(x) paste0("lambda", c(112,222,113,333,123)+(111*(x*3))))
        q.classe.par <- sapply(0:k, function(x) paste0("q", c(13,12,21,31)+(11*(x*3))))
        mu.classe.par <- sapply(0:k, function(x) paste0("mu", c(2,3)+(x*3)))
    }
    classe.par <- rbind(lambda.classe.par, q.classe.par, mu.classe.par)
    mod.par.vec <- as.vector(mod.par)
    classe.par.vec <- as.vector(classe.par)
    res <- cbind(classe.par.vec, mod.par.vec)
    colnames(res) <- c("ClaSSE", "HiGeoSSE")
    return(res)
}

qNamesKeyClaSSEtoGeoSSE <- function(k, area, init){
    ## k: The number of hidden states in the HiGeoSSE model.
    ## area: One of the three areas AB, A or B.
    ## init: The initial number for the areas. This is the order that AB, A and B appears in the model. AB = 1, A = 2, B = 3.
    ## Messing with this might change the order of the states in the simulation and you might be simulating stuff different from what you think.
    if( k >= 3 ){
        classe.id <- sprintf("%02d", seq(from=init, by=3, length.out=k+1))
        ## sprintf is formatting numbers to have fixed two digits.
    } else{
        classe.id <- seq(from=init, by=3, length.out=k+1)
    }
    higeosse.id <- seq(from=0, by=1, length.out=k+1)

    q.hidden.classe <- vector(mode="character", length=((k+1)^2)-(k+1))
    q.hidden.higeosse <- vector(mode="character", length=((k+1)^2)-(k+1))
    ## Length is number of cells in a matrix without the diagonal
    count <- 1
    for( i in 1:(k+1)){
        for(j in 1:(k+1)){
            if( i == j ) next
            q.hidden.classe[count] <- paste0("q", classe.id[i], classe.id[j])
            q.hidden.higeosse[count] <- paste0("q", area, higeosse.id[i], higeosse.id[j])
            count <- count+1
        }
    }
    return( cbind(q.hidden.classe, q.hidden.higeosse) )
}

TranslateParsMakerHiGeoSSE <- function(k){
    ## Returns a table with the names of the parameters for the ClaSSE model that corresponds to the
    ## HiGeoSSE model. All other parameters of the ClaSSE model are set to zero.
    ## k: number of hidden states in the model.

    divpars <- ParNamesKeyClaSSEtoHiGeoSSE(k)
    transpars.AB <- qNamesKeyClaSSEtoGeoSSE(k=k, area="AB", init=1)
    transpars.A <- qNamesKeyClaSSEtoGeoSSE(k=k, area="A", init=2)
    transpars.B <- qNamesKeyClaSSEtoGeoSSE(k=k, area="B", init=3)

    parkey <- rbind( divpars, transpars.AB, transpars.A, transpars.B )
    colnames( parkey ) <- c("classe.pars", "higeosse.pars")
    return( parkey )
}
