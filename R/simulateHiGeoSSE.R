## Function to simulate a tree using the HiGeoSSE process. Depends on 'tree.classe' function of the package diversitree.

## Little issue with an unexported function from diversitree.
## 'diversitree:::default.argnames.classe'
## As an option I will lift/copy the function as an unexported function of the hisse package.
## So we get 'hisse:::default.argnames.classe'

## Entries will be a matrix with exactly 7 rows in the same order as the arguments of the GeoSSE model. The number of columns will inform the number of hidden states.
## Also a squared matrix with dimensions equal to a multiplier of 3 (ordered A, B, AB). Multiples of 3 will inform the number of hidden states.
## Function need to privide a option to just return a list with empty parameters for the own function in the correct format. This will help a lot when using this.
## Function will return the same vectors that the SimulateHisse does. However, returning the $results in the same format will be hard. Will just return a phylogeny with tip states maybe.

SimulateHiGeoSSE <- function(mod.matrix, trans.matrix, max.taxa=Inf, max.time=Inf, max.wall.time=Inf, include.extinct=FALSE, x0="AB0", override.safeties=FALSE, return.empty.pars=FALSE, control.empty.pars=NULL){
    

}

## An additional function that 'SimulateHiGeoSSE' needs.

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
