
######################################################################################################################################
######################################################################################################################################
### TransMatMaker -- Builds transition rate matrix for easy use in the main function
######################################################################################################################################
######################################################################################################################################


TransMatMakerfHiSSE <- function(hidden.traits=0, make.null=FALSE, cat.trans.vary=FALSE){
    
    if(hidden.traits == 0){
        trans.mat <- matrix(0, 2, 2)
        trans.mat[2,1] <- 1
        trans.mat[1,2] <- 2
        diag(trans.mat) <- NA
        rownames(trans.mat) <- colnames(trans.mat) <-  c("(0)","(1)")
    }else{
        trans.mat <- matrix(NA, 2*(hidden.traits+1), 2*(hidden.traits+1))
        Range <- ncol(trans.mat) - 1
        Range <- -Range:Range
        Range <- Range[order(abs(Range))]
        max.par <- 0
        for(diag.index in Range[c(2:(length(Range)-2))]){
            if(make.null == TRUE){
                if(diag.index == -1 | diag.index == 1){
                    max.par <- max.par + 1
                    replace.vec <- rep(NA, 3)
                    tmp.vec <- 1:((hidden.traits+1)*2)
                    replace.vec[tmp.vec[tmp.vec%%2==1]] <- max.par
                    trans.mat[row(trans.mat) == (col(trans.mat) - diag.index)] <- replace.vec
                    max.par.tmp <- max.par + 1
                }else{
                    if(diag.index == -2 | diag.index == 2){
                        if(cat.trans.vary == TRUE){
                            max.par <- max.par + 1
                            if(hidden.traits==1){
                                replace.vec <- c(rep(max.par, 2*hidden.traits))
                            }
                            if(hidden.traits==2){
                                replace.vec <- c(rep(max.par, 2*(hidden.traits-1)), rep(max.par+1, 2*(hidden.traits-1)))
                                max.par <- max.par + 1
                            }
                            if(hidden.traits==3){
                                replace.vec <- c(rep(max.par, 2*(hidden.traits-2)), rep(max.par+1, 2*(hidden.traits-2)), rep(max.par+2, 2*(hidden.traits-2)))
                                max.par <- max.par + 2
                            }
                        }else{
                            max.par <- max.par.tmp
                            replace.vec <- c(rep(max.par, 2*hidden.traits))
                        }
                        trans.mat[row(trans.mat) == (col(trans.mat) - diag.index)] <- replace.vec
                    }
                    if(diag.index == -4 | diag.index == 4){
                        if(cat.trans.vary == TRUE){
                            max.par <- max.par + 1
                            if(hidden.traits == 2){
                                replace.vec <- rep(max.par, 2*(hidden.traits-1))
                            }
                            if(hidden.traits == 3){
                                replace.vec <- c(rep(max.par, 2*(hidden.traits-2)), rep(max.par+1, 2*(hidden.traits-2)))
                                max.par <- max.par + 1
                            }
                        }else{
                            max.par <- max.par.tmp
                            replace.vec <- rep(max.par, 2*(hidden.traits-1))
                        }
                        trans.mat[row(trans.mat) == (col(trans.mat) - diag.index)] <- replace.vec
                    }
                    if(diag.index == -6 | diag.index == 6){
                        if(cat.trans.vary == TRUE){
                            max.par <- max.par + 1
                            replace.vec <- rep(max.par, 2*(hidden.traits-2))
                            
                        }else{
                            max.par <- max.par.tmp
                            replace.vec <- rep(max.par, 2)
                        }
                        trans.mat[row(trans.mat) == (col(trans.mat) - diag.index)] <- replace.vec
                    }
                }
            }else{
                if(diag.index == -1 | diag.index == 1){
                    state.trans <- 1:(2*(hidden.traits+1))
                    replace.vec <- rep(NA, 3)
                    tmp.vec <- 1:((hidden.traits+1)*2)
                    if(diag.index == -1){
                        replace.vec[tmp.vec[tmp.vec%%2==1]] <- state.trans[state.trans%%2 == 1]
                    }else{
                        replace.vec[tmp.vec[tmp.vec%%2==1]] <- state.trans[state.trans%%2 == 0]
                    }
                    trans.mat[row(trans.mat) == (col(trans.mat) - diag.index)] <- replace.vec
                    max.par <- 2*(hidden.traits+1)
                    max.par.tmp <- 2*(hidden.traits+1) + 1
                }else{
                    if(diag.index == -2 | diag.index == 2){
                        if(cat.trans.vary == TRUE){
                            max.par <- max.par + 1
                            if(hidden.traits==1){
                                replace.vec <- c(rep(max.par, 2*hidden.traits))
                            }
                            if(hidden.traits==2){
                                replace.vec <- c(rep(max.par, 2*(hidden.traits-1)), rep(max.par+1, 2*(hidden.traits-1)))
                                max.par <- max.par + 1
                            }
                            if(hidden.traits==3){
                                replace.vec <- c(rep(max.par, 2*(hidden.traits-2)), rep(max.par+1, 2*(hidden.traits-2)), rep(max.par+2, 2*(hidden.traits-2)))
                                max.par <- max.par + 2
                            }
                        }else{
                            max.par <- max.par.tmp
                            replace.vec <- c(rep(max.par, 2*hidden.traits))
                        }
                        trans.mat[row(trans.mat) == (col(trans.mat) - diag.index)] <- replace.vec
                    }
                    if(diag.index == -4 | diag.index == 4){
                        if(cat.trans.vary == TRUE){
                            max.par <- max.par + 1
                            if(hidden.traits == 2){
                                replace.vec <- rep(max.par, 2*(hidden.traits-1))
                            }
                            if(hidden.traits == 3){
                                replace.vec <- c(rep(max.par, 2*(hidden.traits-2)), rep(max.par+1, 2*(hidden.traits-2)))
                                max.par <- max.par + 1
                            }
                        }else{
                            max.par <- max.par.tmp
                            replace.vec <- rep(max.par, 2*(hidden.traits-1))
                        }
                        trans.mat[row(trans.mat) == (col(trans.mat) - diag.index)] <- replace.vec
                    }
                    if(diag.index == -6 | diag.index == 6){
                        if(cat.trans.vary == TRUE){
                            max.par <- max.par + 1
                            replace.vec <- rep(max.par, 2*(hidden.traits-2))
                            
                        }else{
                            max.par <- max.par.tmp
                            replace.vec <- rep(max.par, 2*(hidden.traits-2))
                        }
                        trans.mat[row(trans.mat) == (col(trans.mat) - diag.index)] <- replace.vec
                    }
                }
            }
            rownames(trans.mat) <- colnames(trans.mat) <- c("(0A)","(1A)","(0B)","(1B)","(0C)","(1C)","(0D)","(1D)")[1:((hidden.traits+1)*2)]
        }
    }
    return(trans.mat)
}





GetAllDiags <- function(inmat, sorted = TRUE) {
    Range <- ncol(inmat) - 1
    Range <- -Range:Range
    if (isTRUE(sorted)) Range <- Range[order(abs(Range))]
    lapply(Range, function(x) {
        inmat[row(inmat) == (col(inmat) - x)]
    })
}



