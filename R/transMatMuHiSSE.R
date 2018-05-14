
######################################################################################################################################
######################################################################################################################################
### TransMatMaker -- Builds transition rate matrix for GeoHiSSE and our special MuSSE model
######################################################################################################################################
######################################################################################################################################

TransMatMakerMuHiSSE <- function(hidden.traits=0, make.null=FALSE, include.diagonals=FALSE, make.special.for.teo=FALSE){
    if(hidden.traits == 0){
        trans.mat <- TransMatMuSSEsingle(cat.number=1, include.diagonals=include.diagonals)
    }else{
        if(hidden.traits == 1){
            if(make.null == TRUE){
                sub.mat1 <- TransMatMuSSEsingle(cat.number=1, include.diagonals=include.diagonals)
                sub.mat2 <- TransMatMuSSEsingle(cat.number=1, include.diagonals=include.diagonals)
                max.par <- max(sub.mat2, na.rm=TRUE)
                sub.mat3 <- matrix(NA, 4, 4)
                diag(sub.mat3) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1,sub.mat3), cbind(sub.mat3,sub.mat2))
            }else{
                sub.mat1 <- TransMatMuSSEsingle(cat.number=1, include.diagonals=include.diagonals)
                sub.mat2 <- TransMatMuSSEsingle(cat.number=2, include.diagonals=include.diagonals)
                max.par <- max(sub.mat2, na.rm=TRUE)
                sub.mat3 <- matrix(NA, 4, 4)
                diag(sub.mat3) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat3), cbind(sub.mat3, sub.mat2))
            }
            rownames(trans.mat) <- colnames(trans.mat) <- c("(00A)","(01A)","(10A)","(11A)", "(00B)","(01B)","(10B)","(11B)")
        }
        
        if(hidden.traits == 2){
            if(make.null == TRUE){
                sub.mat1 <- TransMatMuSSEsingle(cat.number=1, include.diagonals=include.diagonals)
                max.par <- max(sub.mat1, na.rm=TRUE)
                sub.mat2 <- matrix(NA, 4, 4)
                diag(sub.mat2) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat1, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat1))
            }else{
                sub.mat1 <- TransMatMuSSEsingle(cat.number=1, include.diagonals=include.diagonals)
                sub.mat2 <- TransMatMuSSEsingle(cat.number=2, include.diagonals=include.diagonals)
                sub.mat3 <- TransMatMuSSEsingle(cat.number=3, include.diagonals=include.diagonals)
                sub.mat4 <- matrix(NA, 4, 4)
                max.par <- max(sub.mat3, na.rm=TRUE)
                diag(sub.mat4) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat4, sub.mat4), cbind(sub.mat4, sub.mat2, sub.mat4), cbind(sub.mat4, sub.mat4, sub.mat3))
            }
            rownames(trans.mat) <- colnames(trans.mat) <- c("(00A)","(01A)","(10A)","(11A)", "(00B)","(01B)","(10B)","(11B)", "(00C)","(01C)","(10C)","(11C)")
        }
        
        if(hidden.traits == 3){
            if(make.special.for.teo == TRUE){
                sub.mat1 <- TransMatMuSSEsingle(cat.number=1, include.diagonals=include.diagonals)
                sub.mat2 <- TransMatMuSSEsingle(cat.number=2, include.diagonals=include.diagonals)
                max.par <- max(sub.mat2, na.rm=TRUE)
                sub.mat3 <- matrix(NA, 4, 4)
                diag(sub.mat3) <- max.par + 1
                max.par <- max(sub.mat3, na.rm=TRUE)
                sub.mat4 <- matrix(NA, 4, 4)
                diag(sub.mat4) <- max.par + 1
                trans.mat.1 <- rbind(cbind(sub.mat1, sub.mat4), cbind(sub.mat3, sub.mat2))
                max.par <- max(sub.mat4, na.rm=TRUE)
                sub.mat5 <- matrix(NA, 8, 8)
                diag(sub.mat5) <- max.par + 1
                trans.mat <- rbind(cbind(trans.mat.1,sub.mat5), cbind(sub.mat5, trans.mat.1))
            }else{
                if(make.null == TRUE){
                    sub.mat1 <- TransMatMuSSEsingle(cat.number=1, include.diagonals=include.diagonals)
                    max.par <- max(sub.mat1, na.rm=TRUE)
                    sub.mat2 <- matrix(NA, 4, 4)
                    diag(sub.mat2) <- max.par + 1
                    trans.mat <- rbind(cbind(sub.mat1, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat1, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat1, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat1))
                }else{
                    sub.mat1 <- TransMatMuSSEsingle(cat.number=1, include.diagonals=include.diagonals)
                    sub.mat2 <- TransMatMuSSEsingle(cat.number=2, include.diagonals=include.diagonals)
                    sub.mat3 <- TransMatMuSSEsingle(cat.number=3, include.diagonals=include.diagonals)
                    sub.mat4 <- TransMatMuSSEsingle(cat.number=4, include.diagonals=include.diagonals)
                    sub.mat5 <- matrix(NA, 4, 4)
                    max.par <- max(sub.mat4, na.rm=TRUE)
                    diag(sub.mat5) <- max.par + 1
                    trans.mat <- rbind(cbind(sub.mat1, sub.mat5, sub.mat5, sub.mat5), cbind(sub.mat5, sub.mat2, sub.mat5, sub.mat5), cbind(sub.mat5, sub.mat5, sub.mat3, sub.mat5), cbind(sub.mat5, sub.mat5, sub.mat5, sub.mat4))
                }
            }
            rownames(trans.mat) <- colnames(trans.mat) <- c("(00A)","(01A)","(10A)","(11A)", "(00B)","(01B)","(10B)","(11B)", "(00C)","(01C)","(10C)","(11C)", "(00D)","(01D)","(10D)","(11D)")
        }
        
        if(hidden.traits == 4){
            if(make.null == TRUE){
                sub.mat1 <- TransMatMuSSEsingle(cat.number=1, include.diagonals=include.diagonals)
                max.par <- max(sub.mat1, na.rm=TRUE)
                sub.mat2 <- matrix(NA, 4, 4)
                diag(sub.mat2) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1))
            }else{
                sub.mat1 <- TransMatMuSSEsingle(cat.number=1, include.diagonals=include.diagonals)
                sub.mat2 <- TransMatMuSSEsingle(cat.number=2, include.diagonals=include.diagonals)
                sub.mat3 <- TransMatMuSSEsingle(cat.number=3, include.diagonals=include.diagonals)
                sub.mat4 <- TransMatMuSSEsingle(cat.number=4, include.diagonals=include.diagonals)
                sub.mat5 <- TransMatMuSSEsingle(cat.number=5, include.diagonals=include.diagonals)
                sub.mat6 <- matrix(NA, 4, 4)
                max.par <- max(sub.mat5, na.rm=TRUE)
                diag(sub.mat6) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat6, sub.mat6, sub.mat6, sub.mat6), cbind(sub.mat6, sub.mat2, sub.mat6, sub.mat6, sub.mat6), cbind(sub.mat6, sub.mat6, sub.mat3, sub.mat6, sub.mat6), cbind(sub.mat6, sub.mat6, sub.mat6, sub.mat4, sub.mat6),  cbind(sub.mat6, sub.mat6, sub.mat6, sub.mat6, sub.mat5))
            }
            rownames(trans.mat) <- colnames(trans.mat) <- c("(00A)","(01A)","(10A)","(11A)", "(00B)","(01B)","(10B)","(11B)", "(00C)","(01C)","(10C)","(11C)", "(00D)","(01D)","(10D)","(11D)", "(00E)","(01E)","(10E)","(11E)")
        }
        
        if(hidden.traits == 5){
            if(make.null == TRUE){
                sub.mat1 <- TransMatMuSSEsingle(cat.number=1, include.diagonals=include.diagonals)
                max.par <- max(sub.mat1, na.rm=TRUE)
                sub.mat2 <- matrix(NA, 4, 4)
                diag(sub.mat2) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1))
            }else{
                sub.mat1 <- TransMatMuSSEsingle(cat.number=1, include.diagonals=include.diagonals)
                sub.mat2 <- TransMatMuSSEsingle(cat.number=2, include.diagonals=include.diagonals)
                sub.mat3 <- TransMatMuSSEsingle(cat.number=3, include.diagonals=include.diagonals)
                sub.mat4 <- TransMatMuSSEsingle(cat.number=4, include.diagonals=include.diagonals)
                sub.mat5 <- TransMatMuSSEsingle(cat.number=5, include.diagonals=include.diagonals)
                sub.mat6 <- TransMatMuSSEsingle(cat.number=6, include.diagonals=include.diagonals)
                sub.mat7 <- matrix(NA, 4, 4)
                max.par <- max(sub.mat6, na.rm=TRUE)
                diag(sub.mat7) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat7, sub.mat7, sub.mat7, sub.mat7, sub.mat7), cbind(sub.mat7, sub.mat2, sub.mat7, sub.mat7, sub.mat7, sub.mat7), cbind(sub.mat7, sub.mat7, sub.mat3, sub.mat7, sub.mat7, sub.mat7), cbind(sub.mat7, sub.mat7, sub.mat7, sub.mat4, sub.mat7, sub.mat7),  cbind(sub.mat7, sub.mat7, sub.mat7, sub.mat7, sub.mat5, sub.mat7), cbind(sub.mat7, sub.mat7, sub.mat7, sub.mat7, sub.mat7, sub.mat6))
            }
            rownames(trans.mat) <- colnames(trans.mat) <- c("(00A)","(01A)","(10A)","(11A)", "(00B)","(01B)","(10B)","(11B)", "(00C)","(01C)","(10C)","(11C)", "(00D)","(01D)","(10D)","(11D)", "(00E)","(01E)","(10E)","(11E)", "(00F)","(01F)","(10F)","(11F)")
        }
        
        if(hidden.traits == 6){
            if(make.null == TRUE){
                sub.mat1 <- TransMatMuSSEsingle(cat.number=1, include.diagonals=include.diagonals)
                max.par <- max(sub.mat1, na.rm=TRUE)
                sub.mat2 <- matrix(NA, 4, 4)
                diag(sub.mat2) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1))
            }else{
                sub.mat1 <- TransMatMuSSEsingle(cat.number=1, include.diagonals=include.diagonals)
                sub.mat2 <- TransMatMuSSEsingle(cat.number=2, include.diagonals=include.diagonals)
                sub.mat3 <- TransMatMuSSEsingle(cat.number=3, include.diagonals=include.diagonals)
                sub.mat4 <- TransMatMuSSEsingle(cat.number=4, include.diagonals=include.diagonals)
                sub.mat5 <- TransMatMuSSEsingle(cat.number=5, include.diagonals=include.diagonals)
                sub.mat6 <- TransMatMuSSEsingle(cat.number=6, include.diagonals=include.diagonals)
                sub.mat7 <- TransMatMuSSEsingle(cat.number=7, include.diagonals=include.diagonals)
                sub.mat8 <- matrix(NA, 4, 4)
                max.par <- max(sub.mat7, na.rm=TRUE)
                diag(sub.mat8) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat8, sub.mat8, sub.mat8, sub.mat8, sub.mat8, sub.mat8), cbind(sub.mat8, sub.mat2, sub.mat8, sub.mat8, sub.mat8, sub.mat8, sub.mat8), cbind(sub.mat8, sub.mat8, sub.mat3, sub.mat8, sub.mat8, sub.mat8, sub.mat8), cbind(sub.mat8, sub.mat8, sub.mat8, sub.mat4, sub.mat8, sub.mat8, sub.mat8),  cbind(sub.mat8, sub.mat8, sub.mat8, sub.mat8, sub.mat5, sub.mat8, sub.mat8), cbind(sub.mat8, sub.mat8, sub.mat8, sub.mat8, sub.mat8, sub.mat6, sub.mat8), cbind(sub.mat8, sub.mat8, sub.mat8, sub.mat8, sub.mat8, sub.mat8, sub.mat7))
            }
            rownames(trans.mat) <- colnames(trans.mat) <- c("(00A)","(01A)","(10A)","(11A)", "(00B)","(01B)","(10B)","(11B)", "(00C)","(01C)","(10C)","(11C)", "(00D)","(01D)","(10D)","(11D)", "(00E)","(01E)","(10E)","(11E)", "(00F)","(01F)","(10F)","(11F)", "(00G)","(01G)","(10G)","(11G)")
        }
        
        if(hidden.traits == 7){
            if(make.null == TRUE){
                sub.mat1 <- TransMatMuSSEsingle(cat.number=1, include.diagonals=include.diagonals)
                max.par <- max(sub.mat1, na.rm=TRUE)
                sub.mat2 <- matrix(NA, 4, 4)
                diag(sub.mat2) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1))
            }else{
                sub.mat1 <- TransMatMuSSEsingle(cat.number=1, include.diagonals=include.diagonals)
                sub.mat2 <- TransMatMuSSEsingle(cat.number=2, include.diagonals=include.diagonals)
                sub.mat3 <- TransMatMuSSEsingle(cat.number=3, include.diagonals=include.diagonals)
                sub.mat4 <- TransMatMuSSEsingle(cat.number=4, include.diagonals=include.diagonals)
                sub.mat5 <- TransMatMuSSEsingle(cat.number=5, include.diagonals=include.diagonals)
                sub.mat6 <- TransMatMuSSEsingle(cat.number=6, include.diagonals=include.diagonals)
                sub.mat7 <- TransMatMuSSEsingle(cat.number=7, include.diagonals=include.diagonals)
                sub.mat8 <- TransMatMuSSEsingle(cat.number=8, include.diagonals=include.diagonals)
                sub.mat9 <- matrix(NA, 4, 4)
                max.par <- max(sub.mat8, na.rm=TRUE)
                diag(sub.mat9) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat9), cbind(sub.mat9, sub.mat2, sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat9), cbind(sub.mat9, sub.mat9, sub.mat3, sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat9), cbind(sub.mat9, sub.mat9, sub.mat9, sub.mat4, sub.mat9, sub.mat9, sub.mat9, sub.mat9),  cbind(sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat5, sub.mat9, sub.mat9, sub.mat9), cbind(sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat6, sub.mat9, sub.mat9), cbind(sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat7, sub.mat9), cbind(sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat8))
            }
            rownames(trans.mat) <- colnames(trans.mat) <- c("(00A)","(01A)","(10A)","(11A)", "(00B)","(01B)","(10B)","(11B)", "(00C)","(01C)","(10C)","(11C)", "(00D)","(01D)","(10D)","(11D)", "(00E)","(01E)","(10E)","(11E)", "(00F)","(01F)","(10F)","(11F)", "(00G)","(01G)","(10G)","(11G)", "(00H)","(01H)","(10H)","(11H)")
        }
    }
    return(trans.mat)
}


######################################################################################################################################
######################################################################################################################################
### Support function for generating matrices within matrices
######################################################################################################################################
######################################################################################################################################

TransMatMuSSEsingle <- function(cat.number=1, include.diagonals=FALSE){
    if(cat.number == 1){
        rate.mat <- matrix(0, 4, 4)
        diag(rate.mat) <- NA
        if(include.diagonals == TRUE){
            rate.mat[!is.na(rate.mat)] = 1
            matFINAL <- rate.mat
            diag(matFINAL) <- 0
            np <- 12
            index <- matFINAL==1
            rate.mat[index] <- 1:np
        }else{
            #Hard-coded for now
            mat1<-matrix(,4,4)
            mat2<-matrix(,4,4)
            vec.tmp1<-c(0,0,1,1)
            vec.tmp2<-c(0,1,0,1)
            for(i in 1:4){
                mat1[i,]<-abs(vec.tmp1-vec.tmp1[i])
                mat2[i,]<-abs(vec.tmp2-vec.tmp2[i])
            }
            matFINAL <- mat1+mat2
            np <- 8
            index <- matFINAL==1
            rate.mat[index] <- 1:np
        }
    }else{
        rate.mat <- matrix(0, 4, 4)
        diag(rate.mat) <- NA
        if(include.diagonals == TRUE){
            rate.mat[!is.na(rate.mat)] = 1
            matFINAL <- rate.mat
            diag(matFINAL) <- 0
            np <- 12
            index <- matFINAL==1
            rate.mat[index] <- (1 + (12 * (cat.number-1))):(np + (12 * (cat.number-1)))
        }else{
            #Hard-coded for now
            mat1<-matrix(,4,4)
            mat2<-matrix(,4,4)
            vec.tmp1<-c(0,0,1,1)
            vec.tmp2<-c(0,1,0,1)
            for(i in 1:4){
                mat1[i,]<-abs(vec.tmp1-vec.tmp1[i])
                mat2[i,]<-abs(vec.tmp2-vec.tmp2[i])
            }
            matFINAL <- mat1+mat2
            np <- 8
            index <- matFINAL==1
            rate.mat[index] <- (1 + (8 * (cat.number-1))):(np + (8 * (cat.number-1)))
        }
    }
    diag(rate.mat) <- NA
    rownames(rate.mat) <- colnames(rate.mat) <-  c("(00)","(01)","(10)", "(11)")
    return(rate.mat)
}



