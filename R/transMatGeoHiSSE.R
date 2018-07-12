
######################################################################################################################################
######################################################################################################################################
### TransMatMaker -- Builds transition rate matrix for GeoHiSSE and our special MuSSE model
######################################################################################################################################
######################################################################################################################################

TransMatMakerGeoHiSSE <- function(hidden.areas=0, make.null=FALSE, include.jumps=FALSE, separate.extirpation=FALSE){
    if(hidden.areas == 0){
        trans.mat <- TransMatGeoSSEsingle(cat.number=1, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
    }else{
        if(hidden.areas == 1){
            if(make.null == TRUE){
                sub.mat1 <- TransMatGeoSSEsingle(cat.number=1, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat2 <- TransMatGeoSSEsingle(cat.number=1, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                max.par <- max(sub.mat2, na.rm=TRUE)
                sub.mat3 <- matrix(NA, 3, 3)
                diag(sub.mat3) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1,sub.mat3), cbind(sub.mat3,sub.mat2))
            }else{
                sub.mat1 <- TransMatGeoSSEsingle(cat.number=1, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat2 <- TransMatGeoSSEsingle(cat.number=2, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                max.par <- max(sub.mat2, na.rm=TRUE)
                sub.mat3 <- matrix(NA, 3, 3)
                diag(sub.mat3) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat3), cbind(sub.mat3, sub.mat2))
            }
            rownames(trans.mat) <- colnames(trans.mat) <- c("(0A)","(1A)","(01A)", "(0B)","(1B)","(01B)")
        }

        if(hidden.areas == 2){
            if(make.null == TRUE){
                sub.mat1 <- TransMatGeoSSEsingle(cat.number=1, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                max.par <- max(sub.mat1, na.rm=TRUE)
                sub.mat2 <- matrix(NA, 3, 3)
                diag(sub.mat2) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat1, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat1))
            }else{
                sub.mat1 <- TransMatGeoSSEsingle(cat.number=1, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat2 <- TransMatGeoSSEsingle(cat.number=2, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat3 <- TransMatGeoSSEsingle(cat.number=3, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat4 <- matrix(NA, 3, 3)
                max.par <- max(sub.mat3, na.rm=TRUE)
                diag(sub.mat4) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat4, sub.mat4), cbind(sub.mat4, sub.mat2, sub.mat4), cbind(sub.mat4, sub.mat4, sub.mat3))
            }
            rownames(trans.mat) <- colnames(trans.mat) <- c("(0A)","(1A)","(01A)", "(0B)","(1B)","(01B)", "(0C)","(1C)","(01C)")
        }

        if(hidden.areas == 3){
            if(make.null == TRUE){
                sub.mat1 <- TransMatGeoSSEsingle(cat.number=1, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                max.par <- max(sub.mat1, na.rm=TRUE)
                sub.mat2 <- matrix(NA, 3, 3)
                diag(sub.mat2) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat1, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat1, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat1))
            }else{
                sub.mat1 <- TransMatGeoSSEsingle(cat.number=1, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat2 <- TransMatGeoSSEsingle(cat.number=2, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat3 <- TransMatGeoSSEsingle(cat.number=3, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat4 <- TransMatGeoSSEsingle(cat.number=4, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat5 <- matrix(NA, 3, 3)
                max.par <- max(sub.mat4, na.rm=TRUE)
                diag(sub.mat5) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat5, sub.mat5, sub.mat5), cbind(sub.mat5, sub.mat2, sub.mat5, sub.mat5), cbind(sub.mat5, sub.mat5, sub.mat3, sub.mat5), cbind(sub.mat5, sub.mat5, sub.mat5, sub.mat4))
            }
            rownames(trans.mat) <- colnames(trans.mat) <- c("(0A)","(1A)","(01A)", "(0B)","(1B)","(01B)", "(0C)","(1C)","(01C)", "(0D)","(1D)","(01D)")
        }
        
        if(hidden.areas == 4){
            if(make.null == TRUE){
                sub.mat1 <- TransMatGeoSSEsingle(cat.number=1, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                max.par <- max(sub.mat1, na.rm=TRUE)
                sub.mat2 <- matrix(NA, 3, 3)
                diag(sub.mat2) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1))
            }else{
                sub.mat1 <- TransMatGeoSSEsingle(cat.number=1, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat2 <- TransMatGeoSSEsingle(cat.number=2, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat3 <- TransMatGeoSSEsingle(cat.number=3, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat4 <- TransMatGeoSSEsingle(cat.number=4, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat5 <- TransMatGeoSSEsingle(cat.number=5, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat6 <- matrix(NA, 3, 3)
                max.par <- max(sub.mat5, na.rm=TRUE)
                diag(sub.mat6) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat6, sub.mat6, sub.mat6, sub.mat6), cbind(sub.mat6, sub.mat2, sub.mat6, sub.mat6, sub.mat6), cbind(sub.mat6, sub.mat6, sub.mat3, sub.mat6, sub.mat6), cbind(sub.mat6, sub.mat6, sub.mat6, sub.mat4, sub.mat6), cbind(sub.mat6, sub.mat6, sub.mat6, sub.mat6, sub.mat5))
            }
            rownames(trans.mat)  <- colnames(trans.mat) <- c("(0A)","(1A)","(01A)", "(0B)","(1B)","(01B)", "(0C)","(1C)","(01C)", "(0D)","(1D)","(01D)", "(0E)","(1E)","(01E)")
        }
    }
    return(trans.mat)
}


######################################################################################################################################
######################################################################################################################################
### Support function for generating matrices within matrices
######################################################################################################################################
######################################################################################################################################

TransMatGeoSSEsingle <- function(cat.number=1, include.jumps=FALSE, separate.extirpation=FALSE){
    if(cat.number == 1){
        rate.mat <- matrix(0, 3, 3)
        diag(rate.mat) <- NA
        if(include.jumps == TRUE){
            if(separate.extirpation == TRUE){
                #jumps
                rate.mat[2,1] <- 1
                rate.mat[1,2] <- 3
                #extirpation
                rate.mat[3,1] <- 2
                rate.mat[3,2] <- 4
                #normal dispersal
                rate.mat[1,3] <- 5
                rate.mat[2,3] <- 6
            }else{
                #jumps
                rate.mat[2,1] <- 1
                rate.mat[1,2] <- 2
                #normal dispersal
                rate.mat[1,3] <- 3
                rate.mat[2,3] <- 4
            }
        }else{
            if(separate.extirpation == TRUE){
                #extirpation
                rate.mat[3,1] <- 1
                rate.mat[3,2] <- 2
                #normal dispersal
                rate.mat[1,3] <- 3
                rate.mat[2,3] <- 4
            }else{
                #normal dispersal
                rate.mat[1,3] <- 1
                rate.mat[2,3] <- 2
            }
        }
    }
    if(cat.number == 2){
        rate.mat <- matrix(0, 3, 3)
        diag(rate.mat) <- NA
        if(include.jumps == TRUE){
            if(separate.extirpation == TRUE){
                #jumps
                rate.mat[2,1] <- 1 + (6 * (cat.number-1))
                rate.mat[1,2] <- 3 + (6 * (cat.number-1))
                #extirpation
                rate.mat[3,1] <- 2 + (6 * (cat.number-1))
                rate.mat[3,2] <- 4 + (6 * (cat.number-1))
                #normal dispersal
                rate.mat[1,3] <- 5 + (6 * (cat.number-1))
                rate.mat[2,3] <- 6 + (6 * (cat.number-1))
            }else{
                #jumps
                rate.mat[2,1] <- 1 + (4 * (cat.number-1))
                rate.mat[1,2] <- 2 + (4 * (cat.number-1))
                #normal dispersal
                rate.mat[1,3] <- 3 + (4 * (cat.number-1))
                rate.mat[2,3] <- 4 + (4 * (cat.number-1))
            }
        }else{
            if(separate.extirpation == TRUE){
                #extirpation
                rate.mat[3,1] <- 1 + (4 * (cat.number-1))
                rate.mat[3,2] <- 2 + (4 * (cat.number-1))
                #normal dispersal
                rate.mat[1,3] <- 3 + (4 * (cat.number-1))
                rate.mat[2,3] <- 4 + (4 * (cat.number-1))
            }else{
                #normal dispersal
                rate.mat[1,3] <- 1 + (2 * (cat.number-1))
                rate.mat[2,3] <- 2 + (2 * (cat.number-1))
            }
        }
    }
    if(cat.number == 3){
        rate.mat <- matrix(0, 3, 3)
        diag(rate.mat) <- NA
        if(include.jumps == TRUE){
            if(separate.extirpation == TRUE){
                #jumps
                rate.mat[2,1] <- 1 + (6 * (cat.number-1))
                rate.mat[1,2] <- 3 + (6 * (cat.number-1))
                #extirpation
                rate.mat[3,1] <- 2 + (6 * (cat.number-1))
                rate.mat[3,2] <- 4 + (6 * (cat.number-1))
                #normal dispersal
                rate.mat[1,3] <- 5 + (6 * (cat.number-1))
                rate.mat[2,3] <- 6 + (6 * (cat.number-1))
            }else{
                #jumps
                rate.mat[2,1] <- 1 + (4 * (cat.number-1))
                rate.mat[1,2] <- 2 + (4 * (cat.number-1))
                #normal dispersal
                rate.mat[1,3] <- 3 + (4 * (cat.number-1))
                rate.mat[2,3] <- 4 + (4 * (cat.number-1))
            }
        }else{
            if(separate.extirpation == TRUE){
                #extirpation
                rate.mat[3,1] <- 1 + (4 * (cat.number-1))
                rate.mat[3,2] <- 2 + (4 * (cat.number-1))
                #normal dispersal
                rate.mat[1,3] <- 3 + (4 * (cat.number-1))
                rate.mat[2,3] <- 4 + (4 * (cat.number-1))
            }else{
                #normal dispersal
                rate.mat[1,3] <- 1 + (2 * (cat.number-1))
                rate.mat[2,3] <- 2 + (2 * (cat.number-1))
            }
        }
    }
    if(cat.number == 4){
        rate.mat <- matrix(0, 3, 3)
        diag(rate.mat) <- NA
        if(include.jumps == TRUE){
            if(separate.extirpation == TRUE){
                #jumps
                rate.mat[2,1] <- 1 + (6 * (cat.number-1))
                rate.mat[1,2] <- 3 + (6 * (cat.number-1))
                #extirpation
                rate.mat[3,1] <- 2 + (6 * (cat.number-1))
                rate.mat[3,2] <- 4 + (6 * (cat.number-1))
                #normal dispersal
                rate.mat[1,3] <- 5 + (6 * (cat.number-1))
                rate.mat[2,3] <- 6 + (6 * (cat.number-1))
            }else{
                #jumps
                rate.mat[2,1] <- 1 + (4 * (cat.number-1))
                rate.mat[1,2] <- 2 + (4 * (cat.number-1))
                #normal dispersal
                rate.mat[1,3] <- 3 + (4 * (cat.number-1))
                rate.mat[2,3] <- 4 + (4 * (cat.number-1))
            }
        }else{
            if(separate.extirpation == TRUE){
                #extirpation
                rate.mat[3,1] <- 1 + (4 * (cat.number-1))
                rate.mat[3,2] <- 2 + (4 * (cat.number-1))
                #normal dispersal
                rate.mat[1,3] <- 3 + (4 * (cat.number-1))
                rate.mat[2,3] <- 4 + (4 * (cat.number-1))
            }else{
                #normal dispersal
                rate.mat[1,3] <- 1 + (2 * (cat.number-1))
                rate.mat[2,3] <- 2 + (2 * (cat.number-1))
            }
        }
    }
    if(cat.number == 5){
        rate.mat <- matrix(0, 3, 3)
        diag(rate.mat) <- NA
        if(include.jumps == TRUE){
            if(separate.extirpation == TRUE){
                #jumps
                rate.mat[2,1] <- 1 + (6 * (cat.number-1))
                rate.mat[1,2] <- 3 + (6 * (cat.number-1))
                #extirpation
                rate.mat[3,1] <- 2 + (6 * (cat.number-1))
                rate.mat[3,2] <- 4 + (6 * (cat.number-1))
                #normal dispersal
                rate.mat[1,3] <- 5 + (6 * (cat.number-1))
                rate.mat[2,3] <- 6 + (6 * (cat.number-1))
            }else{
                #jumps
                rate.mat[2,1] <- 1 + (4 * (cat.number-1))
                rate.mat[1,2] <- 2 + (4 * (cat.number-1))
                #normal dispersal
                rate.mat[1,3] <- 3 + (4 * (cat.number-1))
                rate.mat[2,3] <- 4 + (4 * (cat.number-1))
            }
        }else{
            if(separate.extirpation == TRUE){
                #extirpation
                rate.mat[3,1] <- 1 + (4 * (cat.number-1))
                rate.mat[3,2] <- 2 + (4 * (cat.number-1))
                #normal dispersal
                rate.mat[1,3] <- 3 + (4 * (cat.number-1))
                rate.mat[2,3] <- 4 + (4 * (cat.number-1))
            }else{
                #normal dispersal
                rate.mat[1,3] <- 1 + (2 * (cat.number-1))
                rate.mat[2,3] <- 2 + (2 * (cat.number-1))
            }
        }
    }
    diag(rate.mat) <- NA
    rownames(rate.mat) <- colnames(rate.mat) <-  c("(0)","(1)","(01)")
    return(rate.mat)
}
