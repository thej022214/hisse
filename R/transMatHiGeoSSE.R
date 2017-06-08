

TransMatMakerHiGeoSSE <- function(hidden.states=FALSE, hidden.cats=1, make.null=FALSE){
    if(hidden.states == FALSE){
        trans.mat <- TransMatGeoSSEsingle(cat.number=1)
    }else{
        if(hidden.cats ==1){
            stop("Need to specify more than one hidden category.")
        }
        if(hidden.cats == 2){
            if(make.null == TRUE){
                sub.mat1 <- TransMatGeoSSEsingle(cat.number=1)
                sub.mat2 <- TransMatGeoSSEsingle(cat.number=1)
                sub.mat3 <- matrix(NA, 3, 3)
                diag(sub.mat3) <- 3
                trans.mat <- rbind(cbind(sub.mat1,sub.mat3), cbind(sub.mat4,sub.mat2))
            }else{
                sub.mat1 <- TransMatGeoSSEsingle(cat.number=1)
                sub.mat2 <- TransMatGeoSSEsingle(cat.number=2)
                sub.mat3 <- matrix(NA, 3, 3)
                diag(sub.mat3) <- 5
                trans.mat <- rbind(cbind(sub.mat1, sub.mat3), cbind(sub.mat3, sub.mat2))
            }
            rownames(trans.mat) <- colnames(trans.mat) <- c("(0A)","(1A)","(01A)", "(0B)","(1B)","(01B)")
        }

        if(hidden.cats == 3){
            if(make.null == TRUE){
                sub.mat1 <- TransMatGeoSSEsingle(cat.number=1)
                sub.mat2 <- matrix(NA, 3, 3)
                diag(sub.mat2) <- 3
                trans.mat <- rbind(cbind(sub.mat1, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat1, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat1))
            }else{
                sub.mat1 <- TransMatGeoSSEsingle(cat.number=1)
                sub.mat2 <- TransMatGeoSSEsingle(cat.number=2)
                sub.mat3 <- TransMatGeoSSEsingle(cat.number=3)
                sub.mat4 <- matrix(NA, 3, 3)
                diag(sub.mat4) <- 7
                trans.mat <- rbind(cbind(sub.mat1, sub.mat4, sub.mat4), cbind(sub.mat4, sub.mat2, sub.mat4), cbind(sub.mat4, sub.mat4, sub.mat3))
            }
            rownames(trans.mat) <- colnames(trans.mat) <- c("(0A)","(1A)","(01A)", "(0B)","(1B)","(01B)", "(0C)","(1C)","(01C)")
        }

        if(hidden.cats == 4){
            if(make.null == TRUE){
                sub.mat1 <- TransMatGeoSSEsingle(cat.number=1)
                sub.mat2 <- matrix(NA, 3, 3)
                diag(sub.mat2) <- 3
                trans.mat <- rbind(cbind(sub.mat1, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat1, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat1, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat1))
            }else{
                sub.mat1 <- TransMatGeoSSEsingle(cat.number=1)
                sub.mat2 <- TransMatGeoSSEsingle(cat.number=2)
                sub.mat3 <- TransMatGeoSSEsingle(cat.number=3)
                sub.mat4 <- TransMatGeoSSEsingle(cat.number=4)
                sub.mat5 <- matrix(NA, 3, 3)
                diag(sub.mat5) <- 9
                trans.mat <- rbind(cbind(sub.mat1, sub.mat5, sub.mat5, sub.mat5), cbind(sub.mat5, sub.mat2, sub.mat5, sub.mat5), cbind(sub.mat5, sub.mat5, sub.mat3, sub.mat5), cbind(sub.mat5, sub.mat5, sub.mat5, sub.mat4))
            }
            rownames(trans.mat) <- colnames(trans.mat) <- c("(0A)","(1A)","(01A)", "(0B)","(1B)","(01B)", "(0C)","(1C)","(01C)", "(0D)","(1D)","(01D)")
        }
        
        if(hidden.cats == 5){
            if(make.null == TRUE){
                sub.mat1 <- TransMatGeoSSEsingle(cat.number=1)
                sub.mat2 <- matrix(NA, 3, 3)
                diag(sub.mat2) <- 3
                trans.mat <- rbind(cbind(sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1))
            }else{
                sub.mat1 <- TransMatGeoSSEsingle(cat.number=1)
                sub.mat2 <- TransMatGeoSSEsingle(cat.number=2)
                sub.mat3 <- TransMatGeoSSEsingle(cat.number=3)
                sub.mat4 <- TransMatGeoSSEsingle(cat.number=4)
                sub.mat5 <- TransMatGeoSSEsingle(cat.number=5)
                sub.mat6 <- matrix(NA, 3, 3)
                diag(sub.mat6) <- 11
                trans.mat <- rbind(cbind(sub.mat1, sub.mat6, sub.mat6, sub.mat6, sub.mat6), cbind(sub.mat6, sub.mat2, sub.mat6, sub.mat6, sub.mat6), cbind(sub.mat6, sub.mat6, sub.mat3, sub.mat6, sub.mat6), cbind(sub.mat6, sub.mat6, sub.mat6, sub.mat4, sub.mat6), cbind(sub.mat6, sub.mat6, sub.mat6, sub.mat6, sub.mat5))
            }
            rownames(trans.mat)  <- colnames(trans.mat) <- c("(0A)","(1A)","(01A)", "(0B)","(1B)","(01B)", "(0C)","(1C)","(01C)", "(0D)","(1D)","(01D)", "(0E)","(1E)","(01E)")
        }
    }
    return(trans.mat)
}


TransMatGeoSSEsingle <- function(cat.number=1){
    if(cat.number == 1){
        rate.mat <- matrix(0, 3, 3)
        diag(rate.mat) <- 3
        rate.mat[1,3] <- 1
        rate.mat[2,3] <- 2
    }
    if(cat.number == 2){
        rate.mat <- matrix(0, 3, 3)
        diag(rate.mat) <- 3
        rate.mat[1,3] <- 3
        rate.mat[2,3] <- 4
    }
    if(cat.number == 3){
        rate.mat <- matrix(0, 3, 3)
        diag(rate.mat) <- 3
        rate.mat[1,3] <- 5
        rate.mat[2,3] <- 6
    }
    if(cat.number == 4){
        rate.mat <- matrix(0, 3, 3)
        diag(rate.mat) <- 3
        rate.mat[1,3] <- 7
        rate.mat[2,3] <- 8
    }
    if(cat.number == 5){
        rate.mat <- matrix(0, 3, 3)
        diag(rate.mat) <- 3
        rate.mat[1,3] <- 9
        rate.mat[2,3] <- 10
    }
    diag(rate.mat) <- NA
    rownames(rate.mat) <- colnames(rate.mat) <-  c("(0)","(1)","(01)")
    return(rate.mat)
}



