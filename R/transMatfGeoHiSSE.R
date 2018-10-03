
######################################################################################################################################
######################################################################################################################################
### TransMatMaker -- Builds transition rate matrix for GeoHiSSE and our special MuSSE model
######################################################################################################################################
######################################################################################################################################

TransMatMakerfGeoHiSSE <- function(hidden.areas=0, make.null=FALSE, include.jumps=FALSE, separate.extirpation=FALSE){
    if(hidden.areas == 0){
        trans.mat <- TransMatfGeoSSEsingle(cat.number=1, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
    }else{
        if(hidden.areas == 1){
            if(make.null == TRUE){
                sub.mat1 <- TransMatfGeoSSEsingle(cat.number=1, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat2 <- TransMatfGeoSSEsingle(cat.number=1, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                max.par <- max(sub.mat2, na.rm=TRUE)
                sub.mat3 <- matrix(NA, 3, 3)
                diag(sub.mat3) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1,sub.mat3), cbind(sub.mat3,sub.mat2))
            }else{
                sub.mat1 <- TransMatfGeoSSEsingle(cat.number=1, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat2 <- TransMatfGeoSSEsingle(cat.number=2, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                max.par <- max(sub.mat2, na.rm=TRUE)
                sub.mat3 <- matrix(NA, 3, 3)
                diag(sub.mat3) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat3), cbind(sub.mat3, sub.mat2))
            }
            rownames(trans.mat) <- colnames(trans.mat) <- c("(00A)","(11A)","(01A)", "(00B)","(11B)","(01B)")
        }

        if(hidden.areas == 2){
            if(make.null == TRUE){
                sub.mat1 <- TransMatfGeoSSEsingle(cat.number=1, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                max.par <- max(sub.mat1, na.rm=TRUE)
                sub.mat2 <- matrix(NA, 3, 3)
                diag(sub.mat2) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat1, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat1))
            }else{
                sub.mat1 <- TransMatfGeoSSEsingle(cat.number=1, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat2 <- TransMatfGeoSSEsingle(cat.number=2, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat3 <- TransMatfGeoSSEsingle(cat.number=3, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat4 <- matrix(NA, 3, 3)
                max.par <- max(sub.mat3, na.rm=TRUE)
                diag(sub.mat4) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat4, sub.mat4), cbind(sub.mat4, sub.mat2, sub.mat4), cbind(sub.mat4, sub.mat4, sub.mat3))
            }
            rownames(trans.mat) <- colnames(trans.mat) <- c("(00A)","(11A)","(01A)", "(00B)","(11B)","(01B)", "(00C)","(11C)","(01C)")
        }

        if(hidden.areas == 3){
            if(make.null == TRUE){
                sub.mat1 <- TransMatfGeoSSEsingle(cat.number=1, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                max.par <- max(sub.mat1, na.rm=TRUE)
                sub.mat2 <- matrix(NA, 3, 3)
                diag(sub.mat2) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat1, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat1, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat1))
            }else{
                sub.mat1 <- TransMatfGeoSSEsingle(cat.number=1, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat2 <- TransMatfGeoSSEsingle(cat.number=2, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat3 <- TransMatfGeoSSEsingle(cat.number=3, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat4 <- TransMatfGeoSSEsingle(cat.number=4, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat5 <- matrix(NA, 3, 3)
                max.par <- max(sub.mat4, na.rm=TRUE)
                diag(sub.mat5) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat5, sub.mat5, sub.mat5), cbind(sub.mat5, sub.mat2, sub.mat5, sub.mat5), cbind(sub.mat5, sub.mat5, sub.mat3, sub.mat5), cbind(sub.mat5, sub.mat5, sub.mat5, sub.mat4))
            }
            rownames(trans.mat) <- colnames(trans.mat) <- c("(00A)","(11A)","(01A)", "(00B)","(11B)","(01B)", "(00C)","(11C)","(01C)", "(00D)","(11D)","(01D)")
        }
        
        if(hidden.areas == 4){
            if(make.null == TRUE){
                sub.mat1 <- TransMatfGeoSSEsingle(cat.number=1, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                max.par <- max(sub.mat1, na.rm=TRUE)
                sub.mat2 <- matrix(NA, 3, 3)
                diag(sub.mat2) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1))
            }else{
                sub.mat1 <- TransMatfGeoSSEsingle(cat.number=1, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat2 <- TransMatfGeoSSEsingle(cat.number=2, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat3 <- TransMatfGeoSSEsingle(cat.number=3, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat4 <- TransMatfGeoSSEsingle(cat.number=4, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat5 <- TransMatfGeoSSEsingle(cat.number=5, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat6 <- matrix(NA, 3, 3)
                max.par <- max(sub.mat5, na.rm=TRUE)
                diag(sub.mat6) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat6, sub.mat6, sub.mat6, sub.mat6), cbind(sub.mat6, sub.mat2, sub.mat6, sub.mat6, sub.mat6), cbind(sub.mat6, sub.mat6, sub.mat3, sub.mat6, sub.mat6), cbind(sub.mat6, sub.mat6, sub.mat6, sub.mat4, sub.mat6), cbind(sub.mat6, sub.mat6, sub.mat6, sub.mat6, sub.mat5))
            }
            rownames(trans.mat)  <- colnames(trans.mat) <- c("(00A)","(11A)","(01A)", "(00B)","(11B)","(01B)", "(00C)","(11C)","(01C)", "(00D)","(11D)","(01D)", "(00E)","(11E)","(01E)")
        }
        
        if(hidden.areas == 5){
            if(make.null == TRUE){
                sub.mat1 <- TransMatfGeoSSEsingle(cat.number=1, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                max.par <- max(sub.mat1, na.rm=TRUE)
                sub.mat2 <- matrix(NA, 3, 3)
                diag(sub.mat2) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1))
            }else{
                sub.mat1 <- TransMatfGeoSSEsingle(cat.number=1, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat2 <- TransMatfGeoSSEsingle(cat.number=2, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat3 <- TransMatfGeoSSEsingle(cat.number=3, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat4 <- TransMatfGeoSSEsingle(cat.number=4, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat5 <- TransMatfGeoSSEsingle(cat.number=5, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat6 <- TransMatfGeoSSEsingle(cat.number=6, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat7 <- matrix(NA, 3, 3)
                max.par <- max(sub.mat6, na.rm=TRUE)
                diag(sub.mat7) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat7, sub.mat7, sub.mat7, sub.mat7, sub.mat7), cbind(sub.mat7, sub.mat2, sub.mat7, sub.mat7, sub.mat7, sub.mat7), cbind(sub.mat7, sub.mat7, sub.mat3, sub.mat7, sub.mat7, sub.mat7), cbind(sub.mat7, sub.mat7, sub.mat7, sub.mat4, sub.mat7, sub.mat7), cbind(sub.mat7, sub.mat7, sub.mat7, sub.mat7, sub.mat5, sub.mat7), cbind(sub.mat7, sub.mat7, sub.mat7, sub.mat7, sub.mat7, sub.mat6))
            }
            rownames(trans.mat)  <- colnames(trans.mat) <- c("(00A)","(11A)","(01A)", "(00B)","(11B)","(01B)", "(00C)","(11C)","(01C)", "(00D)","(11D)","(01D)", "(00E)","(11E)","(01E)", "(00F)","(11F)","(01F)")
        }
        
        if(hidden.areas == 6){
            if(make.null == TRUE){
                sub.mat1 <- TransMatfGeoSSEsingle(cat.number=1, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                max.par <- max(sub.mat1, na.rm=TRUE)
                sub.mat2 <- matrix(NA, 3, 3)
                diag(sub.mat2) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1))
            }else{
                sub.mat1 <- TransMatfGeoSSEsingle(cat.number=1, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat2 <- TransMatfGeoSSEsingle(cat.number=2, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat3 <- TransMatfGeoSSEsingle(cat.number=3, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat4 <- TransMatfGeoSSEsingle(cat.number=4, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat5 <- TransMatfGeoSSEsingle(cat.number=5, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat6 <- TransMatfGeoSSEsingle(cat.number=6, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat7 <- TransMatfGeoSSEsingle(cat.number=7, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat8 <- matrix(NA, 3, 3)
                max.par <- max(sub.mat7, na.rm=TRUE)
                diag(sub.mat8) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat8, sub.mat8, sub.mat8, sub.mat8, sub.mat8, sub.mat8), cbind(sub.mat8, sub.mat2, sub.mat8, sub.mat8, sub.mat8, sub.mat8, sub.mat8), cbind(sub.mat8, sub.mat8, sub.mat3, sub.mat8, sub.mat8, sub.mat8, sub.mat8), cbind(sub.mat8, sub.mat8, sub.mat8, sub.mat4, sub.mat8, sub.mat8, sub.mat8), cbind(sub.mat8, sub.mat8, sub.mat8, sub.mat8, sub.mat5, sub.mat8, sub.mat8), cbind(sub.mat8, sub.mat8, sub.mat8, sub.mat8, sub.mat8, sub.mat6, sub.mat8), cbind(sub.mat8, sub.mat8, sub.mat8, sub.mat8, sub.mat8, sub.mat8, sub.mat7))
            }
            rownames(trans.mat)  <- colnames(trans.mat) <- c("(00A)","(11A)","(01A)", "(00B)","(11B)","(01B)", "(00C)","(11C)","(01C)", "(00D)","(11D)","(01D)", "(00E)","(11E)","(01E)", "(00F)","(11F)","(01F)", "(00G)","(11G)","(01G)")
        }
        
        if(hidden.areas == 7){
            if(make.null == TRUE){
                sub.mat1 <- TransMatfGeoSSEsingle(cat.number=1, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                max.par <- max(sub.mat1, na.rm=TRUE)
                sub.mat2 <- matrix(NA, 3, 3)
                diag(sub.mat2) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1))
            }else{
                sub.mat1 <- TransMatfGeoSSEsingle(cat.number=1, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat2 <- TransMatfGeoSSEsingle(cat.number=2, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat3 <- TransMatfGeoSSEsingle(cat.number=3, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat4 <- TransMatfGeoSSEsingle(cat.number=4, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat5 <- TransMatfGeoSSEsingle(cat.number=5, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat6 <- TransMatfGeoSSEsingle(cat.number=6, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat7 <- TransMatfGeoSSEsingle(cat.number=7, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat8 <- TransMatfGeoSSEsingle(cat.number=8, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat9 <- matrix(NA, 3, 3)
                max.par <- max(sub.mat8, na.rm=TRUE)
                diag(sub.mat9) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat9), cbind(sub.mat9, sub.mat2, sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat9), cbind(sub.mat9, sub.mat9, sub.mat3, sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat9), cbind(sub.mat9, sub.mat9, sub.mat9, sub.mat4, sub.mat9, sub.mat9, sub.mat9, sub.mat9), cbind(sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat5, sub.mat9, sub.mat9, sub.mat9), cbind(sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat6, sub.mat9, sub.mat9), cbind(sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat7, sub.mat9), cbind(sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat8))
            }
            rownames(trans.mat)  <- colnames(trans.mat) <- c("(00A)","(11A)","(01A)", "(00B)","(11B)","(01B)", "(00C)","(11C)","(01C)", "(00D)","(11D)","(01D)", "(00E)","(11E)","(01E)", "(00F)","(11F)","(01F)", "(00G)","(11G)","(01G)", "(00H)","(11H)","(01H)")
        }

        if(hidden.areas == 8){
            if(make.null == TRUE){
                sub.mat1 <- TransMatfGeoSSEsingle(cat.number=1, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                max.par <- max(sub.mat1, na.rm=TRUE)
                sub.mat2 <- matrix(NA, 3, 3)
                diag(sub.mat2) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1))
            }else{
                sub.mat1 <- TransMatfGeoSSEsingle(cat.number=1, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat2 <- TransMatfGeoSSEsingle(cat.number=2, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat3 <- TransMatfGeoSSEsingle(cat.number=3, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat4 <- TransMatfGeoSSEsingle(cat.number=4, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat5 <- TransMatfGeoSSEsingle(cat.number=5, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat6 <- TransMatfGeoSSEsingle(cat.number=6, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat7 <- TransMatfGeoSSEsingle(cat.number=7, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat8 <- TransMatfGeoSSEsingle(cat.number=8, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat9 <- TransMatfGeoSSEsingle(cat.number=9, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat10 <- matrix(NA, 3, 3)
                max.par <- max(sub.mat9, na.rm=TRUE)
                diag(sub.mat10) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat10, sub.mat10, sub.mat10, sub.mat10, sub.mat10, sub.mat10, sub.mat10, sub.mat10), cbind(sub.mat10, sub.mat2, sub.mat10, sub.mat10, sub.mat10, sub.mat10, sub.mat10, sub.mat10, sub.mat10), cbind(sub.mat10, sub.mat10, sub.mat3, sub.mat10, sub.mat10, sub.mat10, sub.mat10, sub.mat10, sub.mat10), cbind(sub.mat10, sub.mat10, sub.mat10, sub.mat4, sub.mat10, sub.mat10, sub.mat10, sub.mat10, sub.mat10), cbind(sub.mat10, sub.mat10, sub.mat10, sub.mat10, sub.mat5, sub.mat10, sub.mat10, sub.mat10, sub.mat10), cbind(sub.mat10, sub.mat10, sub.mat10, sub.mat10, sub.mat10, sub.mat6, sub.mat10, sub.mat10, sub.mat10), cbind(sub.mat10, sub.mat10, sub.mat10, sub.mat10, sub.mat10, sub.mat10, sub.mat7, sub.mat10, sub.mat10), cbind(sub.mat10, sub.mat10, sub.mat10, sub.mat10, sub.mat10, sub.mat10, sub.mat10, sub.mat8, sub.mat10), cbind(sub.mat10, sub.mat10, sub.mat10, sub.mat10, sub.mat10, sub.mat10, sub.mat10, sub.mat10, sub.mat9))
            }
            rownames(trans.mat)  <- colnames(trans.mat) <- c("(00A)","(11A)","(01A)", "(00B)","(11B)","(01B)", "(00C)","(11C)","(01C)", "(00D)","(11D)","(01D)", "(00E)","(11E)","(01E)", "(00F)","(11F)","(01F)", "(00G)","(11G)","(01G)", "(00H)","(11H)","(01H)", "(00I)","(11I)","(01I)")
        }
        
        if(hidden.areas == 9){
            if(make.null == TRUE){
                sub.mat1 <- TransMatfGeoSSEsingle(cat.number=1, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                max.par <- max(sub.mat1, na.rm=TRUE)
                sub.mat2 <- matrix(NA, 3, 3)
                diag(sub.mat2) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1))
            }else{
                sub.mat1 <- TransMatfGeoSSEsingle(cat.number=1, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat2 <- TransMatfGeoSSEsingle(cat.number=2, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat3 <- TransMatfGeoSSEsingle(cat.number=3, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat4 <- TransMatfGeoSSEsingle(cat.number=4, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat5 <- TransMatfGeoSSEsingle(cat.number=5, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat6 <- TransMatfGeoSSEsingle(cat.number=6, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat7 <- TransMatfGeoSSEsingle(cat.number=7, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat8 <- TransMatfGeoSSEsingle(cat.number=8, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat9 <- TransMatfGeoSSEsingle(cat.number=9, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat10 <- TransMatfGeoSSEsingle(cat.number=10, include.jumps=include.jumps, separate.extirpation=separate.extirpation)
                sub.mat11 <- matrix(NA, 3, 3)
                max.par <- max(sub.mat10, na.rm=TRUE)
                diag(sub.mat11) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat11, sub.mat11, sub.mat11, sub.mat11, sub.mat11, sub.mat11, sub.mat11, sub.mat11, sub.mat11), cbind(sub.mat11, sub.mat2, sub.mat11, sub.mat11, sub.mat11, sub.mat11, sub.mat11, sub.mat11, sub.mat11, sub.mat11), cbind(sub.mat11, sub.mat11, sub.mat3, sub.mat11, sub.mat11, sub.mat11, sub.mat11, sub.mat11, sub.mat11, sub.mat11), cbind(sub.mat11, sub.mat11, sub.mat11, sub.mat4, sub.mat11, sub.mat11, sub.mat11, sub.mat11, sub.mat11, sub.mat11), cbind(sub.mat11, sub.mat11, sub.mat11, sub.mat11, sub.mat5, sub.mat11, sub.mat11, sub.mat11, sub.mat11, sub.mat11), cbind(sub.mat11, sub.mat11, sub.mat11, sub.mat11, sub.mat11, sub.mat6, sub.mat11, sub.mat11, sub.mat11, sub.mat11), cbind(sub.mat11, sub.mat11, sub.mat11, sub.mat11, sub.mat11, sub.mat11, sub.mat7, sub.mat11, sub.mat11, sub.mat11), cbind(sub.mat11, sub.mat11, sub.mat11, sub.mat11, sub.mat11, sub.mat11, sub.mat11, sub.mat8, sub.mat11, sub.mat11), cbind(sub.mat11, sub.mat11, sub.mat11, sub.mat11, sub.mat11, sub.mat11, sub.mat11, sub.mat11, sub.mat9, sub.mat11), cbind(sub.mat11, sub.mat11, sub.mat11, sub.mat11, sub.mat11, sub.mat11, sub.mat11, sub.mat11, sub.mat11, sub.mat10))
            }
            rownames(trans.mat)  <- colnames(trans.mat) <- c("(00A)","(11A)","(01A)", "(00B)","(11B)","(01B)", "(00C)","(11C)","(01C)", "(00D)","(11D)","(01D)", "(00E)","(11E)","(01E)", "(00F)","(11F)","(01F)", "(00G)","(11G)","(01G)", "(00H)","(11H)","(01H)", "(00I)","(11I)","(01I)", "(00J)","(11J)","(01J)")
        }
    }
    return(trans.mat)
}


######################################################################################################################################
######################################################################################################################################
### Support function for generating matrices within matrices
######################################################################################################################################
######################################################################################################################################

TransMatfGeoSSEsingle <- function(cat.number=1, include.jumps=FALSE, separate.extirpation=FALSE){
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
    }else{
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
    rownames(rate.mat) <- colnames(rate.mat) <-  c("(00)","(11)","(01)")
    return(rate.mat)
}
