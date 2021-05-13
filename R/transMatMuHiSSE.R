
######################################################################################################################################
######################################################################################################################################
### TransMatMaker -- Builds transition rate matrix for GeoHiSSE and our special MuSSE model
######################################################################################################################################
######################################################################################################################################

TransMatMakerMuHiSSE <- function(hidden.traits=0, make.null=FALSE, include.diagonals=FALSE, cat.trans.vary=FALSE){
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
                max.par <- max(sub.mat3, na.rm=TRUE)
                sub.mat4 <- matrix(NA, 4, 4)
                diag(sub.mat4) <- max.par + 1
                if(cat.trans.vary == TRUE){
                    trans.mat <- rbind(cbind(sub.mat1,sub.mat4), cbind(sub.mat3,sub.mat2))
                }else{
                    trans.mat <- rbind(cbind(sub.mat1,sub.mat3), cbind(sub.mat3,sub.mat2))
                }
            }else{
                sub.mat1 <- TransMatMuSSEsingle(cat.number=1, include.diagonals=include.diagonals)
                sub.mat2 <- TransMatMuSSEsingle(cat.number=2, include.diagonals=include.diagonals)
                max.par <- max(sub.mat2, na.rm=TRUE)
                sub.mat3 <- matrix(NA, 4, 4)
                diag(sub.mat3) <- max.par + 1
                max.par <- max(sub.mat3, na.rm=TRUE)
                sub.mat4 <- matrix(NA, 4, 4)
                diag(sub.mat4) <- max.par + 1
                if(cat.trans.vary == TRUE){
                    trans.mat <- rbind(cbind(sub.mat1, sub.mat4), cbind(sub.mat3, sub.mat2))
                }else{
                    trans.mat <- rbind(cbind(sub.mat1, sub.mat3), cbind(sub.mat3, sub.mat2))
                }
            }
            rownames(trans.mat) <- colnames(trans.mat) <- c("(00A)","(01A)","(10A)","(11A)", "(00B)","(01B)","(10B)","(11B)")
        }
        
        if(hidden.traits == 2){
            if(make.null == TRUE){
                sub.mat1 <- TransMatMuSSEsingle(cat.number=1, include.diagonals=include.diagonals)
                max.par <- max(sub.mat1, na.rm=TRUE)
                sub.mat2 <- matrix(NA, 4, 4)
                diag(sub.mat2) <- max.par + 1
                max.par <- max(sub.mat2, na.rm=TRUE)
                sub.mat3 <- matrix(NA, 4, 4)
                diag(sub.mat3) <- max.par + 1
                max.par <- max(sub.mat3, na.rm=TRUE)
                sub.mat4 <- matrix(NA, 4, 4)
                diag(sub.mat4) <- max.par + 1
                max.par <- max(sub.mat4, na.rm=TRUE)
                sub.mat5 <- matrix(NA, 4, 4)
                diag(sub.mat5) <- max.par + 1
                max.par <- max(sub.mat5, na.rm=TRUE)
                sub.mat6 <- matrix(NA, 4, 4)
                diag(sub.mat6) <- max.par + 1
                max.par <- max(sub.mat6, na.rm=TRUE)
                sub.mat7 <- matrix(NA, 4, 4)
                diag(sub.mat7) <- max.par + 1
                if(cat.trans.vary == TRUE){
                    trans.mat <- rbind(cbind(sub.mat1, sub.mat4, sub.mat7), cbind(sub.mat2, sub.mat1, sub.mat5), cbind(sub.mat6, sub.mat3, sub.mat1))
                }else{
                    trans.mat <- rbind(cbind(sub.mat1, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat1, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat1))
                }
            }else{
                sub.mat1 <- TransMatMuSSEsingle(cat.number=1, include.diagonals=include.diagonals)
                sub.mat2 <- TransMatMuSSEsingle(cat.number=2, include.diagonals=include.diagonals)
                sub.mat3 <- TransMatMuSSEsingle(cat.number=3, include.diagonals=include.diagonals)
                sub.mat4 <- matrix(NA, 4, 4)
                max.par <- max(sub.mat3, na.rm=TRUE)
                diag(sub.mat4) <- max.par + 1
                sub.mat5 <- matrix(NA, 4, 4)
                max.par <- max(sub.mat4, na.rm=TRUE)
                diag(sub.mat5) <- max.par + 1
                max.par <- max(sub.mat5, na.rm=TRUE)
                sub.mat6 <- matrix(NA, 4, 4)
                diag(sub.mat6) <- max.par + 1
                max.par <- max(sub.mat6, na.rm=TRUE)
                sub.mat7 <- matrix(NA, 4, 4)
                diag(sub.mat7) <- max.par + 1
                max.par <- max(sub.mat7, na.rm=TRUE)
                sub.mat8 <- matrix(NA, 4, 4)
                diag(sub.mat8) <- max.par + 1
                max.par <- max(sub.mat8, na.rm=TRUE)
                sub.mat9 <- matrix(NA, 4, 4)
                diag(sub.mat9) <- max.par + 1
                if(cat.trans.vary == TRUE){
                    trans.mat <- rbind(cbind(sub.mat1, sub.mat6, sub.mat9), cbind(sub.mat4, sub.mat2, sub.mat7), cbind(sub.mat8, sub.mat5, sub.mat3))
                }else{
                    trans.mat <- rbind(cbind(sub.mat1, sub.mat4, sub.mat4), cbind(sub.mat4, sub.mat2, sub.mat4), cbind(sub.mat4, sub.mat4, sub.mat3))
                }
            }
            rownames(trans.mat) <- colnames(trans.mat) <- c("(00A)","(01A)","(10A)","(11A)", "(00B)","(01B)","(10B)","(11B)", "(00C)","(01C)","(10C)","(11C)")
        }
        
        if(hidden.traits == 3){
            if(make.null == TRUE){
                sub.mat1 <- TransMatMuSSEsingle(cat.number=1, include.diagonals=include.diagonals)
                max.par <- max(sub.mat1, na.rm=TRUE)
                sub.mat2 <- matrix(NA, 4, 4)
                diag(sub.mat2) <- max.par + 1
                max.par <- max(sub.mat2, na.rm=TRUE)
                sub.mat3 <- matrix(NA, 4, 4)
                diag(sub.mat3) <- max.par + 1
                max.par <- max(sub.mat3, na.rm=TRUE)
                sub.mat4 <- matrix(NA, 4, 4)
                diag(sub.mat4) <- max.par + 1
                max.par <- max(sub.mat4, na.rm=TRUE)
                sub.mat5 <- matrix(NA, 4, 4)
                diag(sub.mat5) <- max.par + 1
                max.par <- max(sub.mat5, na.rm=TRUE)
                sub.mat6 <- matrix(NA, 4, 4)
                diag(sub.mat6) <- max.par + 1
                max.par <- max(sub.mat6, na.rm=TRUE)
                sub.mat7 <- matrix(NA, 4, 4)
                diag(sub.mat7) <- max.par + 1
                max.par <- max(sub.mat7, na.rm=TRUE)
                sub.mat8 <- matrix(NA, 4, 4)
                diag(sub.mat8) <- max.par + 1
                max.par <- max(sub.mat8, na.rm=TRUE)
                sub.mat9 <- matrix(NA, 4, 4)
                diag(sub.mat9) <- max.par + 1
                max.par <- max(sub.mat9, na.rm=TRUE)
                sub.mat10 <- matrix(NA, 4, 4)
                diag(sub.mat10) <- max.par + 1
                max.par <- max(sub.mat10, na.rm=TRUE)
                sub.mat11 <- matrix(NA, 4, 4)
                diag(sub.mat11) <- max.par + 1
                max.par <- max(sub.mat11, na.rm=TRUE)
                sub.mat12 <- matrix(NA, 4, 4)
                diag(sub.mat12) <- max.par + 1
                max.par <- max(sub.mat12, na.rm=TRUE)
                sub.mat13 <- matrix(NA, 4, 4)
                diag(sub.mat13) <- max.par + 1
                if(cat.trans.vary == TRUE){
                    trans.mat <- rbind(cbind(sub.mat1, sub.mat5, sub.mat10, sub.mat13), cbind(sub.mat2, sub.mat1, sub.mat6, sub.mat11), cbind(sub.mat8, sub.mat3, sub.mat1, sub.mat7), cbind(sub.mat12, sub.mat9, sub.mat4, sub.mat1))
                }else{
                    trans.mat <- rbind(cbind(sub.mat1, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat1, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat1, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat1))
                }
            }else{
                sub.mat1 <- TransMatMuSSEsingle(cat.number=1, include.diagonals=include.diagonals)
                sub.mat2 <- TransMatMuSSEsingle(cat.number=2, include.diagonals=include.diagonals)
                sub.mat3 <- TransMatMuSSEsingle(cat.number=3, include.diagonals=include.diagonals)
                sub.mat4 <- TransMatMuSSEsingle(cat.number=4, include.diagonals=include.diagonals)
                sub.mat5 <- matrix(NA, 4, 4)
                max.par <- max(sub.mat4, na.rm=TRUE)
                diag(sub.mat5) <- max.par + 1
                max.par <- max(sub.mat5, na.rm=TRUE)
                sub.mat6 <- matrix(NA, 4, 4)
                diag(sub.mat6) <- max.par + 1
                max.par <- max(sub.mat6, na.rm=TRUE)
                sub.mat7 <- matrix(NA, 4, 4)
                diag(sub.mat7) <- max.par + 1
                max.par <- max(sub.mat7, na.rm=TRUE)
                sub.mat8 <- matrix(NA, 4, 4)
                diag(sub.mat8) <- max.par + 1
                max.par <- max(sub.mat8, na.rm=TRUE)
                sub.mat9 <- matrix(NA, 4, 4)
                diag(sub.mat9) <- max.par + 1
                max.par <- max(sub.mat9, na.rm=TRUE)
                sub.mat10 <- matrix(NA, 4, 4)
                diag(sub.mat10) <- max.par + 1
                max.par <- max(sub.mat10, na.rm=TRUE)
                sub.mat11 <- matrix(NA, 4, 4)
                diag(sub.mat11) <- max.par + 1
                max.par <- max(sub.mat11, na.rm=TRUE)
                sub.mat12 <- matrix(NA, 4, 4)
                diag(sub.mat12) <- max.par + 1
                max.par <- max(sub.mat12, na.rm=TRUE)
                sub.mat13 <- matrix(NA, 4, 4)
                diag(sub.mat13) <- max.par + 1
                max.par <- max(sub.mat13, na.rm=TRUE)
                sub.mat14 <- matrix(NA, 4, 4)
                diag(sub.mat14) <- max.par + 1
                max.par <- max(sub.mat14, na.rm=TRUE)
                sub.mat15 <- matrix(NA, 4, 4)
                diag(sub.mat15) <- max.par + 1
                max.par <- max(sub.mat15, na.rm=TRUE)
                sub.mat16 <- matrix(NA, 4, 4)
                diag(sub.mat16) <- max.par + 1
                if(cat.trans.vary == TRUE){
                    trans.mat <- rbind(cbind(sub.mat1, sub.mat8, sub.mat13, sub.mat16), cbind(sub.mat5, sub.mat2, sub.mat9, sub.mat14), cbind(sub.mat11, sub.mat6, sub.mat3, sub.mat10), cbind(sub.mat15, sub.mat12, sub.mat7, sub.mat4))
                }else{
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
                max.par <- max(sub.mat2, na.rm=TRUE)
                sub.mat3 <- matrix(NA, 4, 4)
                diag(sub.mat3) <- max.par + 1
                max.par <- max(sub.mat3, na.rm=TRUE)
                sub.mat4 <- matrix(NA, 4, 4)
                diag(sub.mat4) <- max.par + 1
                max.par <- max(sub.mat4, na.rm=TRUE)
                sub.mat5 <- matrix(NA, 4, 4)
                diag(sub.mat5) <- max.par + 1
                max.par <- max(sub.mat5, na.rm=TRUE)
                sub.mat6 <- matrix(NA, 4, 4)
                diag(sub.mat6) <- max.par + 1
                max.par <- max(sub.mat6, na.rm=TRUE)
                sub.mat7 <- matrix(NA, 4, 4)
                diag(sub.mat7) <- max.par + 1
                max.par <- max(sub.mat7, na.rm=TRUE)
                sub.mat8 <- matrix(NA, 4, 4)
                diag(sub.mat8) <- max.par + 1
                max.par <- max(sub.mat8, na.rm=TRUE)
                sub.mat9 <- matrix(NA, 4, 4)
                diag(sub.mat9) <- max.par + 1
                max.par <- max(sub.mat9, na.rm=TRUE)
                sub.mat10 <- matrix(NA, 4, 4)
                diag(sub.mat10) <- max.par + 1
                max.par <- max(sub.mat10, na.rm=TRUE)
                sub.mat11 <- matrix(NA, 4, 4)
                diag(sub.mat11) <- max.par + 1
                max.par <- max(sub.mat11, na.rm=TRUE)
                sub.mat12 <- matrix(NA, 4, 4)
                diag(sub.mat12) <- max.par + 1
                max.par <- max(sub.mat12, na.rm=TRUE)
                sub.mat13 <- matrix(NA, 4, 4)
                diag(sub.mat13) <- max.par + 1
                max.par <- max(sub.mat13, na.rm=TRUE)
                sub.mat14 <- matrix(NA, 4, 4)
                diag(sub.mat14) <- max.par + 1
                max.par <- max(sub.mat14, na.rm=TRUE)
                sub.mat15 <- matrix(NA, 4, 4)
                diag(sub.mat15) <- max.par + 1
                max.par <- max(sub.mat15, na.rm=TRUE)
                sub.mat16 <- matrix(NA, 4, 4)
                diag(sub.mat16) <- max.par + 1
                max.par <- max(sub.mat16, na.rm=TRUE)
                sub.mat17 <- matrix(NA, 4, 4)
                diag(sub.mat17) <- max.par + 1
                max.par <- max(sub.mat17, na.rm=TRUE)
                sub.mat18 <- matrix(NA, 4, 4)
                diag(sub.mat18) <- max.par + 1
                max.par <- max(sub.mat18, na.rm=TRUE)
                sub.mat19 <- matrix(NA, 4, 4)
                diag(sub.mat19) <- max.par + 1
                max.par <- max(sub.mat19, na.rm=TRUE)
                sub.mat20 <- matrix(NA, 4, 4)
                diag(sub.mat20) <- max.par + 1
                max.par <- max(sub.mat20, na.rm=TRUE)
                sub.mat21 <- matrix(NA, 4, 4)
                diag(sub.mat21) <- max.par + 1
                if(cat.trans.vary == TRUE){
                    trans.mat <- rbind(cbind(sub.mat1, sub.mat6, sub.mat13, sub.mat18, sub.mat21), cbind(sub.mat2, sub.mat1, sub.mat7, sub.mat14, sub.mat19), cbind(sub.mat10, sub.mat3, sub.mat1, sub.mat8, sub.mat15), cbind(sub.mat16, sub.mat11, sub.mat4, sub.mat1, sub.mat9), cbind(sub.mat20, sub.mat17, sub.mat12, sub.mat5, sub.mat1))
                }else{
                    trans.mat <- rbind(cbind(sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1))
                }
            }else{
                sub.mat1 <- TransMatMuSSEsingle(cat.number=1, include.diagonals=include.diagonals)
                sub.mat2 <- TransMatMuSSEsingle(cat.number=2, include.diagonals=include.diagonals)
                sub.mat3 <- TransMatMuSSEsingle(cat.number=3, include.diagonals=include.diagonals)
                sub.mat4 <- TransMatMuSSEsingle(cat.number=4, include.diagonals=include.diagonals)
                sub.mat5 <- TransMatMuSSEsingle(cat.number=5, include.diagonals=include.diagonals)
                sub.mat6 <- matrix(NA, 4, 4)
                max.par <- max(sub.mat5, na.rm=TRUE)
                diag(sub.mat6) <- max.par + 1
                max.par <- max(sub.mat6, na.rm=TRUE)
                sub.mat7 <- matrix(NA, 4, 4)
                diag(sub.mat7) <- max.par + 1
                max.par <- max(sub.mat7, na.rm=TRUE)
                sub.mat8 <- matrix(NA, 4, 4)
                diag(sub.mat8) <- max.par + 1
                max.par <- max(sub.mat8, na.rm=TRUE)
                sub.mat9 <- matrix(NA, 4, 4)
                diag(sub.mat9) <- max.par + 1
                max.par <- max(sub.mat9, na.rm=TRUE)
                sub.mat10 <- matrix(NA, 4, 4)
                diag(sub.mat10) <- max.par + 1
                max.par <- max(sub.mat10, na.rm=TRUE)
                sub.mat11 <- matrix(NA, 4, 4)
                diag(sub.mat11) <- max.par + 1
                max.par <- max(sub.mat11, na.rm=TRUE)
                sub.mat12 <- matrix(NA, 4, 4)
                diag(sub.mat12) <- max.par + 1
                max.par <- max(sub.mat12, na.rm=TRUE)
                sub.mat13 <- matrix(NA, 4, 4)
                diag(sub.mat13) <- max.par + 1
                max.par <- max(sub.mat13, na.rm=TRUE)
                sub.mat14 <- matrix(NA, 4, 4)
                diag(sub.mat14) <- max.par + 1
                max.par <- max(sub.mat14, na.rm=TRUE)
                sub.mat15 <- matrix(NA, 4, 4)
                diag(sub.mat15) <- max.par + 1
                max.par <- max(sub.mat15, na.rm=TRUE)
                sub.mat16 <- matrix(NA, 4, 4)
                diag(sub.mat16) <- max.par + 1
                max.par <- max(sub.mat16, na.rm=TRUE)
                sub.mat17 <- matrix(NA, 4, 4)
                diag(sub.mat17) <- max.par + 1
                max.par <- max(sub.mat17, na.rm=TRUE)
                sub.mat18 <- matrix(NA, 4, 4)
                diag(sub.mat18) <- max.par + 1
                max.par <- max(sub.mat18, na.rm=TRUE)
                sub.mat19 <- matrix(NA, 4, 4)
                diag(sub.mat19) <- max.par + 1
                max.par <- max(sub.mat19, na.rm=TRUE)
                sub.mat20 <- matrix(NA, 4, 4)
                diag(sub.mat20) <- max.par + 1
                max.par <- max(sub.mat20, na.rm=TRUE)
                sub.mat21 <- matrix(NA, 4, 4)
                diag(sub.mat21) <- max.par + 1
                max.par <- max(sub.mat21, na.rm=TRUE)
                sub.mat22 <- matrix(NA, 4, 4)
                diag(sub.mat22) <- max.par + 1
                max.par <- max(sub.mat22, na.rm=TRUE)
                sub.mat23 <- matrix(NA, 4, 4)
                diag(sub.mat23) <- max.par + 1
                max.par <- max(sub.mat23, na.rm=TRUE)
                sub.mat24 <- matrix(NA, 4, 4)
                diag(sub.mat24) <- max.par + 1
                max.par <- max(sub.mat24, na.rm=TRUE)
                sub.mat25 <- matrix(NA, 4, 4)
                diag(sub.mat25) <- max.par + 1
                if(cat.trans.vary == TRUE){
                    trans.mat <- rbind(cbind(sub.mat1, sub.mat10, sub.mat17, sub.mat22, sub.mat25), cbind(sub.mat6, sub.mat2, sub.mat11, sub.mat18, sub.mat23), cbind(sub.mat14, sub.mat7, sub.mat3, sub.mat12, sub.mat19), cbind(sub.mat20, sub.mat15, sub.mat8, sub.mat4, sub.mat13), cbind(sub.mat24, sub.mat21, sub.mat16, sub.mat9, sub.mat5))
                }else{
                    trans.mat <- rbind(cbind(sub.mat1, sub.mat6, sub.mat6, sub.mat6, sub.mat6), cbind(sub.mat6, sub.mat2, sub.mat6, sub.mat6, sub.mat6), cbind(sub.mat6, sub.mat6, sub.mat3, sub.mat6, sub.mat6), cbind(sub.mat6, sub.mat6, sub.mat6, sub.mat4, sub.mat6),  cbind(sub.mat6, sub.mat6, sub.mat6, sub.mat6, sub.mat5))
                }
            }
            rownames(trans.mat) <- colnames(trans.mat) <- c("(00A)","(01A)","(10A)","(11A)", "(00B)","(01B)","(10B)","(11B)", "(00C)","(01C)","(10C)","(11C)", "(00D)","(01D)","(10D)","(11D)", "(00E)","(01E)","(10E)","(11E)")
        }
        
        if(hidden.traits == 5){
            if(make.null == TRUE){
                sub.mat1 <- TransMatMuSSEsingle(cat.number=1, include.diagonals=include.diagonals)
                max.par <- max(sub.mat1, na.rm=TRUE)
                sub.mat2 <- matrix(NA, 4, 4)
                diag(sub.mat2) <- max.par + 1
                max.par <- max(sub.mat2, na.rm=TRUE)
                sub.mat3 <- matrix(NA, 4, 4)
                diag(sub.mat3) <- max.par + 1
                max.par <- max(sub.mat3, na.rm=TRUE)
                sub.mat4 <- matrix(NA, 4, 4)
                diag(sub.mat4) <- max.par + 1
                max.par <- max(sub.mat4, na.rm=TRUE)
                sub.mat5 <- matrix(NA, 4, 4)
                diag(sub.mat5) <- max.par + 1
                max.par <- max(sub.mat5, na.rm=TRUE)
                sub.mat6 <- matrix(NA, 4, 4)
                diag(sub.mat6) <- max.par + 1
                max.par <- max(sub.mat6, na.rm=TRUE)
                sub.mat7 <- matrix(NA, 4, 4)
                diag(sub.mat7) <- max.par + 1
                max.par <- max(sub.mat7, na.rm=TRUE)
                sub.mat8 <- matrix(NA, 4, 4)
                diag(sub.mat8) <- max.par + 1
                max.par <- max(sub.mat8, na.rm=TRUE)
                sub.mat9 <- matrix(NA, 4, 4)
                diag(sub.mat9) <- max.par + 1
                max.par <- max(sub.mat9, na.rm=TRUE)
                sub.mat10 <- matrix(NA, 4, 4)
                diag(sub.mat10) <- max.par + 1
                max.par <- max(sub.mat10, na.rm=TRUE)
                sub.mat11 <- matrix(NA, 4, 4)
                diag(sub.mat11) <- max.par + 1
                max.par <- max(sub.mat11, na.rm=TRUE)
                sub.mat12 <- matrix(NA, 4, 4)
                diag(sub.mat12) <- max.par + 1
                max.par <- max(sub.mat12, na.rm=TRUE)
                sub.mat13 <- matrix(NA, 4, 4)
                diag(sub.mat13) <- max.par + 1
                max.par <- max(sub.mat13, na.rm=TRUE)
                sub.mat14 <- matrix(NA, 4, 4)
                diag(sub.mat14) <- max.par + 1
                max.par <- max(sub.mat14, na.rm=TRUE)
                sub.mat15 <- matrix(NA, 4, 4)
                diag(sub.mat15) <- max.par + 1
                max.par <- max(sub.mat15, na.rm=TRUE)
                sub.mat16 <- matrix(NA, 4, 4)
                diag(sub.mat16) <- max.par + 1
                max.par <- max(sub.mat16, na.rm=TRUE)
                sub.mat17 <- matrix(NA, 4, 4)
                diag(sub.mat17) <- max.par + 1
                max.par <- max(sub.mat17, na.rm=TRUE)
                sub.mat18 <- matrix(NA, 4, 4)
                diag(sub.mat18) <- max.par + 1
                max.par <- max(sub.mat18, na.rm=TRUE)
                sub.mat19 <- matrix(NA, 4, 4)
                diag(sub.mat19) <- max.par + 1
                max.par <- max(sub.mat19, na.rm=TRUE)
                sub.mat20 <- matrix(NA, 4, 4)
                diag(sub.mat20) <- max.par + 1
                max.par <- max(sub.mat20, na.rm=TRUE)
                sub.mat21 <- matrix(NA, 4, 4)
                diag(sub.mat21) <- max.par + 1
                max.par <- max(sub.mat21, na.rm=TRUE)
                sub.mat22 <- matrix(NA, 4, 4)
                diag(sub.mat22) <- max.par + 1
                max.par <- max(sub.mat22, na.rm=TRUE)
                sub.mat23 <- matrix(NA, 4, 4)
                diag(sub.mat23) <- max.par + 1
                max.par <- max(sub.mat23, na.rm=TRUE)
                sub.mat24 <- matrix(NA, 4, 4)
                diag(sub.mat24) <- max.par + 1
                max.par <- max(sub.mat24, na.rm=TRUE)
                sub.mat25 <- matrix(NA, 4, 4)
                diag(sub.mat25) <- max.par + 1
                max.par <- max(sub.mat25, na.rm=TRUE)
                sub.mat26 <- matrix(NA, 4, 4)
                diag(sub.mat26) <- max.par + 1
                max.par <- max(sub.mat26, na.rm=TRUE)
                sub.mat27 <- matrix(NA, 4, 4)
                diag(sub.mat27) <- max.par + 1
                max.par <- max(sub.mat27, na.rm=TRUE)
                sub.mat28 <- matrix(NA, 4, 4)
                diag(sub.mat28) <- max.par + 1
                max.par <- max(sub.mat28, na.rm=TRUE)
                sub.mat29 <- matrix(NA, 4, 4)
                diag(sub.mat29) <- max.par + 1
                max.par <- max(sub.mat29, na.rm=TRUE)
                sub.mat30 <- matrix(NA, 4, 4)
                diag(sub.mat30) <- max.par + 1
                max.par <- max(sub.mat30, na.rm=TRUE)
                sub.mat31 <- matrix(NA, 4, 4)
                diag(sub.mat31) <- max.par + 1
                if(cat.trans.vary == TRUE){
                    trans.mat <- rbind(cbind(sub.mat1, sub.mat7, sub.mat16, sub.mat23, sub.mat28, sub.mat31), cbind(sub.mat2, sub.mat1, sub.mat8, sub.mat17, sub.mat24, sub.mat29), cbind(sub.mat12, sub.mat3, sub.mat1, sub.mat9, sub.mat18, sub.mat25), cbind(sub.mat20, sub.mat13, sub.mat4, sub.mat1, sub.mat10, sub.mat19), cbind(sub.mat26, sub.mat21, sub.mat14, sub.mat5, sub.mat1, sub.mat11), cbind(sub.mat30, sub.mat27, sub.mat22, sub.mat15, sub.mat6, sub.mat1))
                }else{
                    trans.mat <- rbind(cbind(sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1))
                }
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
                max.par <- max(sub.mat7, na.rm=TRUE)
                sub.mat8 <- matrix(NA, 4, 4)
                diag(sub.mat8) <- max.par + 1
                max.par <- max(sub.mat8, na.rm=TRUE)
                sub.mat9 <- matrix(NA, 4, 4)
                diag(sub.mat9) <- max.par + 1
                max.par <- max(sub.mat9, na.rm=TRUE)
                sub.mat10 <- matrix(NA, 4, 4)
                diag(sub.mat10) <- max.par + 1
                max.par <- max(sub.mat10, na.rm=TRUE)
                sub.mat11 <- matrix(NA, 4, 4)
                diag(sub.mat11) <- max.par + 1
                max.par <- max(sub.mat11, na.rm=TRUE)
                sub.mat12 <- matrix(NA, 4, 4)
                diag(sub.mat12) <- max.par + 1
                max.par <- max(sub.mat12, na.rm=TRUE)
                sub.mat13 <- matrix(NA, 4, 4)
                diag(sub.mat13) <- max.par + 1
                max.par <- max(sub.mat13, na.rm=TRUE)
                sub.mat14 <- matrix(NA, 4, 4)
                diag(sub.mat14) <- max.par + 1
                max.par <- max(sub.mat14, na.rm=TRUE)
                sub.mat15 <- matrix(NA, 4, 4)
                diag(sub.mat15) <- max.par + 1
                max.par <- max(sub.mat15, na.rm=TRUE)
                sub.mat16 <- matrix(NA, 4, 4)
                diag(sub.mat16) <- max.par + 1
                max.par <- max(sub.mat16, na.rm=TRUE)
                sub.mat17 <- matrix(NA, 4, 4)
                diag(sub.mat17) <- max.par + 1
                max.par <- max(sub.mat17, na.rm=TRUE)
                sub.mat18 <- matrix(NA, 4, 4)
                diag(sub.mat18) <- max.par + 1
                max.par <- max(sub.mat18, na.rm=TRUE)
                sub.mat19 <- matrix(NA, 4, 4)
                diag(sub.mat19) <- max.par + 1
                max.par <- max(sub.mat19, na.rm=TRUE)
                sub.mat20 <- matrix(NA, 4, 4)
                diag(sub.mat20) <- max.par + 1
                max.par <- max(sub.mat20, na.rm=TRUE)
                sub.mat21 <- matrix(NA, 4, 4)
                diag(sub.mat21) <- max.par + 1
                max.par <- max(sub.mat21, na.rm=TRUE)
                sub.mat22 <- matrix(NA, 4, 4)
                diag(sub.mat22) <- max.par + 1
                max.par <- max(sub.mat22, na.rm=TRUE)
                sub.mat23 <- matrix(NA, 4, 4)
                diag(sub.mat23) <- max.par + 1
                max.par <- max(sub.mat23, na.rm=TRUE)
                sub.mat24 <- matrix(NA, 4, 4)
                diag(sub.mat24) <- max.par + 1
                max.par <- max(sub.mat24, na.rm=TRUE)
                sub.mat25 <- matrix(NA, 4, 4)
                diag(sub.mat25) <- max.par + 1
                max.par <- max(sub.mat25, na.rm=TRUE)
                sub.mat26 <- matrix(NA, 4, 4)
                diag(sub.mat26) <- max.par + 1
                max.par <- max(sub.mat26, na.rm=TRUE)
                sub.mat27 <- matrix(NA, 4, 4)
                diag(sub.mat27) <- max.par + 1
                max.par <- max(sub.mat27, na.rm=TRUE)
                sub.mat28 <- matrix(NA, 4, 4)
                diag(sub.mat28) <- max.par + 1
                max.par <- max(sub.mat28, na.rm=TRUE)
                sub.mat29 <- matrix(NA, 4, 4)
                diag(sub.mat29) <- max.par + 1
                max.par <- max(sub.mat29, na.rm=TRUE)
                sub.mat30 <- matrix(NA, 4, 4)
                diag(sub.mat30) <- max.par + 1
                max.par <- max(sub.mat30, na.rm=TRUE)
                sub.mat31 <- matrix(NA, 4, 4)
                diag(sub.mat31) <- max.par + 1
                max.par <- max(sub.mat31, na.rm=TRUE)
                sub.mat32 <- matrix(NA, 4, 4)
                diag(sub.mat32) <- max.par + 1
                max.par <- max(sub.mat32, na.rm=TRUE)
                sub.mat33 <- matrix(NA, 4, 4)
                diag(sub.mat33) <- max.par + 1
                max.par <- max(sub.mat33, na.rm=TRUE)
                sub.mat34 <- matrix(NA, 4, 4)
                diag(sub.mat34) <- max.par + 1
                max.par <- max(sub.mat34, na.rm=TRUE)
                sub.mat35 <- matrix(NA, 4, 4)
                diag(sub.mat35) <- max.par + 1
                max.par <- max(sub.mat35, na.rm=TRUE)
                sub.mat36 <- matrix(NA, 4, 4)
                diag(sub.mat36) <- max.par + 1
                if(cat.trans.vary == TRUE){
                    trans.mat <- rbind(cbind(sub.mat1, sub.mat12, sub.mat21, sub.mat28, sub.mat33, sub.mat36), cbind(sub.mat7, sub.mat2, sub.mat13, sub.mat22, sub.mat29, sub.mat34), cbind(sub.mat17, sub.mat8, sub.mat3, sub.mat14, sub.mat23, sub.mat30), cbind(sub.mat25, sub.mat18, sub.mat9, sub.mat4, sub.mat15, sub.mat24), cbind(sub.mat31, sub.mat26, sub.mat19, sub.mat10, sub.mat5, sub.mat16), cbind(sub.mat35, sub.mat32, sub.mat27, sub.mat20, sub.mat11, sub.mat6))
                }else{
                    trans.mat <- rbind(cbind(sub.mat1, sub.mat7, sub.mat7, sub.mat7, sub.mat7, sub.mat7), cbind(sub.mat7, sub.mat2, sub.mat7, sub.mat7, sub.mat7, sub.mat7), cbind(sub.mat7, sub.mat7, sub.mat3, sub.mat7, sub.mat7, sub.mat7), cbind(sub.mat7, sub.mat7, sub.mat7, sub.mat4, sub.mat7, sub.mat7),  cbind(sub.mat7, sub.mat7, sub.mat7, sub.mat7, sub.mat5, sub.mat7), cbind(sub.mat7, sub.mat7, sub.mat7, sub.mat7, sub.mat7, sub.mat6))
                }
            }
            rownames(trans.mat) <- colnames(trans.mat) <- c("(00A)","(01A)","(10A)","(11A)", "(00B)","(01B)","(10B)","(11B)", "(00C)","(01C)","(10C)","(11C)", "(00D)","(01D)","(10D)","(11D)", "(00E)","(01E)","(10E)","(11E)", "(00F)","(01F)","(10F)","(11F)")
        }
        
        if(hidden.traits == 6){
            if(make.null == TRUE){
                sub.mat1 <- TransMatMuSSEsingle(cat.number=1, include.diagonals=include.diagonals)
                max.par <- max(sub.mat1, na.rm=TRUE)
                sub.mat2 <- matrix(NA, 4, 4)
                diag(sub.mat2) <- max.par + 1
                max.par <- max(sub.mat2, na.rm=TRUE)
                sub.mat3 <- matrix(NA, 4, 4)
                diag(sub.mat3) <- max.par + 1
                max.par <- max(sub.mat3, na.rm=TRUE)
                sub.mat4 <- matrix(NA, 4, 4)
                diag(sub.mat4) <- max.par + 1
                max.par <- max(sub.mat4, na.rm=TRUE)
                sub.mat5 <- matrix(NA, 4, 4)
                diag(sub.mat5) <- max.par + 1
                max.par <- max(sub.mat5, na.rm=TRUE)
                sub.mat6 <- matrix(NA, 4, 4)
                diag(sub.mat6) <- max.par + 1
                max.par <- max(sub.mat6, na.rm=TRUE)
                sub.mat7 <- matrix(NA, 4, 4)
                diag(sub.mat7) <- max.par + 1
                max.par <- max(sub.mat7, na.rm=TRUE)
                sub.mat8 <- matrix(NA, 4, 4)
                diag(sub.mat8) <- max.par + 1
                max.par <- max(sub.mat8, na.rm=TRUE)
                sub.mat9 <- matrix(NA, 4, 4)
                diag(sub.mat9) <- max.par + 1
                max.par <- max(sub.mat9, na.rm=TRUE)
                sub.mat10 <- matrix(NA, 4, 4)
                diag(sub.mat10) <- max.par + 1
                max.par <- max(sub.mat10, na.rm=TRUE)
                sub.mat11 <- matrix(NA, 4, 4)
                diag(sub.mat11) <- max.par + 1
                max.par <- max(sub.mat11, na.rm=TRUE)
                sub.mat12 <- matrix(NA, 4, 4)
                diag(sub.mat12) <- max.par + 1
                max.par <- max(sub.mat12, na.rm=TRUE)
                sub.mat13 <- matrix(NA, 4, 4)
                diag(sub.mat13) <- max.par + 1
                max.par <- max(sub.mat13, na.rm=TRUE)
                sub.mat14 <- matrix(NA, 4, 4)
                diag(sub.mat14) <- max.par + 1
                max.par <- max(sub.mat14, na.rm=TRUE)
                sub.mat15 <- matrix(NA, 4, 4)
                diag(sub.mat15) <- max.par + 1
                max.par <- max(sub.mat15, na.rm=TRUE)
                sub.mat16 <- matrix(NA, 4, 4)
                diag(sub.mat16) <- max.par + 1
                max.par <- max(sub.mat16, na.rm=TRUE)
                sub.mat17 <- matrix(NA, 4, 4)
                diag(sub.mat17) <- max.par + 1
                max.par <- max(sub.mat17, na.rm=TRUE)
                sub.mat18 <- matrix(NA, 4, 4)
                diag(sub.mat18) <- max.par + 1
                max.par <- max(sub.mat18, na.rm=TRUE)
                sub.mat19 <- matrix(NA, 4, 4)
                diag(sub.mat19) <- max.par + 1
                max.par <- max(sub.mat19, na.rm=TRUE)
                sub.mat20 <- matrix(NA, 4, 4)
                diag(sub.mat20) <- max.par + 1
                max.par <- max(sub.mat20, na.rm=TRUE)
                sub.mat21 <- matrix(NA, 4, 4)
                diag(sub.mat21) <- max.par + 1
                max.par <- max(sub.mat21, na.rm=TRUE)
                sub.mat22 <- matrix(NA, 4, 4)
                diag(sub.mat22) <- max.par + 1
                max.par <- max(sub.mat22, na.rm=TRUE)
                sub.mat23 <- matrix(NA, 4, 4)
                diag(sub.mat23) <- max.par + 1
                max.par <- max(sub.mat23, na.rm=TRUE)
                sub.mat24 <- matrix(NA, 4, 4)
                diag(sub.mat24) <- max.par + 1
                max.par <- max(sub.mat24, na.rm=TRUE)
                sub.mat25 <- matrix(NA, 4, 4)
                diag(sub.mat25) <- max.par + 1
                max.par <- max(sub.mat25, na.rm=TRUE)
                sub.mat26 <- matrix(NA, 4, 4)
                diag(sub.mat26) <- max.par + 1
                max.par <- max(sub.mat26, na.rm=TRUE)
                sub.mat27 <- matrix(NA, 4, 4)
                diag(sub.mat27) <- max.par + 1
                max.par <- max(sub.mat27, na.rm=TRUE)
                sub.mat28 <- matrix(NA, 4, 4)
                diag(sub.mat28) <- max.par + 1
                max.par <- max(sub.mat28, na.rm=TRUE)
                sub.mat29 <- matrix(NA, 4, 4)
                diag(sub.mat29) <- max.par + 1
                max.par <- max(sub.mat29, na.rm=TRUE)
                sub.mat30 <- matrix(NA, 4, 4)
                diag(sub.mat30) <- max.par + 1
                max.par <- max(sub.mat30, na.rm=TRUE)
                sub.mat31 <- matrix(NA, 4, 4)
                diag(sub.mat31) <- max.par + 1
                max.par <- max(sub.mat31, na.rm=TRUE)
                sub.mat32 <- matrix(NA, 4, 4)
                diag(sub.mat32) <- max.par + 1
                max.par <- max(sub.mat32, na.rm=TRUE)
                sub.mat33 <- matrix(NA, 4, 4)
                diag(sub.mat33) <- max.par + 1
                max.par <- max(sub.mat33, na.rm=TRUE)
                sub.mat34 <- matrix(NA, 4, 4)
                diag(sub.mat34) <- max.par + 1
                max.par <- max(sub.mat34, na.rm=TRUE)
                sub.mat35 <- matrix(NA, 4, 4)
                diag(sub.mat35) <- max.par + 1
                max.par <- max(sub.mat35, na.rm=TRUE)
                sub.mat36 <- matrix(NA, 4, 4)
                diag(sub.mat36) <- max.par + 1
                max.par <- max(sub.mat36, na.rm=TRUE)
                sub.mat37 <- matrix(NA, 4, 4)
                diag(sub.mat37) <- max.par + 1
                max.par <- max(sub.mat37, na.rm=TRUE)
                sub.mat38 <- matrix(NA, 4, 4)
                diag(sub.mat38) <- max.par + 1
                max.par <- max(sub.mat38, na.rm=TRUE)
                sub.mat39 <- matrix(NA, 4, 4)
                diag(sub.mat39) <- max.par + 1
                max.par <- max(sub.mat39, na.rm=TRUE)
                sub.mat40 <- matrix(NA, 4, 4)
                diag(sub.mat40) <- max.par + 1
                max.par <- max(sub.mat40, na.rm=TRUE)
                sub.mat41 <- matrix(NA, 4, 4)
                diag(sub.mat41) <- max.par + 1
                max.par <- max(sub.mat41, na.rm=TRUE)
                sub.mat42 <- matrix(NA, 4, 4)
                diag(sub.mat42) <- max.par + 1
                max.par <- max(sub.mat42, na.rm=TRUE)
                sub.mat43 <- matrix(NA, 4, 4)
                diag(sub.mat43) <- max.par + 1
                if(cat.trans.vary == TRUE){
                    trans.mat <- rbind(cbind(sub.mat1, sub.mat8, sub.mat19, sub.mat28, sub.mat35, sub.mat40, sub.mat43), cbind(sub.mat2, sub.mat1, sub.mat9, sub.mat20, sub.mat29, sub.mat36, sub.mat41), cbind(sub.mat14, sub.mat3, sub.mat1, sub.mat10, sub.mat21, sub.mat30, sub.mat37), cbind(sub.mat24, sub.mat15, sub.mat4, sub.mat1, sub.mat11, sub.mat22, sub.mat31), cbind(sub.mat32, sub.mat25, sub.mat16, sub.mat5, sub.mat1, sub.mat12, sub.mat23), cbind(sub.mat38, sub.mat33, sub.mat26, sub.mat17, sub.mat6, sub.mat1, sub.mat13), cbind(sub.mat42, sub.mat39, sub.mat34, sub.mat27, sub.mat18, sub.mat7, sub.mat1))
                }else{
                    trans.mat <- rbind(cbind(sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1))
                }
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
                max.par <- max(sub.mat8, na.rm=TRUE)
                sub.mat9 <- matrix(NA, 4, 4)
                diag(sub.mat9) <- max.par + 1
                max.par <- max(sub.mat9, na.rm=TRUE)
                sub.mat10 <- matrix(NA, 4, 4)
                diag(sub.mat10) <- max.par + 1
                max.par <- max(sub.mat10, na.rm=TRUE)
                sub.mat11 <- matrix(NA, 4, 4)
                diag(sub.mat11) <- max.par + 1
                max.par <- max(sub.mat11, na.rm=TRUE)
                sub.mat12 <- matrix(NA, 4, 4)
                diag(sub.mat12) <- max.par + 1
                max.par <- max(sub.mat12, na.rm=TRUE)
                sub.mat13 <- matrix(NA, 4, 4)
                diag(sub.mat13) <- max.par + 1
                max.par <- max(sub.mat13, na.rm=TRUE)
                sub.mat14 <- matrix(NA, 4, 4)
                diag(sub.mat14) <- max.par + 1
                max.par <- max(sub.mat14, na.rm=TRUE)
                sub.mat15 <- matrix(NA, 4, 4)
                diag(sub.mat15) <- max.par + 1
                max.par <- max(sub.mat15, na.rm=TRUE)
                sub.mat16 <- matrix(NA, 4, 4)
                diag(sub.mat16) <- max.par + 1
                max.par <- max(sub.mat16, na.rm=TRUE)
                sub.mat17 <- matrix(NA, 4, 4)
                diag(sub.mat17) <- max.par + 1
                max.par <- max(sub.mat17, na.rm=TRUE)
                sub.mat18 <- matrix(NA, 4, 4)
                diag(sub.mat18) <- max.par + 1
                max.par <- max(sub.mat18, na.rm=TRUE)
                sub.mat19 <- matrix(NA, 4, 4)
                diag(sub.mat19) <- max.par + 1
                max.par <- max(sub.mat19, na.rm=TRUE)
                sub.mat20 <- matrix(NA, 4, 4)
                diag(sub.mat20) <- max.par + 1
                max.par <- max(sub.mat20, na.rm=TRUE)
                sub.mat21 <- matrix(NA, 4, 4)
                diag(sub.mat21) <- max.par + 1
                max.par <- max(sub.mat21, na.rm=TRUE)
                sub.mat22 <- matrix(NA, 4, 4)
                diag(sub.mat22) <- max.par + 1
                max.par <- max(sub.mat22, na.rm=TRUE)
                sub.mat23 <- matrix(NA, 4, 4)
                diag(sub.mat23) <- max.par + 1
                max.par <- max(sub.mat23, na.rm=TRUE)
                sub.mat24 <- matrix(NA, 4, 4)
                diag(sub.mat24) <- max.par + 1
                max.par <- max(sub.mat24, na.rm=TRUE)
                sub.mat25 <- matrix(NA, 4, 4)
                diag(sub.mat25) <- max.par + 1
                max.par <- max(sub.mat25, na.rm=TRUE)
                sub.mat26 <- matrix(NA, 4, 4)
                diag(sub.mat26) <- max.par + 1
                max.par <- max(sub.mat26, na.rm=TRUE)
                sub.mat27 <- matrix(NA, 4, 4)
                diag(sub.mat27) <- max.par + 1
                max.par <- max(sub.mat27, na.rm=TRUE)
                sub.mat28 <- matrix(NA, 4, 4)
                diag(sub.mat28) <- max.par + 1
                max.par <- max(sub.mat28, na.rm=TRUE)
                sub.mat29 <- matrix(NA, 4, 4)
                diag(sub.mat29) <- max.par + 1
                max.par <- max(sub.mat29, na.rm=TRUE)
                sub.mat30 <- matrix(NA, 4, 4)
                diag(sub.mat30) <- max.par + 1
                max.par <- max(sub.mat30, na.rm=TRUE)
                sub.mat31 <- matrix(NA, 4, 4)
                diag(sub.mat31) <- max.par + 1
                max.par <- max(sub.mat31, na.rm=TRUE)
                sub.mat32 <- matrix(NA, 4, 4)
                diag(sub.mat32) <- max.par + 1
                max.par <- max(sub.mat32, na.rm=TRUE)
                sub.mat33 <- matrix(NA, 4, 4)
                diag(sub.mat33) <- max.par + 1
                max.par <- max(sub.mat33, na.rm=TRUE)
                sub.mat34 <- matrix(NA, 4, 4)
                diag(sub.mat34) <- max.par + 1
                max.par <- max(sub.mat34, na.rm=TRUE)
                sub.mat35 <- matrix(NA, 4, 4)
                diag(sub.mat35) <- max.par + 1
                max.par <- max(sub.mat35, na.rm=TRUE)
                sub.mat36 <- matrix(NA, 4, 4)
                diag(sub.mat36) <- max.par + 1
                max.par <- max(sub.mat36, na.rm=TRUE)
                sub.mat37 <- matrix(NA, 4, 4)
                diag(sub.mat37) <- max.par + 1
                max.par <- max(sub.mat37, na.rm=TRUE)
                sub.mat38 <- matrix(NA, 4, 4)
                diag(sub.mat38) <- max.par + 1
                max.par <- max(sub.mat38, na.rm=TRUE)
                sub.mat39 <- matrix(NA, 4, 4)
                diag(sub.mat39) <- max.par + 1
                max.par <- max(sub.mat39, na.rm=TRUE)
                sub.mat40 <- matrix(NA, 4, 4)
                diag(sub.mat40) <- max.par + 1
                max.par <- max(sub.mat40, na.rm=TRUE)
                sub.mat41 <- matrix(NA, 4, 4)
                diag(sub.mat41) <- max.par + 1
                max.par <- max(sub.mat41, na.rm=TRUE)
                sub.mat42 <- matrix(NA, 4, 4)
                diag(sub.mat42) <- max.par + 1
                max.par <- max(sub.mat42, na.rm=TRUE)
                sub.mat43 <- matrix(NA, 4, 4)
                diag(sub.mat43) <- max.par + 1
                max.par <- max(sub.mat43, na.rm=TRUE)
                sub.mat44 <- matrix(NA, 4, 4)
                diag(sub.mat44) <- max.par + 1
                max.par <- max(sub.mat44, na.rm=TRUE)
                sub.mat44 <- matrix(NA, 4, 4)
                diag(sub.mat44) <- max.par + 1
                max.par <- max(sub.mat44, na.rm=TRUE)
                sub.mat45 <- matrix(NA, 4, 4)
                diag(sub.mat45) <- max.par + 1
                max.par <- max(sub.mat45, na.rm=TRUE)
                sub.mat46 <- matrix(NA, 4, 4)
                diag(sub.mat46) <- max.par + 1
                max.par <- max(sub.mat46, na.rm=TRUE)
                sub.mat47 <- matrix(NA, 4, 4)
                diag(sub.mat47) <- max.par + 1
                max.par <- max(sub.mat47, na.rm=TRUE)
                sub.mat48 <- matrix(NA, 4, 4)
                diag(sub.mat48) <- max.par + 1
                max.par <- max(sub.mat48, na.rm=TRUE)
                sub.mat49 <- matrix(NA, 4, 4)
                diag(sub.mat49) <- max.par + 1
                if(cat.trans.vary == TRUE){
                    trans.mat <- rbind(cbind(sub.mat1, sub.mat14, sub.mat25, sub.mat34, sub.mat41, sub.mat46, sub.mat49), cbind(sub.mat8, sub.mat2, sub.mat15, sub.mat26, sub.mat35, sub.mat42, sub.mat47), cbind(sub.mat20, sub.mat9, sub.mat3, sub.mat16, sub.mat27, sub.mat36, sub.mat43), cbind(sub.mat30, sub.mat21, sub.mat10, sub.mat4, sub.mat17, sub.mat28, sub.mat37),  cbind(sub.mat38, sub.mat31, sub.mat22, sub.mat11, sub.mat5, sub.mat18, sub.mat29), cbind(sub.mat44, sub.mat39, sub.mat32, sub.mat23, sub.mat12, sub.mat6, sub.mat19), cbind(sub.mat48, sub.mat45, sub.mat40, sub.mat33, sub.mat24, sub.mat13, sub.mat7))
                }else{
                    trans.mat <- rbind(cbind(sub.mat1, sub.mat8, sub.mat8, sub.mat8, sub.mat8, sub.mat8, sub.mat8), cbind(sub.mat8, sub.mat2, sub.mat8, sub.mat8, sub.mat8, sub.mat8, sub.mat8), cbind(sub.mat8, sub.mat8, sub.mat3, sub.mat8, sub.mat8, sub.mat8, sub.mat8), cbind(sub.mat8, sub.mat8, sub.mat8, sub.mat4, sub.mat8, sub.mat8, sub.mat8),  cbind(sub.mat8, sub.mat8, sub.mat8, sub.mat8, sub.mat5, sub.mat8, sub.mat8), cbind(sub.mat8, sub.mat8, sub.mat8, sub.mat8, sub.mat8, sub.mat6, sub.mat8), cbind(sub.mat8, sub.mat8, sub.mat8, sub.mat8, sub.mat8, sub.mat8, sub.mat7))
                }
            }
            rownames(trans.mat) <- colnames(trans.mat) <- c("(00A)","(01A)","(10A)","(11A)", "(00B)","(01B)","(10B)","(11B)", "(00C)","(01C)","(10C)","(11C)", "(00D)","(01D)","(10D)","(11D)", "(00E)","(01E)","(10E)","(11E)", "(00F)","(01F)","(10F)","(11F)", "(00G)","(01G)","(10G)","(11G)")
        }
        
        if(hidden.traits == 7){
            if(make.null == TRUE){
                sub.mat1 <- TransMatMuSSEsingle(cat.number=1, include.diagonals=include.diagonals)
                max.par <- max(sub.mat1, na.rm=TRUE)
                sub.mat2 <- matrix(NA, 4, 4)
                diag(sub.mat2) <- max.par + 1
                max.par <- max(sub.mat2, na.rm=TRUE)
                sub.mat3 <- matrix(NA, 4, 4)
                diag(sub.mat3) <- max.par + 1
                max.par <- max(sub.mat3, na.rm=TRUE)
                sub.mat4 <- matrix(NA, 4, 4)
                diag(sub.mat4) <- max.par + 1
                max.par <- max(sub.mat4, na.rm=TRUE)
                sub.mat5 <- matrix(NA, 4, 4)
                diag(sub.mat5) <- max.par + 1
                max.par <- max(sub.mat5, na.rm=TRUE)
                sub.mat6 <- matrix(NA, 4, 4)
                diag(sub.mat6) <- max.par + 1
                max.par <- max(sub.mat6, na.rm=TRUE)
                sub.mat7 <- matrix(NA, 4, 4)
                diag(sub.mat7) <- max.par + 1
                max.par <- max(sub.mat7, na.rm=TRUE)
                sub.mat8 <- matrix(NA, 4, 4)
                diag(sub.mat8) <- max.par + 1
                max.par <- max(sub.mat8, na.rm=TRUE)
                sub.mat9 <- matrix(NA, 4, 4)
                diag(sub.mat9) <- max.par + 1
                max.par <- max(sub.mat9, na.rm=TRUE)
                sub.mat10 <- matrix(NA, 4, 4)
                diag(sub.mat10) <- max.par + 1
                max.par <- max(sub.mat10, na.rm=TRUE)
                sub.mat11 <- matrix(NA, 4, 4)
                diag(sub.mat11) <- max.par + 1
                max.par <- max(sub.mat11, na.rm=TRUE)
                sub.mat12 <- matrix(NA, 4, 4)
                diag(sub.mat12) <- max.par + 1
                max.par <- max(sub.mat12, na.rm=TRUE)
                sub.mat13 <- matrix(NA, 4, 4)
                diag(sub.mat13) <- max.par + 1
                max.par <- max(sub.mat13, na.rm=TRUE)
                sub.mat14 <- matrix(NA, 4, 4)
                diag(sub.mat14) <- max.par + 1
                max.par <- max(sub.mat14, na.rm=TRUE)
                sub.mat15 <- matrix(NA, 4, 4)
                diag(sub.mat15) <- max.par + 1
                max.par <- max(sub.mat15, na.rm=TRUE)
                sub.mat16 <- matrix(NA, 4, 4)
                diag(sub.mat16) <- max.par + 1
                max.par <- max(sub.mat16, na.rm=TRUE)
                sub.mat17 <- matrix(NA, 4, 4)
                diag(sub.mat17) <- max.par + 1
                max.par <- max(sub.mat17, na.rm=TRUE)
                sub.mat18 <- matrix(NA, 4, 4)
                diag(sub.mat18) <- max.par + 1
                max.par <- max(sub.mat18, na.rm=TRUE)
                sub.mat19 <- matrix(NA, 4, 4)
                diag(sub.mat19) <- max.par + 1
                max.par <- max(sub.mat19, na.rm=TRUE)
                sub.mat20 <- matrix(NA, 4, 4)
                diag(sub.mat20) <- max.par + 1
                max.par <- max(sub.mat20, na.rm=TRUE)
                sub.mat21 <- matrix(NA, 4, 4)
                diag(sub.mat21) <- max.par + 1
                max.par <- max(sub.mat21, na.rm=TRUE)
                sub.mat22 <- matrix(NA, 4, 4)
                diag(sub.mat22) <- max.par + 1
                max.par <- max(sub.mat22, na.rm=TRUE)
                sub.mat23 <- matrix(NA, 4, 4)
                diag(sub.mat23) <- max.par + 1
                max.par <- max(sub.mat23, na.rm=TRUE)
                sub.mat24 <- matrix(NA, 4, 4)
                diag(sub.mat24) <- max.par + 1
                max.par <- max(sub.mat24, na.rm=TRUE)
                sub.mat25 <- matrix(NA, 4, 4)
                diag(sub.mat25) <- max.par + 1
                max.par <- max(sub.mat25, na.rm=TRUE)
                sub.mat26 <- matrix(NA, 4, 4)
                diag(sub.mat26) <- max.par + 1
                max.par <- max(sub.mat26, na.rm=TRUE)
                sub.mat27 <- matrix(NA, 4, 4)
                diag(sub.mat27) <- max.par + 1
                max.par <- max(sub.mat27, na.rm=TRUE)
                sub.mat28 <- matrix(NA, 4, 4)
                diag(sub.mat28) <- max.par + 1
                max.par <- max(sub.mat28, na.rm=TRUE)
                sub.mat29 <- matrix(NA, 4, 4)
                diag(sub.mat29) <- max.par + 1
                max.par <- max(sub.mat29, na.rm=TRUE)
                sub.mat30 <- matrix(NA, 4, 4)
                diag(sub.mat30) <- max.par + 1
                max.par <- max(sub.mat30, na.rm=TRUE)
                sub.mat31 <- matrix(NA, 4, 4)
                diag(sub.mat31) <- max.par + 1
                max.par <- max(sub.mat31, na.rm=TRUE)
                sub.mat32 <- matrix(NA, 4, 4)
                diag(sub.mat32) <- max.par + 1
                max.par <- max(sub.mat32, na.rm=TRUE)
                sub.mat33 <- matrix(NA, 4, 4)
                diag(sub.mat33) <- max.par + 1
                max.par <- max(sub.mat33, na.rm=TRUE)
                sub.mat34 <- matrix(NA, 4, 4)
                diag(sub.mat34) <- max.par + 1
                max.par <- max(sub.mat34, na.rm=TRUE)
                sub.mat35 <- matrix(NA, 4, 4)
                diag(sub.mat35) <- max.par + 1
                max.par <- max(sub.mat35, na.rm=TRUE)
                sub.mat36 <- matrix(NA, 4, 4)
                diag(sub.mat36) <- max.par + 1
                max.par <- max(sub.mat36, na.rm=TRUE)
                sub.mat37 <- matrix(NA, 4, 4)
                diag(sub.mat37) <- max.par + 1
                max.par <- max(sub.mat37, na.rm=TRUE)
                sub.mat38 <- matrix(NA, 4, 4)
                diag(sub.mat38) <- max.par + 1
                max.par <- max(sub.mat38, na.rm=TRUE)
                sub.mat39 <- matrix(NA, 4, 4)
                diag(sub.mat39) <- max.par + 1
                max.par <- max(sub.mat39, na.rm=TRUE)
                sub.mat40 <- matrix(NA, 4, 4)
                diag(sub.mat40) <- max.par + 1
                max.par <- max(sub.mat40, na.rm=TRUE)
                sub.mat41 <- matrix(NA, 4, 4)
                diag(sub.mat41) <- max.par + 1
                max.par <- max(sub.mat41, na.rm=TRUE)
                sub.mat42 <- matrix(NA, 4, 4)
                diag(sub.mat42) <- max.par + 1
                max.par <- max(sub.mat42, na.rm=TRUE)
                sub.mat43 <- matrix(NA, 4, 4)
                diag(sub.mat43) <- max.par + 1
                max.par <- max(sub.mat43, na.rm=TRUE)
                sub.mat44 <- matrix(NA, 4, 4)
                diag(sub.mat44) <- max.par + 1
                max.par <- max(sub.mat44, na.rm=TRUE)
                sub.mat44 <- matrix(NA, 4, 4)
                diag(sub.mat44) <- max.par + 1
                max.par <- max(sub.mat44, na.rm=TRUE)
                sub.mat45 <- matrix(NA, 4, 4)
                diag(sub.mat45) <- max.par + 1
                max.par <- max(sub.mat45, na.rm=TRUE)
                sub.mat46 <- matrix(NA, 4, 4)
                diag(sub.mat46) <- max.par + 1
                max.par <- max(sub.mat46, na.rm=TRUE)
                sub.mat47 <- matrix(NA, 4, 4)
                diag(sub.mat47) <- max.par + 1
                max.par <- max(sub.mat47, na.rm=TRUE)
                sub.mat48 <- matrix(NA, 4, 4)
                diag(sub.mat48) <- max.par + 1
                max.par <- max(sub.mat48, na.rm=TRUE)
                sub.mat49 <- matrix(NA, 4, 4)
                diag(sub.mat49) <- max.par + 1
                max.par <- max(sub.mat49, na.rm=TRUE)
                sub.mat50 <- matrix(NA, 4, 4)
                diag(sub.mat50) <- max.par + 1
                max.par <- max(sub.mat50, na.rm=TRUE)
                sub.mat51 <- matrix(NA, 4, 4)
                diag(sub.mat51) <- max.par + 1
                max.par <- max(sub.mat51, na.rm=TRUE)
                sub.mat52 <- matrix(NA, 4, 4)
                diag(sub.mat52) <- max.par + 1
                max.par <- max(sub.mat52, na.rm=TRUE)
                sub.mat53 <- matrix(NA, 4, 4)
                diag(sub.mat53) <- max.par + 1
                max.par <- max(sub.mat53, na.rm=TRUE)
                sub.mat54 <- matrix(NA, 4, 4)
                diag(sub.mat54) <- max.par + 1
                max.par <- max(sub.mat54, na.rm=TRUE)
                sub.mat55 <- matrix(NA, 4, 4)
                diag(sub.mat55) <- max.par + 1
                max.par <- max(sub.mat55, na.rm=TRUE)
                sub.mat56 <- matrix(NA, 4, 4)
                diag(sub.mat56) <- max.par + 1
                max.par <- max(sub.mat56, na.rm=TRUE)
                sub.mat57 <- matrix(NA, 4, 4)
                diag(sub.mat57) <- max.par + 1
                if(cat.trans.vary == TRUE){
                    trans.mat <- rbind(cbind(sub.mat1, sub.mat9, sub.mat22, sub.mat33, sub.mat42, sub.mat49, sub.mat54, sub.mat57), cbind(sub.mat2, sub.mat1, sub.mat10, sub.mat23, sub.mat34, sub.mat43, sub.mat50, sub.mat55), cbind(sub.mat16, sub.mat3, sub.mat1, sub.mat11, sub.mat24, sub.mat35, sub.mat44, sub.mat51), cbind(sub.mat28, sub.mat17, sub.mat4, sub.mat1, sub.mat12, sub.mat25, sub.mat36, sub.mat45), cbind(sub.mat38, sub.mat29, sub.mat18, sub.mat5, sub.mat1, sub.mat13, sub.mat26, sub.mat37), cbind(sub.mat46, sub.mat39, sub.mat30, sub.mat19, sub.mat6, sub.mat1, sub.mat14, sub.mat27), cbind(sub.mat52, sub.mat47, sub.mat40, sub.mat31, sub.mat20, sub.mat7, sub.mat1, sub.mat15), cbind(sub.mat56, sub.mat53, sub.mat48, sub.mat41, sub.mat32, sub.mat21, sub.mat8, sub.mat1))
                }else{
                    trans.mat <- rbind(cbind(sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2), cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1))
                }
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
                max.par <- max(sub.mat9, na.rm=TRUE)
                sub.mat10 <- matrix(NA, 4, 4)
                diag(sub.mat10) <- max.par + 1
                max.par <- max(sub.mat10, na.rm=TRUE)
                sub.mat11 <- matrix(NA, 4, 4)
                diag(sub.mat11) <- max.par + 1
                max.par <- max(sub.mat11, na.rm=TRUE)
                sub.mat12 <- matrix(NA, 4, 4)
                diag(sub.mat12) <- max.par + 1
                max.par <- max(sub.mat12, na.rm=TRUE)
                sub.mat13 <- matrix(NA, 4, 4)
                diag(sub.mat13) <- max.par + 1
                max.par <- max(sub.mat13, na.rm=TRUE)
                sub.mat14 <- matrix(NA, 4, 4)
                diag(sub.mat14) <- max.par + 1
                max.par <- max(sub.mat14, na.rm=TRUE)
                sub.mat15 <- matrix(NA, 4, 4)
                diag(sub.mat15) <- max.par + 1
                max.par <- max(sub.mat15, na.rm=TRUE)
                sub.mat16 <- matrix(NA, 4, 4)
                diag(sub.mat16) <- max.par + 1
                max.par <- max(sub.mat16, na.rm=TRUE)
                sub.mat17 <- matrix(NA, 4, 4)
                diag(sub.mat17) <- max.par + 1
                max.par <- max(sub.mat17, na.rm=TRUE)
                sub.mat18 <- matrix(NA, 4, 4)
                diag(sub.mat18) <- max.par + 1
                max.par <- max(sub.mat18, na.rm=TRUE)
                sub.mat19 <- matrix(NA, 4, 4)
                diag(sub.mat19) <- max.par + 1
                max.par <- max(sub.mat19, na.rm=TRUE)
                sub.mat20 <- matrix(NA, 4, 4)
                diag(sub.mat20) <- max.par + 1
                max.par <- max(sub.mat20, na.rm=TRUE)
                sub.mat21 <- matrix(NA, 4, 4)
                diag(sub.mat21) <- max.par + 1
                max.par <- max(sub.mat21, na.rm=TRUE)
                sub.mat22 <- matrix(NA, 4, 4)
                diag(sub.mat22) <- max.par + 1
                max.par <- max(sub.mat22, na.rm=TRUE)
                sub.mat23 <- matrix(NA, 4, 4)
                diag(sub.mat23) <- max.par + 1
                max.par <- max(sub.mat23, na.rm=TRUE)
                sub.mat24 <- matrix(NA, 4, 4)
                diag(sub.mat24) <- max.par + 1
                max.par <- max(sub.mat24, na.rm=TRUE)
                sub.mat25 <- matrix(NA, 4, 4)
                diag(sub.mat25) <- max.par + 1
                max.par <- max(sub.mat25, na.rm=TRUE)
                sub.mat26 <- matrix(NA, 4, 4)
                diag(sub.mat26) <- max.par + 1
                max.par <- max(sub.mat26, na.rm=TRUE)
                sub.mat27 <- matrix(NA, 4, 4)
                diag(sub.mat27) <- max.par + 1
                max.par <- max(sub.mat27, na.rm=TRUE)
                sub.mat28 <- matrix(NA, 4, 4)
                diag(sub.mat28) <- max.par + 1
                max.par <- max(sub.mat28, na.rm=TRUE)
                sub.mat29 <- matrix(NA, 4, 4)
                diag(sub.mat29) <- max.par + 1
                max.par <- max(sub.mat29, na.rm=TRUE)
                sub.mat30 <- matrix(NA, 4, 4)
                diag(sub.mat30) <- max.par + 1
                max.par <- max(sub.mat30, na.rm=TRUE)
                sub.mat31 <- matrix(NA, 4, 4)
                diag(sub.mat31) <- max.par + 1
                max.par <- max(sub.mat31, na.rm=TRUE)
                sub.mat32 <- matrix(NA, 4, 4)
                diag(sub.mat32) <- max.par + 1
                max.par <- max(sub.mat32, na.rm=TRUE)
                sub.mat33 <- matrix(NA, 4, 4)
                diag(sub.mat33) <- max.par + 1
                max.par <- max(sub.mat33, na.rm=TRUE)
                sub.mat34 <- matrix(NA, 4, 4)
                diag(sub.mat34) <- max.par + 1
                max.par <- max(sub.mat34, na.rm=TRUE)
                sub.mat35 <- matrix(NA, 4, 4)
                diag(sub.mat35) <- max.par + 1
                max.par <- max(sub.mat35, na.rm=TRUE)
                sub.mat36 <- matrix(NA, 4, 4)
                diag(sub.mat36) <- max.par + 1
                max.par <- max(sub.mat36, na.rm=TRUE)
                sub.mat37 <- matrix(NA, 4, 4)
                diag(sub.mat37) <- max.par + 1
                max.par <- max(sub.mat37, na.rm=TRUE)
                sub.mat38 <- matrix(NA, 4, 4)
                diag(sub.mat38) <- max.par + 1
                max.par <- max(sub.mat38, na.rm=TRUE)
                sub.mat39 <- matrix(NA, 4, 4)
                diag(sub.mat39) <- max.par + 1
                max.par <- max(sub.mat39, na.rm=TRUE)
                sub.mat40 <- matrix(NA, 4, 4)
                diag(sub.mat40) <- max.par + 1
                max.par <- max(sub.mat40, na.rm=TRUE)
                sub.mat41 <- matrix(NA, 4, 4)
                diag(sub.mat41) <- max.par + 1
                max.par <- max(sub.mat41, na.rm=TRUE)
                sub.mat42 <- matrix(NA, 4, 4)
                diag(sub.mat42) <- max.par + 1
                max.par <- max(sub.mat42, na.rm=TRUE)
                sub.mat43 <- matrix(NA, 4, 4)
                diag(sub.mat43) <- max.par + 1
                max.par <- max(sub.mat43, na.rm=TRUE)
                sub.mat44 <- matrix(NA, 4, 4)
                diag(sub.mat44) <- max.par + 1
                max.par <- max(sub.mat44, na.rm=TRUE)
                sub.mat44 <- matrix(NA, 4, 4)
                diag(sub.mat44) <- max.par + 1
                max.par <- max(sub.mat44, na.rm=TRUE)
                sub.mat45 <- matrix(NA, 4, 4)
                diag(sub.mat45) <- max.par + 1
                max.par <- max(sub.mat45, na.rm=TRUE)
                sub.mat46 <- matrix(NA, 4, 4)
                diag(sub.mat46) <- max.par + 1
                max.par <- max(sub.mat46, na.rm=TRUE)
                sub.mat47 <- matrix(NA, 4, 4)
                diag(sub.mat47) <- max.par + 1
                max.par <- max(sub.mat47, na.rm=TRUE)
                sub.mat48 <- matrix(NA, 4, 4)
                diag(sub.mat48) <- max.par + 1
                max.par <- max(sub.mat48, na.rm=TRUE)
                sub.mat49 <- matrix(NA, 4, 4)
                diag(sub.mat49) <- max.par + 1
                max.par <- max(sub.mat49, na.rm=TRUE)
                sub.mat50 <- matrix(NA, 4, 4)
                diag(sub.mat50) <- max.par + 1
                max.par <- max(sub.mat50, na.rm=TRUE)
                sub.mat51 <- matrix(NA, 4, 4)
                diag(sub.mat51) <- max.par + 1
                max.par <- max(sub.mat51, na.rm=TRUE)
                sub.mat52 <- matrix(NA, 4, 4)
                diag(sub.mat52) <- max.par + 1
                max.par <- max(sub.mat52, na.rm=TRUE)
                sub.mat53 <- matrix(NA, 4, 4)
                diag(sub.mat53) <- max.par + 1
                max.par <- max(sub.mat53, na.rm=TRUE)
                sub.mat54 <- matrix(NA, 4, 4)
                diag(sub.mat54) <- max.par + 1
                max.par <- max(sub.mat54, na.rm=TRUE)
                sub.mat55 <- matrix(NA, 4, 4)
                diag(sub.mat55) <- max.par + 1
                max.par <- max(sub.mat55, na.rm=TRUE)
                sub.mat56 <- matrix(NA, 4, 4)
                diag(sub.mat56) <- max.par + 1
                max.par <- max(sub.mat56, na.rm=TRUE)
                sub.mat57 <- matrix(NA, 4, 4)
                diag(sub.mat57) <- max.par + 1
                max.par <- max(sub.mat57, na.rm=TRUE)
                sub.mat58 <- matrix(NA, 4, 4)
                diag(sub.mat58) <- max.par + 1
                max.par <- max(sub.mat58, na.rm=TRUE)
                sub.mat59 <- matrix(NA, 4, 4)
                diag(sub.mat59) <- max.par + 1
                max.par <- max(sub.mat59, na.rm=TRUE)
                sub.mat60 <- matrix(NA, 4, 4)
                diag(sub.mat60) <- max.par + 1
                max.par <- max(sub.mat60, na.rm=TRUE)
                sub.mat61 <- matrix(NA, 4, 4)
                diag(sub.mat61) <- max.par + 1
                max.par <- max(sub.mat61, na.rm=TRUE)
                sub.mat62 <- matrix(NA, 4, 4)
                diag(sub.mat62) <- max.par + 1
                max.par <- max(sub.mat62, na.rm=TRUE)
                sub.mat63 <- matrix(NA, 4, 4)
                diag(sub.mat63) <- max.par + 1
                max.par <- max(sub.mat63, na.rm=TRUE)
                sub.mat64 <- matrix(NA, 4, 4)
                diag(sub.mat64) <- max.par + 1
                if(cat.trans.vary == TRUE){
                    trans.mat <- rbind(cbind(sub.mat1, sub.mat16, sub.mat29, sub.mat40, sub.mat49, sub.mat56, sub.mat61, sub.mat64), cbind(sub.mat9, sub.mat2, sub.mat17, sub.mat30, sub.mat41, sub.mat50, sub.mat57, sub.mat62), cbind(sub.mat23, sub.mat10, sub.mat3, sub.mat18, sub.mat31, sub.mat42, sub.mat51, sub.mat58), cbind(sub.mat35, sub.mat24, sub.mat11, sub.mat4, sub.mat19, sub.mat32, sub.mat43, sub.mat52),  cbind(sub.mat45, sub.mat36, sub.mat25, sub.mat12, sub.mat5, sub.mat20, sub.mat33, sub.mat44), cbind(sub.mat53, sub.mat46, sub.mat37, sub.mat26, sub.mat13, sub.mat6, sub.mat21, sub.mat34), cbind(sub.mat59, sub.mat54, sub.mat47, sub.mat38, sub.mat27, sub.mat14, sub.mat7, sub.mat22), cbind(sub.mat63, sub.mat60, sub.mat55, sub.mat48, sub.mat39, sub.mat28, sub.mat15, sub.mat8))
                }else{
                    trans.mat <- rbind(cbind(sub.mat1, sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat9), cbind(sub.mat9, sub.mat2, sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat9), cbind(sub.mat9, sub.mat9, sub.mat3, sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat9), cbind(sub.mat9, sub.mat9, sub.mat9, sub.mat4, sub.mat9, sub.mat9, sub.mat9, sub.mat9),  cbind(sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat5, sub.mat9, sub.mat9, sub.mat9), cbind(sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat6, sub.mat9, sub.mat9), cbind(sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat7, sub.mat9), cbind(sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat9, sub.mat8))
                }
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




