
######################################################################################################################################
######################################################################################################################################
### TransMatMaker -- Builds transition rate matrix for GeoHiSSE and our special MuSSE model
## Daniel: Code updated to deal with GeoHiSSE models with more than two endemic areas.
######################################################################################################################################
######################################################################################################################################

TransMatMakerGeoHiSSE <- function(hidden.areas=0, endemic.areas=2, make.null=FALSE, include.jumps=FALSE, separate.extirpation=FALSE){
    ## Check parameters of the function and return human readable errors:
    if( endemic.areas > 3 ){
        stop("At the moment just up to three endemic areas are supported.")
    }
    if( hidden.areas > 4 & endemic.areas == 2 ){
        stop("The GeoHiSSE model with 2 endemic areas supports up to 4 hidden states.")
    }
    if( hidden.areas > 12 & endemic.areas == 3 ){
        stop("The GeoHiSSE model with 3 endemic areas supports up to 12 hidden states.")
    }

    if( endemic.areas > 2){
        ## Then call the function to deal with the larger model.
        trans.mat <- TransMatGeoSSE_Plus(areas.number=endemic.areas, cat.number=hidden.areas+1
                                       , include.jumps=include.jumps, separate.extirpation=separate.extirpation
                                       , make.null=make.null)
        return( trans.mat )
    }
    
    if(hidden.areas == 0){
        trans.mat <- TransMatGeoSSEsingle(cat.number=1, include.jumps=include.jumps
                                        , separate.extirpation=separate.extirpation)
    }else{
        if(hidden.areas == 1){
            if(make.null == TRUE){
                sub.mat1 <- TransMatGeoSSEsingle(cat.number=1, include.jumps=include.jumps
                                               , separate.extirpation=separate.extirpation)
                sub.mat2 <- TransMatGeoSSEsingle(cat.number=1, include.jumps=include.jumps
                                               , separate.extirpation=separate.extirpation)
                max.par <- max(sub.mat2, na.rm=TRUE)
                sub.mat3 <- matrix(NA, 3, 3)
                diag(sub.mat3) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1,sub.mat3), cbind(sub.mat3,sub.mat2))
            }else{
                sub.mat1 <- TransMatGeoSSEsingle(cat.number=1, include.jumps=include.jumps
                                               , separate.extirpation=separate.extirpation)
                sub.mat2 <- TransMatGeoSSEsingle(cat.number=2, include.jumps=include.jumps
                                               , separate.extirpation=separate.extirpation)
                max.par <- max(sub.mat2, na.rm=TRUE)
                sub.mat3 <- matrix(NA, 3, 3)
                diag(sub.mat3) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat3), cbind(sub.mat3, sub.mat2))
            }
            rownames(trans.mat) <- colnames(trans.mat) <- c("(0A)","(1A)","(01A)", "(0B)","(1B)","(01B)")
        }

        if(hidden.areas == 2){
            if(make.null == TRUE){
                sub.mat1 <- TransMatGeoSSEsingle(cat.number=1, include.jumps=include.jumps
                                               , separate.extirpation=separate.extirpation)
                max.par <- max(sub.mat1, na.rm=TRUE)
                sub.mat2 <- matrix(NA, 3, 3)
                diag(sub.mat2) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat2, sub.mat2), cbind(sub.mat2, sub.mat1, sub.mat2)
                                 , cbind(sub.mat2, sub.mat2, sub.mat1))
            }else{
                sub.mat1 <- TransMatGeoSSEsingle(cat.number=1, include.jumps=include.jumps
                                               , separate.extirpation=separate.extirpation)
                sub.mat2 <- TransMatGeoSSEsingle(cat.number=2, include.jumps=include.jumps
                                               , separate.extirpation=separate.extirpation)
                sub.mat3 <- TransMatGeoSSEsingle(cat.number=3, include.jumps=include.jumps
                                               , separate.extirpation=separate.extirpation)
                sub.mat4 <- matrix(NA, 3, 3)
                max.par <- max(sub.mat3, na.rm=TRUE)
                diag(sub.mat4) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat4, sub.mat4), cbind(sub.mat4, sub.mat2, sub.mat4)
                                 , cbind(sub.mat4, sub.mat4, sub.mat3))
            }
            rownames(trans.mat) <- colnames(trans.mat) <- c("(0A)","(1A)","(01A)", "(0B)","(1B)","(01B)", "(0C)","(1C)","(01C)")
        }

        if(hidden.areas == 3){
            if(make.null == TRUE){
                sub.mat1 <- TransMatGeoSSEsingle(cat.number=1, include.jumps=include.jumps
                                               , separate.extirpation=separate.extirpation)
                max.par <- max(sub.mat1, na.rm=TRUE)
                sub.mat2 <- matrix(NA, 3, 3)
                diag(sub.mat2) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat2, sub.mat2, sub.mat2)
                                 , cbind(sub.mat2, sub.mat1, sub.mat2, sub.mat2)
                                 , cbind(sub.mat2, sub.mat2, sub.mat1, sub.mat2)
                                 , cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat1))
            }else{
                sub.mat1 <- TransMatGeoSSEsingle(cat.number=1, include.jumps=include.jumps
                                               , separate.extirpation=separate.extirpation)
                sub.mat2 <- TransMatGeoSSEsingle(cat.number=2, include.jumps=include.jumps
                                               , separate.extirpation=separate.extirpation)
                sub.mat3 <- TransMatGeoSSEsingle(cat.number=3, include.jumps=include.jumps
                                               , separate.extirpation=separate.extirpation)
                sub.mat4 <- TransMatGeoSSEsingle(cat.number=4, include.jumps=include.jumps
                                               , separate.extirpation=separate.extirpation)
                sub.mat5 <- matrix(NA, 3, 3)
                max.par <- max(sub.mat4, na.rm=TRUE)
                diag(sub.mat5) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat5, sub.mat5, sub.mat5)
                                 , cbind(sub.mat5, sub.mat2, sub.mat5, sub.mat5)
                                 , cbind(sub.mat5, sub.mat5, sub.mat3, sub.mat5)
                                 , cbind(sub.mat5, sub.mat5, sub.mat5, sub.mat4))
            }
            rownames(trans.mat) <- colnames(trans.mat) <- c("(0A)","(1A)","(01A)", "(0B)","(1B)","(01B)", "(0C)","(1C)","(01C)", "(0D)","(1D)","(01D)")
        }
        
        if(hidden.areas == 4){
            if(make.null == TRUE){
                sub.mat1 <- TransMatGeoSSEsingle(cat.number=1, include.jumps=include.jumps
                                               , separate.extirpation=separate.extirpation)
                max.par <- max(sub.mat1, na.rm=TRUE)
                sub.mat2 <- matrix(NA, 3, 3)
                diag(sub.mat2) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat2, sub.mat2, sub.mat2, sub.mat2)
                                 , cbind(sub.mat2, sub.mat1, sub.mat2, sub.mat2, sub.mat2)
                                 , cbind(sub.mat2, sub.mat2, sub.mat1, sub.mat2, sub.mat2)
                                 , cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat1, sub.mat2)
                                 , cbind(sub.mat2, sub.mat2, sub.mat2, sub.mat2, sub.mat1))
            }else{
                sub.mat1 <- TransMatGeoSSEsingle(cat.number=1, include.jumps=include.jumps
                                               , separate.extirpation=separate.extirpation)
                sub.mat2 <- TransMatGeoSSEsingle(cat.number=2, include.jumps=include.jumps
                                               , separate.extirpation=separate.extirpation)
                sub.mat3 <- TransMatGeoSSEsingle(cat.number=3, include.jumps=include.jumps
                                               , separate.extirpation=separate.extirpation)
                sub.mat4 <- TransMatGeoSSEsingle(cat.number=4, include.jumps=include.jumps
                                               , separate.extirpation=separate.extirpation)
                sub.mat5 <- TransMatGeoSSEsingle(cat.number=5, include.jumps=include.jumps
                                               , separate.extirpation=separate.extirpation)
                sub.mat6 <- matrix(NA, 3, 3)
                max.par <- max(sub.mat5, na.rm=TRUE)
                diag(sub.mat6) <- max.par + 1
                trans.mat <- rbind(cbind(sub.mat1, sub.mat6, sub.mat6, sub.mat6, sub.mat6)
                                 , cbind(sub.mat6, sub.mat2, sub.mat6, sub.mat6, sub.mat6)
                                 , cbind(sub.mat6, sub.mat6, sub.mat3, sub.mat6, sub.mat6)
                                 , cbind(sub.mat6, sub.mat6, sub.mat6, sub.mat4, sub.mat6)
                                 , cbind(sub.mat6, sub.mat6, sub.mat6, sub.mat6, sub.mat5))
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

## Make function to produce the transition matrix when more than two endemic areas are present.
## This takes into account if the model is CID or not.
TransMatGeoSSE_Plus <- function(areas.number, cat.number, include.jumps, separate.extirpation, make.null){
    ## At the moment only areas.number == 3 is supported:
    if( areas.number > 3 ){
        stop( "Internal function supports only 3 endemic areas." )
    }

    if( make.null ){
        ## This case transitions have the same index across all layers.
        ## Can just construct the matrix without increasing the indexes.
        trans.mat <- matrix(0, nrow = 6*cat.number, ncol = 6*cat.number)
        diag( trans.mat ) <- NA
        geo.mat.names <- paste0( "(", paste0(c(1,2,3,12,13,23), rep(LETTERS[1:cat.number], each = 6)), ")" )
        rownames( trans.mat ) <- colnames( trans.mat ) <- geo.mat.names
        for( i in seq(from = 0, by = 6, length.out = cat.number) ){
            trans.mat[1+i,(4:5)+i] <- 1:2
            trans.mat[2+i,c(4,6)+i] <- 3:4
            trans.mat[3+i,c(5,6)+i] <- 5:6
        }
        if( cat.number > 1 ){
            np <- max( trans.mat, na.rm = TRUE )
            sub.mat <- diag(np+1, nrow = 6, ncol = 6)
            for( i in 0:(cat.number-1) ){
                for( j in 0:(cat.number-1) ){
                    if( i != j ){
                        trans.mat[(1:6)+(i*6),(1:6)+(j*6)] <- sub.mat
                    }
                }
            }    
        }
        if( include.jumps ){
            ## Transitions between endemic areas are allowed in this case.
            ## For 3 endemic areas, the position of the endemic areas are fixed on the matrix.
            max_id <- max(trans.mat, na.rm = TRUE)
            id_vec_jumps <- 1:6 + max_id
            for( i in seq(from = 0, by = 6, length.out = cat.number) ){
                trans.mat[1+i,c(2,3)+i] <- id_vec_jumps[1:2]
                trans.mat[2+i,c(1,3)+i] <- id_vec_jumps[3:4]
                trans.mat[3+i,c(1,2)+i] <- id_vec_jumps[5:6]
                ## Get the updated max_id value.
            }
        }
        if( separate.extirpation ){
            ## Transitions between endemic areas are allowed in this case.
            ## For 3 endemic areas, the position of the endemic areas are fixed on the matrix.
            max_id <- max(trans.mat, na.rm = TRUE)
            id_vec_ext <- 1:6 + max_id
            for( i in seq(from = 0, by = 6, length.out = cat.number) ){
                trans.mat[4+i,c(1,2)+i] <- id_vec_ext[1:2]
                trans.mat[5+i,c(1,3)+i] <- id_vec_ext[3:4]
                trans.mat[6+i,c(2,3)+i] <- id_vec_ext[5:6]
                ## Get the updated max_id value.
                max_id <- max(trans.mat, na.rm = TRUE)
            }
        }
        return( trans.mat )
    } else{
        ## Here is the standard case when the indexes are different across the models.
        trans.mat <- matrix(0, nrow = 6*cat.number, ncol = 6*cat.number)
        diag( trans.mat ) <- NA
        geo.mat.names <- paste0( "(", paste0(c(1,2,3,12,13,23), rep(LETTERS[1:cat.number], each = 6)), ")" )
        rownames( trans.mat ) <- colnames( trans.mat ) <- geo.mat.names
        for( i in seq(from = 0, by = 6, length.out = cat.number) ){
            trans.mat[1+i,(4:5)+i] <- (1:2)+i
            trans.mat[2+i,c(4,6)+i] <- (3:4)+i
            trans.mat[3+i,c(5,6)+i] <- (5:6)+i
        }
        if( cat.number > 1 ){
            np <- max( trans.mat, na.rm = TRUE )
            sub.mat <- diag(np+1, nrow = 6, ncol = 6)
            for( i in 0:(cat.number-1) ){
                for( j in 0:(cat.number-1) ){
                    if( i != j ){
                        trans.mat[(1:6)+(i*6),(1:6)+(j*6)] <- sub.mat
                    }
                }
            }    
        }
        if( include.jumps ){
            ## Transitions between endemic areas are allowed in this case.
            ## For 3 endemic areas, the position of the endemic areas are fixed on the matrix.
            max_id <- max(trans.mat, na.rm = TRUE)
            for( i in seq(from = 0, by = 6, length.out = cat.number) ){
                id_vec_jumps <- 1:6 + max_id
                trans.mat[1+i,c(2,3)+i] <- id_vec_jumps[1:2]
                trans.mat[2+i,c(1,3)+i] <- id_vec_jumps[3:4]
                trans.mat[3+i,c(1,2)+i] <- id_vec_jumps[5:6]
                ## Get the updated max_id value.
                max_id <- max(trans.mat, na.rm = TRUE)
            }
        }
        if( separate.extirpation ){
            ## Transitions between endemic areas are allowed in this case.
            ## For 3 endemic areas, the position of the endemic areas are fixed on the matrix.
            max_id <- max(trans.mat, na.rm = TRUE)
            for( i in seq(from = 0, by = 6, length.out = cat.number) ){
                id_vec_ext <- 1:6 + max_id
                trans.mat[4+i,c(1,2)+i] <- id_vec_ext[1:2]
                trans.mat[5+i,c(1,3)+i] <- id_vec_ext[3:4]
                trans.mat[6+i,c(2,3)+i] <- id_vec_ext[5:6]
                ## Get the updated max_id value.
                max_id <- max(trans.mat, na.rm = TRUE)
            }
        }
        return( trans.mat )
    }
    
}
