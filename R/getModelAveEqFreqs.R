
######################################################################################################################################
######################################################################################################################################
### EquiFreq -- calculates a model-averaged equilibrium frequency for each state under a set of parameters in HiSSE and GeoHiSSE
######################################################################################################################################
######################################################################################################################################

GetModelAveEqFreqs <- function(x, max.time, model.type="hisse", get.rates=FALSE, rate.type="turnover", get.all.states=FALSE){

    ## Need to first check if max.time is present and return an error if not.
    if( missing(max.time) ){
        stop("argument 'max.time' is missing, with no default")
    }
    
    if(model.type == "hisse"){
        res <- c()
        hisse.results <- x

        ## The object class can be a vector.
        if( inherits(x=hisse.results, what="hisse.fit") ){
            ## we have to make a list so we can run this generally
            tmp.list <- list()
            tmp.list[[1]] <- hisse.results
            hisse.results <- tmp.list
        }
        for(model.index in 1:length(hisse.results)){
            if(inherits(x=hisse.results[[model.index]], what="hisse.fit")){
                ##Modify the data file
                data.new <- data.frame(hisse.results[[model.index]]$data[,2], hisse.results[[model.index]]$data[,2], row.names=hisse.results[[model.index]]$data[,1])
                data.new <- data.new[hisse.results[[model.index]]$phy$tip.label,]
                cache = ParametersToPass(hisse.results[[model.index]]$phy, data.new[,1], model.vec=hisse.results[[model.index]]$solution, f=hisse.results[[model.index]]$f, timeslice=NULL, hidden.states=TRUE)
                transformed.pars <- ParameterTransform(hisse.results[[model.index]]$solution[1:4], hisse.results[[model.index]]$solution[5:8])
                cache$lambda0 <- transformed.pars[1]
                cache$lambda1 <- transformed.pars[2]
                cache$lambdaA <- transformed.pars[3]
                cache$lambdaB <- transformed.pars[4]
                cache$death0 <- transformed.pars[5]
                cache$death1 <- transformed.pars[6]
                cache$deathA <- transformed.pars[7]
                cache$deathB <- transformed.pars[8]
                cache$turnover.beta.factor0 = 1 / dbeta(0.1, hisse.results[[model.index]]$solution[21], hisse.results[[model.index]]$solution[25])
                cache$turnover.beta.factor1 = 1 / dbeta(0.1, hisse.results[[model.index]]$solution[22], hisse.results[[model.index]]$solution[26])
                cache$turnover.beta.factorA = 1 / dbeta(0.1, hisse.results[[model.index]]$solution[23], hisse.results[[model.index]]$solution[27])
                cache$turnover.beta.factorB = 1 / dbeta(0.1, hisse.results[[model.index]]$solution[24], hisse.results[[model.index]]$solution[28])
                cache$eps.beta.factor0 = 1 / dbeta(0.1, hisse.results[[model.index]]$solution[29], hisse.results[[model.index]]$solution[33])
                cache$eps.beta.factor1 = 1 / dbeta(0.1, hisse.results[[model.index]]$solution[30], hisse.results[[model.index]]$solution[34])
                cache$eps.beta.factorA = 1 / dbeta(0.1, hisse.results[[model.index]]$solution[31], hisse.results[[model.index]]$solution[35])
                cache$eps.beta.factorB = 1 / dbeta(0.1, hisse.results[[model.index]]$solution[32], hisse.results[[model.index]]$solution[36])

                if(hisse.results[[model.index]]$root.type=="madfitz"){
                    get.starting.probs <- DownPass(phy=hisse.results[[model.index]]$phy, cache=cache, hidden.states=TRUE, condition.on.survival=hisse.results[[model.index]]$condition.on.survival, root.type=hisse.results[[model.index]]$root.type, root.p=hisse.results[[model.index]]$root.p, get.phi=TRUE, ode.eps=0)$compD.root
                }else{
                    get.starting.probs <- hisse.results[[model.index]]$root.p
                }
                out <- lsoda(c(state0A=get.starting.probs[1],state1A=get.starting.probs[2],state0B=get.starting.probs[3], state1B=get.starting.probs[4]), times=c(0, max.time), func=EqFreqHiSSE, parms=NULL, cache=cache, rtol=1e-8, atol=1e-8)[-1,-1]
                if(get.rates == TRUE){
                    rescaled.probs.0 <- c(out[1],out[3]) / sum(out[1],out[3])
                    rescaled.probs.1 <- c(out[2],out[4]) / sum(out[2],out[4])
                    if(rate.type == "turnover"){
                        state0 <- ((cache$lambda0+cache$death0) * rescaled.probs.0[1]) + ((cache$lambdaA+cache$deathA) * rescaled.probs.0[2])
                        state1 <- ((cache$lambda1+cache$death1) * rescaled.probs.1[1]) + ((cache$lambdaB+cache$deathB) * rescaled.probs.1[2])
                    }
                    if(rate.type == "net.div"){
                        state0 <- ((cache$lambda0-cache$death0) * rescaled.probs.0[1]) + ((cache$lambdaA-cache$deathA) * rescaled.probs.0[2])
                        state1 <- ((cache$lambda1-cache$death1) * rescaled.probs.1[1]) + ((cache$lambdaB-cache$deathB) * rescaled.probs.1[2])
                    }
                    if(rate.type == "speciation"){
                        state0 <- (cache$lambda0 * rescaled.probs.0[1]) + (cache$lambdaA * rescaled.probs.0[2])
                        state1 <- (cache$lambda1 * rescaled.probs.1[1]) + (cache$lambdaB * rescaled.probs.1[2])
                    }
                    out.mat <- matrix(c(state0, state1), 1, 2)
                    colnames(out.mat) <- c("0", "1")
                    res <- rbind(res, out.mat)
                }else{
                    if(get.all.states == TRUE){
                        out.mat <- t(matrix(out, 2, 2))
                        colnames(out.mat) <- c("0", "1")
                        rownames(out.mat) <- c("A", "B")
                        return(out.mat / sum(out.mat))
                    }else{
                        out.mat <- t(matrix(out, 2, 2))
                        colnames(out.mat) <- c("0", "1")
                        res <- rbind(res, colSums(out.mat)/sum(out.mat))
                    }
                }
            }else{
                data.new <- data.frame(hisse.results[[model.index]]$data[,2], hisse.results[[model.index]]$data[,2], row.names=hisse.results[[model.index]]$data[,1])
                data.new <- data.new[hisse.results[[model.index]]$phy$tip.label,]
                cache = ParametersToPassNull(hisse.results[[model.index]]$phy, data.new[,1], model.vec=hisse.results[[model.index]]$solution, f=hisse.results[[model.index]]$f)
                transformed.pars <- ParameterTransform(hisse.results[[model.index]]$solution[1:8], hisse.results[[model.index]]$solution[9:16])
                cache$lambda0A <- transformed.pars[1]
                cache$lambda0B <- transformed.pars[2]
                cache$lambda0C <- transformed.pars[3]
                cache$lambda0D <- transformed.pars[4]
                cache$lambda1A <- transformed.pars[5]
                cache$lambda1B <- transformed.pars[6]
                cache$lambda1C <- transformed.pars[7]
                cache$lambda1D <- transformed.pars[8]
                cache$death0A <- transformed.pars[9]
                cache$death0B <- transformed.pars[10]
                cache$death0C <- transformed.pars[11]
                cache$death0D <- transformed.pars[12]
                cache$death1A <- transformed.pars[13]
                cache$death1B <- transformed.pars[14]
                cache$death1C <- transformed.pars[15]
                cache$death1D <- transformed.pars[16]

                if(hisse.results[[model.index]]$root.type=="madfitz"){
                    get.starting.probs <- DownPassNull(phy=hisse.results[[model.index]]$phy, cache=cache, condition.on.survival=hisse.results[[model.index]]$condition.on.survival, root.type=hisse.results[[model.index]]$root.type, root.p=hisse.results[[model.index]]$root.p, get.phi=TRUE)$compD.root
                }else{
                    get.starting.probs <- hisse.results[[model.index]]$root.p
                }
                out <- lsoda(c(state0A=get.starting.probs[1],state1A=get.starting.probs[2],state0B=get.starting.probs[3],state1B=get.starting.probs[4],state0C=get.starting.probs[5],state1C=get.starting.probs[6],state0D=get.starting.probs[7],state1D=get.starting.probs[8]), times=c(0, max.time), func=EqFreqCID4, parms=NULL, cache=cache, rtol=1e-8, atol=1e-8)[-1,-1]
                if(get.rates == TRUE){
                    rescaled.probs.0 <- c(out[1],out[3],out[5],out[7]) / sum(out[1],out[3],out[5],out[7])
                    rescaled.probs.1 <- c(out[2],out[4],out[6],out[8]) / sum(out[2],out[4],out[6],out[8])
                    if(rate.type == "turnover"){
                        state0 <- ((cache$lambda0A+cache$death0A) * rescaled.probs.0[1]) + ((cache$lambda0B+cache$death0B) * rescaled.probs.0[2]) + ((cache$lambda0C+cache$death0C) * rescaled.probs.0[3]) + ((cache$lambda0D+cache$death0D) * rescaled.probs.0[4])
                        state1 <- ((cache$lambda1A+cache$death1A) * rescaled.probs.1[1]) + ((cache$lambda1B+cache$death1B) * rescaled.probs.1[2]) + ((cache$lambda1C+cache$death1C) * rescaled.probs.1[3]) + ((cache$lambda1D+cache$death1D) * rescaled.probs.1[4])
                    }
                    if(rate.type == "net.div"){
                        state0 <- ((cache$lambda0A-cache$death0A) * rescaled.probs.0[1]) + ((cache$lambda0B-cache$death0B) * rescaled.probs.0[2]) + ((cache$lambda0C-cache$death0C) * rescaled.probs.0[3]) + ((cache$lambda0D-cache$death0D) * rescaled.probs.0[4])
                        state1 <- ((cache$lambda1A-cache$death1A) * rescaled.probs.1[1]) + ((cache$lambda1B-cache$death1B) * rescaled.probs.1[2]) + ((cache$lambda1C-cache$death1C) * rescaled.probs.1[3]) + ((cache$lambda1D-cache$death1D) * rescaled.probs.1[4])
                    }
                    if(rate.type == "speciation"){
                        state0 <- (cache$lambda0A * rescaled.probs.0[1]) + (cache$lambda0B * rescaled.probs.0[2]) + (cache$lambda0C * rescaled.probs.0[3]) + (cache$lambda0D * rescaled.probs.0[4])
                        state1 <- (cache$lambda1A * rescaled.probs.1[1]) + (cache$lambda1B * rescaled.probs.1[2]) + (cache$lambda1C * rescaled.probs.1[3]) + (cache$lambda1D * rescaled.probs.1[4])
                    }
                    out.mat <- matrix(c(state0, state1), 1, 2)
                    colnames(out.mat) <- c("0", "1")
                    res <- rbind(res, out.mat)
                }else{
                    if(get.all.states == TRUE){
                        out.mat <- t(matrix(out, 2, 4))
                        colnames(out.mat) <- c("0", "1")
                        rownames(out.mat) <- c("A", "B", "C", "D")
                        return(out.mat / sum(out.mat))
                    }else{
                        out.mat <- t(matrix(out, 2, 4))
                        colnames(out.mat) <- c("0", "1")
                        res <- rbind(res, colSums(out.mat)/sum(out.mat))
                    }
                }
            }
        }
        AIC.vector <- sapply(hisse.results, "[[", "AIC")
        delta.AIC.vector <- AIC.vector - min(AIC.vector)
        rel.likelihood <- exp(-0.5 * delta.AIC.vector)
        AIC.weight.vector <- rel.likelihood / sum(rel.likelihood)
        final.eq.freq <- apply(res, 2, weighted.mean, w=AIC.weight.vector)
        if(get.rates == TRUE){
            return(final.eq.freq)
        }else{
            return(final.eq.freq /sum(final.eq.freq))
        }
    }
    
    if(model.type=="geohisse"){
        res <- c()
        geohisse.results <- x
        
        ## The object class can be a vector.
        if( inherits(x=geohisse.results, what="geohisse.fit") ) {
            ## we have to make a list so we can run this generally
            tmp.list <- list()
            tmp.list[[1]] <- geohisse.results
            geohisse.results <- tmp.list
        }
        for(model.index in 1:length(geohisse.results)){
            if(geohisse.results[[model.index]]$assume.cladogenetic == TRUE){
                ##Modify the data file
                data.new <- data.frame(geohisse.results[[model.index]]$data[,2], geohisse.results[[model.index]]$data[,2], row.names=geohisse.results[[model.index]]$data[,1])
                data.new <- data.new[geohisse.results[[model.index]]$phy$tip.label,]
                cache = ParametersToPassGeoHiSSE(geohisse.results[[model.index]]$phy, data.new[,1], model.vec=geohisse.results[[model.index]]$solution, f=geohisse.results[[model.index]]$f, hidden.states=TRUE)
                if( !geohisse.results[[model.index]]$root.type %in% c("madfitz","user","equal") ){
                    stop("Option for root.type is not implemented. Check help for GeoHiSSE.")
                }
                if(geohisse.results[[model.index]]$root.type=="madfitz"){
                    get.starting.probs <- DownPassGeoHisse(phy=geohisse.results[[model.index]]$phy, cache=cache, hidden.states=TRUE, condition.on.survival=geohisse.results[[model.index]]$condition.on.survival, root.type=geohisse.results[[model.index]]$root.type, root.p=geohisse.results[[model.index]]$root.p, get.phi=TRUE)$compD.root
                }else{
                    get.starting.probs <- geohisse.results[[model.index]]$root.p
                }
                out <- lsoda(c(state0A=get.starting.probs[1],state1A=get.starting.probs[2],state01A=get.starting.probs[3], state0B=get.starting.probs[4],state1B=get.starting.probs[5],state01B=get.starting.probs[6], state0C=get.starting.probs[7],state1C=get.starting.probs[8],state01C=get.starting.probs[9], state0D=get.starting.probs[10],state1D=get.starting.probs[11],state01D=get.starting.probs[12], state0E=get.starting.probs[13],state1E=get.starting.probs[14],state01E=get.starting.probs[15]), times=c(0, max.time), func=EqFreqsGeoHiSSE, parms=NULL, cache=cache, rtol=1e-8, atol=1e-8)[-1,-1]
            }else{
                data.new <- data.frame(geohisse.results[[model.index]]$data[,2], geohisse.results[[model.index]]$data[,2], row.names=geohisse.results[[model.index]]$data[,1])
                data.new <- data.new[geohisse.results[[model.index]]$phy$tip.label,]
                cache = ParametersToPassMuSSE(geohisse.results[[model.index]]$phy, data.new[,1], model.vec=geohisse.results[[model.index]]$solution, f=geohisse.results[[model.index]]$f, hidden.states=TRUE)
                if( !geohisse.results[[model.index]]$root.type %in% c("madfitz","user","equal") ){
                    stop("Option for root.type is not implemented. Check help for GeoHiSSE.")
                }
                if(geohisse.results[[model.index]]$root.type=="madfitz"){
                    get.starting.probs <- DownPassMusse(phy=geohisse.results[[model.index]]$phy, cache=cache, hidden.states=TRUE, condition.on.survival=geohisse.results[[model.index]]$condition.on.survival, root.type=geohisse.results[[model.index]]$root.type, root.p=geohisse.results[[model.index]]$root.p, get.phi=TRUE)$compD.root
                }else{
                    get.starting.probs <- geohisse.results[[model.index]]$root.p
                }
                out <- lsoda(c(state0A=get.starting.probs[1],state1A=get.starting.probs[2],state01A=get.starting.probs[3], state0B=get.starting.probs[4],state1B=get.starting.probs[5],state01B=get.starting.probs[6], state0C=get.starting.probs[7],state1C=get.starting.probs[8],state01C=get.starting.probs[9], state0D=get.starting.probs[10],state1D=get.starting.probs[11],state01D=get.starting.probs[12], state0E=get.starting.probs[13],state1E=get.starting.probs[14],state01E=get.starting.probs[15]), times=c(0, max.time), func=EqFreqsMuSSE, parms=NULL, cache=cache, rtol=1e-8, atol=1e-8)[-1,-1]
            }
            if(get.rates == TRUE){
                rescaled.probs.0 <- c(out[1],out[4],out[7],out[10],out[13]) / sum(out[1],out[4],out[7],out[10],out[13])
                rescaled.probs.1 <- c(out[2],out[5],out[8],out[11],out[14]) / sum(out[2],out[5],out[8],out[11],out[14])
                rescaled.probs.01 <- c(out[3],out[6],out[9],out[12],out[15]) / sum(out[3],out[6],out[9],out[12],out[15])
                if(rate.type == "turnover"){
                    state0 <- ((cache$s0A+cache$x0A) * rescaled.probs.0[1]) + ((cache$s0B+cache$x0B) * rescaled.probs.0[2]) + ((cache$s0C+cache$x0C) * rescaled.probs.0[3]) + ((cache$s0D+cache$x0D) * rescaled.probs.0[4]) + ((cache$s0E+cache$x0E) * rescaled.probs.0[5])
                    state1 <- ((cache$s1A+cache$x1A) * rescaled.probs.1[1]) + ((cache$s1B+cache$x1B) * rescaled.probs.1[2]) + ((cache$s1C+cache$x1C) * rescaled.probs.1[3]) + ((cache$s1D+cache$x1D) * rescaled.probs.1[4]) + ((cache$s1E+cache$x1E) * rescaled.probs.1[5])
                    state01 <- ((cache$s0A+cache$s1A+cache$s01A) * rescaled.probs.01[1]) + ((cache$s0B+cache$s1B+cache$s01B) * rescaled.probs.01[2]) + ((cache$s0C+cache$s1C+cache$s01C) * rescaled.probs.01[3]) + ((cache$s0D+cache$s1D+cache$s01D) * rescaled.probs.01[4]) + ((cache$s0E+cache$s1E+cache$s01E) * rescaled.probs.01[5])
                }
                if(rate.type == "net.div"){
                    state0 <- ((cache$s0A-cache$x0A) * rescaled.probs.0[1]) + ((cache$s0B-cache$x0B) * rescaled.probs.0[2]) + ((cache$s0C-cache$x0C) * rescaled.probs.0[3]) + ((cache$s0D-cache$x0D) * rescaled.probs.0[4]) + ((cache$s0E-cache$x0E) * rescaled.probs.0[5])
                    state1 <- ((cache$s1A-cache$x1A) * rescaled.probs.1[1]) + ((cache$s1B-cache$x1B) * rescaled.probs.1[2]) + ((cache$s1C-cache$x1C) * rescaled.probs.1[3]) + ((cache$s1D-cache$x1D) * rescaled.probs.1[4]) + ((cache$s1E-cache$x1E) * rescaled.probs.1[5])
                    state01 <- ((cache$s0A+cache$s1A+cache$s01A) * rescaled.probs.01[1]) + ((cache$s0B+cache$s1B+cache$s01B) * rescaled.probs.01[2]) + ((cache$s0C+cache$s1C+cache$s01C) * rescaled.probs.01[3]) + ((cache$s0D+cache$s1D+cache$s01D) * rescaled.probs.01[4]) + ((cache$s0E+cache$s1E+cache$s01E) * rescaled.probs.01[5])
                }
                if(rate.type == "speciation"){
                    state0 <- (cache$s0A * rescaled.probs.0[1]) + (cache$s0B * rescaled.probs.0[2]) + (cache$s0C * rescaled.probs.0[3]) + (cache$s0D * rescaled.probs.0[4]) + (cache$s0E * rescaled.probs.0[5])
                    state1 <- (cache$s1A * rescaled.probs.1[1]) + (cache$s1B * rescaled.probs.1[2]) + (cache$s1C * rescaled.probs.1[3]) + (cache$s1D * rescaled.probs.1[4]) + (cache$s1E * rescaled.probs.1[5])
                    state01 <- (cache$s01A * rescaled.probs.01[1]) + (cache$s01B * rescaled.probs.01[2]) + (cache$s01C * rescaled.probs.01[3]) + (cache$s01D * rescaled.probs.01[4]) + (cache$s01E * rescaled.probs.01[5])
                }
                out.mat <- matrix(c(state0, state1, state01), 1, 3)
                colnames(out.mat) <- c("0", "1", "01")
                res <- rbind(res, out.mat)
            }else{
                if(get.all.states == TRUE){
                    out.mat <- t(matrix(out, 3, 5))
                    colnames(out.mat) <- c("0", "1", "01")
                    rownames(out.mat) <- c("A", "B", "C", "D", "E")
                    return(out.mat / sum(out.mat))
                }else{
                    out.mat <- t(matrix(out, 3, 5))
                    colnames(out.mat) <- c("0", "1", "01")
                    res <- rbind(res, colSums(out.mat)/sum(out.mat))
                }
            }
        }
        AIC.vector <- sapply(geohisse.results, "[[", "AIC")
        delta.AIC.vector <- AIC.vector - min(AIC.vector)
        rel.likelihood <- exp(-0.5 * delta.AIC.vector)
        AIC.weight.vector <- rel.likelihood / sum(rel.likelihood)
        final.eq.freq <- apply(res, 2, weighted.mean, w=AIC.weight.vector)
        if(get.rates == TRUE){
            return(final.eq.freq)
        }else{
            return(final.eq.freq /sum(final.eq.freq))
        }
    }
}


EqFreqHiSSE <- function(t, y, parms, cache){
    dN0AdT = cache$lambda0 * y[1] - cache$death0 * y[1] - cache$q01 * y[1] - cache$q0A * y[1] - cache$q0B * y[1] + cache$q10 * y[2] + cache$qA0 * y[3] + cache$qB0 * y[4]
    dN1AdT = cache$lambda1 * y[2] - cache$death1 * y[2] - cache$q10 * y[2] - cache$q1A * y[2] - cache$q1B * y[2] + cache$q01 * y[1] + cache$qA1 * y[3] + cache$qB1 * y[4]
    dN0BdT = cache$lambdaA * y[3] - cache$deathA * y[3] - cache$qA0 * y[3] - cache$qA1 * y[3] - cache$qAB * y[3] + cache$q0A * y[1] + cache$q1A * y[2] + cache$qBA * y[4]
    dN1BdT = cache$lambdaB * y[4] - cache$deathB * y[4] - cache$qB0 * y[4] - cache$qB1 * y[4] - cache$qBA * y[4] + cache$q0B * y[1] + cache$q1B * y[2] + cache$qAB * y[3]
    
    return(list(c(dN0AdT,dN1AdT,dN0BdT,dN1BdT)))
}


EqFreqCID4 <- function(t, y, parms, cache){
    dN0AdT = cache$lambda0A * y[1] - cache$death0A * y[1] - cache$q0A1A * y[1] - cache$q0A0B * y[1] - cache$q0A0C * y[1] - cache$q0A0D * y[1] + cache$q1A0A * y[2] + cache$q0B0A * y[3] + cache$q0C0A * y[5] + cache$q0D0A * y[7]
    dN1AdT = cache$lambda1A * y[2] - cache$death1A * y[2] - cache$q1A0A * y[2] - cache$q1A1B * y[2] - cache$q1A1C * y[2] - cache$q1A1D * y[2] + cache$q0A1A * y[1] + cache$q1B1A * y[4]  + cache$q1C1A * y[6] + cache$q1D1A * y[8]

    dN0BdT = cache$lambda0B * y[3] - cache$death0B * y[3] - cache$q0B0A * y[3] - cache$q0B1B * y[3] - cache$q0B0C * y[3] - cache$q0B0D * y[3] + cache$q0A0B * y[1] + cache$q1B0B * y[4] + cache$q0C0B * y[5] + cache$q0D0B * y[7]
    dN1BdT = cache$lambda1B * y[4] - cache$death1B * y[4] - cache$q1B1A * y[4] - cache$q1B0B * y[4] - cache$q1B1C * y[4] - cache$q1B1D * y[4] + cache$q1A1B * y[2] + cache$q0B1B * y[3] + cache$q1C1B * y[6] + cache$q1D1B * y[8]

    dN0CdT = cache$lambda0C * y[5] - cache$death0C * y[5] - cache$q0C0A * y[5] - cache$q0C0B * y[5] - cache$q0C1C * y[5] - cache$q0C0D * y[5] + cache$q0A0C * y[1] + cache$q0B0C * y[3] + cache$q1C0C * y[6] + cache$q0D0C * y[7]
    dN1CdT = cache$lambda1C * y[6] - cache$death1C * y[6] - cache$q1C1A * y[6] - cache$q1C0C * y[6] - cache$q1C1B * y[6] - cache$q1C1D * y[6] + cache$q1A1C * y[2] + cache$q1B1C * y[4] + cache$q0C1C * y[5] + cache$q1D1C * y[8]
    
    dN0DdT = cache$lambda0D * y[7] - cache$death0D * y[7] - cache$q0D0A * y[7] - cache$q0D0B * y[7] - cache$q0D0C * y[7] - cache$q0D1D * y[7] + cache$q0A0D * y[1] + cache$q0B0D * y[3] + cache$q0C0D * y[5] + cache$q1D0D * y[8]
    dN1DdT = cache$lambda1D * y[8] - cache$death1D * y[8] - cache$q1D1A * y[8] - cache$q1D1B * y[8] - cache$q1D1C * y[8] - cache$q1D0D * y[8] + cache$q1A1D * y[2] + cache$q1B1D * y[4] + cache$q1C1D * y[6] + cache$q0D1D * y[7]
    
    return(list(c(dN0AdT,dN1AdT, dN0BdT,dN1BdT, dN0CdT,dN1CdT, dN0DdT,dN1DdT)))
}


#GeoHiSSE equilibrium freqs.
EqFreqsGeoHiSSE <- function(t, y, parms, cache){
    #cat A:
    dN0AdT  = cache$s0A * y[1] + cache$d01A_0A * y[3] + cache$s01A * y[3] + cache$s0A * y[3] - cache$x0A * y[1] - cache$d0A_01A * y[1] - cache$d0A_1A * y[1] - cache$d0A_0B * y[1] - cache$d0A_0C * y[1] - cache$d0A_0D * y[1] - cache$d0A_0E * y[1] + cache$d0B_0A * y[4] + cache$d0C_0A * y[7] + cache$d0D_0A * y[10] + cache$d0E_0A * y[13]
    dN1AdT  = cache$s1A * y[2] + cache$d01A_1A * y[3] + cache$s01A * y[3] + cache$s1A * y[3] - cache$x1A * y[2] - cache$d1A_01A * y[2] - cache$d1A_0A * y[2] - cache$d1A_1B * y[2] - cache$d1A_1C * y[2] - cache$d1A_1D * y[2] - cache$d1A_1E * y[2] + cache$d1B_1A * y[5] + cache$d1C_1A * y[8] + cache$d1D_1A * y[11] + cache$d1E_1A * y[14]
    dN01AdT = cache$d0A_01A * y[1] + cache$d1A_01A * y[2] - cache$s01A * y[3] - cache$d01A_1A * y[3] - cache$d01A_0A * y[3] - cache$d01A_01B * y[3] - cache$d01A_01C * y[3] - cache$d01A_01D * y[3] - cache$d01A_01E * y[3] + cache$d01B_01A * y[6] + cache$d01C_01A * y[9] + cache$d01D_01A * y[12] + cache$d01E_01A * y[15]
    #cat B:
    dN0BdT  = cache$s0B * y[4] + cache$d01B_0B * y[6] + cache$s01B * y[6] + cache$s0B * y[6] - cache$x0B * y[4] - cache$d0B_01B * y[4] - cache$d0B_1B * y[4] - cache$d0B_0A * y[4] - cache$d0A_0C * y[4] - cache$d0A_0D * y[4] - cache$d0A_0E * y[4] + cache$d0A_0B * y[1] + cache$d0C_0B * y[7] + cache$d0D_0B * y[10] + cache$d0E_0B * y[13]
    dN1BdT  = cache$s1B * y[5] + cache$d01B_1B * y[6] + cache$s01B * y[6] + cache$s1B * y[6] - cache$x1B * y[5] - cache$d1B_01B * y[5] - cache$d1B_0B * y[5] - cache$d1B_1A * y[5] - cache$d1B_1C * y[5] - cache$d1B_1D * y[5] - cache$d1B_1E * y[5] + cache$d1A_1B * y[2] + cache$d1C_1B * y[8] + cache$d1D_1B * y[11] + cache$d1E_1B * y[14]
    dN01BdT = cache$d0B_01B * y[4] + cache$d1B_01B * y[5] - cache$s01B * y[6] - cache$d01B_1B * y[6] - cache$d01B_0B * y[6] - cache$d01B_01A * y[6] - cache$d01B_01C * y[6] - cache$d01B_01D * y[6] - cache$d01B_01E * y[6] + cache$d01A_01B * y[3] + cache$d01C_01B * y[9] + cache$d01D_01B * y[12] + cache$d01E_01B * y[15]
    #cat C:
    dN0CdT  = cache$s0C * y[7] + cache$d01C_0C * y[9] + cache$s01C * y[9] + cache$s0C * y[9] - cache$x0C * y[7] - cache$d0C_01C * y[7] - cache$d0C_1C * y[7] - cache$d0C_0A * y[7] - cache$d0C_0B * y[7] - cache$d0C_0D * y[7] - cache$d0C_0E * y[7] + cache$d0A_0C * y[1] + cache$d0B_0C * y[4] + cache$d0D_0C * y[10] + cache$d0E_0C * y[13]
    dN1CdT  = cache$s1C * y[8] + cache$d01C_1C * y[9] + cache$s01C * y[9] + cache$s1C * y[9] - cache$x1C * y[8] - cache$d1C_01C * y[8] - cache$d1C_0C * y[8] - cache$d1C_1A * y[8] - cache$d1C_1B * y[8] - cache$d1C_1D * y[8] - cache$d1C_1E * y[8] + cache$d1A_1C * y[2] + cache$d1B_1C * y[5] + cache$d1D_1C * y[11] + cache$d1E_1C * y[14]
    dN01CdT = cache$d0C_01C * y[7] + cache$d1C_01C * y[8] - cache$s01C * y[9] - cache$d01C_1C * y[9] - cache$d01C_0C * y[9] - cache$d01C_01A * y[9] - cache$d01C_01B * y[9] - cache$d01C_01D * y[9] - cache$d01C_01E * y[9] + cache$d01A_01C * y[3] + cache$d01B_01C * y[6] + cache$d01D_01C * y[12] + cache$d01E_01C * y[15]
    #cat D:
    dN0DdT  = cache$s0D * y[10] + cache$d01D_0D * y[12] + cache$s01D * y[12] + cache$s0D * y[12] - cache$x0D * y[10] - cache$d0D_01D * y[10] - cache$d0D_1D * y[10] - cache$d0D_0A * y[10] - cache$d0D_0B * y[10] - cache$d0D_0C * y[10] - cache$d0D_0E * y[10] + cache$d0A_0D * y[1] + cache$d0B_0D * y[4] + cache$d0C_0D * y[7] + cache$d0E_0D * y[13]
    dN1DdT  = cache$s1D * y[11] + cache$d01D_1D * y[12] + cache$s01D * y[12] + cache$s1D * y[12] - cache$x1D * y[11] - cache$d1D_01D * y[11] - cache$d1D_0D * y[11] - cache$d1D_1A * y[11] - cache$d1D_1B * y[11] - cache$d1D_1C * y[11] - cache$d1D_1E * y[11] + cache$d1A_1D * y[2] + cache$d1B_1D * y[5] + cache$d1C_1D * y[8] + cache$d1E_1D * y[14]
    dN01DdT = cache$d0D_01D * y[10] + cache$d1D_01D * y[11] - cache$s01D * y[12] - cache$d01D_1D * y[12] - cache$d01D_0D * y[12] - cache$d01D_01A * y[12] - cache$d01D_01B * y[12] - cache$d01D_01C * y[12] - cache$d01D_01E * y[12] + cache$d01A_01D * y[3] + cache$d01B_01D * y[6] + cache$d01C_01D * y[9] + cache$d01E_01D * y[15]
    #cat E:
    dN0EdT  = cache$s0E * y[13] + cache$d01E_0E * y[15] + cache$s01E * y[15] + cache$s0E * y[15] - cache$x0E * y[13] - cache$d0E_01E * y[13] - cache$d0E_1E * y[13] - cache$d0E_0A * y[13] - cache$d0E_0B * y[13] - cache$d0E_0C * y[13] - cache$d0E_0D * y[13] + cache$d0A_0E * y[1] + cache$d0B_0E * y[4] + cache$d0C_0E * y[7] + cache$d0D_0E * y[10]
    dN1EdT  = cache$s1E * y[14] + cache$d01E_1E * y[15] + cache$s01E * y[15] + cache$s1E * y[15] - cache$x1E * y[14] - cache$d1E_01E * y[14] - cache$d1E_0E * y[14] - cache$d1E_1A * y[14] - cache$d1E_1B * y[14] - cache$d1E_1C * y[14] - cache$d1E_1D * y[14] + cache$d1A_1E * y[2] + cache$d1B_1E * y[5] + cache$d1C_1E * y[8] + cache$d1D_1E * y[11]
    dN01EdT = cache$d0E_01E * y[13] + cache$d1E_01E * y[14] - cache$s01E * y[15] - cache$d01E_1E * y[15] - cache$d01E_0E * y[15] - cache$d01E_01A * y[15] - cache$d01E_01B * y[15] - cache$d01E_01C * y[15] - cache$d01E_01D * y[15] + cache$d01A_01E * y[3] + cache$d01B_01E * y[6] + cache$d01C_01E * y[9] + cache$d01D_01E * y[12]
    
    return(list(c(dN0AdT,dN1AdT,dN01AdT, dN0BdT,dN1BdT,dN01BdT, dN0CdT,dN1CdT,dN01CdT, dN0DdT,dN1DdT,dN01DdT, dN0EdT,dN1EdT,dN01EdT)))
}


#MuSSE equlibrium freqs.
EqFreqsMuSSE <- function(t, y, parms, cache){
    #cat A:
    dN0AdT = cache$s0A * y[1] + cache$d01A_0A * y[3] + cache$d1A_0A * y[2] - cache$x0A * y[1] - cache$d0A_01A * y[1] - cache$d0A_1A * y[1] - cache$d0A_0B * y[1] - cache$d0A_0C * y[1] - cache$d0A_0D * y[1] - cache$d0A_0E * y[1] + cache$d0B_0A * y[4] + cache$d0C_0A * y[7] + cache$d0D_0A * y[10] + cache$d0E_0A * y[13]
    dN1AdT = cache$s1A * y[2] + cache$d01A_1A * y[3] + cache$d0A_1A * y[1] - cache$x1A * y[2] - cache$d1A_01A * y[2] - cache$d1A_0A * y[2] - cache$d1A_1B * y[2] - cache$d1A_1C * y[2] - cache$d1A_1D * y[2] - cache$d1A_1E * y[2] + cache$d1B_1A * y[5] + cache$d1C_1A * y[8] + cache$d1D_1A * y[11] + cache$d1E_1A * y[14]
    dN01AdT = cache$s01A * y[3] + cache$d0A_01A * y[1] + cache$d1A_01A * y[2] - cache$x01A * y[3] - cache$d01A_1A * y[3] - cache$d01A_0A * y[3] - cache$d01A_01B * y[3] - cache$d01A_01C * y[3] - cache$d01A_01D * y[3] - cache$d01A_01E * y[3] + cache$d01B_01A * y[6] + cache$d01C_01A * y[9] + cache$d01D_01A * y[12] + cache$d01E_01A * y[15]
    #cat B:
    dN0BdT = cache$s0B * y[4] + cache$d01B_0B * y[6] + cache$d1B_0B * y[5] - cache$x0B * y[4] - cache$d0B_01B * y[4] - cache$d0B_1B * y[4] - cache$d0B_0A * y[4] - cache$d0A_0C * y[4] - cache$d0A_0D * y[4] - cache$d0A_0E * y[4] + cache$d0A_0B * y[1] + cache$d0C_0B * y[7] + cache$d0D_0B * y[10] + cache$d0E_0B * y[13]
    dN1BdT = cache$s1B * y[5] + cache$d01B_1B * y[6] + cache$d0B_1B * y[4] - cache$x1B * y[5] - cache$d1B_01B * y[5] - cache$d1B_0B * y[5] - cache$d1B_1A * y[5] - cache$d1B_1C * y[5] - cache$d1B_1D * y[5] - cache$d1B_1E * y[5] + cache$d1A_1B * y[2] + cache$d1C_1B * y[8] + cache$d1D_1B * y[11] + cache$d1E_1B * y[14]
    dN01BdT = cache$s01B * y[6] + cache$d0B_01B * y[4] + cache$d1B_01B * y[5] - cache$x01B * y[6] - cache$d01B_1B * y[6] - cache$d01B_0B * y[6] - cache$d01B_01A * y[6] - cache$d01B_01C * y[6] - cache$d01B_01D * y[6] - cache$d01B_01E * y[6] + cache$d01A_01B * y[3] + cache$d01C_01B * y[9] + cache$d01D_01B * y[12] + cache$d01E_01B * y[15]
    #cat C:
    dN0CdT = cache$s0C * y[7] + cache$d01C_0C * y[9] + cache$d1C_0C * y[8] - cache$x0C * y[7] - cache$d0C_01C * y[7] - cache$d0C_1C * y[7] - cache$d0C_0A * y[7] - cache$d0C_0B * y[7] - cache$d0C_0D * y[7] - cache$d0C_0E * y[7] + cache$d0A_0C * y[1] + cache$d0B_0C * y[4] + cache$d0D_0C * y[10] + cache$d0E_0C * y[13]
    dN1CdT = cache$s1C * y[8] + cache$d01C_1C * y[9] + cache$d0C_1C * y[7] - cache$x1C * y[8] - cache$d1C_01C * y[8] - cache$d1C_0C * y[8] - cache$d1C_1A * y[8] - cache$d1C_1B * y[8] - cache$d1C_1D * y[8] - cache$d1C_1E * y[8] + cache$d1A_1C * y[2] + cache$d1B_1C * y[5] + cache$d1D_1C * y[11] + cache$d1E_1C * y[14]
    dN01CdT = cache$s01C * y[9] + cache$d0C_01C * y[7] + cache$d1C_01C * y[8] - cache$x01C * y[9] - cache$d01C_1C * y[9] - cache$d01C_0C * y[9] - cache$d01C_01A * y[9] - cache$d01C_01B * y[9] - cache$d01C_01D * y[9] - cache$d01C_01E * y[9] + cache$d01A_01C * y[3] + cache$d01B_01C * y[6] + cache$d01D_01C * y[12] + cache$d01E_01C * y[15]
    #cat D:
    dN0DdT = cache$s0D * y[10] + cache$d01D_0D * y[12] + cache$d1D_0D * y[11] - cache$x0D * y[10] - cache$d0D_01D * y[10] - cache$d0D_1D * y[10] - cache$d0D_0A * y[10] - cache$d0D_0B * y[10] - cache$d0D_0C * y[10] - cache$d0D_0E * y[10] + cache$d0A_0D * y[1] + cache$d0B_0D * y[4] + cache$d0C_0D * y[7] + cache$d0E_0D * y[13]
    dN1DdT = cache$s1D * y[11] + cache$d01D_1D * y[12] + cache$d0D_1D * y[10] - cache$x1D * y[11] - cache$d1D_01D * y[11] - cache$d1D_0D * y[11] - cache$d1D_1A * y[11] - cache$d1D_1B * y[11] - cache$d1D_1C * y[11] - cache$d1D_1E * y[11] + cache$d1A_1D * y[2] + cache$d1B_1D * y[5] + cache$d1C_1D * y[8] + cache$d1E_1D * y[14]
    dN01DdT = cache$s01D * y[12] + cache$d0D_01D * y[10] + cache$d1D_01D * y[11] - cache$x01D * y[12] - cache$d01D_1D * y[12] - cache$d01D_0D * y[12] - cache$d01D_01A * y[12] - cache$d01D_01B * y[12] - cache$d01D_01C * y[12] - cache$d01D_01E * y[12] + cache$d01A_01D * y[3] + cache$d01B_01D * y[6] + cache$d01C_01D * y[9] + cache$d01E_01D * y[15]
    #cat E:
    dN0EdT = cache$s0E * y[13] + cache$d01E_0E * y[15] + cache$d1E_0E * y[14] - cache$x0E * y[13] - cache$d0E_01E * y[13] - cache$d0E_1E * y[13] - cache$d0E_0A * y[13] - cache$d0E_0B * y[13] - cache$d0E_0C * y[13] - cache$d0E_0D * y[13] + cache$d0A_0E * y[1] + cache$d0B_0E * y[4] + cache$d0C_0E * y[7] + cache$d0D_0E * y[10]
    dN1EdT = cache$s1E * y[14] + cache$d01E_1E * y[15] + cache$d0E_1E * y[13] - cache$x1E * y[14] - cache$d1E_01E * y[14] - cache$d1E_0E * y[14] - cache$d1E_1A * y[14] - cache$d1E_1B * y[14] - cache$d1E_1C * y[14] - cache$d1E_1D * y[14] + cache$d1A_1E * y[2] + cache$d1B_1E * y[5] + cache$d1C_1E * y[8] + cache$d1D_1E * y[11]
    dN01EdT = cache$s01E * y[15] + cache$d0E_01E * y[13] + cache$d1E_01E * y[14] - cache$x01E * y[15] - cache$d01E_1E * y[15] - cache$d01E_0E * y[15] - cache$d01E_01A * y[15] - cache$d01E_01B * y[15] - cache$d01E_01C * y[15] - cache$d01E_01D * y[15] + cache$d01A_01E * y[3] + cache$d01B_01E * y[6] + cache$d01C_01E * y[9] + cache$d01D_01E * y[12]
    
    return(list(c(dN0AdT,dN1AdT,dN01AdT, dN0BdT,dN1BdT,dN01BdT, dN0CdT,dN1CdT,dN01CdT, dN0DdT,dN1DdT,dN01DdT, dN0EdT,dN1EdT,dN01EdT)))
}
