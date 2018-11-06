
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
                
                if(hisse.results[[model.index]]$root.type=="madfitz" | hisse.results[[model.index]]$root.type=="herr_als"){
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
                
                if(hisse.results[[model.index]]$root.type=="madfitz" | hisse.results[[model.index]]$root.type=="herr_als"){
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
                if( !geohisse.results[[model.index]]$root.type %in% c("madfitz","herr_als","user","equal") ){
                    stop("Option for root.type is not implemented. Check help for GeoHiSSE.")
                }
                if(geohisse.results[[model.index]]$root.type=="madfitz" | geohisse.results[[model.index]]$root.type=="herr_als"){
                    get.starting.probs <- DownPassGeoHisse(phy=geohisse.results[[model.index]]$phy, cache=cache, hidden.states=TRUE, condition.on.survival=geohisse.results[[model.index]]$condition.on.survival, root.type=geohisse.results[[model.index]]$root.type, root.p=geohisse.results[[model.index]]$root.p, get.phi=TRUE)$compD.root
                }else{
                    get.starting.probs <- geohisse.results[[model.index]]$root.p
                }
                out <- lsoda(c(state0A=get.starting.probs[1],state1A=get.starting.probs[2],state01A=get.starting.probs[3], state0B=get.starting.probs[4],state1B=get.starting.probs[5],state01B=get.starting.probs[6], state0C=get.starting.probs[7],state1C=get.starting.probs[8],state01C=get.starting.probs[9], state0D=get.starting.probs[10],state1D=get.starting.probs[11],state01D=get.starting.probs[12], state0E=get.starting.probs[13],state1E=get.starting.probs[14],state01E=get.starting.probs[15]), times=c(0, max.time), func=EqFreqsGeoHiSSE, parms=NULL, cache=cache, rtol=1e-8, atol=1e-8)[-1,-1]
            }else{
                data.new <- data.frame(geohisse.results[[model.index]]$data[,2], geohisse.results[[model.index]]$data[,2], row.names=geohisse.results[[model.index]]$data[,1])
                data.new <- data.new[geohisse.results[[model.index]]$phy$tip.label,]
                cache = ParametersToPassMuSSE(geohisse.results[[model.index]]$phy, data.new[,1], model.vec=geohisse.results[[model.index]]$solution, f=geohisse.results[[model.index]]$f, hidden.states=TRUE)
                if( !geohisse.results[[model.index]]$root.type %in% c("madfitz","herr_als","user","equal") ){
                    stop("Option for root.type is not implemented. Check help for GeoHiSSE.")
                }
                if(geohisse.results[[model.index]]$root.type=="madfitz" | geohisse.results[[model.index]]$root.type=="herr_als"){
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
    
    if(model.type=="muhisse"){
        res <- c()
        muhisse.results <- x
        
        ## The object class can be a vector.
        if( inherits(x=muhisse.results, what="muhisse.fit") ) {
            ## we have to make a list so we can run this generally
            tmp.list <- list()
            tmp.list[[1]] <- muhisse.results
            muhisse.results <- tmp.list
        }
        for(model.index in 1:length(muhisse.results)){
            data.new <- data.frame(muhisse.results[[model.index]]$data[,2], muhisse.results[[model.index]]$data[,3], row.names=muhisse.results[[model.index]]$data[,1])
            data.new <- data.new[muhisse.results[[model.index]]$phy$tip.label,]
            # Some new prerequisites #
            gen <- FindGenerations(muhisse.results[[model.index]]$phy)
            dat.tab <- OrganizeData(data=data.new, phy=muhisse.results[[model.index]]$phy, f=muhisse.results[[model.index]]$f, hidden.states=TRUE)
            nb.tip <- Ntip(muhisse.results[[model.index]]$phy)
            nb.node <- muhisse.results[[model.index]]$phy$Nnode
            ##########################
            cache = ParametersToPassMuHiSSE(model.vec=muhisse.results[[model.index]]$solution, hidden.states=TRUE, nb.tip=nb.tip, nb.node=nb.node, bad.likelihood=exp(-500), ode.eps=muhisse.results[[model.index]]$ode.eps)
            if( !muhisse.results[[model.index]]$root.type %in% c("madfitz","herr_als","user","equal") ){
                stop("Option for root.type is not implemented. Check help for MuHiSSE.")
            }
            if(muhisse.results[[model.index]]$root.type=="madfitz" | muhisse.results[[model.index]]$root.type=="herr_als"){
                get.starting.probs <- DownPassMuHisse(dat.tab=dat.tab, cache=cache, gen=gen, condition.on.survival=muhisse.results[[model.index]]$condition.on.survival, root.type=muhisse.results[[model.index]]$root.type, root.p=muhisse.results[[model.index]]$root.p, get.phi=TRUE)$compD.root
            }else{
                get.starting.probs <- muhisse.results[[model.index]]$root.p
            }
            out <- lsoda(c(state00A=get.starting.probs[1],state01A=get.starting.probs[2],state10A=get.starting.probs[3],state11A=get.starting.probs[4], state00B=get.starting.probs[5],state01B=get.starting.probs[6],state10B=get.starting.probs[7],state11B=get.starting.probs[8], state00C=get.starting.probs[9],state01C=get.starting.probs[10],state10C=get.starting.probs[11],state11C=get.starting.probs[12], state00D=get.starting.probs[13],state01D=get.starting.probs[14],state10D=get.starting.probs[15],state11D=get.starting.probs[16], state00E=get.starting.probs[17],state01E=get.starting.probs[18],state10E=get.starting.probs[19],state11E=get.starting.probs[20], state00F=get.starting.probs[21],state01F=get.starting.probs[22],state10F=get.starting.probs[23],state11F=get.starting.probs[24], state00G=get.starting.probs[25],state01G=get.starting.probs[26],state10G=get.starting.probs[27],state11G=get.starting.probs[28], state00H=get.starting.probs[29],state01H=get.starting.probs[30],state10H=get.starting.probs[31],state11H=get.starting.probs[32]), times=c(0, max.time), func=EqFreqsMuHiSSE, parms=NULL, cache=cache, rtol=1e-8, atol=1e-8)[-1,-1]
            if(get.rates == TRUE){
                rescaled.probs.00 <- c(out[1],out[5],out[9],out[13],out[17],out[21],out[25],out[29]) / sum(out[1],out[5],out[9],out[13],out[17],out[21],out[25],out[29])
                rescaled.probs.01 <- c(out[2],out[6],out[10],out[14],out[18],out[22],out[26],out[30]) / sum(out[2],out[6],out[10],out[14],out[18],out[22],out[26],out[30])
                rescaled.probs.10 <- c(out[3],out[7],out[11],out[15],out[19],out[23],out[27],out[31]) / sum(out[3],out[7],out[11],out[15],out[19],out[23],out[27],out[31])
                rescaled.probs.11 <- c(out[4],out[8],out[12],out[16],out[20],out[24],out[28],out[32]) / sum(out[4],out[8],out[12],out[16],out[20],out[24],out[28],out[32])

                if(rate.type == "turnover"){
                    state00 <- ((cache$lambda00A+cache$mu00A) * rescaled.probs.00[1]) + ((cache$lambda00B+cache$mu00B) * rescaled.probs.00[2]) + ((cache$lambda00C+cache$mu00C) * rescaled.probs.00[3]) + ((cache$lambda00D+cache$mu00D) * rescaled.probs.00[4]) + ((cache$lambda00E+cache$mu00E) * rescaled.probs.00[5]) + ((cache$lambda00F+cache$mu00F) * rescaled.probs.00[6]) + ((cache$lambda00G+cache$mu00G) * rescaled.probs.00[7]) + ((cache$lambda00H+cache$mu00H) * rescaled.probs.00[8])
                    state01 <- ((cache$lambda01A+cache$mu01A) * rescaled.probs.01[1]) + ((cache$lambda01B+cache$mu01B) * rescaled.probs.01[2]) + ((cache$lambda01C+cache$mu01C) * rescaled.probs.01[3]) + ((cache$lambda01D+cache$mu01D) * rescaled.probs.01[4]) + ((cache$lambda01E+cache$mu01E) * rescaled.probs.01[5]) + ((cache$lambda01F+cache$mu01F) * rescaled.probs.01[6]) + ((cache$lambda01G+cache$mu01G) * rescaled.probs.01[7]) + ((cache$lambda01H+cache$mu01H) * rescaled.probs.01[8])
                    state10 <- ((cache$lambda10A+cache$mu10A) * rescaled.probs.10[1]) + ((cache$lambda10B+cache$mu10B) * rescaled.probs.10[2]) + ((cache$lambda10C+cache$mu10C) * rescaled.probs.10[3]) + ((cache$lambda10D+cache$mu10D) * rescaled.probs.10[4]) + ((cache$lambda10E+cache$mu10E) * rescaled.probs.10[5]) + ((cache$lambda10F+cache$mu10F) * rescaled.probs.10[6]) + ((cache$lambda10G+cache$mu10G) * rescaled.probs.10[7]) + ((cache$lambda10H+cache$mu10H) * rescaled.probs.10[8])
                    state11 <- ((cache$lambda11A+cache$mu11A) * rescaled.probs.11[1]) + ((cache$lambda11B+cache$mu11B) * rescaled.probs.11[2]) + ((cache$lambda11C+cache$mu11C) * rescaled.probs.11[3]) + ((cache$lambda11D+cache$mu11D) * rescaled.probs.11[4]) + ((cache$lambda11E+cache$mu11E) * rescaled.probs.11[5]) + ((cache$lambda11F+cache$mu11F) * rescaled.probs.11[6]) + ((cache$lambda11G+cache$mu11G) * rescaled.probs.11[7]) + ((cache$lambda11H+cache$mu11H) * rescaled.probs.11[8])
                }
                if(rate.type == "net.div"){
                    state00 <- ((cache$lambda00A-cache$mu00A) * rescaled.probs.00[1]) + ((cache$lambda00B-cache$mu00B) * rescaled.probs.00[2]) + ((cache$lambda00C-cache$mu00C) * rescaled.probs.00[3]) + ((cache$lambda00D-cache$mu00D) * rescaled.probs.00[4]) + ((cache$lambda00E-cache$mu00E) * rescaled.probs.00[5]) + ((cache$lambda00F-cache$mu00F) * rescaled.probs.00[6]) + ((cache$lambda00G-cache$mu00G) * rescaled.probs.00[7]) + ((cache$lambda00H-cache$mu00H) * rescaled.probs.00[8])
                    state01 <- ((cache$lambda01A-cache$mu01A) * rescaled.probs.01[1]) + ((cache$lambda01B-cache$mu01B) * rescaled.probs.01[2]) + ((cache$lambda01C-cache$mu01C) * rescaled.probs.01[3]) + ((cache$lambda01D-cache$mu01D) * rescaled.probs.01[4]) + ((cache$lambda01E-cache$mu01E) * rescaled.probs.01[5]) + ((cache$lambda01F-cache$mu01F) * rescaled.probs.01[6]) + ((cache$lambda01G-cache$mu01G) * rescaled.probs.01[7]) + ((cache$lambda01H-cache$mu01H) * rescaled.probs.01[8])
                    state10 <- ((cache$lambda10A-cache$mu10A) * rescaled.probs.10[1]) + ((cache$lambda10B-cache$mu10B) * rescaled.probs.10[2]) + ((cache$lambda10C-cache$mu10C) * rescaled.probs.10[3]) + ((cache$lambda10D-cache$mu10D) * rescaled.probs.10[4]) + ((cache$lambda10E-cache$mu10E) * rescaled.probs.10[5]) + ((cache$lambda10F-cache$mu10F) * rescaled.probs.10[6]) + ((cache$lambda10G-cache$mu10G) * rescaled.probs.10[7]) + ((cache$lambda10H-cache$mu10H) * rescaled.probs.10[8])
                    state11 <- ((cache$lambda11A-cache$mu11A) * rescaled.probs.11[1]) + ((cache$lambda11B-cache$mu11B) * rescaled.probs.11[2]) + ((cache$lambda11C-cache$mu11C) * rescaled.probs.11[3]) + ((cache$lambda11D-cache$mu11D) * rescaled.probs.11[4]) + ((cache$lambda11E-cache$mu11E) * rescaled.probs.11[5]) + ((cache$lambda11F-cache$mu11F) * rescaled.probs.11[6]) + ((cache$lambda11G-cache$mu11G) * rescaled.probs.11[7]) + ((cache$lambda11H-cache$mu11H) * rescaled.probs.11[8])
                }
                if(rate.type == "speciation"){
                    state00 <- (cache$lambda00A * rescaled.probs.00[1]) + (cache$lambda00B * rescaled.probs.00[2]) + (cache$lambda00C * rescaled.probs.00[3]) + (cache$lambda00D * rescaled.probs.00[4]) + (cache$lambda00E * rescaled.probs.00[5]) + (cache$lambda00F * rescaled.probs.00[6]) + (cache$lambda00G * rescaled.probs.00[7]) + (cache$lambda00H * rescaled.probs.00[8])
                    state01 <- (cache$lambda01A * rescaled.probs.01[1]) + (cache$lambda01B * rescaled.probs.01[2]) + (cache$lambda01C * rescaled.probs.01[3]) + (cache$lambda01D * rescaled.probs.01[4]) + (cache$lambda01E * rescaled.probs.01[5]) + (cache$lambda01F * rescaled.probs.01[6]) + (cache$lambda01G * rescaled.probs.01[7]) + (cache$lambda01H * rescaled.probs.01[8])
                    state10 <- (cache$lambda10A * rescaled.probs.10[1]) + (cache$lambda10B * rescaled.probs.10[2]) + (cache$lambda10C * rescaled.probs.10[3]) + (cache$lambda10D * rescaled.probs.10[4]) + (cache$lambda10E * rescaled.probs.10[5]) + (cache$lambda10F * rescaled.probs.10[6]) + (cache$lambda10G * rescaled.probs.10[7]) + (cache$lambda10H * rescaled.probs.10[8])
                    state11 <- (cache$lambda11A * rescaled.probs.11[1]) + (cache$lambda11B * rescaled.probs.11[2]) + (cache$lambda11C * rescaled.probs.11[3]) + (cache$lambda11D * rescaled.probs.11[4]) + (cache$lambda11E * rescaled.probs.11[5]) + (cache$lambda11F * rescaled.probs.11[6]) + (cache$lambda11G * rescaled.probs.11[7]) + (cache$lambda11H * rescaled.probs.11[8])
                }
                out.mat <- matrix(c(state00, state01, state10, state11), 1, 4)
                colnames(out.mat) <- c("00", "01", "10", "11")
                res <- rbind(res, out.mat)
            }else{
                if(get.all.states == TRUE){
                    out.mat <- t(matrix(out, 4, 8))
                    colnames(out.mat) <- c("00", "01", "10", "11")
                    rownames(out.mat) <- c("A", "B", "C", "D", "E", "F", "G", "H")
                    return(out.mat / sum(out.mat))
                }else{
                    out.mat <- t(matrix(out, 4, 8))
                    colnames(out.mat) <- c("00", "01", "10", "11")
                    res <- rbind(res, colSums(out.mat)/sum(out.mat))
                }
            }
        }
        AIC.vector <- sapply(muhisse.results, "[[", "AIC")
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



#MuHiSSE equlibrium freqs.
EqFreqsMuHiSSE <- function(t, y, parms, cache){
    
    #cat A:
    dN00AdT <- cache$lambda00A * y[1] - cache$mu00A * y[1] - cache$q00A_01A * y[1] - cache$q00A_10A * y[1] - cache$q00A_11A * y[1] + cache$q01A_00A * y[2] + cache$q10A_00A * y[3] + cache$q11A_00A * y[4] - cache$q00A_00B * y[1] - cache$q00A_00C * y[1] - cache$q00A_00D * y[1] - cache$q00A_00E * y[1] - cache$q00A_00F * y[1] - cache$q00A_00G * y[1] - cache$q00A_00H * y[1] + cache$q00B_00A * y[5] + cache$q00C_00A * y[9] + cache$q00D_00A * y[13] + cache$q00E_00A * y[17] + cache$q00F_00A * y[21] + cache$q00G_00A * y[25] + cache$q00H_00A * y[29]
    dN01AdT <- cache$lambda01A * y[2] - cache$mu01A * y[2] - cache$q01A_00A * y[2] - cache$q01A_10A * y[2] - cache$q01A_11A * y[2] + cache$q00A_01A * y[1] + cache$q10A_01A * y[3] + cache$q11A_01A * y[4] - cache$q01A_01B * y[2] - cache$q01A_01C * y[2] - cache$q01A_01D * y[2] - cache$q01A_01E * y[2] - cache$q01A_01F * y[2] - cache$q01A_01G * y[2] - cache$q01A_01H * y[2] + cache$q01B_01A * y[6] + cache$q01C_01A * y[10] + cache$q01D_01A * y[14] + cache$q01E_01A * y[18] + cache$q01F_01A * y[22] + cache$q01G_01A * y[26] + cache$q01H_01A * y[30]
    dN10AdT <- cache$lambda10A * y[3] - cache$mu10A * y[3] - cache$q10A_00A * y[3] - cache$q10A_01A * y[3] - cache$q10A_11A * y[3] + cache$q00A_10A * y[1] + cache$q01A_10A * y[2] + cache$q11A_10A * y[4] - cache$q10A_10B * y[3] - cache$q10A_10C * y[3] - cache$q10A_10D * y[3] - cache$q10A_10E * y[3] - cache$q10A_10F * y[3] - cache$q10A_10G * y[3] - cache$q10A_10H * y[3] + cache$q10B_10A * y[7] + cache$q10C_10A * y[11] + cache$q10D_10A * y[15] + cache$q10E_10A * y[19] + cache$q10F_10A * y[23] + cache$q10G_10A * y[27] + cache$q10H_10A * y[31]
    dN11AdT <- cache$lambda11A * y[4] - cache$mu11A * y[4] - cache$q11A_00A * y[4] - cache$q11A_01A * y[4] - cache$q11A_10A * y[4] + cache$q00A_11A * y[1] + cache$q01A_11A * y[2] + cache$q10A_11A * y[3] - cache$q11A_11B * y[4] - cache$q11A_11C * y[4] - cache$q11A_11D * y[4] - cache$q11A_11E * y[4] - cache$q11A_11F * y[4] - cache$q11A_11G * y[4] - cache$q11A_11H * y[4] + cache$q11B_11A * y[8] + cache$q11C_11A * y[12] + cache$q11D_11A * y[16] + cache$q11E_11A * y[20] + cache$q11F_11A * y[24] + cache$q11G_11A * y[28] + cache$q11H_11A * y[32]
    
    #cat B:
    dN00BdT <- cache$lambda00B * y[5] - cache$mu00B * y[5] - cache$q00B_01B * y[5] - cache$q00B_10B * y[5] - cache$q00B_11B * y[5] + cache$q01B_00B * y[6] + cache$q10B_00B * y[7] + cache$q11B_00B * y[8] - cache$q00B_00A * y[5] - cache$q00B_00C * y[5] - cache$q00B_00D * y[5] - cache$q00B_00E * y[5] - cache$q00B_00F * y[5] - cache$q00B_00G * y[5] - cache$q00B_00H * y[5] + cache$q00A_00B * y[1] + cache$q00C_00B * y[9] + cache$q00D_00B * y[13] + cache$q00E_00B * y[17] + cache$q00F_00B * y[21] + cache$q00G_00B * y[25] + cache$q00H_00B * y[29]
    dN01BdT <- cache$lambda01B * y[6] - cache$mu01B * y[6] - cache$q01B_00B * y[6] - cache$q01B_10B * y[6] - cache$q01B_11B * y[6] + cache$q00B_01B * y[5] + cache$q10B_01B * y[7] + cache$q11B_01B * y[8] - cache$q01B_01A * y[6] - cache$q01B_01C * y[6] - cache$q01B_01D * y[6] - cache$q01B_01E * y[6] - cache$q01B_01F * y[6] - cache$q01B_01G * y[6] - cache$q01B_01H * y[6] + cache$q01A_01B * y[2] + cache$q01C_01B * y[10] + cache$q01D_01B * y[14] + cache$q01E_01B * y[18] + cache$q01F_01B * y[22] + cache$q01G_01B * y[26] + cache$q01H_01B * y[30]
    dN10BdT <- cache$lambda10B * y[7] - cache$mu10B * y[7] - cache$q10B_00B * y[7] - cache$q10B_01B * y[7] - cache$q10B_11B * y[7] + cache$q00B_10B * y[5] + cache$q01B_10B * y[6] + cache$q11B_10B * y[8] - cache$q10B_10A * y[7] - cache$q10B_10C * y[7] - cache$q10B_10D * y[7] - cache$q10B_10E * y[7] - cache$q10B_10F * y[7] - cache$q10B_10G * y[7] - cache$q10B_10H * y[7] + cache$q10A_10B * y[3] + cache$q10C_10B * y[11] + cache$q10D_10B * y[15] + cache$q10E_10B * y[19] + cache$q10F_10B * y[23] + cache$q10G_10B * y[27] + cache$q10H_10B * y[31]
    dN11BdT <- cache$lambda11B * y[8] - cache$mu11B * y[8] - cache$q11B_00B * y[8] - cache$q11B_01B * y[8] - cache$q11B_10B * y[8] + cache$q00B_11B * y[5] + cache$q01B_11B * y[6] + cache$q10B_11B * y[7] - cache$q11B_11A * y[8] - cache$q11B_11C * y[8] - cache$q11B_11D * y[8] - cache$q11B_11E * y[8] - cache$q11B_11F * y[8] - cache$q11B_11G * y[8] - cache$q11B_11H * y[8] + cache$q11A_11B * y[4] + cache$q11C_11B * y[12] + cache$q11D_11B * y[16] + cache$q11E_11B * y[20] + cache$q11F_11B * y[24] + cache$q11G_11B * y[28] + cache$q11H_11B * y[32]
    
    #cat C:
    dN00CdT <- cache$lambda00C * y[9] - cache$mu00C * y[9] - cache$q00C_01C * y[9] - cache$q00C_10C * y[9] - cache$q00C_11C * y[9] + cache$q01C_00C * y[10] + cache$q10C_00C * y[11] + cache$q11C_00C * y[12] - cache$q00C_00A * y[9] - cache$q00C_00B * y[9] - cache$q00C_00D * y[9] - cache$q00C_00E * y[9] - cache$q00C_00F * y[9] - cache$q00C_00G * y[9] - cache$q00C_00H * y[9] + cache$q00A_00C * y[1] + cache$q00B_00C * y[5] + cache$q00D_00C * y[13] + cache$q00E_00C * y[17] + cache$q00F_00C * y[21] + cache$q00G_00C * y[25] + cache$q00H_00C * y[29]
    dN01CdT <- cache$lambda01C * y[10] - cache$mu01C * y[10] - cache$q01C_00C * y[10] - cache$q01C_10C * y[10] - cache$q01C_11C * y[10] + cache$q00C_01C * y[9] + cache$q10C_01C * y[11] + cache$q11C_01C * y[12] - cache$q01C_01A * y[10] - cache$q01C_01B * y[10] - cache$q01C_01D * y[10] - cache$q01C_01E * y[10] - cache$q01C_01F * y[10] - cache$q01C_01G * y[10] - cache$q01C_01H * y[10] + cache$q01A_01C * y[2] + cache$q01B_01C * y[6] + cache$q01D_01C * y[14] + cache$q01E_01C * y[18] + cache$q01F_01C * y[22] + cache$q01G_01C * y[26] + cache$q01H_01C * y[30]
    dN10CdT <- cache$lambda10C * y[11] - cache$mu10C * y[11] - cache$q10C_00C * y[11] - cache$q10C_01C * y[11] - cache$q10C_11C * y[11] + cache$q00C_10C * y[9] + cache$q01C_10C * y[10] + cache$q11C_10C * y[12] - cache$q10C_10A * y[11] - cache$q10C_10B * y[11] - cache$q10C_10D * y[11] - cache$q10C_10E * y[11] - cache$q10C_10F * y[11] - cache$q10C_10G * y[11] - cache$q10C_10H * y[11] + cache$q10A_10C * y[3] + cache$q10B_10C * y[7] + cache$q10D_10C * y[15] + cache$q10E_10C * y[19] + cache$q10F_10C * y[23] + cache$q10G_10C * y[27] + cache$q10H_10C * y[31]
    dN11CdT <- cache$lambda11C * y[12] - cache$mu11C * y[12] - cache$q11C_00C * y[12] - cache$q11C_01C * y[12] - cache$q11C_10C * y[12] + cache$q00C_11C * y[9] + cache$q01C_11C * y[10] + cache$q10C_11C * y[11] - cache$q11C_11A * y[12] - cache$q11C_11B * y[12] - cache$q11C_11D * y[12] - cache$q11C_11E * y[12] - cache$q11C_11F * y[12] - cache$q11C_11G * y[12] - cache$q11C_11H * y[12] + cache$q11A_11C * y[4] + cache$q11B_11C * y[8] + cache$q11D_11C * y[16] + cache$q11E_11C * y[20] + cache$q11F_11C * y[24] + cache$q11G_11C * y[28] + cache$q11H_11C * y[32]
    
    #cat D:
    dN00DdT <- cache$lambda00D * y[13] - cache$mu00D * y[13] - cache$q00D_01D * y[13] - cache$q00D_10D * y[13] - cache$q00D_11D * y[13] + cache$q01D_00D * y[14] + cache$q10D_00D * y[15] + cache$q11D_00D * y[16] - cache$q00D_00A * y[13] - cache$q00D_00B * y[13] - cache$q00D_00C * y[13] - cache$q00D_00E * y[13] - cache$q00D_00F * y[13] - cache$q00D_00G * y[13] - cache$q00D_00H * y[13] + cache$q00A_00D * y[1] + cache$q00B_00D * y[5] + cache$q00C_00D * y[9] + cache$q00E_00D * y[17] + cache$q00F_00D * y[21] + cache$q00G_00D * y[25] + cache$q00H_00D * y[29]
    dN01DdT <- cache$lambda01D * y[14] - cache$mu01D * y[14] - cache$q01D_00D * y[14] - cache$q01D_10D * y[14] - cache$q01D_11D * y[14] + cache$q00D_01D * y[13] + cache$q10D_01D * y[15] + cache$q11D_01D * y[16] - cache$q01D_01A * y[14] - cache$q01D_01B * y[14] - cache$q01D_01C * y[14] - cache$q01D_01E * y[14] - cache$q01D_01F * y[14] - cache$q01D_01G * y[14] - cache$q01D_01H * y[14] + cache$q01A_01D * y[2] + cache$q01B_01D * y[6] + cache$q01C_01D * y[10] + cache$q01E_01D * y[18] + cache$q01F_01D * y[22] + cache$q01G_01D * y[26] + cache$q01H_01D * y[30]
    dN10DdT <- cache$lambda10D * y[15] - cache$mu10D * y[15] - cache$q10D_00D * y[15] - cache$q10D_01D * y[15] - cache$q10D_11D * y[15] + cache$q00D_10D * y[13] + cache$q01D_10D * y[14] + cache$q11D_10D * y[16] - cache$q10D_10A * y[15] - cache$q10D_10B * y[15] - cache$q10D_10C * y[15] - cache$q10D_10E * y[15] - cache$q10D_10F * y[15] - cache$q10D_10G * y[15] - cache$q10D_10H * y[15] + cache$q10A_10D * y[3] + cache$q10B_10D * y[7] + cache$q10C_10D * y[11] + cache$q10E_10D * y[19] + cache$q10F_10D * y[23] + cache$q10G_10D * y[27] + cache$q10H_10D * y[31]
    dN11DdT <- cache$lambda11D * y[16] - cache$mu11D * y[16] - cache$q11D_00D * y[16] - cache$q11D_01D * y[16] - cache$q11D_10D * y[16] + cache$q00D_11D * y[13] + cache$q01D_11D * y[14] + cache$q10D_11D * y[15] - cache$q11D_11A * y[16] - cache$q11D_11B * y[16] - cache$q11D_11C * y[16] - cache$q11D_11E * y[16] - cache$q11D_11F * y[16] - cache$q11D_11G * y[16] - cache$q11D_11H * y[16] + cache$q11A_11D * y[4] + cache$q11B_11D * y[8] + cache$q11C_11D * y[12] + cache$q11E_11D * y[20] + cache$q11F_11D * y[24] + cache$q11G_11D * y[28] + cache$q11H_11D * y[32]
    
    #cat E:
    dN00EdT <- cache$lambda00E * y[17] - cache$mu00E * y[17] - cache$q00E_01E * y[17] - cache$q00E_10E * y[17] - cache$q00E_11E * y[17] + cache$q01E_00E * y[18] + cache$q10E_00E * y[19] + cache$q11E_00E * y[20] - cache$q00E_00A * y[17] - cache$q00E_00B * y[17] - cache$q00E_00C * y[17] - cache$q00E_00D * y[17] - cache$q00E_00F * y[17] - cache$q00E_00G * y[17] - cache$q00E_00H * y[17] + cache$q00A_00E * y[1] + cache$q00B_00E * y[5] + cache$q00C_00E * y[9] + cache$q00D_00E * y[13] + cache$q00F_00E * y[21] + cache$q00G_00E * y[25] + cache$q00H_00E * y[29]
    dN01EdT <- cache$lambda01E * y[18] - cache$mu01E * y[18] - cache$q01E_00E * y[18] - cache$q01E_10E * y[18] - cache$q01E_11E * y[18] + cache$q00E_01E * y[17] + cache$q10E_01E * y[19] + cache$q11E_01E * y[20] - cache$q01E_01A * y[18] - cache$q01E_01B * y[18] - cache$q01E_01C * y[18] - cache$q01E_01D * y[18] - cache$q01E_01F * y[18] - cache$q01E_01G * y[18] - cache$q01E_01H * y[18] + cache$q01A_01E * y[2] + cache$q01B_01E * y[6] + cache$q01C_01E * y[10] + cache$q01D_01E * y[14] + cache$q01F_01E * y[22] + cache$q01G_01E * y[26] + cache$q01H_01E * y[30]
    dN10EdT <- cache$lambda10E * y[19] - cache$mu10E * y[19] - cache$q10E_00E * y[19] - cache$q10E_01E * y[19] - cache$q10E_11E * y[19] + cache$q00E_10E * y[17] + cache$q01E_10E * y[18] + cache$q11E_10E * y[20] - cache$q10E_10A * y[19] - cache$q10E_10B * y[19] - cache$q10E_10C * y[19] - cache$q10E_10D * y[19] - cache$q10E_10F * y[19] - cache$q10E_10G * y[19] - cache$q10E_10H * y[19] + cache$q10A_10E * y[3] + cache$q10B_10E * y[7] + cache$q10C_10E * y[11] + cache$q10D_10E * y[15] + cache$q10F_10E * y[23] + cache$q10G_10E * y[27] + cache$q10H_10E * y[31]
    dN11EdT <- cache$lambda11E * y[20] - cache$mu11E * y[20] - cache$q11E_00E * y[20] - cache$q11E_01E * y[20] - cache$q11E_10E* y[20] + cache$q00E_11E * y[17] + cache$q01E_11E * y[18] + cache$q10E_11E * y[19] - cache$q11E_11A * y[20] - cache$q11E_11B * y[20] - cache$q11E_11C * y[20] - cache$q11E_11D * y[20] - cache$q11E_11F * y[20] - cache$q11E_11G * y[20] - cache$q11E_11H * y[20] + cache$q11A_11E * y[4] + cache$q11B_11E * y[8] + cache$q11C_11E * y[12] + cache$q11D_11E * y[16] + cache$q11F_11E * y[24] + cache$q11G_11E * y[28] + cache$q11H_11E * y[32]

    #cat F:
    dN00FdT <- cache$lambda00F * y[21] - cache$mu00F * y[21] - cache$q00F_01F * y[21] - cache$q00F_10F * y[21] - cache$q00F_11F * y[21] + cache$q01F_00F * y[22] + cache$q10F_00F * y[23] + cache$q11F_00F * y[24] - cache$q00F_00A * y[21] - cache$q00F_00B * y[21] - cache$q00F_00C * y[21] - cache$q00F_00D * y[21] - cache$q00F_00E * y[21] - cache$q00F_00G * y[21] - cache$q00F_00H * y[21] + cache$q00A_00F * y[1] + cache$q00B_00F * y[5] + cache$q00C_00F * y[9] + cache$q00D_00F * y[13] + cache$q00E_00F * y[17] + cache$q00G_00F * y[25] + cache$q00H_00F * y[29]
    dN01FdT <- cache$lambda01F * y[22] - cache$mu01F * y[22] - cache$q01F_00F * y[22] - cache$q01F_10F * y[22] - cache$q01F_11F * y[22] + cache$q00F_01F * y[21] + cache$q10F_01F * y[23] + cache$q11F_01F * y[24] - cache$q01F_01A * y[22] - cache$q01F_01B * y[22] - cache$q01F_01C * y[22] - cache$q01F_01D * y[22] - cache$q01F_01E * y[22] - cache$q01F_01G * y[22] - cache$q01F_01H * y[22] + cache$q01A_01F * y[2] + cache$q01B_01F * y[6] + cache$q01C_01F * y[10] + cache$q01D_01F * y[14] + cache$q01E_01F * y[18] + cache$q01G_01F * y[26] + cache$q01H_01F * y[30]
    dN10FdT <- cache$lambda10F * y[23] - cache$mu10F * y[23] - cache$q10F_00F * y[23] - cache$q10F_01F * y[23] - cache$q10F_11F * y[23] + cache$q00F_10F * y[21] + cache$q01F_10F * y[22] + cache$q11F_10F * y[24] - cache$q10F_10A * y[23] - cache$q10F_10B * y[23] - cache$q10F_10C * y[23] - cache$q10F_10D * y[23] - cache$q10F_10E * y[23] - cache$q10F_10G * y[23] - cache$q10F_10H * y[23] + cache$q10A_10F * y[3] + cache$q10B_10F * y[7] + cache$q10C_10F * y[11] + cache$q10D_10F * y[15] + cache$q10E_10F * y[19] + cache$q10G_10F * y[27] + cache$q10H_10F * y[31]
    dN11FdT <- cache$lambda11F * y[24] - cache$mu11F * y[24] - cache$q11F_00F * y[24] - cache$q11F_01F * y[24] - cache$q11F_10F * y[24] + cache$q00F_11F * y[21] + cache$q01F_11F * y[22] + cache$q10F_11F * y[23] - cache$q11F_11A * y[24] - cache$q11F_11B * y[24] - cache$q11F_11C * y[24] - cache$q11F_11D * y[24] - cache$q11F_11E * y[24] - cache$q11F_11G * y[24] - cache$q11F_11H * y[24] + cache$q11A_11F * y[4] + cache$q11B_11F * y[8] + cache$q11C_11F * y[12] + cache$q11D_11F * y[16] + cache$q11E_11F * y[20] + cache$q11G_11F * y[28] + cache$q11H_11F * y[32]
    
    #cat G:
    dN00GdT <- cache$lambda00G * y[25] - cache$mu00G * y[25] - cache$q00G_01G * y[25] - cache$q00G_10G * y[25] - cache$q00G_11G * y[25] + cache$q01G_00G * y[26] + cache$q10G_00G * y[27] + cache$q11G_00G * y[28] - cache$q00G_00A * y[25] - cache$q00G_00B * y[25] - cache$q00G_00C * y[25] - cache$q00G_00D * y[25] - cache$q00G_00E * y[25] - cache$q00G_00F * y[25] - cache$q00G_00H * y[25] + cache$q00A_00G * y[1] + cache$q00B_00G * y[5] + cache$q00C_00G * y[9] + cache$q00D_00G * y[13] + cache$q00E_00G * y[17] + cache$q00F_00G * y[21] + cache$q00H_00G * y[29]
    dN01GdT <- cache$lambda01G * y[26] - cache$mu01G * y[26] - cache$q01G_00G * y[26] - cache$q01G_10G * y[26] - cache$q01G_11G * y[26] + cache$q00G_01G * y[25] + cache$q10G_01G * y[27] + cache$q11G_01G * y[28] - cache$q01G_01A * y[26] - cache$q01G_01B * y[26] - cache$q01G_01C * y[26] - cache$q01G_01D * y[26] - cache$q01G_01E * y[26] - cache$q01G_01F * y[26] - cache$q01G_01H * y[26] + cache$q01A_01G * y[2] + cache$q01B_01G * y[6] + cache$q01C_01G * y[10] + cache$q01D_01G * y[14] + cache$q01E_01G * y[18] + cache$q01F_01G * y[22] + cache$q01H_01G * y[30]
    dN10GdT <- cache$lambda10G * y[27] - cache$mu10G * y[27] - cache$q10G_00G * y[27] - cache$q10G_01G * y[27] - cache$q10G_11G * y[27] + cache$q00G_10G * y[25] + cache$q01G_10G * y[26] + cache$q11G_10G * y[28] - cache$q10G_10A * y[27] - cache$q10G_10B * y[27] - cache$q10G_10C * y[27] - cache$q10G_10D * y[27] - cache$q10G_10E * y[27] - cache$q10G_10F * y[27] - cache$q10G_10H * y[27] + cache$q10A_10G * y[3] + cache$q10B_10G * y[7] + cache$q10C_10G * y[11] + cache$q10D_10G * y[15] + cache$q10E_10G * y[19] + cache$q10F_10G * y[23] + cache$q10H_10G * y[31]
    dN11GdT <- cache$lambda11G * y[28] - cache$mu11G * y[28] - cache$q11G_00G * y[28] - cache$q11G_01G * y[28] - cache$q11G_10G * y[28] + cache$q00G_11G * y[25] + cache$q01G_11G * y[26] + cache$q10G_11G * y[27] - cache$q11G_11A * y[28] - cache$q11G_11B * y[28] - cache$q11G_11C * y[28] - cache$q11G_11D * y[28] - cache$q11G_11E * y[28] - cache$q11G_11F * y[28] - cache$q11G_11H * y[28] + cache$q11A_11G * y[4] + cache$q11B_11G * y[8] + cache$q11C_11G * y[12] + cache$q11D_11G * y[16] + cache$q11E_11G * y[20] + cache$q11F_11G * y[24] + cache$q11H_11G * y[32]
    
    #cat H:
    dN00HdT <- cache$lambda00H * y[29] - cache$mu00H * y[29] - cache$q00H_01H * y[29] - cache$q00H_10H * y[29] - cache$q00H_11H * y[29] + cache$q01H_00H * y[30] + cache$q10H_00H * y[31] + cache$q11H_00H * y[32] - cache$q00H_00A * y[29] - cache$q00H_00B * y[29] - cache$q00H_00C * y[29] - cache$q00H_00D * y[29] - cache$q00H_00E * y[29] - cache$q00H_00F * y[29] - cache$q00H_00G * y[29] + cache$q00A_00H * y[1] + cache$q00B_00H * y[5] + cache$q00C_00H * y[9] + cache$q00D_00H * y[13] + cache$q00E_00H * y[17] + cache$q00F_00H * y[21] + cache$q00G_00H * y[25]
    dN01HdT <- cache$lambda01H * y[30] - cache$mu01H * y[30] - cache$q01H_00H * y[30] - cache$q01H_10H * y[30] - cache$q01H_11H * y[30] + cache$q00H_01H * y[29] + cache$q10H_01H * y[31] + cache$q11H_01H * y[32] - cache$q01H_01A * y[30] - cache$q01H_01B * y[30] - cache$q01H_01C * y[30] - cache$q01H_01D * y[30] - cache$q01H_01E * y[30] - cache$q01H_01F * y[30] - cache$q01H_01G * y[30] + cache$q01A_01H * y[2] + cache$q01B_01H * y[6] + cache$q01C_01H * y[10] + cache$q01D_01H * y[14] + cache$q01E_01H * y[18] + cache$q01F_01H * y[22] + cache$q01G_01H * y[26]
    dN10HdT <- cache$lambda10H * y[31] - cache$mu10H * y[31] - cache$q10H_00H * y[31] - cache$q10H_01H * y[31] - cache$q10H_11H * y[31] + cache$q00H_10H * y[29] + cache$q01H_10H * y[30] + cache$q11H_10H * y[32] - cache$q10H_10A * y[31] - cache$q10H_10B * y[31] - cache$q10H_10C * y[31] - cache$q10H_10D * y[31] - cache$q10H_10E * y[31] - cache$q10H_10F * y[31] - cache$q10H_10G * y[31] + cache$q10A_10H * y[3] + cache$q10B_10H * y[7] + cache$q10C_10H * y[11] + cache$q10D_10H * y[15] + cache$q10E_10H * y[19] + cache$q10F_10H * y[23] + cache$q10G_10H * y[27]
    dN11HdT <- cache$lambda11H * y[32] - cache$mu11H * y[32] - cache$q11H_00H * y[32] - cache$q11H_01H * y[32] - cache$q11H_10H * y[32] + cache$q00H_11H * y[29] + cache$q01H_11H * y[30] + cache$q10H_11H * y[31] - cache$q11H_11A * y[32] - cache$q11H_11B * y[32] - cache$q11H_11C * y[32] - cache$q11H_11D * y[32] - cache$q11H_11E * y[32] - cache$q11H_11F * y[32] - cache$q11H_11G * y[32] + cache$q11A_11H * y[4] + cache$q11B_11H * y[8] + cache$q11C_11H * y[12] + cache$q11D_11G * y[16] + cache$q11E_11H * y[20] + cache$q11F_11H * y[24] + cache$q11G_11H * y[28]
    
    
    return(list(c(dN00AdT,dN01AdT,dN10AdT,dN11AdT, dN00BdT,dN01BdT,dN10BdT,dN11BdT, dN00CdT,dN01CdT,dN10CdT,dN11CdT, dN00DdT,dN01DdT,dN10DdT,dN11DdT, dN00EdT,dN01EdT,dN10EdT,dN11EdT, dN00FdT,dN01FdT,dN10FdT,dN11FdT, dN00GdT,dN01GdT,dN10GdT,dN11GdT, dN00HdT,dN01HdT,dN10HdT,dN11HdT)))
}





