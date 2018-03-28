
######################################################################################################################################
######################################################################################################################################
### EquiFreq -- calculates a model-averaged equilibrium frequency for each state under a set of parameters in hiGeoSSE
######################################################################################################################################
######################################################################################################################################

GetModelAveEqFreqs <- function(x, max.time){
    res <- c()
    higeosse.results <- x
    if(class(higeosse.results)!="list") { #we have to make a list so we can run this generally
        tmp.list <- list()
        tmp.list[[1]] <- higeosse.results
        higeosse.results <- tmp.list
    }
    for(model.index in 1:length(higeosse.results)){
        if(higeosse.results[[model.index]]$assume.cladogenetic == TRUE){
            ##Modify the data file
            data.new <- data.frame(higeosse.results[[model.index]]$data[,2], higeosse.results[[model.index]]$data[,2], row.names=higeosse.results[[model.index]]$data[,1])
            data.new <- data.new[higeosse.results[[model.index]]$phy$tip.label,]
            cache = ParametersToPassHiGeoSSE(higeosse.results[[model.index]]$phy, data.new[,1], model.vec=higeosse.results[[model.index]]$solution, f=higeosse.results[[model.index]]$f, hidden.states=TRUE)
            if(higeosse.results[[model.index]]$root.type=="madfitz"){
                get.starting.probs <- DownPassHiGeosse(phy=higeosse.results[[model.index]]$phy, cache=cache, hidden.states=TRUE, condition.on.survival=higeosse.results[[model.index]]$condition.on.survival, root.type=higeosse.results[[model.index]]$root.type, root.p=higeosse.results[[model.index]]$root.p, get.phi=TRUE)$compD.root
            }
            
            out <- lsoda(c(state0A=get.starting.probs[1],state1A=get.starting.probs[2],state01A=get.starting.probs[3], state0B=get.starting.probs[4],state1B=get.starting.probs[5],state01B=get.starting.probs[6], state0C=get.starting.probs[7],state1C=get.starting.probs[8],state01C=get.starting.probs[9], state0D=get.starting.probs[10],state1D=get.starting.probs[11],state01D=get.starting.probs[12], state0E=get.starting.probs[13],state1E=get.starting.probs[14],state01E=get.starting.probs[15]), times=c(0, max.time), func=EqFreqsHiGeoSSE, parms=NULL, cache=cache, rtol=1e-8, atol=1e-8)[-1,-1]
        }else{
            data.new <- data.frame(higeosse.results[[model.index]]$data[,2], higeosse.results[[model.index]]$data[,2], row.names=higeosse.results[[model.index]]$data[,1])
            data.new <- data.new[higeosse.results[[model.index]]$phy$tip.label,]
            cache = ParametersToPassMuSSE(higeosse.results[[model.index]]$phy, data.new[,1], model.vec=higeosse.results[[model.index]]$solution, f=higeosse.results[[model.index]]$f, hidden.states=TRUE)
            if(higeosse.results[[model.index]]$root.type=="madfitz"){
                get.starting.probs <- DownPassMusse(phy=higeosse.results[[model.index]]$phy, cache=cache, hidden.states=TRUE, condition.on.survival=higeosse.results[[model.index]]$condition.on.survival, root.type=higeosse.results[[model.index]]$root.type, root.p=higeosse.results[[model.index]]$root.p, get.phi=TRUE)$compD.root
            }
            out <- lsoda(c(state0A=get.starting.probs[1],state1A=get.starting.probs[2],state01A=get.starting.probs[3], state0B=get.starting.probs[4],state1B=get.starting.probs[5],state01B=get.starting.probs[6], state0C=get.starting.probs[7],state1C=get.starting.probs[8],state01C=get.starting.probs[9], state0D=get.starting.probs[10],state1D=get.starting.probs[11],state01D=get.starting.probs[12], state0E=get.starting.probs[13],state1E=get.starting.probs[14],state01E=get.starting.probs[15]), times=c(0, max.time), func=EqFreqsMuSSE, parms=NULL, cache=cache, rtol=1e-8, atol=1e-8)[-1,-1]
        }
        out.mat <- t(matrix(out, 3, 5))
        colnames(out.mat) <- c("0", "1", "01")
        res <- rbind(res, colSums(out.mat)/sum(out.mat))
    }
    AIC.vector <- sapply(higeosse.results, "[[", "AIC")
    delta.AIC.vector <- AIC.vector - min(AIC.vector)
    rel.likelihood <- exp(-0.5 * delta.AIC.vector)
    AIC.weight.vector <- rel.likelihood / sum(rel.likelihood)
    final.eq.freq <- apply(res, 2, weighted.mean, w=AIC.weight.vector)
    return(final.eq.freq /sum(final.eq.freq))
}




#HiGeoSSE equilibrium freqs.
EqFreqsHiGeoSSE <- function(t, y, parms, cache){
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


