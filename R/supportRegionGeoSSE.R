######################################################################################################################################
######################################################################################################################################
### Adaptive Bootstrap GeoHiSSE -- Simulating confidence intervals for parameters estimated in GeoHiSSE
######################################################################################################################################
######################################################################################################################################

SupportRegionGeoSSE <- function(geohisse.obj, n.points=1000, scale.int=0.1, desired.delta=2, min.number.points=10, verbose=TRUE){
    phy <- geohisse.obj$phy
    data <- geohisse.obj$data
    data.new<-data.frame(data[,2], data[,2], row.names=data[,1])
    data.new<-data.new[phy$tip.label,]
    f <- geohisse.obj$f
    assume.cladogenetic <- geohisse.obj$assume.cladogenetic
    np <- max(geohisse.obj$index.par) - 1
    par <- numeric(np)
    free.parameters <- which(geohisse.obj$index.par < max(geohisse.obj$index.par))
    np.sequence <- 1:np
    for(i in np.sequence){
        par[i] <- geohisse.obj$solution[which(geohisse.obj$index.par == np.sequence[i])[1]]
    }
    hidden.states=geohisse.obj$hidden.states
    condition.on.survival=geohisse.obj$condition.on.survival
    root.type=geohisse.obj$root.type
    root.p=geohisse.obj$root.p
    
    lower <- exp(geohisse.obj$lower.bounds)
    upper <- exp(geohisse.obj$upper.bounds)
    
    #Bad Jeremy! Hard-coded column headers...
    #if(assume.cladogenetic==TRUE){
    #    interval.names <- c("lnLik", "s0A", "s1A", "s01A", "x0A", "x1A", "d0A_1A ", "d0A_01A", "d1A_0A", "d1A_01A", "d01A_0A", "d01A_1A", "d0A_0B", "d0A_0C", "d0A_0D", "d0A_0E", "d1A_1B", "d1A_1C", "d1A_1D", "d1A_1E", "d01A_01B", "d01A_01C", "d01A_01D", "d01A_01E", "s0B", "s1B", "s01B", "x0B", "x1B", "d0B_1B ", "d0B_01B", "d1B_0B", "d1B_01B", "d01B_0B", "d01B_1B", "d0B_0A", "d0B_0C", "d0B_0D", "d0B_0E", "d1B_1A", "d1B_1C", "d1B_1D", "d1B_1E", "d01B_01A", "d01B_01C", "d01B_01D", "d01B_01E", "s0C", "s1C", "s01C", "x0C", "x1C", "d0C_1C ", "d0C_01C", "d1C_0C", "d1C_01C", "d01C_0C", "d01C_1C", "d0C_0A", "d0C_0B", "d0C_0D", "d0C_0E", "d1C_1A", "d1C_1B", "d1C_1D", "d1C_1E", "d01C_01A", "d01C_01B", "d01C_01D", "d01C_01E", "s0D", "s1D", "s01D", "x0D", "x1D", "d0D_1D ", "d0D_01D", "d1D_0D", "d1D_01D", "d01D_0D", "d01D_1D", "d0D_0A", "d0D_0B", "d0D_0C", "d0D_0E", "d1D_1A", "d1D_1B", "d1D_1C", "d1D_1E", "d01D_01A", "d01D_01B", "d01D_01C", "d01D_01E", "s0E", "s1E", "s01E", "x0E", "x1E", "d0E_1E ", "d0E_01E", "d1E_0E", "d1E_01E", "d01E_0E", "d01E_1E", "d0E_0A", "d0E_0B", "d0E_0C", "d0E_0D", "d1E_1A", "d1E_1B", "d1E_1C", "d1E_1D", "d01E_01A", "d01E_01B", "d01E_01C", "d01E_01D")
    #}else{
    #    interval.names <- c("lnLik", "s0A", "s1A", "s01A", "x0A", "x1A", "x01A", "d0A_1A ", "d0A_01A", "d1A_0A", "d1A_01A", "d01A_0A", "d01A_1A", "d0A_0B", "d0A_0C", "d0A_0D", "d0A_0E", "d1A_1B", "d1A_1C", "d1A_1D", "d1A_1E", "d01A_01B", "d01A_01C", "d01A_01D", "d01A_01E", "s0B", "s1B", "s01B", "x0B", "x1B", "x01B", "d0B_1B ", "d0B_01B", "d1B_0B", "d1B_01B", "d01B_0B", "d01B_1B", "d0B_0A", "d0B_0C", "d0B_0D", "d0B_0E", "d1B_1A", "d1B_1C", "d1B_1D", "d1B_1E", "d01B_01A", "d01B_01C", "d01B_01D", "d01B_01E", "s0C", "s1C", "s01C", "x0C", "x1C", "x01C", "d0C_1C ", "d0C_01C", "d1C_0C", "d1C_01C", "d01C_0C", "d01C_1C", "d0C_0A", "d0C_0B", "d0C_0D", "d0C_0E", "d1C_1A", "d1C_1B", "d1C_1D", "d1C_1E", "d01C_01A", "d01C_01B", "d01C_01D", "d01C_01E", "s0D", "s1D", "s01D", "x0D", "x1D", "x01D", "d0D_1D ", "d0D_01D", "d1D_0D", "d1D_01D", "d01D_0D", "d01D_1D", "d0D_0A", "d0D_0B", "d0D_0C", "d0D_0E", "d1D_1A", "d1D_1B", "d1D_1C", "d1D_1E", "d01D_01A", "d01D_01B", "d01D_01C", "d01D_01E", "s0E", "s1E", "s01E", "x0E", "x1E", "x01E", "d0E_1E ", "d0E_01E", "d1E_0E", "d1E_01E", "d01E_0E", "d01E_1E", "d0E_0A", "d0E_0B", "d0E_0C", "d0E_0D", "d1E_1A", "d1E_1B", "d1E_1C", "d1E_1D", "d01E_01A", "d01E_01B", "d01E_01C", "d01E_01D")
    #}
    
    interval.names <- c("lnLik", "tau00A","tau11A","tau01A","ef00A","ef11A","d00A_11A","d00A_01A","d11A_00A","d11A_01A","d01A_00A","d01A_11A","d00A_00B","d00A_00C","d00A_00D","d00A_00E","d00A_00F","d00A_00G","d00A_00H","d00A_00I","d00A_00J","d11A_11B","d11A_11C","d11A_11D","d11A_11E","d11A_11F","d11A_11G","d11A_11H","d11A_11I","d11A_11J","d01A_01B","d01A_01C","d01A_01D","d01A_01E","d01A_01F","d01A_01G","d01A_01H","d01A_01I","d01A_01J","tau00B","tau11B","tau01B","ef00B","ef11B","d00B_11B","d00B_01B","d11B_00B","d11B_01B","d01B_00B","d01B_11B","d00B_00A","d00B_00C","d00B_00D","d00B_00E","d00B_00F","d00B_00G","d00B_00H","d00B_00I","d00B_00J","d11B_11A","d11B_11C","d11B_11D","d11B_11E","d11B_11F","d11B_11G","d11B_11H","d11B_11I","d11B_11J","d01B_01A","d01B_01C","d01B_01D","d01B_01E","d01B_01F","d01B_01G","d01B_01H","d01B_01I","d01B_01J","tau00C","tau11C","tau01C","ef00C","ef11C","d00C_11C","d00C_01C","d11C_00C","d11C_01C","d01C_00C","d01C_11C","d00C_00A","d00C_00B","d00C_00D","d00C_00E","d00C_00F","d00C_00G","d00C_00H","d00C_00I","d00C_00J","d11C_11A","d11C_11B","d11C_11D","d11C_11E","d11C_11F","d11C_11G","d11C_11H","d11C_11I","d11C_11J","d01C_01A","d01C_01B","d01C_01D","d01C_01E","d01C_01F","d01C_01G","d01C_01H","d01C_01I","d01C_01J","tau00D","tau11D","tau01D","ef00D","ef11D","d00D_11D","d00D_01D","d11D_00D","d11D_01D","d01D_00D","d01D_11D","d00D_00A","d00D_00B","d00D_00C","d00D_00E","d00D_00F","d00D_00G","d00D_00H","d00D_00I","d00D_00J","d11D_11A","d11D_11B","d11D_11C","d11D_11E","d11D_11F","d11D_11G","d11D_11H","d11D_11I","d11D_11J","d01D_01A","d01D_01B","d01D_01C","d01D_01E","d01D_01F","d01D_01G","d01D_01H","d01D_01I","d01D_01J","tau00E","tau11E","tau01E","ef00E","ef11E","d00E_11E","d00E_01E","d11E_00E","d11E_01E","d01E_00E","d01E_11E","d00E_00A","d00E_00B","d00E_00C","d00E_00D","d00E_00F","d00E_00G","d00E_00H","d00E_00I","d00E_00J","d11E_11A","d11E_11B","d11E_11C","d11E_11D","d11E_11F","d11E_11G","d11E_11H","d11E_11I","d11E_11J","d01E_01A","d01E_01B","d01E_01C","d01E_01D","d01E_01F","d01E_01G","d01E_01H","d01E_01I","d01E_01J","tau00F","tau11F","tau01F","ef00F","ef11F","d00F_11F","d00F_01F","d11F_00F","d11F_01F","d01F_00F","d01F_11F","d00F_00A","d00F_00B","d00F_00C","d00F_00D","d00F_00E","d00F_00G","d00F_00H","d00F_00I","d00F_00J","d11F_11A","d11F_11B","d11F_11C","d11F_11D","d11F_11E","d11F_11G","d11F_11H","d11F_11I","d11F_11J","d01F_01A","d01F_01B","d01F_01C","d01F_01D","d01F_01E","d01F_01G","d01F_01H","d01F_01I","d01F_01J","tau00G","tau11G","tau01G","ef00G","ef11G","d00G_11G","d00G_01G","d11G_00G","d11G_01G","d01G_00G","d01G_11G","d00G_00A","d00G_00B","d00G_00C","d00G_00D","d00G_00E","d00G_00F","d00G_00H","d00G_00I","d00G_00J","d11G_11A","d11G_11B","d11G_11C","d11G_11D","d11G_11E","d11G_11F","d11G_11H","d11G_11I","d11G_11J","d01G_01A","d01G_01B","d01G_01C","d01G_01D","d01G_01E","d01G_01F","d01G_01H","d01G_01I","d01G_01J","tau00H","tau11H","tau01H","ef00H","ef11H","d00H_11H","d00H_01H","d11H_00H","d11H_01H","d01H_00H","d01H_11H","d00H_00A","d00H_00B","d00H_00C","d00H_00D","d00H_00E","d00H_00F","d00H_00G","d00H_00I","d00H_00J","d11H_11A","d11H_11B","d11H_11C","d11H_11D","d11H_11E","d11H_11F","d11H_11G","d11H_11I","d11H_11J","d01H_01A","d01H_01B","d01H_01C","d01H_01D","d01H_01E","d01H_01F","d01H_01G","d01H_01I","d01H_01J","tau00I","tau11I","tau01I","ef00I","ef11I","d00I_11I","d00I_01I","d11I_00I","d11I_01I","d01I_00I","d01I_11I","d00I_00A","d00I_00B","d00I_00C","d00I_00D","d00I_00E","d00I_00F","d00I_00G","d00I_00H","d00I_00J","d11I_11A","d11I_11B","d11I_11C","d11I_11D","d11I_11E","d11I_11F","d11I_11G","d11I_11H","d11I_11J","d01I_01A","d01I_01B","d01I_01C","d01I_01D","d01I_01E","d01I_01F","d01I_01G","d01I_01H","d01I_01J","tau00J","tau11J","tau01J","ef00J","ef11J","d00J_11J","d00J_01J","d11J_00J","d11J_01J","d01J_00J","d01J_11J","d00J_00A","d00J_00B","d00J_00C","d00J_00D","d00J_00E","d00J_00F","d00J_00G","d00J_00H","d00J_00I","d11J_11A","d11J_11B","d11J_11C","d11J_11D","d11J_11E","d11J_11F","d11J_11G","d11J_11H","d11J_11I","d01J_01A","d01J_01B","d01J_01C","d01J_01D","d01J_01E","d01J_01F","d01J_01G","d01J_01H","d01J_01I")
    
    interval.results <- AdaptiveConfidenceIntervalSamplingGeoHiSSE(par, lower=lower, upper=upper, desired.delta = desired.delta, n.points=n.points, verbose=verbose, phy=phy, data=data.new, index.par=geohisse.obj$index.par, f=f, hidden.states=hidden.states, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, scale.int=scale.int, assume.cladogenetic=assume.cladogenetic, min.number.points=min.number.points)
    interval.results.final <- matrix(0, n.points+1, length(geohisse.obj$index.par))
    for(i in 1:(n.points+1)){
        par.rep <- unlist(interval.results[i,-1],use.names=FALSE)
        interval.results.final[i,] <- c(par.rep,0)[geohisse.obj$index.par]
    }
    interval.results.final <- cbind(interval.results[,1], interval.results.final)
    interval.results.in <- interval.results.final[which(interval.results.final[,1] - min(interval.results.final[,1])<=desired.delta),]
    if(class(interval.results.in)=="numeric"){
        stop("Only the MLE is in the desired range. Try reducing scale.int.", call.=FALSE)
    }else{
        ci.interval = apply(interval.results.in, 2, quantile)
        colnames(interval.results.final) <- colnames(interval.results.in) <- colnames(ci.interval) <- interval.names
        obj = NULL
        obj$ci <- ci.interval
        obj$points.within.region = interval.results.in
        obj$all.points = interval.results.final
        class(obj) = "geohisse.support"
        return(obj)
    }
}


AdaptiveConfidenceIntervalSamplingGeoHiSSE <- function(par, lower, upper, desired.delta=2, n.points=5000, verbose=TRUE, phy, data, index.par, f, hidden.states, condition.on.survival, root.type, root.p, scale.int, assume.cladogenetic=TRUE, min.number.points=10) {
    
    # Some new prerequisites #
    gen <- FindGenerations(phy)
    dat.tab <- OrganizeDataGeo(data=data[,1], phy=phy, f=f, hidden.states=hidden.states)
    ##########################
    
    #Wrangle the data so that we can make use of DownPass easily:
    actual.params = which(index.par < max(index.par))
    model.vec <- numeric(length(index.par))
    model.vec[] <- c(par,0)[index.par]
    cache <- ParametersToPassGeoHiSSEfast(model.vec=model.vec, hidden.states=hidden.states, assume.cladogenetic=assume.cladogenetic, nb.tip=Ntip(phy), nb.node=Nnode(phy), bad.likelihood=exp(-300), ode.eps=0)
    phy$node.label <- NULL
    starting <- -DownPassGeoHissefast(dat.tab=dat.tab, gen=gen, cache=cache, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p)

    #Generate the multipliers for feeling the boundaries:
    min.multipliers <- rep(1, length(par))
    max.multipliers <- rep(1, length(par))
    results <- data.frame(data.frame(matrix(nrow=n.points+1, ncol=1+length(par))))
    results[1,] <- unname(c(starting, par))
    for (i in sequence(n.points)) {
        sim.points <- NA
        while(is.na(sim.points[1])) {
            sim.points <- GenerateValues(par, lower=lower, upper=upper, scale.int=scale.int, examined.max=max.multipliers*apply(results[which(results[,1]-min(results[,1], na.rm=TRUE)<=desired.delta),-1], 2, max, na.rm=TRUE), examined.min=min.multipliers*apply(results[which(results[,1]-min(results[,1], na.rm=TRUE)<=desired.delta),-1], 2, min, na.rm=TRUE))
        }
        par <- sim.points
        model.vec <- numeric(length(index.par))
        model.vec[] <- c(sim.points,0)[index.par]
        cache <- ParametersToPassGeoHiSSEfast(model.vec=model.vec, hidden.states=hidden.states, assume.cladogenetic=assume.cladogenetic, nb.tip=Ntip(phy), nb.node=Nnode(phy), bad.likelihood=exp(-300), ode.eps=0)
        phy$node.label <- NULL
        second <- -DownPassGeoHissefast(dat.tab=dat.tab, gen=gen, cache=cache, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p)
        results[i+1,] <- c(second, sim.points)
        if(i%%20==0) {
            for (j in sequence(length(par))) {
                returned.range <- range(results[which((results[,1]-min(results[,1], na.rm=TRUE))<desired.delta), j+1], na.rm=TRUE)
                total.range <- range(results[,j+1], na.rm=TRUE)
                width.ratio <- diff(returned.range)/diff(total.range)
                if(is.na(width.ratio)) {
                    width.ratio=1
                }
                if(width.ratio > 0.5) { #we are not sampling widely enough
                    min.multipliers[j] <- min.multipliers[j] * (1-scale.int)
                    max.multipliers[j] <- max.multipliers[j] * (1+scale.int) #expand the range
                } else {
                    if(width.ratio < 0.02 & i>100) { # we are sampling too widely
                        min.multipliers[j] <- min.multipliers[j] * (1+scale.int)
                        max.multipliers[j] <- max.multipliers[j] * (1-scale.int) #contract the range
                    } else {
                        min.multipliers[j] <- 1
                        max.multipliers[j] <- 1
                    }
                }
            }
        }
        if (verbose && i%%100==0) {
            cat(paste(i, "of", n.points, "points done"), "\n")
        }
    }
    while(length(which((results[,1]-min(results[,1], na.rm=TRUE))<desired.delta))<min.number.points) {
        warning("Did not generate enough points in the region; restarting to create additional points")
        print(paste("Now doing an additional", 2+round(n.points/4), "points to the", dim(results)[1], "ones already done because not enough points in the good enough region were sampled"))
        scale.int <- scale.int*0.5
        new.results <- AdaptiveConfidenceIntervalSamplingGeoHiSSE(par, lower=lower, upper=upper, desired.delta = desired.delta, n.points=n.points, verbose=verbose, phy=phy, data=data, index.par=index.par, f=f, hidden.states=hidden.states, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, scale.int=scale.int, assume.cladogenetic=assume.cladogenetic, min.number.points=0)
        results <- rbind(results, new.results[-1,])
    }
    return(results)
}


print.geohisse.support <- function(x,...){
    
    cat("\nSupport Region\n")
    print(x$ci[,-1])
    cat("\n")
}


