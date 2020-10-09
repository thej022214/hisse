######################################################################################################################################
######################################################################################################################################
### Adaptive Bootstrap GeoHiSSE -- Simulating confidence intervals for parameters estimated in GeoHiSSE
######################################################################################################################################
######################################################################################################################################

SupportRegionGeoSSE.old <- function(geohisse.old.obj, n.points=1000, scale.int=0.1, desired.delta=2, min.number.points=10, verbose=TRUE){
    phy <- geohisse.old.obj$phy
    data <- geohisse.old.obj$data
    data.new<-data.frame(data[,2], data[,2], row.names=data[,1])
    data.new<-data.new[phy$tip.label,]
    f <- geohisse.old.obj$f
    assume.cladogenetic <- geohisse.old.obj$assume.cladogenetic
    np <- max(geohisse.old.obj$index.par) - 1
    par <- numeric(np)
    free.parameters <- which(geohisse.old.obj$index.par < max(geohisse.old.obj$index.par))
    np.sequence <- 1:np
    for(i in np.sequence){
        par[i] <- geohisse.old.obj$solution[which(geohisse.old.obj$index.par == np.sequence[i])[1]]
    }
    hidden.states=geohisse.old.obj$hidden.areas
    condition.on.survival=geohisse.old.obj$condition.on.survival
    root.type=geohisse.old.obj$root.type
    root.p=geohisse.old.obj$root.p

    lower <- exp(geohisse.old.obj$lower.bounds)
    upper <- exp(geohisse.old.obj$upper.bounds)

    #Bad Jeremy! Hard-coded column headers...
    if(assume.cladogenetic==TRUE){
        interval.names <- c("lnLik", "s0A", "s1A", "s01A", "x0A", "x1A", "d0A_1A ", "d0A_01A", "d1A_0A", "d1A_01A", "d01A_0A", "d01A_1A", "d0A_0B", "d0A_0C", "d0A_0D", "d0A_0E", "d1A_1B", "d1A_1C", "d1A_1D", "d1A_1E", "d01A_01B", "d01A_01C", "d01A_01D", "d01A_01E", "s0B", "s1B", "s01B", "x0B", "x1B", "d0B_1B ", "d0B_01B", "d1B_0B", "d1B_01B", "d01B_0B", "d01B_1B", "d0B_0A", "d0B_0C", "d0B_0D", "d0B_0E", "d1B_1A", "d1B_1C", "d1B_1D", "d1B_1E", "d01B_01A", "d01B_01C", "d01B_01D", "d01B_01E", "s0C", "s1C", "s01C", "x0C", "x1C", "d0C_1C ", "d0C_01C", "d1C_0C", "d1C_01C", "d01C_0C", "d01C_1C", "d0C_0A", "d0C_0B", "d0C_0D", "d0C_0E", "d1C_1A", "d1C_1B", "d1C_1D", "d1C_1E", "d01C_01A", "d01C_01B", "d01C_01D", "d01C_01E", "s0D", "s1D", "s01D", "x0D", "x1D", "d0D_1D ", "d0D_01D", "d1D_0D", "d1D_01D", "d01D_0D", "d01D_1D", "d0D_0A", "d0D_0B", "d0D_0C", "d0D_0E", "d1D_1A", "d1D_1B", "d1D_1C", "d1D_1E", "d01D_01A", "d01D_01B", "d01D_01C", "d01D_01E", "s0E", "s1E", "s01E", "x0E", "x1E", "d0E_1E ", "d0E_01E", "d1E_0E", "d1E_01E", "d01E_0E", "d01E_1E", "d0E_0A", "d0E_0B", "d0E_0C", "d0E_0D", "d1E_1A", "d1E_1B", "d1E_1C", "d1E_1D", "d01E_01A", "d01E_01B", "d01E_01C", "d01E_01D")
    }else{
        interval.names <- c("lnLik", "s0A", "s1A", "s01A", "x0A", "x1A", "x01A", "d0A_1A ", "d0A_01A", "d1A_0A", "d1A_01A", "d01A_0A", "d01A_1A", "d0A_0B", "d0A_0C", "d0A_0D", "d0A_0E", "d1A_1B", "d1A_1C", "d1A_1D", "d1A_1E", "d01A_01B", "d01A_01C", "d01A_01D", "d01A_01E", "s0B", "s1B", "s01B", "x0B", "x1B", "x01B", "d0B_1B ", "d0B_01B", "d1B_0B", "d1B_01B", "d01B_0B", "d01B_1B", "d0B_0A", "d0B_0C", "d0B_0D", "d0B_0E", "d1B_1A", "d1B_1C", "d1B_1D", "d1B_1E", "d01B_01A", "d01B_01C", "d01B_01D", "d01B_01E", "s0C", "s1C", "s01C", "x0C", "x1C", "x01C", "d0C_1C ", "d0C_01C", "d1C_0C", "d1C_01C", "d01C_0C", "d01C_1C", "d0C_0A", "d0C_0B", "d0C_0D", "d0C_0E", "d1C_1A", "d1C_1B", "d1C_1D", "d1C_1E", "d01C_01A", "d01C_01B", "d01C_01D", "d01C_01E", "s0D", "s1D", "s01D", "x0D", "x1D", "x01D", "d0D_1D ", "d0D_01D", "d1D_0D", "d1D_01D", "d01D_0D", "d01D_1D", "d0D_0A", "d0D_0B", "d0D_0C", "d0D_0E", "d1D_1A", "d1D_1B", "d1D_1C", "d1D_1E", "d01D_01A", "d01D_01B", "d01D_01C", "d01D_01E", "s0E", "s1E", "s01E", "x0E", "x1E", "x01E", "d0E_1E ", "d0E_01E", "d1E_0E", "d1E_01E", "d01E_0E", "d01E_1E", "d0E_0A", "d0E_0B", "d0E_0C", "d0E_0D", "d1E_1A", "d1E_1B", "d1E_1C", "d1E_1D", "d01E_01A", "d01E_01B", "d01E_01C", "d01E_01D")
    }

    interval.results <- AdaptiveConfidenceIntervalSamplingGeoHiSSE.old(par, lower=lower, upper=upper, desired.delta = desired.delta, n.points=n.points, verbose=verbose, phy=phy, data=data.new, index.par=geohisse.old.obj$index.par, f=f, hidden.states=hidden.states, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, scale.int=scale.int, assume.cladogenetic=assume.cladogenetic, min.number.points=min.number.points)
    interval.results.final <- matrix(0, n.points+1, length(geohisse.old.obj$index.par))
    for(i in 1:(n.points+1)){
        par.rep <- unlist(interval.results[i,-1],use.names=FALSE)
        interval.results.final[i,] <- c(par.rep,0)[geohisse.old.obj$index.par]
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




AdaptiveConfidenceIntervalSamplingGeoHiSSE.old <- function(par, lower, upper, desired.delta=2, n.points=5000, verbose=TRUE, phy, data, index.par, f, hidden.states, condition.on.survival, root.type, root.p, scale.int, assume.cladogenetic=TRUE, min.number.points=10) {
    #Wrangle the data so that we can make use of DownPass easily:
    actual.params = which(index.par < max(index.par))
    model.vec <- numeric(length(index.par))
    model.vec[] <- c(par,0)[index.par]
    if(assume.cladogenetic == TRUE){
        cache = ParametersToPassGeoHiSSE(phy, data[,1], f, model.vec=model.vec, hidden.states=hidden.states)
        phy$node.label <- NULL
        starting <- -DownPassGeoHisse(phy=phy, cache=cache, hidden.states=TRUE, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p)
    }else{
        cache = ParametersToPassMuSSE(phy, data[,1], f, model.vec=model.vec, hidden.states=hidden.states)
        phy$node.label <- NULL
        starting <- -DownPassMusse(phy=phy, cache=cache, hidden.states=TRUE, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p)
    }
    #Generate the multipliers for feeling the boundaries:
    min.multipliers <- rep(1, length(par))
    max.multipliers <- rep(1, length(par))
    results <- data.frame(data.frame(matrix(nrow=n.points+1, ncol=1+length(par))))
    results[1,] <- unname(c(starting, par))
    for (i in sequence(n.points)) {
        sim.points <- NA
        while(is.na(sim.points[1])) {
            sim.points <- GenerateValuesGeo(par, lower=lower, upper=upper, scale.int=scale.int, examined.max=max.multipliers*apply(results[which(results[,1]-min(results[,1], na.rm=TRUE)<=desired.delta),-1], 2, max, na.rm=TRUE), examined.min=min.multipliers*apply(results[which(results[,1]-min(results[,1], na.rm=TRUE)<=desired.delta),-1], 2, min, na.rm=TRUE))
        }
        par <- sim.points
        model.vec <- numeric(length(index.par))
        model.vec[] <- c(sim.points,0)[index.par]
        if(assume.cladogenetic == TRUE){
            cache = ParametersToPassGeoHiSSE(phy, data[,1], f, model.vec=model.vec, hidden.states=hidden.states)
            phy$node.label <- NULL
            second <- -DownPassGeoHisse(phy=phy, cache=cache, hidden.states=TRUE, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p)
        }else{
            cache = ParametersToPassMuSSE(phy, data[,1], f, model.vec=model.vec, hidden.states=hidden.states)
            phy$node.label <- NULL
            second <- -DownPassMusse(phy=phy, cache=cache, hidden.states=TRUE, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p)
        }
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
        new.results <- AdaptiveConfidenceIntervalSamplingGeoHiSSE.old(par, lower=lower, upper=upper, desired.delta = desired.delta, n.points=n.points, verbose=verbose, phy=phy, data=data, index.par=index.par, f=f, hidden.states=hidden.states, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, scale.int=scale.int, assume.cladogenetic=assume.cladogenetic, min.number.points=0)
        results <- rbind(results, new.results[-1,])
    }
    return(results)
}

## Strangely this function is not present in the development branch, only in the master branch.
## DANIEL: I can see that this function is used once in this file, so I will keep it here.
GenerateValuesGeo <- function(par, lower, upper, scale.int, max.tries=100, expand.prob=0, examined.max, examined.min) {
    pass=FALSE
    tries=0
    while(!pass && tries<=max.tries) {
        tries <- tries+1
        pass=TRUE
        new.vals <- rep(NA, length(par))
        for(i in sequence(length(par))) {
            examined.max[i] <- max(0.001, examined.max[i])
            min.val <- min(max(lower[i], (1-scale.int)*examined.min[i]), examined.max[i]) #just in case min is greater than max
            max.val <- max(min(upper[i], (1+scale.int)*examined.max[i]), examined.min[i])
            if(isTRUE(all.equal(min.val, max.val))) {
                min.val <- min.val * 0.9999
                max.val <- max.val * 1.0001
            }
            new.vals[i] <- runif(1, min.val, max.val)
            if(new.vals[i]<lower[i]) {
                pass=FALSE
            }
            if(new.vals[i]>upper[i]) {
                pass=FALSE
            }
        }
    }
    if(tries>max.tries) {
        return(NA)
    }
    return(new.vals)
}

print.geohisse.support <- function(x,...){

    cat("\nSupport Region\n")
    print(x$ci[,-1])
    cat("\n")
}
