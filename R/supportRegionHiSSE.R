######################################################################################################################################
######################################################################################################################################
### Adaptive Bootstrap fHiSSE -- Simulating confidence intervals for parameters estimated in the fast version of HiSSE
######################################################################################################################################
######################################################################################################################################

SupportRegionHiSSE <- function(hisse.obj, n.points=1000, scale.int=0.1, desired.delta=2, min.number.points=10, verbose=TRUE){
    phy <- hisse.obj$phy
    data <- hisse.obj$data
    data.new <- data.frame(data[,2], data[,2], row.names=data[,1])
    data.new <- data.new[phy$tip.label,]
    
    f <- hisse.obj$f
    np <- max(hisse.obj$index.par) - 1
    par <- numeric(np)
    free.parameters <- which(hisse.obj$index.par < max(hisse.obj$index.par))
    np.sequence <- 1:np
    
    for(i in np.sequence){
        par[i] <- hisse.obj$solution[which(hisse.obj$index.par == np.sequence[i])[1]]
    }
    hidden.states=hisse.obj$hidden.states
    condition.on.survival=hisse.obj$condition.on.survival
    root.type=hisse.obj$root.type
    root.p=hisse.obj$root.p
    includes.fossils = hisse.obj$includes.fossils
    fix.type = hisse.obj$fix.type
    
    lower <- exp(hisse.obj$lower.bounds)
    upper <- exp(hisse.obj$upper.bounds)
    
    #Bad Jeremy! Hard-coded column headers...
    interval.names <- c("lnLik", "turnover0A","turnover1A","eps0A","eps1A","q0A1A","q1A0A","q0A0B","q0A0C","q0A0D","q1A1B","q1A1C","q1A1D","turnover0B","turnover1B","eps0B","eps1B","q0B1B","q1B0B","q0B0A","q0B0C","q0B0D","q1B1A","q1B1C","q1B1D","turnover0C","turnover1C","eps0C","eps1C","q0C1C","q1C0C","q0C0A","q0C0B","q0C0D","q1C1A","q1C1B","q1C1D","turnover0D","turnover1D","eps0D","eps1D","q0D1D","q1D0D","q0D0A","q0D0B","q0D0C","q1D1A","q1D1B","q1D1C","psi")
    
    interval.results <- AdaptiveConfidenceIntervalSamplingfHiSSE(par, lower=lower, upper=upper, desired.delta = desired.delta, n.points=n.points, verbose=verbose, phy=phy, data=data.new, index.par=hisse.obj$index.par, f=f, hidden.states=hidden.states, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, includes.fossils=includes.fossils, fix.type=fix.type, scale.int=scale.int, min.number.points=min.number.points)
    interval.results.final <- matrix(0, n.points+1, length(hisse.obj$index.par))
    for(i in 1:(n.points+1)){
        par.rep <- unlist(interval.results[i,-1],use.names=FALSE)
        interval.results.final[i,] <- c(par.rep,0)[hisse.obj$index.par]
    }
    interval.results.final <- cbind(interval.results[,1], interval.results.final)
    interval.results.in <- interval.results.final[which(interval.results.final[,1] - min(interval.results.final[,1])<=desired.delta),]
    if(inherits(interval.results.in[1], what="numeric")){
        ci.interval = apply(interval.results.in, 2, quantile)
        colnames(interval.results.final) <- colnames(interval.results.in) <- colnames(ci.interval) <- interval.names
        obj = NULL
        obj$ci <- ci.interval
        obj$points.within.region = interval.results.in
        obj$all.points = interval.results.final
        class(obj) = "hisse.support"
        return(obj)
    }else{
        stop("Only the MLE is in the desired range. Try reducing scale.int.", call.=FALSE)
    }
}




AdaptiveConfidenceIntervalSamplingfHiSSE <- function(par, lower, upper, desired.delta=2, n.points=5000, verbose=TRUE, phy, data, index.par, f, hidden.states, condition.on.survival, root.type, root.p, includes.fossils, fix.type, scale.int, min.number.points=10) {
    
    # Some new prerequisites #
    gen <- FindGenerations(phy)
    dat.tab <- OrganizeDataHiSSE(data=data, phy=phy, f=f, hidden.states=hidden.states)
    if(includes.fossils == TRUE){
        fossil.taxa <- which(dat.tab$branch.type == 1)
    }else{
        fossil.taxa <- NULL
    }
    ##########################
    
    #Wrangle the data so that we can make use of DownPass easily:
    actual.params = which(index.par < max(index.par))
    model.vec <- numeric(length(index.par))
    model.vec[] <- c(par,0)[index.par]
    cache <- ParametersToPassfHiSSE(model.vec=model.vec, hidden.states=hidden.states, nb.tip=Ntip(phy), nb.node=Nnode(phy), bad.likelihood=exp(-300), f=f, ode.eps=0)
    
    if(includes.fossils == TRUE){
        starting <- -DownPassHiSSE(dat.tab=dat.tab, gen=gen, cache=cache, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, node=fix.type[,1], state=NULL, fossil.taxa=fossil.taxa, fix.type=fix.type[,2])
    }else{
        starting <- -DownPassHiSSE(dat.tab=dat.tab, gen=gen, cache=cache, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, node=NULL, state=NULL, fossil.taxa=fossil.taxa, fix.type=NULL)
    }
    
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
        cache <- ParametersToPassfHiSSE(model.vec=model.vec, hidden.states=hidden.states, nb.tip=Ntip(phy), nb.node=Nnode(phy), bad.likelihood=exp(-300), f=f, ode.eps=0)
        if(includes.fossils == TRUE){
            second <- -DownPassHiSSE(dat.tab=dat.tab, gen=gen, cache=cache, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, node=fix.type[,1], state=NULL, fossil.taxa=fossil.taxa, fix.type=fix.type[,2])
        }else{
            second <- -DownPassHiSSE(dat.tab=dat.tab, gen=gen, cache=cache, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, node=NULL, state=NULL, fossil.taxa=fossil.taxa, fix.type=NULL)
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
        scale.int <- scale.int*0.5
        new.results <- AdaptiveConfidenceIntervalSamplingfHiSSE(par, lower=lower, upper=upper, desired.delta = desired.delta, n.points=n.points, verbose=verbose, phy=phy, data=data, index.par=index.par, f=f, hidden.states=hidden.states, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, includes.fossils=includes.fossils, fix.type=fix.type, scale.int=scale.int, min.number.points=0)
        results <- rbind(results, new.results[-1,])
    }
    return(results)
}


print.hisse.support <- function(x,...){
    cat("\nSupport Region\n")
    print(x$ci[,-1])
    cat("\n")
}

