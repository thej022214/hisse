######################################################################################################################################
######################################################################################################################################
### Adaptive Bootstrap MiSSE -- Simulating confidence intervals for parameters estimated in MiSSE
######################################################################################################################################
######################################################################################################################################

SupportRegionMiSSE <- function(misse.obj, n.points=1000, scale.int=0.1, desired.delta=2, min.number.points=10, verbose=TRUE){
    
    #If no fossils should be fine, but necessary to call this tree if fossils are present:
    phy <- misse.obj$phy.w.k
    
    f <- misse.obj$f
    np <- max(misse.obj$index.par) - 1
    par <- numeric(np)
    free.parameters <- which(misse.obj$index.par < max(misse.obj$index.par))
    np.sequence <- 1:np
    
    for(i in np.sequence){
        par[i] <- misse.obj$solution[which(misse.obj$index.par == np.sequence[i])[1]]
    }
    hidden.states=misse.obj$hidden.states
    condition.on.survival=misse.obj$condition.on.survival
    root.type=misse.obj$root.type
    root.p=misse.obj$root.p
    includes.fossils = misse.obj$includes.fossils
    fix.type = misse.obj$fix.type
    strat.intervals = misse.obj$strat.intervals
    
    lower <- exp(misse.obj$lower.bounds)
    upper <- exp(misse.obj$upper.bounds)
    
    #Bad Jeremy! Hard-coded column headers...
    interval.names <- c("lnLik", "turnover0A","eps0A", "turnover0B","eps0B", "turnover0C","eps0C", "turnover0D","eps0D", "turnover0E","eps0E", "turnover0F","eps0F", "turnover0G","eps0G", "turnover0H","eps0H", "turnover0I","eps0I", "turnover0J","eps0J", "turnover0K","eps0K", "turnover0L","eps0L", "turnover0M","eps0M", "turnover0N","eps0N", "turnover0O","eps0O", "turnover0P","eps0P", "turnover0Q","eps0Q", "turnover0R","eps0R", "turnover0S","eps0S", "turnover0T","eps0T", "turnover0U","eps0U", "turnover0V","eps0V","turnover0W","eps0W","turnover0X","eps0X", "turnover0Y","eps0Y", "turnover0Z","eps0Z", "q0", "psi")
    
    interval.results <- AdaptiveConfidenceIntervalSamplingMiSSE(par, lower=lower, upper=upper, desired.delta = desired.delta, n.points=n.points, verbose=verbose, phy=phy, index.par=misse.obj$index.par, f=f, hidden.states=hidden.states, fixed.eps=misse.obj$fixed.eps, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, includes.fossils=includes.fossils, strat.intervals=strat.intervals, fix.type=fix.type, scale.int=scale.int, min.number.points=min.number.points)
    interval.results.final <- matrix(0, n.points+1, length(misse.obj$index.par))
    for(i in 1:(n.points+1)){
        par.rep <- unlist(interval.results[i,-1],use.names=FALSE)
        interval.results.final[i,] <- c(par.rep,0)[misse.obj$index.par]
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
        class(obj) = "misse.support"
        return(obj)
    }else{
        stop("Only the MLE is in the desired range. Try reducing scale.int.", call.=FALSE)
    }
}


AdaptiveConfidenceIntervalSamplingMiSSE <- function(par, lower, upper, desired.delta=2, n.points=5000, verbose=TRUE, phy, index.par, f, hidden.states, fixed.eps, condition.on.survival, root.type, root.p, includes.fossils, strat.intervals, fix.type, scale.int, min.number.points=10) {
    
    # Some new prerequisites #
    if(includes.fossils == TRUE){
        if(!is.null(strat.intervals)){
            strat.cache <- GetStratInfo(strat.intervals=strat.intervals)
            gen <- FindGenerations(phy)
            dat.tab <- OrganizeDataMiSSE(phy=phy, f=f, hidden.states=hidden.states, includes.intervals=TRUE, intervening.intervals=strat.cache$intervening.intervals)
            #These are all inputs for generating starting values:
            edge_details <- GetEdgeDetails(phy, includes.intervals=TRUE, intervening.intervals=strat.cache$intervening.intervals)
            fossil.taxa <- edge_details$tipward_node[which(edge_details$type == "extinct_tip" | edge_details$type == "k_extinct_interval")]
        }else{
            gen <- FindGenerations(phy)
            dat.tab <- OrganizeDataMiSSE(phy=phy, f=f, hidden.states=hidden.states, includes.intervals=FALSE, intervening.intervals=NULL)
            #These are all inputs for generating starting values:
            edge_details <- GetEdgeDetails(phy, includes.intervals=FALSE, intervening.intervals=NULL)
            fossil.taxa <- edge_details$tipward_node[which(edge_details$type == "extinct_tip")]
        }
    }else{
        fossil.taxa <- NULL
        gen <- FindGenerations(phy)
        dat.tab <- OrganizeDataMiSSE(phy=phy, f=f, hidden.states=hidden.states)
    }
    ##########################
    
    #Wrangle the data so that we can make use of DownPass easily:
    actual.params <- which(index.par < max(index.par))
    model.vec <- numeric(length(index.par))
    model.vec[] <- c(par,0)[index.par]
    cache <- ParametersToPassMiSSE(model.vec=model.vec, hidden.states=hidden.states, fixed.eps=fixed.eps, nb.tip=Ntip(phy), nb.node=Nnode(phy), bad.likelihood=exp(-300), ode.eps=0)
    
    if(includes.fossils == TRUE){
        if(!is.null(strat.intervals)){
            starting <- -(DownPassMisse(dat.tab=dat.tab, gen=gen, cache=cache, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, node=fix.type[,1], state=NULL, fossil.taxa=fossil.taxa, fix.type=fix.type[,2]) + (strat.cache$k*log(cache$psi)) + (cache$psi*strat.cache$l_s))
        }else{
            starting <- -DownPassMisse(dat.tab=dat.tab, gen=gen, cache=cache, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, node=fix.type[,1], state=NULL, fossil.taxa=fossil.taxa, fix.type=fix.type[,2])
        }
    }else{
        starting <- -DownPassMisse(dat.tab=dat.tab, gen=gen, cache=cache, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, node=NULL, fossil.taxa=fossil.taxa, fix.type=NULL)
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
        cache <- ParametersToPassMiSSE(model.vec=model.vec, hidden.states=hidden.states, fixed.eps=fixed.eps, nb.tip=Ntip(phy), nb.node=Nnode(phy), bad.likelihood=exp(-300), ode.eps=0)
        if(includes.fossils == TRUE){
            if(!is.null(strat.intervals)){
                second <- -(DownPassMisse(dat.tab=dat.tab, gen=gen, cache=cache, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, node=fix.type[,1], state=NULL, fossil.taxa=fossil.taxa, fix.type=fix.type[,2]) + (strat.cache$k*log(cache$psi)) + (cache$psi*strat.cache$l_s))
            }else{
                second <- -DownPassMisse(dat.tab=dat.tab, gen=gen, cache=cache, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, node=fix.type[,1], state=NULL, fossil.taxa=fossil.taxa, fix.type=fix.type[,2])
            }
        }else{
            second <- -DownPassMisse(dat.tab=dat.tab, gen=gen, cache=cache, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, node=NULL, fossil.taxa=fossil.taxa, fix.type=NULL)
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
        new.results <- AdaptiveConfidenceIntervalSamplingMiSSE(par, lower=lower, upper=upper, desired.delta = desired.delta, n.points=n.points, verbose=verbose, phy=phy, index.par=index.par, f=f, hidden.states=hidden.states, fixed.eps=fixed.eps, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, includes.fossils=includes.fossils, fix.type=fix.type, scale.int=scale.int, min.number.points=0)
        results <- rbind(results, new.results[-1,])
    }
    return(results)
}


print.misse.support <- function(x,...){
    
    cat("\nSupport Region\n")
    print(x$ci[,-1])
    cat("\n")
}


