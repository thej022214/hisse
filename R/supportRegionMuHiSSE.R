######################################################################################################################################
######################################################################################################################################
### Adaptive Bootstrap MuHiSSE -- Simulating confidence intervals for parameters estimated in MuHiSSE
######################################################################################################################################
######################################################################################################################################

SupportRegionMuHiSSE <- function(muhisse.obj, n.points=1000, scale.int=0.1, desired.delta=2, min.number.points=10, verbose=TRUE){
    phy <- muhisse.obj$phy
    data <- muhisse.obj$data
    data.new <- data.frame(data[,2], data[,3], row.names=data[,1])
    data.new <- data.new[phy$tip.label,]

    f <- muhisse.obj$f
    np <- max(muhisse.obj$index.par) - 1
    par <- numeric(np)
    free.parameters <- which(muhisse.obj$index.par < max(muhisse.obj$index.par))
    np.sequence <- 1:np
    
    for(i in np.sequence){
        par[i] <- muhisse.obj$solution[which(muhisse.obj$index.par == np.sequence[i])[1]]
    }
    hidden.states=muhisse.obj$hidden.states
    condition.on.survival=muhisse.obj$condition.on.survival
    root.type=muhisse.obj$root.type
    root.p=muhisse.obj$root.p
    
    lower <- exp(muhisse.obj$lower.bounds)
    upper <- exp(muhisse.obj$upper.bounds)
    
    #Bad Jeremy! Hard-coded column headers...
    interval.names <- c("lnLik", "turnover00A","turnover01A","turnover10A","turnover11A","eps00A","eps01A","eps10A","eps11A","q00A_01A","q00A_10A","q00A_11A","q01A_00A","q01A_10A","q01A_11A","q10A_00A","q10A_01A","q10A_11A","q11A_00A","q11A_01A","q11A_10A","q00A_00B","q00A_00C","q00A_00D","q00A_00E","q00A_00F","q00A_00G","q00A_00H","q01A_01B","q01A_01C","q01A_01D","q01A_01E","q01A_01F","q01A_01G","q01A_01H","q10A_10B","q10A_10C","q10A_10D","q10A_10E","q10A_10F","q10A_10G","q10A_10H","q11A_11B","q11A_11C","q11A_11D","q11A_11E","q11A_11F","q11A_11G","q11A_11H","turnover00B","turnover01B","turnover10B","turnover11B","eps00B","eps01B","eps10B","eps11B","q00B_01B","q00B_10B","q00B_11B","q01B_00B","q01B_10B","q01B_11B","q10B_00B","q10B_01B","q10B_11B","q11B_00B","q11B_01B","q11B_10B","q00B_00A","q00B_00C","q00B_00D","q00B_00E","q00B_00F","q00B_00G","q00B_00H","q01B_01A","q01B_01C","q01B_01D","q01B_01E","q01B_01F","q01B_01G","q01B_01H","q10B_10A","q10B_10C","q10B_10D","q10B_10E","q10B_10F","q10B_10G","q10B_10H","q11B_11A","q11B_11C","q11B_11D","q11B_11E","q11B_11F","q11B_11G","q11B_11H","turnover00C","turnover01C","turnover10C","turnover11C","eps00C","eps01C","eps10C","eps11C","q00C_01C","q00C_10C","q00C_11C","q01C_00C","q01C_10C","q01C_11C","q10C_00C","q10C_01C","q10C_11C","q11C_00C","q11C_01C","q11C_10C","q00C_00A","q00C_00B","q00C_00D","q00C_00E","q00C_00F","q00C_00G","q00C_00H","q01C_01A","q01C_01B","q01C_01D","q01C_01E","q01C_01F","q01C_01G","q01C_01H","q10C_10A","q10C_10B","q10C_10D","q10C_10E","q10C_10F","q10C_10G","q10C_10H","q11C_11A","q11C_11B","q11C_11D","q11C_11E","q11C_11F","q11C_11G","q11C_11H","turnover00D","turnover01D","turnover10D","turnover11D","eps00D","eps01D","eps10D","eps11D","q00D_01D","q00D_10D","q00D_11D","q01D_00D","q01D_10D","q01D_11D","q10D_00D","q10D_01D","q10D_11D","q11D_00D","q11D_01D","q11D_10D","q00D_00A","q00D_00B","q00D_00C","q00D_00E","q00D_00F","q00D_00G","q00D_00H","q01D_01A","q01D_01B","q01D_01C","q01D_01E","q01D_01F","q01D_01G","q01D_01H","q10D_10A","q10D_10B","q10D_10C","q10D_10E","q10D_10F","q10D_10G","q10D_10H","q11D_11A","q11D_11B","q11D_11C","q11D_11E","q11D_11F","q11D_11G","q11D_11H","turnover00E","turnover01E","turnover10E","turnover11E","eps00E","eps01E","eps10E","eps11E","q00E_01E","q00E_10E","q00E_11E","q01E_00E","q01E_10E","q01E_11E","q10E_00E","q10E_01E","q10E_11E","q11E_00E","q11E_01E","q11E_10E","q00E_00A","q00E_00B","q00E_00C","q00E_00D","q00E_00F","q00E_00G","q00E_00H","q01E_01A","q01E_01B","q01E_01C","q01E_01D","q01E_01F","q01E_01G","q01E_01H","q10E_10A","q10E_10B","q10E_10C","q10E_10D","q10E_10F","q10E_10G","q10E_10H","q11E_11A","q11E_11B","q11E_11C","q11E_11D","q11E_11F","q11E_11G","q11E_11H","turnover00F","turnover01F","turnover10F","turnover11F","eps00F","eps01F","eps10F","eps11F","q00F_01F","q00F_10F","q00F_11F","q01F_00F","q01F_10F","q01F_11F","q10F_00F","q10F_01F","q10F_11F","q11F_00F","q11F_01F","q11F_10F","q00F_00A","q00F_00B","q00F_00C","q00F_00D","q00F_00E","q00F_00G","q00F_00H","q01F_01A","q01F_01B","q01F_01C","q01F_01D","q01F_01E","q01F_01G","q01F_01H","q10F_10A","q10F_10B","q10F_10C","q10F_10D","q10F_10E","q10F_10G","q10F_10H","q11F_11A","q11F_11B","q11F_11C","q11F_11D","q11F_11E","q11F_11G","q11F_11H","turnover00G","turnover01G","turnover10G","turnover11G","eps00G","eps01G","eps10G","eps11G","q00G_01G","q00G_10G","q00G_11G","q01G_00G","q01G_10G","q01G_11G","q10G_00G","q10G_01G","q10G_11G","q11G_00G","q11G_01G","q11G_10G","q00G_00A","q00G_00B","q00G_00C","q00G_00D","q00G_00E","q00G_00F","q00G_00H","q01G_01A","q01G_01B","q01G_01C","q01G_01D","q01G_01E","q01G_01F","q01G_01H","q10G_10A","q10G_10B","q10G_10C","q10G_10D","q10G_10E","q10G_10F","q10G_10H","q11G_11A","q11G_11B","q11G_11C","q11G_11D","q11G_11E","q11G_11F","q11G_11H","turnover00H","turnover01H","turnover10H","turnover11H","eps00H","eps01H","eps10H","eps11H","q00H_01H","q00H_10H","q00H_11H","q01H_00H","q01H_10H","q01H_11H","q10H_00H","q10H_01H","q10H_11H","q11H_00H","q11H_01H","q11H_10H","q00H_00A","q00H_00B","q00H_00C","q00H_00D","q00H_00E","q00H_00F","q00H_00G","q01H_01A","q01H_01B","q01H_01C","q01H_01D","q01H_01E","q01H_01F","q01H_01G","q10H_10A","q10H_10B","q10H_10C","q10H_10D","q10H_10E","q10H_10F","q10H_10G","q11H_11A","q11H_11B","q11H_11C","q11H_11D","q11H_11E","q11H_11F","q11H_11G")

    interval.results <- AdaptiveConfidenceIntervalSamplingMuHiSSE(par, lower=lower, upper=upper, desired.delta = desired.delta, n.points=n.points, verbose=verbose, phy=phy, data=data.new, index.par=muhisse.obj$index.par, f=f, hidden.states=hidden.states, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, scale.int=scale.int, min.number.points=min.number.points)
    interval.results.final <- matrix(0, n.points+1, length(muhisse.obj$index.par))
    for(i in 1:(n.points+1)){
        par.rep <- unlist(interval.results[i,-1],use.names=FALSE)
        interval.results.final[i,] <- c(par.rep,0)[muhisse.obj$index.par]
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
        class(obj) = "muhisse.support"
        return(obj)
    }
}



AdaptiveConfidenceIntervalSamplingMuHiSSE <- function(par, lower, upper, desired.delta=2, n.points=5000, verbose=TRUE, phy, data, index.par, f, hidden.states, condition.on.survival, root.type, root.p, scale.int, min.number.points=10) {
    
    # Some new prerequisites #
    gen <- FindGenerations(phy)
    dat.tab <- OrganizeData(data=data, phy=phy, f=f, hidden.states=hidden.states)
    ##########################

    #Wrangle the data so that we can make use of DownPass easily:
    actual.params = which(index.par < max(index.par))
    model.vec <- numeric(length(index.par))
    model.vec[] <- c(par,0)[index.par]
    cache <- ParametersToPassMuHiSSE(model.vec=model.vec, hidden.states=hidden.states, nb.tip=Ntip(phy), nb.node=Nnode(phy), bad.likelihood=exp(-300), ode.eps=0)
    phy$node.label <- NULL
    starting <- -DownPassMuHisse(dat.tab=dat.tab, gen=gen, cache=cache, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p)
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
        cache = ParametersToPassMuHiSSE(model.vec=model.vec, hidden.states=hidden.states, nb.tip=Ntip(phy), nb.node=Nnode(phy), bad.likelihood=exp(-300), ode.eps=0)
        phy$node.label <- NULL
        second <- -DownPassMuHisse(dat.tab=dat.tab, gen=gen, cache=cache, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p)
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
        new.results <- AdaptiveConfidenceIntervalSamplingMuHiSSE(par, lower=lower, upper=upper, desired.delta = desired.delta, n.points=n.points, verbose=verbose, phy=phy, data=data, index.par=index.par, f=f, hidden.states=hidden.states, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, scale.int=scale.int, min.number.points=0)
        results <- rbind(results, new.results[-1,])
    }
    return(results)
}


print.muhisse.support <- function(x,...){
    
    cat("\nSupport Region\n")
    print(x$ci[,-1])
    cat("\n")
}


