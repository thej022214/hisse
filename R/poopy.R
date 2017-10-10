
######################################################################################################################################
######################################################################################################################################
### Adaptive Bootstrap -- Simulating confidence intervals for parameters estimated in HiSSE
######################################################################################################################################
######################################################################################################################################

SupportRegion <- function(hisse.obj, n.points=1000, scale.int=0.1, desired.delta=2, output.type="turnover", hidden.states=TRUE, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, verbose=TRUE){
    if(class(hisse.obj) == "hisse.null4.fit"){
        phy <- hisse.obj$phy
        data <- hisse.obj$data
        data.new<-data.frame(data[,2], data[,2], row.names=data[,1])
        data.new<-data.new[phy$tip.label,]
        f <- hisse.obj$f
        np <- max(hisse.obj$index.par)
        par <- numeric(np)
        free.parameters <- which(hisse.obj$index.par < max(hisse.obj$index.par))
        np.sequence <- 1:np
        for(i in np.sequence){
            par[i] <- hisse.obj$solution[which(hisse.obj$index.par == np.sequence[i])[1]]
        }

        lower <- hisse.obj$lower.bounds
        upper <- hisse.obj$upper.bounds

        #Bad Jeremy! Hard-coded column headers...
        if(output.type == "turnover"){
            interval.names <- c("lnLik", "turn.0A", "turn.0B", "turn.0C", "turn.0D", "turn.1A", "turn.1B", "turn.1C", "turn.1D", "eps.0A", "eps.0B", "eps.0C", "eps.0D", "eps.1A", "eps.1B", "eps.1C", "eps.1D", "q0B0A", "q0C0A", "q0D0A", "q1A0A", "q0A0B", "q0C0B", "q0D0B", "q1B0B", "q0A0C", "q0B0C", "q0D0C", "q1C0C", "q0A0D", "q0B0D", "q0C0D", "q1D0D", "q0A1A", "q1B1A", "q1C1A", "q1D1A", "q0B1B", "q1A1B", "q1C1B", "q1D1B", "q0C1C", "q1A1C", "q1B1C", "q1D1C", "q0D1D", "q1A1D", "q1B1D", "q1C1D")
        }
        if(output.type == "net.div"){
            interval.names <- c("lnLik", "netdiv.0A", "netdiv.0B", "netdiv.0C", "netdiv.0D", "netdiv.1A", "netdiv.1B", "netdiv.1C", "netdiv.1D", "eps.0A", "eps.0B", "eps.0C", "eps.0D", "eps.1A", "eps.1B", "eps.1C", "eps.1D", "q0B0A", "q0C0A", "q0D0A", "q1A0A", "q0A0B", "q0C0B", "q0D0B", "q1B0B", "q0A0C", "q0B0C", "q0D0C", "q1C0C", "q0A0D", "q0B0D", "q0C0D", "q1D0D", "q0A1A", "q1B1A", "q1C1A", "q1D1A", "q0B1B", "q1A1B", "q1C1B", "q1D1B", "q0C1C", "q1A1C", "q1B1C", "q1D1C", "q0D1D", "q1A1D", "q1B1D", "q1C1D")
        }
        if(output.type == "raw"){
            interval.names <- c("lnLik", "lambda.0A", "lambda.0B", "lambda.0C", "lambda.0D", "lambda.1A", "lambda.1B", "lambda.1C", "lambda.1D", "mu.0A", "mu.0B", "mu.0C", "mu.0D", "mu.1A", "mu.1B", "mu.1C", "mu.1D", "q0B0A", "q0C0A", "q0D0A", "q1A0A", "q0A0B", "q0C0B", "q0D0B", "q1B0B", "q0A0C", "q0B0C", "q0D0C", "q1C0C", "q0A0D", "q0B0D", "q0C0D", "q1D0D", "q0A1A", "q1B1A", "q1C1A", "q1D1A", "q0B1B", "q1A1B", "q1C1B", "q1D1B", "q0C1C", "q1A1C", "q1B1C", "q1D1C", "q0D1D", "q1A1D", "q1B1D", "q1C1D")
        }

        interval.results <- AdaptiveConfidenceIntervalSampling(par, lower=lower, upper=upper, desired.delta = desired.delta, n.points=n.points, verbose=verbose, phy=phy, data=data.new, index.par=hisse.obj$index.par, f=f, hidden.states=hidden.states, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, scale.int=scale.int, hisse.null.four=TRUE)
        interval.results.final <- matrix(0, n.points+1, length(hisse.obj$index.par))
        for(i in 1:(n.points+1)){
            par.rep <- unlist(interval.results[i,-1],use.names=FALSE)
            interval.results.final[i,] <- c(par.rep,0)[hisse.obj$index.par]
        }
        interval.results.final <- cbind(interval.results[,1], interval.results.final)
        if(output.type == "net.div"){
            lambda.0A <- interval.results.final[,2] / (1 + interval.results.final[,10])
            lambda.0B <- interval.results.final[,3] / (1 + interval.results.final[,11])
            lambda.0C <- interval.results.final[,4] / (1 + interval.results.final[,12])
            lambda.0D <- interval.results.final[,5] / (1 + interval.results.final[,13])
            lambda.1A <- interval.results.final[,6] / (1 + interval.results.final[,14])
            lambda.1B <- interval.results.final[,7] / (1 + interval.results.final[,15])
            lambda.1C <- interval.results.final[,8] / (1 + interval.results.final[,16])
            lambda.1D <- interval.results.final[,9] / (1 + interval.results.final[,17])
            mu.0A <- (interval.results.final[,2] * interval.results.final[,10]) / (1 + interval.results.final[,10])
            mu.0B <- (interval.results.final[,3] * interval.results.final[,11]) / (1 + interval.results.final[,11])
            mu.0C <- (interval.results.final[,4] * interval.results.final[,12]) / (1 + interval.results.final[,12])
            mu.0D <- (interval.results.final[,5] * interval.results.final[,13]) / (1 + interval.results.final[,13])
            mu.1A <- (interval.results.final[,6] * interval.results.final[,14]) / (1 + interval.results.final[,14])
            mu.1B <- (interval.results.final[,7] * interval.results.final[,15]) / (1 + interval.results.final[,15])
            mu.1C <- (interval.results.final[,8] * interval.results.final[,16]) / (1 + interval.results.final[,16])
            mu.1D <- (interval.results.final[,9] * interval.results.final[,17]) / (1 + interval.results.final[,17])
            interval.results.final[,2] <- lambda.0A - mu.0A
            interval.results.final[,3] <- lambda.0B - mu.0B
            interval.results.final[,4] <- lambda.0C - mu.0C
            interval.results.final[,5] <- lambda.0D - mu.0D
            interval.results.final[,6] <- lambda.1A - mu.1A
            interval.results.final[,7] <- lambda.1B - mu.1B
            interval.results.final[,8] <- lambda.1C - mu.1C
            interval.results.final[,9] <- lambda.1D - mu.1D
        }
        if(output.type == "raw"){
            lambda.0A <- interval.results.final[,2] / (1 + interval.results.final[,10])
            lambda.0B <- interval.results.final[,3] / (1 + interval.results.final[,11])
            lambda.0C <- interval.results.final[,4] / (1 + interval.results.final[,12])
            lambda.0D <- interval.results.final[,5] / (1 + interval.results.final[,13])
            lambda.1A <- interval.results.final[,6] / (1 + interval.results.final[,14])
            lambda.1B <- interval.results.final[,7] / (1 + interval.results.final[,15])
            lambda.1C <- interval.results.final[,8] / (1 + interval.results.final[,16])
            lambda.1D <- interval.results.final[,9] / (1 + interval.results.final[,17])
            mu.0A <- (interval.results.final[,2] * interval.results.final[,10]) / (1 + interval.results.final[,10])
            mu.0B <- (interval.results.final[,3] * interval.results.final[,11]) / (1 + interval.results.final[,11])
            mu.0C <- (interval.results.final[,4] * interval.results.final[,12]) / (1 + interval.results.final[,12])
            mu.0D <- (interval.results.final[,5] * interval.results.final[,13]) / (1 + interval.results.final[,13])
            mu.1A <- (interval.results.final[,6] * interval.results.final[,14]) / (1 + interval.results.final[,14])
            mu.1B <- (interval.results.final[,7] * interval.results.final[,15]) / (1 + interval.results.final[,15])
            mu.1C <- (interval.results.final[,8] * interval.results.final[,16]) / (1 + interval.results.final[,16])
            mu.1D <- (interval.results.final[,9] * interval.results.final[,17]) / (1 + interval.results.final[,17])
            interval.results.final[,2] <- lambda.0A
            interval.results.final[,3] <- lambda.0B
            interval.results.final[,4] <- lambda.0C
            interval.results.final[,5] <- lambda.0D
            interval.results.final[,6] <- lambda.1A
            interval.results.final[,7] <- lambda.1B
            interval.results.final[,8] <- lambda.1C
            interval.results.final[,9] <- lambda.1D
            interval.results.final[,10] <- mu.0A
            interval.results.final[,11] <- mu.0B
            interval.results.final[,12] <- mu.0C
            interval.results.final[,13] <- mu.0D
            interval.results.final[,14] <- mu.1A
            interval.results.final[,15] <- mu.1B
            interval.results.final[,16] <- mu.1C
            interval.results.final[,17] <- mu.1D
        }
        interval.results.in <- interval.results.final[which(interval.results.final[,1] - min(interval.results.final[,1])<=desired.delta),]
        ci.interval = apply(interval.results.in, 2, quantile)
        colnames(interval.results.final) <- colnames(interval.results.in) <- colnames(ci.interval) <- interval.names
        obj = NULL
        obj$ci <- ci.interval
        obj$points.within.region = interval.results.in
        obj$all.points = interval.results.final
    }else{
        phy <- hisse.obj$phy
        data <- hisse.obj$data
        data.new<-data.frame(data[,2], data[,2], row.names=data[,1])
        data.new<-data.new[phy$tip.label,]
        f <- hisse.obj$f
        np <- max(hisse.obj$index.par)-1
        par <- numeric(np)
        free.parameters <- which(hisse.obj$index.par < max(hisse.obj$index.par))
        np.sequence <- 1:np
        for(i in np.sequence){
            par[i] <- hisse.obj$solution[which(hisse.obj$index.par == np.sequence[i])[1]]
        }

        lower <- hisse.obj$lower.bounds
        upper <- hisse.obj$upper.bounds

        #Bad Jeremy! Hard-coded column headers...
        if(output.type == "turnover"){
            interval.names <- c("lnLik", "turn.0A", "turn.1A", "turn.0B", "turn.1B", "eps.0A", "eps.1A", "eps.0B", "eps.1B","q1A0A","q0B0A","q1B0A","q0A1A","q0B1A","q1B1A","q0A0B","q1A0B","q1B0B","q0A1B","q1A1B","q0A1B","turn.alpha.0A","turn.alpha.1A", "turn.alpha.0B", "turn.alpha.1B", "turn.beta.0A","turn.beta.1A", "turn.beta.0B", "turn.beta.1B", "eps.alpha.0A","eps.alpha.1A", "eps.alpha.0B", "eps.alpha.1B", "eps.beta.0A","eps.beta.1A", "eps.beta.0B", "eps.beta.1B", "turn.slice.0A","turn.slice.1A", "turn.slice.0B", "turn.slice.1B", "eps.slice.0A","eps.slice.1A", "eps.slice.0B", "eps.slice.1B", "q0A1A.slice","q1A0A.slice","q0A0B.slice","q0B0A.slice","q1A1B.slice","q1B1A.slice","q0A1B.slice","q1B0A.slice","q1A0B.slice","q0B1A.slice","q1B0B.slice","q0B1B.slice")
        }
        if(output.type == "net.div"){
            interval.names <- c("lnLik", "netdiv.0A", "netdiv.1A", "netdiv.0B", "netdiv.1B", "eps.0A", "eps.1A", "eps.0B", "eps.1B","q1A0A","q0B0A","q1B0A","q0A1A","q0B1A","q1B1A","q0A0B","q1A0B","q1B0B","q0A1B","q1A1B","q0A1B","turn.alpha.0A","turn.alpha.1A", "turn.alpha.0B", "turn.alpha.1B", "turn.beta.0A","turn.beta.1A", "turn.beta.0B", "turn.beta.1B", "eps.alpha.0A","eps.alpha.1A", "eps.alpha.0B", "eps.alpha.1B", "eps.beta.0A","eps.beta.1A", "eps.beta.0B", "eps.beta.1B", "turn.slice.0A","turn.slice.1A", "turn.slice.0B", "turn.slice.1B", "eps.slice.0A","eps.slice.1A", "eps.slice.0B", "eps.slice.1B", "q0A1A.slice","q1A0A.slice","q0A0B.slice","q0B0A.slice","q1A1B.slice","q1B1A.slice","q0A1B.slice","q1B0A.slice","q1A0B.slice","q0B1A.slice","q1B0B.slice","q0B1B.slice")
        }
        if(output.type == "raw"){
            interval.names <- c("lnLik", "lambda.0A", "lambda.1A", "lambda.0B", "lambda.1B", "mu.0A", "mu.1A", "mu.0B", "mu.1B","q1A0A","q0B0A","q1B0A","q0A1A","q0B1A","q1B1A","q0A0B","q1A0B","q1B0B","q0A1B","q1A1B","q0A1B","turn.alpha.0A","turn.alpha.1A", "turn.alpha.0B", "turn.alpha.1B", "turn.beta.0A","turn.beta.1A", "turn.beta.0B", "turn.beta.1B", "eps.alpha.0A","eps.alpha.1A", "eps.alpha.0B", "eps.alpha.1B", "eps.beta.0A","eps.beta.1A", "eps.beta.0B", "eps.beta.1B", "turn.slice.0A","turn.slice.1A", "turn.slice.0B", "turn.slice.1B", "eps.slice.0A","eps.slice.1A", "eps.slice.0B", "eps.slice.1B", "q0A1A.slice","q1A0A.slice","q0A0B.slice","q0B0A.slice","q1A1B.slice","q1B1A.slice","q0A1B.slice","q1B0A.slice","q1A0B.slice","q0B1A.slice","q1B0B.slice","q0B1B.slice")
        }

        interval.results <- AdaptiveConfidenceIntervalSampling(par, lower=lower, upper=upper, desired.delta = desired.delta, n.points=n.points, verbose=verbose, phy=phy, data=data.new, index.par=hisse.obj$index.par, f=f, hidden.states=hidden.states, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, scale.int=scale.int, hisse.null.four=FALSE)
        interval.results.final <- matrix(0, n.points+1, length(hisse.obj$index.par))
        for(i in 1:(n.points+1)){
            par.rep <- unlist(interval.results[i,-1],use.names=FALSE)
            interval.results.final[i,] <- c(par.rep,0)[hisse.obj$index.par]
        }
        interval.results.final[,21:56] = 1
        interval.results.final <- cbind(interval.results[,1], interval.results.final)
        if(output.type == "net.div"){
            lambda.0A <- interval.results.final[,2] / (1 + interval.results.final[,6])
            lambda.1A <- interval.results.final[,3] / (1 + interval.results.final[,7])
            lambda.0B <- interval.results.final[,4] / (1 + interval.results.final[,8])
            lambda.1B <- interval.results.final[,5] / (1 + interval.results.final[,9])
            mu.0A <- (interval.results.final[,2] * interval.results.final[,6]) / (1 + interval.results.final[,6])
            mu.1A <- (interval.results.final[,3] * interval.results.final[,7]) / (1 + interval.results.final[,7])
            mu.0B <- (interval.results.final[,4] * interval.results.final[,8]) / (1 + interval.results.final[,8])
            mu.1B <- (interval.results.final[,5] * interval.results.final[,9]) / (1 + interval.results.final[,9])
            interval.results.final[,2] <- lambda.0A - mu.0A
            interval.results.final[,3] <- lambda.1A - mu.1A
            interval.results.final[,4] <- lambda.0B - mu.0B
            interval.results.final[,5] <- lambda.1B - mu.1B
        }
        if(output.type == "raw"){
            lambda.0A <- interval.results.final[,2] / (1 + interval.results.final[,6])
            lambda.1A <- interval.results.final[,3] / (1 + interval.results.final[,7])
            lambda.0B <- interval.results.final[,4] / (1 + interval.results.final[,8])
            lambda.1B <- interval.results.final[,5] / (1 + interval.results.final[,9])
            mu.0A <- (interval.results.final[,2] * interval.results.final[,6]) / (1 + interval.results.final[,6])
            mu.1A <- (interval.results.final[,3] * interval.results.final[,7]) / (1 + interval.results.final[,7])
            mu.0B <- (interval.results.final[,4] * interval.results.final[,8]) / (1 + interval.results.final[,8])
            mu.1B <- (interval.results.final[,5] * interval.results.final[,9]) / (1 + interval.results.final[,9])
            interval.results.final[,2] <- lambda.0A
            interval.results.final[,3] <- lambda.1A
            interval.results.final[,4] <- lambda.0B
            interval.results.final[,5] <- lambda.1B
            interval.results.final[,6] <- mu.0A
            interval.results.final[,7] <- mu.1A
            interval.results.final[,8] <- mu.0B
            interval.results.final[,9] <- mu.1B
        }
        interval.results.in <- interval.results.final[which(interval.results.final[,1] - min(interval.results.final[,1])<=desired.delta),]
        ci.interval = apply(interval.results.in, 2, quantile)
        colnames(interval.results.final) <- colnames(interval.results.in) <- colnames(ci.interval) <- interval.names

        obj = NULL
        obj$ci <- ci.interval[,1:21]
        obj$points.within.region = interval.results.in[,1:21]
        obj$all.points = interval.results.final[,1:21]
    }
    class(obj) = "hisse.support"
    return(obj)
}


AdaptiveConfidenceIntervalSampling <- function(par, lower, upper, desired.delta=2, n.points=5000, verbose=TRUE, phy, data, index.par, f, hidden.states, condition.on.survival, root.type, root.p, scale.int, hisse.null.four=FALSE, min.number.points=10) {
    #Wrangle the data so that we can make use of DownPass easily:
    actual.params = which(index.par < max(index.par))
    model.vec <- numeric(length(index.par))
    model.vec[] <- c(par,0)[index.par]
    if(hisse.null.four == TRUE){
        cache = ParametersToPassNull(phy, data[,1], f=f, model.vec)
        phy$node.label <- NULL
        starting <- -DownPassNull(phy, cache, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p)
    }else{
        model.vec.tmp = model.vec[21:56]
        model.vec.tmp[model.vec.tmp==0] = 1
        model.vec[21:56] = model.vec.tmp
        cache = ParametersToPass(phy, data[,1], model.vec, f=f, timeslice=NULL, hidden.states=hidden.states)
        cache$turnover.beta.factor0 = 1 / dbeta(0.1, model.vec[21], model.vec[25])
        cache$turnover.beta.factor1 = 1 / dbeta(0.1, model.vec[22], model.vec[26])
        cache$turnover.beta.factorA = 1 / dbeta(0.1, model.vec[23], model.vec[27])
        cache$turnover.beta.factorB = 1 / dbeta(0.1, model.vec[24], model.vec[28])
        cache$eps.beta.factor0 = 1 / dbeta(0.1, model.vec[29], model.vec[33])
        cache$eps.beta.factor1 = 1 / dbeta(0.1, model.vec[30], model.vec[34])
        cache$eps.beta.factorA = 1 / dbeta(0.1, model.vec[31], model.vec[35])
        cache$eps.beta.factorB = 1 / dbeta(0.1, model.vec[32], model.vec[36])
        phy$node.label <- NULL
        #############################################################
        #Now assess the likelihood at the MLE:
        starting <- -DownPass(phy, cache, hidden.states=hidden.states, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p)
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
        if(hisse.null.four == TRUE){
            cache = ParametersToPassNull(phy, data[,1], f=f, model.vec)
            second <- -DownPassNull(phy, cache, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p)
        }else{
            model.vec.tmp = model.vec[21:56]
            model.vec.tmp[model.vec.tmp==0] = 1
            model.vec[21:56] = model.vec.tmp
            cache = ParametersToPass(phy, data[,1], model.vec, f=f, timeslice=NULL, hidden.states=hidden.states)
            cache$turnover.beta.factor0 = 1 / dbeta(0.1, model.vec[21], model.vec[25])
            cache$turnover.beta.factor1 = 1 / dbeta(0.1, model.vec[22], model.vec[26])
            cache$turnover.beta.factorA = 1 / dbeta(0.1, model.vec[23], model.vec[27])
            cache$turnover.beta.factorB = 1 / dbeta(0.1, model.vec[24], model.vec[28])
            cache$eps.beta.factor0 = 1 / dbeta(0.1, model.vec[29], model.vec[33])
            cache$eps.beta.factor1 = 1 / dbeta(0.1, model.vec[30], model.vec[34])
            cache$eps.beta.factorA = 1 / dbeta(0.1, model.vec[31], model.vec[35])
            cache$eps.beta.factorB = 1 / dbeta(0.1, model.vec[32], model.vec[36])
            second <- -DownPass(phy, cache, hidden.states=hidden.states, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p)
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
      new.results <- AdaptiveConfidenceIntervalSampling(par=par, lower=lower, upper=upper, desired.delta=desired.delta, n.points=2+round(n.points/4), verbose=verbose, phy=phy, data=data, index.par=index.par, f=f, hidden.states=hidden.states, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, scale.int=scale.int, hisse.null.four=hisse.null.four, min.number.points=0)
      results <- rbind(results, new.results[-1,])
    }
    return(results)
}



GenerateValues <- function(par, lower, upper, scale.int, max.tries=100, expand.prob=0, examined.max, examined.min) {
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


print.hisse.support <- function(x,...){

    cat("\nSupport Region\n")
    print(x$ci[,-1])
    cat("\n")
}
