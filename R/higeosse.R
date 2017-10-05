
#library(deSolve)
#library(diversitree)
#library(nloptr)
#dyn.load("../src/higeosse-ext-derivs.so")
#dyn.load("../src/canonical_geosse-ext-derivs.so")


######################################################################################################################################
######################################################################################################################################
### HiGeoSSE -- Expanded set of GeoSSE models for examining diversification in relation to geographic range evolution
######################################################################################################################################
######################################################################################################################################

higeosse <- function(phy, data, f=c(1,1,1), speciation=c(1,2,3), extinction=c(1,2), hidden.states=FALSE, trans.rate=NULL, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, sann=FALSE, sann.its=10000, bounded.search=TRUE, max.tol=.Machine$double.eps^.25, mag.san.start=0.5, starting.vals=NULL, speciation.upper=1000, extinction.upper=1000, trans.upper=100, ode.eps=0){
    
    if(!is.null(root.p)) {
        root.type="user"
        root.p <- root.p / sum(root.p)
        if(hidden.states ==TRUE & length(root.p)==2){
            root.p <- rep(root.p, 2)
            root.p <- root.p / sum(root.p)
            warning("For hidden states, you need to specify the root.p for all four hidden states. We have adjusted it so that there's equal chance for 0A as 0B, and for 1A as 1B")
        }
    }
    
    if(is.null(trans.rate)){
        stop("Rate matrix needed. See TransMatMakerHiGeoSSE() to create one.")
    }
    
    if(hidden.states == TRUE & dim(trans.rate)[1]<4){
        stop("You chose a hidden state but this is not reflected in the transition matrix")
    }

    pars <- numeric(95)

    if(dim(trans.rate)[2]==3){
        rate.cats <- 1
        pars.tmp <- speciation
        extinction.tmp <- extinction
        extinction.tmp[which(extinction.tmp > 0)] = (extinction.tmp[which( extinction.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, extinction.tmp)
        trans.tmp <- c(trans.rate["(0)", "(01)"], trans.rate["(1)", "(01)"])
        trans.tmp[which(trans.tmp > 0)] = (trans.tmp[which(trans.tmp > 0)] + max(pars.tmp))
        category.rates.unique <- 0
        pars.tmp <- c(pars.tmp, trans.tmp)
        pars[1:7] <- pars.tmp
    }
    
    if(dim(trans.rate)[2]==6){
        rate.cats <- 2
        pars.tmp <- speciation
        extinction.tmp <- extinction
        extinction.tmp[which(extinction.tmp > 0)] = (extinction.tmp[which( extinction.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, extinction.tmp)
        rows <- c("(0A)", "(1A)", "(0B)", "(1B)")
        cols <- c("(01A)", "(01A)", "(01B)", "(01B)")
        trans.tmp <- trans.rate[cbind(rows,cols)]
        trans.tmp[which(trans.tmp > 0)] = (trans.tmp[which(trans.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, trans.tmp)
        category.tmp <- trans.rate[which(trans.rate==max(trans.rate, na.rm=TRUE))]
        category.rate.shift <- rep(max(pars.tmp)+1, length(category.tmp))
        category.rate.shiftA <- c(category.rate.shift[1], rep(0,3), category.rate.shift[2], rep(0,3), category.rate.shift[3], rep(0,3))
        category.rate.shiftB <- c(category.rate.shift[4], rep(0,3), category.rate.shift[5], rep(0,3), category.rate.shift[6], rep(0,3))
        pars.tmp <- c(speciation[1:3], extinction.tmp[1:2], trans.tmp[1:2], category.rate.shiftA, speciation[4:6], extinction.tmp[3:4], trans.tmp[3:4], category.rate.shiftB)
        pars[1:length(pars.tmp)] <- pars.tmp
    }
    
    if(dim(trans.rate)[2]==9){
        rate.cats <- 3
        pars.tmp <- speciation
        extinction.tmp <- extinction
        extinction.tmp[which(extinction.tmp > 0)] = (extinction.tmp[which( extinction.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, extinction.tmp)
        rows <- c("(0A)", "(1A)", "(0B)", "(1B)", "(0C)", "(1C)")
        cols <- c("(01A)", "(01A)", "(01B)", "(01B)", "(01C)", "(01C)")
        trans.tmp <- trans.rate[cbind(rows,cols)]
        trans.tmp[which(trans.tmp > 0)] = (trans.tmp[which(trans.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, trans.tmp)
        category.tmp <- trans.rate[which(trans.rate==max(trans.rate, na.rm=TRUE))]
        category.rate.shift <- rep(max(pars.tmp)+1, length(category.tmp))
        category.rate.shiftA <- c(category.rate.shift[1:2], rep(0,2), category.rate.shift[3:4], rep(0,2), category.rate.shift[5:6], rep(0,2))
        category.rate.shiftB <- c(category.rate.shift[7:8], rep(0,2), category.rate.shift[9:10], rep(0,2), category.rate.shift[11:12], rep(0,2))
        category.rate.shiftC <- c(category.rate.shift[13:14], rep(0,2), category.rate.shift[15:16], rep(0,2), category.rate.shift[17:18], rep(0,2))
        pars.tmp <- c(speciation[1:3], extinction.tmp[1:2], trans.tmp[1:2], category.rate.shiftA, speciation[4:6], extinction.tmp[3:4], trans.tmp[3:4], category.rate.shiftB, speciation[7:9], extinction.tmp[5:6], trans.tmp[5:6], category.rate.shiftC)
        pars[1:length(pars.tmp)] <- pars.tmp
    }
    
    if(dim(trans.rate)[2]==12){
        rate.cats <- 4
        pars.tmp <- speciation
        extinction.tmp <- extinction
        extinction.tmp[which(extinction.tmp > 0)] = (extinction.tmp[which( extinction.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, extinction.tmp)
        rows <- c("(0A)", "(1A)", "(0B)", "(1B)", "(0C)", "(1C)", "(0D)", "(1D)")
        cols <- c("(01A)", "(01A)", "(01B)", "(01B)", "(01C)", "(01C)", "(01D)", "(01D)")
        trans.tmp <- trans.rate[cbind(rows,cols)]
        trans.tmp[which(trans.tmp > 0)] = (trans.tmp[which(trans.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, trans.tmp)
        category.tmp <- trans.rate[which(trans.rate==max(trans.rate, na.rm=TRUE))]
        category.rate.shift <- rep(max(pars.tmp)+1, length(category.tmp))
        category.rate.shiftA <- c(category.rate.shift[1:3], rep(0,1), category.rate.shift[4:6], rep(0,1), category.rate.shift[7:9], rep(0,1))
        category.rate.shiftB <- c(category.rate.shift[10:12], rep(0,1), category.rate.shift[13:15], rep(0,1), category.rate.shift[16:18], rep(0,1))
        category.rate.shiftC <- c(category.rate.shift[19:21], rep(0,1), category.rate.shift[22:24], rep(0,1), category.rate.shift[25:27], rep(0,1))
        category.rate.shiftD <- c(category.rate.shift[28:30], rep(0,1), category.rate.shift[31:33], rep(0,1), category.rate.shift[34:36], rep(0,1))
        category.rates.all <- c(category.rate.shiftA, category.rate.shiftB, category.rate.shiftC, category.rate.shiftD)
        category.rates.unique <- length(unique(category.rates.all[category.rates.all>0]))
        pars.tmp <- c(speciation[1:3], extinction.tmp[1:2], trans.tmp[1:2], category.rate.shiftA, speciation[4:6], extinction.tmp[3:4], trans.tmp[3:4], category.rate.shiftB, speciation[7:9], extinction.tmp[5:6], trans.tmp[5:6], category.rate.shiftC, speciation[10:12], extinction.tmp[7:8], trans.tmp[7:8], category.rate.shiftD)
        pars[1:length(pars.tmp)] <- pars.tmp
    }

    if(dim(trans.rate)[2]==15){
        rate.cats <- 5
        pars.tmp <- speciation
        extinction.tmp <- extinction
        extinction.tmp[which(extinction.tmp > 0)] = (extinction.tmp[which(extinction.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, extinction.tmp)
        rows <- c("(0A)", "(1A)", "(0B)", "(1B)", "(0C)", "(1C)", "(0D)", "(1D)", "(0E)", "(1E)")
        cols <- c("(01A)", "(01A)", "(01B)", "(01B)", "(01C)", "(01C)", "(01D)", "(01D)", "(01E)", "(01E)")
        trans.tmp <- trans.rate[cbind(rows,cols)]
        trans.tmp[which(trans.tmp > 0)] = (trans.tmp[which(trans.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, trans.tmp)
        category.tmp <- trans.rate[which(trans.rate==max(trans.rate, na.rm=TRUE))]
        category.rate.shift <- rep(max(pars.tmp)+1, length(category.tmp))
        category.rate.shiftA <- category.rate.shift[1:12]
        category.rate.shiftB <- category.rate.shift[13:24]
        category.rate.shiftC <- category.rate.shift[25:36]
        category.rate.shiftD <- category.rate.shift[37:48]
        category.rate.shiftE <- category.rate.shift[49:60]
        category.rates.all <- c(category.rate.shiftA, category.rate.shiftB, category.rate.shiftC, category.rate.shiftD, category.rate.shiftE)
        category.rates.unique <- length(unique(category.rates.all[category.rates.all>0]))
        pars.tmp <- c(speciation[1:3], extinction.tmp[1:2], trans.tmp[1:2], category.rate.shiftA, speciation[4:6], extinction.tmp[3:4], trans.tmp[3:4], category.rate.shiftB, speciation[7:9], extinction.tmp[5:6], trans.tmp[5:6], category.rate.shiftC, speciation[10:12], extinction.tmp[7:8], trans.tmp[7:8], category.rate.shiftD, speciation[13:15], extinction.tmp[9:10], trans.tmp[9:10], category.rate.shiftE)
        pars[1:length(pars.tmp)] <- pars.tmp
    }
    
    np <- max(pars)
    pars[pars==0] <- np+1

    cat("Initializing...", "\n")

    data.new <- data.frame(data[,2], data[,2], row.names=data[,1])
    data.new <- data.new[phy$tip.label,]
   
   #This is used to scale starting values to account for sampling:
    if(length(f) == 3){
        samp.freq.tree <- Ntip(phy) / sum(table(data.new[,1]) / f)
    }else{
        if(length(f) == Ntip(phy)){
            samp.freq.tree <- Ntip(phy) / sum(table(data.new[,1]) / mean(f))
        }else{
            stop("The vector of sampling frequencies does not match the number of tips in the tree.")
        }
    }

    if(sum(extinction)==0){
        init.pars <- starting.point.geosse(phy, eps=0, samp.freq.tree=samp.freq.tree)
    }else{
        init.pars <- starting.point.geosse(phy, eps=mag.san.start, samp.freq.tree=samp.freq.tree)
    }
    names(init.pars) <- NULL

    if(is.null(starting.vals)){
        def.set.pars <- rep(c(log(init.pars[1:3]), log(init.pars[4:5]), log(init.pars[6:7]), rep(log(.01), 12)), rate.cats)
    }else{
        def.set.pars <- rep(c(log(starting.vals[1:3]), log(starting.vals[4:5]), log(starting.vals[6:7]), rep(log(0.01), 12)), rate.cats)
    }
    if(bounded.search == TRUE){
        upper.full <- rep(c(rep(log(speciation.upper),3), rep(log(extinction.upper),2), rep(log(trans.upper),2), rep(log(10), 12)), rate.cats)
    }else{
        upper.full <- rep(21,length(def.set.pars))
    }

    np.sequence <- 1:np
    ip <- numeric(np)
    upper <- numeric(np)
    for(i in np.sequence){
        ip[i] <- def.set.pars[which(pars == np.sequence[i])[1]]
        upper[i] <- upper.full[which(pars == np.sequence[i])[1]]
    }
    lower <- rep(-20, length(ip))

    if(sann == FALSE){
        if(bounded.search == TRUE){
            cat("Finished. Beginning bounded subplex routine...", "\n")
            opts <- list("algorithm" = "NLOPT_LN_SBPLX", "maxeval" = 100000, "ftol_rel" = max.tol)
            out = nloptr(x0=ip, eval_f=DevOptimizeHiGeoSSE, ub=upper, lb=lower, opts=opts, pars=pars, phy=phy, data=data.new[,1], f=f, hidden.states=hidden.states, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, np=np, ode.eps=ode.eps)
            solution <- numeric(length(pars))
            solution[] <- c(exp(out$solution), 0)[pars]
            loglik = -out$objective
        }else{
            cat("Finished. Beginning subplex routine...", "\n")
            out = subplex(ip, fn=DevOptimizeHiGeoSSE, control=list(reltol=max.tol, parscale=rep(0.1, length(ip))), pars=pars, phy=phy, data=data.new[,1], f=f, hidden.states=hidden.states, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, np=np, ode.eps=ode.eps)
            solution <- numeric(length(pars))
            solution[] <- c(exp(out$par), 0)[pars]
            loglik = -out$value
        }
    }else{
        cat("Finished. Beginning simulated annealing...", "\n")
        out.sann = GenSA(ip, fn=DevOptimizeHiGeoSSE, lower=lower, upper=upper, control=list(max.call=sann.its), pars=pars, phy=phy, data=data.new[,1], f=f, hidden.states=hidden.states, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, np=np, ode.eps=ode.eps)
        cat("Finished. Refining using subplex routine...", "\n")
        out = nloptr(x0=out.sann$par, eval_f=DevOptimizeHiGeoSSE, ub=upper, lb=lower, opts=opts, pars=pars, phy=phy, data=data.new[,1], f=f, hidden.states=hidden.states, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, np=np, ode.eps=ode.eps)
        solution <- numeric(length(pars))
        solution[] <- c(exp(out$solution), 0)[pars]
        loglik = -out$objective
    }

    cat("Finished. Summarizing results...", "\n")

    obj = list(loglik = loglik, AIC = -2*loglik+2*np, AICc = -2*loglik+(2*np*(Ntip(phy)/(Ntip(phy)-np-1))), solution=solution, index.par=pars, f=f, hidden.states=hidden.states, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, phy=phy, data=data, trans.matrix=trans.rate, max.tol=max.tol, starting.vals=ip, upper.bounds=upper, lower.bounds=lower, ode.eps=ode.eps)
    class(obj) = "higeosse.fit"

    return(obj)

}

######################################################################################################################################
######################################################################################################################################
### The function used to optimize parameters:
######################################################################################################################################
######################################################################################################################################


DevOptimizeHiGeoSSE <- function(p, pars, phy, data, f, hidden.states, condition.on.survival, root.type, root.p, np, ode.eps) {
    #Generates the final vector with the appropriate parameter estimates in the right place:
    p.new <- exp(p)
    print(p.new)
    model.vec <- numeric(length(pars))
    model.vec[] <- c(p.new, 0)[pars]

    cache = ParametersToPassHiGeoSSE(phy=phy, data=data, f=f, model.vec=model.vec, hidden.states=hidden.states)
    
    logl <- DownPassHiGeosse(phy, cache, hidden.states=hidden.states, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, ode.eps=ode.eps)
    return(-logl)
}


#Taken from the GeoSSE code modified to account for sampling -- credit goes to Emma Goldberg:
starting.point.geosse <- function(tree, eps=0.5, samp.freq.tree) {
    if (eps == 0) {
        s <- (log(Ntip(tree)) - log(2)) / max(branching.times(tree))
        s <- s/samp.freq.tree
        x <- 0
        d <- s/10
    } else {
        n <- Ntip(tree)
        r <- ( log( (n/2) * (1 - eps*eps) + 2*eps + (1 - eps)/2 *
        sqrt( n * (n*eps*eps - 8*eps + 2*n*eps + n))) - log(2)
        ) / max(branching.times(tree))
        s <- r / (1 - eps)
        x <- s * eps
        s <- s/samp.freq.tree
        x <- x - ((s * samp.freq.tree) * (1 - 1/samp.freq.tree))
        d <- x
    }
    p <- c(s, s, s, x, x, d, d)
    names(p) <- c("sA",  "sB",  "sAB", "xA" , "xB"  ,"dA"  ,"dB")
    p
}


######################################################################################################################################
######################################################################################################################################
### The HiGeoSSE down pass that carries out the integration and returns the likelihood:
######################################################################################################################################
######################################################################################################################################

DownPassHiGeosse <- function(phy, cache, hidden.states, bad.likelihood=-10000000000, condition.on.survival, root.type, root.p, get.phi=FALSE, node=NULL, state=NULL, ode.eps=0) {
	#Some preliminaries:
	nb.tip <- length(phy$tip.label)
	nb.node <- phy$Nnode
	phy <- reorder(phy, "pruningwise")
	anc <- unique(phy$edge[,1])		
	TIPS <- 1:nb.tip

	if(hidden.states == FALSE){
		compD <- matrix(0, nrow=nb.tip + nb.node, ncol=3)
		compE <- matrix(0, nrow=nb.tip + nb.node, ncol=3)
	}else{
		compD <- matrix(0, nrow=nb.tip + nb.node, ncol=15)
		compE <- matrix(0, nrow=nb.tip + nb.node, ncol=15)
	}
    
	#Initializes the tip sampling and sets internal nodes to be zero:
	ncols = dim(compD)[2]
	if(length(cache$f) == 3){
		for(i in 1:(nb.tip)){
			compD[i,] <- cache$f * cache$states[i,]
			compE[i,] <- rep((1-cache$f), ncols/3)
		}
	}else{
		for(i in 1:(nb.tip)){
			compD[i,] <- cache$f[i] * cache$states[i,]
			compE[i,] <- rep((1-cache$f[i]), ncols/3)
		}
	}
	logcomp <- c()
	#Start the postorder traversal indexing lists by node number: 
	for (i in seq(from = 1, length.out = nb.node)) {
		#A vector of all the internal nodes:
		focal <- anc[i]
		desRows <- which(phy$edge[,1]==focal)
		desNodes <- phy$edge[desRows,2]
		#Note: when the tree has been reordered branching.times are no longer valid. Fortunately, we extract this information in the initial cache setup. Also focal is the rootward node, whereas desNodes represent a vector of all descendant nodes:
		cache$rootward.age <- cache$split.times[which(names(cache$split.times)==focal)]
		
        v <- c()
		phi <- c()
		for (desIndex in sequence(length(desRows))){
			cache$focal.edge.length <- phy$edge.length[desRows[desIndex]]
			cache$tipward.age <- cache$rootward.age - cache$focal.edge.length
			#Strange rounding errors. A tip age should be zero. This ensures that:
			if(cache$tipward.age < .Machine$double.eps^0.5){
				cache$tipward.age = 0
			}
			cache$node.D <- compD[desNodes[desIndex],]
			cache$node.E <- compE[desNodes[desIndex],]
			##Call to lsoda that utilizes C code. Requires a lot of inputs. Note that for now we hardcode the NUMELEMENTS arguments. The reason for this is because with lsoda we can only pass a vector of parameters.
			if(hidden.states == FALSE){
				pars <- list(cache$s0A, cache$s1A, cache$s01A, cache$x0A, cache$x1A, cache$d0A_01A, cache$d1A_01A)
                NUMELEMENTS <- 7 #needed for passing in vector to C
				padded.pars <- rep(0, NUMELEMENTS)
				pars <- c(unlist(pars))
				stopifnot(length(padded.pars)<=NUMELEMENTS)
				padded.pars[sequence(length(pars))]<-pars
				yini <-c(E_0=cache$node.E[1], E_1=cache$node.E[2], E_01=cache$node.E[3], D_N0=cache$node.D[1], D_N1=cache$node.D[2], D_N2=cache$node.D[3])
				times=c(cache$tipward.age, cache$rootward.age)
                
                runSilent <- function() {
                    options(warn = -1)
                    on.exit(options(warn = 0))
                    capture.output(res <- lsoda(yini, times, func = "classe_geosse_equivalent_derivs", padded.pars, initfunc="initmod_geosse", dllname = "hisse", rtol=1e-8, atol=1e-8))
                    #capture.output(res <- lsoda(yini, times, func = "classe_geosse_equivalent_derivs", padded.pars, initfunc="initmod_geosse", dll = "canonical_geosse-ext-derivs", rtol=1e-8, atol=1e-8))
                    res
                }
                prob.subtree.cal.full <- runSilent()
			}else{
				pars <- list(cache$s0A, cache$s1A, cache$s01A, cache$x0A, cache$x1A, cache$d0A_01A, cache$d1A_01A,  cache$d0A_0B, cache$d0A_0C, cache$d0A_0D, cache$d0A_0E, cache$d1A_1B, cache$d1A_1C, cache$d1A_1D, cache$d1A_1E, cache$d01A_01B, cache$d01A_01C, cache$d01A_01D, cache$d01A_01E, cache$s0B, cache$s1B, cache$s01B, cache$x0B, cache$x1B, cache$d0B_01B, cache$d1B_01B, cache$d0B_0A, cache$d0B_0C, cache$d0B_0D, cache$d0B_0E, cache$d1B_1A, cache$d1B_1C, cache$d1B_1D, cache$d1B_1E, cache$d01B_01A, cache$d01B_01C, cache$d01B_01D, cache$d01B_01E,  cache$s0C, cache$s1C, cache$s01C, cache$x0C, cache$x1C, cache$d0C_01C, cache$d1C_01C,  cache$d0C_0A, cache$d0C_0B, cache$d0C_0D, cache$d0C_0E, cache$d1C_1A, cache$d1C_1B, cache$d1C_1D, cache$d1C_1E, cache$d01C_01A, cache$d01C_01B, cache$d01C_01D, cache$d01C_01E,  cache$s0D, cache$s1D, cache$s01D, cache$x0D, cache$x1D, cache$d0D_01D, cache$d1D_01D,  cache$d0D_0A, cache$d0D_0B, cache$d0D_0C, cache$d0D_0E, cache$d1D_1A, cache$d1D_1B, cache$d1D_1C, cache$d1D_1E, cache$d01D_01A, cache$d01D_01B, cache$d01D_01C, cache$d01D_01E,  cache$s0E, cache$s1E, cache$s01E, cache$x0E, cache$x1E, cache$d0E_01E, cache$d1E_01E,  cache$d0E_0A, cache$d0E_0B, cache$d0E_0C, cache$d0E_0D, cache$d1E_1A, cache$d1E_1B, cache$d1E_1C, cache$d1E_1ED, cache$d01E_01A, cache$d01E_01B, cache$d01E_01C, cache$d01E_01D)
				NUMELEMENTS <- 95 #needed for passing in vector to C
				padded.pars <- rep(0, NUMELEMENTS)
				pars <- c(unlist(pars))
				stopifnot(length(padded.pars)<=NUMELEMENTS)
				padded.pars[sequence(length(pars))]<-pars
                yini <- c(E_0A=cache$node.E[1], E_1A=cache$node.E[2], E_01A=cache$node.E[3], E_0B=cache$node.E[4], E_1B=cache$node.E[5], E_01B=cache$node.E[6], E_0C=cache$node.E[7], E_1C=cache$node.E[8], E_01C=cache$node.E[9], E_0D=cache$node.E[10], E_1D=cache$node.E[11], E_01D=cache$node.E[12], E_0E=cache$node.E[13], E_1E=cache$node.E[14], E_01E=cache$node.E[15], D_N0A=cache$node.D[1], D_N1A=cache$node.D[2], D_N01A=cache$node.D[3], D_N0B=cache$node.D[4], D_N1B=cache$node.D[5], D_N01B=cache$node.D[6], D_N0C=cache$node.D[7], D_N1C=cache$node.D[8], D_N01C=cache$node.D[9], D_N0D=cache$node.D[10], D_N1D=cache$node.D[11], D_N01D=cache$node.D[12], D_N0E=cache$node.D[13], D_N1E=cache$node.D[14], D_N01E=cache$node.D[15])
				times=c(cache$tipward.age, cache$rootward.age)
                
                runSilent <- function() {
                    options(warn = -1)
                    on.exit(options(warn = 0))
                    capture.output(res <- lsoda(yini, times, func = "higeosse_derivs", padded.pars, initfunc="initmod_higeosse", dllname = "hisse", rtol=1e-8, atol=1e-8))
                    #capture.output(res <- lsoda(yini, times, func = "higeosse_derivs", padded.pars, initfunc="initmod_higeosse", dll = "higeosse-ext-derivs", rtol=1e-8, atol=1e-8))
                    res
                }
                prob.subtree.cal.full <- runSilent()
			}

			######## THIS CHECKS TO ENSURE THAT THE INTEGRATION WAS SUCCESSFUL ###########
            if(attributes(prob.subtree.cal.full)$istate[1] < 0){
				return(bad.likelihood)
			}else{
				prob.subtree.cal <- prob.subtree.cal.full[-1,-1]
			}
			##############################################################################
			
            if(hidden.states == FALSE){
				if(is.nan(prob.subtree.cal[3]) | is.nan(prob.subtree.cal[4]) | is.nan(prob.subtree.cal[5])){
					return(bad.likelihood)
				}
                #This is default and cannot change, but if we get a negative probability, discard the results:
				if(prob.subtree.cal[4]<0 | prob.subtree.cal[5]<0 | prob.subtree.cal[6]<0){
					return(bad.likelihood)
				}
                #This can be modified at the input, but if the sum of the D's at the end of a branch are less than some value, then discard the results. A little more stringent than diversitree, but with difficult problems, this stabilizes things immensely.
                if(sum(prob.subtree.cal[4:6]) < ode.eps){
                    return(bad.likelihood)
                }

				#Designating phi here because of its relation to Morlon et al (2011) and using "e" would be confusing:
				phi <- c(phi, prob.subtree.cal[1:3])
				v <- rbind(v, prob.subtree.cal[4:6])
			}else{
                if(any(is.nan(prob.subtree.cal[16:30]))){
					return(bad.likelihood)
				}
                #This is default and cannot change, but if we get a negative probability, discard the results:
                if(any(prob.subtree.cal[16:30] < 0)){
                    return(bad.likelihood)
				}
                #This can be modified at the input, but if the sum of the D's at the end of a branch are less than some value, then discard the results. A little more stringent than diversitree, but with difficult problems, this stabilizes things immensely.
                if(sum(prob.subtree.cal[10:18]) < ode.eps){
                    return(bad.likelihood)
                }
                
				#Designating phi here because of its relation to Morlon et al (2011) and using "e" would be confusing:
                phi <- c(phi, prob.subtree.cal[1:15])
                v <- rbind(v, prob.subtree.cal[16:30])
			}
		}
		if(hidden.states == TRUE){
			compD[focal,1] <- v[1,1] * v[2,1] * cache$s0A
			compD[focal,2] <- v[1,2] * v[2,2] * cache$s1A
            compD[focal,3] <- 0.5 * (v[1,3] * v[2,1] + v[1,1] * v[2,3]) * cache$s0A + 0.5 * (v[1,3] * v[2,2] + v[1,2] * v[2,3]) * cache$s1A + 0.5 * (v[1,1] * v[2,2] + v[1,2] * v[2,1]) * cache$s01A
            v <- v[,-c(1:3)]
            compD[focal,4] <- v[1,1] * v[2,1] * cache$s0B
            compD[focal,5] <- v[1,2] * v[2,2] * cache$s1B
            compD[focal,6] <- 0.5 * (v[1,3] * v[2,1] + v[1,1] * v[2,3]) * cache$s0B + 0.5 * (v[1,3] * v[2,2] + v[1,2] * v[2,3]) * cache$s1B + 0.5 * (v[1,1] * v[2,2] + v[1,2] * v[2,1]) * cache$s01B
            v <- v[,-c(1:3)]
            compD[focal,7] <- v[1,1] * v[2,1] * cache$s0C
            compD[focal,8] <- v[1,2] * v[2,2] * cache$s1C
            compD[focal,9] <- 0.5 * (v[1,3] * v[2,1] + v[1,1] * v[2,3]) * cache$s0C + 0.5 * (v[1,3] * v[2,2] + v[1,2] * v[2,3]) * cache$s1C + 0.5 * (v[1,1] * v[2,2] + v[1,2] * v[2,1]) * cache$s01C
            v <- v[,-c(1:3)]
            compD[focal,10] <- v[1,1] * v[2,1] * cache$s0D
            compD[focal,11] <- v[1,2] * v[2,2] * cache$s1D
            compD[focal,12] <- 0.5 * (v[1,3] * v[2,1] + v[1,1] * v[2,3]) * cache$s0D + 0.5 * (v[1,3] * v[2,2] + v[1,2] * v[2,3]) * cache$s1D + 0.5 * (v[1,1] * v[2,2] + v[1,2] * v[2,1]) * cache$s01D
            v <- v[,-c(1:3)]
            compD[focal,13] <- v[1,1] * v[2,1] * cache$s0E
            compD[focal,14] <- v[1,2] * v[2,2] * cache$s1E
            compD[focal,15] <- 0.5 * (v[1,3] * v[2,1] + v[1,1] * v[2,3]) * cache$s0E + 0.5 * (v[1,3] * v[2,2] + v[1,2] * v[2,3]) * cache$s1E + 0.5 * (v[1,1] * v[2,2] + v[1,2] * v[2,1]) * cache$s01E
            compE[focal,] <- phi[1:15]
			if(!is.null(node)){
				fixer = numeric(15)
				fixer[state] = 1
				if(node == focal){
					compD[focal,] <- compD[focal,] * fixer
				}
				#compE[focal,] <- compE[focal,] * fixer
			}
		}else{
            compD[focal,1] <- v[1,1] * v[2,1] * cache$s0A
            compD[focal,2] <- v[1,2] * v[2,2] * cache$s1A
            compD[focal,3] <- 0.5 * (v[1,3] * v[2,1] + v[1,1] * v[2,3]) * cache$s0A + 0.5 * (v[1,3] * v[2,2] + v[1,2] * v[2,3]) * cache$s1A + 0.5 * (v[1,1] * v[2,2] + v[1,2] * v[2,1]) * cache$s01A
            compE[focal,] <- phi[1:3]
            if(!is.null(node)){
				if(node == focal){
					fixer = c(0,0,0)
					fixer[state] = 1
					compD[focal,] <- compD[focal,] * fixer
				}
			}			
		}
		###########################
		#Logcompensation bit for dealing with underflow issues. Need to give a necessary shoutout to Rich FitzJohn -- we follow his diversitree approach. VERIFIED that it works properly:
        tmp <- sum(compD[focal,])
        compD[focal,] <- compD[focal,] / tmp
        logcomp <- c(logcomp, log(tmp))
        
	}
	root.node <- nb.tip + 1L
	if (is.na(sum(log(compD[root.node,]))) || is.na(log(sum(1-compE[root.node,])))){
		return(bad.likelihood)
	}else{
		if(root.type == "madfitz"){
            root.p = compD[root.node,]/sum(compD[root.node,])
            root.p[which(is.na(root.p))] = 0
		}
		if(root.type == "equal"){
			root.p = c(rep(1/length(which(compD[root.node,] > 0)), length(compD[root.node,])))
			root.p[which(!compD[root.node,] > 0)] = 0
		}
		if(root.type == "user"){
			root.p = root.p
		}
		if(condition.on.survival == TRUE){
			if(hidden.states == FALSE){
                lambda <- c(cache$s0A, cache$s1A, sum(c(cache$s0A, cache$s1A, cache$s01A)))
				compD[root.node,] <- compD[root.node,] / sum(root.p * lambda * (1 - compE[root.node,])^2)
				#Corrects for possibility that you have 0/0:
				compD[root.node,which(is.na(compD[root.node,]))] = 0
			}else{
                lambda <- c(cache$s0A, cache$s1A, sum(c(cache$s0A, cache$s1A, cache$s01A)), cache$s0B, cache$s1B, sum(c(cache$s0B, cache$s1B, cache$s01B)), cache$s0C, cache$s1C, sum(c(cache$s0C, cache$s1C, cache$s01C)), cache$s0D, cache$s1D, sum(c(cache$s0D, cache$s1D, cache$s01D)), cache$s0E, cache$s1E, sum(c(cache$s0E, cache$s1E, cache$s01E)))
                compD[root.node,] <- compD[root.node,] / sum(root.p * lambda * (1 - compE[root.node,])^2)
                #Corrects for possibility that you have 0/0:
				compD[root.node,which(is.na(compD[root.node,]))] = 0
			}
		}
		loglik <- log(sum(compD[root.node,] * root.p)) + sum(logcomp)
	}
	if(get.phi==TRUE){
		obj = NULL
		obj$compD.root = compD[root.node,]
		obj$compE = compE
		return(obj)
	}else{
		return(loglik)
	}
}


######################################################################################################################################
######################################################################################################################################
### Cache object for storing parameters that are used throughout HiGeoSSE:
######################################################################################################################################
######################################################################################################################################

ParametersToPassHiGeoSSE <- function(phy, data, f, model.vec, hidden.states){
    #Provides an initial object that contains all the parameters to be passed among functions. This will also be used to pass other things are we move down the tree (see DownPassGeoSSE):
    obj <- NULL
    obj$phy = phy
    
    if(hidden.states == FALSE){
        states = matrix(0,Ntip(phy),3)
        for(i in 1:Ntip(phy)){
            if(data[i]==1){states[i,1]=1}
            if(data[i]==2){states[i,2]=1}
            if(data[i]==0){states[i,3]=1}
        }
    }
    if(hidden.states == TRUE){
        states = matrix(0,Ntip(phy),15)
        for(i in 1:Ntip(phy)){
            if(data[i]==1){states[i,c(1,4,7,10,13)]=1}
            if(data[i]==2){states[i,c(2,5,8,11,14)]=1}
            if(data[i]==0){states[i,c(3,6,9,12,15)]=1}
        }
    }
    obj$states = states
    obj$tot_time = max(branching.times(phy))
    obj$f = f

    obj$s0A = model.vec[1]
    obj$s1A = model.vec[2]
    obj$s01A = model.vec[3]
    obj$x0A = model.vec[4]
    obj$x1A = model.vec[5]
    obj$d0A_01A = model.vec[6]
    obj$d1A_01A = model.vec[7]
    
    obj$d0A_0B = model.vec[8]
    obj$d0A_0C = model.vec[9]
    obj$d0A_0D = model.vec[10]
    obj$d0A_0E = model.vec[11]
    obj$d1A_1B = model.vec[12]
    obj$d1A_1C = model.vec[13]
    obj$d1A_1D = model.vec[14]
    obj$d1A_1E = model.vec[15]
    obj$d01A_01B = model.vec[16]
    obj$d01A_01C = model.vec[17]
    obj$d01A_01D = model.vec[18]
    obj$d01A_01E = model.vec[19]

    obj$s0B = model.vec[20]
    obj$s1B = model.vec[21]
    obj$s01B = model.vec[22]
    obj$x0B = model.vec[23]
    obj$x1B = model.vec[24]
    obj$d0B_01B = model.vec[25]
    obj$d1B_01B = model.vec[26]
    
    obj$d0B_0A = model.vec[27]
    obj$d0B_0C = model.vec[28]
    obj$d0B_0D = model.vec[29]
    obj$d0B_0E = model.vec[30]
    obj$d1B_1A = model.vec[31]
    obj$d1B_1C = model.vec[32]
    obj$d1B_1D = model.vec[33]
    obj$d1B_1E = model.vec[34]
    obj$d01B_01A = model.vec[35]
    obj$d01B_01C = model.vec[36]
    obj$d01B_01D = model.vec[37]
    obj$d01B_01E = model.vec[38]

    obj$s0C = model.vec[39]
    obj$s1C = model.vec[40]
    obj$s01C = model.vec[41]
    obj$x0C = model.vec[42]
    obj$x1C = model.vec[43]
    obj$d0C_01C = model.vec[44]
    obj$d1C_01C = model.vec[45]
    
    obj$d0C_0A = model.vec[46]
    obj$d0C_0B = model.vec[47]
    obj$d0C_0D = model.vec[48]
    obj$d0C_0E = model.vec[49]
    obj$d1C_1A = model.vec[50]
    obj$d1C_1B = model.vec[51]
    obj$d1C_1D = model.vec[52]
    obj$d1C_1E = model.vec[53]
    obj$d01C_01A = model.vec[54]
    obj$d01C_01B = model.vec[55]
    obj$d01C_01D = model.vec[56]
    obj$d01C_01E = model.vec[57]
    
    obj$s0D = model.vec[58]
    obj$s1D = model.vec[59]
    obj$s01D = model.vec[60]
    obj$x0D = model.vec[61]
    obj$x1D = model.vec[62]
    obj$d0D_01D = model.vec[63]
    obj$d1D_01D = model.vec[64]
    
    obj$d0D_0A = model.vec[65]
    obj$d0D_0B = model.vec[66]
    obj$d0D_0C = model.vec[67]
    obj$d0D_0E = model.vec[68]
    obj$d1D_1A = model.vec[69]
    obj$d1D_1B = model.vec[70]
    obj$d1D_1C = model.vec[71]
    obj$d1D_1E = model.vec[72]
    obj$d01D_01A = model.vec[73]
    obj$d01D_01B = model.vec[74]
    obj$d01D_01C = model.vec[75]
    obj$d01D_01E = model.vec[76]
    
    obj$s0E = model.vec[77]
    obj$s1E = model.vec[78]
    obj$s01E = model.vec[79]
    obj$x0E = model.vec[80]
    obj$x1E = model.vec[81]
    obj$d0E_01E = model.vec[82]
    obj$d1E_01E = model.vec[83]
    
    obj$d0E_0A = model.vec[84]
    obj$d0E_0B = model.vec[85]
    obj$d0E_0C = model.vec[86]
    obj$d0E_0D = model.vec[87]
    obj$d1E_1A = model.vec[88]
    obj$d1E_1B = model.vec[89]
    obj$d1E_1C = model.vec[90]
    obj$d1E_1ED = model.vec[91]
    obj$d01E_01A = model.vec[92]
    obj$d01E_01B = model.vec[93]
    obj$d01E_01C = model.vec[94]
    obj$d01E_01D = model.vec[95]

    obj$split.times = sort(branching.times(phy), decreasing=TRUE)
    
    return(obj)
}


######################################################################################################################################
######################################################################################################################################
### Print function for our diversity class:
######################################################################################################################################
######################################################################################################################################

#PUT HERE






######################################################################################################################################
######################################################################################################################################
### Code for testing purposes:
######################################################################################################################################
######################################################################################################################################

#pars <- c(1.5, 0.5, 1.0, 0.7, 0.7, 2.5, 0.5)
#names(pars) <- diversitree:::default.argnames.geosse()

## Simulate a tree
#set.seed(5)
#phy <- tree.geosse(pars, max.t=4, x0=0)
#data <- data.frame(g_s = names(phy$tip.state), states = phy$tip.state)
## See the data
#statecols <- c("AB"="purple", "A"="blue", "B"="red")

## The likelihood function
#lik <- make.geosse(phy, phy$tip.state)
#lik(pars)

#states <- data.frame(phy$tip.state, phy$tip.state, row.names=names(phy$tip.state))

#states <- states[phy$tip.label,]
#names(pars) <- NULL
#model.vec <- numeric(95)
#tests against GeoSSE
#model.vec <- rep(c(pars[1:3], pars[4:5], pars[6:7], rep(0,12)), 5)
#model.vec <- numeric(95)
#model.vec[1:7] <- c(pars[1:3], pars[4:5], pars[6:7])
#phy$node.label <- NULL
#cache <- ParametersToPassHiGeoSSE(phy, states[,1], f=c(1,1,1), model.vec, hidden.states=FALSE)
#ll.geosse <- DownPassHiGeosse(phy=phy, cache=cache, hidden.states=FALSE, bad.likelihood=-10000000000, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL)
#cache <- ParametersToPassHiGeoSSE(phy, states[,1], f=c(1,1,1), model.vec, hidden.states=TRUE)
#ll.higeosse <- DownPassHiGeosse(phy=phy, cache=cache, hidden.states=TRUE, bad.likelihood=-10000000000, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL)

#tests against NULL
#lambda_A <- pars[1]
#lambda_B <- pars[2]
#lambda_C <- pars[3]
#model.vec <- c(rep(lambda_A, 3), pars[4:5], pars[6:7], rep(1,6), rep(lambda_B, 3), pars[4:5], pars[6:7], rep(lambda_C, 3), pars[1:3], pars[4:5], pars[6:7], rep(1,6))

#phy$node.label <- NULL

#cache <- ParametersToPassGeoSSE(phy, states[,1], f=c(1,1,1), model.vec, hidden.states=TRUE)

#ll <- DownPassGeosse(phy=phy, cache=cache, hidden.states=TRUE, bad.likelihood=-10000000000, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL)




