

#library(ape)
#library(deSolve)
#library(subplex)
#library(phytools)
#library(nloptr)
#library(GenSA)
#library(data.table)
#dyn.load("misse-ext-derivs.so")
#dyn.load("birthdeath-ext-derivs.so")

######################################################################################################################################
######################################################################################################################################
### MiSSE -- Examines shifts in diversification in relation to hidden states ONLY
######################################################################################################################################
######################################################################################################################################

MiSSE <- function(phy, f=1, turnover=c(1,2), eps=c(1,2), condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, sann=FALSE, sann.its=10000, bounded.search=TRUE, max.tol=.Machine$double.eps^.50, starting.vals=NULL, turnover.upper=10000, eps.upper=3, trans.upper=100, restart.obj=NULL, ode.eps=0){
    
    ## Temporary fix for the current BUG:
    if( !is.null(phy$node.label) ) phy$node.label <- NULL
    
    if(!is.null(root.p)) {
        root.p <- root.p / sum(root.p)
        if(hidden.states ==TRUE){
            if(length(root.p) < 26){
                root.p.new <- numeric(26)
                root.p.new[1:length(root.p)] <- root.p
                root.p <- root.p.new
                root.p <- root.p / sum(root.p)
            }
        }
    }

    if(!root.type == "madfitz" & !root.type == "herr_als"){
        stop("Check that you specified a proper root.type option. Options are 'madfitz' or 'herr_als'. See help for more details.", call.=FALSE)
    }
    
    if(length(turnover) != length(eps)){
        stop("The number of turnover parameters need to match the number of extinction fraction parameters.", call.=FALSE)
    }


    ntips <- Ntip(phy)
    param.count <- sum(c(length(unique(turnover)), length(unique(eps)), 1))
    if(param.count > (ntips/10)){
        warning("You might not have enough data to fit this model well", call.=FALSE, immediate.=TRUE)
    }

    pars <- numeric(53)
    rate.cats <- hidden.states <- length(turnover)
    turnover.tmp <- numeric(26)
    turnover.tmp[1:length(turnover)] <- turnover
    pars.tmp <- turnover
    eps.tmp <- numeric(26)
    eps.tmp[1:length(eps)] <- eps
    eps.tmp[which(eps.tmp > 0)] = (eps.tmp[which( eps.tmp > 0)] + max(pars.tmp))
    pars.tmp <- c(pars.tmp, eps.tmp)
    if(rate.cats > 1){
        trans.tmp <- 1
        trans.tmp <- trans.tmp + max(pars.tmp)
    }else{
        trans.tmp <- 0
    }
    pars.tmp <- c(turnover.tmp[1], eps.tmp[1], turnover.tmp[2], eps.tmp[2], turnover.tmp[3], eps.tmp[3], turnover.tmp[4], eps.tmp[4], turnover.tmp[5], eps.tmp[5], turnover.tmp[6], eps.tmp[6], turnover.tmp[7], eps.tmp[7], turnover.tmp[8], eps.tmp[8], turnover.tmp[9], eps.tmp[9], turnover.tmp[10], eps.tmp[10], turnover.tmp[11], eps.tmp[11], turnover.tmp[12], eps.tmp[12], turnover.tmp[13], eps.tmp[13], turnover.tmp[14], eps.tmp[14], turnover.tmp[15], eps.tmp[15], turnover.tmp[16], eps.tmp[16], turnover.tmp[17], eps.tmp[17], turnover.tmp[18], eps.tmp[18], turnover.tmp[19], eps.tmp[19], turnover.tmp[20], eps.tmp[20], turnover.tmp[21], eps.tmp[21], turnover.tmp[22], eps.tmp[22], turnover.tmp[23], eps.tmp[23], turnover.tmp[24], eps.tmp[24], turnover.tmp[25], eps.tmp[25], turnover.tmp[26], eps.tmp[26], trans.tmp[1])
    pars[1:length(pars.tmp)] <- pars.tmp
    np <- max(pars)
    pars[pars==0] <- np + 1

    cat("Initializing...", "\n")

    #This is used to scale starting values to account for sampling:
    if(length(f) == 1){
            samp.freq.tree <- f
    }else{
        if(length(f) == Ntip(phy)){
            stop("This functionality has been temporarily removed.")
        }else{
            stop("The vector of sampling frequencies does not match the number of tips in the tree.")
        }
    }

    if(is.null(restart.obj)){
        if(sum(eps)==0){
            init.pars <- starting.point.generator(phy, 1, samp.freq.tree, yule=TRUE)
        }else{
            init.pars <- starting.point.generator(phy, 1, samp.freq.tree, yule=FALSE)
            if(any(init.pars[2] == 0)){
                init.pars[2] = 1e-6
            }
        }
        names(init.pars) <- NULL

        if(is.null(starting.vals)){
            def.set.pars <- rep(c(log(init.pars[1]+init.pars[2]), log(init.pars[2]/init.pars[1])), rate.cats)
            trans.start <- log(rate.cats/sum(phy$edge.length))
        }else{
            ## Check the length of the stating vector:
            if( length( starting.vals ) != 3 ){
                stop("Incorrect length for starting.vals vector.")
            }
            def.set.pars <- rep(c(log(starting.vals[1]), log(starting.vals[2])), rate.cats)
            trans.start <- log(starting.vals[3])
        }
        if(bounded.search == TRUE){
            upper.full <- rep(c(log(turnover.upper), log(eps.upper)), rate.cats)
        }else{
            upper.full <- rep(21, length(def.set.pars))
        }

        if(rate.cats > 1){
            np.sequence <- 1:(np-1)
            ip <- numeric(np-1)
            upper <- numeric(np-1)
        }else{
            np.sequence <- 1:np
            ip <- numeric(np)
            upper <- numeric(np)
        }

        for(i in np.sequence){
            ip[i] <- def.set.pars[which(pars == np.sequence[i])[1]]
            upper[i] <- upper.full[which(pars == np.sequence[i])[1]]
        }
        if(rate.cats > 1){
            ip <- c(ip, trans.start)
            upper <- c(upper, log(trans.upper))
        }
        lower <- rep(-20, length(ip))
    }else{
        upper <- restart.obj$upper.bounds
        lower <- restart.obj$lower.bounds
        pars <- restart.obj$index.par
        ip <- numeric(length(unique(restart.obj$index.par))-1)
        for(k in 1:length(ip)){
            ip[k] <- log(restart.obj$solution[which(restart.obj$index.par==k)][1])
        }
    }

    # Some new prerequisites #
    gen <- FindGenerations(phy)
    dat.tab <- OrganizeDataMiSSE(phy=phy, f=f, hidden.states=hidden.states)
    nb.tip <- Ntip(phy)
    nb.node <- phy$Nnode
    ##########################

    if(sann == FALSE){
        if(bounded.search == TRUE){
            cat("Finished. Beginning bounded subplex routine...", "\n")
            opts <- list("algorithm" = "NLOPT_LN_SBPLX", "maxeval" = 100000, "ftol_rel" = max.tol)
            out = nloptr(x0=ip, eval_f=DevOptimizeMiSSE, ub=upper, lb=lower, opts=opts, pars=pars, dat.tab=dat.tab, gen=gen, hidden.states=hidden.states, nb.tip=nb.tip, nb.node=nb.node, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, np=np, ode.eps=ode.eps)
            solution <- numeric(length(pars))
            solution[] <- c(exp(out$solution), 0)[pars]
            loglik = -out$objective
        }else{
            cat("Finished. Beginning subplex routine...", "\n")
            out = subplex(ip, fn=DevOptimizeMiSSE, control=list(reltol=max.tol, parscale=rep(0.1, length(ip))), pars=pars, dat.tab=dat.tab, gen=gen, hidden.states=hidden.states, nb.tip=nb.tip, nb.node=nb.node, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, np=np, ode.eps=ode.eps)
            solution <- numeric(length(pars))
            solution[] <- c(exp(out$par), 0)[pars]
            loglik = -out$value
        }
    }else{
        cat("Finished. Beginning simulated annealing...", "\n")
        out.sann = GenSA(ip, fn=DevOptimizeMiSSE, lower=lower, upper=upper, control=list(max.call=sann.its), pars=pars, dat.tab=dat.tab, gen=gen, hidden.states=hidden.states, nb.tip=nb.tip, nb.node=nb.node, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, np=np, ode.eps=ode.eps)
        cat("Finished. Refining using subplex routine...", "\n")
        opts <- list("algorithm" = "NLOPT_LN_SBPLX", "maxeval" = 100000, "ftol_rel" = max.tol)
        out <- nloptr(x0=out.sann$par, eval_f=DevOptimizeMiSSE, ub=upper, lb=lower, opts=opts, pars=pars, dat.tab=dat.tab, gen=gen, hidden.states=hidden.states, nb.tip=nb.tip, nb.node=nb.node, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, np=np, ode.eps=ode.eps)
        solution <- numeric(length(pars))
        solution[] <- c(exp(out$solution), 0)[pars]

        loglik = -out$objective
    }

    names(solution) <- c("turnover0A","eps0A", "turnover0B","eps0B", "turnover0C","eps0C", "turnover0D","eps0D", "turnover0E","eps0E", "turnover0F","eps0F", "turnover0G","eps0G", "turnover0H","eps0H", "turnover0I","eps0I", "turnover0J","eps0J", "turnover0K","eps0K", "turnover0L","eps0L", "turnover0M","eps0M", "turnover0N","eps0N", "turnover0O","eps0O", "turnover0P","eps0P", "turnover0Q","eps0Q", "turnover0R","eps0R", "turnover0S","eps0S", "turnover0T","eps0T", "turnover0U","eps0U", "turnover0V","eps0V","turnover0W","eps0W","turnover0X","eps0X", "turnover0Y","eps0Y", "turnover0Z","eps0Z", "q0")

    cat("Finished. Summarizing results...", "\n")

    obj = list(loglik = loglik, AIC = -2*loglik+2*np, AICc = -2*loglik+(2*np*(Ntip(phy)/(Ntip(phy)-np-1))), solution=solution, index.par=pars, f=f, hidden.states=hidden.states, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, phy=phy, max.tol=max.tol, starting.vals=ip, upper.bounds=upper, lower.bounds=lower, ode.eps=ode.eps, turnover=turnover, eps=eps)
    class(obj) <- append(class(obj), "misse.fit")
    return(obj)
}


# Note it'd be tempting to make this parallel. Go for it, but note that MiSSE itself already sprawls parallel (could be data.table) so there might not be the speedups you want -- BCO
MiSSEGreedy <- function(phy, f=1, turnover.tries=sequence(26), eps.constant=c(TRUE,FALSE), stop.count=2, stop.deltaAICc=10, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, sann=FALSE, sann.its=10000, bounded.search=TRUE, max.tol=.Machine$double.eps^.50, starting.vals=NULL, turnover.upper=10000, eps.upper=3, trans.upper=100, restart.obj=NULL, ode.eps=0) {
    misse.list <- list()
    first.AICc <- Inf
    for (eps.index in seq_along(eps.constant)) {
        best.AICc <- first.AICc #reset so we start over at one turnover parameter and work our way back up
        times.since.close.enough <- 0
        for (turnover.index in seq_along(turnover.tries)) {
            if(turnover.index>1 | eps.index==1) { # don't do the 1 turnover, 1 eps model twice -- overcounts it
                turnover <- sequence(turnover.tries[turnover.index])
                eps <- turnover
                if(eps.constant[eps.index]) {
                    eps <- rep(1, length(turnover))
                }
                starting.time <- Sys.time()
                cat("\nNow starting run with", turnover.tries[turnover.index], "turnover categories and", length(unique(eps)), "extinction fraction categories", "\n")
                cat("Starting at ", as.character(starting.time), "\n", sep="")
                current.run <- MiSSE(phy, f=f, turnover=turnover, eps=eps, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, sann=sann, sann.its=sann.its, bounded.search=bounded.search, max.tol=max.tol, starting.vals=starting.vals, turnover.upper=turnover.upper, eps.upper=eps.upper, trans.upper=trans.upper, restart.obj=restart.obj, ode.eps=ode.eps)
                misse.list[[length(misse.list)+1]] <- current.run
                ending.time <- Sys.time()
                cat("Finished at ", as.character(ending.time),"; this model took ", round(difftime(ending.time,starting.time, units="mins"),2), " minutes\n",sep="")
                
                cat("Current AICc is ", current.run$AICc, "\n", sep="")
                if(is.infinite(first.AICc)) {
                    first.AICc <- current.run$AICc # so when we reset, use the 1,1 AICc, not Inf
                }
                if(current.run$AICc < best.AICc) {
                    cat("Found better AICc by ",  best.AICc - current.run$AICc, "\n", sep="")
                    best.AICc <- current.run$AICc
                    times.since.close.enough <- 0
                } else if ((current.run$AICc - best.AICc ) < stop.deltaAICc) {
                    cat("Found worse AICc by ",  current.run$AICc - best.AICc , ", but this is still within ", stop.deltaAICc, " of the best", "\n", sep="")
                    times.since.close.enough <- 0
                } else {
                    times.since.close.enough <- times.since.close.enough + 1
                    cat("Found worse AICc by ",  current.run$AICc - best.AICc , ", it has been ", times.since.close.enough, " models since finding one within ", stop.deltaAICc, " of the best", "\n", sep="")
                    
                    if(times.since.close.enough > stop.count) {
                        break()
                    }
                }
            }
        }
    }
    return(misse.list)
}


######################################################################################################################################
######################################################################################################################################
### The function used to optimize parameters:
######################################################################################################################################
######################################################################################################################################


DevOptimizeMiSSE <- function(p, pars, dat.tab, gen, hidden.states, nb.tip=nb.tip, nb.node=nb.node, condition.on.survival, root.type, root.p, np, ode.eps) {
    #Generates the final vector with the appropriate parameter estimates in the right place:
    p.new <- exp(p)
    model.vec <- numeric(length(pars))
    model.vec[] <- c(p.new, 0)[pars]
    cache <- ParametersToPassMiSSE(model.vec=model.vec, hidden.states=hidden.states, nb.tip=nb.tip, nb.node=nb.node, bad.likelihood=exp(-500), ode.eps=ode.eps)
    logl <- DownPassMisse(dat.tab=dat.tab, cache=cache, gen=gen, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p)
    return(-logl)
}



######################################################################################################################################
######################################################################################################################################
### The various utility functions used
######################################################################################################################################
######################################################################################################################################


OrganizeDataMiSSE <- function(phy, f, hidden.states){
    ### Ughy McUgherson. This is a must in order to pass CRAN checks: http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
    DesNode = NULL

    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode

    compD <- matrix(0, nrow=nb.tip, ncol=26)
    compE <- matrix(0, nrow=nb.tip, ncol=26)
    #Initializes the tip sampling and sets internal nodes to be zero:
    ncols = hidden.states
    if(length(f) == 1){
        for(i in 1:(nb.tip)){
            compD[i,1:hidden.states] <- f
            compE[i,1:hidden.states] <- rep((1-f), ncols)
        }
    }else{
        for(i in 1:(nb.tip)){
            compD[i,] <- f[i]
            compE[i,] <- rep((1-f[i]), ncols)
        }
    }
    node.ages <- c(rep(0, Ntip(phy)), branching.times(phy))
    table.ages <- matrix(0, dim(phy$edge)[1], 2)
    for(row.index in 1:dim(phy$edge)[1]){
        table.ages[row.index,1] <- node.ages[phy$edge[row.index,1]]
        table.ages[row.index,2] <- node.ages[phy$edge[row.index,2]]
    }
    tmp.df <- cbind(table.ages, phy$edge.length, phy$edge[,1], phy$edge[,2], 0, matrix(0, nrow(table.ages), ncol(compD)), matrix(0, nrow(table.ages), ncol(compE)))
    colnames(tmp.df) <- c("RootwardAge", "TipwardAge", "BranchLength", "FocalNode", "DesNode", "comp", paste("compD", 1:ncol(compD), sep="_"), paste("compE", 1:ncol(compE), sep="_"))
    dat.tab <- as.data.table(tmp.df)
    setkey(dat.tab, DesNode)
    for (j in 1:(dim(compD)[2])){
        dat.tab[data.table(c(1:nb.tip)), paste("compD", j, sep="_") := compD[,j]]
        dat.tab[data.table(c(1:nb.tip)), paste("compE", j, sep="_") := compE[,j]]
    }
    return(dat.tab)
}


SingleChildProbMiSSE <- function(cache, compD, compE, start.time, end.time, x){
    pars <- c(cache$lambda0A, cache$mu0A, cache$lambda0B, cache$mu0B, cache$lambda0C, cache$mu0C, cache$lambda0D, cache$mu0D, cache$lambda0E, cache$mu0E, cache$lambda0F, cache$mu0F, cache$lambda0G, cache$mu0G, cache$lambda0H, cache$mu0H, cache$lambda0I, cache$mu0I, cache$lambda0J, cache$mu0J, cache$lambda0K, cache$mu0K, cache$lambda0L, cache$mu0L, cache$lambda0M, cache$mu0M, cache$lambda0N, cache$mu0N, cache$lambda0O, cache$mu0O, cache$lambda0P, cache$mu0P, cache$lambda0Q, cache$mu0Q, cache$lambda0R, cache$mu0R, cache$lambda0S, cache$mu0S, cache$lambda0T, cache$mu0T, cache$lambda0U, cache$mu0U, cache$lambda0V, cache$mu0V, cache$lambda0W, cache$mu0W, cache$lambda0X, cache$mu0X, cache$lambda0Y, cache$mu0Y, cache$lambda0Z, cache$mu0Z, cache$q0, cache$hidden.states)
    yini <- c(E0A = compE[1], E0B = compE[2], E0C = compE[3], E0D = compE[4], E0E = compE[5], E0F = compE[6], E0G = compE[7], E0H = compE[8], E0I = compE[9], E0J = compE[10], E0K = compE[11], E0L = compE[12], E0M = compE[13], E0N = compE[14], E0O = compE[15], E0P = compE[16], E0Q = compE[17], E0R = compE[18], E0S = compE[19], E0T = compE[20], E0U = compE[21], E0V = compE[22], E0W = compE[23], E0X = compE[24], E0Y = compE[25], E0Z = compE[26], D0A = compD[1], D0B = compD[2], D0C = compD[3], D0D = compD[4], D0E = compD[5], D0F = compD[6], D0G = compD[7], D0H = compD[8], D0I = compD[9], D0J = compD[10], D0K = compD[11], D0L = compD[12], D0M = compD[13], D0N = compD[14], D0O = compD[15], D0P = compD[16], D0Q = compD[17], D0R = compD[18], D0S = compD[19], D0T = compD[20], D0U = compD[21], D0V = compD[22], D0W = compD[23], D0X = compD[24], D0Y = compD[25], D0Z = compD[26])
    times=c(start.time, end.time)

    #prob.subtree.cal.full <- lsoda(yini, times, func = "misse_derivs", pars, initfunc="initmod_misse", dll = "misse-ext-derivs", rtol=1e-8, atol=1e-8)
    prob.subtree.cal.full <- lsoda(yini, times, func = "misse_derivs", pars, initfunc="initmod_misse", dllname = "hisse", rtol=1e-8, atol=1e-8)

    ######## THIS CHECKS TO ENSURE THAT THE INTEGRATION WAS SUCCESSFUL ###########
    if(attributes(prob.subtree.cal.full)$istate[1] < 0){
        prob.subtree.cal <- prob.subtree.cal.full[-1,-1]
        if(cache$hidden.states == TRUE){
            prob.subtree.cal[27:52] <- cache$bad.likelihood
            return(prob.subtree.cal)
        }
    }else{
        prob.subtree.cal <- prob.subtree.cal.full[-1,-1]
    }
    ##############################################################################

    if(any(is.nan(prob.subtree.cal[27:52]))){
        prob.subtree.cal[27:52] <- cache$bad.likelihood
        return(prob.subtree.cal)
    }
    #This is default and cannot change, but if we get a negative probability, discard the results:
    if(any(prob.subtree.cal[27:52] < 0)){
        prob.subtree.cal[27:52] <- cache$bad.likelihood
        return(prob.subtree.cal)
    }
    if(sum(prob.subtree.cal[27:52]) < cache$ode.eps){
        prob.subtree.cal[27:52] <- cache$bad.likelihood
        return(prob.subtree.cal)
    }

    return(prob.subtree.cal)
}


FocalNodeProbMiSSE <- function(cache, dat.tab, generations){
    ### Ughy McUgherson. This is a must in order to pass CRAN checks: http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
    DesNode = NULL
    FocalNode = NULL

    setkey(dat.tab, FocalNode)
    CurrentGenData <- dat.tab[data.table(generations)]
    tmp <- t(apply(CurrentGenData, 1, function(z) SingleChildProbMiSSE(cache, z[7:32], z[33:58],  z[2], z[1])))
    v.mat <- matrix(tmp[seq(1,nrow(tmp)-1,2),27:52] * tmp[seq(2,nrow(tmp),2),27:52], length(unique(CurrentGenData$FocalNode)), 26)
    v.mat <- v.mat * matrix(c(cache$lambda0A, cache$lambda0B, cache$lambda0C, cache$lambda0D, cache$lambda0E, cache$lambda0F, cache$lambda0G, cache$lambda0H, cache$lambda0I, cache$lambda0J, cache$lambda0K, cache$lambda0L, cache$lambda0M, cache$lambda0N, cache$lambda0O, cache$lambda0P, cache$lambda0Q, cache$lambda0R, cache$lambda0S, cache$lambda0T, cache$lambda0U, cache$lambda0V, cache$lambda0W, cache$lambda0X, cache$lambda0Y, cache$lambda0Z), length(unique(CurrentGenData$FocalNode)), 26, byrow=TRUE)
    phi.mat <- matrix(tmp[seq(1,nrow(tmp)-1,2),1:26], length(unique(CurrentGenData$FocalNode)), 26)
    if(!is.null(cache$node)){
        if(which(generations == cache$node)){
            fixer = numeric(26)
            fixer[cache$state] = 1
            v.mat[which(generations == cache$node),] <- v.mat[which(generations == cache$node),] * fixer
        }
    }
    tmp.comp <- rowSums(v.mat)
    tmp.probs <- v.mat / tmp.comp
    setkey(dat.tab, DesNode)
    for (j in 1:(dim(tmp.probs)[2])){
        dat.tab[data.table(c(generations)), paste("compD", j, sep="_") := tmp.probs[,j]]
        dat.tab[data.table(c(generations)), paste("compE", j, sep="_") := phi.mat[,j]]
    }
    dat.tab[data.table(c(generations)), "comp" := tmp.comp]
    return(dat.tab)
}


#Have to calculate root prob separately because it is not a descendant in our table. Could add it, but I worry about the NA that is required.
GetRootProbMiSSE <- function(cache, dat.tab, generations){
    ### Ughy McUgherson. This is a must in order to pass CRAN checks: http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
    FocalNode = NULL

    setkey(dat.tab, FocalNode)
    CurrentGenData <- dat.tab[data.table(generations)]
    tmp <- t(apply(CurrentGenData, 1, function(z) SingleChildProbMiSSE(cache, z[7:32], z[33:58],  z[2], z[1])))
    v.mat <- matrix(tmp[seq(1,nrow(tmp)-1,2),27:52] * tmp[seq(2,nrow(tmp),2),27:52], length(unique(CurrentGenData$FocalNode)), 26)
    v.mat <- v.mat * matrix(c(cache$lambda0A, cache$lambda0B, cache$lambda0C, cache$lambda0D, cache$lambda0E, cache$lambda0F, cache$lambda0G, cache$lambda0H, cache$lambda0I, cache$lambda0J, cache$lambda0K, cache$lambda0L, cache$lambda0M, cache$lambda0N, cache$lambda0O, cache$lambda0P, cache$lambda0Q, cache$lambda0R, cache$lambda0S, cache$lambda0T, cache$lambda0U, cache$lambda0V, cache$lambda0W, cache$lambda0X, cache$lambda0Y, cache$lambda0Z), length(unique(CurrentGenData$FocalNode)), 26, byrow=TRUE)
    phi.mat <- matrix(tmp[seq(1,nrow(tmp)-1,2),1:26], length(unique(CurrentGenData$FocalNode)), 26)
    if(!is.null(cache$node)){
        if(which(generations == cache$node)){
            fixer = numeric(26)
            fixer[cache$state] = 1
            v.mat[which(generations == cache$node),] <- v.mat[which(generations == cache$node),] * fixer
        }
    }

    tmp.comp <- rowSums(v.mat)
    tmp.probs <- v.mat / tmp.comp

    return(cbind(tmp.comp, phi.mat, tmp.probs))
}



######################################################################################################################################
######################################################################################################################################
### The MiSSE type down pass that carries out the integration and returns the likelihood:
######################################################################################################################################
######################################################################################################################################

DownPassMisse <- function(dat.tab, gen, cache, condition.on.survival, root.type, root.p, get.phi=FALSE, node=NULL, state=NULL) {

    ### Ughy McUgherson. This is a must in order to pass CRAN checks: http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
    DesNode = NULL
    compE = NULL

    nb.tip <- cache$nb.tip
    nb.node <- cache$nb.node
    TIPS <- 1:nb.tip
    for(i in 1:length(gen)){
        if(i == length(gen)){
            if(!is.null(node)){
                if(node %in% gen[[i]]){
                    cache$node <- node
                    cache$state <- state
                    res.tmp <- GetRootProbMiSSE(cache=cache, dat.tab=dat.tab, generations=gen[[i]])
                    cache$node <- NULL
                    cache$state <- NULL
                }else{
                    res.tmp <- GetRootProbMiSSE(cache=cache, dat.tab=dat.tab, generations=gen[[i]])
                }
            }else{
                res.tmp <- GetRootProbMiSSE(cache=cache, dat.tab=dat.tab, generations=gen[[i]])
            }
            compD.root <- res.tmp[c(28:53)]
            compE.root <- res.tmp[c(2:27)]
            setkey(dat.tab, DesNode)
            comp <- dat.tab[["comp"]]
            comp <- c(comp[-TIPS], res.tmp[1])
        }else{
            if(!is.null(node)){
                if(node %in% gen[[i]]){
                    cache$node <- node
                    cache$state <- state
                    dat.tab <- FocalNodeProbMiSSE(cache, dat.tab, gen[[i]])
                    cache$node <- NULL
                    cache$state <- NULL
                }else{
                    dat.tab <- FocalNodeProbMiSSE(cache, dat.tab, gen[[i]])
                }
            }else{
                dat.tab <- FocalNodeProbMiSSE(cache, dat.tab, gen[[i]])
            }
        }
    }
    if (is.na(sum(log(compD.root))) || is.na(log(sum(1-compE.root)))){
        return(log(cache$bad.likelihood)^13)
    }else{
        if(root.type == "madfitz" | root.type == "herr_als"){
            if(is.null(root.p)){
                root.p = compD.root/sum(compD.root)
                root.p[which(is.na(root.p))] = 0
            }
        }
        if(condition.on.survival == TRUE){
            if(root.type == "madfitz"){
                lambda <- c(cache$lambda0A, cache$lambda0B, cache$lambda0C, cache$lambda0D, cache$lambda0E, cache$lambda0F, cache$lambda0G, cache$lambda0H, cache$lambda0I, cache$lambda0J, cache$lambda0K, cache$lambda0L, cache$lambda0M, cache$lambda0N, cache$lambda0O, cache$lambda0P, cache$lambda0Q, cache$lambda0R, cache$lambda0S, cache$lambda0T, cache$lambda0U, cache$lambda0V, cache$lambda0W, cache$lambda0X, cache$lambda0Y, cache$lambda0Z)
                compD.root <- compD.root / sum(root.p * lambda * (1 - compE.root)^2)
                #Corrects for possibility that you have 0/0:
                compD.root[which(is.na(compD.root))] = 0
                loglik <- log(sum(compD.root * root.p)) + sum(log(comp))
            }else{
                lambda <- c(cache$lambda0A, cache$lambda0B, cache$lambda0C, cache$lambda0D, cache$lambda0E, cache$lambda0F, cache$lambda0G, cache$lambda0H, cache$lambda0I, cache$lambda0J, cache$lambda0K, cache$lambda0L, cache$lambda0M, cache$lambda0N, cache$lambda0O, cache$lambda0P, cache$lambda0Q, cache$lambda0R, cache$lambda0S, cache$lambda0T, cache$lambda0U, cache$lambda0V, cache$lambda0W, cache$lambda0X, cache$lambda0Y, cache$lambda0Z)
                compD.root <- (compD.root * root.p) / (lambda * (1 - compE.root)^2)
                #Corrects for possibility that you have 0/0:
                compD.root[which(is.na(compD.root))] = 0
                loglik <- log(sum(compD.root)) + sum(log(comp))
            }
        }
        if(!is.finite(loglik)){
            return(log(cache$bad.likelihood)^7)
        }
    }
    if(get.phi==TRUE){
        obj = NULL
        obj$compD.root = compD.root/sum(compD.root)
        obj$compE = compE
        obj$root.p = root.p
        return(obj)
    }else{
        return(loglik)
    }
}


######################################################################################################################################
######################################################################################################################################
### Cache object for storing parameters that are used throughout MiSSE:
######################################################################################################################################
######################################################################################################################################

ParametersToPassMiSSE <- function(model.vec, hidden.states, nb.tip, nb.node, bad.likelihood, ode.eps){
    #Provides an initial object that contains all the parameters to be passed among functions. This will also be used to pass other things are we move down the tree (see DownPassGeoSSE):
    obj <- NULL

    obj$hidden.states <- hidden.states
    obj$nb.tip <- nb.tip
    obj$nb.node <- nb.node
    obj$bad.likelihood <- bad.likelihood
    obj$ode.eps <- ode.eps

    ##Hidden State A
    obj$lambda0A = model.vec[1] / (1 + model.vec[2])
    obj$mu0A = (model.vec[2] * model.vec[1]) / (1 + model.vec[2])

    ##Hidden State B
    obj$lambda0B = model.vec[3] / (1 + model.vec[4])
    obj$mu0B = (model.vec[4] * model.vec[3]) / (1 + model.vec[4])

    ##Hidden State C
    obj$lambda0C = model.vec[5] / (1 + model.vec[6])
    obj$mu0C = (model.vec[6] * model.vec[5]) / (1 + model.vec[6])

    ##Hidden State D
    obj$lambda0D = model.vec[7] / (1 + model.vec[8])
    obj$mu0D = (model.vec[8] * model.vec[7]) / (1 + model.vec[8])

    ##Hidden State E
    obj$lambda0E = model.vec[9] / (1 + model.vec[10])
    obj$mu0E = (model.vec[10] * model.vec[9]) / (1 + model.vec[10])

    ##Hidden State F
    obj$lambda0F = model.vec[11] / (1 + model.vec[12])
    obj$mu0F = (model.vec[12] * model.vec[11]) / (1 + model.vec[12])

    ##Hidden State G
    obj$lambda0G = model.vec[13] / (1 + model.vec[14])
    obj$mu0G = (model.vec[14] * model.vec[13]) / (1 + model.vec[14])

    ##Hidden State H
    obj$lambda0H = model.vec[15] / (1 + model.vec[16])
    obj$mu0H = (model.vec[16] * model.vec[15]) / (1 + model.vec[16])

    ##Hidden State I
    obj$lambda0I = model.vec[17] / (1 + model.vec[18])
    obj$mu0I = (model.vec[18] * model.vec[17]) / (1 + model.vec[18])

    ##Hidden State J
    obj$lambda0J = model.vec[19] / (1 + model.vec[20])
    obj$mu0J = (model.vec[20] * model.vec[19]) / (1 + model.vec[20])

    ##Hidden State K
    obj$lambda0K = model.vec[21] / (1 + model.vec[22])
    obj$mu0K = (model.vec[22] * model.vec[21]) / (1 + model.vec[22])

    ##Hidden State L
    obj$lambda0L = model.vec[23] / (1 + model.vec[24])
    obj$mu0L = (model.vec[24] * model.vec[23]) / (1 + model.vec[24])

    ##Hidden State M
    obj$lambda0M = model.vec[25] / (1 + model.vec[26])
    obj$mu0M = (model.vec[26] * model.vec[25]) / (1 + model.vec[26])

    ##Hidden State N
    obj$lambda0N = model.vec[27] / (1 + model.vec[28])
    obj$mu0N = (model.vec[28] * model.vec[27]) / (1 + model.vec[28])

    ##Hidden State O
    obj$lambda0O = model.vec[29] / (1 + model.vec[30])
    obj$mu0O = (model.vec[30] * model.vec[29]) / (1 + model.vec[30])

    ##Hidden State P
    obj$lambda0P = model.vec[31] / (1 + model.vec[32])
    obj$mu0P = (model.vec[32] * model.vec[31]) / (1 + model.vec[32])

    ##Hidden State Q
    obj$lambda0Q = model.vec[33] / (1 + model.vec[34])
    obj$mu0Q = (model.vec[34] * model.vec[33]) / (1 + model.vec[34])

    ##Hidden State R
    obj$lambda0R = model.vec[35] / (1 + model.vec[36])
    obj$mu0R = (model.vec[36] * model.vec[35]) / (1 + model.vec[36])

    ##Hidden State S
    obj$lambda0S = model.vec[37] / (1 + model.vec[38])
    obj$mu0S = (model.vec[38] * model.vec[37]) / (1 + model.vec[38])

    ##Hidden State T
    obj$lambda0T = model.vec[39] / (1 + model.vec[40])
    obj$mu0T = (model.vec[40] * model.vec[39]) / (1 + model.vec[40])

    ##Hidden State U
    obj$lambda0U = model.vec[41] / (1 + model.vec[42])
    obj$mu0U = (model.vec[42] * model.vec[41]) / (1 + model.vec[42])

    ##Hidden State V
    obj$lambda0V = model.vec[43] / (1 + model.vec[44])
    obj$mu0V = (model.vec[44] * model.vec[43]) / (1 + model.vec[44])

    ##Hidden State W
    obj$lambda0W = model.vec[45] / (1 + model.vec[46])
    obj$mu0W = (model.vec[46] * model.vec[45]) / (1 + model.vec[46])

    ##Hidden State X
    obj$lambda0X = model.vec[47] / (1 + model.vec[48])
    obj$mu0X = (model.vec[48] * model.vec[47]) / (1 + model.vec[48])

    ##Hidden State Y
    obj$lambda0Y = model.vec[49] / (1 + model.vec[50])
    obj$mu0Y = (model.vec[50] * model.vec[49]) / (1 + model.vec[50])

    ##Hidden State Z
    obj$lambda0Z = model.vec[51] / (1 + model.vec[52])
    obj$mu0Z = (model.vec[52] * model.vec[51]) / (1 + model.vec[52])

    obj$q0 = model.vec[53]
    return(obj)
}


print.misse.fit <- function(x,...){
    ## Function to print a "muhisse.fit" object.
    set.zero <- max( x$index.par )
    ## Keep only the parameters estimated:
    par.list <- x$solution[ !x$index.par == set.zero ]
    ntips <- Ntip( x$phy )
    nstates <- x$hidden.states
    output <- c(x$loglik, x$AIC, x$AICc, ntips, nstates)
    names(output) <- c("lnL", "AIC", "AICc", "n.taxa", "n.hidden.states")
    cat("\n")
    cat("Fit \n")
    print(output)
    cat("\n")
    cat("Model parameters: \n")
    cat("\n")
    print(par.list)
    cat("\n")
}




#phy <- read.tree("../vignettes/whales_Steemanetal2009.tre")
#phy <- read.tree("whales_Slateretal2010.tre")
## print(p.new)
#gen <- hisse:::FindGenerations(phy)
#dat.tab <- hisse:::OrganizeDataMiSSE(phy=phy, f=1, hidden.states=1)
#nb.tip <- Ntip(phy)
#nb.node <- phy$Nnode
#model.vec <- c(0.103624, 5.207178e-09, rep(0,52), 1)

#cache = hisse:::ParametersToPassMiSSE(model.vec=model.vec, hidden.states=1, nb.tip=nb.tip, nb.node=nb.node, bad.likelihood=exp(-500), ode.eps=0)#
#logl <- hisse:::DownPassMisse(dat.tab=dat.tab, cache=cache, gen=gen, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL)
#right.logl <- -277.6942
#round(logl,4) == round(right.logl,4)
