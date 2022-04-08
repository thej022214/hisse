
######################################################################################################################################
######################################################################################################################################
### HiSSE -- a faster version assumes four possible hidden states, associated with or without particular character states
######################################################################################################################################
######################################################################################################################################

hisse <- function(phy, data, f=c(1,1), turnover=c(1,2), eps=c(1,2), hidden.states=FALSE, trans.rate=NULL, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, includes.fossils=FALSE, k.samples=NULL, strat.intervals=NULL, sann=TRUE, sann.its=5000, bounded.search=TRUE, max.tol=.Machine$double.eps^.50, starting.vals=NULL, turnover.upper=10000, eps.upper=3, trans.upper=100, restart.obj=NULL, ode.eps=0, dt.threads=1){
    
    ## Temporary fix for the current BUG:
    if( !is.null(phy$node.label) ) phy$node.label <- NULL
    
    if(!is.ultrametric(phy) & includes.fossils == FALSE){
        warning("Tree is not ultrametric. Used force.ultrametric() function to coerce the tree to be ultrametric - see note above.")
        edge_details <- GetEdgeDetails(phy, includes.intervals=FALSE, intervening.intervals=NULL)
        if(any(edge_details$type == "extinct_tip")){
            phy <- force.ultrametric(phy)
        }
    }

    if(!is.null(root.p)) {
        if(hidden.states ==TRUE){
            if( length( root.p ) == 2 ){
                root.p <- rep(root.p, 4)
                root.p <- root.p / sum(root.p)
                warning("For hidden states, you need to specify the root.p for all hidden states. We have adjusted it so that there's equal chance for among all hidden states.")
            } else{
                root.p.new <- numeric(8)
                root.p.new[1:length(root.p)] <- root.p
                root.p <- root.p.new
                root.p <- root.p / sum(root.p)
            }
        }else{
            ## All good:
            root.p <- root.p / sum(root.p)
        }
    }
    
    setDTthreads(threads=dt.threads)
    
    if(sann == FALSE & is.null(starting.vals)){
        warning("You have chosen to rely on the internal starting points that generally work but does not guarantee finding the MLE.")
    }
    
    if(!root.type == "madfitz" & !root.type == "herr_als"){
        stop("Check that you specified a proper root.type option. Options are 'madfitz' or 'herr_als'. See help for more details.", call.=FALSE)
    }
    
    if(is.null(trans.rate)){
        stop("Rate matrix needed. See TransMatMakerHiSSE() to create one.")
    }
    
    if(hidden.states == TRUE & dim(trans.rate)[1]<4){
        stop("You chose a hidden state but this is not reflected in the transition matrix")
    }
    
    ## Return error message if the data is not in the correct format.
    if( !inherits(data, what = c("matrix", "data.frame")) ){
        stop("'data' needs to be a matrix or data.frame with 2 columns. See help.")
    }
    if( !ncol( data ) == 2 ){
        stop("'data' needs to be a matrix or data.frame with 2 columns. See help.")
    }
    
    ## Check if 'hidden.states' parameter is congruent with the turnover vector:
    if( length(turnover) > 2 & !hidden.states ){
        stop("Turnover has more than 2 elements but 'hidden.states' was set to FALSE. Please set 'hidden.states' to TRUE if the model include more than one rate class.")
    }
    if( length(turnover) == 2 & hidden.states ){
        stop("Turnover has only 2 elements but 'hidden.states' was set to TRUE. Please set 'hidden.states' to FALSE if the model does not include hidden rate classes.")
    }
    
    pars <- numeric(48)
    
    if(dim(trans.rate)[2]==2){
        rate.cats <- 1
        pars.tmp <- turnover
        eps.tmp <- eps
        eps.tmp[which(eps.tmp > 0)] = (eps.tmp[which( eps.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, eps.tmp)
        trans.tmp <- c(trans.rate["(0)", "(1)"], trans.rate["(1)", "(0)"])
        trans.tmp[which(trans.tmp > 0)] = (trans.tmp[which(trans.tmp > 0)] + max(pars.tmp))
        category.rates.unique <- 0
        pars.tmp <- c(pars.tmp, trans.tmp)
        pars[1:6] <- pars.tmp
    }
    
    if(dim(trans.rate)[2]==4){
        rate.cats <- 2
        pars.tmp <- turnover
        eps.tmp <- eps
        eps.tmp[which(eps.tmp > 0)] = (eps.tmp[which( eps.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, eps.tmp)
        for.late.adjust <- max(pars.tmp)
        rows <- c("(0A)", "(1A)", "(0B)", "(1B)")
        cols <- c("(1A)", "(0A)", "(1B)", "(0B)")
        trans.tmp <- trans.rate[cbind(rows,cols)]
        trans.tmp[which(trans.tmp > 0)] = (trans.tmp[which(trans.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, trans.tmp)
        rows <- c("(0A)", "(1A)", "(0B)", "(1B)")
        cols <- c("(0B)", "(1B)", "(0A)", "(1A)")
        category.tmp <- trans.rate[cbind(rows,cols)]
        category.tmp[category.tmp==0] <- NA
        category.rate.shift <- category.tmp + for.late.adjust
        category.rate.shift[is.na(category.rate.shift)] <- 0
        category.rate.shiftA <- c(category.rate.shift[1], rep(0,2), category.rate.shift[2], rep(0,2))
        category.rate.shiftB <- c(category.rate.shift[3], rep(0,2), category.rate.shift[4], rep(0,2))
        pars.tmp <- c(turnover[1:2], eps.tmp[1:2], trans.tmp[1:2], category.rate.shiftA, turnover[3:4], eps.tmp[3:4], trans.tmp[3:4], category.rate.shiftB)
        pars[1:length(pars.tmp)] <- pars.tmp
    }
    
    if(dim(trans.rate)[2]==6){
        rate.cats <- 3
        pars.tmp <- turnover
        eps.tmp <- eps
        eps.tmp[which(eps.tmp > 0)] = (eps.tmp[which( eps.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, eps.tmp)
        for.late.adjust <- max(pars.tmp)
        rows <- c("(0A)", "(1A)", "(0B)", "(1B)", "(0C)", "(1C)")
        cols <- c("(1A)", "(0A)", "(1B)", "(0B)", "(1C)", "(0C)")
        trans.tmp <- trans.rate[cbind(rows,cols)]
        trans.tmp[which(trans.tmp > 0)] = (trans.tmp[which(trans.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, trans.tmp)
        rows <- c("(0A)", "(0A)", "(1A)", "(1A)", "(0B)", "(0B)", "(1B)", "(1B)", "(0C)", "(0C)", "(1C)", "(1C)")
        cols <- c("(0B)", "(0C)", "(1B)", "(1C)", "(0A)", "(0C)", "(1A)", "(1C)", "(0A)", "(0B)", "(1A)", "(1B)")
        category.tmp <- trans.rate[cbind(rows,cols)]
        category.tmp[category.tmp==0] <- NA
        category.rate.shift <- category.tmp + for.late.adjust
        category.rate.shift[is.na(category.rate.shift)] <- 0
        category.rate.shiftA <- c(category.rate.shift[1:2], rep(0,1), category.rate.shift[3:4], rep(0,1))
        category.rate.shiftB <- c(category.rate.shift[5:6], rep(0,1), category.rate.shift[7:8], rep(0,1))
        category.rate.shiftC <- c(category.rate.shift[9:10], rep(0,1), category.rate.shift[11:12], rep(0,1))
        pars.tmp <- c(turnover[1:2], eps.tmp[1:2], trans.tmp[1:2], category.rate.shiftA, turnover[3:4], eps.tmp[3:4], trans.tmp[3:4], category.rate.shiftB, turnover[5:6], eps.tmp[5:6], trans.tmp[5:6], category.rate.shiftC)
        pars[1:length(pars.tmp)] <- pars.tmp
    }
    
    if(dim(trans.rate)[2]==8){
        rate.cats <- 4
        pars.tmp <- turnover
        eps.tmp <- eps
        eps.tmp[which(eps.tmp > 0)] = (eps.tmp[which( eps.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, eps.tmp)
        for.late.adjust <- max(pars.tmp)
        rows <- c("(0A)", "(1A)", "(0B)", "(1B)", "(0C)", "(1C)", "(0D)", "(1D)")
        cols <- c("(1A)", "(0A)", "(1B)", "(0B)", "(1C)", "(0C)", "(1D)", "(0D)")
        trans.tmp <- trans.rate[cbind(rows,cols)]
        trans.tmp[which(trans.tmp > 0)] = (trans.tmp[which(trans.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, trans.tmp)
        rows <- c("(0A)", "(0A)", "(0A)", "(1A)", "(1A)", "(1A)", "(0B)", "(0B)", "(0B)", "(1B)", "(1B)", "(1B)", "(0C)", "(0C)", "(0C)", "(1C)", "(1C)", "(1C)", "(0D)", "(0D)", "(0D)", "(1D)", "(1D)", "(1D)")
        cols <- c("(0B)", "(0C)", "(0D)", "(1B)", "(1C)", "(1D)", "(0A)", "(0C)", "(0D)", "(1A)", "(1C)", "(1D)", "(0A)", "(0B)", "(0D)", "(1A)", "(1B)", "(1D)", "(0A)", "(0B)", "(0C)", "(1A)", "(1B)", "(1C)")
        category.tmp <- trans.rate[cbind(rows,cols)]
        category.tmp[category.tmp==0] <- NA
        category.rate.shift <- category.tmp + for.late.adjust
        category.rate.shift[is.na(category.rate.shift)] <- 0
        category.rate.shiftA <- c(category.rate.shift[1:3], category.rate.shift[4:6])
        category.rate.shiftB <- c(category.rate.shift[7:9], category.rate.shift[10:12])
        category.rate.shiftC <- c(category.rate.shift[13:15], category.rate.shift[16:18])
        category.rate.shiftD <- c(category.rate.shift[19:21], category.rate.shift[22:24])
        pars.tmp <- c(turnover[1:2], eps.tmp[1:2], trans.tmp[1:2], category.rate.shiftA, turnover[3:4], eps.tmp[3:4], trans.tmp[3:4], category.rate.shiftB, turnover[5:6], eps.tmp[5:6], trans.tmp[5:6], category.rate.shiftC, turnover[7:8], eps.tmp[7:8], trans.tmp[7:8], category.rate.shiftD)
        pars[1:length(pars.tmp)] <- pars.tmp
    }
    
    if(includes.fossils == TRUE){
        pars <- c(pars, max(pars.tmp)+1)
    }else{
        pars <- c(pars, 0)
    }
    
    np <- max(pars)
    pars[pars==0] <- np+1
    
    cat("Initializing...", "\n")
    
    # Some new prerequisites #
    if(includes.fossils == TRUE){
        if(!is.null(k.samples)){
            phy.og <- phy
            psi.type <- "m+k"
            split.times <- dateNodes(phy, rootAge=max(node.depth.edgelength(phy)))[-c(1:Ntip(phy))]
            strat.cache <- NULL
            k.samples <- k.samples[order(as.numeric(k.samples[,3]), decreasing=FALSE),]
            phy <- AddKNodes(phy, k.samples)
            fix.type <- GetKSampleMRCA(phy, k.samples)
            no.k.samples <- length(k.samples[,1])
            gen <- FindGenerations(phy)
            dat.tab <- OrganizeDataHiSSE(phy=phy, f=f, hidden.states=hidden.states, includes.intervals=FALSE, intervening.intervals=NULL)
            #These are all inputs for generating starting values:
            edge_details <- GetEdgeDetails(phy, includes.intervals=FALSE, intervening.intervals=NULL)
            fossil.taxa <- edge_details$tipward_node[which(edge_details$type == "extinct_tip")]
            fossil.ages <- dat.tab$TipwardAge[which(dat.tab$DesNode %in% fossil.taxa)]
            n.tax.starting <- Ntip(phy)-length(fossil.taxa)-no.k.samples
            data <- AddKData(data, k.samples)
        }else{
            if(!is.null(strat.intervals)){
                phy.og <- phy
                psi.type <- "m+int"
                split.times.plus.tips <- dateNodes(phy, rootAge=max(node.depth.edgelength(phy)))
                split.times <- split.times.plus.tips[-c(1:Ntip(phy))]
                strat.cache <- GetStratInfo(strat.intervals=strat.intervals)
                k.samples <- GetIntervalToK(strat.intervals, intervening.intervals=strat.cache$intervening.intervals)
                extinct.tips <- which(round(k.samples$timefrompresent,8) %in% round(split.times.plus.tips[c(1:Ntip(phy))],8))
                if(length(extinct.tips > 0)){
                    k.samples <- k.samples[-extinct.tips,]
                }
                phy <- AddKNodes(phy, k.samples)
                fix.type <- GetKSampleMRCA(phy, k.samples, strat.intervals=TRUE)
                gen <- FindGenerations(phy)
                dat.tab <- OrganizeDataHiSSE(phy=phy, f=f, hidden.states=hidden.states, includes.intervals=TRUE, intervening.intervals=strat.cache$intervening.intervals)
                #These are all inputs for generating starting values:
                edge_details <- GetEdgeDetails(phy, includes.intervals=TRUE, intervening.intervals=strat.cache$intervening.intervals)
                fossil.taxa <- edge_details$tipward_node[which(edge_details$type == "extinct_tip" | edge_details$type == "k_extinct_interval")]
                fossil.ages <- dat.tab$TipwardAge[which(dat.tab$DesNode %in% fossil.taxa)]
                data <- AddKData(data, k.samples)
                n.tax.starting <- Ntip(phy)-length(fossil.taxa)-dim(fix.type)[1]
            }else{
                phy.og <- phy
                psi.type <- "m_only"
                fix.type <- NULL
                split.times <- dateNodes(phy, rootAge=max(node.depth.edgelength(phy)))[-c(1:Ntip(phy))]
                strat.cache <- NULL
                no.k.samples <- 0
                gen <- FindGenerations(phy)
                dat.tab <- OrganizeDataMiSSE(phy=phy, f=f, hidden.states=hidden.states, includes.intervals=FALSE, intervening.intervals=NULL)
                #These are all inputs for generating starting values:
                edge_details <- GetEdgeDetails(phy, includes.intervals=FALSE, intervening.intervals=NULL)
                fossil.taxa <- edge_details$tipward_node[which(edge_details$type == "extinct_tip")]
                fossil.ages <- dat.tab$TipwardAge[which(dat.tab$DesNode %in% fossil.taxa)]
                n.tax.starting <- Ntip(phy)-length(fossil.taxa)-no.k.samples
            }
        }
    }else{
        phy.og <- phy
        gen <- FindGenerations(phy)
        data.new <- data.frame(data[,2], data[,2], row.names=data[,1])
        data.new <- data.new[phy$tip.label,]
        dat.tab <- OrganizeDataHiSSE(data=data.new, phy=phy, f=f, hidden.states=hidden.states, includes.intervals=FALSE, intervening.intervals=NULL)
        fossil.taxa <- NULL
        fix.type <- NULL
        psi.type <- NULL
        strat.cache <- NULL
    }
    nb.tip <- Ntip(phy)
    nb.node <- phy$Nnode
    ##########################
    
    #This is used to scale starting values to account for sampling:
    if(length(f) == 2){
        freqs <- table(apply(data.new, 1, function(x) switch(paste0(x, collapse=""), "00" = 1, "11" = 2, "22"=1)))
        ## Fixing the structure of the freqs vector.
        freqs.vec <- rep(0, times = 2)
        freqs.vec[as.numeric(names(freqs))] <- as.numeric( freqs )
        samp.freq.tree <- Ntip(phy) / sum(freqs.vec / f)
    }else{
        if(length(f) == Ntip(phy)){
            stop("This functionality has been temporarily removed.")
            #samp.freq.tree <- Ntip(phy) / sum(table(data.new[,1]) / mean(f[as.numeric(names(freqs))+1]))
        }else{
            stop("The vector of sampling frequencies does not match the number of tips in the tree.")
        }
    }
    if(is.null(restart.obj)){
        if(sum(eps)==0){
            if(includes.fossils == TRUE){
                stop("You input a tree that contains fossils but you are assuming no extinction. Check your function call.")
            }else{
                init.pars <- starting.point.generator(phy, 2, samp.freq.tree, yule=TRUE)
            }
        }else{
            if(includes.fossils == TRUE){
                if(!is.null(strat.intervals)){
                    branch.type = NULL
                    cols <- c("FocalNode","DesNode", "RootwardAge", "TipwardAge", "branch.type")
                    seg.map <- dat.tab[, cols, with=FALSE]
                    #remove k tips -- we do not do anything with them.
                    setkey(seg.map, branch.type)
                    #drop the k.tips because we do not do calculation on these zero length edges:
                    seg.map <- seg.map[branch.type != 2]
                    init.pars <- starting.point.generator.intervals(k=2, samp.freq.tree, n.tax=n.tax.starting, seg_map=seg.map, split.times=split.times, fossil.ages=fossil.ages, strat.cache=strat.cache)
                }else{
                    init.pars <- starting.point.generator.fossils(n.tax=n.tax.starting, k=2, samp.freq.tree, fossil.taxa=fossil.taxa, fossil.ages=fossil.ages, no.k.samples=no.k.samples, split.times=split.times)
                }
                psi.start <- init.pars[length(init.pars)]
            }else{
                init.pars <- starting.point.generator(phy, k=2, samp.freq.tree, yule=FALSE)
            }
            if(any(init.pars[3:4] == 0)){
                init.pars[3:4] = 1e-6
            }
        }
        names(init.pars) <- NULL
        
        if(is.null(starting.vals)){
            def.set.pars <- rep(c(log(init.pars[1:2]+init.pars[3:4]), log(init.pars[3:4]/init.pars[1:2]), rep(log(init.pars[5]),8)), rate.cats)
        }else{
            ## Check if 'starting.vals' has the correct format.
            if( !length(starting.vals) %in% c(3,12) ){
                stop("Wrong length of starting.vals")
            }
            if( length(starting.vals) == 5 ){
                cat("Using developer mode for starting.vals.", "\n")
                def.set.pars <- rep(c(log(starting.vals[1:2]), log(starting.vals[3:4]), rep(log(starting.vals[5]), 8)), rate.cats)
                
            } else{
                def.set.pars <- rep(c(log( rep(starting.vals[1],2) ), log( rep(starting.vals[2],2) ), log( rep(starting.vals[3],8))), rate.cats)
            }
        }
        
        if(bounded.search == TRUE){
            upper.full <- rep(c(rep(log(turnover.upper),2), rep(log(eps.upper),2), rep(log(trans.upper),8)), rate.cats)
        }else{
            upper.full <- rep(21,length(def.set.pars))
        }
        
        if(includes.fossils == TRUE){
            np.sequence <- 1:(np-1)
            ip <- numeric(np-1)
            upper <- numeric(np-1)
            for(i in np.sequence){
                ip[i] <- def.set.pars[which(pars == np.sequence[i])[1]]
                upper[i] <- upper.full[which(pars == np.sequence[i])[1]]
            }
            ip <- c(ip, log(psi.start))
            upper <- c(upper, log(trans.upper))
            lower <- rep(-20, length(ip))
        }else{
            np.sequence <- 1:(np)
            ip <- numeric(np)
            upper <- numeric(np)
            for(i in np.sequence){
                ip[i] <- def.set.pars[which(pars == np.sequence[i])[1]]
                upper[i] <- upper.full[which(pars == np.sequence[i])[1]]
            }
            lower <- rep(-20, length(ip))
        }
    }else{
        upper <- restart.obj$upper.bounds
        lower <- restart.obj$lower.bounds
        pars <- restart.obj$index.par
        ip <- numeric(length(unique(restart.obj$index.par))-1)
        for(k in 1:length(ip)){
            ip[k] <- log(restart.obj$solution[which(restart.obj$index.par==k)][1])
        }
    }
    
    if(sann == FALSE){
        if(bounded.search == TRUE){
            cat("Finished. Beginning bounded subplex routine...", "\n")
            opts <- list("algorithm" = "NLOPT_LN_SBPLX", "maxeval" = 100000, "ftol_rel" = max.tol)
            out = nloptr(x0=ip, eval_f=DevOptimizefHiSSE, ub=upper, lb=lower, opts=opts, pars=pars, dat.tab=dat.tab, gen=gen, hidden.states=hidden.states, nb.tip=nb.tip, nb.node=nb.node, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, np=np, f=f, ode.eps=ode.eps, fossil.taxa=fossil.taxa, fix.type=fix.type, strat.cache=strat.cache)
            sann.counts <- NULL
            solution <- numeric(length(pars))
            solution[] <- c(exp(out$solution), 0)[pars]
            loglik = -out$objective
        }else{
            cat("Finished. Beginning subplex routine...", "\n")
            out = subplex(ip, fn=DevOptimizefHiSSE, control=list(reltol=max.tol, parscale=rep(0.1, length(ip))), pars=pars, dat.tab=dat.tab, gen=gen, hidden.states=hidden.states, nb.tip=nb.tip, nb.node=nb.node, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, np=np, f=f, ode.eps=ode.eps, fossil.taxa=fossil.taxa, fix.type=fix.type, strat.cache=strat.cache)
            sann.counts <- NULL
            solution <- numeric(length(pars))
            solution[] <- c(exp(out$par), 0)[pars]
            loglik = -out$value
        }
    }else{
        cat("Finished. Beginning simulated annealing...", "\n")
        opts <- list("algorithm"="NLOPT_GD_STOGO", "maxeval" = 100000, "ftol_rel" = max.tol)
        out.sann = GenSA(ip, fn=DevOptimizefHiSSE, lower=lower, upper=upper, control=list(max.call=sann.its), pars=pars, dat.tab=dat.tab, gen=gen, hidden.states=hidden.states, nb.tip=nb.tip, nb.node=nb.node, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, np=np, f=f, ode.eps=ode.eps, fossil.taxa=fossil.taxa, fix.type=fix.type, strat.cache=strat.cache)
        sann.counts <- out.sann$counts
        cat("Finished. Refining using subplex routine...", "\n")
        opts <- list("algorithm" = "NLOPT_LN_SBPLX", "maxeval" = 100000, "ftol_rel" = max.tol)
        out <- nloptr(x0=out.sann$par, eval_f=DevOptimizefHiSSE, ub=upper, lb=lower, opts=opts, pars=pars, dat.tab=dat.tab, gen=gen, hidden.states=hidden.states, nb.tip=nb.tip, nb.node=nb.node, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, np=np, ode.eps=ode.eps, f=f, fossil.taxa=fossil.taxa, fix.type=fix.type, strat.cache=strat.cache)
        solution <- numeric(length(pars))
        solution[] <- c(exp(out$solution), 0)[pars]
        
        loglik = -out$objective
    }
    
    names(solution) <- c("turnover0A","turnover1A","eps0A","eps1A","q0A1A","q1A0A","q0A0B","q0A0C","q0A0D","q1A1B","q1A1C","q1A1D","turnover0B","turnover1B","eps0B","eps1B","q0B1B","q1B0B","q0B0A","q0B0C","q0B0D","q1B1A","q1B1C","q1B1D","turnover0C","turnover1C","eps0C","eps1C","q0C1C","q1C0C","q0C0A","q0C0B","q0C0D","q1C1A","q1C1B","q1C1D","turnover0D","turnover1D","eps0D","eps1D","q0D1D","q1D0D","q0D0A","q0D0B","q0D0C","q1D1A","q1D1B","q1D1C","psi")
    
    cat("Finished. Summarizing results...", "\n")
    
    obj = list(loglik = loglik, AIC = -2*loglik+2*np, AICc = -2*loglik+(2*np*(Ntip(phy)/(Ntip(phy)-np-1))), solution=solution, index.par=pars, f=f, hidden.states=hidden.states, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, phy=phy.og, phy.w.k=phy, data=data, trans.matrix=trans.rate, max.tol=max.tol, starting.vals=ip, upper.bounds=upper, lower.bounds=lower, ode.eps=ode.eps, includes.fossils=includes.fossils, k.samples=k.samples, strat.intervals=strat.intervals, fix.type=fix.type, psi.type=psi.type, sann.counts=sann.counts)
    class(obj) <- append(class(obj), "hisse.fit")
    return(obj)
    
    return(obj)
}


######################################################################################################################################
######################################################################################################################################
### The function used to optimize parameters:
######################################################################################################################################
######################################################################################################################################

#Function used for optimizing parameters:
DevOptimizefHiSSE <- function(p, pars, dat.tab, gen, hidden.states, nb.tip=nb.tip, nb.node=nb.node, condition.on.survival, root.type, root.p, np, f, ode.eps, fossil.taxa, fix.type, strat.cache) {
    #Generates the final vector with the appropriate parameter estimates in the right place:
    p.new <- exp(p)
    ## print(p.new)
    model.vec <- numeric(length(pars))
    model.vec[] <- c(p.new, 0)[pars]
    cache <- ParametersToPassfHiSSE(model.vec=model.vec, hidden.states=hidden.states, nb.tip=nb.tip, nb.node=nb.node, bad.likelihood=exp(-300), f=f, ode.eps=ode.eps)
    if(!is.null(fix.type)){
        if(!is.null(strat.cache)){
            logl <- DownPassHiSSE(dat.tab=dat.tab, cache=cache, gen=gen, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, node=fix.type[,1], state=NULL, fossil.taxa=fossil.taxa, fix.type=fix.type[,2]) + (strat.cache$k*log(cache$psi)) + (cache$psi*strat.cache$l_s)
        }else{
            logl <- DownPassHiSSE(dat.tab=dat.tab, cache=cache, gen=gen, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, node=fix.type[,1], state=NULL, fossil.taxa=fossil.taxa, fix.type=fix.type[,2])
        }
    }else{
        logl <- DownPassHiSSE(dat.tab=dat.tab, cache=cache, gen=gen, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, node=NULL, state=NULL, fossil.taxa=fossil.taxa, fix.type=NULL)
    }
    return(-logl)
}


######################################################################################################################################
######################################################################################################################################
### The various utility functions used
######################################################################################################################################
######################################################################################################################################

OrganizeDataHiSSE <- function(data, phy, f, hidden.states, includes.intervals=FALSE, intervening.intervals=NULL){
    ### Ughy McUgherson. This is a must in order to pass CRAN checks: http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
    DesNode = NULL
    
    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
    
    if(hidden.states == FALSE){
        states = matrix(0,Ntip(phy),2)
        x <- data[,1]
        for(i in 1:Ntip(phy)){
            if(x[i]==0){states[i,1]=1}
            if(x[i]==1){states[i,2]=1}
            if(x[i]==2){states[i,c(1,2)]=1}
        }
        compD <- matrix(0, nrow=nb.tip, ncol=2)
        compE <- matrix(0, nrow=nb.tip, ncol=2)
    }
    
    if(hidden.states == TRUE){
        states = matrix(0,Ntip(phy),8)
        x <- data[,1]
        for(i in 1:Ntip(phy)){
            if(x[i]==0){states[i,c(1,3,5,7)]=1}
            if(x[i]==1){states[i,c(2,4,6,8)]=1}
            if(x[i]==2){states[i,c(1,3,5,7, 2,4,6,8)]=1}
        }
        compD <- matrix(0, nrow=nb.tip, ncol=8)
        compE <- matrix(0, nrow=nb.tip, ncol=8)
    }
    
    #Initializes the tip sampling and sets internal nodes to be zero:
    ncols = dim(compD)[2]
    if(length(f) == 2){
        for(i in 1:(nb.tip)){
            compD[i,] <- f * states[i,]
            compE[i,] <- rep((1-f), ncols/2)
        }
    }else{
        for(i in 1:(nb.tip)){
            compD[i,] <- f[i] * states[i,]
            compE[i,] <- rep((1-f[i]), ncols/2)
        }
    }
    
    #This seems stupid but I cannot figure out how to get data.table to not make this column a factor. When a factor this is not right. For posterity, let it be known Jeremy would rather just retain the character instead of this mess:
    edge_details <- GetEdgeDetails(phy, includes.intervals=includes.intervals, intervening.intervals=intervening.intervals)
    branch.type <- edge_details$type
    branch.type[which(branch.type == "extant_tip")] <- 0
    branch.type[which(branch.type == "internal")] <- 0
    branch.type[which(branch.type == "extinct_tip")] <- 1
    branch.type[which(branch.type == "k_tip")] <- 2
    branch.type[which(branch.type == "k_k_interval")] <- 3
    branch.type[which(branch.type == "k_extinct_interval")] <- 3
    branch.type[which(branch.type == "k_extant_interval")] <- 3
    branch.type[which(branch.type == "intervening_interval")] <- 4

    tmp.df <- cbind(edge_details[,1:5], 0, matrix(0, nrow(edge_details), ncol(compD)), matrix(0, nrow(edge_details), ncol(compE)), as.numeric(branch.type))
    colnames(tmp.df) <- c("RootwardAge", "TipwardAge", "BranchLength", "FocalNode", "DesNode", "comp", paste("compD", 1:ncol(compD), sep="_"), paste("compE", 1:ncol(compE), sep="_"), "branch.type")
    dat.tab <- as.data.table(tmp.df)
    setkey(dat.tab, DesNode)
    cols <- names(dat.tab)
    for (j in 1:(dim(compD)[2])){
        dat.tab[data.table(c(1:nb.tip)), paste("compD", j, sep="_") := compD[,j]]
        #set(dat.tab, 1:nb.tip, cols[6+j], compD[,j])
        dat.tab[data.table(c(1:nb.tip)), paste("compE", j, sep="_") := compE[,j]]
        #set(dat.tab, 1:nb.tip, cols[38+j], compE[,j])
    }
    return(dat.tab)
}


SingleChildProbHiSSE <- function(cache, pars, compD, compE, start.time, end.time, branch.type){
    if(any(!is.finite(c(compD, compE)))) { # something went awry at a previous step. Bail!
        prob.subtree.cal <- rep(0, ifelse(cache$hidden.states == TRUE, 16, 4))
        prob.subtree.cal[(1+length(prob.subtree.cal)/2):length(prob.subtree.cal)] <- cache$bad.likelihood
        return(prob.subtree.cal)
    }
    if(cache$hidden.states == TRUE){
        yini <-c(E0A=compE[1], E1A=compE[2], E0B=compE[3], E1B=compE[4], E0C=compE[5], E1C=compE[6], E0D=compE[7], E1D=compE[8], D0A=compD[1], D1A=compD[2], D0B=compD[3], D1B=compD[4], D0C=compD[5], D1C=compD[6], D0D=compD[7], D1D=compD[8])
        if(branch.type == 2){
            times <- c(0, end.time)
        }else{
            times=c(start.time, end.time)
        }
        runSilent <- function(branch.type) {
            options(warn = -1)
            on.exit(options(warn = 0))
            if(branch.type == 3){
                capture.output(res <- lsoda(yini, times, func = "maddison_DE_strat_fhisse", pars, initfunc="initmod_fbisse", dllname = "hisse", rtol=1e-8, atol=1e-8))
            }else{
                capture.output(res <- lsoda(yini, times, func = "maddison_DE_fhisse", pars, initfunc="initmod_fhisse", dllname = "hisse", rtol=1e-8, atol=1e-8))
            }
            res
        }
        #prob.subtree.cal.full <- lsoda(yini, times, func = "maddison_DE_fhisse", pars, initfunc="initmod_fhisse", dll = "fhisse-ext-derivs", rtol=1e-8, atol=1e-8)
        #prob.subtree.cal.full <- lsoda(yini, times, func = "maddison_DE_fhisse", pars, initfunc="initmod_fhisse", dllname = "hisse", rtol=1e-8, atol=1e-8)
        prob.subtree.cal.full <- runSilent(branch.type=branch.type)
    }else{
        yini <-c(E0=compE[1], E1=compE[2], D0=compD[1], D1=compD[2])
        if(branch.type == 2){
            times <- c(0, end.time)
        }else{
            times=c(start.time, end.time)
        }
        runSilent <- function(branch.type) {
            options(warn = -1)
            on.exit(options(warn = 0))
            if(branch.type == 3){
                capture.output(res <- lsoda(yini, times, func = "maddison_DE_strat_fbisse", pars, initfunc="initmod_fbisse", dllname = "hisse", rtol=1e-8, atol=1e-8))
            }else{
                capture.output(res <- lsoda(yini, times, func = "maddison_DE_fbisse", pars, initfunc="initmod_fbisse", dllname = "hisse", rtol=1e-8, atol=1e-8))
            }
            res
        }
        #prob.subtree.cal.full <- lsoda(yini, times, func = "maddison_DE_fbisse", pars, initfunc="initmod_fbisse", dll = "fbisse-ext-derivs", rtol=1e-8, atol=1e-8)
        #prob.subtree.cal.full <- lsoda(yini, times, func = "maddison_DE_fbisse", pars, initfunc="initmod_fbisse", dllname = "hisse", rtol=1e-8, atol=1e-8)
        prob.subtree.cal.full <- runSilent(branch.type=branch.type)
    }
    ######## THIS CHECKS TO ENSURE THAT THE INTEGRATION WAS SUCCESSFUL ###########
    if(attributes(prob.subtree.cal.full)$istate[1] < 0){
        prob.subtree.cal <- prob.subtree.cal.full[-1,-1]
        if(cache$hidden.states == TRUE){
            prob.subtree.cal[9:16] <- cache$bad.likelihood
            return(prob.subtree.cal)
        }else{
            prob.subtree.cal[3:4] <- cache$bad.likelihood
            return(prob.subtree.cal)
        }
    }else{
        if(branch.type == 2){
            prob.subtree.cal <- prob.subtree.cal.full[-1,-1]
            if(cache$hidden.states == TRUE){
                prob.subtree.cal[9:16] <- compD
            }else{
                prob.subtree.cal[3:4] <- compD
            }
        }else{
            if(branch.type == 1){
                if(cache$hidden.states == TRUE){
                    prob.subtree.cal <- prob.subtree.cal.full[-1,-1]
                    prob.subtree.cal[9:16] <- prob.subtree.cal[1:8] * compD
                }else{
                    prob.subtree.cal <- prob.subtree.cal.full[-1,-1]
                    prob.subtree.cal[3:4] <- prob.subtree.cal[1:2] * compD
                }
            }else{
                if(branch.type == 4){
                    if(cache$hidden.states == TRUE){
                        prob.subtree.cal.wevents <- prob.subtree.cal.full[-1,-1]
                        prob.subtree.cal.full <- runSilent(branch.type=3)
                        prob.subtree.cal.noevents <- prob.subtree.cal.full[-1,-1]
                        unobs.spec.probs <- 1 - (prob.subtree.cal.noevents[9:16]/prob.subtree.cal.wevents[9:16])
                        prob.subtree.cal.wevents[9:16] <- prob.subtree.cal.wevents[9:16] * unobs.spec.probs
                        prob.subtree.cal.wevents[which(is.na(prob.subtree.cal.wevents))] <- 0
                        prob.subtree.cal <- prob.subtree.cal.wevents
                    }else{
                        prob.subtree.cal.wevents <- prob.subtree.cal.full[-1,-1]
                        prob.subtree.cal.full <- runSilent(branch.type=3)
                        prob.subtree.cal.noevents <- prob.subtree.cal.full[-1,-1]
                        unobs.spec.probs <- 1 - (prob.subtree.cal.noevents[3:4]/prob.subtree.cal.wevents[3:4])
                        prob.subtree.cal.wevents[3:4] <- prob.subtree.cal.wevents[3:4] * unobs.spec.probs
                        prob.subtree.cal.wevents[which(is.na(prob.subtree.cal.wevents))] <- 0
                        prob.subtree.cal <- prob.subtree.cal.wevents
                    }
                }else{
                    prob.subtree.cal <- prob.subtree.cal.full[-1,-1]
                }
            }
        }
    }
    ##############################################################################
    
    if(cache$hidden.states == TRUE){
        if(any(is.nan(prob.subtree.cal[9:16]))){
            prob.subtree.cal[9:16] <- cache$bad.likelihood
            return(prob.subtree.cal)
        }
        #This is default and cannot change, but if we get a negative probability, discard the results:
        if(any(prob.subtree.cal[9:16] < 0)){
            prob.subtree.cal[9:16] <- cache$bad.likelihood
            return(prob.subtree.cal)
        }
        if(sum(prob.subtree.cal[9:16]) < cache$ode.eps){
            prob.subtree.cal[9:16] <- cache$bad.likelihood
            return(prob.subtree.cal)
        }
    }else{
        if(any(is.nan(prob.subtree.cal[3:4]))){
            prob.subtree.cal[3:4] <- cache$bad.likelihood
            return(prob.subtree.cal)
        }
        #This is default and cannot change, but if we get a negative probability, discard the results:
        if(any(prob.subtree.cal[3:4] < 0)){
            prob.subtree.cal[3:4] <- cache$bad.likelihood
            return(prob.subtree.cal)
        }
        if(sum(prob.subtree.cal[3:4]) < cache$ode.eps){
            prob.subtree.cal[3:4] <- cache$bad.likelihood
            return(prob.subtree.cal)
        }
    }
    return(prob.subtree.cal)
}


FocalNodeProbHiSSE <- function(cache, pars, lambdas, dat.tab, generations){
    ### Ughy McUgherson. This is a must in order to pass CRAN checks: http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
    DesNode = NULL
    FocalNode = NULL
    . = NULL
    
    gens <- data.table(c(generations))
    setkey(dat.tab, FocalNode)
    CurrentGenData <- dat.tab[gens]
    if(cache$hidden.states == TRUE){
        tmp <- t(apply(CurrentGenData, 1, function(z) SingleChildProbHiSSE(cache, pars, z[7:14], z[15:22],  z[2], z[1], z[23])))
        v.mat <- matrix(tmp[seq(1,nrow(tmp)-1,2),9:16] * tmp[seq(2,nrow(tmp),2),9:16], length(unique(CurrentGenData$FocalNode)), 8)
        v.mat <- v.mat * matrix(lambdas, length(unique(CurrentGenData$FocalNode)), 8, byrow=TRUE)
        phi.mat <- matrix(tmp[seq(1,nrow(tmp)-1,2),1:8], length(unique(CurrentGenData$FocalNode)), 8)
        if(!is.null(cache$node)){
            if(any(cache$node %in% generations)){
                for(fix.index in 1:length(cache$node)){
                    if(cache$fix.type[fix.index] == "event"){
                        #basically we are using the node to fix the state along a branch, but we do not want to assume a true speciation event occurred here.
                        lambdas.check <- lambdas
                        lambdas.check[which(lambdas==0)] <- 1
                        #The initial condition for a k.sample is D(t)*psi
                        v.mat[which(generations == cache$node[fix.index]),] <- (v.mat[which(generations == cache$node[fix.index]),] / lambdas.check) * cache$psi
                    }else{
                        if(cache$fix.type[fix.index] == "interval"){
                            lambdas.check <- lambdas
                            lambdas.check[which(lambdas==0)] <- 1
                            v.mat[which(generations == cache$node[fix.index]),] <- v.mat[which(generations == cache$node[fix.index]),] / lambdas.check
                        }else{
                            fixer = numeric(8)
                            fixer[cache$state[fix.index]] = 1
                            v.mat[which(generations == cache$node[fix.index]),] <- v.mat[which(generations == cache$node[fix.index]),] * fixer
                        }
                    }
                }
            }
        }
        tmp.comp <- rowSums(v.mat)
        tmp.probs <- v.mat / tmp.comp
        setkey(dat.tab, DesNode)
        rows <- dat.tab[.(generations), which=TRUE]
        cols <- names(dat.tab)
        for (j in 1:(dim(tmp.probs)[2])){
            #dat.tab[gens, cols[6+j] := tmp.probs[,j]]
            set(dat.tab, rows, cols[6+j], tmp.probs[,j])
            #dat.tab[gens, cols[38+j] := phi.mat[,j]]
            set(dat.tab, rows, cols[14+j], phi.mat[,j])
        }
        dat.tab[gens, "comp" := tmp.comp]
    }else{
        tmp <- t(apply(CurrentGenData, 1, function(z) SingleChildProbHiSSE(cache, pars, z[7:8], z[9:10], z[2], z[1], z[11])))
        v.mat <- matrix(tmp[seq(1,nrow(tmp)-1,2),3:4] * tmp[seq(2,nrow(tmp),2),3:4], length(unique(CurrentGenData$FocalNode)), 2)
        v.mat <- v.mat * matrix(lambdas, length(unique(CurrentGenData$FocalNode)), 2, byrow=TRUE)
        phi.mat <- matrix(tmp[seq(1,nrow(tmp)-1,2),1:2], length(unique(CurrentGenData$FocalNode)), 2)
        if(!is.null(cache$node)){
            if(any(cache$node %in% generations)){
                for(fix.index in 1:length(cache$node)){
                    if(cache$fix.type[fix.index] == "event"){
                        #basically we are using the node to fix the state along a branch, but we do not want to assume a true speciation event occurred here.
                        lambdas.check <- lambdas
                        lambdas.check[which(lambdas==0)] <- 1
                        #The initial condition for a k.sample is D(t)*psi
                        v.mat[which(generations == cache$node[fix.index]),] <- (v.mat[which(generations == cache$node[fix.index]),] / lambdas.check) * cache$psi
                    }else{
                        if(cache$fix.type[fix.index] == "interval"){
                            lambdas.check <- lambdas
                            lambdas.check[which(lambdas==0)] <- 1
                            v.mat[which(generations == cache$node[fix.index]),] <- v.mat[which(generations == cache$node[fix.index]),] / lambdas.check
                        }else{
                            fixer = numeric(2)
                            fixer[cache$state[fix.index]] = 1
                            v.mat[which(generations == cache$node[fix.index]),] <- v.mat[which(generations == cache$node[fix.index]),] * fixer
                        }
                    }
                }
            }
        }
        tmp.comp <- rowSums(v.mat)
        tmp.probs <- v.mat / tmp.comp
        setkey(dat.tab, DesNode)
        #gens <- data.table(c(generations))
        rows <- dat.tab[.(generations), which=TRUE]
        cols <- names(dat.tab)
        for (j in 1:(dim(tmp.probs)[2])){
            #dat.tab[gens, cols[6+j] := tmp.probs[,j]]
            set(dat.tab, rows, cols[6+j], tmp.probs[,j])
            #dat.tab[gens, cols[10+j] := phi.mat[,j]]
            set(dat.tab, rows, cols[8+j], phi.mat[,j])
        }
        dat.tab[gens, "comp" := tmp.comp]
    }
    return(dat.tab)
}


#Have to calculate root prob separately because it is not a descendant in our table. Could add it, but I worry about the NA that is required.
GetRootProbHiSSE <- function(cache, pars, lambdas, dat.tab, generations){
    ### Ughy McUgherson. This is a must in order to pass CRAN checks: http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
    FocalNode = NULL
    . = NULL
    
    gens <- data.table(c(generations))
    setkey(dat.tab, FocalNode)
    CurrentGenData <- dat.tab[gens]
    if(cache$hidden.states == TRUE){
        tmp <- t(apply(CurrentGenData, 1, function(z) SingleChildProbHiSSE(cache, pars, z[7:14], z[15:22], z[2], z[1], z[23])))
        v.mat <- matrix(tmp[seq(1,nrow(tmp)-1,2),9:16] * tmp[seq(2,nrow(tmp),2),9:16], length(unique(CurrentGenData$FocalNode)), 8)
        v.mat <- v.mat * matrix(lambdas, length(unique(CurrentGenData$FocalNode)), 8, byrow=TRUE)
        phi.mat <- matrix(tmp[seq(1,nrow(tmp)-1,2),1:8], length(unique(CurrentGenData$FocalNode)), 8)
        if(!is.null(cache$node)){
            if(any(cache$node %in% generations)){
                for(fix.index in 1:length(cache$node)){
                    if(cache$fix.type[fix.index] == "event"){
                        #basically we are using the node to fix the state along a branch, but we do not want to assume a true speciation event occurred here.
                        lambdas.check <- lambdas
                        lambdas.check[which(lambdas==0)] <- 1
                        #The initial condition for a k.sample is D(t)*psi
                        v.mat[which(generations == cache$node[fix.index]),] <- (v.mat[which(generations == cache$node[fix.index]),] / lambdas.check) * cache$psi
                    }else{
                        if(cache$fix.type[fix.index] == "interval"){
                            lambdas.check <- lambdas
                            lambdas.check[which(lambdas==0)] <- 1
                            v.mat[which(generations == cache$node[fix.index]),] <- v.mat[which(generations == cache$node[fix.index]),] / lambdas.check
                        }else{
                            fixer = numeric(8)
                            fixer[cache$state[fix.index]] = 1
                            v.mat[which(generations == cache$node[fix.index]),] <- v.mat[which(generations == cache$node[fix.index]),] * fixer
                        }
                    }
                }
            }
        }
    }else{
        tmp <- t(apply(CurrentGenData, 1, function(z) SingleChildProbHiSSE(cache, pars, z[7:8], z[9:10], z[2], z[1], z[11])))
        v.mat <- matrix(tmp[seq(1,nrow(tmp)-1,2),3:4] * tmp[seq(2,nrow(tmp),2),3:4], length(unique(CurrentGenData$FocalNode)), 2)
        v.mat <- v.mat * matrix(lambdas, length(unique(CurrentGenData$FocalNode)), 2, byrow=TRUE)
        phi.mat <- matrix(tmp[seq(1,nrow(tmp)-1,2),1:2], length(unique(CurrentGenData$FocalNode)), 2)
        if(!is.null(cache$node)){
            if(any(cache$node %in% generations)){
                for(fix.index in 1:length(cache$node)){
                    if(cache$fix.type[fix.index] == "event"){
                        #basically we are using the node to fix the state along a branch, but we do not want to assume a true speciation event occurred here.
                        lambdas.check <- lambdas
                        lambdas.check[which(lambdas==0)] <- 1
                        #The initial condition for a k.sample is D(t)*psi
                        v.mat[which(generations == cache$node[fix.index]),] <- (v.mat[which(generations == cache$node[fix.index]),] / lambdas.check) * cache$psi
                    }else{
                        if(cache$fix.type[fix.index] == "interval"){
                            lambdas.check <- lambdas
                            lambdas.check[which(lambdas==0)] <- 1
                            v.mat[which(generations == cache$node[fix.index]),] <- v.mat[which(generations == cache$node[fix.index]),] / lambdas.check
                        }else{
                            fixer = numeric(2)
                            fixer[cache$state[fix.index]] = 1
                            v.mat[which(generations == cache$node[fix.index]),] <- v.mat[which(generations == cache$node[fix.index]),] * fixer
                        }
                    }
                }
            }
        }
    }
    
    tmp.comp <- rowSums(v.mat)
    tmp.probs <- v.mat / tmp.comp
    
    return(cbind(tmp.comp, phi.mat, tmp.probs))
}


#Calculates the initial conditions for fossil taxa in the tree.
GetFossilInitialsHiSSE <- function(cache, pars, lambdas, dat.tab, fossil.taxa){
    ### Ughy McUgherson. This is a must in order to pass CRAN checks: http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
    DesNode = NULL
    . = NULL
    if(cache$hidden.states == TRUE){
        fossils <- data.table(c(fossil.taxa))
        setkey(dat.tab, DesNode)
        CurrentGenData <- dat.tab[fossils]
        tmp <- t(apply(CurrentGenData, 1, function(z) SingleChildProbHiSSE(cache, pars, z[7:14], z[15:22], 0, z[2], 1)))
        tmp.probs <- matrix(tmp[,9:16], length(fossil.taxa), 8) * cache$psi
        phi.mat <- matrix(tmp[,1:8], length(fossil.taxa), 8)
        setkey(dat.tab, DesNode)
        rows <- dat.tab[.(fossils), which=TRUE]
        cols <- names(dat.tab)
        for (j in 1:(dim(tmp.probs)[2])){
            #dat.tab[data.table(c(generations)), paste("compD", j, sep="_") := tmp.probs[,j]]
            set(dat.tab, rows, cols[6+j], tmp.probs[,j])
            #dat.tab[data.table(c(generations)), paste("compE", j, sep="_") := phi.mat[,j]]
            set(dat.tab, rows, cols[14+j], phi.mat[,j])
        }
    }else{
        fossils <- data.table(c(fossil.taxa))
        setkey(dat.tab, DesNode)
        CurrentGenData <- dat.tab[fossils]
        tmp <- t(apply(CurrentGenData, 1, function(z) SingleChildProbHiSSE(cache, pars, z[7:8], z[9:10], 0, z[2], 1)))
        tmp.probs <- matrix(tmp[,3:4], length(fossil.taxa), 2) * cache$psi
        phi.mat <- matrix(tmp[,1:2], length(fossil.taxa), 2)
        
        setkey(dat.tab, DesNode)
        rows <- dat.tab[.(fossils), which=TRUE]
        cols <- names(dat.tab)
        for (j in 1:(dim(tmp.probs)[2])){
            #dat.tab[data.table(c(generations)), paste("compD", j, sep="_") := tmp.probs[,j]]
            set(dat.tab, rows, cols[6+j], tmp.probs[,j])
            #dat.tab[gens, cols[10+j] := phi.mat[,j]]
            set(dat.tab, rows, cols[8+j], phi.mat[,j])
        }
    }
    dat.tab[fossils,"branch.type" := 0]
    return(dat.tab)
}



######################################################################################################################################
######################################################################################################################################
### The fast HiSSE type down pass that carries out the integration and returns the likelihood:
######################################################################################################################################
######################################################################################################################################

DownPassHiSSE <- function(dat.tab, gen, cache, condition.on.survival, root.type, root.p, get.phi=FALSE, node=NULL, state=NULL, fossil.taxa=NULL, fix.type=NULL) {
    
    ### Ughy McUgherson. This is a must in order to pass CRAN checks: http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
    DesNode = NULL
    compE = NULL
    . = NULL
    
    ## Moved this here instead of above because it actually significantly improved the speed.
    if(cache$hidden.states == FALSE){
        pars <- c(cache$lambda0A,cache$lambda1A,cache$mu0A,cache$mu1A,cache$q0A1A,cache$q1A0A,cache$psi)
        lambda <- c(cache$lambda0A, cache$lambda1A)
    }else{
        pars <- c(cache$lambda0A,cache$lambda1A,cache$mu0A,cache$mu1A,cache$q0A1A,cache$q1A0A,cache$q0A0B,cache$q0A0C,cache$q0A0D,cache$q1A1B,cache$q1A1C,cache$q1A1D,cache$lambda0B,cache$lambda1B ,cache$mu0B,cache$mu1B,cache$q0B1B,cache$q1B0B,cache$q0B0A,cache$q0B0C,cache$q0B0D,cache$q1B1A,cache$q1B1C,cache$q1B1D,cache$lambda0C ,cache$lambda1C ,cache$mu0C,cache$mu1C,cache$q0C1C,cache$q1C0C,cache$q0C0A,cache$q0C0B,cache$q0C0D,cache$q1C1A,cache$q1C1B,cache$q1C1D,cache$lambda0D,cache$lambda1D,cache$mu0D,cache$mu1D,cache$q0D1D,cache$q1D0D,cache$q0D0A,cache$q0D0B,cache$q0D0C,cache$q1D1A,cache$q1D1B,cache$q1D1C,cache$psi)
        lambda <- c(cache$lambda0A,cache$lambda1A,cache$lambda0B,cache$lambda1B, cache$lambda0C,cache$lambda1C, cache$lambda0D,cache$lambda1D)
    }
    
    if(!is.null(fossil.taxa)){
        #Gets initial conditions for fossil taxa:
        dat.tab.copy <- copy(dat.tab)
        dat.tab.copy <- GetFossilInitialsHiSSE(cache=cache, pars=pars, lambdas=lambda, dat.tab=dat.tab.copy, fossil.taxa=fossil.taxa)
    }else{
        dat.tab.copy <- copy(dat.tab)
    }
    
    TIPS <- 1:cache$nb.tip
    for(i in 1:length(gen)){
        if(i == length(gen)){
            if(!is.null(node)){
                if(any(node %in% gen[[i]])){
                    cache$node <- node
                    cache$state <- state
                    cache$fix.type <- fix.type
                    res.tmp <- GetRootProbHiSSE(cache=cache, pars=pars, lambdas=lambda, dat.tab=dat.tab.copy, generations=gen[[i]])
                    cache$node <- NULL
                    cache$state <- NULL
                    cache$fix.type <- NULL
                }else{
                    res.tmp <- GetRootProbHiSSE(cache=cache, pars=pars, lambdas=lambda, dat.tab=dat.tab.copy, generations=gen[[i]])
                }
            }else{
                res.tmp <- GetRootProbHiSSE(cache=cache, pars=pars, lambdas=lambda, dat.tab=dat.tab.copy, generations=gen[[i]])
            }
            if(cache$hidden.states == TRUE){
                compD.root <- res.tmp[c(10:17)]
                compE.root <- res.tmp[c(2:9)]
            }else{
                compD.root <- res.tmp[c(4:5)]
                compE.root <- res.tmp[c(2:3)]
            }
            setkey(dat.tab.copy, DesNode)
            comp <- dat.tab.copy[["comp"]]
            comp <- c(comp[-TIPS], res.tmp[1])
        }else{
            if(!is.null(node)){
                if(any(node %in% gen[[i]])){
                    cache$node <- node
                    cache$state <- state
                    cache$fix.type <- fix.type
                    dat.tab.copy <- FocalNodeProbHiSSE(cache, pars=pars, lambdas=lambda, dat.tab.copy, gen[[i]])
                    cache$node <- NULL
                    cache$state <- NULL
                    cache$fix.type <- NULL
                }else{
                    dat.tab.copy <- FocalNodeProbHiSSE(cache, pars=pars, lambdas=lambda, dat.tab.copy, gen[[i]])
                }
            }else{
                dat.tab.copy <- FocalNodeProbHiSSE(cache, pars=pars, lambdas=lambda, dat.tab.copy, gen[[i]])
            }
        }
    }
    
    if(!is.null(fossil.taxa)){
        if(all(fix.type == "event")){
            if(cache$hidden.states == TRUE){
                pars[length(pars)] <- 0
                cache$psi <- 0
                init.d <- rep(cache$f, 4)
                init.e <- rep(1-cache$f, 4)
                phi.mat <- SingleChildProbHiSSE(cache, pars, init.d, init.e, 0, max(dat.tab$RootwardAge), 0)
                compE.root <- matrix(phi.mat[1:8], 1, 8)
            }else{
                pars[length(pars)] <- 0
                cache$psi <- 0
                init.d <- cache$f
                init.e <- 1-cache$f
                phi.mat <- SingleChildProbHiSSE(cache, pars, init.d, init.e, 0, max(dat.tab$RootwardAge), 0)
                compE.root <- matrix(phi.mat[1:2], 1, 2)
            }
        }
    }
    if(is.na(sum(log(compD.root))) || is.na(log(sum(1-compE.root)))){
        return(log(cache$bad.likelihood)^13)
    }else{
        if(root.type == "madfitz" | root.type == "herr_als"){
            if(is.null(root.p)){
                root.p = compD.root/sum(compD.root)
                root.p[which(is.na(root.p))] = 0
            }
        }
        if(condition.on.survival == TRUE){
            if(cache$hidden.states == FALSE){
                if(root.type == "madfitz"){
                    compD.root <- compD.root / sum(root.p * lambda * (1 - compE.root)^2)
                    #Corrects for possibility that you have 0/0:
                    compD.root[which(is.na(compD.root))] = 0
                    loglik <- log(sum(compD.root * root.p)) + sum(log(comp))
                }else{
                    compD.root <- (compD.root * root.p) / (lambda * (1 - compE.root)^2)
                    #Corrects for possibility that you have 0/0:
                    compD.root[which(is.na(compD.root))] = 0
                    loglik <- log(sum(compD.root)) + sum(log(comp))
                }
            }else{
                if(root.type == "madfitz"){
                    compD.root <- compD.root / sum(root.p * lambda * (1 - compE.root)^2)
                    #Corrects for possibility that you have 0/0:
                    compD.root[which(is.na(compD.root))] = 0
                    loglik <- log(sum(compD.root * root.p)) + sum(log(comp))
                }else{
                    compD.root <- (compD.root * root.p) / (lambda * (1 - compE.root)^2)
                    #Corrects for possibility that you have 0/0:
                    compD.root[which(is.na(compD.root))] = 0
                    loglik <- log(sum(compD.root)) + sum(log(comp))
                }
            }
        }
        if(!is.finite(loglik)){
            return(log(cache$bad.likelihood)^13)
        }
    }
    if(get.phi==TRUE){
        obj = NULL
        obj$compD.root = compD.root/sum(compD.root)
        obj$compE = compE.root
        obj$root.p = root.p
        return(obj)
    }else{
        return(loglik)
    }
}


######################################################################################################################################
######################################################################################################################################
### Cache object for storing parameters that are used throughout MuHiSSE:
######################################################################################################################################
######################################################################################################################################

ParametersToPassfHiSSE <- function(model.vec, hidden.states, nb.tip, nb.node, bad.likelihood, f, ode.eps){
    #Provides an initial object that contains all the parameters to be passed among functions. This will also be used to pass other things are we move down the tree (see DownPassGeoSSE):
    obj <- NULL
    
    obj$hidden.states <- hidden.states
    obj$nb.tip <- nb.tip
    obj$nb.node <- nb.node
    obj$bad.likelihood <- bad.likelihood
    obj$ode.eps <- ode.eps
    obj$f <- f
    
    ##Hidden State A
    obj$lambda0A <- model.vec[1] / (1 + model.vec[3])
    obj$lambda1A <- model.vec[2] / (1 + model.vec[4])
    obj$mu0A <- (model.vec[3] * model.vec[1]) / (1 + model.vec[3])
    obj$mu1A <- (model.vec[4] * model.vec[2]) / (1 + model.vec[4])
    obj$q0A1A <- model.vec[5]
    obj$q1A0A <- model.vec[6]
    
    obj$q0A0B <- model.vec[7]
    obj$q0A0C <- model.vec[8]
    obj$q0A0D <- model.vec[9]
    obj$q1A1B <- model.vec[10]
    obj$q1A1C <- model.vec[11]
    obj$q1A1D <- model.vec[12]
    
    ##Hidden State B
    obj$lambda0B <- model.vec[13] / (1 + model.vec[15])
    obj$lambda1B <- model.vec[14] / (1 + model.vec[16])
    obj$mu0B <- (model.vec[15] * model.vec[13]) / (1 + model.vec[15])
    obj$mu1B <- (model.vec[16] * model.vec[14]) / (1 + model.vec[16])
    obj$q0B1B <- model.vec[17]
    obj$q1B0B <- model.vec[18]
    
    obj$q0B0A <- model.vec[19]
    obj$q0B0C <- model.vec[20]
    obj$q0B0D <- model.vec[21]
    obj$q1B1A <- model.vec[22]
    obj$q1B1C <- model.vec[23]
    obj$q1B1D <- model.vec[24]
    
    ##Hidden State C
    obj$lambda0C <- model.vec[25] / (1 + model.vec[27])
    obj$lambda1C <- model.vec[26] / (1 + model.vec[28])
    obj$mu0C <- (model.vec[27] * model.vec[25]) / (1 + model.vec[27])
    obj$mu1C <- (model.vec[28] * model.vec[26]) / (1 + model.vec[28])
    obj$q0C1C <- model.vec[29]
    obj$q1C0C <- model.vec[30]
    obj$q0C0A <- model.vec[31]
    
    obj$q0C0B <- model.vec[32]
    obj$q0C0D <- model.vec[33]
    obj$q1C1A <- model.vec[34]
    obj$q1C1B <- model.vec[35]
    obj$q1C1D <- model.vec[36]
    
    ##Hidden State D
    obj$lambda0D <- model.vec[37] / (1 + model.vec[39])
    obj$lambda1D <- model.vec[38] / (1 + model.vec[40])
    obj$mu0D <- (model.vec[39] * model.vec[37]) / (1 + model.vec[39])
    obj$mu1D <- (model.vec[40] * model.vec[38]) / (1 + model.vec[40])
    obj$q0D1D <- model.vec[41]
    obj$q1D0D <- model.vec[42]
    
    obj$q0D0A <- model.vec[43]
    obj$q0D0B <- model.vec[44]
    obj$q0D0C <- model.vec[45]
    obj$q1D1A <- model.vec[46]
    obj$q1D1B <- model.vec[47]
    obj$q1D1C <- model.vec[48]
    
    obj$psi <- model.vec[49]
    
    return(obj)
}


print.hisse.fit <- function(x,...){
    ## Function to print a "muhisse.fit" object.
    set.zero <- max( x$index.par )
    ## Keep only the parameters estimated:
    par.list <- x$solution[ !x$index.par == set.zero ]
    ntips <- Ntip( x$phy )
    nstates <- ncol( x$trans.matrix )/2
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
    if(!is.null(x$psi.type)){
        if(x$psi.type == "m_only"){
            cat("psi estimate reflects fossil tip sampling only. \n")
        }else{
            cat("psi estimate reflects both fossil edge and tip sampling. \n")
        }
    }
}

