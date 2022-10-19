

######################################################################################################################################
######################################################################################################################################
### MiSSE -- Examines shifts in diversification in relation to hidden states ONLY
######################################################################################################################################
######################################################################################################################################

MiSSE <- function(phy, f=1, turnover=c(1,2), eps=c(1,2), fixed.eps=NULL, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, includes.fossils=FALSE, k.samples=NULL, strat.intervals=NULL, sann=TRUE, sann.its=5000, sann.temp=5230, sann.seed=-100377, bounded.search=TRUE, max.tol=.Machine$double.eps^.50, starting.vals=NULL, turnover.upper=10000, eps.upper=3, trans.upper=100, restart.obj=NULL, ode.eps=0, dt.threads=1, expand.mode=FALSE){
    
    misse_start_time <- Sys.time()
    
    #This makes it easier to handle missegreedy with fixed values
    if(length(fixed.eps)>0) {
        if(is.na(fixed.eps)) {
            fixed.eps <- NULL
        }else{
            eps <- numeric(length(turnover))
        }
    }
    if(expand.mode) {
        if(length(turnover)==1) {
            turnover <- sequence(turnover)
        }
        
        if(length(eps)==1 & eps[1]!=0) {
            eps <- sequence(eps)
        }
        if(length(eps)<length(turnover)) {
            eps <- rep(eps, length(turnover))[sequence(length(turnover))] # i.e., c(1,2,1,2,1)
        }
        if(length(turnover)<length(eps)) {
            turnover <- rep(turnover, length(eps))[sequence(length(eps))]
        }
    }
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
    
    setDTthreads(threads=dt.threads)
    
    #if(!is.ultrametric(phy) & includes.fossils == FALSE){
    #    warning("Tree is not ultrametric. Used force.ultrametric() function to coerce the tree to be ultrametric - see note above.")
    #    edge_details <- GetEdgeDetails(phy, includes.intervals=FALSE, intervening.intervals=NULL)
    #    if(any(edge_details$type == "extinct_tip")){
    #        phy <- force.ultrametric(phy)
    #    }
    #}
    
    if(sann == FALSE & is.null(starting.vals)){
        warning("You have chosen to rely on the internal starting points that generally work but does not guarantee finding the MLE.")
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
        pars.tmp <- c(pars.tmp, trans.tmp)
    }else{
        trans.tmp <- 0
    }
    pars.tmp <- c(turnover.tmp[1], eps.tmp[1], turnover.tmp[2], eps.tmp[2], turnover.tmp[3], eps.tmp[3], turnover.tmp[4], eps.tmp[4], turnover.tmp[5], eps.tmp[5], turnover.tmp[6], eps.tmp[6], turnover.tmp[7], eps.tmp[7], turnover.tmp[8], eps.tmp[8], turnover.tmp[9], eps.tmp[9], turnover.tmp[10], eps.tmp[10], turnover.tmp[11], eps.tmp[11], turnover.tmp[12], eps.tmp[12], turnover.tmp[13], eps.tmp[13], turnover.tmp[14], eps.tmp[14], turnover.tmp[15], eps.tmp[15], turnover.tmp[16], eps.tmp[16], turnover.tmp[17], eps.tmp[17], turnover.tmp[18], eps.tmp[18], turnover.tmp[19], eps.tmp[19], turnover.tmp[20], eps.tmp[20], turnover.tmp[21], eps.tmp[21], turnover.tmp[22], eps.tmp[22], turnover.tmp[23], eps.tmp[23], turnover.tmp[24], eps.tmp[24], turnover.tmp[25], eps.tmp[25], turnover.tmp[26], eps.tmp[26], trans.tmp[1])
    pars[1:length(pars.tmp)] <- pars.tmp
    
    if(includes.fossils == TRUE){
        pars <- c(pars, max(pars.tmp)+1)
    }else{
        pars <- c(pars, 0)
    }
    np <- max(pars)
    pars[pars==0] <- np + 1
    
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
            dat.tab <- OrganizeDataMiSSE(phy=phy, f=f, hidden.states=hidden.states, includes.intervals=FALSE, intervening.intervals=NULL, includes.fossils=includes.fossils)
            #These are all inputs for generating starting values:
            edge_details <- GetEdgeDetails(phy, includes.intervals=FALSE, intervening.intervals=NULL)
            fossil.taxa <- edge_details$tipward_node[which(edge_details$type == "extinct_tip")]
            fossil.ages <- dat.tab$TipwardAge[which(dat.tab$DesNode %in% fossil.taxa)]
            n.tax.starting <- Ntip(phy)-length(fossil.taxa)-no.k.samples
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
                dat.tab <- OrganizeDataMiSSE(phy=phy, f=f, hidden.states=hidden.states, includes.intervals=TRUE, intervening.intervals=strat.cache$intervening.intervals, includes.fossils=includes.fossils)
                #These are all inputs for generating starting values:
                edge_details <- GetEdgeDetails(phy, includes.intervals=TRUE, intervening.intervals=strat.cache$intervening.intervals)
                fossil.taxa <- edge_details$tipward_node[which(edge_details$type == "extinct_tip" | edge_details$type == "k_extinct_interval")]
                fossil.ages <- dat.tab$TipwardAge[which(dat.tab$DesNode %in% fossil.taxa)]
                n.tax.starting <- Ntip(phy)-length(fossil.taxa)-dim(fix.type)[1]
            }else{
                phy.og <- phy
                psi.type <- "m_only"
                fix.type <- NULL
                split.times <- dateNodes(phy, rootAge=max(node.depth.edgelength(phy)))[-c(1:Ntip(phy))]
                strat.cache <- NULL
                no.k.samples <- 0
                gen <- FindGenerations(phy)
                dat.tab <- OrganizeDataMiSSE(phy=phy, f=f, hidden.states=hidden.states, includes.intervals=FALSE, intervening.intervals=NULL, includes.fossils=includes.fossils)
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
        dat.tab <- OrganizeDataMiSSE(phy=phy, f=f, hidden.states=hidden.states, includes.intervals=FALSE, intervening.intervals=NULL, includes.fossils=includes.fossils)
        fossil.taxa <- NULL
        fix.type <- NULL
        psi.type <- NULL
        strat.cache <- NULL
    }
    nb.tip <- Ntip(phy)
    nb.node <- phy$Nnode
    ##########################
    
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
            if(includes.fossils == TRUE){
                stop("You input a tree that contains fossils but you are assuming no extinction. Check your function call.")
            }else{
                init.pars <- starting.point.generator(phy, k=2, samp.freq.tree, yule=TRUE)
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
            if(init.pars[3] == 0){
                init.pars[3] = 1e-6
            }
        }
        names(init.pars) <- NULL
        
        if(is.null(starting.vals)){
            def.set.pars <- rep(NA, 2*rate.cats)
            scaling=seq(from=1.2, to=0.8, length.out=rate.cats)
            for(cat.index in sequence(rate.cats)) {
                def.set.pars[1+(cat.index-1)*2] <- log((init.pars[1]+init.pars[3]) * ifelse(length(unique(turnover))==1, 1, scaling[cat.index]))
                def.set.pars[2+(cat.index-1)*2] <- log(ifelse(length(unique(eps))==1, 1, scaling[cat.index]) * init.pars[3]/init.pars[1])
            }
            #def.set.pars <- rep(c(log(init.pars[1]+init.pars[2]), log(init.pars[2]/init.pars[1])), rate.cats)
            #trans.start <- log(rate.cats/sum(phy$edge.length))
            trans.start <- log(init.pars[5])
        }else{
            if(includes.fossils == TRUE){
                ## Check the length of the stating vector:
                if( length( starting.vals ) != 4 ){
                    stop("Incorrect length for starting.vals vector.")
                }
                def.set.pars <- rep(c(log(starting.vals[1]), log(starting.vals[2])), rate.cats)
                trans.start <- log(starting.vals[3])
                psi.start <- log(starting.vals[4])
            }else{
                ## Check the length of the stating vector:
                if( length( starting.vals ) != 3 ){
                    stop("Incorrect length for starting.vals vector.")
                }
                def.set.pars <- rep(c(log(starting.vals[1]), log(starting.vals[2])), rate.cats)
                trans.start <- log(starting.vals[3])
            }
        }
        
        if(bounded.search == TRUE){
            upper.full <- rep(c(log(turnover.upper), log(eps.upper)), rate.cats)
        }else{
            upper.full <- rep(21, length(def.set.pars))
        }
        
        if(rate.cats > 1){
            if(includes.fossils == TRUE){
                np.sequence <- 1:(np-2)
                ip <- numeric(np-2)
                upper <- numeric(np-2)
            }else{
                np.sequence <- 1:(np-1)
                ip <- numeric(np-1)
                upper <- numeric(np-1)
            }
        }else{
            if(includes.fossils == TRUE){
                np.sequence <- 1:(np-1)
                ip <- numeric(np-1)
                upper <- numeric(np-1)
            }else{
                np.sequence <- 1:np
                ip <- numeric(np)
                upper <- numeric(np)
            }
        }
        
        for(i in np.sequence){
            ip[i] <- def.set.pars[which(pars == np.sequence[i])[1]]
            upper[i] <- upper.full[which(pars == np.sequence[i])[1]]
        }
        if(rate.cats > 1){
            if(includes.fossils == TRUE){
                ip <- c(ip, trans.start, log(psi.start))
                upper <- c(upper, log(trans.upper), log(trans.upper))
            }else{
                ip <- c(ip, trans.start)
                upper <- c(upper, log(trans.upper))
            }
        }else{
            if(includes.fossils == TRUE){
                ip <- c(ip, log(psi.start))
                upper <- c(upper, log(trans.upper))
            }
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
    if(sann == FALSE){
        if(bounded.search == TRUE){
            cat("Finished. Beginning bounded subplex routine...", "\n")
            opts <- list("algorithm" = "NLOPT_LN_SBPLX", "maxeval" = 100000, "ftol_rel" = max.tol)
            out = nloptr(x0=ip, eval_f=DevOptimizeMiSSE, ub=upper, lb=lower, opts=opts, pars=pars, dat.tab=dat.tab, gen=gen, hidden.states=hidden.states, fixed.eps=fixed.eps, nb.tip=nb.tip, nb.node=nb.node, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, np=np, ode.eps=ode.eps, fossil.taxa=fossil.taxa, fix.type=fix.type, strat.cache=strat.cache)
            sann.counts <- NULL
            solution <- numeric(length(pars))
            solution[] <- c(exp(out$solution), 0)[pars]
            loglik = -out$objective
        }else{
            cat("Finished. Beginning subplex routine...", "\n")
            out = subplex(ip, fn=DevOptimizeMiSSE, control=list(reltol=max.tol, parscale=rep(0.1, length(ip))), pars=pars, dat.tab=dat.tab, gen=gen, hidden.states=hidden.states, fixed.eps=fixed.eps, nb.tip=nb.tip, nb.node=nb.node, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, np=np, ode.eps=ode.eps, fossil.taxa=fossil.taxa, fix.type=fix.type, strat.cache=strat.cache)
            sann.counts <- NULL
            solution <- numeric(length(pars))
            solution[] <- c(exp(out$par), 0)[pars]
            loglik = -out$value
        }
    }else{
        cat("Finished. Beginning simulated annealing...", "\n")
        out.sann = GenSA(ip, fn=DevOptimizeMiSSE, lower=log(exp(lower)+0.0000000001), upper=log(exp(upper)-0.0000000001), control=list(max.call=sann.its, temperature=sann.temp, seed=sann.seed), pars=pars, dat.tab=dat.tab, gen=gen, hidden.states=hidden.states, fixed.eps=fixed.eps, nb.tip=nb.tip, nb.node=nb.node, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, np=np, ode.eps=ode.eps, fossil.taxa=fossil.taxa, fix.type=fix.type, strat.cache=strat.cache)
        sann.counts <- out.sann$counts
        cat("Finished. Refining using subplex routine...", "\n")
        opts <- list("algorithm" = "NLOPT_LN_SBPLX", "maxeval" = 100000, "ftol_rel" = max.tol)
        out <- nloptr(x0=out.sann$par, eval_f=DevOptimizeMiSSE, ub=upper, lb=lower, opts=opts, pars=pars, dat.tab=dat.tab, gen=gen, hidden.states=hidden.states, fixed.eps=fixed.eps, nb.tip=nb.tip, nb.node=nb.node, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, np=np, ode.eps=ode.eps, fossil.taxa=fossil.taxa, fix.type=fix.type, strat.cache=strat.cache)
        solution <- numeric(length(pars))
        solution[] <- c(exp(out$solution), 0)[pars]
        
        loglik = -out$objective
    }
    
    names(solution) <- c("turnover0A","eps0A", "turnover0B","eps0B", "turnover0C","eps0C", "turnover0D","eps0D", "turnover0E","eps0E", "turnover0F","eps0F", "turnover0G","eps0G", "turnover0H","eps0H", "turnover0I","eps0I", "turnover0J","eps0J", "turnover0K","eps0K", "turnover0L","eps0L", "turnover0M","eps0M", "turnover0N","eps0N", "turnover0O","eps0O", "turnover0P","eps0P", "turnover0Q","eps0Q", "turnover0R","eps0R", "turnover0S","eps0S", "turnover0T","eps0T", "turnover0U","eps0U", "turnover0V","eps0V","turnover0W","eps0W","turnover0X","eps0X", "turnover0Y","eps0Y", "turnover0Z","eps0Z", "q0", "psi")
    
    cat("Finished. Summarizing results...", "\n")
    misse_end_time <- Sys.time()
    obj = list(loglik = loglik, AIC = -2*loglik+2*np, AICc = -2*loglik+(2*np*(Ntip(phy)/(Ntip(phy)-np-1))), solution=solution, index.par=pars, f=f, hidden.states=hidden.states, fixed.eps=fixed.eps, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, phy=phy.og, phy.w.k=phy, max.tol=max.tol, starting.vals=ip, upper.bounds=upper, lower.bounds=lower, ode.eps=ode.eps, turnover=turnover, eps=eps, elapsed.minutes=difftime(misse_end_time, misse_start_time, units="min"), includes.fossils=includes.fossils, k.samples=k.samples, strat.intervals=strat.intervals, fix.type=fix.type, psi.type=psi.type, sann.counts=sann.counts)
    class(obj) <- append(class(obj), "misse.fit")
    return(obj)
}



######################################################################################################################################
######################################################################################################################################
### MiSSEGreedy -- Automated algorithm for running MiSSE of varying complexity
######################################################################################################################################
######################################################################################################################################

# options(error = utils::recover)
# a <- hisse:::MiSSEGreedyNew(ape::rcoal(50), possible.combos=hisse:::generateMiSSEGreedyCombinations(4), n.cores=4, save.file='~/Downloads/greedy.rda')
MiSSEGreedy <- function(phy, f=1, possible.combos = generateMiSSEGreedyCombinations(shuffle.start=TRUE), stop.deltaAICc=10, save.file=NULL, n.cores=NULL, chunk.size=10, check.fits=FALSE, remove.bad=FALSE, n.tries=2, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, includes.fossils=FALSE, k.samples=NULL, strat.intervals=NULL, sann=TRUE, sann.its=5000, sann.temp=5230, sann.seed=-100377, bounded.search=TRUE, max.tol=.Machine$double.eps^.50, starting.vals=NULL, turnover.upper=10000, eps.upper=3, trans.upper=100, restart.obj=NULL, ode.eps=0) {
    
    misse.list <- list()
    chunk.size <- ifelse(is.null(chunk.size),ifelse(is.null(n.cores),1,n.cores), chunk.size)
    total.chunks <- ceiling(nrow(possible.combos)/chunk.size)
    possible.combos$runorder <- NA
    possible.combos$deltaAICc <- NA
    possible.combos$AICc <- NA
	possible.combos$predictedAICc <- NA
	possible.combos$lnL <- NA
    possible.combos$AIC <- NA
    possible.combos$elapsedMinutes <- NA
    possible.combos$predictedMinutes <- NA
	final.combos <- possible.combos
	all.examined <- c()
    
    for (batch_index in sequence(total.chunks)) { # So, if we can do parallel, we do it in chunks so all cores are busy
        starting.time <- Sys.time()
		focal.models <- sequence(min(chunk.size, nrow(possible.combos)))
		if(batch_index>1) {
            #save(final.combos, possible.combos, file="stufffordebugging.Rsave")
            #final.combos$predictedAICc <- stats::predict(stats::glm(AICc ~ turnover + eps + turnover*eps, data=subset(final.combos, !is.na(final.combos$AICc))), newdata=final.combos) #Idea here is to focus on the best candidate models
			
			model.distances <- as.matrix(dist(final.combos[,c("turnover", "eps")], diag=TRUE, upper=TRUE))
			unknown.indices <- which(is.na(final.combos$AICc))
			for(unknown.index in unknown.indices) {
				model.distances[,unknown.index] <- 1e10	
			}
			for (model_index in sequence(nrow(final.combos))) {
				best <- which(model.distances[model_index,] == min(model.distances[model_index,]))
				final.combos$predictedAICc[model_index] <- min(final.combos[best, "AICc"]+model.distances[model_index, best], na.rm=TRUE)		
			}
			
            best.ones <- base::order(final.combos$predictedAICc)
			best.ones <- best.ones[!(best.ones %in% which(!is.na(final.combos$AICc)))]
			focal.models <- best.ones[1:min(chunk.size, length(best.ones))]
		}
        all.examined <- append(all.examined, focal.models)
		
        #local.combos <- possible.combos[(1+chunk.size*(batch_index-1)):min(nrow(possible.combos), chunk.size*batch_index) ,]
		local.combos <- final.combos[focal.models,]
        
        #cat("\nNow starting run with", paste(range(local.combos$turnover), collapse="-"), "turnover categories and", paste(range(local.combos$eps), collapse="-"), "extinction fraction categories", "\n")
        cat("Starting at ", as.character(starting.time), "\n running on ", n.cores, " cores.", sep="")
        cat("\n")
        print(local.combos[,1:3])
        cat("\n")
        
        misse.list <- append(misse.list, parallel::mcmapply(
        MiSSE,
        eps=local.combos$eps,
        turnover=local.combos$turnover,
        fixed.eps=local.combos$fixed.eps,
        MoreArgs=list(
        phy=phy,
        f=f,
        condition.on.survival=condition.on.survival,
        root.type=root.type,
        root.p=root.p,
        includes.fossils=includes.fossils,
        k.samples=k.samples,
        strat.intervals=strat.intervals,
        sann=sann,
        sann.its=sann.its,
        sann.temp=sann.temp,
        sann.seed=sann.seed,
        bounded.search=bounded.search,
        max.tol=max.tol,
        starting.vals=starting.vals,
        turnover.upper=turnover.upper,
        eps.upper=eps.upper,
        trans.upper=trans.upper,
        restart.obj=restart.obj,
        ode.eps=ode.eps,
        expand.mode=TRUE
        ),
        mc.cores=ifelse(is.null(n.cores),1,n.cores),
        SIMPLIFY=FALSE
        ))
        
        AICc <- unlist(lapply(misse.list, "[[", "AICc"))
        deltaAICc <- AICc-min(AICc)
        min.deltaAICc.this.chunk <- min(deltaAICc[(1+chunk.size*(batch_index-1)):min(nrow(final.combos), chunk.size*batch_index)])
        
        # final.combos$lnL[1:min(nrow(final.combos), chunk.size*batch_index)] <- unlist(lapply(misse.list, "[[", "loglik"))
        # final.combos$AIC[1:min(nrow(final.combos), chunk.size*batch_index)] <- unlist(lapply(misse.list, "[[", "AIC"))
        # final.combos$AICc[1:min(nrow(final.combos), chunk.size*batch_index)] <- AICc
        # final.combos$deltaAICc[1:min(nrow(final.combos), chunk.size*batch_index)] <- deltaAICc
        # final.combos$elapsedMinutes[1:min(nrow(final.combos), chunk.size*batch_index)] <- unlist(lapply(misse.list, "[[", "elapsed.minutes"))
        
        # data.for.fit <- data.frame(nparam=(final.combos$eps+possible.combos$turnover)[1:min(nrow(final.combos), chunk.size*batch_index)], logmin=log(final.combos$elapsedMinutes[1:min(nrow(final.combos), chunk.size*batch_index)]))
        # data.for.prediction <- data.frame(nparam=(final.combos$eps+final.combos$turnover))

        final.combos$lnL[all.examined] <- unlist(lapply(misse.list, "[[", "loglik"))
        final.combos$AIC[all.examined] <- unlist(lapply(misse.list, "[[", "AIC"))
        final.combos$AICc[all.examined] <- AICc
        final.combos$deltaAICc[all.examined] <- deltaAICc
        final.combos$elapsedMinutes[all.examined] <- unlist(lapply(misse.list, "[[", "elapsed.minutes"))
		final.combos$runorder[focal.models] <- batch_index
        
        data.for.fit <- data.frame(nparam=(final.combos$eps+final.combos$turnover)[all.examined], logmin=log(final.combos$elapsedMinutes[all.examined]))
        data.for.prediction <- data.frame(nparam=(final.combos$eps+final.combos$turnover))

        suppressWarnings(final.combos$predictedMinutes <- exp(predict(lm(logmin ~ nparam, data=data.for.fit), newdata=data.for.prediction)))
        
        cat("\nResults so far\n")
        final.combos[] <- lapply(final.combos, function(x) { if(!is.numeric(x)) as.numeric(x) else x })
        print(round(final.combos,2))
        
        if(!is.null(save.file)) {
            save(misse.list, final.combos, file=save.file)
        }
        
        if(batch_index<total.chunks) {
            if(stop.deltaAICc>min.deltaAICc.this.chunk) {
                print(paste0("Best AICc in this set of parallel runs was ", round(min.deltaAICc.this.chunk,2), " which is less than the cutoff to stop running (",stop.deltaAICc,"), so starting another set of parallel runs"))
            } else {
                print(paste0("Best AICc in this set of parallel runs was ", round(min.deltaAICc.this.chunk,2), " which is greater than the cutoff to stop running (",stop.deltaAICc,"), so stopping here"))
                break()
            }
        }
        
        # print("\n")
        # print(local.combos)
        # misse.list <- append(misse.list, MiSSE(
        #     eps=local.combos$eps[1],
        #     turnover=local.combos$turnover[1],
        #     fixed.eps=local.combos$fixed.eps[1],
        #
        #         phy=phy,
        #         f=f,
        #         condition.on.survival=condition.on.survival,
        #         root.type=root.type,
        #         root.p=root.p,
        #         sann=sann,
        #         sann.its=sann.its,
        #         bounded.search=bounded.search,
        #         max.tol=max.tol,
        #         starting.vals=starting.vals,
        #         turnover.upper=turnover.upper,
        #         eps.upper=eps.upper,
        #         trans.upper=trans.upper,
        #         restart.obj=restart.obj,
        #         ode.eps=ode.eps,
        #         expand.mode=TRUE
        # ))
    }

    if(check.fits == TRUE){
        cat("Checking model fits...", "\n")
        misse.list.updated <- MiSSENet(misse.list=misse.list, n.tries=n.tries, remove.bad=remove.bad, dont.rerun=FALSE, save.file=save.file, n.cores=ifelse(is.null(n.cores),1,n.cores), sann=sann, sann.its=sann.its, sann.temp=sann.temp, bounded.search=bounded.search, starting.vals=starting.vals, turnover.upper=turnover.upper, eps.upper=eps.upper, trans.upper=trans.upper, restart.obj=restart.obj)
        return(misse.list.updated)
    }else{
        if(remove.bad == TRUE){
            misse.list.updated <- MiSSENet(misse.list=misse.list, n.tries=n.tries, remove.bad=remove.bad, dont.rerun=TRUE, save.file=save.file, n.cores=ifelse(is.null(n.cores),1,n.cores), sann=sann, sann.its=sann.its, sann.temp=sann.temp, bounded.search=bounded.search, starting.vals=starting.vals, turnover.upper=turnover.upper, eps.upper=eps.upper, trans.upper=trans.upper, restart.obj=restart.obj)
            return(misse.list.updated)
        }else{
            return(misse.list)
        }
    }
}

SummarizeMiSSEGreedy <- function(greedy.result, min.weight=0.01, n.cores=1, recon=TRUE) {
	parameters <- c("turnover", "net.div", "speciation", "extinction", "extinction.fraction")
	AICc_weights <- GetAICWeights(greedy.result, criterion="AICc")
	good_enough <- which(AICc_weights>min.weight)
	best_model <- which.max(AICc_weights)
	rates <- NA
	if(recon) {
		recons <- list()
		rates <- array(dim=c(ape::Ntip(greedy.result[[1]]$phy) + ape::Nnode(greedy.result[[1]]$phy), length(parameters), length(good_enough)+2), dimnames=list(c(greedy.result[[1]]$phy$tip.label, (1+ape::Ntip(greedy.result[[1]]$phy)):(ape::Ntip(greedy.result[[1]]$phy) + ape::Nnode(greedy.result[[1]]$phy))), parameters, c(good_enough, "model_average", "best")))
		for(i in sequence(length(good_enough))){
			local_model <- greedy.result[[good_enough[i]]]
			
			recons[[i]] <- MarginReconMiSSE(phy=local_model$phy, f=local_model$f, pars=local_model$solution[local_model$index.par], hidden.states=local_model$hidden.states, fixed.eps=local_model$fixed.eps, condition.on.survival=local_model$condition.on.survival, root.type=local_model$root.type, root.p=local_model$root.p, includes.fossils=local_model$includes.fossils, k.samples=local_model$k.samples, strat.intervals=local_model$strat.intervals, AIC=local_model$AICc, get.tips.only=FALSE, verbose=FALSE, n.cores=n.cores)
		}

		for(param_index in sequence(length(parameters))) {
			rate.param <- parameters[param_index]
			for(i in sequence(length(good_enough))){
				rates.tips <- ConvertManyToRate(recons[i], rate.param, "tip.mat")
				rates.internal <- ConvertManyToRate(recons[i], rate.param, "node.mat")
				rates[,param_index,i] <- c(rates.tips, rates.internal)
				if(good_enough[i]==best_model) {
					rates[,param_index,2+length(good_enough)] <- c(rates.tips, rates.internal)
				}
				
			}
			rates.tips.avg <- ConvertManyToRate(recons, rate.param, "tip.mat")
			rates.internal.avg <- ConvertManyToRate(recons, rate.param, "node.mat")
			rates[,param_index,1+length(good_enough)] <- c(rates.tips.avg, rates.internal.avg)
		}
	}
	overview <- data.frame(model_number=sequence(length(greedy.result)))
	overview$lnL <- unlist(lapply(greedy.result, "[[", "loglik"))
	overview$AICc <- unlist(lapply(greedy.result, "[[", "AICc"))
	overview$deltaAICc <- overview$AICc-min(overview$AICc)
	overview$AIC <- unlist(lapply(greedy.result, "[[", "AIC"))
	overview$hidden.states <- unlist(lapply(greedy.result, "[[", "hidden.states"))
	overview$n.turnover <- NA
	overview$n.eps <- NA
	overview$n.transition_rate <- 1
	overview$n.freeparam <- NA
	for(i in sequence(length(greedy.result))) {
		overview$n.turnover[i] <- max(greedy.result[[i]]$turnover)
		overview$n.eps[i] <- max(greedy.result[[i]]$eps)
		overview$n.freeparam[i] <- overview$n.turnover[i]+overview$n.eps[i]+1
	}
	overview$AICc_weights <- AICc_weights
	overview$Included_In_Model_Average <- FALSE
	overview$Included_In_Model_Average[good_enough] <- TRUE
	overview$elapsed.minutes <- unlist(lapply(greedy.result, "[[", "elapsed.minutes"))

	return(list(overview=overview, rates=rates))
}



generateMiSSEGreedyCombinations <- function(max.param=52, turnover.tries=sequence(26), eps.tries=sequence(26), fixed.eps.tries=NA, vary.both=TRUE, shuffle.start=TRUE) {
    fixed.eps <- eps <- "CRAN wants this declared somehow"
    if(vary.both) {
        combos <- expand.grid(turnover=turnover.tries, eps=eps.tries, fixed.eps=fixed.eps.tries)
    } else {
        combos <- rbind(
        expand.grid(turnover=turnover.tries, eps=1, fixed.eps=fixed.eps.tries),
        expand.grid(turnover=1, eps=eps.tries, fixed.eps=fixed.eps.tries)
        )
    }
    combos$eps[which(!is.na(combos$fixed.eps))] <- 0
    rownames(combos) <- NULL
    combos <- subset(combos, eps==0 | is.na(fixed.eps)) # Don't estimate multiple eps while also fixing eps
    combos <- combos[!duplicated(combos),]
    combos <- combos[which(combos$turnover + combos$eps <= max.param),]
    combos <- combos[base::order(combos$eps, combos$turnover, decreasing=FALSE),]
    if(shuffle.start) {
        model_pref <- 1/(combos$turnover + 2*combos$eps)^4 #reorder with weight on simplest ones first
        combos <- combos[sample.int(n=nrow(combos), prob=model_pref, replace=FALSE), ]
    }
    rownames(combos) <- NULL
    return(combos)
}


#Original version -- now defunct but kept for posterity.
MiSSEGreedyOLD <- function(phy, f=1, turnover.tries=sequence(26), eps.constant=c(TRUE,FALSE), stop.count=2, stop.deltaAICc=10, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, sann=FALSE, sann.its=10000, bounded.search=TRUE, max.tol=.Machine$double.eps^.50, starting.vals=NULL, turnover.upper=10000, eps.upper=3, trans.upper=100, restart.obj=NULL, ode.eps=0, n.cores=NULL) {
    misse.list <- list()
    first.AICc <- Inf
    for (eps.index in seq_along(eps.constant)) {
        best.AICc <- first.AICc #reset so we start over at one turnover parameter and work our way back up
        times.since.close.enough <- 0
        for (turnover.index in seq_along(turnover.tries)) {
            if(turnover.index>1 | eps.index==1) { # don't do the 1 turnover, 1 eps model twice -- overcounts it
                starting.time <- Sys.time()
                turnover <- sequence(turnover.tries[turnover.index])
                eps <- turnover
                if(eps.constant[eps.index]) {
                    eps <- rep(1, length(turnover))
                }
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

DevOptimizeMiSSE <- function(p, pars, dat.tab, gen, hidden.states, fixed.eps, nb.tip, nb.node, condition.on.survival, root.type, root.p, np, ode.eps, fossil.taxa, fix.type, strat.cache) {
    #Generates the final vector with the appropriate parameter estimates in the right place:
    p.new <- exp(p)
    model.vec <- numeric(length(pars))
    model.vec[] <- c(p.new, 0)[pars]
    cache <- ParametersToPassMiSSE(model.vec=model.vec, hidden.states=hidden.states, fixed.eps=fixed.eps, nb.tip=nb.tip, nb.node=nb.node, bad.likelihood=exp(-300), ode.eps=ode.eps)
    if(!is.null(fix.type)){
        if(!is.null(strat.cache)){
            logl <- DownPassMisse(dat.tab=dat.tab, gen=gen, cache=cache, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, node=fix.type[,1], state=NULL, fossil.taxa=fossil.taxa, fix.type=fix.type[,2]) + (strat.cache$k*log(cache$psi)) + (cache$psi*strat.cache$l_s)
        }else{
            logl <- DownPassMisse(dat.tab=dat.tab, gen=gen, cache=cache, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, node=fix.type[,1], state=NULL, fossil.taxa=fossil.taxa, fix.type=fix.type[,2])
        }
    }else{
        logl <- DownPassMisse(dat.tab=dat.tab, gen=gen, cache=cache, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, node=NULL, fossil.taxa=fossil.taxa, fix.type=NULL)
    }
    return(-logl)
}



######################################################################################################################################
######################################################################################################################################
### The various utility functions used
######################################################################################################################################
######################################################################################################################################

GetTreeTable <- function(phy, root.age=NULL){
    if(is.null(root.age)){
        node.ages <- dateNodes(phy, rootAge=max(node.depth.edgelength(phy)))
    }else{
        node.ages <- dateNodes(phy, rootAge=root.age)
    }
    max.age <- max(node.ages)
    table.ages <- matrix(0, dim(phy$edge)[1], 2)
    for(row.index in 1:dim(phy$edge)[1]){
        table.ages[row.index,1] <- node.ages[phy$edge[row.index,1]]
        table.ages[row.index,2] <- node.ages[phy$edge[row.index,2]]
    }
    full.table <- cbind(table.ages, phy$edge.length, phy$edge)
    return(full.table)
}


GetEdgeDetails <- function(phy, includes.intervals=FALSE, intervening.intervals=NULL){
    table.info <- GetTreeTable(phy, root.age=NULL)
    edges_detail <- as.data.frame(table.info)
    colnames(edges_detail) <- c("rootward_age", "tipward_age", "edge.length", "rootward_node", "tipward_node")
    edges_detail$type <- "internal"
    edges_detail$type[which(edges_detail$tipward_node<=ape::Ntip(phy))] <- "extant_tip"
    for(row.index in 1:dim(table.info)[1]){
        if(table.info[row.index,5]<=ape::Ntip(phy)){
            if(table.info[row.index,2] > .Machine$double.eps^.50){
                edges_detail$type[row.index] <- "extinct_tip"
            }
        }
    }
    edges_detail$type[which(edges_detail$tipward_node %in% which(grepl("Ksamp_",phy$tip.label)))] <- "k_tip"

    if(includes.intervals == TRUE){
        Ksamples_rootward_nodes <- edges_detail$rootward_node[which(edges_detail$type=="k_tip")]
        Extinct_tipward_nodes <- edges_detail$tipward_node[which(edges_detail$type=="extinct_tip")]
        Extant_tipward_nodes <- edges_detail$tipward_node[which(edges_detail$type=="extant_tip")]
        tipwards_nodes_of_rootwards_nodes <- edges_detail$tipward_node[edges_detail$rootward_node %in% Ksamples_rootward_nodes]
        for(rootward_focal in Ksamples_rootward_nodes) {
            focal_rows <- which(edges_detail$rootward_node == rootward_focal)
            for(row_index in focal_rows) {
                ###For k --> k
                if(edges_detail$tipward_node[row_index] %in% Ksamples_rootward_nodes) {
                    edges_detail$type[row_index] <- "k_k_interval"
                }
                ###For k --> extinct
                if(edges_detail$tipward_node[row_index] %in% Extinct_tipward_nodes) {
                    edges_detail$type[row_index] <- "k_extinct_interval"
                }
                ###For k --> extant
                if(edges_detail$tipward_node[row_index] %in% Extant_tipward_nodes) {
                    edges_detail$type[row_index] <- "k_extant_interval"
                }
            }
        }
        
        #First go and find the intervening intervals. Should be k_k_intervals, so use intervening interval table to find and replace:
        if(!is.null(intervening.intervals)){
            tmp <- which(round(edges_detail$tipward_age,10) %in% round(intervening.intervals$o_i,10))
            for(match.index in 1:length(tmp)){
                if(edges_detail[tmp[match.index],]$type == "k_k_interval"){
                    edges_detail[tmp[match.index],]$type <- "intervening_interval"
                }
            }
        }
        
        #Next, check any tip intervals that are actually subtended by a k_k interval:
        for(check.index in 1:dim(edges_detail)[1]){
            if(edges_detail$type[check.index] == "k_extant_interval"){
                sister.taxa <- GetSister(phy, edges_detail$rootward_node[check.index])
                if(edges_detail$type[which(edges_detail$tipward_node==sister.taxa)] == "k_tip"){
                    rootward.of.sister.taxa <- edges_detail$rootward_node[which(edges_detail$tipward_node==sister.taxa)]
                    tmp <- edges_detail[which(edges_detail$rootward_node==rootward.of.sister.taxa),]
                    if(any(tmp$type == "k_k_interval")){
                        edges_detail$type[check.index] <- "extant_tip"
                    }
                }
            }
        }

        #Finally, check any extinct tip intervals that are actually subtended by a k_k interval:
        for(check.index in 1:dim(edges_detail)[1]){
            if(edges_detail$type[check.index] == "k_extinct_interval"){
                sister.taxa <- GetSister(phy, edges_detail$rootward_node[check.index])
                if(edges_detail$type[which(edges_detail$tipward_node==sister.taxa)] == "k_tip"){
                    rootward.of.sister.taxa <- edges_detail$rootward_node[which(edges_detail$tipward_node==sister.taxa)]
                    tmp <- edges_detail[which(edges_detail$rootward_node==rootward.of.sister.taxa),]
                    if(any(tmp$type == "k_k_interval")){
                        edges_detail$type[check.index] <- "extinct_tip"
                    }
                }
            }
        }
    }
    return(edges_detail)
}


OrganizeDataMiSSE <- function(phy, f, hidden.states, includes.intervals=FALSE, intervening.intervals=NULL, includes.fossils=FALSE){
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
    
    #This seems stupid but I cannot figure out how to get data.table to not make this column a factor. When a factor this is not right. For posterity, let it be known Jeremy would rather just retain the character instead of this mess:
    edge_details <- GetEdgeDetails(phy, includes.intervals=includes.intervals, intervening.intervals=intervening.intervals)
    if(includes.fossils == TRUE){
        branch.type <- edge_details$type
        branch.type[which(branch.type == "extant_tip")] <- 0
        branch.type[which(branch.type == "internal")] <- 0
        branch.type[which(branch.type == "extinct_tip")] <- 1
        branch.type[which(branch.type == "k_tip")] <- 2
        branch.type[which(branch.type == "k_k_interval")] <- 3
        branch.type[which(branch.type == "k_extinct_interval")] <- 3
        branch.type[which(branch.type == "k_extant_interval")] <- 3
        branch.type[which(branch.type == "intervening_interval")] <- 4
    }else{
        branch.type <- rep(0, length(edge_details$type))
    }
    
    tmp.df <- cbind(edge_details[,1:5], 0, matrix(0, nrow(edge_details), ncol(compD)), matrix(0, nrow(edge_details), ncol(compE)), as.numeric(branch.type))
    colnames(tmp.df) <- c("RootwardAge", "TipwardAge", "BranchLength", "FocalNode", "DesNode", "comp", paste("compD", 1:ncol(compD), sep="_"), paste("compE", 1:ncol(compE), sep="_"), "branch.type")
    dat.tab <- as.data.table(tmp.df)
    setkey(dat.tab, DesNode)
    cols <- names(dat.tab)
    for (j in 1:(dim(compD)[2])){
        #dat.tab[data.table(c(1:nb.tip)), paste("compD", j, sep="_") := compD[,j]]
        set(dat.tab, 1:nb.tip, cols[6+j], compD[,j])
        #dat.tab[data.table(c(1:nb.tip)), paste("compE", j, sep="_") := compE[,j]]
        set(dat.tab, 1:nb.tip, cols[32+j], compE[,j])
    }
    return(dat.tab)
}


SingleChildProbMiSSE <- function(cache, pars, compD, compE, start.time, end.time, branch.type){
    
    if(any(!is.finite(c(compD, compE)))) { # something went awry at a previous step. Bail!
        prob.subtree.cal <- rep(0, 26*2)
        prob.subtree.cal[(1+length(prob.subtree.cal)/2):length(prob.subtree.cal)] <- cache$bad.likelihood
        return(prob.subtree.cal)
    }
    
    yini <- c(E0A = compE[1], E0B = compE[2], E0C = compE[3], E0D = compE[4], E0E = compE[5], E0F = compE[6], E0G = compE[7], E0H = compE[8], E0I = compE[9], E0J = compE[10], E0K = compE[11], E0L = compE[12], E0M = compE[13], E0N = compE[14], E0O = compE[15], E0P = compE[16], E0Q = compE[17], E0R = compE[18], E0S = compE[19], E0T = compE[20], E0U = compE[21], E0V = compE[22], E0W = compE[23], E0X = compE[24], E0Y = compE[25], E0Z = compE[26], D0A = compD[1], D0B = compD[2], D0C = compD[3], D0D = compD[4], D0E = compD[5], D0F = compD[6], D0G = compD[7], D0H = compD[8], D0I = compD[9], D0J = compD[10], D0K = compD[11], D0L = compD[12], D0M = compD[13], D0N = compD[14], D0O = compD[15], D0P = compD[16], D0Q = compD[17], D0R = compD[18], D0S = compD[19], D0T = compD[20], D0U = compD[21], D0V = compD[22], D0W = compD[23], D0X = compD[24], D0Y = compD[25], D0Z = compD[26])
    if(branch.type == 2){
        times=c(0, end.time)
    }else{
        times=c(start.time, end.time)
    }
    runSilent <- function(branch.type) {
        options(warn = -1)
        on.exit(options(warn = 0))
        if(branch.type == 3){
            capture.output(res <- lsoda(yini, times, func = "misse_strat_derivs", pars, initfunc="initmod_misse", dllname = "hisse", rtol=1e-8, atol=1e-8))
        }else{
            capture.output(res <- lsoda(yini, times, func = "misse_derivs", pars, initfunc="initmod_misse", dllname = "hisse", rtol=1e-8, atol=1e-8))
        }
        res
    }
    #prob.subtree.cal.full <- lsoda(yini, times, func = "misse_derivs", pars, initfunc="initmod_misse", dll = "misse-ext-derivs", rtol=1e-8, atol=1e-8)
    #prob.subtree.cal.full <- lsoda(yini, times, func = "misse_derivs", pars, initfunc="initmod_misse", dllname = "hisse", rtol=1e-8, atol=1e-8)
    prob.subtree.cal.full <- runSilent(branch.type=branch.type)
    
    ######## THIS CHECKS TO ENSURE THAT THE INTEGRATION WAS SUCCESSFUL ###########
    if(attributes(prob.subtree.cal.full)$istate[1] < 0){
        prob.subtree.cal <- prob.subtree.cal.full[-1,-1]
        prob.subtree.cal[27:52] <- cache$bad.likelihood
        return(prob.subtree.cal)
    }else{
        if(branch.type == 2){
            prob.subtree.cal <- prob.subtree.cal.full[-1,-1]
            #So, what this does is simply maintains the original input which here should just be a bunch of ones.
            prob.subtree.cal[27:52] <- compD
        }else{
            if(branch.type == 4){
                prob.subtree.cal.wevents <- prob.subtree.cal.full[-1,-1]
                prob.subtree.cal.full <- runSilent(branch.type=3)
                prob.subtree.cal.noevents <- prob.subtree.cal.full[-1,-1]
                unobs.spec.probs <- 1 - (prob.subtree.cal.noevents[27:52]/prob.subtree.cal.wevents[27:52])
                prob.subtree.cal.wevents[27:52] <- prob.subtree.cal.wevents[27:52] * unobs.spec.probs
                prob.subtree.cal.wevents[which(is.na(prob.subtree.cal.wevents))] <- 0
                prob.subtree.cal <- prob.subtree.cal.wevents
            }else{
                prob.subtree.cal <- prob.subtree.cal.full[-1,-1]
            }
        }
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


FocalNodeProbMiSSE <- function(cache, pars, lambdas, dat.tab, generations){
    ### Ughy McUgherson. This is a must in order to pass CRAN checks: http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
    DesNode = NULL
    FocalNode = NULL
    . = NULL
    
    gens <- data.table(c(generations))
    setkey(dat.tab, FocalNode)
    CurrentGenData <- dat.tab[gens]
    tmp <- t(apply(CurrentGenData, 1, function(z) SingleChildProbMiSSE(cache, pars, z[7:32], z[33:58],  z[2], z[1], z[59])))
    v.mat <- matrix(tmp[seq(1,nrow(tmp)-1,2),27:52] * tmp[seq(2,nrow(tmp),2),27:52], length(unique(CurrentGenData$FocalNode)), 26)
    v.mat <- v.mat * matrix(lambdas, length(unique(CurrentGenData$FocalNode)), 26, byrow=TRUE)
    phi.mat <- matrix(tmp[seq(1,nrow(tmp)-1,2),1:26], length(unique(CurrentGenData$FocalNode)), 26)
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
                        fixer = numeric(26)
                        fixer[cache$state] = 1
                        v.mat[which(generations == cache$node[fix.index]),] <- v.mat[which(generations == cache$node[fix.index]),] * fixer
                    }
                }
            }
        }
    }
    tmp.comp <- rowSums(v.mat)
    tmp.probs <- v.mat / tmp.comp
    #tmp.probs <- v.mat
    setkey(dat.tab, DesNode)
    rows <- dat.tab[.(generations), which=TRUE]
    cols <- names(dat.tab)
    for (j in 1:(dim(tmp.probs)[2])){
        #dat.tab[data.table(c(generations)), paste("compD", j, sep="_") := tmp.probs[,j]]
        set(dat.tab, rows, cols[6+j], tmp.probs[,j])
        #dat.tab[data.table(c(generations)), paste("compE", j, sep="_") := phi.mat[,j]]
        set(dat.tab, rows, cols[32+j], phi.mat[,j])
    }
    dat.tab[data.table(c(generations)), "comp" := tmp.comp]
    return(dat.tab)
}


#Have to calculate root prob separately because it is not a descendant in our table. Could add it, but I worry about the NA that is required.
GetRootProbMiSSE <- function(cache, pars, lambdas, dat.tab, generations){
    ### Ughy McUgherson. This is a must in order to pass CRAN checks: http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
    FocalNode = NULL
    gens <- data.table(c(generations))
    setkey(dat.tab, FocalNode)
    CurrentGenData <- dat.tab[gens]
    tmp <- t(apply(CurrentGenData, 1, function(z) SingleChildProbMiSSE(cache, pars, z[7:32], z[33:58],  z[2], z[1], z[59])))
    v.mat <- matrix(tmp[seq(1,nrow(tmp)-1,2),27:52] * tmp[seq(2,nrow(tmp),2),27:52], length(unique(CurrentGenData$FocalNode)), 26)
    v.mat <- v.mat * matrix(lambdas, length(unique(CurrentGenData$FocalNode)), 26, byrow=TRUE)
    phi.mat <- matrix(tmp[seq(1,nrow(tmp)-1,2),1:26], length(unique(CurrentGenData$FocalNode)), 26)
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
                        #basically we are using the node to fix the state along a branch, but we do not want to assume a true speciation event occurred here.
                        lambdas.check <- lambdas
                        lambdas.check[which(lambdas==0)] <- 1
                        v.mat[which(generations == cache$node[fix.index]),] <- v.mat[which(generations == cache$node[fix.index]),] / lambdas.check
                    }else{
                        fixer = numeric(26)
                        fixer[cache$state] = 1
                        v.mat[which(generations == cache$node[fix.index]),] <- v.mat[which(generations == cache$node[fix.index]),] * fixer
                    }
                }
            }
        }
    }
    
    tmp.comp <- rowSums(v.mat)
    tmp.probs <- v.mat / tmp.comp
    #tmp.probs <- v.mat

    return(cbind(tmp.comp, phi.mat, tmp.probs))
}


#Calculates the initial conditions for fossil taxa in the tree.
GetFossilInitialsMiSSE <- function(cache, pars, lambdas, dat.tab, fossil.taxa){
    ### Ughy McUgherson. This is a must in order to pass CRAN checks: http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
    DesNode = NULL
    . = NULL
    
    fossils <- data.table(c(fossil.taxa))
    setkey(dat.tab, DesNode)
    CurrentGenData <- dat.tab[fossils]
    tmp <- t(apply(CurrentGenData, 1, function(z) SingleChildProbMiSSE(cache, pars, z[7:32], z[33:58], 0, z[2], 1)))
    tmp.probs <- matrix(tmp[,1:26], length(fossil.taxa), 26) * as.matrix(CurrentGenData[,7:32]) * cache$psi
    phi.mat <- matrix(tmp[,1:26], length(fossil.taxa), 26)
    setkey(dat.tab, DesNode)
    rows <- dat.tab[.(fossils), which=TRUE]
    cols <- names(dat.tab)
    for (j in 1:(dim(tmp.probs)[2])){
        #dat.tab[data.table(c(generations)), paste("compD", j, sep="_") := tmp.probs[,j]]
        set(dat.tab, rows, cols[6+j], tmp.probs[,j])
        #dat.tab[data.table(c(generations)), paste("compE", j, sep="_") := phi.mat[,j]]
        set(dat.tab, rows, cols[32+j], phi.mat[,j])
    }
    return(dat.tab)
}


######################################################################################################################################
######################################################################################################################################
### The MiSSE type down pass that carries out the integration and returns the likelihood:
######################################################################################################################################
######################################################################################################################################

DownPassMisse <- function(dat.tab, gen, cache, condition.on.survival, root.type, root.p, get.phi=FALSE, node=NULL, state=NULL, fossil.taxa=NULL, fix.type=NULL) {
    
    ### Ughy McUgherson. This is a must in order to pass CRAN checks: http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
    DesNode = NULL
    compE = NULL
    . = NULL
    
    pars <- c(cache$lambda0A, cache$mu0A, cache$lambda0B, cache$mu0B, cache$lambda0C, cache$mu0C, cache$lambda0D, cache$mu0D, cache$lambda0E, cache$mu0E, cache$lambda0F, cache$mu0F, cache$lambda0G, cache$mu0G, cache$lambda0H, cache$mu0H, cache$lambda0I, cache$mu0I, cache$lambda0J, cache$mu0J, cache$lambda0K, cache$mu0K, cache$lambda0L, cache$mu0L, cache$lambda0M, cache$mu0M, cache$lambda0N, cache$mu0N, cache$lambda0O, cache$mu0O, cache$lambda0P, cache$mu0P, cache$lambda0Q, cache$mu0Q, cache$lambda0R, cache$mu0R, cache$lambda0S, cache$mu0S, cache$lambda0T, cache$mu0T, cache$lambda0U, cache$mu0U, cache$lambda0V, cache$mu0V, cache$lambda0W, cache$mu0W, cache$lambda0X, cache$mu0X, cache$lambda0Y, cache$mu0Y, cache$lambda0Z, cache$mu0Z, cache$q0, cache$hidden.states, cache$psi)
    
    lambda <- c(cache$lambda0A, cache$lambda0B, cache$lambda0C, cache$lambda0D, cache$lambda0E, cache$lambda0F, cache$lambda0G, cache$lambda0H, cache$lambda0I, cache$lambda0J, cache$lambda0K, cache$lambda0L, cache$lambda0M, cache$lambda0N, cache$lambda0O, cache$lambda0P, cache$lambda0Q, cache$lambda0R, cache$lambda0S, cache$lambda0T, cache$lambda0U, cache$lambda0V, cache$lambda0W, cache$lambda0X, cache$lambda0Y, cache$lambda0Z)
    
    if(!is.null(fossil.taxa)){
        #Gets initial conditions for fossil taxa:
        dat.tab.copy <- copy(dat.tab)
        dat.tab.copy <- GetFossilInitialsMiSSE(cache=cache, pars=pars, lambdas=lambda, dat.tab=dat.tab.copy, fossil.taxa=fossil.taxa)
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
                    res.tmp <- GetRootProbMiSSE(cache=cache, pars=pars, lambdas=lambda, dat.tab=dat.tab.copy, generations=gen[[i]])
                    cache$node <- NULL
                    cache$state <- NULL
                    cache$fix.type <- NULL
                }else{
                    res.tmp <- GetRootProbMiSSE(cache=cache, pars=pars, lambdas=lambda, dat.tab=dat.tab.copy, generations=gen[[i]])
                }
            }else{
                res.tmp <- GetRootProbMiSSE(cache=cache, pars=pars, lambdas=lambda, dat.tab=dat.tab.copy, generations=gen[[i]])
            }
            compD.root <- res.tmp[c(28:53)]
            compE.root <- res.tmp[c(2:27)]
            setkey(dat.tab.copy, DesNode)
            comp <- dat.tab.copy[["comp"]]
            comp <- c(comp[-TIPS], res.tmp[1])
        }else{
            if(!is.null(node)){
                if(any(node %in% gen[[i]])){
                    cache$node <- node
                    cache$state <- state
                    cache$fix.type <- fix.type
                    dat.tab.copy <- FocalNodeProbMiSSE(cache, pars=pars, lambdas=lambda, dat.tab=dat.tab.copy, gen[[i]])
                    cache$node <- NULL
                    cache$state <- NULL
                    cache$fix.type <- NULL
                }else{
                    dat.tab.copy <- FocalNodeProbMiSSE(cache, pars=pars, lambdas=lambda, dat.tab=dat.tab.copy, gen[[i]])
                }
            }else{
                dat.tab.copy <- FocalNodeProbMiSSE(cache, pars=pars, lambdas=lambda, dat.tab=dat.tab.copy, gen[[i]])
            }
        }
    }
        
    if(!is.null(fossil.taxa)){
        #this overrides the post order traversal E and recalculates assuming psi=0. See pg. 400 Stadler 2010.
        #if(all(fix.type == "event")){ ## not sure yet for intervals. This is unclear in Stadler et al 2018.
            pars[length(pars)] <- 0
            cache$psi <- 0
            phi.mat <- SingleChildProbMiSSE(cache, pars, c(as.matrix(dat.tab[1,7:32])), c(as.matrix(dat.tab[1, 33:58])),  0, max(dat.tab$RootwardAge), 0)
            compE.root <- matrix(phi.mat[1:26], 1, 26)
        #}
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
                #lambda <- c(cache$lambda0A, cache$lambda0B, cache$lambda0C, cache$lambda0D, cache$lambda0E, cache$lambda0F, cache$lambda0G, cache$lambda0H, cache$lambda0I, cache$lambda0J, cache$lambda0K, cache$lambda0L, cache$lambda0M, cache$lambda0N, cache$lambda0O, cache$lambda0P, cache$lambda0Q, cache$lambda0R, cache$lambda0S, cache$lambda0T, cache$lambda0U, cache$lambda0V, cache$lambda0W, cache$lambda0X, cache$lambda0Y, cache$lambda0Z)
                compD.root <- compD.root / sum(root.p * lambda * (1 - compE.root)^2)
                #Corrects for possibility that you have 0/0:
                compD.root[which(is.na(compD.root))] = 0
                loglik <- log(sum(compD.root * root.p)) + sum(log(comp))
                #loglik <- log(sum(compD.root * root.p))
            }else{
                #lambda <- c(cache$lambda0A, cache$lambda0B, cache$lambda0C, cache$lambda0D, cache$lambda0E, cache$lambda0F, cache$lambda0G, cache$lambda0H, cache$lambda0I, cache$lambda0J, cache$lambda0K, cache$lambda0L, cache$lambda0M, cache$lambda0N, cache$lambda0O, cache$lambda0P, cache$lambda0Q, cache$lambda0R, cache$lambda0S, cache$lambda0T, cache$lambda0U, cache$lambda0V, cache$lambda0W, cache$lambda0X, cache$lambda0Y, cache$lambda0Z)
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
        obj$compE = compE.root
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

ParametersToPassMiSSE <- function(model.vec, hidden.states, fixed.eps, nb.tip, nb.node, bad.likelihood, ode.eps){
    #Provides an initial object that contains all the parameters to be passed among functions. This will also be used to pass other things are we move down the tree (see DownPassGeoSSE):
    obj <- NULL
    
    obj$hidden.states <- hidden.states
    obj$nb.tip <- nb.tip
    obj$nb.node <- nb.node
    obj$bad.likelihood <- bad.likelihood
    obj$ode.eps <- ode.eps
    
    if(is.null(fixed.eps)){
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
    }else{
        ##Hidden State A
        obj$lambda0A = model.vec[1] / (1 + fixed.eps)
        obj$mu0A = (fixed.eps * model.vec[1]) / (1 + fixed.eps)
        
        ##Hidden State B
        obj$lambda0B = model.vec[3] / (1 + fixed.eps)
        obj$mu0B = (fixed.eps * model.vec[3]) / (1 + fixed.eps)
        
        ##Hidden State C
        obj$lambda0C = model.vec[5] / (1 + fixed.eps)
        obj$mu0C = (fixed.eps * model.vec[5]) / (1 + fixed.eps)
        
        ##Hidden State D
        obj$lambda0D = model.vec[7] / (1 + fixed.eps)
        obj$mu0D = (fixed.eps * model.vec[7]) / (1 + fixed.eps)
        
        ##Hidden State E
        obj$lambda0E = model.vec[9] / (1 + fixed.eps)
        obj$mu0E = (fixed.eps * model.vec[9]) / (1 + fixed.eps)
        
        ##Hidden State F
        obj$lambda0F = model.vec[11] / (1 + fixed.eps)
        obj$mu0F = (fixed.eps * model.vec[11]) / (1 + fixed.eps)
        
        ##Hidden State G
        obj$lambda0G = model.vec[13] / (1 + fixed.eps)
        obj$mu0G = (fixed.eps * model.vec[13]) / (1 + fixed.eps)
        
        ##Hidden State H
        obj$lambda0H = model.vec[15] / (1 + fixed.eps)
        obj$mu0H = (fixed.eps * model.vec[15]) / (1 + fixed.eps)
        
        ##Hidden State I
        obj$lambda0I = model.vec[17] / (1 + fixed.eps)
        obj$mu0I = (fixed.eps * model.vec[17]) / (1 + fixed.eps)
        
        ##Hidden State J
        obj$lambda0J = model.vec[19] / (1 + fixed.eps)
        obj$mu0J = (fixed.eps * model.vec[19]) / (1 + fixed.eps)
        
        ##Hidden State K
        obj$lambda0K = model.vec[21] / (1 + fixed.eps)
        obj$mu0K = (fixed.eps * model.vec[21]) / (1 + fixed.eps)
        
        ##Hidden State L
        obj$lambda0L = model.vec[23] / (1 + fixed.eps)
        obj$mu0L = (fixed.eps * model.vec[23]) / (1 + fixed.eps)
        
        ##Hidden State M
        obj$lambda0M = model.vec[25] / (1 + fixed.eps)
        obj$mu0M = (fixed.eps * model.vec[25]) / (1 + fixed.eps)
        
        ##Hidden State N
        obj$lambda0N = model.vec[27] / (1 + fixed.eps)
        obj$mu0N = (fixed.eps * model.vec[27]) / (1 + fixed.eps)
        
        ##Hidden State O
        obj$lambda0O = model.vec[29] / (1 + fixed.eps)
        obj$mu0O = (fixed.eps * model.vec[29]) / (1 + fixed.eps)
        
        ##Hidden State P
        obj$lambda0P = model.vec[31] / (1 + fixed.eps)
        obj$mu0P = (fixed.eps * model.vec[31]) / (1 + fixed.eps)
        
        ##Hidden State Q
        obj$lambda0Q = model.vec[33] / (1 + fixed.eps)
        obj$mu0Q = (fixed.eps * model.vec[33]) / (1 + fixed.eps)
        
        ##Hidden State R
        obj$lambda0R = model.vec[35] / (1 + fixed.eps)
        obj$mu0R = (fixed.eps * model.vec[35]) / (1 + fixed.eps)
        
        ##Hidden State S
        obj$lambda0S = model.vec[37] / (1 + fixed.eps)
        obj$mu0S = (fixed.eps * model.vec[37]) / (1 + fixed.eps)
        
        ##Hidden State T
        obj$lambda0T = model.vec[39] / (1 + fixed.eps)
        obj$mu0T = (fixed.eps * model.vec[39]) / (1 + fixed.eps)
        
        ##Hidden State U
        obj$lambda0U = model.vec[41] / (1 + fixed.eps)
        obj$mu0U = (fixed.eps * model.vec[41]) / (1 + fixed.eps)
        
        ##Hidden State V
        obj$lambda0V = model.vec[43] / (1 + fixed.eps)
        obj$mu0V = (fixed.eps * model.vec[43]) / (1 + fixed.eps)
        
        ##Hidden State W
        obj$lambda0W = model.vec[45] / (1 + fixed.eps)
        obj$mu0W = (fixed.eps * model.vec[45]) / (1 + fixed.eps)
        
        ##Hidden State X
        obj$lambda0X = model.vec[47] / (1 + fixed.eps)
        obj$mu0X = (fixed.eps * model.vec[47]) / (1 + fixed.eps)
        
        ##Hidden State Y
        obj$lambda0Y = model.vec[49] / (1 + fixed.eps)
        obj$mu0Y = (fixed.eps * model.vec[49]) / (1 + fixed.eps)
        
        ##Hidden State Z
        obj$lambda0Z = model.vec[51] / (1 + fixed.eps)
        obj$mu0Z = (fixed.eps * model.vec[51]) / (1 + fixed.eps)
    }
    obj$q0 = model.vec[53]
    obj$psi = model.vec[54]
    
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
    if(!is.null(x$fixed.eps)){
        cat("Fixed eps used: \n")
        cat("\n")
        print(x$fixed.eps)
        cat("\n")
    }
    if(!is.null(x$psi.type)){
        if(x$psi.type == "m_only"){
            cat("psi estimate reflects fossil tip sampling only. \n")
        }
        if(x$psi.type == "m+k"){
            cat("psi estimate reflects both fossil edge and tip sampling. \n")
        }
        if(x$psi.type == "m+int"){
            cat("psi estimate reflects stratigraphic interval sampling. \n")
        }
    }
}




#phy <- read.tree("../vignettes/whales_Steemanetal2009.tre")
#phy <- read.tree("whales_Slateretal2010.tre")
## print(p.new)
#gen <- hisse:::FindGenerations(phy)
#dat.tab <- hisse:::OrganizeDataMiSSE(phy=phy, f=1, hidden.states=1)
#nb.tip <- Ntip(phy)
#nb.node <- phy$Nnode
#model.vec <- c(0.103624, 5.207178e-09, rep(0,52), 1)

#cache = hisse:::ParametersToPassMiSSE(model.vec=model.vec, hidden.states=1, nb.tip=nb.tip, nb.node=nb.node, bad.likelihood=exp(-250), ode.eps=0)#
#logl <- hisse:::DownPassMisse(dat.tab=dat.tab, cache=cache, gen=gen, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL)
#right.logl <- -277.6942
#round(logl,4) == round(right.logl,4)

