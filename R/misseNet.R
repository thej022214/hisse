######################################################################################################################################
######################################################################################################################################
### Assessses model optimization and triages those that fail
######################################################################################################################################
######################################################################################################################################

#Takes a vector of free parameters and find the distance between them and gives back an
#edge matrix for models 1 freepar away from each other
GetEdges <- function(turn.free.p, eps.free.p, nodes){
    #Step 1 get model distance based on turnover
    distances.turn <- dist(turn.free.p)
    dist.mat.turn <- as.matrix(distances.turn)
    dist.mat.turn[lower.tri(dist.mat.turn)] <- max(dist.mat.turn)
    
    #Step 2: Get model distance based on eps
    distances.eps <- dist(eps.free.p)
    dist.mat.eps <- as.matrix(distances.eps)
    dist.mat.eps[lower.tri(dist.mat.eps)] <- max(dist.mat.eps)
    
    #Step 3: Get model diff of 1 in turn direction, zero for eps
    one.par.diff.turn <- which(dist.mat.turn == 1 & dist.mat.eps == 0, arr.ind=TRUE)
    #Step 4: Get model diff of 1 in eps direction, zero for turn
    one.par.diff.eps <- which(dist.mat.eps == 1 & dist.mat.turn == 0, arr.ind=TRUE)
    #Step 5: Combine Step 3 and 4
    one.par.diff <- rbind(one.par.diff.turn, one.par.diff.eps)
    #Step 5: Loop through table to assign proper node names for use in igraph
    tmp.edges.mat <- c()
    for(index in 1:dim(one.par.diff)[1]){
        tmp.edges.mat <- rbind(tmp.edges.mat, c(nodes[one.par.diff[index,1]], nodes[one.par.diff[index,2]]))
    }
    final.edge.mat <- data.frame(x=tmp.edges.mat[,1], y=tmp.edges.mat[,2], weight=rep(1,dim(tmp.edges.mat)[1]))
    #Fin.
    return(final.edge.mat)
}


FindBadOptim <- function(graph, nodes, loglik.vec){
    bad.models <- c()
    for(node.index in 1:length(nodes)){
        desc.nodes <- igraph::neighbors(graph, igraph::V(graph)[node.index])
        for(neighbor.index in 1:length(desc.nodes)){
            if(length(desc.nodes) > 0){
                if(round(loglik.vec[node.index],3) - round(loglik.vec[desc.nodes[neighbor.index]],3) > 0.1){
                    bad.models <- c(bad.models, desc.nodes[neighbor.index])
                }
            }
        }
    }
    if(length(bad.models) == 0){
        return(0)
    }else{
        return(unique(bad.models))
    }
}


RerunBadOptim <- function(bad.fits, misse.list, sann, sann.its, sann.temp, bounded.search, starting.vals, turnover.upper, eps.upper, trans.upper, restart.obj, retries){

    if(sann == FALSE){
        sann = TRUE
        bounded.search = TRUE
    }
    current.loglik <- misse.list[[bad.fits]]$loglik
    current.fit <- misse.list[[bad.fits]]
    print(paste("current loglik:", current.loglik))
    sann.its <- sann.its * 2
    sann.temp <- sann.temp
    sann.seed <- c(-487281, -391945, -149149, -919193, -682017)[retries[[bad.fits]]+1]
    redo.fit <- MiSSE(misse.list[[bad.fits]]$phy, f=misse.list[[bad.fits]]$f, turnover=misse.list[[bad.fits]]$turnover, eps=misse.list[[bad.fits]]$eps, condition.on.survival=misse.list[[bad.fits]]$condition.on.survival, root.type=misse.list[[bad.fits]]$root.type, root.p=misse.list[[bad.fits]]$root.p, includes.fossils=misse.list[[bad.fits]]$includes.fossils, k.samples=misse.list[[bad.fits]]$k.samples, sann=sann, sann.its=sann.its, sann.temp=sann.temp, sann.seed=sann.seed, bounded.search=bounded.search, max.tol=misse.list[[bad.fits]]$max.tol, starting.vals=starting.vals, turnover.upper=turnover.upper, eps.upper=eps.upper, trans.upper=trans.upper, restart.obj=misse.list[[bad.fits]], ode.eps=misse.list[[bad.fits]]$ode.eps)
    print(paste("new loglik:", redo.fit$loglik))
    if(redo.fit$loglik > current.loglik){
        misse.list[[bad.fits]] <- redo.fit
        return(redo.fit)
    }
    return(misse.list[[bad.fits]])
}


MiSSENet <- function(misse.list, n.tries=2, remove.bad=TRUE, dont.rerun=FALSE, save.file=NULL, n.cores=1, sann=TRUE, sann.its=5000, sann.temp=5230, bounded.search=TRUE, starting.vals=NULL, turnover.upper=10000, eps.upper=3, trans.upper=100, restart.obj=NULL){
    
    #Step 1: Make igraph of examined model space
    tmp <- c()
    for(model.index in 1:length(misse.list)){
        tmp <- rbind(tmp, c(length(unique(misse.list[[model.index]]$turnover)), length(unique(misse.list[[model.index]]$eps)), misse.list[[model.index]]$loglik, misse.list[[model.index]]$AICc))
    }
    model.space <- data.frame(turnover=tmp[,1], eps=tmp[,2], loglik=tmp[,3], aic=tmp[,4])
    nodes <- paste("T", model.space$turnover, "_", "E", model.space$eps, sep="")
    edges <- GetEdges(model.space$turnover, model.space$eps, nodes)
    graph.df <- igraph::graph.data.frame(edges, vertices=nodes)
    
    #Step 2: Locate the bad fits by traversing graph in step 1
    bad.fits <- FindBadOptim(graph=graph.df, nodes=nodes, loglik.vec=model.space$loglik)
    
    if(dont.rerun == TRUE){
        if(remove.bad == TRUE){
            for(bad.index in 1:length(bad.fits)){
                misse.list[[bad.fits[bad.index]]] <- NA
            }
            misse.list <- Filter(Negate(anyNA), misse.list)
            return(misse.list)
        }else{
            return(bad.fits)
        }
    }else{
        if(bad.fits[1] != 0){
            model.retries <- as.list(as.list(numeric(length(misse.list))))
            done.enough <- FALSE
            while(done.enough == FALSE){
                cat("Current set of model fits identified as poorly optimized:", bad.fits, "\n")
                #Step 3: Rerun the bad fits by increasing simulated annealing temperature? Still trying to figure out
                bad.fits <- bad.fits[which(unlist(model.retries)[bad.fits] < n.tries)]
                rerun.fits <- parallel::mclapply(bad.fits, RerunBadOptim, misse.list=misse.list, sann=sann, sann.its=sann.its, sann.temp=sann.temp, bounded.search=bounded.search, starting.vals=starting.vals, turnover.upper=turnover.upper, eps.upper=eps.upper, trans.upper=trans.upper, restart.obj=restart.obj, retries=model.retries, mc.cores=ifelse(is.null(n.cores),1,n.cores))
                ### Keeping this here for debugging purposes ###
                #RerunBadOptim(bad.fits=bad.fits[1], misse.list=misse.list, sann=sann, sann.its=sann.its, sann.temp=sann.temp, bounded.search=bounded.search, starting.vals=starting.vals, turnover.upper=turnover.upper, eps.upper=eps.upper, trans.upper=trans.upper, restart.obj=restart.obj)
                ################################################
                for(rerun.index in 1:length(rerun.fits)){
                    misse.list[[bad.fits[rerun.index]]] <- rerun.fits[[rerun.index]]
                    model.retries[[bad.fits[rerun.index]]] <- model.retries[[bad.fits[rerun.index]]] + 1
                }
                
                #Repeat step 1:
                tmp <- c()
                for(model.index in 1:length(misse.list)){
                    tmp <- rbind(tmp, c(length(unique(misse.list[[model.index]]$turnover)), length(unique(misse.list[[model.index]]$eps)), misse.list[[model.index]]$loglik, misse.list[[model.index]]$AICc))
                }
                model.space <- data.frame(turnover=tmp[,1], eps=tmp[,2], loglik=tmp[,3], aic=tmp[,4])
                nodes <- paste("T", model.space$turnover, "_", "E", model.space$eps, sep="")
                edges <- GetEdges(model.space$turnover, model.space$eps, nodes)
                graph.df <- igraph::graph.data.frame(edges, vertices=nodes)
                bad.fits <- FindBadOptim(graph=graph.df, nodes=nodes, loglik.vec=model.space$loglik)
                if(all(unlist(model.retries)[bad.fits] == n.tries) | bad.fits[1] == 0){
                    done.enough = TRUE
                }else{
                    cat("One or more models still identified as poorly optimized.", "\n")
                }
            }
            if(bad.fits[1] != 0){
                cat("Maximum number of reassessments reached. These models should be discarded:", bad.fits, "\n")
                if(remove.bad == TRUE){
                    if(!is.null(save.file)) {
                        save(misse.list, file=save.file)
                    }
                    if(length(bad.fits)){
                        for(bad.index in 1:length(bad.fits)){
                            misse.list[[bad.fits[bad.index]]] <- NA
                        }
                        misse.list <- Filter(Negate(anyNA), misse.list)
                    }
                }
            }
            return(misse.list)
        }else{
            cat("All models appeared to have optimized well.", "\n")
            return(misse.list)
        }
    }
}

