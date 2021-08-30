

######################################################################################################################################
######################################################################################################################################
### Simulator for birth-death models with fossils and sampled lineages
######################################################################################################################################
######################################################################################################################################

GetSister <- function(tree,node,mode=c("number","label")){
    mode <- mode[1]
    if(is.character(node)) node <- match(node,c(tree$tip.label,tree$node.label))
    sisters <- tree$edge[which(tree$edge[,1]==tree$edge[which(tree$edge[,2]==node),1]),2]
    sisters <- setdiff(sisters,node)
    if(mode=="number") return(sisters)
    else if(mode=="label"){
        result <- list()
        n <- length(tree$tip.label)
        if(is.null(tree$node.label)&&any(sisters>n)) result$nodes<-sisters[which(sisters>n)]
        else if(any(sisters>n)) result$nodes<-tree$node.label[sisters[which(sisters>n)]-n]
        if(any(sisters<=n)) result$tips<-tree$tip.label[sisters[which(sisters<=n)]]
        return(result)
    }
}


GetEdgeCombined <- function(phy) {
    edge_combined <- as.data.frame(cbind(phy$edge, phy$edge.length))
    extinct_taxa <- geiger::is.extinct(phy, tol=.Machine$double.eps^.50)
    extinct_taxa_indices <- which((phy$tip.label %in% extinct_taxa))
    living_taxa_indices <- which(!(phy$tip.label %in% extinct_taxa))
    edge_combined$livingdescendants <- FALSE
    edge_combined$livingdescendants[edge_combined[,2] %in% living_taxa_indices] <- TRUE
    for (i in seq_along(living_taxa_indices)) {
        ancestral_nodes <- phangorn::Ancestors(phy, node=living_taxa_indices[i])
        edge_combined$livingdescendants[edge_combined[,2] %in% ancestral_nodes] <- TRUE
    }
    if(length(living_taxa_indices)>0) {
        edge_combined$livingdescendants[edge_combined[,1]==max(edge_combined[,1])] <- TRUE
    }
    edge_combined <- cbind(edge_combined, phytools::nodeHeights(phy))
    colnames(edge_combined) <- c("rootwardnode", "tipwardnode", "edge.length",  "livingdescendants", "rootwardheight", "tipwardheight")
    edge_combined$istip <- FALSE
    edge_combined$istip[edge_combined[,2] <= ape::Ntip(phy)] <- TRUE
    return(edge_combined)
}


GetFossils <- function(phy, psi=0.1) {
    edge_combined=NULL
    if(is.null(edge_combined)) {
        edge_combined <- GetEdgeCombined(phy)
    }
    maxheight = max(edge_combined$tipwardheight)
    fossils <- as.data.frame(matrix(ncol=9, nrow=0))
    colnames(fossils) <- c("rootwardnode", "tipwardnode", "timefromroot", "timefromrootwardnode", "timefrompresent", "fossiltype_short", "fossiltype_long", "order_on_branch", "has_sampled_descendant")
    for(i in sequence(nrow(edge_combined))) {
        total_time <- edge_combined$edge.length[i]
        starting_time <- edge_combined$rootwardheight[i]
        elapsed_time <- rexp(1, rate=psi)
        order_on_branch <- 1
        while(elapsed_time < total_time) {
            fossiltype_long <- paste0(ifelse(edge_combined$livingdescendants[i], "surviving_", "extinct_"), ifelse(edge_combined$istip[i], "terminal", "internal"))
            fossils <- rbind(fossils, data.frame(rootwardnode=edge_combined$rootwardnode[i], tipwardnode=edge_combined$tipwardnode[i], timefromroot=elapsed_time+starting_time, timefromrootwardnode=elapsed_time, timefrompresent=maxheight-(elapsed_time+starting_time), fossiltype_short=ifelse(edge_combined$livingdescendants[i], "from_surviving_lineage", "from_extinct"), fossiltype_long=fossiltype_long, order_on_branch=order_on_branch))
            elapsed_time <- elapsed_time + rexp(1, rate=psi)
            order_on_branch <- order_on_branch + 1
        }
    }
    fossils$has_sampled_descendant <- FALSE
    fossils$has_sampled_descendant[fossils$fossiltype_short=="from_surviving_lineage"] <- TRUE # clearly true b/c it has an extant tip
    fossils <- fossils[order(fossils$tipwardnode, fossils$order_on_branch, decreasing=TRUE),]
    fossils$has_sampled_descendant[duplicated(fossils$tipwardnode)] <- TRUE  # there's a more recent fossil on the same branch

	# loop to find k fossils that have sampled descendants from branches that are extinct
	for (i in sequence(nrow(fossils))) {
		all_descendants <- phangorn::Descendants(phy, fossils$tipwardnode[i], type="all")
#		all_descendants <- all_descendants[which(all_descendants>ape::Ntip(phy))] # so we only look at the internal nodes
		other_fossils <- subset(fossils, fossils$tipwardnode!=fossils$tipwardnode[i])
		if(length(all_descendants)>0) {
			if(any(all_descendants %in% other_fossils$rootwardnode)) {
				fossils$has_sampled_descendant[i] <- TRUE
			}
            if(any(fossils$tipwardnode[i] %in% other_fossils$rootwardnode)) {
                fossils$has_sampled_descendant[i] <- TRUE
            }
		}
	}
	
    #fossils$fossiltype_mk <- NA
    #fossils$fossiltype_mk[which(fossils$fossiltype_long=="extinct_terminal")] <- "m"
    #fossils$fossiltype_mk[which(fossils$fossiltype_long=="extinct_internal" & !fossils$has_sampled_descendant)] <- "m"
    #fossils$fossiltype_mk[which(fossils$fossiltype_long=="extinct_internal" & fossils$has_sampled_descendant)] <- "k"
    
    fossils <- fossils[order(fossils$timefrompresent),]
    
    return(fossils)
}


GetStratigraphicIntervals <- function(phy, f) {
    f$edge_root_tip <- paste0(f$rootwardnode, "_", f$tipwardnode)
    strat_intervals <- data.frame(edge_root_tip=unique(f$edge_root_tip), rootwardnode=NA, tipwardnode=NA, startingtimefromroot=NA, endingtimefromroot=NA, intervallength=0, tipwardendisextant=FALSE)
    for (i in sequence(nrow(strat_intervals))) {
        matching_f <- subset(f,edge_root_tip== strat_intervals$edge_root_tip[i])
        strat_intervals$rootwardnode[i] <- matching_f$rootwardnode[1]
        strat_intervals$tipwardnode[i] <- matching_f$tipwardnode[1]
        if(nrow(matching_f)==1 && matching_f$fossiltype_long[1]=="surviving_terminal") {
            strat_intervals$tipwardendisextant[i] <- TRUE
        }
        strat_intervals$startingtimefromroot[i] <- min(matching_f$timefromroot)
        strat_intervals$endingtimefromroot[i] <- ifelse(strat_intervals$tipwardendisextant[i], max(node.depth.edgelength(phy)), max(matching_f$timefromroot))
        strat_intervals$intervallength[i] <- strat_intervals$endingtimefromroot[i] - strat_intervals$startingtimefromroot[i]
    }
    return(strat_intervals)
}


ProcessSimSample <- function(phy, f, strat.intervals=FALSE){
    
    #Step 1: Get MRCA of K samples for table. We want MRCA1, MRCA2, TIMEFROMPRESENT, Trait
    k.samples <- f[which(f$has_sampled_descendant == TRUE),]
    tmp.k.samples <- as.list(1:length(k.samples$tipwardnode))
    for(sample.index in sequence(length(k.samples$tipwardnode))){
        all.desc <- phytools::getDescendants(phy, node=k.samples$tipwardnode[sample.index])
        tips.only <- all.desc[all.desc<(Ntip(phy)+1)]
        tmp.k.samples[[sample.index]] <- phy$tip.label[tips.only]
    }

    #Step 2: Place sampled extinct taxa in the tree.
    extinct.samples <- f[which(f$fossiltype_long=="extinct_terminal" | f$fossiltype_long=="extinct_internal"),]
    extinct.samples <- extinct.samples[which(extinct.samples$has_sampled_descendant == FALSE),]
    
    ntip <- Ntip(phy)
    tip.desc.ext.list <- as.list(1:length(extinct.samples$tipwardnode))
    for(og.node.index in sequence(length(extinct.samples$tipwardnode))){
        all.desc <- phytools::getDescendants(phy, node=extinct.samples[og.node.index,2])
        tip.desc.ext.list[[og.node.index]] <- all.desc[all.desc<(Ntip(phy)+1)]
    }
    
    new.fossil.names <- c()
    for(sample.index in sequence(length(extinct.samples$tipwardnode))){
        tip.name <- paste("ex_t", sample.index, sep="")
        new.fossil.names <- c(new.fossil.names, tip.name)
        tip <- list(edge=matrix(c(2,1),1,2), tip.label=tip.name, edge.length=0, Nnode=1)
        #So, since the fossil taxa are added after the tips, they should just be Ntip+i as their tip index:
        tip.desc.ext.list[[sample.index]] <- c(tip.desc.ext.list[[sample.index]], ntip + sample.index)
        class(tip) <- "phylo"
        position <- phy$edge.length[which(phy$edge[,2]==extinct.samples$tipwardnode[sample.index])] - extinct.samples$timefromrootwardnode[sample.index]
        phy <- bind.tree(phy, tip, extinct.samples$tipwardnode[sample.index], position=position)
        for(node.index in sequence(length(extinct.samples$tipwardnode))){
            tmp <- getMRCA(phy, tip.desc.ext.list[[node.index]])
            if(is.null(tmp)){
                extinct.samples$tipwardnode[node.index] <- tip.desc.ext.list[[node.index]]
            }else{
                extinct.samples$tipwardnode[node.index] <- tmp
            }
        }
    }
    
    if(strat.intervals == TRUE){
        strat <- hisse:::GetStratigraphicIntervals(phy, f)
    }else{
        #Step 3: Now process the k.sample table so that any added extinct taxa get listed as taxa whose MRCA is where subtending k.sample will be placed:
        for(sample.index in sequence(length(k.samples$tipwardnode))){
            samples.taxa <- tmp.k.samples[[sample.index]]
            mrca <- getMRCA(phy, samples.taxa)
            if(!is.null(mrca)){
                if(k.samples$fossiltype_long[sample.index] == "extinct_terminal" | k.samples$fossiltype_long[sample.index] == "extinct_internal"){
                    sister.taxa <- GetSister(phy, mrca)
                    if(sister.taxa > Ntip(phy)){
                        mrca <- mrca
                    }else{
                        mrca <- getMRCA(phy, c(phy$tip.label[sister.taxa], samples.taxa[1]))
                    }
                    all.desc <- phytools::getDescendants(phy, node=mrca)
                    tips.only <- all.desc[all.desc<(Ntip(phy)+1)]
                    tmp.k.samples[[sample.index]] <- phy$tip.label[tips.only]
                }else{
                    all.desc <- phytools::getDescendants(phy, node=mrca)
                    tips.only <- all.desc[all.desc<(Ntip(phy)+1)]
                    tmp.k.samples[[sample.index]] <- phy$tip.label[tips.only]
                }
            }else{
                if(k.samples$fossiltype_long[sample.index] == "extinct_terminal" | k.samples$fossiltype_long[sample.index] == "extinct_internal"){
                    new.extinct.sister <- GetSister(phy, samples.taxa)
                    tmp.k.samples[[sample.index]] <- c(samples.taxa, phy$tip.label[new.extinct.sister])
                }else{
                    tmp.k.samples[[sample.index]] <- samples.taxa
                }
            }
        }
        
        #Step 4: Now prune down to extant and sampled fossils
        fossil.tips <- geiger::is.extinct(phy,tol=.Machine$double.eps^.50)
        extant.tips <- phy$tip.label[!(phy$tip.label %in% fossil.tips)]
        sampled.tips <- c(extant.tips, new.fossil.names)
        unsampled.tips <- fossil.tips[!(fossil.tips %in% sampled.tips)]
        
        new.tree <- ape::drop.tip(phy, unsampled.tips)
        
        #Step 5: Now make k table by taking the list of taxon names, removing those that were culled, and finding the two that define MRCA of subtending sample.
        tmp <- c()
        for(sample.index in sequence(length(k.samples$tipwardnode))){
            samples.taxa <- tmp.k.samples[[sample.index]]
            samples.taxa <- samples.taxa[which(samples.taxa %in% new.tree$tip.label)]
            if(length(samples.taxa)==1){
                tmp <- rbind(tmp, c(samples.taxa[1], samples.taxa[1], k.samples$timefrompresent[sample.index]))
            }else{
                tmp <- rbind(tmp, c(samples.taxa[1], samples.taxa[length(samples.taxa)], k.samples$timefrompresent[sample.index]))
            }
        }
        
        k.samples <- data.frame(taxon1=tmp[,1], taxon2=tmp[,2], timefrompresent=tmp[,3], stringsAsFactors=FALSE)
        
        keep.root.samples=FALSE
        if(keep.root.samples == FALSE){
            root.edge.samples <- which(as.numeric(k.samples[,3]) > max(node.depth.edgelength(new.tree)))
            if(length(root.edge.samples)>0){
                k.samples <- k.samples[-root.edge.samples,]
            }
        }
    }
    obj <- list(phy=new.tree, k.samples=k.samples)
    return(obj)
}



######################################################################################################################################
######################################################################################################################################
### Utility functions for dealing with fossil data
######################################################################################################################################
######################################################################################################################################

GetKSampleMRCA <- function(phy, k.samples){
    k.sample.tip.no <- grep("Ksamp*", x=phy$tip.label)
    mrca.ksamp <- phy$edge[which(phy$edge[,2] %in% k.sample.tip.no),1]
    fix.info <- data.frame(node=mrca.ksamp, type="event", state=rep(0, length(mrca.ksamp)), stringsAsFactors=FALSE)
    return(fix.info)
}


AddKNodes <- function(phy, k.info){
    #Step 1: Get Descendants for a the MRCA of input taxa
    ntip <- Ntip(phy)
    split.times <- dateNodes(phy, rootAge=max(node.depth.edgelength(phy)))
    tip.desc.list <- as.list(1:length(k.info[,1]))
    mrca.list <- as.list(1:length(k.info[,1]))
    for(og.index in sequence(length(k.info[,1]))){
        if(k.info[og.index,1] == k.info[og.index,2]){
            tip.desc.list[[og.index]] <- k.info[og.index,1]
            mrca.list[[og.index]] <- which(phy$tip.label == k.info[og.index,1])
        }else{
            mrca <- getMRCA(phy, c(k.info[og.index,1], k.info[og.index,2]))
            all.desc <- phytools::getDescendants(phy, node=mrca)
            tip.desc.list[[og.index]] <- phy$tip.label[all.desc[all.desc<(Ntip(phy)+1)]]
            mrca.list[[og.index]] <- mrca
        }
    }
    
    #Step 2:
    for(sample.index in sequence(length(k.info[,1]))){
        tip.name <- paste("Ksamp", sample.index, sep="_")
        tip <- list(edge=matrix(c(2,1),1,2), tip.label=tip.name, edge.length=0, Nnode=1)
        #So, since the fossil taxa are added after the tips, they should just be Ntip+i as their tip index:
        tip.desc.list[[sample.index]] <- c(tip.desc.list[[sample.index]], tip.name)
        class(tip) <- "phylo"
        timefromtipward <- (as.numeric(k.info$timefrompresent[sample.index]) - split.times[which(names(split.times) == mrca.list[[sample.index]])])
        position <- timefromtipward
        phy <- bind.tree(phy, tip, mrca.list[[sample.index]], position=position)
        split.times <- dateNodes(phy, rootAge=max(node.depth.edgelength(phy)))
        for(node.index in sequence(length(k.info[,1]))){
            if(sample.index != node.index){
                tmp <- getMRCA(phy, tip.desc.list[[node.index]])
                if(!is.null(tmp)){
                    sister.taxa <- GetSister(phy, tmp)
                    if(length(sister.taxa) != 0){
                        is.ksample.sis <- strsplit(phy$tip.label[sister.taxa], "_")[[1]]
                        if(!is.na(is.ksample.sis[1])){
                            if(is.ksample.sis[1] == "Ksamp"){
                                samples.taxa <- tip.desc.list[[node.index]]
                                tip.desc.list[[node.index]] <- c(samples.taxa, phy$tip.label[sister.taxa])
                                mrca.list[[node.index]] <- getMRCA(phy, tip.desc.list[[node.index]])
                            }else{
                                samples.taxa <- tip.desc.list[[node.index]]
                                tip.desc.list[[node.index]] <- samples.taxa
                                mrca.list[[node.index]] <- getMRCA(phy, tip.desc.list[[node.index]])
                            }
                        }else{
                            samples.taxa <- tip.desc.list[[node.index]]
                            tip.desc.list[[node.index]] <- samples.taxa
                            mrca.list[[node.index]] <- getMRCA(phy, tip.desc.list[[node.index]])
                        }
                    }else{
                        samples.taxa <- tip.desc.list[[node.index]]
                        tip.desc.list[[node.index]] <- samples.taxa
                        mrca.list[[node.index]] <- getMRCA(phy, tip.desc.list[[node.index]])
                    }
                }else{
                    sister.taxa <- GetSister(phy, tip.desc.list[[node.index]])
                    if(length(sister.taxa) != 0){
                        is.ksample.sis <- strsplit(phy$tip.label[sister.taxa], "_")[[1]]
                        if(!is.na(is.ksample.sis[1])){
                            if(is.ksample.sis[1] == "Ksamp"){
                                samples.taxa <- tip.desc.list[[node.index]]
                                tip.desc.list[[node.index]] <- c(samples.taxa, phy$tip.label[sister.taxa])
                                mrca.list[[node.index]] <- getMRCA(phy, tip.desc.list[[node.index]])
                            }else{
                                samples.taxa <- tip.desc.list[[node.index]]
                                tip.desc.list[[node.index]] <- samples.taxa
                                tmp <- getMRCA(phy, tip.desc.list[[node.index]])
                                if(!is.null(tmp)){
                                    mrca.list[[node.index]] <- tmp
                                }
                            }
                        }else{
                            samples.taxa <- tip.desc.list[[node.index]]
                            tip.desc.list[[node.index]] <- samples.taxa
                            tmp <- getMRCA(phy, tip.desc.list[[node.index]])
                            if(!is.null(tmp)){
                                mrca.list[[node.index]] <- tmp
                            }
                        }
                    }else{
                        samples.taxa <- tip.desc.list[[node.index]]
                        tip.desc.list[[node.index]] <- samples.taxa
                        if(!is.null(tmp)){
                            mrca.list[[node.index]] <- tmp
                        }
                    }
                }
            }else{
                tmp <- getMRCA(phy, tip.desc.list[[node.index]])
                mrca.list[[node.index]] <- tmp
            }
        }
    }
    return(phy)
}


AddKData <- function(data, k.samples, muhisse=FALSE){
    for(event.index in sequence(length(k.samples$taxon1))){
        if(muhisse == TRUE){
            new.data <- c(paste("Ksamp", event.index, sep="_"), k.samples$state1[event.index], k.samples$state2[event.index])
        }else{
            new.data <- c(paste("Ksamp", event.index, sep="_"), rep(k.samples$state[event.index],2))
        }
        data <- rbind(data, new.data)
    }
    return(data)
}



#data to play with:
# ntax=200
# true.psi=0.1
# try(sim.tab <- hisse::SimulateHisse(turnover=c(0.25,0.25), eps=rep(0.25,2), max.taxa=ntax, x0=0, transition.rates=matrix(c(NA, 0.005, 0.005, NA), nrow=2), nstart=2))
# phy <- SimToPhylo(sim.tab, include.extinct=TRUE)
# f <- hisse:::GetFossils(phy, psi=true.psi)


#CheckPlot <- function(phy, extinct.samples, k.samples){
#    split.times <- dateNodes(phy, rootAge=max(node.depth.edgelength(phy)))
#    plot(phy, cex=.5)
#    abline(v=as.numeric(extinct.samples[,5]), col="blue")
#    abline(v=c(max(split.times) - as.numeric(k.samples[,3])), col="red")
#}

#phy <- pp$phy
#k.info <- pp$k.samples
#test <- AddKNodes(pp$phy, pp$k.samples)
#extinct.samples <- f[which(f$fossiltype_long=="extinct_terminal" | f$fossiltype_long=="extinct_internal"),]
#extinct.samples <- extinct.samples[which(extinct.samples$has_sampled_descendant == FALSE),]
#CheckPlot(test, extinct.samples=extinct.samples, k.samples=k.info)

#dat.tab <- hisse:::OrganizeDataMiSSE(phy=pp$phy, f=1, hidden.states=1)
#These are all inputs for generating starting values:
#fossil.taxa <- which(dat.tab$branch.type == 1)
#fossil.ages <- dat.tab$TipwardAge[which(dat.tab$branch.type == 1)]





