

######################################################################################################################################
######################################################################################################################################
### Simulator for birth-death models with fossils and sampled lineages
######################################################################################################################################
######################################################################################################################################

GetSister <- function(tree, node, mode=c("number","label")){
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
    
    fossils <- fossils[order(fossils$timefrompresent),]
    
    return(fossils)
}


GetStratigraphicIntervals <- function(phy, f) {
    f$edge_root_tip <- paste0(f$rootwardnode, "_", f$tipwardnode)
    strat_intervals <- data.frame(edge_root_tip=unique(f$edge_root_tip), rootwardnode=NA, tipwardnode=NA, startingtimefromroot=NA, endingtimefromroot=NA, intervallength=0, tipwardendisextant=FALSE)
    for (i in sequence(nrow(strat_intervals))) {
        matching_f <- subset(f, f$edge_root_tip == strat_intervals$edge_root_tip[i])
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


ProcessSimSample <- function(phy, f){
    
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

    split.times <- paleotree::dateNodes(phy, rootAge=max(node.depth.edgelength(phy)))
    #Step 3: Now process the k.sample table so that any added extinct taxa get listed as taxa whose MRCA is where subtending k.sample will be placed:
    for(sample.index in sequence(length(k.samples$tipwardnode))){
        #obtains the list of taxa whose subtending branch is where the k sample is located:
        samples.taxa <- tmp.k.samples[[sample.index]]
        #now get the MRCA of of these samples:
        mrca <- getMRCA(phy, samples.taxa)
        #if not null we are in the tree:
        if(!is.null(mrca)){
            #Are they on an eventually extinct branch?
            if(k.samples$fossiltype_long[sample.index] == "extinct_terminal" | k.samples$fossiltype_long[sample.index] == "extinct_internal"){
                #Get the sister taxa for this set of taxa:
                sister.taxa <- GetSister(phy, mrca)
                #Is the sister taxon a node or a tip?
                if(sister.taxa > Ntip(phy)){
                    mrca <- mrca
                }else{
                    mrca.new <- getMRCA(phy, c(phy$tip.label[sister.taxa], samples.taxa[1]))
                    if(k.samples$timefrompresent[sample.index] < unname(split.times[mrca.new])){
                        mrca <- mrca
                    }else{
                        mrca <- mrca.new
                    }
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
    
    obj <- list(phy=new.tree, k.samples=k.samples)
    return(obj)
}


ProcessSimStrat <- function(phy, f){
    ## STEP 1: Prune f table to only include those that represent intervals:
    f$edge_root_tip <- paste0(f$rootwardnode, "_", f$tipwardnode)
    unique_edge_root_tip <- unique(f$edge_root_tip)
    f.red <- as.data.frame(matrix(ncol=9, nrow=0))
    colnames(f.red) <- c("rootwardnode", "tipwardnode", "timefromroot", "timefromrootwardnode", "timefrompresent", "fossiltype_short", "fossiltype_long", "order_on_branch", "has_sampled_descendant")
    for (i in sequence(length(unique_edge_root_tip))) {
        matching_f <- subset(f, f$edge_root_tip == unique_edge_root_tip[i])
        if(dim(matching_f)[1] > 1){
            #ok get y_i:
            tmp.f.red <- matching_f[which.min(as.numeric(matching_f$timefromroot)),]
            #now get o_i:
            tmp.f.red <- rbind(tmp.f.red, matching_f[which.max(as.numeric(matching_f$timefromroot)),])
        }else{
            #so here, y_i is either a tip or this is a single k sample on a single edge.
            tmp.f.red <- matching_f
        }
        f.red <- rbind(f.red, tmp.f.red)
    }
    f.red <- f.red[order(as.numeric(f.red$timefromroot),decreasing=TRUE),]
    
    ## STEP 2: Now process the remaining stuff:
    points.added <- ProcessSimSample(phy, f.red)
    split.times <- paleotree::dateNodes(points.added$phy, rootAge=max(node.depth.edgelength(points.added$phy)))
    
    ## STEP 3: Now take this information and format for stat range evolution:
    # taxon1 taxon2 timefrompresentrootward timefrompresenttipward type
    strat.intervals <- as.data.frame(matrix(ncol=5, nrow=0))
    colnames(strat.intervals) <- c("taxon1", "taxon2", "timefrompresentroot", "timefrompresenttip", "type")
    
    k.samples <- points.added$k.samples
    k.samples$mrca_taxa <- paste0(k.samples$taxon1, "_", k.samples$taxon2)
    unique_mrca_taxa <- unique(k.samples$mrca_taxa)
    for(mrca.index in 1:length(unique_mrca_taxa)){
        matching_rows <- subset(k.samples, k.samples$mrca_taxa == unique_mrca_taxa[mrca.index])
        if(dim(matching_rows)[1] > 1){
            if(dim(matching_rows)[1] > 2){
                #here means type is interval is a range and there are multiple ranges on this edge:
                info.from.previous <- as.data.frame(matrix(ncol=2, nrow=dim(matching_rows)[1]))
                for(tmp.index in 1:dim(matching_rows)[1]){
                    tmp <- f.red[which(f.red$timefrompresent == matching_rows$timefrompresent[tmp.index]),]
                    info.from.previous[tmp.index,] <- c(tmp$edge_root_tip, tmp$timefrompresent)
                }
                unique_edge_root_tip <- unique(info.from.previous[,1])
                for(unique.index in 1:length(unique_edge_root_tip)){
                    matching_tmp <- subset(info.from.previous, info.from.previous[,1] == unique_edge_root_tip[unique.index])
                    if(dim(matching_tmp)[1] > 1){
                        tmp <- as.data.frame(matrix(ncol=5, nrow=1))
                        colnames(tmp) <- c("taxon1", "taxon2", "timefrompresentroot", "timefrompresenttip", "type")
                        tmp[1,] <- c(matching_rows$taxon1[1], matching_rows$taxon2[1], max(as.numeric(matching_tmp[,2])), min(as.numeric(matching_tmp[,2])), "R")
                        strat.intervals <- rbind(strat.intervals, tmp)
                    }else{
                        #Singleton or Tip range
                        if(unique.index == 1){
                            if(matching_rows$taxon1[1] == matching_rows$taxon2[1]){
                                #Tip Range
                                tmp <- as.data.frame(matrix(ncol=5, nrow=1))
                                colnames(tmp) <- c("taxon1", "taxon2", "timefrompresentroot", "timefrompresenttip", "type")
                                tmp[1,] <- c(matching_rows$taxon1[1], matching_rows$taxon2[1], as.numeric(matching_tmp[2]), split.times[which(points.added$phy$tip.label == matching_rows$taxon1[1])], "R")
                                strat.intervals <- rbind(strat.intervals, tmp)
                            }else{
                                #Singleton
                                tmp <- as.data.frame(matrix(ncol=5, nrow=1))
                                colnames(tmp) <- c("taxon1", "taxon2", "timefrompresentroot", "timefrompresenttip", "type")
                                tmp[1,] <- c(matching_rows$taxon1[1], matching_rows$taxon2[1], as.numeric(matching_tmp[2]), as.numeric(matching_tmp[2]), "S")
                                strat.intervals <- rbind(strat.intervals, tmp)
                            }
                        }else{
                            #Singleton behind another interval
                            tmp <- as.data.frame(matrix(ncol=5, nrow=1))
                            colnames(tmp) <- c("taxon1", "taxon2", "timefrompresentroot", "timefrompresenttip", "type")
                            tmp[1,] <- c(matching_rows$taxon1[1], matching_rows$taxon2[1], as.numeric(matching_tmp[2]), as.numeric(matching_tmp[2]), "S")
                            strat.intervals <- rbind(strat.intervals, tmp)
                        }
                    }
                }
            }else{
                #Interval is a range and there is only a single range on this edge:
                tmp <- as.data.frame(matrix(ncol=5, nrow=1))
                colnames(tmp) <- c("taxon1", "taxon2", "timefrompresentroot", "timefrompresenttip", "type")
                tmp[1,] <- c(matching_rows$taxon1[1], matching_rows$taxon2[1], max(as.numeric(matching_rows$timefrompresent)), min(as.numeric(matching_rows$timefrompresent)), "R")
                strat.intervals <- rbind(strat.intervals, tmp)
            }
        }else{
            #Interval type is a singleton fossil or tip range:
            if(matching_rows$taxon1[1] == matching_rows$taxon2[1]){
                #Tip range:
                tmp <- as.data.frame(matrix(ncol=5, nrow=1))
                colnames(tmp) <- c("taxon1", "taxon2", "timefrompresentroot", "timefrompresenttip", "type")
                tmp[1,] <- c(matching_rows$taxon1[1], matching_rows$taxon2[1], max(as.numeric(matching_rows$timefrompresent)), split.times[which(points.added$phy$tip.label == matching_rows$taxon1[1])], "R")
                strat.intervals <- rbind(strat.intervals, tmp)
            }else{
                #Singleton fossil:
                tmp <- as.data.frame(matrix(ncol=5, nrow=1))
                colnames(tmp) <- c("taxon1", "taxon2", "timefrompresentroot", "timefrompresenttip", "type")
                tmp[1,] <- c(matching_rows$taxon1[1], matching_rows$taxon2[1], max(as.numeric(matching_rows$timefrompresent)), min(as.numeric(matching_rows$timefrompresent)), "S")
                strat.intervals <- rbind(strat.intervals, tmp)
            }
        }
    }
    obj <- list(phy=points.added$phy, strat.intervals=strat.intervals)
    
    return(obj)
}


######################################################################################################################################
######################################################################################################################################
### Utility functions for dealing with fossil data
######################################################################################################################################
######################################################################################################################################

GetKSampleMRCA <- function(phy, k.samples, strat.intervals=FALSE){
    k.sample.tip.no <- grep("Ksamp*", x=phy$tip.label)
    mrca.ksamp <- phy$edge[which(phy$edge[,2] %in% k.sample.tip.no),1]
    if(strat.intervals == TRUE){
        fix.info <- data.frame(node=mrca.ksamp, type="interval", state=rep(0, length(mrca.ksamp)), stringsAsFactors=FALSE)
    }else{
        fix.info <- data.frame(node=mrca.ksamp, type="event", state=rep(0, length(mrca.ksamp)), stringsAsFactors=FALSE)
    }
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


######################################################################################################################################
######################################################################################################################################
### Utility functions for dealing with stratigraphic intervals
######################################################################################################################################
######################################################################################################################################

GetStratInfo <- function(strat.intervals){
    #Step 1: Get the easy stuff first...
    obj <- NULL
    obj$k <- (length(strat.intervals$type[which(strat.intervals$type == "R")]) * 2) + length(strat.intervals$type[which(strat.intervals$type == "S")]) - length(which(as.numeric(strat.intervals$timefrompresenttip) < .Machine$double.eps^.50))
    obj$l_s <- sum(as.numeric(strat.intervals$timefrompresentroot) - as.numeric(strat.intervals$timefrompresenttip))
    
    #Step 2: Get the hard stuff last...
    intervening.intervals.tmp <- NULL
    strat.intervals$mrca_taxa <- paste0(strat.intervals$taxon1, "_", strat.intervals$taxon2)
    unique_mrca_taxa <- unique(strat.intervals$mrca_taxa)
    for(unique.index in 1:length(unique_mrca_taxa)){
        matching_rows <- subset(strat.intervals, strat.intervals$mrca_taxa == unique_mrca_taxa[unique.index])
        if(dim(matching_rows)[1] > 1){
            #Order so that they are decreasing -- in theory could have several intervals down a long edge:
            matching_rows <- matching_rows[order(as.numeric(matching_rows$timefrompresenttip),decreasing=FALSE),]
            intervening.intervals <- c()
            for(index in 1:(dim(matching_rows)[1]-1)){
                intervening.intervals.tmp <- rbind(intervening.intervals.tmp, c(matching_rows$taxon1[index], matching_rows$taxon2[index], as.numeric(matching_rows$timefrompresenttip[index+1]), as.numeric(matching_rows$timefrompresentroot[index]), matching_rows$type[index+1]))
            }
        }
    }
    if(!is.null(intervening.intervals.tmp)){
        intervening.intervals <- data.frame(taxon1=intervening.intervals.tmp[,1], taxon2=intervening.intervals.tmp[,2], ya=as.numeric(intervening.intervals.tmp[,3]), o_i=as.numeric(intervening.intervals.tmp[,4]), ya_type=intervening.intervals.tmp[,5], stringsAsFactors=FALSE)
        obj$intervening.intervals <- intervening.intervals
    }
    return(obj)
}


GetIntervalToK <- function(strat.intervals, intervening.intervals=NULL){
    #Step 1: Reduce table down to just R type fossils, R meaning ranges.
    intervals.only <- subset(strat.intervals, strat.intervals$type == "R")
    #Step 2: Loop down table and organize everything to match k sample input:
    tmp <- c()
    for(row.index in 1:dim(intervals.only)[1]){
        tmp <- rbind(tmp, c(intervals.only$taxon1[row.index], intervals.only$taxon2[row.index], as.numeric(intervals.only$timefrompresenttip[row.index])))
        tmp <- rbind(tmp, c(intervals.only$taxon1[row.index], intervals.only$taxon2[row.index], as.numeric(intervals.only$timefrompresentroot[row.index])))
    }
    
    if(!is.null(intervening.intervals)){
        #check to see if any S types in intervening.intervals should be retained so we can demarcate intervening intervals:
        extra.k.samples <- subset(intervening.intervals, intervening.intervals$ya_type == "S")
        #If zero it will just fill row with NAs.
        if(dim(extra.k.samples)[1] > 0){
            #loop over the the set of S types:
            for(extra.row.index in 1:dim(extra.k.samples)[1]){
                #Match the tipward time first. If matches that means the o_i of intervening interval is an S fossil (rare, but possible).
                time.match <- which(round(as.numeric(tmp[,3]),10) %in% round(extra.k.samples$o_i[extra.row.index],10))
                if(length(time.match) > 0){
                    #Means that the o_i is an R type range.
                    tmp <- rbind(tmp, c(extra.k.samples$taxon1[extra.row.index], extra.k.samples$taxon2[extra.row.index], as.numeric(extra.k.samples$ya[extra.row.index])))
                }else{
                    #Means that the o_i is an S type range. So retain both.
                    tmp <- rbind(tmp, c(extra.k.samples$taxon1[extra.row.index], extra.k.samples$taxon2[extra.row.index], as.numeric(extra.k.samples$o_i[extra.row.index])))
                    tmp <- rbind(tmp, c(extra.k.samples$taxon1[extra.row.index], extra.k.samples$taxon2[extra.row.index], as.numeric(extra.k.samples$ya[extra.row.index])))
                }
            }
        }
        
        #Now let us check some of the ya R intervening intervals to check if the o_i is an S type fossil:
        extra.k.samples2 <- subset(intervening.intervals, intervening.intervals$ya_type == "R")
        if(dim(extra.k.samples2)[1] > 0){
            for(extra.row.index in 1:dim(extra.k.samples2)[1]){
                #Match the tipward time first. If matches that means front part of intervening interval is already in our set of retained S fossils
                time.match <- which(round(as.numeric(tmp[,3]),10) %in% round(extra.k.samples2$o_i[extra.row.index],10))
                if(length(time.match) == 0){
                    #if it does not match then it is not included and we should include it.
                    tmp <- rbind(tmp, c(extra.k.samples2$taxon1[extra.row.index], extra.k.samples2$taxon2[extra.row.index], as.numeric(extra.k.samples2$o_i[extra.row.index])))
                }
            }
        }
    }
    
    k.samples <- data.frame(taxon1=as.character(tmp[,1]), taxon2=as.character(tmp[,2]), timefrompresent=as.numeric(tmp[,3]), stringsAsFactors=FALSE)
    #Do not need end of interval if it is a tip.
    tip.samples <- which(as.numeric(k.samples$timefrompresent) < .Machine$double.eps^.50)
    if(length(tip.samples) > 0){
        k.samples <- k.samples[-tip.samples,]
    }
    k.samples <- k.samples[order(as.numeric(k.samples[,3]),decreasing=FALSE),]
    return(k.samples)
}


