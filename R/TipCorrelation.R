
TipCorrelation <- function(phy, tip.rate, trait, log=TRUE, remove.cherries=TRUE, scaled=TRUE, positivise=TRUE, use.lmorigin=TRUE) {
  if(!inherits(c(tip.rate, trait), "numeric")) {
    stop("Arguments 'tip.rate' and 'trait' should be of class numeric.")
  }
  if(!inherits(phy, "phylo")) {
    stop("Arguments 'phy' should be of class phylo.")
  }
  if(!setequal(phy$tip.label, names(tip.rate))) {
    stop("Names in tree and tip.rate don't match.")
  }
  if(!setequal(phy$tip.label, names(trait))) {
    stop("Names in tree and trait don't match.")
  }
  phy$node.label <- NULL
  # Few lines to guarantee data is in the same order #-------------
  is_tip <- phy$edge[,2] <= length(phy$tip.label)
  ordered_tips <- phy$edge[is_tip, 2]
  right_order <- as.character(phy$tip.label[ordered_tips])
  trait <- trait[match(as.character(right_order), as.character(names(trait)))]
  tip.rate <- tip.rate[match(as.character(right_order), as.character(names(tip.rate)))]
  #-------
  t0 = trait
  r0 = tip.rate
  if(log) {
    t0 <- log(t0)
    r0 <- log(r0)    
  } 
  # PICs #-------------
  t0Pic <- ape::pic(t0, phy, scaled = scaled)
  r0Pic <- ape::pic(r0, phy, scaled = scaled)
  if(remove.cherries) {
  # Identifying cherries #-------------
  internals <- ape::Ntip(phy) + sequence(ape::Nnode(phy))
  ndescendants <- rep(NA, length(internals))
  for (i in seq_along(internals)) {
    ndescendants[i] <- length(phytools::getDescendants(phy, node=internals[i]))
  }
  cherry_internals <- internals[ndescendants==2]
  # Remove cherries #-------------
  t0Pic <- t0Pic[!(as.numeric(names(t0Pic)) %in% cherry_internals)]
  r0Pic <- r0Pic[!(as.numeric(names(r0Pic)) %in% cherry_internals)]
  }
  if(positivise) {
  # Positivizing PICs #-------------
  r0Pic[which(t0Pic < 0)] <- -1 * r0Pic[which(t0Pic < 0)]
  t0Pic <- abs(t0Pic)
  }
  trait = t0Pic
  tip.rate = r0Pic
  if(use.lmorigin){
    picModel <- lm(tip.rate ~ 0 + trait)
    result <- ape::lmorigin(picModel, picModel$model, origin=TRUE, nperm=999)
  } else {
    result <- lm(tip.rate ~ trait)
  }
  all_results <- list(result, tip.rate, trait)
  names(all_results) <- c("correlation","tip.rate PIC","trait PIC")
  return(all_results)
}


