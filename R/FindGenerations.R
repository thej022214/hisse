## Common function used by multiple methods.

#Note that James Boyko came up with the FindGenerations code, so the insanity of this Rumsfeldian speak is all him....
FindGenerations <- function(phy){
    generation <- list()
    known <- 1:Ntip(phy)
    unknown <- phy$edge[,1]
    needed <- phy$edge[,2]
    root <- min(unknown)
    i <- 1
    repeat{
        knowable <- unknown[needed %in% known]
        knowable <- knowable[duplicated(knowable)]
        generation[[i]] <-  knowable
        
        known <- c(known, knowable)
        needed <- needed[!unknown %in% knowable]
        unknown <- unknown[!unknown %in% knowable]
        i <- i + 1
        if (any(root == knowable)) break
    }
    res <- generation
    return(res)
}
