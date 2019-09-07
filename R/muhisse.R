

######################################################################################################################################
######################################################################################################################################
### MuHiSSE -- Expanded set of MuSSE models for examining diversification in relation to multiple trait evolution
######################################################################################################################################
######################################################################################################################################

MuHiSSE <- function(phy, data, f=c(1,1,1,1), turnover=c(1,2,3,4), eps=c(1,2,3,4), hidden.states=FALSE, trans.rate=NULL, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, sann=FALSE, sann.its=10000, bounded.search=TRUE, max.tol=.Machine$double.eps^.50, starting.vals=NULL, turnover.upper=10000, eps.upper=3, trans.upper=100, restart.obj=NULL, ode.eps=0, dt.threads=1){
    
    ## Temporary fix for the current BUG:
    if( !is.null(phy$node.label) ) phy$node.label <- NULL
    
    if(!is.null(root.p)) {
        if(hidden.states ==TRUE){
            if( length( root.p ) == 4 ){
                root.p <- rep(root.p, 4)
                root.p <- root.p / sum(root.p)
                warning("For hidden states, you need to specify the root.p for all hidden states. We have adjusted it so that there's equal chance for among all hidden states.")
            } else{
                root.p.new <- numeric(32)
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
    
    if(any(f == 0)){
        f[which(f==0)] <- 1
    }
    
    if(!root.type == "madfitz" & !root.type == "herr_als"){
        stop("Check that you specified a proper root.type option. Options are 'madfitz' or 'herr_als'. See help for more details.", call.=FALSE)
    }
    
    if(is.null(trans.rate)){
        stop("Rate matrix needed. See TransMatMakerMuHiSSE() to create one.")
    }
    
    if(hidden.states == TRUE & dim(trans.rate)[1]<4){
        stop("You chose a hidden state but this is not reflected in the transition matrix")
    }
    
    ## Return error message if the data is not in the correct format.
    if( !inherits(data, what = c("matrix", "data.frame")) ){
        stop("'data' needs to be a matrix or data.frame with 3 columns. See help.")
    }
    if( !ncol( data ) == 3 ){
        stop("'data' needs to be a matrix or data.frame with 3 columns. See help.")
    }
    
    ## Check if 'hidden.states' parameter is congruent with the turnover vector:
    if( length(turnover) > 4 & !hidden.states ){
        stop("Turnover has more than 4 elements but 'hidden.states' was set to FALSE. Please set 'hidden.states' to TRUE if the model include more than one rate class.")
    }
    if( length(turnover) == 4 & hidden.states ){
        stop("Turnover has only 4 elements but 'hidden.states' was set to TRUE. Please set 'hidden.states' to FALSE if the model does not include hidden rate classes.")
    }
    
    pars <- numeric(384)
    
    if(dim(trans.rate)[2]==4){
        rate.cats <- 1
        pars.tmp <- turnover
        eps.tmp <- eps
        eps.tmp[which(eps.tmp > 0)] = (eps.tmp[which( eps.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, eps.tmp)
        trans.tmp <- c(trans.rate["(00)", "(01)"], trans.rate["(00)", "(10)"], trans.rate["(00)", "(11)"], trans.rate["(01)", "(00)"], trans.rate["(01)", "(10)"], trans.rate["(01)", "(11)"], trans.rate["(10)", "(00)"], trans.rate["(10)", "(01)"], trans.rate["(10)", "(11)"], trans.rate["(11)", "(00)"], trans.rate["(11)", "(01)"], trans.rate["(11)", "(10)"])
        trans.tmp[which(trans.tmp > 0)] = (trans.tmp[which(trans.tmp > 0)] + max(pars.tmp))
        category.rates.unique <- 0
        pars.tmp <- c(pars.tmp, trans.tmp)
        pars[1:20] <- pars.tmp
    }
    
    if(dim(trans.rate)[2]==8){
        rate.cats <- 2
        pars.tmp <- turnover
        eps.tmp <- eps
        eps.tmp[which(eps.tmp > 0)] = (eps.tmp[which( eps.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, eps.tmp)
        for.late.adjust <- max(pars.tmp)
        rows <- c("(00A)", "(00A)", "(00A)", "(01A)", "(01A)", "(01A)", "(10A)", "(10A)", "(10A)", "(11A)", "(11A)", "(11A)", "(00B)", "(00B)", "(00B)", "(01B)", "(01B)", "(01B)", "(10B)", "(10B)", "(10B)", "(11B)", "(11B)", "(11B)")
        cols <- c("(01A)", "(10A)", "(11A)", "(00A)", "(10A)", "(11A)", "(00A)", "(01A)", "(11A)", "(00A)", "(01A)", "(10A)", "(01B)", "(10B)", "(11B)", "(00B)", "(10B)", "(11B)", "(00B)", "(01B)", "(11B)", "(00B)", "(01B)", "(10B)")
        trans.tmp <- trans.rate[cbind(rows,cols)]
        trans.tmp[which(trans.tmp > 0)] = (trans.tmp[which(trans.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, trans.tmp)
        rows <- c("(00A)", "(01A)", "(10A)", "(11A)", "(00B)", "(01B)", "(10B)", "(11B)")
        cols <- c("(00B)", "(01B)", "(10B)", "(11B)", "(00A)", "(01A)", "(10A)", "(11A)")
        category.tmp <- trans.rate[cbind(rows,cols)]
        category.tmp[category.tmp==0] <- NA
        category.rate.shift <- category.tmp + for.late.adjust
        category.rate.shift[is.na(category.rate.shift)] <- 0
        category.rate.shiftA <- c(category.rate.shift[1], rep(0,6), category.rate.shift[2], rep(0,6), category.rate.shift[3], rep(0,6), category.rate.shift[4], rep(0,6))
        category.rate.shiftB <- c(category.rate.shift[5], rep(0,6), category.rate.shift[6], rep(0,6), category.rate.shift[7], rep(0,6), category.rate.shift[8], rep(0,6))
        pars.tmp <- c(turnover[1:4], eps.tmp[1:4], trans.tmp[1:12], category.rate.shiftA, turnover[5:8], eps.tmp[5:8], trans.tmp[13:24], category.rate.shiftB)
        pars[1:length(pars.tmp)] <- pars.tmp
    }
    
    if(dim(trans.rate)[2]==12){
        rate.cats <- 3
        pars.tmp <- turnover
        eps.tmp <- eps
        eps.tmp[which(eps.tmp > 0)] = (eps.tmp[which( eps.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, eps.tmp)
        for.late.adjust <- max(pars.tmp)
        rows <- c("(00A)", "(00A)", "(00A)", "(01A)", "(01A)", "(01A)", "(10A)", "(10A)", "(10A)", "(11A)", "(11A)", "(11A)", "(00B)", "(00B)", "(00B)", "(01B)", "(01B)", "(01B)", "(10B)", "(10B)", "(10B)", "(11B)", "(11B)", "(11B)", "(00C)", "(00C)", "(00C)", "(01C)", "(01C)", "(01C)", "(10C)", "(10C)", "(10C)", "(11C)", "(11C)", "(11C)")
        cols <- c("(01A)", "(10A)", "(11A)", "(00A)", "(10A)", "(11A)", "(00A)", "(01A)", "(11A)", "(00A)", "(01A)", "(10A)", "(01B)", "(10B)", "(11B)", "(00B)", "(10B)", "(11B)", "(00B)", "(01B)", "(11B)", "(00B)", "(01B)", "(10B)", "(01C)", "(10C)", "(11C)", "(00C)", "(10C)", "(11C)", "(00C)", "(01C)", "(11C)", "(00C)", "(01C)", "(10C)")
        trans.tmp <- trans.rate[cbind(rows,cols)]
        trans.tmp[which(trans.tmp > 0)] = (trans.tmp[which(trans.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, trans.tmp)
        rows <- c("(00A)", "(00A)", "(01A)", "(01A)", "(10A)", "(10A)", "(11A)", "(11A)", "(00B)", "(00B)", "(01B)", "(01B)", "(10B)", "(10B)", "(11B)", "(11B)", "(00C)", "(00C)", "(01C)", "(01C)", "(10C)", "(10C)", "(11C)", "(11C)")
        cols <- c("(00B)", "(00C)", "(01B)", "(01C)", "(10B)", "(10C)", "(11B)", "(11C)", "(00A)", "(00C)", "(01A)", "(01C)", "(10A)", "(10C)", "(11A)", "(11C)", "(00A)", "(00B)", "(01A)", "(01B)", "(10A)", "(10B)", "(11A)", "(11B)")
        category.tmp <- trans.rate[cbind(rows,cols)]
        category.tmp[category.tmp==0] <- NA
        category.rate.shift <- category.tmp + for.late.adjust
        category.rate.shift[is.na(category.rate.shift)] <- 0
        category.rate.shiftA <- c(category.rate.shift[1:2], rep(0,5), category.rate.shift[3:4], rep(0,5), category.rate.shift[5:6], rep(0,5), category.rate.shift[7:8], rep(0,5))
        category.rate.shiftB <- c(category.rate.shift[9:10], rep(0,5), category.rate.shift[11:12], rep(0,5), category.rate.shift[13:14], rep(0,5), category.rate.shift[15:16], rep(0,5))
        category.rate.shiftC <- c(category.rate.shift[17:18], rep(0,5), category.rate.shift[19:20], rep(0,5), category.rate.shift[21:22], rep(0,5), category.rate.shift[23:24], rep(0,5))
        pars.tmp <- c(turnover[1:4], eps.tmp[1:4], trans.tmp[1:12], category.rate.shiftA, turnover[5:8], eps.tmp[5:8], trans.tmp[13:24], category.rate.shiftB, turnover[9:12], eps.tmp[9:12], trans.tmp[25:36], category.rate.shiftC)
        pars[1:length(pars.tmp)] <- pars.tmp
    }
    
    if(dim(trans.rate)[2]==16){
        rate.cats <- 4
        pars.tmp <- turnover
        eps.tmp <- eps
        eps.tmp[which(eps.tmp > 0)] = (eps.tmp[which( eps.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, eps.tmp)
        for.late.adjust <- max(pars.tmp)
        rows <- c("(00A)", "(00A)", "(00A)", "(01A)", "(01A)", "(01A)", "(10A)", "(10A)", "(10A)", "(11A)", "(11A)", "(11A)", "(00B)", "(00B)", "(00B)", "(01B)", "(01B)", "(01B)", "(10B)", "(10B)", "(10B)", "(11B)", "(11B)", "(11B)", "(00C)", "(00C)", "(00C)", "(01C)", "(01C)", "(01C)", "(10C)", "(10C)", "(10C)", "(11C)", "(11C)", "(11C)", "(00D)", "(00D)", "(00D)", "(01D)", "(01D)", "(01D)", "(10D)", "(10D)", "(10D)", "(11D)", "(11D)", "(11D)")
        cols <- c("(01A)", "(10A)", "(11A)", "(00A)", "(10A)", "(11A)", "(00A)", "(01A)", "(11A)", "(00A)", "(01A)", "(10A)", "(01B)", "(10B)", "(11B)", "(00B)", "(10B)", "(11B)", "(00B)", "(01B)", "(11B)", "(00B)", "(01B)", "(10B)", "(01C)", "(10C)", "(11C)", "(00C)", "(10C)", "(11C)", "(00C)", "(01C)", "(11C)", "(00C)", "(01C)", "(10C)", "(01D)", "(10D)", "(11D)", "(00D)", "(10D)", "(11D)", "(00D)", "(01D)", "(11D)", "(00D)", "(01D)", "(10D)")
        trans.tmp <- trans.rate[cbind(rows,cols)]
        trans.tmp[which(trans.tmp > 0)] = (trans.tmp[which(trans.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, trans.tmp)
        rows <- c("(00A)", "(00A)", "(00A)", "(01A)", "(01A)", "(01A)", "(10A)", "(10A)", "(10A)", "(11A)", "(11A)", "(11A)", "(00B)", "(00B)", "(00B)", "(01B)", "(01B)", "(01B)", "(10B)", "(10B)", "(10B)", "(11B)", "(11B)", "(11B)", "(00C)", "(00C)", "(00C)", "(01C)", "(01C)", "(01C)", "(10C)", "(10C)", "(10C)", "(11C)", "(11C)", "(11C)", "(00D)", "(00D)", "(00D)", "(01D)", "(01D)", "(01D)", "(10D)", "(10D)", "(10D)", "(11D)", "(11D)", "(11D)")
        cols <- c("(00B)", "(00C)", "(00D)", "(01B)", "(01C)", "(01D)", "(10B)", "(10C)", "(10D)", "(11B)", "(11C)", "(11D)", "(00A)", "(00C)", "(00D)", "(01A)", "(01C)", "(01D)", "(10A)", "(10C)", "(10D)", "(11A)", "(11C)", "(11D)", "(00A)", "(00B)", "(00D)", "(01A)", "(01B)", "(01D)", "(10A)", "(10B)", "(10D)", "(11A)", "(11B)", "(11D)", "(00A)", "(00B)", "(00C)", "(01A)", "(01B)", "(01C)", "(10A)", "(10B)", "(10C)", "(11A)", "(11B)", "(11C)")
        category.tmp <- trans.rate[cbind(rows,cols)]
        category.tmp[category.tmp==0] <- NA
        category.rate.shift <- category.tmp + for.late.adjust
        category.rate.shift[is.na(category.rate.shift)] <- 0
        category.rate.shiftA <- c(category.rate.shift[1:3], rep(0,4), category.rate.shift[4:6], rep(0,4), category.rate.shift[7:9], rep(0,4), category.rate.shift[10:12], rep(0,4))
        category.rate.shiftB <- c(category.rate.shift[13:15], rep(0,4), category.rate.shift[16:18], rep(0,4), category.rate.shift[19:21], rep(0,4), category.rate.shift[22:24], rep(0,4))
        category.rate.shiftC <- c(category.rate.shift[25:27], rep(0,4), category.rate.shift[28:30], rep(0,4), category.rate.shift[31:33], rep(0,4), category.rate.shift[34:36], rep(0,4))
        category.rate.shiftD <- c(category.rate.shift[37:39], rep(0,4), category.rate.shift[40:42], rep(0,4), category.rate.shift[43:45], rep(0,4), category.rate.shift[46:48], rep(0,4))
        pars.tmp <- c(turnover[1:4], eps.tmp[1:4], trans.tmp[1:12], category.rate.shiftA, turnover[5:8], eps.tmp[5:8], trans.tmp[13:24], category.rate.shiftB, turnover[9:12], eps.tmp[9:12], trans.tmp[25:36], category.rate.shiftC, turnover[13:16], eps.tmp[13:16], trans.tmp[37:48], category.rate.shiftD)
        pars[1:length(pars.tmp)] <- pars.tmp
    }
    
    if(dim(trans.rate)[2]==20){
        rate.cats <- 5
        pars.tmp <- turnover
        eps.tmp <- eps
        eps.tmp[which(eps.tmp > 0)] = (eps.tmp[which( eps.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, eps.tmp)
        for.late.adjust <- max(pars.tmp)
        rows <- c("(00A)", "(00A)", "(00A)", "(01A)", "(01A)", "(01A)", "(10A)", "(10A)", "(10A)", "(11A)", "(11A)", "(11A)", "(00B)", "(00B)", "(00B)", "(01B)", "(01B)", "(01B)", "(10B)", "(10B)", "(10B)", "(11B)", "(11B)", "(11B)", "(00C)", "(00C)", "(00C)", "(01C)", "(01C)", "(01C)", "(10C)", "(10C)", "(10C)", "(11C)", "(11C)", "(11C)", "(00D)", "(00D)", "(00D)", "(01D)", "(01D)", "(01D)", "(10D)", "(10D)", "(10D)", "(11D)", "(11D)", "(11D)","(00E)", "(00E)", "(00E)", "(01E)", "(01E)", "(01E)", "(10E)", "(10E)", "(10E)", "(11E)", "(11E)", "(11E)")
        cols <- c("(01A)", "(10A)", "(11A)", "(00A)", "(10A)", "(11A)", "(00A)", "(01A)", "(11A)", "(00A)", "(01A)", "(10A)", "(01B)", "(10B)", "(11B)", "(00B)", "(10B)", "(11B)", "(00B)", "(01B)", "(11B)", "(00B)", "(01B)", "(10B)", "(01C)", "(10C)", "(11C)", "(00C)", "(10C)", "(11C)", "(00C)", "(01C)", "(11C)", "(00C)", "(01C)", "(10C)", "(01D)", "(10D)", "(11D)", "(00D)", "(10D)", "(11D)", "(00D)", "(01D)", "(11D)", "(00D)", "(01D)", "(10D)", "(01E)", "(10E)", "(11E)", "(00E)", "(10E)", "(11E)", "(00E)", "(01E)", "(11E)", "(00E)", "(01E)", "(10E)")
        trans.tmp <- trans.rate[cbind(rows,cols)]
        trans.tmp[which(trans.tmp > 0)] = (trans.tmp[which(trans.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, trans.tmp)
        rows <- c("(00A)","(00A)","(00A)","(00A)","(01A)","(01A)","(01A)","(01A)","(10A)","(10A)","(10A)","(10A)","(11A)","(11A)","(11A)","(11A)","(00B)","(00B)","(00B)","(00B)","(01B)","(01B)","(01B)","(01B)","(10B)","(10B)","(10B)","(10B)","(11B)","(11B)","(11B)","(11B)","(00C)","(00C)","(00C)","(00C)","(01C)","(01C)","(01C)","(01C)","(10C)","(10C)","(10C)","(10C)","(11C)","(11C)","(11C)","(11C)","(00D)","(00D)","(00D)","(00D)","(01D)","(01D)","(01D)","(01D)","(10D)","(10D)","(10D)","(10D)","(11D)","(11D)","(11D)","(11D)","(00E)","(00E)","(00E)","(00E)","(01E)","(01E)","(01E)","(01E)","(10E)","(10E)","(10E)","(10E)","(11E)","(11E)","(11E)","(11E)")
        cols <- c("(00B)","(00C)","(00D)","(00E)","(01B)","(01C)","(01D)","(01E)","(10B)","(10C)","(10D)","(10E)","(11B)","(11C)","(11D)","(11E)","(00A)","(00C)","(00D)","(00E)","(01A)","(01C)","(01D)","(01E)","(10A)","(10C)","(10D)","(10E)","(11A)","(11C)","(11D)","(11E)","(00A)","(00B)","(00D)","(00E)","(01A)","(01B)","(01D)","(01E)","(10A)","(10B)","(10D)","(10E)","(11A)","(11B)","(11D)","(11E)","(00A)","(00B)","(00C)","(00E)","(01A)","(01B)","(01C)","(01E)","(10A)","(10B)","(10C)","(10E)","(11A)","(11B)","(11C)","(11E)","(00A)","(00B)","(00C)","(00D)","(01A)","(01B)","(01C)","(01D)","(10A)","(10B)","(10C)","(10D)","(11A)","(11B)","(11C)","(11D)")
        category.tmp <- trans.rate[cbind(rows,cols)]
        category.tmp[category.tmp==0] <- NA
        category.rate.shift <- category.tmp + for.late.adjust
        category.rate.shift[is.na(category.rate.shift)] <- 0
        category.rate.shiftA <- c(category.rate.shift[1:4], rep(0,3), category.rate.shift[5:8], rep(0,3), category.rate.shift[9:12], rep(0,3), category.rate.shift[13:16], rep(0,3))
        category.rate.shiftB <- c(category.rate.shift[17:20], rep(0,3), category.rate.shift[21:24], rep(0,3), category.rate.shift[25:28], rep(0,3), category.rate.shift[29:32], rep(0,3))
        category.rate.shiftC <- c(category.rate.shift[33:36], rep(0,3), category.rate.shift[37:40], rep(0,3), category.rate.shift[41:44], rep(0,3), category.rate.shift[45:48], rep(0,3))
        category.rate.shiftD <- c(category.rate.shift[49:52], rep(0,3), category.rate.shift[53:56], rep(0,3), category.rate.shift[57:60], rep(0,3), category.rate.shift[61:64], rep(0,3))
        category.rate.shiftE <- c(category.rate.shift[65:68], rep(0,3), category.rate.shift[69:72], rep(0,3), category.rate.shift[73:76], rep(0,3), category.rate.shift[77:80], rep(0,3))
        pars.tmp <- c(turnover[1:4], eps.tmp[1:4], trans.tmp[1:12], category.rate.shiftA, turnover[5:8], eps.tmp[5:8], trans.tmp[13:24], category.rate.shiftB, turnover[9:12], eps.tmp[9:12], trans.tmp[25:36], category.rate.shiftC, turnover[13:16], eps.tmp[13:16], trans.tmp[37:48], category.rate.shiftD, turnover[17:20], eps.tmp[17:20], trans.tmp[49:60], category.rate.shiftE)
        pars[1:length(pars.tmp)] <- pars.tmp
    }
    
    if(dim(trans.rate)[2]==24){
        rate.cats <- 6
        pars.tmp <- turnover
        eps.tmp <- eps
        eps.tmp[which(eps.tmp > 0)] = (eps.tmp[which( eps.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, eps.tmp)
        for.late.adjust <- max(pars.tmp)
        rows <- c("(00A)", "(00A)", "(00A)", "(01A)", "(01A)", "(01A)", "(10A)", "(10A)", "(10A)", "(11A)", "(11A)", "(11A)", "(00B)", "(00B)", "(00B)", "(01B)", "(01B)", "(01B)", "(10B)", "(10B)", "(10B)", "(11B)", "(11B)", "(11B)", "(00C)", "(00C)", "(00C)", "(01C)", "(01C)", "(01C)", "(10C)", "(10C)", "(10C)", "(11C)", "(11C)", "(11C)", "(00D)", "(00D)", "(00D)", "(01D)", "(01D)", "(01D)", "(10D)", "(10D)", "(10D)", "(11D)", "(11D)", "(11D)", "(00E)", "(00E)", "(00E)", "(01E)", "(01E)", "(01E)", "(10E)", "(10E)", "(10E)", "(11E)", "(11E)", "(11E)", "(00F)", "(00F)", "(00F)", "(01F)", "(01F)", "(01F)", "(10F)", "(10F)", "(10F)", "(11F)", "(11F)", "(11F)")
        cols <- c("(01A)", "(10A)", "(11A)", "(00A)", "(10A)", "(11A)", "(00A)", "(01A)", "(11A)", "(00A)", "(01A)", "(10A)", "(01B)", "(10B)", "(11B)", "(00B)", "(10B)", "(11B)", "(00B)", "(01B)", "(11B)", "(00B)", "(01B)", "(10B)", "(01C)", "(10C)", "(11C)", "(00C)", "(10C)", "(11C)", "(00C)", "(01C)", "(11C)", "(00C)", "(01C)", "(10C)", "(01D)", "(10D)", "(11D)", "(00D)", "(10D)", "(11D)", "(00D)", "(01D)", "(11D)", "(00D)", "(01D)", "(10D)", "(01E)", "(10E)", "(11E)", "(00E)", "(10E)", "(11E)", "(00E)", "(01E)", "(11E)", "(00E)", "(01E)", "(10E)", "(01F)", "(10F)", "(11F)", "(00F)", "(10F)", "(11F)", "(00F)", "(01F)", "(11F)", "(00F)", "(01F)", "(10F)")
        trans.tmp <- trans.rate[cbind(rows,cols)]
        trans.tmp[which(trans.tmp > 0)] = (trans.tmp[which(trans.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, trans.tmp)
        rows <- c("(00A)","(00A)","(00A)","(00A)","(00A)","(01A)","(01A)","(01A)","(01A)","(01A)","(10A)","(10A)","(10A)","(10A)","(10A)","(11A)","(11A)","(11A)","(11A)","(11A)","(00B)","(00B)","(00B)","(00B)","(00B)","(01B)","(01B)","(01B)","(01B)","(01B)","(10B)","(10B)","(10B)","(10B)","(10B)","(11B)","(11B)","(11B)","(11B)","(11B)","(00C)","(00C)","(00C)","(00C)","(00C)","(01C)","(01C)","(01C)","(01C)","(01C)","(10C)","(10C)","(10C)","(10C)","(10C)","(11C)","(11C)","(11C)","(11C)","(11C)","(00D)","(00D)","(00D)","(00D)","(00D)","(01D)","(01D)","(01D)","(01D)","(01D)","(10D)","(10D)","(10D)","(10D)","(10D)","(11D)","(11D)","(11D)","(11D)","(11D)","(00E)","(00E)","(00E)","(00E)","(00E)","(01E)","(01E)","(01E)","(01E)","(01E)","(10E)","(10E)","(10E)","(10E)","(10E)","(11E)","(11E)","(11E)","(11E)","(11E)","(00F)","(00F)","(00F)","(00F)","(00F)","(01F)","(01F)","(01F)","(01F)","(01F)","(10F)","(10F)","(10F)","(10F)","(10F)","(11F)","(11F)","(11F)","(11F)","(11F)")
        cols <- c("(00B)","(00C)","(00D)","(00E)","(00F)","(01B)","(01C)","(01D)","(01E)","(01F)","(10B)","(10C)","(10D)","(10E)","(10F)","(11B)","(11C)","(11D)","(11E)","(11F)","(00A)","(00C)","(00D)","(00E)","(00F)","(01A)","(01C)","(01D)","(01E)","(01F)","(10A)","(10C)","(10D)","(10E)","(10F)","(11A)","(11C)","(11D)","(11E)","(11F)","(00A)","(00B)","(00D)","(00E)","(00F)","(01A)","(01B)","(01D)","(01E)","(01F)","(10A)","(10B)","(10D)","(10E)","(10F)","(11A)","(11B)","(11D)","(11E)","(11F)","(00A)","(00B)","(00C)","(00E)","(00F)","(01A)","(01B)","(01C)","(01E)","(01F)","(10A)","(10B)","(10C)","(10E)","(10F)","(11A)","(11B)","(11C)","(11E)","(11F)","(00A)","(00B)","(00C)","(00D)","(00F)","(01A)","(01B)","(01C)","(01D)","(01F)","(10A)","(10B)","(10C)","(10D)","(10F)","(11A)","(11B)","(11C)","(11D)","(11F)","(00A)","(00B)","(00C)","(00D)","(00E)","(01A)","(01B)","(01C)","(01D)","(01E)","(10A)","(10B)","(10C)","(10D)","(10E)","(11A)","(11B)","(11C)","(11D)","(11E)")
        category.tmp <- trans.rate[cbind(rows,cols)]
        category.tmp[category.tmp==0] <- NA
        category.rate.shift <- category.tmp + for.late.adjust
        category.rate.shift[is.na(category.rate.shift)] <- 0
        category.rate.shiftA <- c(category.rate.shift[1:5], rep(0,2), category.rate.shift[6:10], rep(0,2), category.rate.shift[11:15], rep(0,2), category.rate.shift[16:20], rep(0,2))
        category.rate.shiftB <- c(category.rate.shift[21:25], rep(0,2), category.rate.shift[26:30], rep(0,2), category.rate.shift[31:35], rep(0,2), category.rate.shift[36:40], rep(0,2))
        category.rate.shiftC <- c(category.rate.shift[41:45], rep(0,2), category.rate.shift[46:50], rep(0,2), category.rate.shift[51:55], rep(0,2), category.rate.shift[56:60], rep(0,2))
        category.rate.shiftD <- c(category.rate.shift[61:65], rep(0,2), category.rate.shift[66:70], rep(0,2), category.rate.shift[71:75], rep(0,2), category.rate.shift[76:80], rep(0,2))
        category.rate.shiftE <- c(category.rate.shift[81:85], rep(0,2), category.rate.shift[86:90], rep(0,2), category.rate.shift[91:95], rep(0,2), category.rate.shift[96:100], rep(0,2))
        category.rate.shiftF <- c(category.rate.shift[101:105], rep(0,2), category.rate.shift[106:110], rep(0,2), category.rate.shift[111:115], rep(0,2), category.rate.shift[116:120], rep(0,2))
        pars.tmp <- c(turnover[1:4], eps.tmp[1:4], trans.tmp[1:12], category.rate.shiftA, turnover[5:8], eps.tmp[5:8], trans.tmp[13:24], category.rate.shiftB, turnover[9:12], eps.tmp[9:12], trans.tmp[25:36], category.rate.shiftC, turnover[13:16], eps.tmp[13:16], trans.tmp[37:48], category.rate.shiftD, turnover[17:20], eps.tmp[17:20], trans.tmp[49:60], category.rate.shiftE, turnover[21:24], eps.tmp[21:24], trans.tmp[61:72], category.rate.shiftF)
        pars[1:length(pars.tmp)] <- pars.tmp
    }
    
    if(dim(trans.rate)[2]==28){
        rate.cats <- 7
        pars.tmp <- turnover
        eps.tmp <- eps
        eps.tmp[which(eps.tmp > 0)] = (eps.tmp[which( eps.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, eps.tmp)
        for.late.adjust <- max(pars.tmp)
        rows <- c("(00A)", "(00A)", "(00A)", "(01A)", "(01A)", "(01A)", "(10A)", "(10A)", "(10A)", "(11A)", "(11A)", "(11A)", "(00B)", "(00B)", "(00B)", "(01B)", "(01B)", "(01B)", "(10B)", "(10B)", "(10B)", "(11B)", "(11B)", "(11B)", "(00C)", "(00C)", "(00C)", "(01C)", "(01C)", "(01C)", "(10C)", "(10C)", "(10C)", "(11C)", "(11C)", "(11C)", "(00D)", "(00D)", "(00D)", "(01D)", "(01D)", "(01D)", "(10D)", "(10D)", "(10D)", "(11D)", "(11D)", "(11D)", "(00E)", "(00E)", "(00E)", "(01E)", "(01E)", "(01E)", "(10E)", "(10E)", "(10E)", "(11E)", "(11E)", "(11E)", "(00F)", "(00F)", "(00F)", "(01F)", "(01F)", "(01F)", "(10F)", "(10F)", "(10F)", "(11F)", "(11F)", "(11F)", "(00G)", "(00G)", "(00G)", "(01G)", "(01G)", "(01G)", "(10G)", "(10G)", "(10G)", "(11G)", "(11G)", "(11G)")
        cols <- c("(01A)", "(10A)", "(11A)", "(00A)", "(10A)", "(11A)", "(00A)", "(01A)", "(11A)", "(00A)", "(01A)", "(10A)", "(01B)", "(10B)", "(11B)", "(00B)", "(10B)", "(11B)", "(00B)", "(01B)", "(11B)", "(00B)", "(01B)", "(10B)", "(01C)", "(10C)", "(11C)", "(00C)", "(10C)", "(11C)", "(00C)", "(01C)", "(11C)", "(00C)", "(01C)", "(10C)", "(01D)", "(10D)", "(11D)", "(00D)", "(10D)", "(11D)", "(00D)", "(01D)", "(11D)", "(00D)", "(01D)", "(10D)", "(01E)", "(10E)", "(11E)", "(00E)", "(10E)", "(11E)", "(00E)", "(01E)", "(11E)", "(00E)", "(01E)", "(10E)", "(01F)", "(10F)", "(11F)", "(00F)", "(10F)", "(11F)", "(00F)", "(01F)", "(11F)", "(00F)", "(01F)", "(10F)", "(01G)", "(10G)", "(11G)", "(00G)", "(10G)", "(11G)", "(00G)", "(01G)", "(11G)", "(00G)", "(01G)", "(10G)")
        trans.tmp <- trans.rate[cbind(rows,cols)]
        trans.tmp[which(trans.tmp > 0)] = (trans.tmp[which(trans.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, trans.tmp)
        rows <- c("(00A)","(00A)","(00A)","(00A)","(00A)","(00A)","(01A)","(01A)","(01A)","(01A)","(01A)","(01A)","(10A)","(10A)","(10A)","(10A)","(10A)","(10A)","(11A)","(11A)","(11A)","(11A)","(11A)","(11A)","(00B)","(00B)","(00B)","(00B)","(00B)","(00B)","(01B)","(01B)","(01B)","(01B)","(01B)","(01B)","(10B)","(10B)","(10B)","(10B)","(10B)","(10B)","(11B)","(11B)","(11B)","(11B)","(11B)","(11B)","(00C)","(00C)","(00C)","(00C)","(00C)","(00C)","(01C)","(01C)","(01C)","(01C)","(01C)","(01C)","(10C)","(10C)","(10C)","(10C)","(10C)","(10C)","(11C)","(11C)","(11C)","(11C)","(11C)","(11C)","(00D)","(00D)","(00D)","(00D)","(00D)","(00D)","(01D)","(01D)","(01D)","(01D)","(01D)","(01D)","(10D)","(10D)","(10D)","(10D)","(10D)","(10D)","(11D)","(11D)","(11D)","(11D)","(11D)","(11D)","(00E)","(00E)","(00E)","(00E)","(00E)","(00E)","(01E)","(01E)","(01E)","(01E)","(01E)","(01E)","(10E)","(10E)","(10E)","(10E)","(10E)","(10E)","(11E)","(11E)","(11E)","(11E)","(11E)","(11E)","(00F)","(00F)","(00F)","(00F)","(00F)","(00F)","(01F)","(01F)","(01F)","(01F)","(01F)","(01F)","(10F)","(10F)","(10F)","(10F)","(10F)","(10F)","(11F)","(11F)","(11F)","(11F)","(11F)","(11F)","(00G)","(00G)","(00G)","(00G)","(00G)","(00G)","(01G)","(01G)","(01G)","(01G)","(01G)","(01G)","(10G)","(10G)","(10G)","(10G)","(10G)","(10G)","(11G)","(11G)","(11G)","(11G)","(11G)","(11G)")
        cols <- c("(00B)","(00C)","(00D)","(00E)","(00F)","(00G)","(01B)","(01C)","(01D)","(01E)","(01F)","(01G)","(10B)","(10C)","(10D)","(10E)","(10F)","(10G)","(11B)","(11C)","(11D)","(11E)","(11F)","(11G)","(00A)","(00C)","(00D)","(00E)","(00F)","(00G)","(01A)","(01C)","(01D)","(01E)","(01F)","(01G)","(10A)","(10C)","(10D)","(10E)","(10F)","(10G)","(11A)","(11C)","(11D)","(11E)","(11F)","(11G)","(00A)","(00B)","(00D)","(00E)","(00F)","(00G)","(01A)","(01B)","(01D)","(01E)","(01F)","(01G)","(10A)","(10B)","(10D)","(10E)","(10F)","(10G)","(11A)","(11B)","(11D)","(11E)","(11F)","(11G)","(00A)","(00B)","(00C)","(00E)","(00F)","(00G)","(01A)","(01B)","(01C)","(01E)","(01F)","(01G)","(10A)","(10B)","(10C)","(10E)","(10F)","(10G)","(11A)","(11B)","(11C)","(11E)","(11F)","(11G)","(00A)","(00B)","(00C)","(00D)","(00F)","(00G)","(01A)","(01B)","(01C)","(01D)","(01F)","(01G)","(10A)","(10B)","(10C)","(10D)","(10F)","(10G)","(11A)","(11B)","(11C)","(11D)","(11F)","(11G)","(00A)","(00B)","(00C)","(00D)","(00E)","(00G)","(01A)","(01B)","(01C)","(01D)","(01E)","(01G)","(10A)","(10B)","(10C)","(10D)","(10E)","(10G)","(11A)","(11B)","(11C)","(11D)","(11E)","(11G)","(00A)","(00B)","(00C)","(00D)","(00E)","(00F)","(01A)","(01B)","(01C)","(01D)","(01E)","(01F)","(10A)","(10B)","(10C)","(10D)","(10E)","(10F)","(11A)","(11B)","(11C)","(11D)","(11E)","(11F)")
        category.tmp <- trans.rate[cbind(rows,cols)]
        category.tmp[category.tmp==0] <- NA
        category.rate.shift <- category.tmp + for.late.adjust
        category.rate.shift[is.na(category.rate.shift)] <- 0
        category.rate.shiftA <- c(category.rate.shift[1:6], rep(0,1), category.rate.shift[7:12], rep(0,1), category.rate.shift[13:18], rep(0,1), category.rate.shift[19:24], rep(0,1))
        category.rate.shiftB <- c(category.rate.shift[25:30], rep(0,1), category.rate.shift[31:36], rep(0,1), category.rate.shift[37:42], rep(0,1), category.rate.shift[43:48], rep(0,1))
        category.rate.shiftC <- c(category.rate.shift[49:54], rep(0,1), category.rate.shift[55:60], rep(0,1), category.rate.shift[61:66], rep(0,1), category.rate.shift[67:72], rep(0,1))
        category.rate.shiftD <- c(category.rate.shift[73:78], rep(0,1), category.rate.shift[79:84], rep(0,1), category.rate.shift[85:90], rep(0,1), category.rate.shift[91:96], rep(0,1))
        category.rate.shiftE <- c(category.rate.shift[97:102], rep(0,1), category.rate.shift[103:108], rep(0,1), category.rate.shift[109:114], rep(0,1), category.rate.shift[115:120], rep(0,1))
        category.rate.shiftF <- c(category.rate.shift[121:126], rep(0,1), category.rate.shift[127:132], rep(0,1), category.rate.shift[133:138], rep(0,1), category.rate.shift[139:144], rep(0,1))
        category.rate.shiftG <- c(category.rate.shift[145:150], rep(0,1), category.rate.shift[151:156], rep(0,1), category.rate.shift[157:162], rep(0,1), category.rate.shift[163:168], rep(0,1))
        pars.tmp <- c(turnover[1:4], eps.tmp[1:4], trans.tmp[1:12], category.rate.shiftA, turnover[5:8], eps.tmp[5:8], trans.tmp[13:24], category.rate.shiftB, turnover[9:12], eps.tmp[9:12], trans.tmp[25:36], category.rate.shiftC, turnover[13:16], eps.tmp[13:16], trans.tmp[37:48], category.rate.shiftD, turnover[17:20], eps.tmp[17:20], trans.tmp[49:60], category.rate.shiftE, turnover[21:24], eps.tmp[21:24], trans.tmp[61:72], category.rate.shiftF, turnover[25:28], eps.tmp[25:28], trans.tmp[73:84], category.rate.shiftG)
        pars[1:length(pars.tmp)] <- pars.tmp
    }
    
    if(dim(trans.rate)[2]==32){
        rate.cats <- 8
        pars.tmp <- turnover
        eps.tmp <- eps
        eps.tmp[which(eps.tmp > 0)] = (eps.tmp[which( eps.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, eps.tmp)
        for.late.adjust <- max(pars.tmp)
        rows <- c("(00A)", "(00A)", "(00A)", "(01A)", "(01A)", "(01A)", "(10A)", "(10A)", "(10A)", "(11A)", "(11A)", "(11A)", "(00B)", "(00B)", "(00B)", "(01B)", "(01B)", "(01B)", "(10B)", "(10B)", "(10B)", "(11B)", "(11B)", "(11B)", "(00C)", "(00C)", "(00C)", "(01C)", "(01C)", "(01C)", "(10C)", "(10C)", "(10C)", "(11C)", "(11C)", "(11C)", "(00D)", "(00D)", "(00D)", "(01D)", "(01D)", "(01D)", "(10D)", "(10D)", "(10D)", "(11D)", "(11D)", "(11D)", "(00E)", "(00E)", "(00E)", "(01E)", "(01E)", "(01E)", "(10E)", "(10E)", "(10E)", "(11E)", "(11E)", "(11E)", "(00F)", "(00F)", "(00F)", "(01F)", "(01F)", "(01F)", "(10F)", "(10F)", "(10F)", "(11F)", "(11F)", "(11F)", "(00G)", "(00G)", "(00G)", "(01G)", "(01G)", "(01G)", "(10G)", "(10G)", "(10G)", "(11G)", "(11G)", "(11G)", "(00H)", "(00H)", "(00H)", "(01H)", "(01H)", "(01H)", "(10H)", "(10H)", "(10H)", "(11H)", "(11H)", "(11H)")
        cols <- c("(01A)", "(10A)", "(11A)", "(00A)", "(10A)", "(11A)", "(00A)", "(01A)", "(11A)", "(00A)", "(01A)", "(10A)", "(01B)", "(10B)", "(11B)", "(00B)", "(10B)", "(11B)", "(00B)", "(01B)", "(11B)", "(00B)", "(01B)", "(10B)", "(01C)", "(10C)", "(11C)", "(00C)", "(10C)", "(11C)", "(00C)", "(01C)", "(11C)", "(00C)", "(01C)", "(10C)", "(01D)", "(10D)", "(11D)", "(00D)", "(10D)", "(11D)", "(00D)", "(01D)", "(11D)", "(00D)", "(01D)", "(10D)", "(01E)", "(10E)", "(11E)", "(00E)", "(10E)", "(11E)", "(00E)", "(01E)", "(11E)", "(00E)", "(01E)", "(10E)", "(01F)", "(10F)", "(11F)", "(00F)", "(10F)", "(11F)", "(00F)", "(01F)", "(11F)", "(00F)", "(01F)", "(10F)", "(01G)", "(10G)", "(11G)", "(00G)", "(10G)", "(11G)", "(00G)", "(01G)", "(11G)", "(00G)", "(01G)", "(10G)", "(01H)", "(10H)", "(11H)", "(00H)", "(10H)", "(11H)", "(00H)", "(01H)", "(11H)", "(00H)", "(01H)", "(10H)")
        trans.tmp <- trans.rate[cbind(rows,cols)]
        trans.tmp[which(trans.tmp > 0)] = (trans.tmp[which(trans.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, trans.tmp)
        rows <- c("(00A)","(00A)","(00A)","(00A)","(00A)","(00A)","(00A)","(01A)","(01A)","(01A)","(01A)","(01A)","(01A)","(01A)","(10A)","(10A)","(10A)","(10A)","(10A)","(10A)","(10A)","(11A)","(11A)","(11A)","(11A)","(11A)","(11A)","(11A)","(00B)","(00B)","(00B)","(00B)","(00B)","(00B)","(00B)","(01B)","(01B)","(01B)","(01B)","(01B)","(01B)","(01B)","(10B)","(10B)","(10B)","(10B)","(10B)","(10B)","(10B)","(11B)","(11B)","(11B)","(11B)","(11B)","(11B)","(11B)","(00C)","(00C)","(00C)","(00C)","(00C)","(00C)","(00C)","(01C)","(01C)","(01C)","(01C)","(01C)","(01C)","(01C)","(10C)","(10C)","(10C)","(10C)","(10C)","(10C)","(10C)","(11C)","(11C)","(11C)","(11C)","(11C)","(11C)","(11C)","(00D)","(00D)","(00D)","(00D)","(00D)","(00D)","(00D)","(01D)","(01D)","(01D)","(01D)","(01D)","(01D)","(01D)","(10D)","(10D)","(10D)","(10D)","(10D)","(10D)","(10D)","(11D)","(11D)","(11D)","(11D)","(11D)","(11D)","(11D)","(00E)","(00E)","(00E)","(00E)","(00E)","(00E)","(00E)","(01E)","(01E)","(01E)","(01E)","(01E)","(01E)","(01E)","(10E)","(10E)","(10E)","(10E)","(10E)","(10E)","(10E)","(11E)","(11E)","(11E)","(11E)","(11E)","(11E)","(11E)","(00F)","(00F)","(00F)","(00F)","(00F)","(00F)","(00F)","(01F)","(01F)","(01F)","(01F)","(01F)","(01F)","(01F)","(10F)","(10F)","(10F)","(10F)","(10F)","(10F)","(10F)","(11F)","(11F)","(11F)","(11F)","(11F)","(11F)","(11F)","(00G)","(00G)","(00G)","(00G)","(00G)","(00G)","(00G)","(01G)","(01G)","(01G)","(01G)","(01G)","(01G)","(01G)","(10G)","(10G)","(10G)","(10G)","(10G)","(10G)","(10G)","(11G)","(11G)","(11G)","(11G)","(11G)","(11G)","(11G)","(00H)","(00H)","(00H)","(00H)","(00H)","(00H)","(00H)","(01H)","(01H)","(01H)","(01H)","(01H)","(01H)","(01H)","(10H)","(10H)","(10H)","(10H)","(10H)","(10H)","(10H)","(11H)","(11H)","(11H)","(11H)","(11H)","(11H)","(11H)")
        cols <- c("(00B)","(00C)","(00D)","(00E)","(00F)","(00G)","(00H)","(01B)","(01C)","(01D)","(01E)","(01F)","(01G)","(01H)","(10B)","(10C)","(10D)","(10E)","(10F)","(10G)","(10H)","(11B)","(11C)","(11D)","(11E)","(11F)","(11G)","(11H)","(00A)","(00C)","(00D)","(00E)","(00F)","(00G)","(00H)","(01A)","(01C)","(01D)","(01E)","(01F)","(01G)","(01H)","(10A)","(10C)","(10D)","(10E)","(10F)","(10G)","(10H)","(11A)","(11C)","(11D)","(11E)","(11F)","(11G)","(11H)","(00A)","(00B)","(00D)","(00E)","(00F)","(00G)","(00H)","(01A)","(01B)","(01D)","(01E)","(01F)","(01G)","(01H)","(10A)","(10B)","(10D)","(10E)","(10F)","(10G)","(10H)","(11A)","(11B)","(11D)","(11E)","(11F)","(11G)","(11H)","(00A)","(00B)","(00C)","(00E)","(00F)","(00G)","(00H)","(01A)","(01B)","(01C)","(01E)","(01F)","(01G)","(01H)","(10A)","(10B)","(10C)","(10E)","(10F)","(10G)","(10H)","(11A)","(11B)","(11C)","(11E)","(11F)","(11G)","(11H)","(00A)","(00B)","(00C)","(00D)","(00F)","(00G)","(00H)","(01A)","(01B)","(01C)","(01D)","(01F)","(01G)","(01H)","(10A)","(10B)","(10C)","(10D)","(10F)","(10G)","(10H)","(11A)","(11B)","(11C)","(11D)","(11F)","(11G)","(11H)","(00A)","(00B)","(00C)","(00D)","(00E)","(00G)","(00H)","(01A)","(01B)","(01C)","(01D)","(01E)","(01G)","(01H)","(10A)","(10B)","(10C)","(10D)","(10E)","(10G)","(10H)","(11A)","(11B)","(11C)","(11D)","(11E)","(11G)","(11H)","(00A)","(00B)","(00C)","(00D)","(00E)","(00F)","(00H)","(01A)","(01B)","(01C)","(01D)","(01E)","(01F)","(01H)","(10A)","(10B)","(10C)","(10D)","(10E)","(10F)","(10H)","(11A)","(11B)","(11C)","(11D)","(11E)","(11F)","(11H)","(00A)","(00B)","(00C)","(00D)","(00E)","(00F)","(00G)","(01A)","(01B)","(01C)","(01D)","(01E)","(01F)","(01G)","(10A)","(10B)","(10C)","(10D)","(10E)","(10F)","(10G)","(11A)","(11B)","(11C)","(11D)","(11E)","(11F)","(11G)")
        category.tmp <- trans.rate[cbind(rows,cols)]
        category.tmp[category.tmp==0] <- NA
        category.rate.shift <- category.tmp + for.late.adjust
        category.rate.shift[is.na(category.rate.shift)] <- 0
        category.rate.shiftA <- c(category.rate.shift[1:28])
        category.rate.shiftB <- c(category.rate.shift[29:56])
        category.rate.shiftC <- c(category.rate.shift[57:84])
        category.rate.shiftD <- c(category.rate.shift[85:112])
        category.rate.shiftE <- c(category.rate.shift[113:140])
        category.rate.shiftF <- c(category.rate.shift[141:168])
        category.rate.shiftG <- c(category.rate.shift[169:196])
        category.rate.shiftH <- c(category.rate.shift[197:224])
        pars.tmp <- c(turnover[1:4], eps.tmp[1:4], trans.tmp[1:12], category.rate.shiftA, turnover[5:8], eps.tmp[5:8], trans.tmp[13:24], category.rate.shiftB, turnover[9:12], eps.tmp[9:12], trans.tmp[25:36], category.rate.shiftC, turnover[13:16], eps.tmp[13:16], trans.tmp[37:48], category.rate.shiftD, turnover[17:20], eps.tmp[17:20], trans.tmp[49:60], category.rate.shiftE, turnover[21:24], eps.tmp[21:24], trans.tmp[61:72], category.rate.shiftF, turnover[25:28], eps.tmp[25:28], trans.tmp[73:84], category.rate.shiftG, turnover[29:32], eps.tmp[29:32], trans.tmp[85:96], category.rate.shiftH)
        pars[1:length(pars.tmp)] <- pars.tmp
    }
    
    np <- max(pars)
    pars[pars==0] <- np+1
    
    cat("Initializing...", "\n")
    
    data.new <- data.frame(data[,2], data[,3], row.names=data[,1])
    data.new <- data.new[phy$tip.label,]
    
    #This is used to scale starting values to account for sampling:
    if(length(f) == 4){
        freqs <- table(apply(data.new, 1, function(x) switch(paste0(x, collapse=""), "00" = 1, "01" = 2, "10" = 3, "11" = 4, "02"=1, "20"=3, "21"=2, "12"=4, "22"=4)))
        
        ## if(length(freqs == 4)){
        ##     freqs[which(!c(1:4) %in% names(freqs))] <- 0
        ##     samp.freq.tree <- Ntip(phy) / sum(freqs / f)
        ## }else{
        ##     samp.freq.tree <- Ntip(phy) / sum(freqs / f)
        ## }
        
        ## Fixing the structure of the freqs vector.
        freqs.vec <- rep(0, times = 4)
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
            init.pars <- starting.point.generator(phy, 4, samp.freq.tree, yule=TRUE)
        }else{
            init.pars <- starting.point.generator(phy, 4, samp.freq.tree, yule=FALSE)
            if(any(init.pars[5:8] == 0)){
                init.pars[5:8] = 1e-6
            }
        }
        names(init.pars) <- NULL
        
        if(is.null(starting.vals)){
            def.set.pars <- rep(c(log(init.pars[1:4]+init.pars[5:8]), log(init.pars[5:8]/init.pars[1:4]), log(init.pars[9:20]), rep(log(.01), 28)), rate.cats)
        }else{
            ## Check if 'starting.vals' has the correct format.
            if( !length(starting.vals) %in% c(3,20) ){
                stop("Wrong length of starting.vals")
            }
            if( length(starting.vals) == 20 ){
                cat("Using developer mode for starting.vals.", "\n")
                def.set.pars <- rep(c(log(starting.vals[1:4]), log(starting.vals[5:8]), log(starting.vals[9:20]), rep(log(0.01), 28)), rate.cats)
            } else{
                def.set.pars <- rep(c(log( rep(starting.vals[1],4) ), log( rep(starting.vals[2],4) ), log( rep(starting.vals[3],12) ), rep(log(0.01), 28)), rate.cats)
            }
        }
        
        if(bounded.search == TRUE){
            upper.full <- rep(c(rep(log(turnover.upper),4), rep(log(eps.upper),4), rep(log(trans.upper),12), rep(log(10), 28)), rate.cats)
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
    dat.tab <- OrganizeData(data=data.new, phy=phy, f=f, hidden.states=hidden.states)
    nb.tip <- Ntip(phy)
    nb.node <- phy$Nnode
    ##########################
    
    if(sann == FALSE){
        if(bounded.search == TRUE){
            cat("Finished. Beginning bounded subplex routine...", "\n")
            opts <- list("algorithm" = "NLOPT_LN_SBPLX", "maxeval" = 100000, "ftol_rel" = max.tol)
            out = nloptr(x0=ip, eval_f=DevOptimizeMuHiSSE, ub=upper, lb=lower, opts=opts, pars=pars, dat.tab=dat.tab, gen=gen, hidden.states=hidden.states, nb.tip=nb.tip, nb.node=nb.node, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, np=np, ode.eps=ode.eps)
            solution <- numeric(length(pars))
            solution[] <- c(exp(out$solution), 0)[pars]
            loglik = -out$objective
        }else{
            cat("Finished. Beginning subplex routine...", "\n")
            out = subplex(ip, fn=DevOptimizeMuHiSSE, control=list(reltol=max.tol, parscale=rep(0.1, length(ip))), pars=pars, dat.tab=dat.tab, gen=gen, hidden.states=hidden.states, nb.tip=nb.tip, nb.node=nb.node, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, np=np, ode.eps=ode.eps)
            solution <- numeric(length(pars))
            solution[] <- c(exp(out$par), 0)[pars]
            loglik = -out$value
        }
    }else{
        cat("Finished. Beginning simulated annealing...", "\n")
        opts <- list("algorithm"="NLOPT_GD_STOGO", "maxeval" = 100000, "ftol_rel" = max.tol)
        out.sann = GenSA(ip, fn=DevOptimizeMuHiSSE, lower=lower, upper=upper, control=list(max.call=sann.its), pars=pars, dat.tab=dat.tab, gen=gen, hidden.states=hidden.states, nb.tip=nb.tip, nb.node=nb.node, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, np=np, ode.eps=ode.eps)
        #out.sann <- stogo(x0=ip, fn=DevOptimizeMuHiSSE, gr=NULL, upper=upper, lower=lower, pars=pars, phy=phy, data=data.new, f=f, hidden.states=hidden.states, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, np=np, ode.eps=ode.eps)
        cat("Finished. Refining using subplex routine...", "\n")
        opts <- list("algorithm" = "NLOPT_LN_SBPLX", "maxeval" = 100000, "ftol_rel" = max.tol)
        out <- nloptr(x0=out.sann$par, eval_f=DevOptimizeMuHiSSE, ub=upper, lb=lower, opts=opts, pars=pars, dat.tab=dat.tab, gen=gen, hidden.states=hidden.states, nb.tip=nb.tip, nb.node=nb.node, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, np=np, ode.eps=ode.eps)
        solution <- numeric(length(pars))
        solution[] <- c(exp(out$solution), 0)[pars]
        
        loglik = -out$objective
    }
    
    names(solution) <- c("turnover00A","turnover01A","turnover10A","turnover11A","eps00A","eps01A","eps10A","eps11A","q00A_01A","q00A_10A","q00A_11A","q01A_00A","q01A_10A","q01A_11A","q10A_00A","q10A_01A","q10A_11A","q11A_00A","q11A_01A","q11A_10A","q00A_00B","q00A_00C","q00A_00D","q00A_00E","q00A_00F","q00A_00G","q00A_00H","q01A_01B","q01A_01C","q01A_01D","q01A_01E","q01A_01F","q01A_01G","q01A_01H","q10A_10B","q10A_10C","q10A_10D","q10A_10E","q10A_10F","q10A_10G","q10A_10H","q11A_11B","q11A_11C","q11A_11D","q11A_11E","q11A_11F","q11A_11G","q11A_11H","turnover00B","turnover01B","turnover10B","turnover11B","eps00B","eps01B","eps10B","eps11B","q00B_01B","q00B_10B","q00B_11B","q01B_00B","q01B_10B","q01B_11B","q10B_00B","q10B_01B","q10B_11B","q11B_00B","q11B_01B","q11B_10B","q00B_00A","q00B_00C","q00B_00D","q00B_00E","q00B_00F","q00B_00G","q00B_00H","q01B_01A","q01B_01C","q01B_01D","q01B_01E","q01B_01F","q01B_01G","q01B_01H","q10B_10A","q10B_10C","q10B_10D","q10B_10E","q10B_10F","q10B_10G","q10B_10H","q11B_11A","q11B_11C","q11B_11D","q11B_11E","q11B_11F","q11B_11G","q11B_11H","turnover00C","turnover01C","turnover10C","turnover11C","eps00C","eps01C","eps10C","eps11C","q00C_01C","q00C_10C","q00C_11C","q01C_00C","q01C_10C","q01C_11C","q10C_00C","q10C_01C","q10C_11C","q11C_00C","q11C_01C","q11C_10C","q00C_00A","q00C_00B","q00C_00D","q00C_00E","q00C_00F","q00C_00G","q00C_00H","q01C_01A","q01C_01B","q01C_01D","q01C_01E","q01C_01F","q01C_01G","q01C_01H","q10C_10A","q10C_10B","q10C_10D","q10C_10E","q10C_10F","q10C_10G","q10C_10H","q11C_11A","q11C_11B","q11C_11D","q11C_11E","q11C_11F","q11C_11G","q11C_11H","turnover00D","turnover01D","turnover10D","turnover11D","eps00D","eps01D","eps10D","eps11D","q00D_01D","q00D_10D","q00D_11D","q01D_00D","q01D_10D","q01D_11D","q10D_00D","q10D_01D","q10D_11D","q11D_00D","q11D_01D","q11D_10D","q00D_00A","q00D_00B","q00D_00C","q00D_00E","q00D_00F","q00D_00G","q00D_00H","q01D_01A","q01D_01B","q01D_01C","q01D_01E","q01D_01F","q01D_01G","q01D_01H","q10D_10A","q10D_10B","q10D_10C","q10D_10E","q10D_10F","q10D_10G","q10D_10H","q11D_11A","q11D_11B","q11D_11C","q11D_11E","q11D_11F","q11D_11G","q11D_11H","turnover00E","turnover01E","turnover10E","turnover11E","eps00E","eps01E","eps10E","eps11E","q00E_01E","q00E_10E","q00E_11E","q01E_00E","q01E_10E","q01E_11E","q10E_00E","q10E_01E","q10E_11E","q11E_00E","q11E_01E","q11E_10E","q00E_00A","q00E_00B","q00E_00C","q00E_00D","q00E_00F","q00E_00G","q00E_00H","q01E_01A","q01E_01B","q01E_01C","q01E_01D","q01E_01F","q01E_01G","q01E_01H","q10E_10A","q10E_10B","q10E_10C","q10E_10D","q10E_10F","q10E_10G","q10E_10H","q11E_11A","q11E_11B","q11E_11C","q11E_11D","q11E_11F","q11E_11G","q11E_11H","turnover00F","turnover01F","turnover10F","turnover11F","eps00F","eps01F","eps10F","eps11F","q00F_01F","q00F_10F","q00F_11F","q01F_00F","q01F_10F","q01F_11F","q10F_00F","q10F_01F","q10F_11F","q11F_00F","q11F_01F","q11F_10F","q00F_00A","q00F_00B","q00F_00C","q00F_00D","q00F_00E","q00F_00G","q00F_00H","q01F_01A","q01F_01B","q01F_01C","q01F_01D","q01F_01E","q01F_01G","q01F_01H","q10F_10A","q10F_10B","q10F_10C","q10F_10D","q10F_10E","q10F_10G","q10F_10H","q11F_11A","q11F_11B","q11F_11C","q11F_11D","q11F_11E","q11F_11G","q11F_11H","turnover00G","turnover01G","turnover10G","turnover11G","eps00G","eps01G","eps10G","eps11G","q00G_01G","q00G_10G","q00G_11G","q01G_00G","q01G_10G","q01G_11G","q10G_00G","q10G_01G","q10G_11G","q11G_00G","q11G_01G","q11G_10G","q00G_00A","q00G_00B","q00G_00C","q00G_00D","q00G_00E","q00G_00F","q00G_00H","q01G_01A","q01G_01B","q01G_01C","q01G_01D","q01G_01E","q01G_01F","q01G_01H","q10G_10A","q10G_10B","q10G_10C","q10G_10D","q10G_10E","q10G_10F","q10G_10H","q11G_11A","q11G_11B","q11G_11C","q11G_11D","q11G_11E","q11G_11F","q11G_11H","turnover00H","turnover01H","turnover10H","turnover11H","eps00H","eps01H","eps10H","eps11H","q00H_01H","q00H_10H","q00H_11H","q01H_00H","q01H_10H","q01H_11H","q10H_00H","q10H_01H","q10H_11H","q11H_00H","q11H_01H","q11H_10H","q00H_00A","q00H_00B","q00H_00C","q00H_00D","q00H_00E","q00H_00F","q00H_00G","q01H_01A","q01H_01B","q01H_01C","q01H_01D","q01H_01E","q01H_01F","q01H_01G","q10H_10A","q10H_10B","q10H_10C","q10H_10D","q10H_10E","q10H_10F","q10H_10G","q11H_11A","q11H_11B","q11H_11C","q11H_11D","q11H_11E","q11H_11F","q11H_11G")
    
    cat("Finished. Summarizing results...", "\n")
    
    obj = list(loglik = loglik, AIC = -2*loglik+2*np, AICc = -2*loglik+(2*np*(Ntip(phy)/(Ntip(phy)-np-1))), solution=solution, index.par=pars, f=f, hidden.states=hidden.states, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, phy=phy, data=data, trans.matrix=trans.rate, max.tol=max.tol, starting.vals=ip, upper.bounds=upper, lower.bounds=lower, ode.eps=ode.eps)
    class(obj) <- append(class(obj), "muhisse.fit")
    return(obj)
}

######################################################################################################################################
######################################################################################################################################
### The function used to optimize parameters:
######################################################################################################################################
######################################################################################################################################


DevOptimizeMuHiSSE <- function(p, pars, dat.tab, gen, hidden.states, nb.tip=nb.tip, nb.node=nb.node, condition.on.survival, root.type, root.p, np, ode.eps) {
    #Generates the final vector with the appropriate parameter estimates in the right place:
    p.new <- exp(p)
    ## print(p.new)
    model.vec <- numeric(length(pars))
    model.vec[] <- c(p.new, 0)[pars]
    cache = ParametersToPassMuHiSSE(model.vec=model.vec, hidden.states=hidden.states, nb.tip=nb.tip, nb.node=nb.node, bad.likelihood=exp(-500), ode.eps=ode.eps)
    logl <- DownPassMuHisse(dat.tab=dat.tab, cache=cache, gen=gen, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p)

    return(-logl)
}



######################################################################################################################################
######################################################################################################################################
### The various utility functions used
######################################################################################################################################
######################################################################################################################################

OrganizeData <- function(data, phy, f, hidden.states){
    ### Ughy McUgherson. This is a must in order to pass CRAN checks: http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
    DesNode = NULL

    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
    
    if(hidden.states == FALSE){
        states = matrix(0,Ntip(phy),4)
        x <- data[,1]
        y <- data[,2]
        for(i in 1:Ntip(phy)){
            if(x[i]==0 & y[i]==0){states[i,1]=1}
            if(x[i]==0 & y[i]==1){states[i,2]=1}
            if(x[i]==1 & y[i]==0){states[i,3]=1}
            if(x[i]==1 & y[i]==1){states[i,4]=1}
            if(x[i]==2 & y[i]==0){states[i,c(1,3)]=1}
            if(x[i]==2 & y[i]==1){states[i,c(2,4)]=1}
            if(x[i]==0 & y[i]==2){states[i,c(1,2)]=1}
            if(x[i]==1 & y[i]==2){states[i,c(3,4)]=1}
            if(x[i]==2 & y[i]==2){states[i,1:4]=1}
        }
        compD <- matrix(0, nrow=nb.tip, ncol=4)
        compE <- matrix(0, nrow=nb.tip, ncol=4)
    }
    
    if(hidden.states == TRUE){
        states = matrix(0,Ntip(phy),32)
        x <- data[,1]
        y <- data[,2]
        for(i in 1:Ntip(phy)){
            if(x[i]==0 & y[i]==0){states[i,c(1,5,9,13,17,21,25,29)]=1}
            if(x[i]==0 & y[i]==1){states[i,c(2,6,10,14,18,22,26,30)]=1}
            if(x[i]==1 & y[i]==0){states[i,c(3,7,11,15,19,23,27,31)]=1}
            if(x[i]==1 & y[i]==1){states[i,c(4,8,12,16,20,24,28,32)]=1}
            if(x[i]==2 & y[i]==0){states[i,c(1,5,9,13,17,21,25,29, 3,7,11,15,19,23,27,31)]=1}
            if(x[i]==2 & y[i]==1){states[i,c(2,6,10,14,18,22,26,30, 4,8,12,16,20,24,28,32)]=1}
            if(x[i]==0 & y[i]==2){states[i,c(1,5,9,13,17,21,25,29, 2,6,10,14,18,22,26,30)]=1}
            if(x[i]==1 & y[i]==2){states[i,c(3,7,11,15,19,23,27,31, 4,8,12,16,20,24,28,32)]=1}
            if(x[i]==2 & y[i]==2){states[i,1:16]=1}
        }
        compD <- matrix(0, nrow=nb.tip, ncol=32)
        compE <- matrix(0, nrow=nb.tip, ncol=32)
    }
    
    if(hidden.states == "TEST1"){
        states = matrix(0,Ntip(phy),32)
        for(i in 1:Ntip(phy)){
            if(data[i]==1){states[i,c(1,5,9,13,17,21,25,29)]=1}
            if(data[i]==2){states[i,c(2,6,10,14,18,22,26,30)]=1}
            if(data[i]==3){states[i,c(3,7,11,15,19,23,27,31)]=1}
            if(data[i]==4){states[i,c(4,8,12,16,20,24,28,32)]=1}
        }
        compD <- matrix(0, nrow=nb.tip, ncol=32)
        compE <- matrix(0, nrow=nb.tip, ncol=32)
    }
    
    if(hidden.states == "TEST2"){
        states = matrix(0,Ntip(phy),4)
        for(i in 1:Ntip(phy)){
            if(data[i]==1){states[i,1]=1}
            if(data[i]==2){states[i,2]=1}
            if(data[i]==3){states[i,3]=1}
            if(data[i]==4){states[i,4]=1}
        }
        compD <- matrix(0, nrow=nb.tip, ncol=4)
        compE <- matrix(0, nrow=nb.tip, ncol=4)
        
    }
    
    #Initializes the tip sampling and sets internal nodes to be zero:
    ncols = dim(compD)[2]
    if(length(f) == 4){
        for(i in 1:(nb.tip)){
            compD[i,] <- f * states[i,]
            compE[i,] <- rep((1-f), ncols/4)
        }
    }else{
        for(i in 1:(nb.tip)){
            compD[i,] <- f[i] * states[i,]
            compE[i,] <- rep((1-f[i]), ncols/4)
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
    cols <- names(dat.tab)
    for (j in 1:(dim(compD)[2])){
        #dat.tab[data.table(c(1:nb.tip)), paste("compD", j, sep="_") := compD[,j]]
        set(dat.tab, 1:nb.tip, cols[6+j], compD[,j])
        #dat.tab[data.table(c(1:nb.tip)), paste("compE", j, sep="_") := compE[,j]]
        set(dat.tab, 1:nb.tip, cols[38+j], compE[,j])
    }
    return(dat.tab)
}


SingleChildProb <- function(cache, pars, compD, compE, start.time, end.time){
    if(cache$hidden.states == TRUE){
        yini <- c(E00A = compE[1], E01A = compE[2], E10A = compE[3], E11A = compE[4], E00B = compE[5], E01B = compE[6], E10B = compE[7], E11B = compE[8], E00C = compE[9], E01C = compE[10], E10C = compE[11], E11C = compE[12], E00D = compE[13], E01D = compE[14], E10D = compE[15], E11D = compE[16], E00E = compE[17], E01E = compE[18], E10E = compE[19], E11E = compE[20], E00F = compE[21], E01F = compE[22], E10F = compE[23], E11F = compE[24], E00G = compE[25], E01G = compE[26], E10G = compE[27], E11G = compE[28], E00H = compE[29], E01H = compE[30], E10H = compE[31], E11H = compE[32], D00A = compD[1], D01A = compD[2], D10A = compD[3], D11A = compD[4], D00B = compD[5], D01B = compD[6], D10B = compD[7], D11B = compD[8], D00C = compD[9], D01C = compD[10], D10C = compD[11], D11C = compD[12], D00D = compD[13], D01D = compD[14], D10D = compD[15], D11D = compD[16], D00E = compD[17], D01E = compD[18], D10E = compD[19], D11E = compD[20], D00F = compD[21], D01F = compD[22], D10F = compD[23], D11F = compD[24], D00G = compD[25], D01G = compD[26], D10G = compD[27], D11G = compD[28], D00H = compD[29], D01H = compD[30], D10H = compD[31], D11H = compD[32])
        times=c(start.time, end.time)
        #prob.subtree.cal.full <- lsoda(yini, times, func = "muhisse_derivs", pars, initfunc="initmod_muhisse", dll = "muhisse-ext-derivs", rtol=1e-8, atol=1e-8)
        prob.subtree.cal.full <- lsoda(yini, times, func = "muhisse_derivs", pars, initfunc="initmod_muhisse", dllname = "hisse", rtol=1e-8, atol=1e-8)
    }else{
        yini <-c(E00=compE[1], E01=compE[2], E10=compE[3], E11=compE[4], D00=compD[1], D01=compD[2], D10=compD[3], D11=compD[4])
        times=c(start.time, end.time)
        #prob.subtree.cal.full <- lsoda(yini, times, func = "musse_derivs", pars, initfunc="initmod_musse", dll = "canonical-musse-ext-derivs", rtol=1e-8, atol=1e-8)
        prob.subtree.cal.full <- lsoda(yini, times, func = "musse_derivs", pars, initfunc="initmod_musse", dllname = "hisse", rtol=1e-8, atol=1e-8)
    }
    ######## THIS CHECKS TO ENSURE THAT THE INTEGRATION WAS SUCCESSFUL ###########
    if(attributes(prob.subtree.cal.full)$istate[1] < 0){
        prob.subtree.cal <- prob.subtree.cal.full[-1,-1]
        if(cache$hidden.states == TRUE){
            prob.subtree.cal[33:64] <- cache$bad.likelihood
            return(prob.subtree.cal)
        }else{
            prob.subtree.cal[5:8] <- cache$bad.likelihood
            return(prob.subtree.cal)
        }
    }else{
        prob.subtree.cal <- prob.subtree.cal.full[-1,-1]
    }
    ##############################################################################
    
    if(cache$hidden.states == TRUE){
        if(any(is.nan(prob.subtree.cal[33:64]))){
            prob.subtree.cal[33:64] <- cache$bad.likelihood
            return(prob.subtree.cal)
        }
        #This is default and cannot change, but if we get a negative probability, discard the results:
        if(any(prob.subtree.cal[33:64] < 0)){
            prob.subtree.cal[33:64] <- cache$bad.likelihood
            return(prob.subtree.cal)
        }
        if(sum(prob.subtree.cal[33:64]) < cache$ode.eps){
            prob.subtree.cal[33:64] <- cache$bad.likelihood
            return(prob.subtree.cal)
        }
    }else{
        if(any(is.nan(prob.subtree.cal[5:8]))){
            prob.subtree.cal[5:8] <- cache$bad.likelihood
            return(prob.subtree.cal)
        }
        #This is default and cannot change, but if we get a negative probability, discard the results:
        if(any(prob.subtree.cal[5:8] < 0)){
            prob.subtree.cal[5:8] <- cache$bad.likelihood
            return(prob.subtree.cal)
        }
        if(sum(prob.subtree.cal[5:8]) < cache$ode.eps){
            prob.subtree.cal[5:8] <- cache$bad.likelihood
            return(prob.subtree.cal)
        }
    }
    return(prob.subtree.cal)
}


FocalNodeProb <- function(cache, pars, lambdas, dat.tab, generations){
    ### Ughy McUgherson. This is a must in order to pass CRAN checks: http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
    DesNode = NULL
    FocalNode = NULL
    . = NULL

    gens <- data.table(c(generations))
    #gens <- dat.tab[.(generations), which=TRUE]
    setkey(dat.tab, FocalNode)
    CurrentGenData <- dat.tab[gens]
    if(cache$hidden.states == TRUE){
        tmp <- t(apply(CurrentGenData, 1, function(z) SingleChildProb(cache, pars, z[7:38], z[39:70],  z[2], z[1])))
        v.mat <- matrix(tmp[seq(1,nrow(tmp)-1,2),33:64] * tmp[seq(2,nrow(tmp),2),33:64], length(unique(CurrentGenData$FocalNode)), 32)
        v.mat <- v.mat * matrix(lambdas, length(unique(CurrentGenData$FocalNode)), 32, byrow=TRUE)
        phi.mat <- matrix(tmp[seq(1,nrow(tmp)-1,2),1:32], length(unique(CurrentGenData$FocalNode)), 32)
        if(!is.null(cache$node)){
            if(any(cache$node %in% generations)){
                for(fix.index in 1:length(cache$node)){
                    if(cache$fix.type[fix.index] == "event"){
                        fixer.tmp = numeric(4)
                        fixer.tmp[cache$state[fix.index]] = 1
                        fixer = rep(fixer.tmp, 8)
                    }else{
                        fixer = numeric(32)
                        fixer[cache$state[fix.index]] = 1
                    }
                    v.mat[which(generations == cache$node[fix.index]),] <- v.mat[which(generations == cache$node[fix.index]),] * fixer
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
            #dat.tab[gens, cols[38+j] := phi.mat[,j]]
            set(dat.tab, rows, cols[38+j], phi.mat[,j])
        }
        dat.tab[gens, "comp" := tmp.comp]
    }else{
        tmp <- t(apply(CurrentGenData, 1, function(z) SingleChildProb(cache, pars, z[7:10], z[11:14],  z[2], z[1])))
        v.mat <- matrix(tmp[seq(1,nrow(tmp)-1,2),5:8] * tmp[seq(2,nrow(tmp),2),5:8], length(unique(CurrentGenData$FocalNode)), 4)
        v.mat <- v.mat * matrix(lambdas, length(unique(CurrentGenData$FocalNode)), 4, byrow=TRUE)
        phi.mat <- matrix(tmp[seq(1,nrow(tmp)-1,2),1:4], length(unique(CurrentGenData$FocalNode)), 4)
        if(!is.null(cache$node)){
            if(any(cache$node %in% generations)){
                for(fix.index in 1:length(cache$node)){
                    fixer = numeric(4)
                    fixer[cache$state[fix.index]] = 1
                    v.mat[which(generations == cache$node[fix.index]),] <- v.mat[which(generations == cache$node[fix.index]),] * fixer
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
            set(dat.tab, rows, cols[10+j], phi.mat[,j])
        }
        dat.tab[gens, "comp" := tmp.comp]
    }
    return(dat.tab)
}


#Have to calculate root prob separately because it is not a descendant in our table. Could add it, but I worry about the NA that is required.
GetRootProb <- function(cache, pars, lambdas, dat.tab, generations){
    ### Ughy McUgherson. This is a must in order to pass CRAN checks: http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
    FocalNode = NULL
    . = NULL

    gens <- data.table(c(generations))
    setkey(dat.tab, FocalNode)
    CurrentGenData <- dat.tab[gens]
    if(cache$hidden.states == TRUE){
        tmp <- t(apply(CurrentGenData, 1, function(z) SingleChildProb(cache, pars, z[7:38], z[39:70],  z[2], z[1])))
        v.mat <- matrix(tmp[seq(1,nrow(tmp)-1,2),33:64] * tmp[seq(2,nrow(tmp),2),33:64], length(unique(CurrentGenData$FocalNode)), 32)
        v.mat <- v.mat * matrix(lambdas, length(unique(CurrentGenData$FocalNode)), 32, byrow=TRUE)
        phi.mat <- matrix(tmp[seq(1,nrow(tmp)-1,2),1:32], length(unique(CurrentGenData$FocalNode)), 32)
        if(!is.null(cache$node)){
            if(any(cache$node %in% generations)){
                for(fix.index in 1:length(cache$node)){
                    if(cache$fix.type[fix.index] == "event"){
                        fixer.tmp = numeric(4)
                        fixer.tmp[cache$state[fix.index]] = 1
                        fixer = rep(fixer.tmp, 8)
                    }else{
                        fixer = numeric(32)
                        fixer[cache$state[fix.index]] = 1
                    }
                    v.mat[which(generations == cache$node[fix.index]),] <- v.mat[which(generations == cache$node[fix.index]),] * fixer
                }
            }
        }
    }else{
        tmp <- t(apply(CurrentGenData, 1, function(z) SingleChildProb(cache, pars, z[7:10], z[11:14],  z[2], z[1])))
        v.mat <- matrix(tmp[seq(1,nrow(tmp)-1,2),5:8] * tmp[seq(2,nrow(tmp),2),5:8], length(unique(CurrentGenData$FocalNode)), 4)
        v.mat <- v.mat * matrix(lambdas, length(unique(CurrentGenData$FocalNode)), 4, byrow=TRUE)
        phi.mat <- matrix(tmp[seq(1,nrow(tmp)-1,2),1:4], length(unique(CurrentGenData$FocalNode)), 4)
        if(!is.null(cache$node)){
            if(any(cache$node %in% generations)){
                for(fix.index in 1:length(cache$node)){
                    fixer = numeric(4)
                    fixer[cache$state[fix.index]] = 1
                    v.mat[which(generations == cache$node[fix.index]),] <- v.mat[which(generations == cache$node[fix.index]),] * fixer
                }
            }
        }
    }
    
    tmp.comp <- rowSums(v.mat)
    tmp.probs <- v.mat / tmp.comp
    
    return(cbind(tmp.comp, phi.mat, tmp.probs))
}



######################################################################################################################################
######################################################################################################################################
### The MuHiSSE type down pass that carries out the integration and returns the likelihood:
######################################################################################################################################
######################################################################################################################################

DownPassMuHisse <- function(dat.tab, gen, cache, condition.on.survival, root.type, root.p, get.phi=FALSE, node=NULL, state=NULL, fix.type=NULL) {
    
    ### Ughy McUgherson. This is a must in order to pass CRAN checks: http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
    DesNode = NULL
    compE = NULL

    ## Moved this here instead of above because it actually significantly improved the speed.
    if(cache$hidden.states == FALSE){
        pars <- c(cache$lambda00A, cache$lambda01A, cache$lambda10A, cache$lambda11A, cache$mu00A, cache$mu01A, cache$mu10A, cache$mu11A, cache$q00A_01A, cache$q00A_10A, cache$q00A_11A, cache$q01A_00A, cache$q01A_10A, cache$q01A_11A, cache$q10A_00A, cache$q10A_01A, cache$q10A_11A, cache$q11A_00A, cache$q11A_01A, cache$q11A_10A)
        lambda <- c(cache$lambda00A, cache$lambda01A, cache$lambda10A, cache$lambda11A)
    }else{
        pars <- c(cache$lambda00A, cache$lambda01A, cache$lambda10A, cache$lambda11A, cache$mu00A, cache$mu01A, cache$mu10A, cache$mu11A, cache$q00A_01A, cache$q00A_10A, cache$q00A_11A, cache$q01A_00A, cache$q01A_10A, cache$q01A_11A, cache$q10A_00A, cache$q10A_01A, cache$q10A_11A, cache$q11A_00A, cache$q11A_01A, cache$q11A_10A, cache$q00A_00B, cache$q00A_00C, cache$q00A_00D, cache$q00A_00E, cache$q00A_00F, cache$q00A_00G, cache$q00A_00H, cache$q01A_01B, cache$q01A_01C, cache$q01A_01D, cache$q01A_01E, cache$q01A_01F, cache$q01A_01G, cache$q01A_01H, cache$q10A_10B, cache$q10A_10C, cache$q10A_10D, cache$q10A_10E, cache$q10A_10F, cache$q10A_10G, cache$q10A_10H, cache$q11A_11B, cache$q11A_11C, cache$q11A_11D, cache$q11A_11E, cache$q11A_11F, cache$q11A_11G, cache$q11A_11H, cache$lambda00B, cache$lambda01B, cache$lambda10B, cache$lambda11B, cache$mu00B, cache$mu01B, cache$mu10B, cache$mu11B, cache$q00B_01B, cache$q00B_10B, cache$q00B_11B, cache$q01B_00B, cache$q01B_10B, cache$q01B_11B, cache$q10B_00B, cache$q10B_01B, cache$q10B_11B, cache$q11B_00B, cache$q11B_01B, cache$q11B_10B, cache$q00B_00A, cache$q00B_00C, cache$q00B_00D, cache$q00B_00E, cache$q00B_00F, cache$q00B_00G, cache$q00B_00H, cache$q01B_01A, cache$q01B_01C, cache$q01B_01D, cache$q01B_01E, cache$q01B_01F, cache$q01B_01G, cache$q01B_01H, cache$q10B_10A, cache$q10B_10C, cache$q10B_10D, cache$q10B_10E, cache$q10B_10F, cache$q10B_10G, cache$q10B_10H, cache$q11B_11A, cache$q11B_11C, cache$q11B_11D, cache$q11B_11E, cache$q11B_11F, cache$q11B_11G, cache$q11B_11H, cache$lambda00C, cache$lambda01C, cache$lambda10C, cache$lambda11C, cache$mu00C, cache$mu01C, cache$mu10C, cache$mu11C, cache$q00C_01C, cache$q00C_10C, cache$q00C_11C, cache$q01C_00C, cache$q01C_10C, cache$q01C_11C, cache$q10C_00C, cache$q10C_01C, cache$q10C_11C, cache$q11C_00C, cache$q11C_01C, cache$q11C_10C, cache$q00C_00A, cache$q00C_00B, cache$q00C_00D, cache$q00C_00E, cache$q00C_00F, cache$q00C_00G, cache$q00C_00H, cache$q01C_01A, cache$q01C_01B, cache$q01C_01D, cache$q01C_01E, cache$q01C_01F, cache$q01C_01G, cache$q01C_01H, cache$q10C_10A, cache$q10C_10B, cache$q10C_10D, cache$q10C_10E, cache$q10C_10F, cache$q10C_10G, cache$q10C_10H, cache$q11C_11A, cache$q11C_11B, cache$q11C_11D, cache$q11C_11E, cache$q11C_11F, cache$q11C_11G, cache$q11C_11H, cache$lambda00D, cache$lambda01D, cache$lambda10D, cache$lambda11D, cache$mu00D, cache$mu01D, cache$mu10D, cache$mu11D, cache$q00D_01D, cache$q00D_10D, cache$q00D_11D, cache$q01D_00D, cache$q01D_10D, cache$q01D_11D, cache$q10D_00D, cache$q10D_01D, cache$q10D_11D, cache$q11D_00D, cache$q11D_01D, cache$q11D_10D, cache$q00D_00A, cache$q00D_00B, cache$q00D_00C, cache$q00D_00E, cache$q00D_00F, cache$q00D_00G, cache$q00D_00H, cache$q01D_01A, cache$q01D_01B, cache$q01D_01C, cache$q01D_01E, cache$q01D_01F, cache$q01D_01G, cache$q01D_01H, cache$q10D_10A, cache$q10D_10B, cache$q10D_10C, cache$q10D_10E, cache$q10D_10F, cache$q10D_10G, cache$q10D_10H, cache$q11D_11A, cache$q11D_11B, cache$q11D_11C, cache$q11D_11E, cache$q11D_11F, cache$q11D_11G, cache$q11D_11H, cache$lambda00E, cache$lambda01E, cache$lambda10E, cache$lambda11E, cache$mu00E, cache$mu01E, cache$mu10E, cache$mu11E, cache$q00E_01E, cache$q00E_10E, cache$q00E_11E, cache$q01E_00E, cache$q01E_10E, cache$q01E_11E, cache$q10E_00E, cache$q10E_01E, cache$q10E_11E, cache$q11E_00E, cache$q11E_01E, cache$q11E_10E, cache$q00E_00A, cache$q00E_00B, cache$q00E_00C, cache$q00E_00D, cache$q00E_00F, cache$q00E_00G, cache$q00E_00H, cache$q01E_01A, cache$q01E_01B, cache$q01E_01C, cache$q01E_01D, cache$q01E_01F, cache$q01E_01G, cache$q01E_01H, cache$q10E_10A, cache$q10E_10B, cache$q10E_10C, cache$q10E_10D, cache$q10E_10F, cache$q10E_10G, cache$q10E_10H, cache$q11E_11A, cache$q11E_11B, cache$q11E_11C, cache$q11E_11D, cache$q11E_11F, cache$q11E_11G, cache$q11E_11H, cache$lambda00F, cache$lambda01F, cache$lambda10F, cache$lambda11F, cache$mu00F, cache$mu01F, cache$mu10F, cache$mu11F, cache$q00F_01F, cache$q00F_10F, cache$q00F_11F, cache$q01F_00F, cache$q01F_10F, cache$q01F_11F, cache$q10F_00F, cache$q10F_01F, cache$q10F_11F, cache$q11F_00F, cache$q11F_01F, cache$q11F_10F, cache$q00F_00A, cache$q00F_00B, cache$q00F_00C, cache$q00F_00D, cache$q00F_00E, cache$q00F_00G, cache$q00F_00H, cache$q01F_01A, cache$q01F_01B, cache$q01F_01C, cache$q01F_01D, cache$q01F_01E, cache$q01F_01G, cache$q01F_01H, cache$q10F_10A, cache$q10F_10B, cache$q10F_10C, cache$q10F_10D, cache$q10F_10E, cache$q10F_10G, cache$q10F_10H, cache$q11F_11A, cache$q11F_11B, cache$q11F_11C, cache$q11F_11D, cache$q11F_11E, cache$q11F_11G, cache$q11F_11H, cache$lambda00G, cache$lambda01G, cache$lambda10G, cache$lambda11G, cache$mu00G, cache$mu01G, cache$mu10G, cache$mu11G, cache$q00G_01G, cache$q00G_10G, cache$q00G_11G, cache$q01G_00G, cache$q01G_10G, cache$q01G_11G, cache$q10G_00G, cache$q10G_01G, cache$q10G_11G, cache$q11G_00G, cache$q11G_01G, cache$q11G_10G, cache$q00G_00A, cache$q00G_00B, cache$q00G_00C, cache$q00G_00D, cache$q00G_00E, cache$q00G_00F, cache$q00G_00H, cache$q01G_01A, cache$q01G_01B, cache$q01G_01C, cache$q01G_01D, cache$q01G_01E, cache$q01G_01F, cache$q01G_01H, cache$q10G_10A, cache$q10G_10B, cache$q10G_10C, cache$q10G_10D, cache$q10G_10E, cache$q10G_10F, cache$q10G_10H, cache$q11G_11A, cache$q11G_11B, cache$q11G_11C, cache$q11G_11D, cache$q11G_11E, cache$q11G_11F, cache$q11G_11H, cache$lambda00H, cache$lambda01H, cache$lambda10H, cache$lambda11H, cache$mu00H, cache$mu01H, cache$mu10H, cache$mu11H, cache$q00H_01H, cache$q00H_10H, cache$q00H_11H, cache$q01H_00H, cache$q01H_10H, cache$q01H_11H, cache$q10H_00H, cache$q10H_01H, cache$q10H_11H, cache$q11H_00H, cache$q11H_01H, cache$q11H_10H, cache$q00H_00A, cache$q00H_00B, cache$q00H_00C, cache$q00H_00D, cache$q00H_00E, cache$q00H_00F, cache$q00H_00G, cache$q01H_01A, cache$q01H_01B, cache$q01H_01C, cache$q01H_01D, cache$q01H_01E, cache$q01H_01F, cache$q01H_01G, cache$q10H_10A, cache$q10H_10B, cache$q10H_10C, cache$q10H_10D, cache$q10H_10E, cache$q10H_10F, cache$q10H_10G, cache$q11H_11A, cache$q11H_11B, cache$q11H_11C, cache$q11H_11D, cache$q11H_11E, cache$q11H_11F, cache$q11H_11G)
        lambda <- c(cache$lambda00A,cache$lambda01A,cache$lambda10A,cache$lambda11A,cache$lambda00B,cache$lambda01B,cache$lambda10B,cache$lambda11B,cache$lambda00C,cache$lambda01C,cache$lambda10C,cache$lambda11C,cache$lambda00D,cache$lambda01D,cache$lambda10D,cache$lambda11D,cache$lambda00E,cache$lambda01E,cache$lambda10E,cache$lambda11E,cache$lambda00F,cache$lambda01F,cache$lambda10F,cache$lambda11F,cache$lambda00G,cache$lambda01G,cache$lambda10G,cache$lambda11G,cache$lambda00H,cache$lambda01H,cache$lambda10H,cache$lambda11H)
    }
    
    nb.tip <- cache$nb.tip
    nb.node <- cache$nb.node
    TIPS <- 1:nb.tip
    for(i in 1:length(gen)){
        if(i == length(gen)){
            if(!is.null(node)){
                if(any(node %in% gen[[i]])){
                    cache$node <- node
                    cache$state <- state
                    cache$fix.type <- fix.type
                    res.tmp <- GetRootProb(cache=cache, pars=pars, lambdas=lambda, dat.tab=dat.tab, generations=gen[[i]])
                    cache$node <- NULL
                    cache$state <- NULL
                    cache$fix.type <- NULL
                }else{
                    res.tmp <- GetRootProb(cache=cache, pars=pars, lambdas=lambda, dat.tab=dat.tab, generations=gen[[i]])
                }
            }else{
                res.tmp <- GetRootProb(cache=cache, pars=pars, lambdas=lambda, dat.tab=dat.tab, generations=gen[[i]])
            }
            if(cache$hidden.states == TRUE){
                compD.root <- res.tmp[c(34:65)]
                compE.root <- res.tmp[c(2:33)]
            }else{
                compD.root <- res.tmp[c(6:9)]
                compE.root <- res.tmp[c(2:5)]
            }
            setkey(dat.tab, DesNode)
            comp <- dat.tab[["comp"]]
            comp <- c(comp[-TIPS], res.tmp[1])
        }else{
            if(!is.null(node)){
                if(any(node %in% gen[[i]])){
                    cache$node <- node
                    cache$state <- state
                    cache$fix.type <- fix.type
                    dat.tab <- FocalNodeProb(cache, pars=pars, lambdas=lambda, dat.tab, gen[[i]])
                    cache$node <- NULL
                    cache$state <- NULL
                    cache$fix.type <- NULL
                }else{
                    dat.tab <- FocalNodeProb(cache, pars=pars, lambdas=lambda, dat.tab, gen[[i]])
                }
            }else{
                dat.tab <- FocalNodeProb(cache, pars=pars, lambdas=lambda, dat.tab, gen[[i]])
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
            if(cache$hidden.states == FALSE){
                if(root.type == "madfitz"){
                    #lambda <- c(cache$lambda00A, cache$lambda01A, cache$lambda10A, cache$lambda11A)
                    compD.root <- compD.root / sum(root.p * lambda * (1 - compE.root)^2)
                    #Corrects for possibility that you have 0/0:
                    compD.root[which(is.na(compD.root))] = 0
                    loglik <- log(sum(compD.root * root.p)) + sum(log(comp))
                }else{
                    #lambda <- c(cache$lambda00A, cache$lambda01A, cache$lambda10A, cache$lambda11A)
                    compD.root <- (compD.root * root.p) / (lambda * (1 - compE.root)^2)
                    #Corrects for possibility that you have 0/0:
                    compD.root[which(is.na(compD.root))] = 0
                    loglik <- log(sum(compD.root)) + sum(log(comp))
                }
            }else{
                if(root.type == "madfitz"){
                    #lambda <- c(cache$lambda00A, cache$lambda01A, cache$lambda10A, cache$lambda11A, cache$lambda00B, cache$lambda01B, cache$lambda10B, cache$lambda11B, cache$lambda00C, cache$lambda01C, cache$lambda10C, cache$lambda11C, cache$lambda00D, cache$lambda01D, cache$lambda10D, cache$lambda11D, cache$lambda00E, cache$lambda01E, cache$lambda10E, cache$lambda11E, cache$lambda00F, cache$lambda01F, cache$lambda10F, cache$lambda11F, cache$lambda00G, cache$lambda01G, cache$lambda10G, cache$lambda11G, cache$lambda00H, cache$lambda01H, cache$lambda10H, cache$lambda11H)
                    compD.root <- compD.root / sum(root.p * lambda * (1 - compE.root)^2)
                    #Corrects for possibility that you have 0/0:
                    compD.root[which(is.na(compD.root))] = 0
                    loglik <- log(sum(compD.root * root.p)) + sum(log(comp))
                }else{
                    #lambda <- c(cache$lambda00A, cache$lambda01A, cache$lambda10A, cache$lambda11A, cache$lambda00B, cache$lambda01B, cache$lambda10B, cache$lambda11B, cache$lambda00C, cache$lambda01C, cache$lambda10C, cache$lambda11C, cache$lambda00D, cache$lambda01D, cache$lambda10D, cache$lambda11D, cache$lambda00E, cache$lambda01E, cache$lambda10E, cache$lambda11E, cache$lambda00F, cache$lambda01F, cache$lambda10F, cache$lambda11F, cache$lambda00G, cache$lambda01G, cache$lambda10G, cache$lambda11G, cache$lambda00H, cache$lambda01H, cache$lambda10H, cache$lambda11H)
                    compD.root <- (compD.root * root.p) / (lambda * (1 - compE.root)^2)
                    #Corrects for possibility that you have 0/0:
                    compD.root[which(is.na(compD.root))] = 0
                    loglik <- log(sum(compD.root)) + sum(log(comp))
                }
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
### Cache object for storing parameters that are used throughout MuHiSSE:
######################################################################################################################################
######################################################################################################################################

ParametersToPassMuHiSSE <- function(model.vec, hidden.states, nb.tip, nb.node, bad.likelihood, ode.eps){
    #Provides an initial object that contains all the parameters to be passed among functions. This will also be used to pass other things are we move down the tree (see DownPassGeoSSE):
    obj <- NULL
    
    obj$hidden.states <- hidden.states
    obj$nb.tip <- nb.tip
    obj$nb.node <- nb.node
    obj$bad.likelihood <- bad.likelihood
    obj$ode.eps <- ode.eps
    
    ##Hidden State A
    obj$lambda00A = model.vec[1] / (1 + model.vec[5])
    obj$lambda01A = model.vec[2] / (1 + model.vec[6])
    obj$lambda10A = model.vec[3] / (1 + model.vec[7])
    obj$lambda11A = model.vec[4] / (1 + model.vec[8])
    obj$mu00A = (model.vec[5] * model.vec[1]) / (1 + model.vec[5])
    obj$mu01A = (model.vec[6] * model.vec[2]) / (1 + model.vec[6])
    obj$mu10A = (model.vec[7] * model.vec[3]) / (1 + model.vec[7])
    obj$mu11A = (model.vec[8] * model.vec[4]) / (1 + model.vec[8])
    obj$q00A_01A = model.vec[9]
    obj$q00A_10A = model.vec[10]
    obj$q00A_11A = model.vec[11]
    obj$q01A_00A = model.vec[12]
    obj$q01A_10A = model.vec[13]
    obj$q01A_11A = model.vec[14]
    obj$q10A_00A = model.vec[15]
    obj$q10A_01A = model.vec[16]
    obj$q10A_11A = model.vec[17]
    obj$q11A_00A = model.vec[18]
    obj$q11A_01A = model.vec[19]
    obj$q11A_10A = model.vec[20]
    
    obj$q00A_00B = model.vec[21]
    obj$q00A_00C = model.vec[22]
    obj$q00A_00D = model.vec[23]
    obj$q00A_00E = model.vec[24]
    obj$q00A_00F = model.vec[25]
    obj$q00A_00G = model.vec[26]
    obj$q00A_00H = model.vec[27]
    obj$q01A_01B = model.vec[28]
    obj$q01A_01C = model.vec[29]
    obj$q01A_01D = model.vec[30]
    obj$q01A_01E = model.vec[31]
    obj$q01A_01F = model.vec[32]
    obj$q01A_01G = model.vec[33]
    obj$q01A_01H = model.vec[34]
    obj$q10A_10B = model.vec[35]
    obj$q10A_10C = model.vec[36]
    obj$q10A_10D = model.vec[37]
    obj$q10A_10E = model.vec[38]
    obj$q10A_10F = model.vec[39]
    obj$q10A_10G = model.vec[40]
    obj$q10A_10H = model.vec[41]
    obj$q11A_11B = model.vec[42]
    obj$q11A_11C = model.vec[43]
    obj$q11A_11D = model.vec[44]
    obj$q11A_11E = model.vec[45]
    obj$q11A_11F = model.vec[46]
    obj$q11A_11G = model.vec[47]
    obj$q11A_11H = model.vec[48]
    
    ##Hidden State B
    obj$lambda00B = model.vec[49] / (1 + model.vec[53])
    obj$lambda01B = model.vec[50] / (1 + model.vec[54])
    obj$lambda10B = model.vec[51] / (1 + model.vec[55])
    obj$lambda11B = model.vec[52] / (1 + model.vec[56])
    obj$mu00B = (model.vec[53] * model.vec[49]) / (1 + model.vec[53])
    obj$mu01B = (model.vec[54] * model.vec[50]) / (1 + model.vec[54])
    obj$mu10B = (model.vec[55] * model.vec[51]) / (1 + model.vec[55])
    obj$mu11B = (model.vec[56] * model.vec[52]) / (1 + model.vec[56])
    obj$q00B_01B = model.vec[57]
    obj$q00B_10B = model.vec[58]
    obj$q00B_11B = model.vec[59]
    obj$q01B_00B = model.vec[60]
    obj$q01B_10B = model.vec[61]
    obj$q01B_11B = model.vec[62]
    obj$q10B_00B = model.vec[63]
    obj$q10B_01B = model.vec[64]
    obj$q10B_11B = model.vec[65]
    obj$q11B_00B = model.vec[66]
    obj$q11B_01B = model.vec[67]
    obj$q11B_10B = model.vec[68]
    
    obj$q00B_00A = model.vec[69]
    obj$q00B_00C = model.vec[70]
    obj$q00B_00D = model.vec[71]
    obj$q00B_00E = model.vec[72]
    obj$q00B_00F = model.vec[73]
    obj$q00B_00G = model.vec[74]
    obj$q00B_00H = model.vec[75]
    obj$q01B_01A = model.vec[76]
    obj$q01B_01C = model.vec[77]
    obj$q01B_01D = model.vec[78]
    obj$q01B_01E = model.vec[79]
    obj$q01B_01F = model.vec[80]
    obj$q01B_01G = model.vec[81]
    obj$q01B_01H = model.vec[82]
    obj$q10B_10A = model.vec[83]
    obj$q10B_10C = model.vec[84]
    obj$q10B_10D = model.vec[85]
    obj$q10B_10E = model.vec[86]
    obj$q10B_10F = model.vec[87]
    obj$q10B_10G = model.vec[88]
    obj$q10B_10H = model.vec[89]
    obj$q11B_11A = model.vec[90]
    obj$q11B_11C = model.vec[91]
    obj$q11B_11D = model.vec[92]
    obj$q11B_11E = model.vec[93]
    obj$q11B_11F = model.vec[94]
    obj$q11B_11G = model.vec[95]
    obj$q11B_11H = model.vec[96]
    
    ##Hidden State C
    obj$lambda00C = model.vec[97] / (1 + model.vec[101])
    obj$lambda01C = model.vec[98] / (1 + model.vec[102])
    obj$lambda10C = model.vec[99] / (1 + model.vec[103])
    obj$lambda11C = model.vec[100] / (1 + model.vec[104])
    obj$mu00C = (model.vec[101] * model.vec[97]) / (1 + model.vec[101])
    obj$mu01C = (model.vec[102] * model.vec[98]) / (1 + model.vec[102])
    obj$mu10C = (model.vec[103] * model.vec[99]) / (1 + model.vec[103])
    obj$mu11C = (model.vec[104] * model.vec[100]) / (1 + model.vec[104])
    obj$q00C_01C = model.vec[105]
    obj$q00C_10C = model.vec[106]
    obj$q00C_11C = model.vec[107]
    obj$q01C_00C = model.vec[108]
    obj$q01C_10C = model.vec[109]
    obj$q01C_11C = model.vec[110]
    obj$q10C_00C = model.vec[111]
    obj$q10C_01C = model.vec[112]
    obj$q10C_11C = model.vec[113]
    obj$q11C_00C = model.vec[114]
    obj$q11C_01C = model.vec[115]
    obj$q11C_10C = model.vec[116]
    
    obj$q00C_00A = model.vec[117]
    obj$q00C_00B = model.vec[118]
    obj$q00C_00D = model.vec[119]
    obj$q00C_00E = model.vec[120]
    obj$q00C_00F = model.vec[121]
    obj$q00C_00G = model.vec[122]
    obj$q00C_00H = model.vec[123]
    obj$q01C_01A = model.vec[124]
    obj$q01C_01B = model.vec[125]
    obj$q01C_01D = model.vec[126]
    obj$q01C_01E = model.vec[127]
    obj$q01C_01F = model.vec[128]
    obj$q01C_01G = model.vec[129]
    obj$q01C_01H = model.vec[130]
    obj$q10C_10A = model.vec[131]
    obj$q10C_10B = model.vec[132]
    obj$q10C_10D = model.vec[133]
    obj$q10C_10E = model.vec[134]
    obj$q10C_10F = model.vec[135]
    obj$q10C_10G = model.vec[136]
    obj$q10C_10H = model.vec[137]
    obj$q11C_11A = model.vec[138]
    obj$q11C_11B = model.vec[139]
    obj$q11C_11D = model.vec[140]
    obj$q11C_11E = model.vec[141]
    obj$q11C_11F = model.vec[142]
    obj$q11C_11G = model.vec[143]
    obj$q11C_11H = model.vec[144]
    
    ##Hidden State D
    obj$lambda00D = model.vec[145] / (1 + model.vec[149])
    obj$lambda01D = model.vec[146] / (1 + model.vec[150])
    obj$lambda10D = model.vec[147] / (1 + model.vec[151])
    obj$lambda11D = model.vec[148] / (1 + model.vec[152])
    obj$mu00D = (model.vec[149] * model.vec[145]) / (1 + model.vec[149])
    obj$mu01D = (model.vec[150] * model.vec[146]) / (1 + model.vec[150])
    obj$mu10D = (model.vec[151] * model.vec[147]) / (1 + model.vec[151])
    obj$mu11D = (model.vec[152] * model.vec[148]) / (1 + model.vec[152])
    obj$q00D_01D = model.vec[153]
    obj$q00D_10D = model.vec[154]
    obj$q00D_11D = model.vec[155]
    obj$q01D_00D = model.vec[156]
    obj$q01D_10D = model.vec[157]
    obj$q01D_11D = model.vec[158]
    obj$q10D_00D = model.vec[159]
    obj$q10D_01D = model.vec[160]
    obj$q10D_11D = model.vec[161]
    obj$q11D_00D = model.vec[162]
    obj$q11D_01D = model.vec[163]
    obj$q11D_10D = model.vec[164]
    
    obj$q00D_00A = model.vec[165]
    obj$q00D_00B = model.vec[166]
    obj$q00D_00C = model.vec[167]
    obj$q00D_00E = model.vec[168]
    obj$q00D_00F = model.vec[169]
    obj$q00D_00G = model.vec[170]
    obj$q00D_00H = model.vec[171]
    obj$q01D_01A = model.vec[172]
    obj$q01D_01B = model.vec[173]
    obj$q01D_01C = model.vec[174]
    obj$q01D_01E = model.vec[175]
    obj$q01D_01F = model.vec[176]
    obj$q01D_01G = model.vec[177]
    obj$q01D_01H = model.vec[178]
    obj$q10D_10A = model.vec[179]
    obj$q10D_10B = model.vec[180]
    obj$q10D_10C = model.vec[181]
    obj$q10D_10E = model.vec[182]
    obj$q10D_10F = model.vec[183]
    obj$q10D_10G = model.vec[184]
    obj$q10D_10H = model.vec[185]
    obj$q11D_11A = model.vec[186]
    obj$q11D_11B = model.vec[187]
    obj$q11D_11C = model.vec[188]
    obj$q11D_11E = model.vec[189]
    obj$q11D_11F = model.vec[190]
    obj$q11D_11G = model.vec[191]
    obj$q11D_11H = model.vec[192]
    
    ##Hidden State E
    obj$lambda00E = model.vec[193] / (1 + model.vec[197])
    obj$lambda01E = model.vec[194] / (1 + model.vec[198])
    obj$lambda10E = model.vec[195] / (1 + model.vec[199])
    obj$lambda11E = model.vec[196] / (1 + model.vec[200])
    obj$mu00E = (model.vec[197] * model.vec[193]) / (1 + model.vec[197])
    obj$mu01E = (model.vec[198] * model.vec[194]) / (1 + model.vec[198])
    obj$mu10E = (model.vec[199] * model.vec[195]) / (1 + model.vec[199])
    obj$mu11E = (model.vec[200] * model.vec[196]) / (1 + model.vec[200])
    obj$q00E_01E = model.vec[201]
    obj$q00E_10E = model.vec[202]
    obj$q00E_11E = model.vec[203]
    obj$q01E_00E = model.vec[204]
    obj$q01E_10E = model.vec[205]
    obj$q01E_11E = model.vec[206]
    obj$q10E_00E = model.vec[207]
    obj$q10E_01E = model.vec[208]
    obj$q10E_11E = model.vec[209]
    obj$q11E_00E = model.vec[210]
    obj$q11E_01E = model.vec[211]
    obj$q11E_10E = model.vec[212]
    
    obj$q00E_00A = model.vec[213]
    obj$q00E_00B = model.vec[214]
    obj$q00E_00C = model.vec[215]
    obj$q00E_00D = model.vec[216]
    obj$q00E_00F = model.vec[217]
    obj$q00E_00G = model.vec[218]
    obj$q00E_00H = model.vec[219]
    obj$q01E_01A = model.vec[220]
    obj$q01E_01B = model.vec[221]
    obj$q01E_01C = model.vec[222]
    obj$q01E_01D = model.vec[223]
    obj$q01E_01F = model.vec[224]
    obj$q01E_01G = model.vec[225]
    obj$q01E_01H = model.vec[226]
    obj$q10E_10A = model.vec[227]
    obj$q10E_10B = model.vec[228]
    obj$q10E_10C = model.vec[229]
    obj$q10E_10D = model.vec[230]
    obj$q10E_10F = model.vec[231]
    obj$q10E_10G = model.vec[232]
    obj$q10E_10H = model.vec[233]
    obj$q11E_11A = model.vec[234]
    obj$q11E_11B = model.vec[235]
    obj$q11E_11C = model.vec[236]
    obj$q11E_11D = model.vec[237]
    obj$q11E_11F = model.vec[238]
    obj$q11E_11G = model.vec[239]
    obj$q11E_11H = model.vec[240]
    
    ##Hidden State F
    obj$lambda00F = model.vec[241] / (1 + model.vec[245])
    obj$lambda01F = model.vec[242] / (1 + model.vec[246])
    obj$lambda10F = model.vec[243] / (1 + model.vec[247])
    obj$lambda11F = model.vec[244] / (1 + model.vec[248])
    obj$mu00F = (model.vec[245] * model.vec[241]) / (1 + model.vec[245])
    obj$mu01F = (model.vec[246] * model.vec[242]) / (1 + model.vec[246])
    obj$mu10F = (model.vec[247] * model.vec[243]) / (1 + model.vec[247])
    obj$mu11F = (model.vec[248] * model.vec[244]) / (1 + model.vec[248])
    obj$q00F_01F = model.vec[249]
    obj$q00F_10F = model.vec[250]
    obj$q00F_11F = model.vec[251]
    obj$q01F_00F = model.vec[252]
    obj$q01F_10F = model.vec[253]
    obj$q01F_11F = model.vec[254]
    obj$q10F_00F = model.vec[255]
    obj$q10F_01F = model.vec[256]
    obj$q10F_11F = model.vec[257]
    obj$q11F_00F = model.vec[258]
    obj$q11F_01F = model.vec[259]
    obj$q11F_10F = model.vec[260]
    
    obj$q00F_00A = model.vec[261]
    obj$q00F_00B = model.vec[262]
    obj$q00F_00C = model.vec[263]
    obj$q00F_00D = model.vec[264]
    obj$q00F_00E = model.vec[265]
    obj$q00F_00G = model.vec[266]
    obj$q00F_00H = model.vec[267]
    obj$q01F_01A = model.vec[268]
    obj$q01F_01B = model.vec[269]
    obj$q01F_01C = model.vec[270]
    obj$q01F_01D = model.vec[271]
    obj$q01F_01E = model.vec[272]
    obj$q01F_01G = model.vec[273]
    obj$q01F_01H = model.vec[274]
    obj$q10F_10A = model.vec[275]
    obj$q10F_10B = model.vec[276]
    obj$q10F_10C = model.vec[277]
    obj$q10F_10D = model.vec[278]
    obj$q10F_10E = model.vec[279]
    obj$q10F_10G = model.vec[280]
    obj$q10F_10H = model.vec[281]
    obj$q11F_11A = model.vec[282]
    obj$q11F_11B = model.vec[283]
    obj$q11F_11C = model.vec[284]
    obj$q11F_11D = model.vec[285]
    obj$q11F_11E = model.vec[286]
    obj$q11F_11G = model.vec[287]
    obj$q11F_11H = model.vec[288]
    
    ##Hidden State G
    obj$lambda00G = model.vec[289] / (1 + model.vec[293])
    obj$lambda01G = model.vec[290] / (1 + model.vec[294])
    obj$lambda10G = model.vec[291] / (1 + model.vec[295])
    obj$lambda11G = model.vec[292] / (1 + model.vec[296])
    obj$mu00G = (model.vec[293] * model.vec[289]) / (1 + model.vec[293])
    obj$mu01G = (model.vec[294] * model.vec[290]) / (1 + model.vec[294])
    obj$mu10G = (model.vec[295] * model.vec[291]) / (1 + model.vec[295])
    obj$mu11G = (model.vec[296] * model.vec[292]) / (1 + model.vec[296])
    obj$q00G_01G = model.vec[297]
    obj$q00G_10G = model.vec[298]
    obj$q00G_11G = model.vec[299]
    obj$q01G_00G = model.vec[300]
    obj$q01G_10G = model.vec[301]
    obj$q01G_11G = model.vec[302]
    obj$q10G_00G = model.vec[303]
    obj$q10G_01G = model.vec[304]
    obj$q10G_11G = model.vec[305]
    obj$q11G_00G = model.vec[306]
    obj$q11G_01G = model.vec[307]
    obj$q11G_10G = model.vec[308]
    
    obj$q00G_00A = model.vec[309]
    obj$q00G_00B = model.vec[310]
    obj$q00G_00C = model.vec[311]
    obj$q00G_00D = model.vec[312]
    obj$q00G_00E = model.vec[313]
    obj$q00G_00F = model.vec[314]
    obj$q00G_00H = model.vec[315]
    obj$q01G_01A = model.vec[316]
    obj$q01G_01B = model.vec[317]
    obj$q01G_01C = model.vec[318]
    obj$q01G_01D = model.vec[319]
    obj$q01G_01E = model.vec[320]
    obj$q01G_01F = model.vec[321]
    obj$q01G_01H = model.vec[322]
    obj$q10G_10A = model.vec[323]
    obj$q10G_10B = model.vec[324]
    obj$q10G_10C = model.vec[325]
    obj$q10G_10D = model.vec[326]
    obj$q10G_10E = model.vec[327]
    obj$q10G_10F = model.vec[328]
    obj$q10G_10H = model.vec[329]
    obj$q11G_11A = model.vec[330]
    obj$q11G_11B = model.vec[331]
    obj$q11G_11C = model.vec[332]
    obj$q11G_11D = model.vec[333]
    obj$q11G_11E = model.vec[334]
    obj$q11G_11F = model.vec[335]
    obj$q11G_11H = model.vec[336]
    
    ##Hidden State H
    obj$lambda00H = model.vec[337] / (1 + model.vec[341])
    obj$lambda01H = model.vec[338] / (1 + model.vec[342])
    obj$lambda10H = model.vec[339] / (1 + model.vec[343])
    obj$lambda11H = model.vec[340] / (1 + model.vec[344])
    obj$mu00H = (model.vec[341] * model.vec[337]) / (1 + model.vec[341])
    obj$mu01H = (model.vec[342] * model.vec[338]) / (1 + model.vec[342])
    obj$mu10H = (model.vec[343] * model.vec[339]) / (1 + model.vec[343])
    obj$mu11H = (model.vec[344] * model.vec[340]) / (1 + model.vec[344])
    obj$q00H_01H = model.vec[345]
    obj$q00H_10H = model.vec[346]
    obj$q00H_11H = model.vec[347]
    obj$q01H_00H = model.vec[348]
    obj$q01H_10H = model.vec[349]
    obj$q01H_11H = model.vec[350]
    obj$q10H_00H = model.vec[351]
    obj$q10H_01H = model.vec[352]
    obj$q10H_11H = model.vec[353]
    obj$q11H_00H = model.vec[354]
    obj$q11H_01H = model.vec[355]
    obj$q11H_10H = model.vec[356]
    
    obj$q00H_00A = model.vec[357]
    obj$q00H_00B = model.vec[358]
    obj$q00H_00C = model.vec[359]
    obj$q00H_00D = model.vec[360]
    obj$q00H_00E = model.vec[361]
    obj$q00H_00F = model.vec[362]
    obj$q00H_00G = model.vec[363]
    obj$q01H_01A = model.vec[364]
    obj$q01H_01B = model.vec[365]
    obj$q01H_01C = model.vec[366]
    obj$q01H_01D = model.vec[367]
    obj$q01H_01E = model.vec[368]
    obj$q01H_01F = model.vec[369]
    obj$q01H_01G = model.vec[370]
    obj$q10H_10A = model.vec[371]
    obj$q10H_10B = model.vec[372]
    obj$q10H_10C = model.vec[373]
    obj$q10H_10D = model.vec[374]
    obj$q10H_10E = model.vec[375]
    obj$q10H_10F = model.vec[376]
    obj$q10H_10G = model.vec[377]
    obj$q11H_11A = model.vec[378]
    obj$q11H_11B = model.vec[379]
    obj$q11H_11C = model.vec[380]
    obj$q11H_11D = model.vec[381]
    obj$q11H_11E = model.vec[382]
    obj$q11H_11F = model.vec[383]
    obj$q11H_11G = model.vec[384]
    
    return(obj)
}


print.muhisse.fit <- function(x,...){
    ## Function to print a "muhisse.fit" object.
    set.zero <- max( x$index.par )
    ## Keep only the parameters estimated:
    par.list <- x$solution[ !x$index.par == set.zero ]
    ntips <- Ntip( x$phy )
    nstates <- ncol( x$trans.matrix )/4
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


