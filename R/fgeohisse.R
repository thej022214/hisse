
#library(ape)
#library(deSolve)
#library(subplex)
#library(phytools)
#library(nloptr)
#library(GenSA)
#library(data.table)
#dyn.load("fcanonical_geosse-ext-derivs.so")
#dyn.load("fgeohisse-ext-derivs.so")
#dyn.load("fnotclasse-more-ext-derivs.so")
#dyn.load("fnotclasse-ext-derivs.so")


######################################################################################################################################
######################################################################################################################################
### fGeoHiSSE -- Very fast version of the expanded set of GeoSSE models
######################################################################################################################################
######################################################################################################################################

fGeoHiSSE <- function(phy, data, f=c(1,1,1), turnover=c(1,2,3), extinct.frac=c(1,2), hidden.areas=FALSE, trans.rate=NULL, assume.cladogenetic=TRUE, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, sann=FALSE, sann.its=10000, bounded.search=TRUE, max.tol=.Machine$double.eps^.50, mag.san.start=0.5, starting.vals=NULL, turnover.upper=1000, extinct.frac.upper=1000, trans.upper=100, restart.obj=NULL, ode.eps=0){
    
    ## Temporary fix for the current BUG:
    if( !is.null(phy$node.label) ) phy$node.label <- NULL
    
    if(!is.null(root.p)) {
        root.p <- root.p / sum(root.p)
        if(hidden.areas ==TRUE & length(root.p)==2){
            root.p <- rep(root.p, 2)
            root.p <- root.p / sum(root.p)
            warning("For hidden areas, you need to specify the root.p for all four hidden states. We have adjusted it so that there's equal chance for 00A as 00B, and for 11A as 11B")
        }
    }
    
    if(!root.type == "madfitz" & !root.type == "herr_als"){
        stop("Check that you specified a proper root.type option. Options are 'madfitz' or 'herr_als'. See help for more details.", call.=FALSE)
    }
    
    if(is.null(trans.rate)){
        stop("Rate matrix needed. See TransMatMakerGeoHiSSE() to create one.")
    }
    
    if(hidden.areas == TRUE & dim(trans.rate)[1]<3){
        stop("You chose a hidden state but this is not reflected in the transition matrix")
    }
        
    ## Return error message if the data is not in the correct format.
    if( !inherits(data, what = c("matrix","data.frame")) ){
        stop("'data' needs to be a matrix or data.frame with 2 columns. See help.")
    }
    if( !ncol( data ) == 2 ){
        stop("'data' needs to be a matrix or data.frame with 2 columns. See help.")
    }
    ## Check if the states are %in% 0:2:
    states.check <- all( as.numeric(data[,2]) %in% 0:2 )
    if( !states.check ){
        stop("states need to be one of 0, 1, or 2. See help.")
    }
    
    ## Check if 'hidden.areas' parameter is congruent with the turnover vector:
    if( length(turnover) > 3 & !hidden.areas ){
        stop("Turnover has more than 3 elements but 'hidden.areas' was set to FALSE. Please set 'hidden.areas' to TRUE if the model include more than one rate class.")
    }
    if( length(turnover) == 3 & hidden.areas ){
        stop("Turnover has only 3 elements but 'hidden.areas' was set to TRUE. Please set 'hidden.areas' to FALSE if the model does not include hidden rate classes.")
    }
    
    pars <- numeric(380)
    
    if(dim(trans.rate)[2]==3){
        rate.cats <- 1
        pars.tmp <- turnover
        extinct.frac.tmp <- extinct.frac
        extinct.frac.tmp[which(extinct.frac.tmp > 0)] = (extinct.frac.tmp[which( extinct.frac.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, extinct.frac.tmp)
        trans.tmp <- c(trans.rate["(00)", "(11)"], trans.rate["(00)", "(01)"], trans.rate["(11)", "(00)"], trans.rate["(11)", "(01)"],  trans.rate["(01)", "(00)"],  trans.rate["(01)", "(11)"])
        trans.tmp[which(trans.tmp > 0)] = (trans.tmp[which(trans.tmp > 0)] + max(pars.tmp))
        category.rates.unique <- 0
        pars.tmp <- c(pars.tmp, trans.tmp)
        pars[1:11] <- pars.tmp
    }
    
    if(dim(trans.rate)[2]==6){
        rate.cats <- 2
        pars.tmp <- turnover
        extinct.frac.tmp <- extinct.frac
        extinct.frac.tmp[which(extinct.frac.tmp > 0)] = (extinct.frac.tmp[which( extinct.frac.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, extinct.frac.tmp)
        for.late.adjust <- max(pars.tmp)
        rows <- c("(00A)", "(00A)",  "(11A)", "(11A)",  "(01A)", "(01A)", "(00B)", "(00B)",  "(11B)", "(11B)",  "(01B)", "(01B)")
        cols <- c("(11A)", "(01A)", "(00A)", "(01A)", "(00A)",  "(11A)",  "(11B)", "(01B)", "(00B)", "(01B)", "(00B)",  "(11B)")
        trans.tmp <- trans.rate[cbind(rows,cols)]
        trans.tmp[which(trans.tmp > 0)] = (trans.tmp[which(trans.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, trans.tmp)
        rows <- c(rep("(00A)", 1), rep("(11A)", 1), rep("(01A)", 1), rep("(00B)", 1), rep("(11B)", 1), rep("(01B)", 1))
        cols <- c("(00B)", "(11B)", "(01B)", "(00A)", "(11A)", "(01A)")
        category.tmp <- trans.rate[cbind(rows,cols)]
        category.rate.shift <- category.tmp + for.late.adjust
        category.rate.shift[is.na(category.rate.shift)] <- 0
        category.rate.shiftA <- c(category.rate.shift[1], rep(0,8), category.rate.shift[2], rep(0,8), category.rate.shift[3], rep(0,8))
        category.rate.shiftB <- c(category.rate.shift[4], rep(0,8), category.rate.shift[5], rep(0,8), category.rate.shift[6], rep(0,8))
        category.rates.all <- c(category.rate.shiftA, category.rate.shiftB)
        category.rates.unique <- length(unique(category.rates.all[category.rates.all>0]))
        pars.tmp <- c(turnover[1:3], extinct.frac.tmp[1:2], trans.tmp[1:6], category.rate.shiftA, turnover[4:6], extinct.frac.tmp[3:4], trans.tmp[7:12], category.rate.shiftB)
        pars[1:length(pars.tmp)] <- pars.tmp
    }
    
    if(dim(trans.rate)[2]==9){
        rate.cats <- 3
        pars.tmp <- turnover
        extinct.frac.tmp <- extinct.frac
        extinct.frac.tmp[which(extinct.frac.tmp > 0)] = (extinct.frac.tmp[which( extinct.frac.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, extinct.frac.tmp)
        for.late.adjust <- max(pars.tmp)
        rows <- c("(00A)", "(00A)",  "(11A)", "(11A)",  "(01A)", "(01A)", "(00B)", "(00B)",  "(11B)", "(11B)",  "(01B)", "(01B)", "(00C)", "(00C)",  "(11C)", "(11C)",  "(01C)", "(01C)")
        cols <- c("(11A)", "(01A)", "(00A)", "(01A)", "(00A)",  "(11A)",  "(11B)", "(01B)", "(00B)", "(01B)", "(00B)",  "(11B)",  "(11C)", "(01C)", "(00C)", "(01C)", "(00C)",  "(11C)")
        trans.tmp <- trans.rate[cbind(rows,cols)]
        trans.tmp[which(trans.tmp > 0)] = (trans.tmp[which(trans.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, trans.tmp)
        rows <- c(rep("(00A)", 2), rep("(11A)", 2), rep("(01A)", 2), rep("(00B)", 2), rep("(11B)", 2), rep("(01B)", 2), rep("(00C)", 2), rep("(11C)", 2), rep("(01C)", 2))
        cols <- c("(00B)", "(00C)", "(11B)", "(11C)", "(01B)", "(01C)", "(00A)", "(00C)", "(11A)", "(11C)", "(01A)", "(01C)", "(00A)", "(00B)", "(11A)", "(11B)", "(01A)", "(01B)")
        category.tmp <- trans.rate[cbind(rows,cols)]
        category.rate.shift <- category.tmp + for.late.adjust
        category.rate.shift[is.na(category.rate.shift)] <- 0
        category.rate.shiftA <- c(category.rate.shift[1:2], rep(0,7), category.rate.shift[3:4], rep(0,7), category.rate.shift[5:6], rep(0,7))
        category.rate.shiftB <- c(category.rate.shift[7:8], rep(0,7), category.rate.shift[9:10], rep(0,7), category.rate.shift[11:12], rep(0,7))
        category.rate.shiftC <- c(category.rate.shift[13:14], rep(0,7), category.rate.shift[15:16], rep(0,7), category.rate.shift[17:18], rep(0,7))
        category.rates.all <- c(category.rate.shiftA, category.rate.shiftB, category.rate.shiftC)
        category.rates.unique <- length(unique(category.rates.all[category.rates.all>0]))
        pars.tmp <- c(turnover[1:3], extinct.frac.tmp[1:2], trans.tmp[1:6], category.rate.shiftA, turnover[4:6], extinct.frac.tmp[3:4], trans.tmp[7:12], category.rate.shiftB, turnover[7:9], extinct.frac.tmp[5:6], trans.tmp[13:18], category.rate.shiftC)
        pars[1:length(pars.tmp)] <- pars.tmp
    }
    
    if(dim(trans.rate)[2]==12){
        rate.cats <- 4
        pars.tmp <- turnover
        extinct.frac.tmp <- extinct.frac
        extinct.frac.tmp[which(extinct.frac.tmp > 0)] = (extinct.frac.tmp[which( extinct.frac.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, extinct.frac.tmp)
        for.late.adjust <- max(pars.tmp)
        rows <- c("(00A)", "(00A)",  "(11A)", "(11A)",  "(01A)", "(01A)", "(00B)", "(00B)",  "(11B)", "(11B)",  "(01B)", "(01B)", "(00C)", "(00C)",  "(11C)", "(11C)",  "(01C)", "(01C)", "(00D)", "(00D)",  "(11D)", "(11D)",  "(01D)", "(01D)")
        cols <- c("(11A)", "(01A)", "(00A)", "(01A)", "(00A)",  "(11A)",  "(11B)", "(01B)", "(00B)", "(01B)", "(00B)",  "(11B)",  "(11C)", "(01C)", "(00C)", "(01C)", "(00C)",  "(11C)",  "(11D)", "(01D)", "(00D)", "(01D)", "(00D)",  "(11D)")
        trans.tmp <- trans.rate[cbind(rows,cols)]
        trans.tmp[which(trans.tmp > 0)] = (trans.tmp[which(trans.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, trans.tmp)
        rows <- c(rep("(00A)", 3), rep("(11A)", 3), rep("(01A)", 3), rep("(00B)", 3), rep("(11B)", 3), rep("(01B)", 3), rep("(00C)", 3), rep("(11C)", 3), rep("(01C)", 3), rep("(00D)", 3), rep("(11D)", 3), rep("(01D)", 3))
        cols <- c("(00B)", "(00C)", "(00D)", "(11B)", "(11C)", "(11D)", "(01B)", "(01C)", "(01D)", "(00A)", "(00C)", "(00D)", "(11A)", "(11C)", "(11D)", "(01A)", "(01C)", "(01D)", "(00A)", "(00B)", "(00D)", "(11A)", "(11B)", "(11D)", "(01A)", "(01B)", "(01D)", "(00A)", "(00B)", "(00C)", "(11A)", "(11B)", "(11C)", "(01A)", "(01B)", "(01C)")
        category.tmp <- trans.rate[cbind(rows,cols)]
        category.rate.shift <- category.tmp + for.late.adjust
        category.rate.shift[is.na(category.rate.shift)] <- 0
        category.rate.shiftA <- c(category.rate.shift[1:3], rep(0,6), category.rate.shift[4:6], rep(0,6), category.rate.shift[7:9], rep(0,6))
        category.rate.shiftB <- c(category.rate.shift[10:12], rep(0,6), category.rate.shift[13:15], rep(0,6), category.rate.shift[16:18], rep(0,6))
        category.rate.shiftC <- c(category.rate.shift[19:21], rep(0,6), category.rate.shift[22:24], rep(0,6), category.rate.shift[25:27], rep(0,6))
        category.rate.shiftD <- c(category.rate.shift[28:30], rep(0,6), category.rate.shift[31:33], rep(0,6), category.rate.shift[34:36], rep(0,6))
        category.rates.all <- c(category.rate.shiftA, category.rate.shiftB, category.rate.shiftC, category.rate.shiftD)
        category.rates.unique <- length(unique(category.rates.all[category.rates.all>0]))
        pars.tmp <- c(turnover[1:3], extinct.frac.tmp[1:2], trans.tmp[1:6], category.rate.shiftA, turnover[4:6], extinct.frac.tmp[3:4], trans.tmp[7:12], category.rate.shiftB, turnover[7:9], extinct.frac.tmp[5:6], trans.tmp[13:18], category.rate.shiftC, turnover[10:12], extinct.frac.tmp[7:8], trans.tmp[19:24], category.rate.shiftD)
        pars[1:length(pars.tmp)] <- pars.tmp
    }
    
    if(dim(trans.rate)[2]==15){
        rate.cats <- 5
        pars.tmp <- turnover
        extinct.frac.tmp <- extinct.frac
        extinct.frac.tmp[which(extinct.frac.tmp > 0)] = (extinct.frac.tmp[which( extinct.frac.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, extinct.frac.tmp)
        for.late.adjust <- max(pars.tmp)
        rows <- c("(00A)", "(00A)",  "(11A)", "(11A)",  "(01A)", "(01A)", "(00B)", "(00B)",  "(11B)", "(11B)",  "(01B)", "(01B)", "(00C)", "(00C)",  "(11C)", "(11C)",  "(01C)", "(01C)", "(00D)", "(00D)",  "(11D)", "(11D)",  "(01D)", "(01D)", "(00E)", "(00E)",  "(11E)", "(11E)",  "(01E)", "(01E)")
        cols <- c("(11A)", "(01A)", "(00A)", "(01A)", "(00A)",  "(11A)",  "(11B)", "(01B)", "(00B)", "(01B)", "(00B)",  "(11B)",  "(11C)", "(01C)", "(00C)", "(01C)", "(00C)",  "(11C)",  "(11D)", "(01D)", "(00D)", "(01D)", "(00D)",  "(11D)",  "(11E)", "(01E)", "(00E)", "(01E)", "(00E)",  "(11E)")
        trans.tmp <- trans.rate[cbind(rows,cols)]
        trans.tmp[which(trans.tmp > 0)] = (trans.tmp[which(trans.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, trans.tmp)
        rows <- c(rep("(00A)", 4), rep("(11A)", 4), rep("(01A)", 4), rep("(00B)", 4), rep("(11B)", 4), rep("(01B)", 4), rep("(00C)", 4), rep("(11C)", 4), rep("(01C)", 4), rep("(00D)", 4), rep("(11D)", 4), rep("(01D)", 4), rep("(00E)", 4), rep("(11E)", 4), rep("(01E)", 4))
        cols <- c("(00B)", "(00C)", "(00D)", "(00E)", "(11B)", "(11C)", "(11D)", "(11E)", "(01B)", "(01C)", "(01D)", "(01E)", "(00A)", "(00C)", "(00D)", "(00E)", "(11A)", "(11C)", "(11D)", "(11E)", "(01A)", "(01C)", "(01D)", "(01E)", "(00A)", "(00B)", "(00D)", "(00E)", "(11A)", "(11B)", "(11D)", "(11E)", "(01A)", "(01B)", "(01D)", "(01E)", "(00A)", "(00B)", "(00C)", "(00E)", "(11A)", "(11B)", "(11C)", "(11E)", "(01A)", "(01B)", "(01C)", "(01E)", "(00A)", "(00B)", "(00C)", "(00D)", "(11A)", "(11B)", "(11C)", "(11D)", "(01A)", "(01B)", "(01C)", "(01D)")
        category.tmp <- trans.rate[cbind(rows,cols)]
        category.rate.shift <- category.tmp + for.late.adjust
        category.rate.shift[is.na(category.rate.shift)] <- 0
        category.rate.shiftA <- c(category.rate.shift[1:4], rep(0,5), category.rate.shift[5:8], rep(0,5), category.rate.shift[9:12], rep(0,5))
        category.rate.shiftB <- c(category.rate.shift[13:16], rep(0,5), category.rate.shift[17:20], rep(0,5), category.rate.shift[21:24], rep(0,5))
        category.rate.shiftC <- c(category.rate.shift[25:28], rep(0,5), category.rate.shift[29:32], rep(0,5), category.rate.shift[33:36], rep(0,5))
        category.rate.shiftD <- c(category.rate.shift[37:40], rep(0,5), category.rate.shift[41:44], rep(0,5), category.rate.shift[45:48], rep(0,5))
        category.rate.shiftE <- c(category.rate.shift[49:52], rep(0,5), category.rate.shift[53:56], rep(0,5), category.rate.shift[57:60], rep(0,5))
        category.rates.all <- c(category.rate.shiftA, category.rate.shiftB, category.rate.shiftC, category.rate.shiftD, category.rate.shiftE)
        category.rates.unique <- length(unique(category.rates.all[category.rates.all>0]))
        pars.tmp <- c(turnover[1:3], extinct.frac.tmp[1:2], trans.tmp[1:6], category.rate.shiftA, turnover[4:6], extinct.frac.tmp[3:4], trans.tmp[7:12], category.rate.shiftB, turnover[7:9], extinct.frac.tmp[5:6], trans.tmp[13:18], category.rate.shiftC, turnover[10:12], extinct.frac.tmp[7:8], trans.tmp[19:24], category.rate.shiftD, turnover[13:15], extinct.frac.tmp[9:10], trans.tmp[25:30], category.rate.shiftE)
        pars[1:length(pars.tmp)] <- pars.tmp
    }
    
    if(dim(trans.rate)[2]==18){
        rate.cats <- 6
        pars.tmp <- turnover
        extinct.frac.tmp <- extinct.frac
        extinct.frac.tmp[which(extinct.frac.tmp > 0)] = (extinct.frac.tmp[which( extinct.frac.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, extinct.frac.tmp)
        for.late.adjust <- max(pars.tmp)
        rows <- c("(00A)", "(00A)",  "(11A)", "(11A)",  "(01A)", "(01A)", "(00B)", "(00B)",  "(11B)", "(11B)",  "(01B)", "(01B)", "(00C)", "(00C)",  "(11C)", "(11C)",  "(01C)", "(01C)", "(00D)", "(00D)",  "(11D)", "(11D)",  "(01D)", "(01D)", "(00E)", "(00E)",  "(11E)", "(11E)",  "(01E)", "(01E)", "(00F)", "(00F)",  "(11F)", "(11F)",  "(01F)", "(01F)")
        cols <- c("(11A)", "(01A)", "(00A)", "(01A)", "(00A)",  "(11A)",  "(11B)", "(01B)", "(00B)", "(01B)", "(00B)",  "(11B)",  "(11C)", "(01C)", "(00C)", "(01C)", "(00C)",  "(11C)",  "(11D)", "(01D)", "(00D)", "(01D)", "(00D)",  "(11D)",  "(11E)", "(01E)", "(00E)", "(01E)", "(00E)",  "(11E)",  "(11F)", "(01F)", "(00F)", "(01F)", "(00F)",  "(11F)")
        trans.tmp <- trans.rate[cbind(rows,cols)]
        trans.tmp[which(trans.tmp > 0)] = (trans.tmp[which(trans.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, trans.tmp)
        rows <- c(rep("(00A)", 5), rep("(11A)", 5), rep("(01A)", 5), rep("(00B)", 5), rep("(11B)", 5), rep("(01B)", 5), rep("(00C)", 5), rep("(11C)", 5), rep("(01C)", 5), rep("(00D)", 5), rep("(11D)", 5), rep("(01D)", 5), rep("(00E)", 5), rep("(11E)", 5), rep("(01E)", 5), rep("(00F)", 5), rep("(11F)", 5), rep("(01F)", 5))
        cols <- c("(00B)", "(00C)", "(00D)", "(00E)", "(00F)", "(11B)", "(11C)", "(11D)", "(11E)", "(11F)", "(01B)", "(01C)", "(01D)", "(01E)", "(01F)", "(00A)", "(00C)", "(00D)", "(00E)", "(00F)", "(11A)", "(11C)", "(11D)", "(11E)", "(11F)", "(01A)", "(01C)", "(01D)", "(01E)", "(01F)", "(00A)", "(00B)", "(00D)", "(00E)", "(00F)", "(11A)", "(11B)", "(11D)", "(11E)", "(11F)", "(01A)", "(01B)", "(01D)", "(01E)", "(01F)", "(00A)", "(00B)", "(00C)", "(00E)", "(00F)", "(11A)", "(11B)", "(11C)", "(11E)", "(11F)", "(01A)", "(01B)", "(01C)", "(01E)", "(01F)", "(00A)", "(00B)", "(00C)", "(00D)", "(00F)", "(11A)", "(11B)", "(11C)", "(11D)", "(11F)", "(01A)", "(01B)", "(01C)", "(01D)", "(01F)", "(00A)", "(00B)", "(00C)", "(00D)", "(00E)", "(11A)", "(11B)", "(11C)", "(11D)", "(11E)", "(01A)", "(01B)", "(01C)", "(01D)", "(01E)")
        category.tmp <- trans.rate[cbind(rows,cols)]
        category.rate.shift <- category.tmp + for.late.adjust
        category.rate.shift[is.na(category.rate.shift)] <- 0
        category.rate.shiftA <- c(category.rate.shift[1:5], rep(0,4), category.rate.shift[6:10], rep(0,4), category.rate.shift[11:15], rep(0,4))
        category.rate.shiftB <- c(category.rate.shift[16:20], rep(0,4), category.rate.shift[21:25], rep(0,4), category.rate.shift[26:30], rep(0,4))
        category.rate.shiftC <- c(category.rate.shift[31:35], rep(0,4), category.rate.shift[36:40], rep(0,4), category.rate.shift[41:45], rep(0,4))
        category.rate.shiftD <- c(category.rate.shift[46:50], rep(0,4), category.rate.shift[51:55], rep(0,4), category.rate.shift[56:60], rep(0,4))
        category.rate.shiftE <- c(category.rate.shift[61:65], rep(0,4), category.rate.shift[66:70], rep(0,4), category.rate.shift[71:75], rep(0,4))
        category.rate.shiftF <- c(category.rate.shift[76:80], rep(0,4), category.rate.shift[81:85], rep(0,4), category.rate.shift[86:90], rep(0,4))
        category.rates.all <- c(category.rate.shiftA, category.rate.shiftB, category.rate.shiftC, category.rate.shiftD, category.rate.shiftE, category.rate.shiftF)
        category.rates.unique <- length(unique(category.rates.all[category.rates.all>0]))
        pars.tmp <- c(turnover[1:3], extinct.frac.tmp[1:2], trans.tmp[1:6], category.rate.shiftA, turnover[4:6], extinct.frac.tmp[3:4], trans.tmp[7:12], category.rate.shiftB, turnover[7:9], extinct.frac.tmp[5:6], trans.tmp[13:18], category.rate.shiftC, turnover[10:12], extinct.frac.tmp[7:8], trans.tmp[19:24], category.rate.shiftD, turnover[13:15], extinct.frac.tmp[9:10], trans.tmp[25:30], category.rate.shiftE, turnover[16:18], extinct.frac.tmp[11:12], trans.tmp[31:36], category.rate.shiftF)
        pars[1:length(pars.tmp)] <- pars.tmp
    }
    
    
    if(dim(trans.rate)[2]==21){
        rate.cats <- 7
        pars.tmp <- turnover
        extinct.frac.tmp <- extinct.frac
        extinct.frac.tmp[which(extinct.frac.tmp > 0)] = (extinct.frac.tmp[which( extinct.frac.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, extinct.frac.tmp)
        for.late.adjust <- max(pars.tmp)
        rows <- c("(00A)", "(00A)",  "(11A)", "(11A)",  "(01A)", "(01A)", "(00B)", "(00B)",  "(11B)", "(11B)",  "(01B)", "(01B)", "(00C)", "(00C)",  "(11C)", "(11C)",  "(01C)", "(01C)", "(00D)", "(00D)",  "(11D)", "(11D)",  "(01D)", "(01D)", "(00E)", "(00E)",  "(11E)", "(11E)",  "(01E)", "(01E)", "(00F)", "(00F)",  "(11F)", "(11F)",  "(01F)", "(01F)", "(00G)", "(00G)",  "(11G)", "(11G)",  "(01G)", "(01G)")
        cols <- c("(11A)", "(01A)", "(00A)", "(01A)", "(00A)",  "(11A)",  "(11B)", "(01B)", "(00B)", "(01B)", "(00B)",  "(11B)",  "(11C)", "(01C)", "(00C)", "(01C)", "(00C)",  "(11C)",  "(11D)", "(01D)", "(00D)", "(01D)", "(00D)",  "(11D)",  "(11E)", "(01E)", "(00E)", "(01E)", "(00E)",  "(11E)",  "(11F)", "(01F)", "(00F)", "(01F)", "(00F)",  "(11F)",  "(11G)", "(01G)", "(00G)", "(01G)", "(00G)",  "(11G)")
        trans.tmp <- trans.rate[cbind(rows,cols)]
        trans.tmp[which(trans.tmp > 0)] = (trans.tmp[which(trans.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, trans.tmp)
        rows <- c(rep("(00A)", 6), rep("(11A)", 6), rep("(01A)", 6), rep("(00B)", 6), rep("(11B)", 6), rep("(01B)", 6), rep("(00C)", 6), rep("(11C)", 6), rep("(01C)", 6), rep("(00D)", 6), rep("(11D)", 6), rep("(01D)", 6), rep("(00E)", 6), rep("(11E)", 6), rep("(01E)", 6), rep("(00F)", 6), rep("(11F)", 6), rep("(01F)", 6), rep("(00G)", 6), rep("(11G)", 6), rep("(01G)", 6))
        cols <- c("(00B)", "(00C)", "(00D)", "(00E)", "(00F)", "(00G)", "(11B)", "(11C)", "(11D)", "(11E)", "(11F)", "(11G)", "(01B)", "(01C)", "(01D)", "(01E)", "(01F)", "(01G)", "(00A)", "(00C)", "(00D)", "(00E)", "(00F)", "(00G)", "(11A)", "(11C)", "(11D)", "(11E)", "(11F)", "(11G)", "(01A)", "(01C)", "(01D)", "(01E)", "(01F)", "(01G)", "(00A)", "(00B)", "(00D)", "(00E)", "(00F)", "(00G)", "(11A)", "(11B)", "(11D)", "(11E)", "(11F)", "(11G)", "(01A)", "(01B)", "(01D)", "(01E)", "(01F)", "(01G)", "(00A)", "(00B)", "(00C)", "(00E)", "(00F)", "(00G)", "(11A)", "(11B)", "(11C)", "(11E)", "(11F)", "(11G)", "(01A)", "(01B)", "(01C)", "(01E)", "(01F)", "(01G)", "(00A)", "(00B)", "(00C)", "(00D)", "(00F)", "(00G)", "(11A)", "(11B)", "(11C)", "(11D)", "(11F)", "(11G)", "(01A)", "(01B)", "(01C)", "(01D)", "(01F)", "(01G)", "(00A)", "(00B)", "(00C)", "(00D)", "(00E)", "(00G)", "(11A)", "(11B)", "(11C)", "(11D)", "(11E)", "(11G)", "(01A)", "(01B)", "(01C)", "(01D)", "(01E)", "(01G)", "(00A)", "(00B)", "(00C)", "(00D)", "(00E)", "(00F)", "(11A)", "(11B)", "(11C)", "(11D)", "(11E)", "(11F)", "(01A)", "(01B)", "(01C)", "(01D)", "(01E)", "(01F)")
        category.tmp <- trans.rate[cbind(rows,cols)]
        category.rate.shift <- category.tmp + for.late.adjust
        category.rate.shift[is.na(category.rate.shift)] <- 0
        category.rate.shiftA <- c(category.rate.shift[1:6], rep(0,3), category.rate.shift[7:12], rep(0,3), category.rate.shift[13:18], rep(0,3))
        category.rate.shiftB <- c(category.rate.shift[19:24], rep(0,3), category.rate.shift[25:30], rep(0,3), category.rate.shift[31:36], rep(0,3))
        category.rate.shiftC <- c(category.rate.shift[37:42], rep(0,3), category.rate.shift[43:48], rep(0,3), category.rate.shift[49:54], rep(0,3))
        category.rate.shiftD <- c(category.rate.shift[55:60], rep(0,3), category.rate.shift[61:66], rep(0,3), category.rate.shift[67:72], rep(0,3))
        category.rate.shiftE <- c(category.rate.shift[73:78], rep(0,3), category.rate.shift[79:84], rep(0,3), category.rate.shift[85:90], rep(0,3))
        category.rate.shiftF <- c(category.rate.shift[91:96], rep(0,3), category.rate.shift[97:102], rep(0,3), category.rate.shift[103:108], rep(0,3))
        category.rate.shiftG <- c(category.rate.shift[109:114], rep(0,3), category.rate.shift[115:120], rep(0,3), category.rate.shift[121:126], rep(0,3))
        category.rates.all <- c(category.rate.shiftA, category.rate.shiftB, category.rate.shiftC, category.rate.shiftD, category.rate.shiftE, category.rate.shiftF, category.rate.shiftG)
        category.rates.unique <- length(unique(category.rates.all[category.rates.all>0]))
        pars.tmp <- c(turnover[1:3], extinct.frac.tmp[1:2], trans.tmp[1:6], category.rate.shiftA, turnover[4:6], extinct.frac.tmp[3:4], trans.tmp[7:12], category.rate.shiftB, turnover[7:9], extinct.frac.tmp[5:6], trans.tmp[13:18], category.rate.shiftC, turnover[10:12], extinct.frac.tmp[7:8], trans.tmp[19:24], category.rate.shiftD, turnover[13:15], extinct.frac.tmp[9:10], trans.tmp[25:30], category.rate.shiftE, turnover[16:18], extinct.frac.tmp[11:12], trans.tmp[31:36], category.rate.shiftF, turnover[19:21], extinct.frac.tmp[13:14], trans.tmp[37:42], category.rate.shiftG)
        pars[1:length(pars.tmp)] <- pars.tmp
    }
    
    
    if(dim(trans.rate)[2]==24){
        rate.cats <- 8
        pars.tmp <- turnover
        extinct.frac.tmp <- extinct.frac
        extinct.frac.tmp[which(extinct.frac.tmp > 0)] = (extinct.frac.tmp[which( extinct.frac.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, extinct.frac.tmp)
        for.late.adjust <- max(pars.tmp)
        rows <- c("(00A)", "(00A)",  "(11A)", "(11A)",  "(01A)", "(01A)", "(00B)", "(00B)",  "(11B)", "(11B)",  "(01B)", "(01B)", "(00C)", "(00C)",  "(11C)", "(11C)",  "(01C)", "(01C)", "(00D)", "(00D)",  "(11D)", "(11D)",  "(01D)", "(01D)", "(00E)", "(00E)",  "(11E)", "(11E)",  "(01E)", "(01E)", "(00F)", "(00F)",  "(11F)", "(11F)",  "(01F)", "(01F)", "(00G)", "(00G)",  "(11G)", "(11G)",  "(01G)", "(01G)", "(00H)", "(00H)",  "(11H)", "(11H)",  "(01H)", "(01H)")
        cols <- c("(11A)", "(01A)", "(00A)", "(01A)", "(00A)",  "(11A)",  "(11B)", "(01B)", "(00B)", "(01B)", "(00B)",  "(11B)",  "(11C)", "(01C)", "(00C)", "(01C)", "(00C)",  "(11C)",  "(11D)", "(01D)", "(00D)", "(01D)", "(00D)",  "(11D)",  "(11E)", "(01E)", "(00E)", "(01E)", "(00E)",  "(11E)",  "(11F)", "(01F)", "(00F)", "(01F)", "(00F)",  "(11F)",  "(11G)", "(01G)", "(00G)", "(01G)", "(00G)",  "(11G)",  "(11H)", "(01H)", "(00H)", "(01H)", "(00H)",  "(11H)")
        trans.tmp <- trans.rate[cbind(rows,cols)]
        trans.tmp[which(trans.tmp > 0)] = (trans.tmp[which(trans.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, trans.tmp)
        rows <- c(rep("(00A)", 7), rep("(11A)", 7), rep("(01A)", 7), rep("(00B)", 7), rep("(11B)", 7), rep("(01B)", 7), rep("(00C)", 7), rep("(11C)", 7), rep("(01C)", 7), rep("(00D)", 7), rep("(11D)", 7), rep("(01D)", 7), rep("(00E)", 7), rep("(11E)", 7), rep("(01E)", 7), rep("(00F)", 7), rep("(11F)", 7), rep("(01F)", 7), rep("(00G)", 7), rep("(11G)", 7), rep("(01G)", 7), rep("(00H)", 7), rep("(11H)", 7), rep("(01H)", 7))
        cols <- c("(00B)", "(00C)", "(00D)", "(00E)", "(00F)", "(00G)", "(00H)", "(11B)", "(11C)", "(11D)", "(11E)", "(11F)", "(11G)", "(11H)", "(01B)", "(01C)", "(01D)", "(01E)", "(01F)", "(01G)", "(01H)", "(00A)", "(00C)", "(00D)", "(00E)", "(00F)", "(00G)", "(00H)", "(11A)", "(11C)", "(11D)", "(11E)", "(11F)", "(11G)", "(11H)", "(01A)", "(01C)", "(01D)", "(01E)", "(01F)", "(01G)", "(01H)", "(00A)", "(00B)", "(00D)", "(00E)", "(00F)", "(00G)", "(00H)", "(11A)", "(11B)", "(11D)", "(11E)", "(11F)", "(11G)", "(11H)", "(01A)", "(01B)", "(01D)", "(01E)", "(01F)", "(01G)", "(01H)", "(00A)", "(00B)", "(00C)", "(00E)", "(00F)", "(00G)", "(00H)", "(11A)", "(11B)", "(11C)", "(11E)", "(11F)", "(11G)", "(11H)", "(01A)", "(01B)", "(01C)", "(01E)", "(01F)", "(01G)", "(01H)", "(00A)", "(00B)", "(00C)", "(00D)", "(00F)", "(00G)", "(00H)", "(11A)", "(11B)", "(11C)", "(11D)", "(11F)", "(11G)", "(11H)", "(01A)", "(01B)", "(01C)", "(01D)", "(01F)", "(01G)", "(01H)", "(00A)", "(00B)", "(00C)", "(00D)", "(00E)", "(00G)", "(00H)", "(11A)", "(11B)", "(11C)", "(11D)", "(11E)", "(11G)", "(11H)", "(01A)", "(01B)", "(01C)", "(01D)", "(01E)", "(01G)", "(01H)", "(00A)", "(00B)", "(00C)", "(00D)", "(00E)", "(00F)", "(00H)", "(11A)", "(11B)", "(11C)", "(11D)", "(11E)", "(11F)", "(11H)", "(01A)", "(01B)", "(01C)", "(01D)", "(01E)", "(01F)", "(01H)", "(00A)", "(00B)", "(00C)", "(00D)", "(00E)", "(00F)", "(00G)", "(11A)", "(11B)", "(11C)", "(11D)", "(11E)", "(11F)", "(11G)", "(01A)", "(01B)", "(01C)", "(01D)", "(01E)", "(01F)", "(01G)")
        category.tmp <- trans.rate[cbind(rows,cols)]
        category.rate.shift <- category.tmp + for.late.adjust
        category.rate.shift[is.na(category.rate.shift)] <- 0
        category.rate.shiftA <- c(category.rate.shift[1:7], rep(0,2), category.rate.shift[8:14], rep(0,2), category.rate.shift[15:21], rep(0,2))
        category.rate.shiftB <- c(category.rate.shift[22:28], rep(0,2), category.rate.shift[29:35], rep(0,2), category.rate.shift[36:42], rep(0,2))
        category.rate.shiftC <- c(category.rate.shift[43:49], rep(0,2), category.rate.shift[50:56], rep(0,2), category.rate.shift[57:63], rep(0,2))
        category.rate.shiftD <- c(category.rate.shift[64:70], rep(0,2), category.rate.shift[71:77], rep(0,2), category.rate.shift[78:84], rep(0,2))
        category.rate.shiftE <- c(category.rate.shift[85:91], rep(0,2), category.rate.shift[92:98], rep(0,2), category.rate.shift[99:105], rep(0,2))
        category.rate.shiftF <- c(category.rate.shift[106:112], rep(0,2), category.rate.shift[113:119], rep(0,2), category.rate.shift[120:126], rep(0,2))
        category.rate.shiftG <- c(category.rate.shift[127:133], rep(0,2), category.rate.shift[134:140], rep(0,2), category.rate.shift[141:147], rep(0,2))
        category.rate.shiftH <- c(category.rate.shift[148:154], rep(0,2), category.rate.shift[155:161], rep(0,2), category.rate.shift[162:168], rep(0,2))
        category.rates.all <- c(category.rate.shiftA, category.rate.shiftB, category.rate.shiftC, category.rate.shiftD, category.rate.shiftE, category.rate.shiftF, category.rate.shiftG, category.rate.shiftH)
        category.rates.unique <- length(unique(category.rates.all[category.rates.all>0]))
        pars.tmp <- c(turnover[1:3], extinct.frac.tmp[1:2], trans.tmp[1:6], category.rate.shiftA, turnover[4:6], extinct.frac.tmp[3:4], trans.tmp[7:12], category.rate.shiftB, turnover[7:9], extinct.frac.tmp[5:6], trans.tmp[13:18], category.rate.shiftC, turnover[10:12], extinct.frac.tmp[7:8], trans.tmp[19:24], category.rate.shiftD, turnover[13:15], extinct.frac.tmp[9:10], trans.tmp[25:30], category.rate.shiftE, turnover[16:18], extinct.frac.tmp[11:12], trans.tmp[31:36], category.rate.shiftF, turnover[19:21], extinct.frac.tmp[13:14], trans.tmp[37:42], category.rate.shiftG, turnover[22:24], extinct.frac.tmp[15:16], trans.tmp[43:48], category.rate.shiftH)
        pars[1:length(pars.tmp)] <- pars.tmp
    }
    
    
    if(dim(trans.rate)[2]==27){
        rate.cats <- 9
        pars.tmp <- turnover
        extinct.frac.tmp <- extinct.frac
        extinct.frac.tmp[which(extinct.frac.tmp > 0)] = (extinct.frac.tmp[which( extinct.frac.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, extinct.frac.tmp)
        for.late.adjust <- max(pars.tmp)
        rows <- c("(00A)", "(00A)",  "(11A)", "(11A)",  "(01A)", "(01A)", "(00B)", "(00B)",  "(11B)", "(11B)",  "(01B)", "(01B)", "(00C)", "(00C)",  "(11C)", "(11C)",  "(01C)", "(01C)", "(00D)", "(00D)",  "(11D)", "(11D)",  "(01D)", "(01D)", "(00E)", "(00E)",  "(11E)", "(11E)",  "(01E)", "(01E)", "(00F)", "(00F)",  "(11F)", "(11F)",  "(01F)", "(01F)", "(00G)", "(00G)",  "(11G)", "(11G)",  "(01G)", "(01G)", "(00H)", "(00H)",  "(11H)", "(11H)",  "(01H)", "(01H)", "(00I)", "(00I)",  "(11I)", "(11I)",  "(01I)", "(01I)")
        cols <- c("(11A)", "(01A)", "(00A)", "(01A)", "(00A)",  "(11A)",  "(11B)", "(01B)", "(00B)", "(01B)", "(00B)",  "(11B)",  "(11C)", "(01C)", "(00C)", "(01C)", "(00C)",  "(11C)",  "(11D)", "(01D)", "(00D)", "(01D)", "(00D)",  "(11D)",  "(11E)", "(01E)", "(00E)", "(01E)", "(00E)",  "(11E)",  "(11F)", "(01F)", "(00F)", "(01F)", "(00F)",  "(11F)",  "(11G)", "(01G)", "(00G)", "(01G)", "(00G)",  "(11G)",  "(11H)", "(01H)", "(00H)", "(01H)", "(00H)",  "(11H)",  "(11I)", "(01I)", "(00I)", "(01I)", "(00I)",  "(11I)")
        trans.tmp <- trans.rate[cbind(rows,cols)]
        trans.tmp[which(trans.tmp > 0)] = (trans.tmp[which(trans.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, trans.tmp)
        rows <- c(rep("(00A)", 8), rep("(11A)", 8), rep("(01A)", 8), rep("(00B)", 8), rep("(11B)", 8), rep("(01B)", 8), rep("(00C)", 8), rep("(11C)", 8), rep("(01C)", 8), rep("(00D)", 8), rep("(11D)", 8), rep("(01D)", 8), rep("(00E)", 8), rep("(11E)", 8), rep("(01E)", 8), rep("(00F)", 8), rep("(11F)", 8), rep("(01F)", 8), rep("(00G)", 8), rep("(11G)", 8), rep("(01G)", 8), rep("(00H)", 8), rep("(11H)", 8), rep("(01H)", 8), rep("(00I)", 8), rep("(11I)", 8), rep("(01I)", 8))
        cols <- c("(00B)", "(00C)", "(00D)", "(00E)", "(00F)", "(00G)", "(00H)", "(00I)", "(11B)", "(11C)", "(11D)", "(11E)", "(11F)", "(11G)", "(11H)", "(11I)", "(01B)", "(01C)", "(01D)", "(01E)", "(01F)", "(01G)", "(01H)", "(01I)", "(00A)", "(00C)", "(00D)", "(00E)", "(00F)", "(00G)", "(00H)", "(00I)", "(11A)", "(11C)", "(11D)", "(11E)", "(11F)", "(11G)", "(11H)", "(11I)", "(01A)", "(01C)", "(01D)", "(01E)", "(01F)", "(01G)", "(01H)", "(01I)", "(00A)", "(00B)", "(00D)", "(00E)", "(00F)", "(00G)", "(00H)", "(00I)", "(11A)", "(11B)", "(11D)", "(11E)", "(11F)", "(11G)", "(11H)", "(11I)", "(01A)", "(01B)", "(01D)", "(01E)", "(01F)", "(01G)", "(01H)", "(01I)", "(00A)", "(00B)", "(00C)", "(00E)", "(00F)", "(00G)", "(00H)", "(00I)", "(11A)", "(11B)", "(11C)", "(11E)", "(11F)", "(11G)", "(11H)", "(11I)", "(01A)", "(01B)", "(01C)", "(01E)", "(01F)", "(01G)", "(01H)", "(01I)", "(00A)", "(00B)", "(00C)", "(00D)", "(00F)", "(00G)", "(00H)", "(00I)", "(11A)", "(11B)", "(11C)", "(11D)", "(11F)", "(11G)", "(11H)", "(11I)", "(01A)", "(01B)", "(01C)", "(01D)", "(01F)", "(01G)", "(01H)", "(01I)", "(00A)", "(00B)", "(00C)", "(00D)", "(00E)", "(00G)", "(00H)", "(00I)", "(11A)", "(11B)", "(11C)", "(11D)", "(11E)", "(11G)", "(11H)", "(11I)", "(01A)", "(01B)", "(01C)", "(01D)", "(01E)", "(01G)", "(01H)", "(01I)", "(00A)", "(00B)", "(00C)", "(00D)", "(00E)", "(00F)", "(00H)", "(00I)", "(11A)", "(11B)", "(11C)", "(11D)", "(11E)", "(11F)", "(11H)", "(11I)", "(01A)", "(01B)", "(01C)", "(01D)", "(01E)", "(01F)", "(01H)", "(01I)", "(00A)", "(00B)", "(00C)", "(00D)", "(00E)", "(00F)", "(00G)", "(00I)", "(11A)", "(11B)", "(11C)", "(11D)", "(11E)", "(11F)", "(11G)", "(11I)", "(01A)", "(01B)", "(01C)", "(01D)", "(01E)", "(01F)", "(01G)", "(01I)", "(00A)", "(00B)", "(00C)", "(00D)", "(00E)", "(00F)", "(00G)", "(00H)", "(11A)", "(11B)", "(11C)", "(11D)", "(11E)", "(11F)", "(11G)", "(11H)", "(01A)", "(01B)", "(01C)", "(01D)", "(01E)", "(01F)", "(01G)", "(01H)")
        category.tmp <- trans.rate[cbind(rows,cols)]
        category.rate.shift <- category.tmp + for.late.adjust
        category.rate.shift[is.na(category.rate.shift)] <- 0
        category.rate.shiftA <- c(category.rate.shift[1:8], rep(0,1), category.rate.shift[9:16], rep(0,1), category.rate.shift[17:24], rep(0,1))
        category.rate.shiftB <- c(category.rate.shift[25:32], rep(0,1), category.rate.shift[33:40], rep(0,1), category.rate.shift[41:48], rep(0,1))
        category.rate.shiftC <- c(category.rate.shift[49:56], rep(0,1), category.rate.shift[57:64], rep(0,1), category.rate.shift[65:72], rep(0,1))
        category.rate.shiftD <- c(category.rate.shift[73:80], rep(0,1), category.rate.shift[81:88], rep(0,1), category.rate.shift[89:96], rep(0,1))
        category.rate.shiftE <- c(category.rate.shift[97:104], rep(0,1), category.rate.shift[105:112], rep(0,1), category.rate.shift[113:120], rep(0,1))
        category.rate.shiftF <- c(category.rate.shift[121:128], rep(0,1), category.rate.shift[129:136], rep(0,1), category.rate.shift[137:144], rep(0,1))
        category.rate.shiftG <- c(category.rate.shift[145:152], rep(0,1), category.rate.shift[153:160], rep(0,1), category.rate.shift[161:168], rep(0,1))
        category.rate.shiftH <- c(category.rate.shift[169:176], rep(0,1), category.rate.shift[177:184], rep(0,1), category.rate.shift[185:192], rep(0,1))
        category.rate.shiftI <- c(category.rate.shift[193:200], rep(0,1), category.rate.shift[201:208], rep(0,1), category.rate.shift[209:216], rep(0,1))
        category.rates.all <- c(category.rate.shiftA, category.rate.shiftB, category.rate.shiftC, category.rate.shiftD, category.rate.shiftE, category.rate.shiftF, category.rate.shiftG, category.rate.shiftH, category.rate.shiftI)
        category.rates.unique <- length(unique(category.rates.all[category.rates.all>0]))
        pars.tmp <- c(turnover[1:3], extinct.frac.tmp[1:2], trans.tmp[1:6], category.rate.shiftA, turnover[4:6], extinct.frac.tmp[3:4], trans.tmp[7:12], category.rate.shiftB, turnover[7:9], extinct.frac.tmp[5:6], trans.tmp[13:18], category.rate.shiftC, turnover[10:12], extinct.frac.tmp[7:8], trans.tmp[19:24], category.rate.shiftD, turnover[13:15], extinct.frac.tmp[9:10], trans.tmp[25:30], category.rate.shiftE, turnover[16:18], extinct.frac.tmp[11:12], trans.tmp[31:36], category.rate.shiftF, turnover[19:21], extinct.frac.tmp[13:14], trans.tmp[37:42], category.rate.shiftG, turnover[22:24], extinct.frac.tmp[15:16], trans.tmp[43:48], category.rate.shiftH, turnover[25:27], extinct.frac.tmp[17:18], trans.tmp[49:54], category.rate.shiftI)
        pars[1:length(pars.tmp)] <- pars.tmp
    }
    
    if(dim(trans.rate)[2]==30){
        rate.cats <- 10
        pars.tmp <- turnover
        extinct.frac.tmp <- extinct.frac
        extinct.frac.tmp[which(extinct.frac.tmp > 0)] = (extinct.frac.tmp[which( extinct.frac.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, extinct.frac.tmp)
        for.late.adjust <- max(pars.tmp)
        rows <- c("(00A)", "(00A)",  "(11A)", "(11A)",  "(01A)", "(01A)", "(00B)", "(00B)",  "(11B)", "(11B)",  "(01B)", "(01B)", "(00C)", "(00C)",  "(11C)", "(11C)",  "(01C)", "(01C)", "(00D)", "(00D)",  "(11D)", "(11D)",  "(01D)", "(01D)", "(00E)", "(00E)",  "(11E)", "(11E)",  "(01E)", "(01E)", "(00F)", "(00F)",  "(11F)", "(11F)",  "(01F)", "(01F)", "(00G)", "(00G)",  "(11G)", "(11G)",  "(01G)", "(01G)", "(00H)", "(00H)",  "(11H)", "(11H)",  "(01H)", "(01H)", "(00I)", "(00I)",  "(11I)", "(11I)",  "(01I)", "(01I)", "(00J)", "(00J)",  "(11J)", "(11J)",  "(01J)", "(01J)")
        cols <- c("(11A)", "(01A)", "(00A)", "(01A)", "(00A)",  "(11A)",  "(11B)", "(01B)", "(00B)", "(01B)", "(00B)",  "(11B)",  "(11C)", "(01C)", "(00C)", "(01C)", "(00C)",  "(11C)",  "(11D)", "(01D)", "(00D)", "(01D)", "(00D)",  "(11D)",  "(11E)", "(01E)", "(00E)", "(01E)", "(00E)",  "(11E)",  "(11F)", "(01F)", "(00F)", "(01F)", "(00F)",  "(11F)",  "(11G)", "(01G)", "(00G)", "(01G)", "(00G)",  "(11G)",  "(11H)", "(01H)", "(00H)", "(01H)", "(00H)",  "(11H)",  "(11I)", "(01I)", "(00I)", "(01I)", "(00I)",  "(11I)",  "(11J)", "(01J)", "(00J)", "(01J)", "(00J)", "(11J)")
        trans.tmp <- trans.rate[cbind(rows,cols)]
        trans.tmp[which(trans.tmp > 0)] = (trans.tmp[which(trans.tmp > 0)] + max(pars.tmp))
        pars.tmp <- c(pars.tmp, trans.tmp)
        rows <- c(rep("(00A)", 9), rep("(11A)", 9), rep("(01A)", 9), rep("(00B)", 9), rep("(11B)", 9), rep("(01B)", 9), rep("(00C)", 9), rep("(11C)", 9), rep("(01C)", 9), rep("(00D)", 9), rep("(11D)", 9), rep("(01D)", 9), rep("(00E)", 9), rep("(11E)", 9), rep("(01E)", 9), rep("(00F)", 9), rep("(11F)", 9), rep("(01F)", 9), rep("(00G)", 9), rep("(11G)", 9), rep("(01G)", 9), rep("(00H)", 9), rep("(11H)", 9), rep("(01H)", 9), rep("(00I)", 9), rep("(11I)", 9), rep("(01I)", 9), rep("(00J)", 9), rep("(11J)", 9), rep("(01J)", 9))
        cols <- c("(00B)", "(00C)", "(00D)", "(00E)", "(00F)", "(00G)", "(00H)", "(00I)", "(00J)", "(11B)", "(11C)", "(11D)", "(11E)", "(11F)", "(11G)", "(11H)", "(11I)", "(11J)", "(01B)", "(01C)", "(01D)", "(01E)", "(01F)", "(01G)", "(01H)", "(01I)", "(01J)", "(00A)", "(00C)", "(00D)", "(00E)", "(00F)", "(00G)", "(00H)", "(00I)", "(00J)", "(11A)", "(11C)", "(11D)", "(11E)", "(11F)", "(11G)", "(11H)", "(11I)", "(11J)", "(01A)", "(01C)", "(01D)", "(01E)", "(01F)", "(01G)", "(01H)", "(01I)", "(01J)", "(00A)", "(00B)", "(00D)", "(00E)", "(00F)", "(00G)", "(00H)", "(00I)", "(00J)", "(11A)", "(11B)", "(11D)", "(11E)", "(11F)", "(11G)", "(11H)", "(11I)", "(11J)", "(01A)", "(01B)", "(01D)", "(01E)", "(01F)", "(01G)", "(01H)", "(01I)", "(01J)", "(00A)", "(00B)", "(00C)", "(00E)", "(00F)", "(00G)", "(00H)", "(00I)", "(00J)", "(11A)", "(11B)", "(11C)", "(11E)", "(11F)", "(11G)", "(11H)", "(11I)", "(11J)", "(01A)", "(01B)", "(01C)", "(01E)", "(01F)", "(01G)", "(01H)", "(01I)", "(01J)", "(00A)", "(00B)", "(00C)", "(00D)", "(00F)", "(00G)", "(00H)", "(00I)", "(00J)", "(11A)", "(11B)", "(11C)", "(11D)", "(11F)", "(11G)", "(11H)", "(11I)", "(11J)", "(01A)", "(01B)", "(01C)", "(01D)", "(01F)", "(01G)", "(01H)", "(01I)", "(01J)", "(00A)", "(00B)", "(00C)", "(00D)", "(00E)", "(00G)", "(00H)", "(00I)", "(00J)", "(11A)", "(11B)", "(11C)", "(11D)", "(11E)", "(11G)", "(11H)", "(11I)", "(11J)", "(01A)", "(01B)", "(01C)", "(01D)", "(01E)", "(01G)", "(01H)", "(01I)", "(01J)", "(00A)", "(00B)", "(00C)", "(00D)", "(00E)", "(00F)", "(00H)", "(00I)", "(00J)", "(11A)", "(11B)", "(11C)", "(11D)", "(11E)", "(11F)", "(11H)", "(11I)", "(11J)", "(01A)", "(01B)", "(01C)", "(01D)", "(01E)", "(01F)", "(01H)", "(01I)", "(01J)", "(00A)", "(00B)", "(00C)", "(00D)", "(00E)", "(00F)", "(00G)", "(00I)", "(00J)", "(11A)", "(11B)", "(11C)", "(11D)", "(11E)", "(11F)", "(11G)", "(11I)", "(11J)", "(01A)", "(01B)", "(01C)", "(01D)", "(01E)", "(01F)", "(01G)", "(01I)", "(01J)",  "(00A)", "(00B)", "(00C)", "(00D)", "(00E)", "(00F)", "(00G)", "(00H)", "(00J)", "(11A)", "(11B)", "(11C)", "(11D)", "(11E)", "(11F)", "(11G)", "(11H)", "(11J)", "(01A)", "(01B)", "(01C)", "(01D)", "(01E)", "(01F)", "(01G)", "(01H)", "(01J)", "(00A)", "(00B)", "(00C)", "(00D)", "(00E)", "(00F)", "(00G)", "(00H)", "(00I)", "(11A)", "(11B)", "(11C)", "(11D)", "(11E)", "(11F)", "(11G)", "(11H)", "(11I)", "(01A)", "(01B)", "(01C)", "(01D)", "(01E)", "(01F)", "(01G)", "(01H)", "(01I)")
        category.tmp <- trans.rate[cbind(rows,cols)]
        category.rate.shift <- category.tmp + for.late.adjust
        category.rate.shift[is.na(category.rate.shift)] <- 0
        category.rate.shiftA <- category.rate.shift[1:27]
        category.rate.shiftB <- category.rate.shift[28:54]
        category.rate.shiftC <- category.rate.shift[55:81]
        category.rate.shiftD <- category.rate.shift[82:108]
        category.rate.shiftE <- category.rate.shift[109:135]
        category.rate.shiftF <- category.rate.shift[136:162]
        category.rate.shiftG <- category.rate.shift[163:189]
        category.rate.shiftH <- category.rate.shift[190:216]
        category.rate.shiftI <- category.rate.shift[217:243]
        category.rate.shiftJ <- category.rate.shift[244:270]
        category.rates.all <- c(category.rate.shiftA, category.rate.shiftB, category.rate.shiftC, category.rate.shiftD, category.rate.shiftE, category.rate.shiftF, category.rate.shiftG, category.rate.shiftH, category.rate.shiftI, category.rate.shiftJ)
        category.rates.unique <- length(unique(category.rates.all[category.rates.all>0]))
        pars.tmp <- c(turnover[1:3], extinct.frac.tmp[1:2], trans.tmp[1:6], category.rate.shiftA, turnover[4:6], extinct.frac.tmp[3:4], trans.tmp[7:12], category.rate.shiftB, turnover[7:9], extinct.frac.tmp[5:6], trans.tmp[13:18], category.rate.shiftC, turnover[10:12], extinct.frac.tmp[7:8], trans.tmp[19:24], category.rate.shiftD, turnover[13:15], extinct.frac.tmp[9:10], trans.tmp[25:30], category.rate.shiftE, turnover[16:18], extinct.frac.tmp[11:12], trans.tmp[31:36], category.rate.shiftF, turnover[19:21], extinct.frac.tmp[13:14], trans.tmp[37:42], category.rate.shiftG, turnover[22:24], extinct.frac.tmp[15:16], trans.tmp[43:48], category.rate.shiftH, turnover[25:27], extinct.frac.tmp[17:18], trans.tmp[49:54], category.rate.shiftI, turnover[28:30], extinct.frac.tmp[19:20], trans.tmp[55:60], category.rate.shiftJ)
        pars[1:length(pars.tmp)] <- pars.tmp
    }

    np <- max(pars)
    pars[pars==0] <- np+1
    
    cat("Initializing...", "\n")
    
    data.new <- data.frame(data[,2], data[,2], row.names=data[,1])
    data.new <- data.new[phy$tip.label,]
    
    #This is used to scale starting values to account for sampling:
    if(length(f) == 3){
        freqs <- table(data.new[,1])
        if(length(freqs == 2)){
            samp.freq.tree <- Ntip(phy) / sum(table(data.new[,1]) / f[as.numeric(names(freqs))+1])
        }else{
            samp.freq.tree <- Ntip(phy) / sum(table(data.new[,1]) / f[as.numeric(names(freqs))+1])
        }
    }else{
        if(length(f) == Ntip(phy)){
            stop("This is functionality has been temporarily removed.")
            #samp.freq.tree <- Ntip(phy) / sum(table(data.new[,1]) / mean(f[as.numeric(names(freqs))+1]))
        }else{
            stop("The vector of sampling frequencies does not match the number of tips in the tree.")
        }
    }
    
    if(is.null(restart.obj)){
        if(sum(extinct.frac)==0){
            init.pars <- starting.point.geosse(phy, eps=0, samp.freq.tree=samp.freq.tree)
        }else{
            #init.pars <- starting.point.geosse(phy, eps=mag.san.start, samp.freq.tree=samp.freq.tree)
            init.pars <- starting.point.generator(phy, k=3, samp.freq.tree=samp.freq.tree, q.div=5, yule=FALSE)
        }
        names(init.pars) <- NULL
        
        if(is.null(starting.vals)){
            #def.set.pars <- rep(c(log(init.pars[1:3]), log(init.pars[4:5]), rep(log(init.pars[6:7]*.1),3), rep(log(.01), 27)), rate.cats)
            #def.set.pars <- rep(c(log(init.pars[1:3]), log(init.pars[4:5]), rep(log(init.pars[7:8]*.1),3), rep(log(.01), 27)), rate.cats)
            def.set.pars <- rep(c(log(init.pars[1]+init.pars[4]), log(init.pars[2]+init.pars[5]), log(sum(init.pars[1:3])), log(init.pars[4]/init.pars[1]),  log(init.pars[5]/init.pars[2]), rep(log(init.pars[7:8]),3), rep(log(0.001), 27)), rate.cats)
        }else{
            def.set.pars <- rep(c(log(starting.vals[1:3]), log(starting.vals[4:5]), rep(log(starting.vals[6:7]),3), rep(log(0.001), 27)), rate.cats)
            #def.set.pars <- rep(c(log(starting.vals[1:3]), log(starting.vals[4:5]), rep(log(init.pars[7:8]),3), rep(log(0.01), 27)), rate.cats)
        }
        if(bounded.search == TRUE){
            #upper.full <- rep(c(rep(log(turnover.upper),3), rep(log(extinct.frac.upper),2), rep(log(trans.upper),2*3), rep(log(10), 27)), rate.cats)
            upper.full <- rep(c(rep(log(10000),3), rep(log(3),2), rep(log(trans.upper), 2*3), rep(log(10), 27)), rate.cats)
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
    dat.tab <- OrganizeDataGeo(data=data.new[,1], phy=phy, f=f, hidden.states=hidden.areas)
    nb.tip <- Ntip(phy)
    nb.node <- phy$Nnode
    ##########################
    
    if(sann == FALSE){
        if(bounded.search == TRUE){
            cat("Finished. Beginning bounded subplex routine...", "\n")
            opts <- list("algorithm" = "NLOPT_LN_SBPLX", "maxeval" = 100000, "ftol_rel" = max.tol)
            out = nloptr(x0=ip, eval_f=DevOptimizeGeoHiSSEfast, ub=upper, lb=lower, opts=opts, pars=pars, dat.tab=dat.tab, gen=gen, hidden.states=hidden.areas, assume.cladogenetic=assume.cladogenetic, nb.tip=nb.tip, nb.node=nb.node, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, np=np, ode.eps=ode.eps)
            solution <- numeric(length(pars))
            solution[] <- c(exp(out$solution), 0)[pars]
            loglik = -out$objective
        }else{
            cat("Finished. Beginning subplex routine...", "\n")
            out = subplex(ip, fn=DevOptimizeGeoHiSSEfast, control=list(reltol=max.tol, parscale=rep(0.1, length(ip))), pars=pars, dat.tab=dat.tab, gen=gen, hidden.states=hidden.areas, assume.cladogenetic=assume.cladogenetic, nb.tip=nb.tip, nb.node=nb.node, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, np=np, ode.eps=ode.eps)
            solution <- numeric(length(pars))
            solution[] <- c(exp(out$par), 0)[pars]
            loglik = -out$value
        }
    }else{
        cat("Finished. Beginning simulated annealing...", "\n")
        out.sann = GenSA(ip, fn=DevOptimizeGeoHiSSEfast, lower=lower, upper=upper, control=list(max.call=sann.its), pars=pars, dat.tab=dat.tab, gen=gen, hidden.states=hidden.areas, assume.cladogenetic=assume.cladogenetic, nb.tip=nb.tip, nb.node=nb.node, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, np=np, ode.eps=ode.eps)
        cat("Finished. Refining using subplex routine...", "\n")
        out = nloptr(x0=out.sann$par, eval_f=DevOptimizeGeoHiSSE, ub=upper, lb=lower, opts=opts, pars=pars, dat.tab=dat.tab, gen=gen, hidden.states=hidden.areas, assume.cladogenetic=assume.cladogenetic, nb.tip=nb.tip, nb.node=nb.node, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, np=np, ode.eps=ode.eps)
        solution <- numeric(length(pars))
        solution[] <- c(exp(out$solution), 0)[pars]
        
        loglik = -out$objective
    }
    
    names(solution) <- c("tau00A","tau11A","tau01A","ef00A","ef11A","d00A_11A","d00A_01A","d11A_00A","d11A_01A","d01A_00A","d01A_11A","d00A_00B","d00A_00C","d00A_00D","d00A_00E","d00A_00F","d00A_00G","d00A_00H","d00A_00I","d00A_00J","d11A_11B","d11A_11C","d11A_11D","d11A_11E","d11A_11F","d11A_11G","d11A_11H","d11A_11I","d11A_11J","d01A_01B","d01A_01C","d01A_01D","d01A_01E","d01A_01F","d01A_01G","d01A_01H","d01A_01I","d01A_01J","tau00B","tau11B","tau01B","ef00B","ef11B","d00B_11B","d00B_01B","d11B_00B","d11B_01B","d01B_00B","d01B_11B","d00B_00A","d00B_00C","d00B_00D","d00B_00E","d00B_00F","d00B_00G","d00B_00H","d00B_00I","d00B_00J","d11B_11A","d11B_11C","d11B_11D","d11B_11E","d11B_11F","d11B_11G","d11B_11H","d11B_11I","d11B_11J","d01B_01A","d01B_01C","d01B_01D","d01B_01E","d01B_01F","d01B_01G","d01B_01H","d01B_01I","d01B_01J","tau00C","tau11C","tau01C","ef00C","ef11C","d00C_11C","d00C_01C","d11C_00C","d11C_01C","d01C_00C","d01C_11C","d00C_00A","d00C_00B","d00C_00D","d00C_00E","d00C_00F","d00C_00G","d00C_00H","d00C_00I","d00C_00J","d11C_11A","d11C_11B","d11C_11D","d11C_11E","d11C_11F","d11C_11G","d11C_11H","d11C_11I","d11C_11J","d01C_01A","d01C_01B","d01C_01D","d01C_01E","d01C_01F","d01C_01G","d01C_01H","d01C_01I","d01C_01J","tau00D","tau11D","tau01D","ef00D","ef11D","d00D_11D","d00D_01D","d11D_00D","d11D_01D","d01D_00D","d01D_11D","d00D_00A","d00D_00B","d00D_00C","d00D_00E","d00D_00F","d00D_00G","d00D_00H","d00D_00I","d00D_00J","d11D_11A","d11D_11B","d11D_11C","d11D_11E","d11D_11F","d11D_11G","d11D_11H","d11D_11I","d11D_11J","d01D_01A","d01D_01B","d01D_01C","d01D_01E","d01D_01F","d01D_01G","d01D_01H","d01D_01I","d01D_01J","tau00E","tau11E","tau01E","ef00E","ef11E","d00E_11E","d00E_01E","d11E_00E","d11E_01E","d01E_00E","d01E_11E","d00E_00A","d00E_00B","d00E_00C","d00E_00D","d00E_00F","d00E_00G","d00E_00H","d00E_00I","d00E_00J","d11E_11A","d11E_11B","d11E_11C","d11E_11D","d11E_11F","d11E_11G","d11E_11H","d11E_11I","d11E_11J","d01E_01A","d01E_01B","d01E_01C","d01E_01D","d01E_01F","d01E_01G","d01E_01H","d01E_01I","d01E_01J","tau00F","tau11F","tau01F","ef00F","ef11F","d00F_11F","d00F_01F","d11F_00F","d11F_01F","d01F_00F","d01F_11F","d00F_00A","d00F_00B","d00F_00C","d00F_00D","d00F_00E","d00F_00G","d00F_00H","d00F_00I","d00F_00J","d11F_11A","d11F_11B","d11F_11C","d11F_11D","d11F_11E","d11F_11G","d11F_11H","d11F_11I","d11F_11J","d01F_01A","d01F_01B","d01F_01C","d01F_01D","d01F_01E","d01F_01G","d01F_01H","d01F_01I","d01F_01J","tau00G","tau11G","tau01G","ef00G","ef11G","d00G_11G","d00G_01G","d11G_00G","d11G_01G","d01G_00G","d01G_11G","d00G_00A","d00G_00B","d00G_00C","d00G_00D","d00G_00E","d00G_00F","d00G_00H","d00G_00I","d00G_00J","d11G_11A","d11G_11B","d11G_11C","d11G_11D","d11G_11E","d11G_11F","d11G_11H","d11G_11I","d11G_11J","d01G_01A","d01G_01B","d01G_01C","d01G_01D","d01G_01E","d01G_01F","d01G_01H","d01G_01I","d01G_01J","tau00H","tau11H","tau01H","ef00H","ef11H","d00H_11H","d00H_01H","d11H_00H","d11H_01H","d01H_00H","d01H_11H","d00H_00A","d00H_00B","d00H_00C","d00H_00D","d00H_00E","d00H_00F","d00H_00G","d00H_00I","d00H_00J","d11H_11A","d11H_11B","d11H_11C","d11H_11D","d11H_11E","d11H_11F","d11H_11G","d11H_11I","d11H_11J","d01H_01A","d01H_01B","d01H_01C","d01H_01D","d01H_01E","d01H_01F","d01H_01G","d01H_01I","d01H_01J","tau00I","tau11I","tau01I","ef00I","ef11I","d00I_11I","d00I_01I","d11I_00I","d11I_01I","d01I_00I","d01I_11I","d00I_00A","d00I_00B","d00I_00C","d00I_00D","d00I_00E","d00I_00F","d00I_00G","d00I_00H","d00I_00J","d11I_11A","d11I_11B","d11I_11C","d11I_11D","d11I_11E","d11I_11F","d11I_11G","d11I_11H","d11I_11J","d01I_01A","d01I_01B","d01I_01C","d01I_01D","d01I_01E","d01I_01F","d01I_01G","d01I_01H","d01I_01J","tau00J","tau11J","tau01J","ef00J","ef11J","d00J_11J","d00J_01J","d11J_00J","d11J_01J","d01J_00J","d01J_11J","d00J_00A","d00J_00B","d00J_00C","d00J_00D","d00J_00E","d00J_00F","d00J_00G","d00J_00H","d00J_00I","d11J_11A","d11J_11B","d11J_11C","d11J_11D","d11J_11E","d11J_11F","d11J_11G","d11J_11H","d11J_11I","d01J_01A","d01J_01B","d01J_01C","d01J_01D","d01J_01E","d01J_01F","d01J_01G","d01J_01H","d01J_01I")

    cat("Finished. Summarizing results...", "\n")
    
    obj = list(loglik = loglik, AIC = -2*loglik+2*np, AICc = -2*loglik+(2*np*(Ntip(phy)/(Ntip(phy)-np-1))), solution=solution, index.par=pars, f=f, hidden.areas=hidden.areas, assume.cladogenetic=assume.cladogenetic, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p, phy=phy, data=data, trans.matrix=trans.rate, max.tol=max.tol, starting.vals=ip, upper.bounds=upper, lower.bounds=lower, ode.eps=ode.eps)
    class(obj) <- append(class(obj), "geohisse.fit")
    return(obj)
}



######################################################################################################################################
######################################################################################################################################
### The function used to optimize parameters:
######################################################################################################################################
######################################################################################################################################


DevOptimizeGeoHiSSEfast <- function(p, pars, dat.tab, gen, hidden.states, assume.cladogenetic=assume.cladogenetic, nb.tip=nb.tip, nb.node=nb.node, condition.on.survival, root.type, root.p, np, ode.eps) {
    #Generates the final vector with the appropriate parameter estimates in the right place:
    p.new <- exp(p)
    #print(p.new)
    model.vec <- numeric(length(pars))
    model.vec[] <- c(p.new, 0)[pars]
    cache = ParametersToPassGeoHiSSEfast(model.vec=model.vec, hidden.states=hidden.states, assume.cladogenetic=assume.cladogenetic, nb.tip=nb.tip, nb.node=nb.node, bad.likelihood=exp(-500), ode.eps=ode.eps)
    #save(cache, file="cache.Rsave")
    if(any(cache$s01A<0, cache$s01B<0, cache$s01C<0, cache$s01D<0, cache$s01E<0, cache$s01F<0, cache$s01G<0, cache$s01H<0, cache$s01I<0, cache$s01J<0)){
        #print("bad")
        #print(c(cache$s01A, cache$s01B, cache$s01C))
        return(-log(cache$bad.likelihood)^13)
    }else{
        logl <- DownPassGeoHissefast(dat.tab=dat.tab, cache=cache, gen=gen, condition.on.survival=condition.on.survival, root.type=root.type, root.p=root.p)
        #print(logl)
        return(-logl)
    }
}



######################################################################################################################################
######################################################################################################################################
### The various utility functions used
######################################################################################################################################
######################################################################################################################################

OrganizeDataGeo <- function(data, phy, f, hidden.states){
    ### Ughy McUgherson. This is a must in order to pass CRAN checks: http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
    DesNode = NULL
    
    nb.tip <- length(phy$tip.label)
    nb.node <- phy$Nnode
    
    if(hidden.states == FALSE){
        states = matrix(0,Ntip(phy),3)
        for(i in 1:Ntip(phy)){
            if(data[i]==1){states[i,1]=1}
            if(data[i]==2){states[i,2]=1}
            if(data[i]==0){states[i,3]=1}
        }
        compD <- matrix(0, nrow=nb.tip, ncol=3)
        compE <- matrix(0, nrow=nb.tip, ncol=3)
    }
    if(hidden.states == TRUE){
        states = matrix(0,Ntip(phy),30)
        for(i in 1:Ntip(phy)){
            if(data[i]==1){states[i,c(1,4,7,10,13,16,19,22,25,28)]=1}
            if(data[i]==2){states[i,c(2,5,8,11,14,17,20,23,26,29)]=1}
            if(data[i]==0){states[i,c(3,6,9,12,15,18,21,24,27,30)]=1}
        }
        compD <- matrix(0, nrow=nb.tip, ncol=30)
        compE <- matrix(0, nrow=nb.tip, ncol=30)
    }
    if(hidden.states == "TEST"){
        states = matrix(0,Ntip(phy),30)
        for(i in 1:Ntip(phy)){
            if(data[i]==1){states[i,1]=1}
            if(data[i]==2){states[i,2]=1}
            if(data[i]==0){states[i,3]=1}
            if(data[i]==4){states[i,4]=1}
            if(data[i]==5){states[i,5]=1}
            if(data[i]==3){states[i,6]=1}
        }
        compD <- matrix(0, nrow=nb.tip, ncol=30)
        compE <- matrix(0, nrow=nb.tip, ncol=30)
    }
    #Initializes the tip sampling and sets internal nodes to be zero:
    ncols = dim(compD)[2]
    if(length(f) == 3){
        for(i in 1:(nb.tip)){
            compD[i,] <- f * states[i,]
            compE[i,] <- rep((1-f), ncols/3)
        }
    }else{
        for(i in 1:(nb.tip)){
            compD[i,] <- f[i] * states[i,]
            compE[i,] <- rep((1-f[i]), ncols/3)
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


SingleChildProbGeo <- function(cache, compD, compE, start.time, end.time, x){
    if(cache$hidden.states == TRUE){
        pars <- c(cache$s00A, cache$s11A, cache$s01A, cache$x00A, cache$x11A, cache$d00A_11A, cache$d00A_01A, cache$d11A_00A, cache$d11A_01A, cache$d01A_00A, cache$d01A_11A, cache$d00A_00B, cache$d00A_00C, cache$d00A_00D, cache$d00A_00E, cache$d00A_00F, cache$d00A_00G, cache$d00A_00H, cache$d00A_00I, cache$d00A_00J, cache$d11A_11B, cache$d11A_11C, cache$d11A_11D, cache$d11A_11E, cache$d11A_11F, cache$d11A_11G, cache$d11A_11H, cache$d11A_11I, cache$d11A_11J, cache$d01A_01B, cache$d01A_01C, cache$d01A_01D, cache$d01A_01E, cache$d01A_01F, cache$d01A_01G, cache$d01A_01H, cache$d01A_01I, cache$d01A_01J, cache$s00B, cache$s11B, cache$s01B, cache$x00B, cache$x11B, cache$d00B_11B, cache$d00B_01B, cache$d11B_00B, cache$d11B_01B, cache$d01B_00B, cache$d01B_11B, cache$d00B_00A, cache$d00B_00C, cache$d00B_00D, cache$d00B_00E, cache$d00B_00F, cache$d00B_00G, cache$d00B_00H, cache$d00B_00I, cache$d00B_00J, cache$d11B_11A, cache$d11B_11C, cache$d11B_11D, cache$d11B_11E, cache$d11B_11F, cache$d11B_11G, cache$d11B_11H, cache$d11B_11I, cache$d11B_11J, cache$d01B_01A, cache$d01B_01C, cache$d01B_01D, cache$d01B_01E, cache$d01B_01F, cache$d01B_01G, cache$d01B_01H, cache$d01B_01I, cache$d01B_01J, cache$s00C, cache$s11C, cache$s01C, cache$x00C, cache$x11C, cache$d00C_11C, cache$d00C_01C, cache$d11C_00C, cache$d11C_01C, cache$d01C_00C, cache$d01C_11C, cache$d00C_00A, cache$d00C_00B, cache$d00C_00D, cache$d00C_00E, cache$d00C_00F, cache$d00C_00G, cache$d00C_00H, cache$d00C_00I, cache$d00C_00J, cache$d11C_11A, cache$d11C_11B, cache$d11C_11D, cache$d11C_11E, cache$d11C_11F, cache$d11C_11G, cache$d11C_11H, cache$d11C_11I, cache$d11C_11J, cache$d01C_01A, cache$d01C_01B, cache$d01C_01D, cache$d01C_01E, cache$d01C_01F, cache$d01C_01G, cache$d01C_01H, cache$d01C_01I, cache$d01C_01J, cache$s00D, cache$s11D, cache$s01D, cache$x00D, cache$x11D, cache$d00D_11D, cache$d00D_01D, cache$d11D_00D, cache$d11D_01D, cache$d01D_00D, cache$d01D_11D, cache$d00D_00A, cache$d00D_00B, cache$d00D_00C, cache$d00D_00E, cache$d00D_00F, cache$d00D_00G, cache$d00D_00H, cache$d00D_00I, cache$d00D_00J, cache$d11D_11A, cache$d11D_11B, cache$d11D_11C, cache$d11D_11E, cache$d11D_11F, cache$d11D_11G, cache$d11D_11H, cache$d11D_11I, cache$d11D_11J, cache$d01D_01A, cache$d01D_01B, cache$d01D_01C, cache$d01D_01E, cache$d01D_01F, cache$d01D_01G, cache$d01D_01H, cache$d01D_01I, cache$d01D_01J, cache$s00E, cache$s11E, cache$s01E, cache$x00E, cache$x11E, cache$d00E_11E, cache$d00E_01E, cache$d11E_00E, cache$d11E_01E, cache$d01E_00E, cache$d01E_11E, cache$d00E_00A, cache$d00E_00B, cache$d00E_00C, cache$d00E_00D, cache$d00E_00F, cache$d00E_00G, cache$d00E_00H, cache$d00E_00I, cache$d00E_00J, cache$d11E_11A, cache$d11E_11B, cache$d11E_11C, cache$d11E_11D, cache$d11E_11F, cache$d11E_11G, cache$d11E_11H, cache$d11E_11I, cache$d11E_11J, cache$d01E_01A, cache$d01E_01B, cache$d01E_01C, cache$d01E_01D, cache$d01E_01F, cache$d01E_01G, cache$d01E_01H, cache$d01E_01I, cache$d01E_01J, cache$s00F, cache$s11F, cache$s01F, cache$x00F, cache$x11F, cache$d00F_11F, cache$d00F_01F, cache$d11F_00F, cache$d11F_01F, cache$d01F_00F, cache$d01F_11F, cache$d00F_00A, cache$d00F_00B, cache$d00F_00C, cache$d00F_00D, cache$d00F_00E, cache$d00F_00G, cache$d00F_00H, cache$d00F_00I, cache$d00F_00J, cache$d11F_11A, cache$d11F_11B, cache$d11F_11C, cache$d11F_11D, cache$d11F_11E, cache$d11F_11G, cache$d11F_11H, cache$d11F_11I, cache$d11F_11J, cache$d01F_01A, cache$d01F_01B, cache$d01F_01C, cache$d01F_01D, cache$d01F_01E, cache$d01F_01G, cache$d01F_01H, cache$d01F_01I, cache$d01F_01J, cache$s00G, cache$s11G, cache$s01G, cache$x00G, cache$x11G, cache$d00G_11G, cache$d00G_01G, cache$d11G_00G, cache$d11G_01G, cache$d01G_00G, cache$d01G_11G, cache$d00G_00A, cache$d00G_00B, cache$d00G_00C, cache$d00G_00D, cache$d00G_00E, cache$d00G_00F, cache$d00G_00H, cache$d00G_00I, cache$d00G_00J, cache$d11G_11A, cache$d11G_11B, cache$d11G_11C, cache$d11G_11D, cache$d11G_11E, cache$d11G_11F, cache$d11G_11H, cache$d11G_11I, cache$d11G_11J, cache$d01G_01A, cache$d01G_01B, cache$d01G_01C, cache$d01G_01D, cache$d01G_01E, cache$d01G_01F, cache$d01G_01H, cache$d01G_01I, cache$d01G_01J, cache$s00H, cache$s11H, cache$s01H, cache$x00H, cache$x11H, cache$d00H_11H, cache$d00H_01H, cache$d11H_00H, cache$d11H_01H, cache$d01H_00H, cache$d01H_11H, cache$d00H_00A, cache$d00H_00B, cache$d00H_00C, cache$d00H_00D, cache$d00H_00E, cache$d00H_00F, cache$d00H_00G, cache$d00H_00I, cache$d00H_00J, cache$d11H_11A, cache$d11H_11B, cache$d11H_11C, cache$d11H_11D, cache$d11H_11E, cache$d11H_11F, cache$d11H_11G, cache$d11H_11I, cache$d11H_11J, cache$d01H_01A, cache$d01H_01B, cache$d01H_01C, cache$d01H_01D, cache$d01H_01E, cache$d01H_01F, cache$d01H_01G, cache$d01H_01I, cache$d01H_01J, cache$s00I, cache$s11I, cache$s01I, cache$x00I, cache$x11I, cache$d00I_11I, cache$d00I_01I, cache$d11I_00I, cache$d11I_01I, cache$d01I_00I, cache$d01I_11I, cache$d00I_00A, cache$d00I_00B, cache$d00I_00C, cache$d00I_00D, cache$d00I_00E, cache$d00I_00F, cache$d00I_00G, cache$d00I_00H, cache$d00I_00J, cache$d11I_11A, cache$d11I_11B, cache$d11I_11C, cache$d11I_11D, cache$d11I_11E, cache$d11I_11F, cache$d11I_11G, cache$d11I_11H, cache$d11I_11J, cache$d01I_01A, cache$d01I_01B, cache$d01I_01C, cache$d01I_01D, cache$d01I_01E, cache$d01I_01F, cache$d01I_01G, cache$d01I_01H, cache$d01I_01J, cache$s00J, cache$s11J, cache$s01J, cache$x00J, cache$x11J, cache$d00J_11J, cache$d00J_01J, cache$d11J_00J, cache$d11J_01J, cache$d01J_00J, cache$d01J_11J, cache$d00J_00A, cache$d00J_00B, cache$d00J_00C, cache$d00J_00D, cache$d00J_00E, cache$d00J_00F, cache$d00J_00G, cache$d00J_00H, cache$d00J_00I, cache$d11J_11A, cache$d11J_11B, cache$d11J_11C, cache$d11J_11D, cache$d11J_11E, cache$d11J_11F, cache$d11J_11G, cache$d11J_11H, cache$d11J_11I, cache$d01J_01A, cache$d01J_01B, cache$d01J_01C, cache$d01J_01D, cache$d01J_01E, cache$d01J_01F, cache$d01J_01G, cache$d01J_01H, cache$d01J_01I)
        yini <- c(E00A = compE[1], E11A = compE[2], E01A = compE[3], E00B = compE[4], E11B = compE[5], E01B = compE[6], E00C = compE[7], E11C = compE[8], E01C = compE[9], E00D = compE[10], E11D = compE[11], E01D = compE[12], E00E = compE[13], E11E = compE[14], E01E = compE[15], E00F = compE[16], E11F = compE[17], E01F = compE[18], E00G = compE[19], E11G = compE[20], E01G = compE[21], E00H = compE[22], E11H = compE[23], E01H = compE[24], E00I = compE[25], E11I = compE[26], E01I = compE[27], E00J = compE[28], E11J = compE[29], E01J = compE[30], D00A = compD[1], D11A = compD[2], D01A = compD[3], D00B = compD[4], D11B = compD[5], D01B = compD[6], D00C = compD[7], D11C = compD[8], D01C = compD[9], D00D = compD[10], D11D = compD[11], D01D = compD[12], D00E = compD[13], D11E = compD[14], D01E = compD[15], D00F = compD[16], D11F = compD[17], D01F = compD[18], D00G = compD[19], D11G = compD[20], D01G = compD[21], D00H = compD[22], D11H = compD[23], D01H = compD[24], D00I = compD[25], D11I = compD[26], D01I = compD[27], D00J = compD[28], D11J = compD[29], D01J = compD[30])
        times=c(start.time, end.time)
        if(cache$assume.cladogenetic == TRUE){
            #prob.subtree.cal.full <-lsoda(yini, times, func = "fgeohisse_derivs", pars, initfunc="initmod_fgeohisse", dll = "fgeohisse-ext-derivs", rtol=1e-8, atol=1e-8)
            prob.subtree.cal.full <- lsoda(yini, times, func = "fgeohisse_derivs", pars, initfunc="initmod_fgeohisse", dllname = "hisse", rtol=1e-8, atol=1e-8)
        }else{
            #prob.subtree.cal.full <- lsoda(yini, times, func = "fnotclasse_more_derivs", pars, initfunc="initmod_fhinoclass", dll = "fnotclasse-more-ext-derivs", rtol=1e-8, atol=1e-8)
            prob.subtree.cal.full <- lsoda(yini, times, func = "fnotclasse_more_derivs", pars, initfunc="initmod_fhinoclass", dllname = "hisse", rtol=1e-8, atol=1e-8)
        }
    }else{
        pars <- c(cache$s00A, cache$s11A, cache$s01A, cache$x00A, cache$x11A, cache$d00A_11A, cache$d00A_01A, cache$d11A_00A, cache$d11A_01A, cache$d01A_00A, cache$d01A_11A)
        yini <-c(E_0=compE[1], E_1=compE[2], E_01=compE[3], D_N0=compD[1], D_N1=compD[2], D_N2=compD[3])
        times=c(start.time, end.time)
        if(cache$assume.cladogenetic == TRUE){
            #prob.subtree.cal.full <- lsoda(yini, times, func = "fclasse_geosse_equivalent_derivs", pars, initfunc="initmod_fgeosse", dll = "fcanonical_geosse-ext-derivs", rtol=1e-8, atol=1e-8)
            prob.subtree.cal.full <- lsoda(yini, times, func = "fclasse_geosse_equivalent_derivs", pars, initfunc="initmod_fgeosse", dllname = "hisse", rtol=1e-8, atol=1e-8)
        }else{
            #prob.subtree.cal.full <-lsoda(yini, times, func = "fnotclasse_derivs", pars, initfunc="initmod_fnoclass", dll = "fnotclasse-ext-derivs", rtol=1e-8, atol=1e-8)
            prob.subtree.cal.full <- lsoda(yini, times, func = "fnotclasse_derivs", pars, initfunc="initmod_fnoclass", dllname = "hisse", rtol=1e-8, atol=1e-8)
        }
    }
    ######## THIS CHECKS TO ENSURE THAT THE INTEGRATION WAS SUCCESSFUL ###########
    if(attributes(prob.subtree.cal.full)$istate[1] < 0){
        prob.subtree.cal <- prob.subtree.cal.full[-1,-1]
        if(cache$hidden.states == TRUE){
            prob.subtree.cal[31:60] <- cache$bad.likelihood
            return(prob.subtree.cal)
        }else{
            prob.subtree.cal[4:6] <- cache$bad.likelihood
            return(prob.subtree.cal)
        }
    }else{
        prob.subtree.cal <- prob.subtree.cal.full[-1,-1]
    }
    ##############################################################################
    
    if(cache$hidden.states == TRUE){
        if(any(is.nan(prob.subtree.cal[31:60]))){
            prob.subtree.cal[31:60] <- cache$bad.likelihood
            return(prob.subtree.cal)
        }
        #This is default and cannot change, but if we get a negative probability, discard the results:
        if(any(prob.subtree.cal[31:60] < 0)){
            prob.subtree.cal[31:60] <- cache$bad.likelihood
            return(prob.subtree.cal)
        }
        if(sum(prob.subtree.cal[31:60]) < cache$ode.eps){
            prob.subtree.cal[31:60] <- cache$bad.likelihood
            return(prob.subtree.cal)
        }
    }else{
        if(is.nan(prob.subtree.cal[4]) | is.nan(prob.subtree.cal[5]) | is.nan(prob.subtree.cal[6])){
            prob.subtree.cal[4:6] <- cache$bad.likelihood
            return(prob.subtree.cal)
        }
        #This is default and cannot change, but if we get a negative probability, discard the results:
        if(prob.subtree.cal[4]<0 | prob.subtree.cal[5]<0 | prob.subtree.cal[6]<0){
            prob.subtree.cal[4:6] <- cache$bad.likelihood
            return(prob.subtree.cal)
        }
        if(sum(prob.subtree.cal[4:6]) < cache$ode.eps){
            prob.subtree.cal[4:6] <- cache$bad.likelihood
            return(prob.subtree.cal)
        }
    }
    return(prob.subtree.cal)
}


FocalNodeProbGeo <- function(cache, dat.tab, generations){
    ### Ughy McUgherson. This is a must in order to pass CRAN checks: http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
    DesNode = NULL
    FocalNode = NULL
    setkey(dat.tab, FocalNode)
    CurrentGenData <- dat.tab[data.table(generations)]
    if(cache$hidden.states == TRUE){
        tmp <- t(apply(CurrentGenData, 1, function(z) SingleChildProbGeo(cache, z[7:36], z[37:66],  z[2], z[1])))
        if(cache$assume.cladogenetic == TRUE){
            #I imagine this is slower than it needs to be. But this is the best I can come up with
            DM <- matrix(tmp[seq(1,nrow(tmp)-1,2),31:60], length(unique(CurrentGenData$FocalNode)), 30)
            DN <- matrix(tmp[seq(2,nrow(tmp),2),31:60], length(unique(CurrentGenData$FocalNode)), 30)
            Ss <- matrix(c(cache$s00A, cache$s11A, cache$s01A, cache$s00B, cache$s11B, cache$s01B, cache$s00C, cache$s11C, cache$s01C, cache$s00D, cache$s11D, cache$s01D, cache$s00E, cache$s11E, cache$s01E, cache$s00F, cache$s11F, cache$s01F, cache$s00G, cache$s11G, cache$s01G, cache$s00H, cache$s11H, cache$s01H, cache$s00I, cache$s11I, cache$s01I, cache$s00J, cache$s11J, cache$s01J), length(unique(CurrentGenData$FocalNode)), 30, byrow=TRUE)
            v.mat <- matrix(0,length(unique(CurrentGenData$FocalNode)), 30)
            for(i in 1:30){
                if(i %% 3 != 0){
                    v.mat[,i] <- DM[,i] * DN[,i] * Ss[,i]
                }else{
                    v.mat[,i] <- 0.5 * (DN[,i] * DM[,i-2] + DN[,i-2] * DM[,i]) * Ss[,i-2] + 0.5 * (DN[,i] * DM[,i-1] + DN[,i-1] * DM[,i]) * Ss[,i-1] + 0.5 * (DN[,i-2] * DM[,i-1] + DN[,i-1] * DM[,i-2]) * Ss[,i]
                }
            }
        }else{
            v.mat <- matrix(tmp[seq(1,nrow(tmp)-1,2),31:60] * tmp[seq(2,nrow(tmp),2),31:60], length(unique(CurrentGenData$FocalNode)), 30)
            v.mat <- v.mat * matrix(c(cache$s00A, cache$s11A, cache$s01A, cache$s00B, cache$s11B, cache$s01B, cache$s00C, cache$s11C, cache$s01C, cache$s00D, cache$s11D, cache$s01D, cache$s00E, cache$s11E, cache$s01E, cache$s00F, cache$s11F, cache$s01F, cache$s00G, cache$s11G, cache$s01G, cache$s00H, cache$s11H, cache$s01H, cache$s00I, cache$s11I, cache$s01I, cache$s00J, cache$s11J, cache$s01J), length(unique(CurrentGenData$FocalNode)), 30, byrow=TRUE)
        }
        phi.mat <- matrix(tmp[seq(1,nrow(tmp)-1,2),1:30], length(unique(CurrentGenData$FocalNode)), 30)
        if(!is.null(cache$node)){
            if(which(generations == cache$node)){
                fixer = numeric(30)
                fixer[cache$state] = 1
                v.mat[which(generations == cache$node),] <- v.mat[which(generations == cache$node),] * fixer
            }
        }
    }else{
        tmp <- t(apply(CurrentGenData, 1, function(z) SingleChildProbGeo(cache, z[7:9], z[10:13],  z[2], z[1])))
        if(cache$assume.cladogenetic == TRUE){
            #I imagine this is slower than it needs to be. But this is the best I can come up with
            DM <- matrix(tmp[seq(1,nrow(tmp)-1,2),4:6], length(unique(CurrentGenData$FocalNode)), 3)
            DN <- matrix(tmp[seq(2,nrow(tmp),2),4:6], length(unique(CurrentGenData$FocalNode)), 3)
            Ss <- matrix(c(cache$s00A, cache$s11A, cache$s01A), length(unique(CurrentGenData$FocalNode)), 3, byrow=TRUE)
            v.mat <- matrix(0,length(unique(CurrentGenData$FocalNode)), 3)
            for(i in 1:3){
                if(i %% 3 != 0){
                    v.mat[,i] <- DM[,i] * DN[,i] * Ss[,i]
                }else{
                    v.mat[,i] <- 0.5 * (DN[,i] * DM[,i-2] + DN[,i-2] * DM[,i]) * Ss[,i-2] + 0.5 * (DN[,i] * DM[,i-1] + DN[,i-1] * DM[,i]) * Ss[,i-1] + 0.5 * (DN[,i-2] * DM[,i-1] + DN[,i-1] * DM[,i-2]) * Ss[,i]
                }
            }
        }else{
            v.mat <- matrix(tmp[seq(1,nrow(tmp)-1,2),4:6] * tmp[seq(2,nrow(tmp),2),4:6], length(unique(CurrentGenData$FocalNode)), 3)
            v.mat <- v.mat * matrix(c(cache$s00A, cache$s11A, cache$s01A), length(unique(CurrentGenData$FocalNode)), 3, byrow=TRUE)
        }
        phi.mat <- matrix(tmp[seq(1,nrow(tmp)-1,2),1:3], length(unique(CurrentGenData$FocalNode)), 3)
        if(!is.null(cache$node)){
            if(which(generations == cache$node)){
                fixer = numeric(3)
                fixer[cache$state] = 1
                v.mat[which(generations == cache$node),] <- v.mat[which(generations == cache$node),] * fixer
            }
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
GetRootProbGeo <- function(cache, dat.tab, generations){
    ### Ughy McUgherson. This is a must in order to pass CRAN checks: http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
    FocalNode = NULL
    
    setkey(dat.tab, FocalNode)
    CurrentGenData <- dat.tab[data.table(generations)]
    if(cache$hidden.states == TRUE){
        tmp <- t(apply(CurrentGenData, 1, function(z) SingleChildProbGeo(cache, z[7:36], z[37:66],  z[2], z[1])))
        if(cache$assume.cladogenetic == TRUE){
            #I imagine this is slower than it needs to be. But this is the best I can come up with
            DM <- matrix(tmp[seq(1,nrow(tmp)-1,2),31:60], length(unique(CurrentGenData$FocalNode)), 30)
            DN <- matrix(tmp[seq(2,nrow(tmp),2),31:60], length(unique(CurrentGenData$FocalNode)), 30)
            Ss <- matrix(c(cache$s00A, cache$s11A, cache$s01A, cache$s00B, cache$s11B, cache$s01B, cache$s00C, cache$s11C, cache$s01C, cache$s00D, cache$s11D, cache$s01D, cache$s00E, cache$s11E, cache$s01E, cache$s00F, cache$s11F, cache$s01F, cache$s00G, cache$s11G, cache$s01G, cache$s00H, cache$s11H, cache$s01H, cache$s00I, cache$s11I, cache$s01I, cache$s00J, cache$s11J, cache$s01J), length(unique(CurrentGenData$FocalNode)), 30, byrow=TRUE)
            v.mat <- matrix(0,length(unique(CurrentGenData$FocalNode)), 30)
            for(i in 1:30){
                if(i %% 3 != 0){
                    v.mat[,i] <- DM[,i] * DN[,i] * Ss[,i]
                }else{
                    v.mat[,i] <- 0.5 * (DN[,i] * DM[,i-2] + DN[,i-2] * DM[,i]) * Ss[,i-2] + 0.5 * (DN[,i] * DM[,i-1] + DN[,i-1] * DM[,i]) * Ss[,i-1] + 0.5 * (DN[,i-2] * DM[,i-1] + DN[,i-1] * DM[,i-2]) * Ss[,i]
                }
            }
        }else{
            v.mat <- matrix(tmp[seq(1,nrow(tmp)-1,2),31:60] * tmp[seq(2,nrow(tmp),2),31:60], length(unique(CurrentGenData$FocalNode)), 30)
            v.mat <- v.mat * matrix(c(cache$s00A, cache$s11A, cache$s01A, cache$s00B, cache$s11B, cache$s01B, cache$s00C, cache$s11C, cache$s01C, cache$s00D, cache$s11D, cache$s01D, cache$s00E, cache$s11E, cache$s01E, cache$s00F, cache$s11F, cache$s01F, cache$s00G, cache$s11G, cache$s01G, cache$s00H, cache$s11H, cache$s01H, cache$s00I, cache$s11I, cache$s01I, cache$s00J, cache$s11J, cache$s01J), length(unique(CurrentGenData$FocalNode)), 30, byrow=TRUE)
        }
        phi.mat <- matrix(tmp[seq(1,nrow(tmp)-1,2),1:30], length(unique(CurrentGenData$FocalNode)), 30)
        if(!is.null(cache$node)){
            if(which(generations == cache$node)){
                fixer = numeric(30)
                fixer[cache$state] = 1
                v.mat[which(generations == cache$node),] <- v.mat[which(generations == cache$node),] * fixer
            }
        }
    }else{
        tmp <- t(apply(CurrentGenData, 1, function(z) SingleChildProbGeo(cache, z[7:9], z[10:13],  z[2], z[1])))
        if(cache$assume.cladogenetic == TRUE){
            #I imagine this is slower than it needs to be. But this is the best I can come up with
            DM <- matrix(tmp[seq(1,nrow(tmp)-1,2),4:6], length(unique(CurrentGenData$FocalNode)), 3)
            DN <- matrix(tmp[seq(2,nrow(tmp),2),4:6], length(unique(CurrentGenData$FocalNode)), 3)
            Ss <- matrix(c(cache$s00A, cache$s11A, cache$s01A), length(unique(CurrentGenData$FocalNode)), 3, byrow=TRUE)
            v.mat <- matrix(0,length(unique(CurrentGenData$FocalNode)), 3)
            for(i in 1:3){
                if(i %% 3 != 0){
                    v.mat[,i] <- DM[,i] * DN[,i] * Ss[,i]
                }else{
                    v.mat[,i] <- 0.5 * (DN[,i] * DM[,i-2] + DN[,i-2] * DM[,i]) * Ss[,i-2] + 0.5 * (DN[,i] * DM[,i-1] + DN[,i-1] * DM[,i]) * Ss[,i-1] + 0.5 * (DN[,i-2] * DM[,i-1] + DN[,i-1] * DM[,i-2]) * Ss[,i]
                }
            }
        }else{
            v.mat <- matrix(tmp[seq(1,nrow(tmp)-1,2),4:6] * tmp[seq(2,nrow(tmp),2),4:6], length(unique(CurrentGenData$FocalNode)), 3)
            v.mat <- v.mat * matrix(c(cache$s00A, cache$s11A, cache$s01A), length(unique(CurrentGenData$FocalNode)), 3, byrow=TRUE)
        }
        phi.mat <- matrix(tmp[seq(1,nrow(tmp)-1,2),1:3], length(unique(CurrentGenData$FocalNode)), 3)
        if(!is.null(cache$node)){
            if(which(generations == cache$node)){
                fixer = numeric(3)
                fixer[cache$state] = 1
                v.mat[which(generations == cache$node),] <- v.mat[which(generations == cache$node),] * fixer
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

DownPassGeoHissefast <- function(dat.tab, gen, cache, condition.on.survival, root.type, root.p, get.phi=FALSE, node=NULL, state=NULL) {
    
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
                    res.tmp <- GetRootProbGeo(cache=cache, dat.tab=dat.tab, generations=gen[[i]])
                    cache$node <- NULL
                    cache$state <- NULL
                }else{
                    res.tmp <- GetRootProbGeo(cache=cache, dat.tab=dat.tab, generations=gen[[i]])
                }
            }else{
                res.tmp <- GetRootProbGeo(cache=cache, dat.tab=dat.tab, generations=gen[[i]])
            }
            if(cache$hidden.states == TRUE){
                compD.root <- res.tmp[c(32:61)]
                compE.root <- res.tmp[c(2:31)]
            }else{
                compD.root <- res.tmp[c(5:7)]
                compE.root <- res.tmp[c(2:4)]
            }
            setkey(dat.tab, DesNode)
            comp <- dat.tab[["comp"]]
            comp <- c(comp[-TIPS], res.tmp[1])
        }else{
            if(!is.null(node)){
                if(node %in% gen[[i]]){
                    cache$node <- node
                    cache$state <- state
                    dat.tab <- FocalNodeProbGeo(cache, dat.tab, gen[[i]])
                    cache$node <- NULL
                    cache$state <- NULL
                }else{
                    dat.tab <- FocalNodeProbGeo(cache, dat.tab, gen[[i]])
                }
            }else{
                dat.tab <- FocalNodeProbGeo(cache, dat.tab, gen[[i]])
            }
        }
    }
    if (is.na(sum(log(compD.root))) || is.na(log(sum(1 - compE.root)))){
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
                    if(cache$assume.cladogenetic == TRUE){
                        lambda <- c(cache$s00A, cache$s11A, sum(c(cache$s00A, cache$s11A, cache$s01A)))
                    }else{
                        lambda <- c(cache$s00A, cache$s11A, cache$s01A)
                    }
                    compD.root <- compD.root / sum(root.p * lambda * (1 - compE.root)^2)
                    #Corrects for possibility that you have 0/0:
                    compD.root[which(is.na(compD.root))] = 0
                    loglik <- log(sum(compD.root * root.p)) + sum(log(comp))
                }else{
                    if(cache$assume.cladogenetic == TRUE){
                        lambda <- c(cache$s00A, cache$s11A, sum(c(cache$s00A, cache$s11A, cache$s01A)))
                    }else{
                        lambda <- c(cache$s00A, cache$s11A, cache$s01A)
                    }
                    compD.root <- compD.root / (lambda * (1 - compE.root)^2)
                    #Corrects for possibility that you have 0/0:
                    compD.root[which(is.na(compD.root))] = 0
                    loglik <- log(sum(compD.root)) + sum(log(comp))
                }
            }else{
                if(root.type == "madfitz"){
                    if(cache$assume.cladogenetic == TRUE) {
                        lambda <- c(cache$s00A, cache$s11A, sum(c(cache$s00A, cache$s11A, cache$s01A)), cache$s00B, cache$s11B, sum(c(cache$s00B, cache$s11B, cache$s01B)), cache$s00C, cache$s11C, sum(c(cache$s00C, cache$s11C, cache$s01C)), cache$s00D, cache$s11D, sum(c(cache$s00D, cache$s11D, cache$s01D)), cache$s00E, cache$s11E, sum(c(cache$s00E, cache$s11E, cache$s01E)), cache$s00F, cache$s11F, sum(c(cache$s00F, cache$s11F, cache$s01F)), cache$s00G, cache$s11G, sum(c(cache$s00G, cache$s11G, cache$s01G)), cache$s00H, cache$s11H, sum(c(cache$s00H, cache$s11H, cache$s01H)), cache$s00I, cache$s11I, sum(c(cache$s00I, cache$s11I, cache$s01I)), cache$s00J, cache$s11J, sum(c(cache$s00J, cache$s11J, cache$s01J)))
                    }else{
                        lambda <- c(cache$s00A, cache$s11A, cache$s01A, cache$s00B, cache$s11B, cache$s01B, cache$s00C, cache$s11C, cache$s01C, cache$s00D, cache$s11D, cache$s01D, cache$s00E, cache$s11E, cache$s01E, cache$s00F, cache$s11F, cache$s01F, cache$s00G, cache$s11G, cache$s01G, cache$s00H, cache$s11H, cache$s01H, cache$s00I, cache$s11I, cache$s01I, cache$s00J, cache$s11J, cache$s01J)
                    }
                    compD.root <- compD.root / sum(root.p * lambda * (1 - compE.root)^2)
                    #Corrects for possibility that you have 0/0:
                    compD.root[which(is.na(compD.root))] = 0
                    loglik <- log(sum(compD.root * root.p)) + sum(log(comp))
                }else{
                    if(cache$assume.cladogenetic == TRUE) {
                        lambda <- c(cache$s00A, cache$s11A, sum(c(cache$s00A, cache$s11A, cache$s01A)), cache$s00B, cache$s11B, sum(c(cache$s00B, cache$s11B, cache$s01B)), cache$s00C, cache$s11C, sum(c(cache$s00C, cache$s11C, cache$s01C)), cache$s00D, cache$s11D, sum(c(cache$s00D, cache$s11D, cache$s01D)), cache$s00E, cache$s11E, sum(c(cache$s00E, cache$s11E, cache$s01E)), cache$s00F, cache$s11F, sum(c(cache$s00F, cache$s11F, cache$s01F)), cache$s00G, cache$s11G, sum(c(cache$s00G, cache$s11G, cache$s01G)), cache$s00H, cache$s11H, sum(c(cache$s00H, cache$s11H, cache$s01H)), cache$s00I, cache$s11I, sum(c(cache$s00I, cache$s11I, cache$s01I)), cache$s00J, cache$s11J, sum(c(cache$s00J, cache$s11J, cache$s01J)))
                    }else{
                        lambda <- c(cache$s00A, cache$s11A, cache$s01A, cache$s00B, cache$s11B, cache$s01B, cache$s00C, cache$s11C, cache$s01C, cache$s00D, cache$s11D, cache$s01D, cache$s00E, cache$s11E, cache$s01E, cache$s00F, cache$s11F, cache$s01F, cache$s00G, cache$s11G, cache$s01G, cache$s00H, cache$s11H, cache$s01H, cache$s00I, cache$s11I, cache$s01I, cache$s00J, cache$s11J, cache$s01J)
                    }
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

ParametersToPassGeoHiSSEfast <- function(model.vec, hidden.states, assume.cladogenetic=TRUE, nb.tip, nb.node, bad.likelihood, ode.eps){
    #Provides an initial object that contains all the parameters to be passed among functions. This will also be used to pass other things are we move down the tree (see DownPassGeoSSE):
    obj <- NULL
    
    obj$hidden.states <- hidden.states
    obj$assume.cladogenetic <- assume.cladogenetic
    obj$nb.tip <- nb.tip
    obj$nb.node <- nb.node
    obj$bad.likelihood <- bad.likelihood
    obj$ode.eps <- ode.eps
    
    
    ######REMAINING ISSUE -- OPTIMIZING FIXED PARAMETERS of TURNOVER AND EPS. Can be done, just a bit trickier than usual.
    ##Hidden State A
    #obj$s00A = model.vec[1]
    #obj$s11A = model.vec[2]
    #obj$s01A = model.vec[3]
    #obj$x00A = model.vec[4]
    #obj$x11A = model.vec[5]
    obj$s00A = model.vec[1] / (1 + model.vec[4])
    obj$s11A = model.vec[2] / (1 + model.vec[5])
    if(model.vec[3] == 0 & assume.cladogenetic==TRUE){
        #For the independent models
        obj$s01A = obj$s00A
    }else{
        obj$s01A = model.vec[3] - obj$s00A - obj$s11A
    }
    obj$x00A = (model.vec[4] * model.vec[1]) / (1 + model.vec[4])
    obj$x11A = (model.vec[5] * model.vec[2]) / (1 + model.vec[5])
    
    obj$d00A_11A = model.vec[6]
    obj$d00A_01A = model.vec[7]
    obj$d11A_00A = model.vec[8]
    obj$d11A_01A = model.vec[9]
    
    #This sets the extinct.frac if necessary
    if(model.vec[10]==0){
        #obj$d01A_00A = model.vec[5]
        obj$d01A_00A = obj$x11A
    }else{
        obj$d01A_00A = model.vec[10]
    }
    if(model.vec[11]==0){
        #obj$d01A_11A = model.vec[4]
        obj$d01A_11A = obj$x00A
    }else{
        obj$d01A_11A = model.vec[11]
    }
    
    obj$d00A_00B = model.vec[12]
    obj$d00A_00C = model.vec[13]
    obj$d00A_00D = model.vec[14]
    obj$d00A_00E = model.vec[15]
    obj$d00A_00F = model.vec[16]
    obj$d00A_00G = model.vec[17]
    obj$d00A_00H = model.vec[18]
    obj$d00A_00I = model.vec[19]
    obj$d00A_00J = model.vec[20]
    obj$d11A_11B = model.vec[21]
    obj$d11A_11C = model.vec[22]
    obj$d11A_11D = model.vec[23]
    obj$d11A_11E = model.vec[24]
    obj$d11A_11F = model.vec[25]
    obj$d11A_11G = model.vec[26]
    obj$d11A_11H = model.vec[27]
    obj$d11A_11I = model.vec[28]
    obj$d11A_11J = model.vec[29]
    obj$d01A_01B = model.vec[30]
    obj$d01A_01C = model.vec[31]
    obj$d01A_01D = model.vec[32]
    obj$d01A_01E = model.vec[33]
    obj$d01A_01F = model.vec[34]
    obj$d01A_01G = model.vec[35]
    obj$d01A_01H = model.vec[36]
    obj$d01A_01I = model.vec[37]
    obj$d01A_01J = model.vec[38]
    
    ##Hidden State B
    #obj$s00B = model.vec[39]
    #obj$s11B = model.vec[40]
    #obj$s01B = model.vec[41]
    #obj$x00B = model.vec[42]
    #obj$x11B = model.vec[43]
    obj$s00B = model.vec[39] / (1 + model.vec[42])
    obj$s11B = model.vec[40] / (1 + model.vec[43])
    if(model.vec[41] == 0 & assume.cladogenetic==TRUE){
        obj$s01B = obj$s00B
    }else{
        obj$s01B = model.vec[41] - obj$s00B - obj$s11B
    }
    obj$x00B = (model.vec[42] * model.vec[39]) / (1 + model.vec[42])
    obj$x11B = (model.vec[43] * model.vec[40]) / (1 + model.vec[43])

    obj$d00B_11B = model.vec[44]
    obj$d00B_01B = model.vec[45]
    obj$d11B_00B = model.vec[46]
    obj$d11B_01B = model.vec[47]
    
    #This sets the extinct.frac if necessary
    if(model.vec[48]==0){
        #obj$d01B_00B = model.vec[43]
        obj$d01B_00B = obj$x11B
    }else{
        obj$d01B_00B = model.vec[48]
    }
    if(model.vec[49]==0){
        #obj$d01B_11B = model.vec[42]
        obj$d01B_11B = obj$x00B
    }else{
        obj$d01B_11B = model.vec[49]
    }
    
    obj$d00B_00A = model.vec[50]
    obj$d00B_00C = model.vec[51]
    obj$d00B_00D = model.vec[52]
    obj$d00B_00E = model.vec[53]
    obj$d00B_00F = model.vec[54]
    obj$d00B_00G = model.vec[55]
    obj$d00B_00H = model.vec[56]
    obj$d00B_00I = model.vec[57]
    obj$d00B_00J = model.vec[58]
    obj$d11B_11A = model.vec[59]
    obj$d11B_11C = model.vec[60]
    obj$d11B_11D = model.vec[61]
    obj$d11B_11E = model.vec[62]
    obj$d11B_11F = model.vec[63]
    obj$d11B_11G = model.vec[64]
    obj$d11B_11H = model.vec[65]
    obj$d11B_11I = model.vec[66]
    obj$d11B_11J = model.vec[67]
    obj$d01B_01A = model.vec[68]
    obj$d01B_01C = model.vec[69]
    obj$d01B_01D = model.vec[70]
    obj$d01B_01E = model.vec[71]
    obj$d01B_01F = model.vec[72]
    obj$d01B_01G = model.vec[73]
    obj$d01B_01H = model.vec[74]
    obj$d01B_01I = model.vec[75]
    obj$d01B_01J = model.vec[76]
    
    ##Hidden State C
    #obj$s00C = model.vec[77]
    #obj$s11C = model.vec[78]
    #obj$s01C = model.vec[79]
    #obj$x00C = model.vec[80]
    #obj$x11C = model.vec[81]
    obj$s00C = model.vec[77] / (1 + model.vec[80])
    obj$s11C = model.vec[78] / (1 + model.vec[81])
    if(model.vec[79] == 0 & assume.cladogenetic==TRUE){
        obj$s01C = obj$s00C
    }else{
        obj$s01C = model.vec[79] - obj$s00C - obj$s11C
    }
    obj$x00C = (model.vec[80] * model.vec[77]) / (1 + model.vec[80])
    obj$x11C = (model.vec[81] * model.vec[78]) / (1 + model.vec[81])

    obj$d00C_11C = model.vec[82]
    obj$d00C_01C = model.vec[83]
    obj$d11C_00C = model.vec[84]
    obj$d11C_01C = model.vec[85]
    
    #This sets the extinct.frac if necessary
    if(model.vec[86]==0){
        #obj$d01C_00C = model.vec[81]
        obj$d01C_00C = obj$x11C
    }else{
        obj$d01C_00C = model.vec[86]
    }
    if(model.vec[87]==0){
        #obj$d01C_11C = model.vec[80]
        obj$d01C_11C = obj$x00C
    }else{
        obj$d01C_11C = model.vec[87]
    }
    
    obj$d00C_00A = model.vec[88]
    obj$d00C_00B = model.vec[89]
    obj$d00C_00D = model.vec[90]
    obj$d00C_00E = model.vec[91]
    obj$d00C_00F = model.vec[92]
    obj$d00C_00G = model.vec[93]
    obj$d00C_00H = model.vec[94]
    obj$d00C_00I = model.vec[95]
    obj$d00C_00J = model.vec[96]
    obj$d11C_11A = model.vec[97]
    obj$d11C_11B = model.vec[98]
    obj$d11C_11D = model.vec[99]
    obj$d11C_11E = model.vec[100]
    obj$d11C_11F = model.vec[101]
    obj$d11C_11G = model.vec[102]
    obj$d11C_11H = model.vec[103]
    obj$d11C_11I = model.vec[104]
    obj$d11C_11J = model.vec[105]
    obj$d01C_01A = model.vec[106]
    obj$d01C_01B = model.vec[107]
    obj$d01C_01D = model.vec[108]
    obj$d01C_01E = model.vec[109]
    obj$d01C_01F = model.vec[110]
    obj$d01C_01G = model.vec[111]
    obj$d01C_01H = model.vec[112]
    obj$d01C_01I = model.vec[113]
    obj$d01C_01J = model.vec[114]
    
    ##Hidden State D
    #obj$s00D = model.vec[115]
    #obj$s11D = model.vec[116]
    #obj$s01D = model.vec[117]
    #obj$x00D = model.vec[118]
    #obj$x11D = model.vec[119]
    obj$s00D = model.vec[115] / (1 + model.vec[118])
    obj$s11D = model.vec[116] / (1 + model.vec[119])
    if(model.vec[117] == 0 & assume.cladogenetic==TRUE){
        obj$s01D = obj$s00D
    }else{
        obj$s01D = model.vec[117] - obj$s00D - obj$s11D
    }
    obj$x00D = (model.vec[118] * model.vec[115]) / (1 + model.vec[118])
    obj$x11D = (model.vec[119] * model.vec[116]) / (1 + model.vec[119])

    obj$d00D_11D = model.vec[120]
    obj$d00D_01D = model.vec[121]
    obj$d11D_00D = model.vec[122]
    obj$d11D_01D = model.vec[123]
    
    #This sets the extinct.frac if necessary
    if(model.vec[124]==0){
        #obj$d01D_00D = model.vec[119]
        obj$d01D_00D = obj$x11D
    }else{
        obj$d01D_00D = model.vec[124]
    }
    if(model.vec[125]==0){
        #obj$d01D_11D = model.vec[118]
        obj$d01D_11D = obj$x00D
    }else{
        obj$d01D_11D = model.vec[125]
    }
    
    obj$d00D_00A = model.vec[126]
    obj$d00D_00B = model.vec[127]
    obj$d00D_00C = model.vec[128]
    obj$d00D_00E = model.vec[129]
    obj$d00D_00F = model.vec[130]
    obj$d00D_00G = model.vec[131]
    obj$d00D_00H = model.vec[132]
    obj$d00D_00I = model.vec[133]
    obj$d00D_00J = model.vec[134]
    obj$d11D_11A = model.vec[135]
    obj$d11D_11B = model.vec[136]
    obj$d11D_11C = model.vec[137]
    obj$d11D_11E = model.vec[138]
    obj$d11D_11F = model.vec[139]
    obj$d11D_11G = model.vec[140]
    obj$d11D_11H = model.vec[141]
    obj$d11D_11I = model.vec[142]
    obj$d11D_11J = model.vec[143]
    obj$d01D_01A = model.vec[144]
    obj$d01D_01B = model.vec[145]
    obj$d01D_01C = model.vec[146]
    obj$d01D_01E = model.vec[147]
    obj$d01D_01F = model.vec[148]
    obj$d01D_01G = model.vec[149]
    obj$d01D_01H = model.vec[150]
    obj$d01D_01I = model.vec[151]
    obj$d01D_01J = model.vec[152]
    
    ##Hidden State E
    #obj$s00E = model.vec[153]
    #obj$s11E = model.vec[154]
    #obj$s01E = model.vec[155]
    #obj$x00E = model.vec[156]
    #obj$x11E = model.vec[157]
    obj$s00E = model.vec[153] / (1 + model.vec[156])
    obj$s11E = model.vec[154] / (1 + model.vec[157])
    if(model.vec[155] == 0 & assume.cladogenetic==TRUE){
        obj$s01E = obj$s00E
    }else{
        obj$s01E = model.vec[155] - obj$s00E - obj$s11E
    }
    obj$x00E = (model.vec[156] * model.vec[153]) / (1 + model.vec[156])
    obj$x11E = (model.vec[157] * model.vec[154]) / (1 + model.vec[157])

    obj$d00E_11E = model.vec[158]
    obj$d00E_01E = model.vec[159]
    obj$d11E_00E = model.vec[160]
    obj$d11E_01E = model.vec[161]
    
    #This sets the extinct.frac if necessary
    if(model.vec[162]==0){
        #obj$d01E_00E = model.vec[157]
        obj$d01E_00E = obj$x11E
    }else{
        obj$d01E_00E = model.vec[162]
    }
    if(model.vec[163]==0){
        #obj$d01E_11E = model.vec[156]
        obj$d01E_11E = obj$x00E
    }else{
        obj$d01E_11E = model.vec[163]
    }
    
    obj$d00E_00A = model.vec[164]
    obj$d00E_00B = model.vec[165]
    obj$d00E_00C = model.vec[166]
    obj$d00E_00D = model.vec[167]
    obj$d00E_00F = model.vec[168]
    obj$d00E_00G = model.vec[169]
    obj$d00E_00H = model.vec[170]
    obj$d00E_00I = model.vec[171]
    obj$d00E_00J = model.vec[172]
    obj$d11E_11A = model.vec[173]
    obj$d11E_11B = model.vec[174]
    obj$d11E_11C = model.vec[175]
    obj$d11E_11D = model.vec[176]
    obj$d11E_11F = model.vec[177]
    obj$d11E_11G = model.vec[178]
    obj$d11E_11H = model.vec[179]
    obj$d11E_11I = model.vec[180]
    obj$d11E_11J = model.vec[181]
    obj$d01E_01A = model.vec[182]
    obj$d01E_01B = model.vec[183]
    obj$d01E_01C = model.vec[184]
    obj$d01E_01D = model.vec[185]
    obj$d01E_01F = model.vec[186]
    obj$d01E_01G = model.vec[187]
    obj$d01E_01H = model.vec[188]
    obj$d01E_01I = model.vec[189]
    obj$d01E_01J = model.vec[190]
    
    ##Hidden State F
    #obj$s00F = model.vec[191]
    #obj$s11F = model.vec[192]
    #obj$s01F = model.vec[193]
    #obj$x00F = model.vec[194]
    #obj$x11F = model.vec[195]
    obj$s00F = model.vec[191] / (1 + model.vec[194])
    obj$s11F = model.vec[192] / (1 + model.vec[195])
    if(model.vec[193] == 0 & assume.cladogenetic==TRUE){
        obj$s01F = obj$s00F
    }else{
        obj$s01F = model.vec[193] - obj$s00F - obj$s11F
    }
    obj$x00F = (model.vec[194] * model.vec[191]) / (1 + model.vec[194])
    obj$x11F = (model.vec[195] * model.vec[192]) / (1 + model.vec[195])

    obj$d00F_11F = model.vec[196]
    obj$d00F_01F = model.vec[197]
    obj$d11F_00F = model.vec[198]
    obj$d11F_01F = model.vec[199]
    
    #This sets the extinct.frac if necessary
    if(model.vec[200]==0){
        #obj$d01F_00F = model.vec[195]
        obj$d01F_00F = obj$x11F
    }else{
        obj$d01F_00F = model.vec[200]
    }
    if(model.vec[201]==0){
        #obj$d01F_11F = model.vec[194]
        obj$d01F_11F = obj$x00F
    }else{
        obj$d01F_11F = model.vec[201]
    }
    
    obj$d00F_00A = model.vec[202]
    obj$d00F_00B = model.vec[203]
    obj$d00F_00C = model.vec[204]
    obj$d00F_00D = model.vec[205]
    obj$d00F_00E = model.vec[206]
    obj$d00F_00G = model.vec[207]
    obj$d00F_00H = model.vec[208]
    obj$d00F_00I = model.vec[209]
    obj$d00F_00J = model.vec[210]
    obj$d11F_11A = model.vec[211]
    obj$d11F_11B = model.vec[212]
    obj$d11F_11C = model.vec[213]
    obj$d11F_11D = model.vec[214]
    obj$d11F_11E = model.vec[215]
    obj$d11F_11G = model.vec[216]
    obj$d11F_11H = model.vec[217]
    obj$d11F_11I = model.vec[218]
    obj$d11F_11J = model.vec[219]
    obj$d01F_01A = model.vec[220]
    obj$d01F_01B = model.vec[221]
    obj$d01F_01C = model.vec[222]
    obj$d01F_01D = model.vec[223]
    obj$d01F_01E = model.vec[224]
    obj$d01F_01G = model.vec[225]
    obj$d01F_01H = model.vec[226]
    obj$d01F_01I = model.vec[227]
    obj$d01F_01J = model.vec[228]
    
    ##Hidden State G
    #obj$s00G = model.vec[229]
    #obj$s11G = model.vec[230]
    #obj$s01G = model.vec[231]
    #obj$x00G = model.vec[232]
    #obj$x11G = model.vec[233]
    obj$s00G = model.vec[229] / (1 + model.vec[232])
    obj$s11G = model.vec[230] / (1 + model.vec[233])
    if(model.vec[231] == 0 & assume.cladogenetic==TRUE){
        obj$s01G = obj$s00G
    }else{
        obj$s01G = model.vec[231] - obj$s00G - obj$s11G
    }
    obj$x00G = (model.vec[232] * model.vec[229]) / (1 + model.vec[232])
    obj$x11G = (model.vec[233] * model.vec[230]) / (1 + model.vec[233])

    obj$d00G_11G = model.vec[234]
    obj$d00G_01G = model.vec[235]
    obj$d11G_00G = model.vec[236]
    obj$d11G_01G = model.vec[237]
    
    #This sets the extinct.frac if necessary
    if(model.vec[238]==0){
        #obj$d01G_00G = model.vec[233]
        obj$d01G_00G = obj$x11G
    }else{
        obj$d01G_00G = model.vec[238]
    }
    if(model.vec[239]==0){
        #obj$d01G_11G = model.vec[232]
        obj$d01G_11G = obj$x00G
    }else{
        obj$d01G_11G = model.vec[239]
    }
    
    obj$d00G_00A = model.vec[240]
    obj$d00G_00B = model.vec[241]
    obj$d00G_00C = model.vec[242]
    obj$d00G_00D = model.vec[243]
    obj$d00G_00E = model.vec[244]
    obj$d00G_00F = model.vec[245]
    obj$d00G_00H = model.vec[246]
    obj$d00G_00I = model.vec[247]
    obj$d00G_00J = model.vec[248]
    obj$d11G_11A = model.vec[249]
    obj$d11G_11B = model.vec[250]
    obj$d11G_11C = model.vec[251]
    obj$d11G_11D = model.vec[252]
    obj$d11G_11E = model.vec[253]
    obj$d11G_11F = model.vec[254]
    obj$d11G_11H = model.vec[255]
    obj$d11G_11I = model.vec[256]
    obj$d11G_11J = model.vec[257]
    obj$d01G_01A = model.vec[258]
    obj$d01G_01B = model.vec[259]
    obj$d01G_01C = model.vec[260]
    obj$d01G_01D = model.vec[261]
    obj$d01G_01E = model.vec[262]
    obj$d01G_01F = model.vec[263]
    obj$d01G_01H = model.vec[264]
    obj$d01G_01I = model.vec[265]
    obj$d01G_01J = model.vec[266]
    
    ##Hidden State H
    #obj$s00H = model.vec[267]
    #obj$s11H = model.vec[268]
    #obj$s01H = model.vec[269]
    #obj$x00H = model.vec[270]
    #obj$x11H = model.vec[271]
    obj$s00H = model.vec[267] / (1 + model.vec[270])
    obj$s11H = model.vec[268] / (1 + model.vec[271])
    if(model.vec[269] == 0 & assume.cladogenetic==TRUE){
        obj$s01H = obj$s00H
    }else{
        obj$s01H = model.vec[269] - obj$s00H - obj$s11H
    }
    obj$x00H = (model.vec[270] * model.vec[267]) / (1 + model.vec[270])
    obj$x11H = (model.vec[271] * model.vec[268]) / (1 + model.vec[271])

    obj$d00H_11H = model.vec[272]
    obj$d00H_01H = model.vec[273]
    obj$d11H_00H = model.vec[274]
    obj$d11H_01H = model.vec[275]
    
    #This sets the extinct.frac if necessary
    if(model.vec[276]==0){
        #obj$d01H_00H = model.vec[271]
        obj$d01H_00H = obj$x11H
    }else{
        obj$d01H_00H = model.vec[276]
    }
    if(model.vec[277]==0){
        #obj$d01H_11H = model.vec[270]
        obj$d01H_11H = obj$x00H
    }else{
        obj$d01H_11H = model.vec[277]
    }
    
    obj$d00H_00A = model.vec[278]
    obj$d00H_00B = model.vec[279]
    obj$d00H_00C = model.vec[280]
    obj$d00H_00D = model.vec[281]
    obj$d00H_00E = model.vec[282]
    obj$d00H_00F = model.vec[283]
    obj$d00H_00G = model.vec[284]
    obj$d00H_00I = model.vec[285]
    obj$d00H_00J = model.vec[286]
    obj$d11H_11A = model.vec[287]
    obj$d11H_11B = model.vec[288]
    obj$d11H_11C = model.vec[289]
    obj$d11H_11D = model.vec[290]
    obj$d11H_11E = model.vec[291]
    obj$d11H_11F = model.vec[292]
    obj$d11H_11G = model.vec[293]
    obj$d11H_11I = model.vec[294]
    obj$d11H_11J = model.vec[295]
    obj$d01H_01A = model.vec[296]
    obj$d01H_01B = model.vec[297]
    obj$d01H_01C = model.vec[298]
    obj$d01H_01D = model.vec[299]
    obj$d01H_01E = model.vec[300]
    obj$d01H_01F = model.vec[301]
    obj$d01H_01G = model.vec[302]
    obj$d01H_01I = model.vec[303]
    obj$d01H_01J = model.vec[304]
    
    ##Hidden State I
    #obj$s00I = model.vec[305]
    #obj$s11I = model.vec[306]
    #obj$s01I = model.vec[307]
    #obj$x00I = model.vec[308]
    #obj$x11I = model.vec[309]
    obj$s00I = model.vec[305] / (1 + model.vec[308])
    obj$s11I = model.vec[306] / (1 + model.vec[309])
    if(model.vec[307] == 0 & assume.cladogenetic==TRUE){
        obj$s01I = obj$s00I
    }else{
        obj$s01I = model.vec[307] - obj$s00I - obj$s11I
    }
    obj$x00I = (model.vec[308] * model.vec[305]) / (1 + model.vec[308])
    obj$x11I = (model.vec[309] * model.vec[306]) / (1 + model.vec[309])

    obj$d00I_11I = model.vec[310]
    obj$d00I_01I = model.vec[311]
    obj$d11I_00I = model.vec[312]
    obj$d11I_01I = model.vec[313]
    
    #This sets the extinct.frac if necessary
    if(model.vec[314]==0){
        #obj$d01I_00I = model.vec[309]
        obj$d01I_00I = obj$x11I
    }else{
        obj$d01I_00I = model.vec[314]
    }
    if(model.vec[315]==0){
        #obj$d01I_11I = model.vec[308]
        obj$d01I_11I = obj$x00I
    }else{
        obj$d01I_11I = model.vec[315]
    }
    
    obj$d00I_00A = model.vec[316]
    obj$d00I_00B = model.vec[317]
    obj$d00I_00C = model.vec[318]
    obj$d00I_00D = model.vec[319]
    obj$d00I_00E = model.vec[320]
    obj$d00I_00F = model.vec[321]
    obj$d00I_00G = model.vec[322]
    obj$d00I_00H = model.vec[323]
    obj$d00I_00J = model.vec[324]
    obj$d11I_11A = model.vec[325]
    obj$d11I_11B = model.vec[326]
    obj$d11I_11C = model.vec[327]
    obj$d11I_11D = model.vec[328]
    obj$d11I_11E = model.vec[329]
    obj$d11I_11F = model.vec[330]
    obj$d11I_11G = model.vec[331]
    obj$d11I_11H = model.vec[332]
    obj$d11I_11J = model.vec[333]
    obj$d01I_01A = model.vec[334]
    obj$d01I_01B = model.vec[335]
    obj$d01I_01C = model.vec[336]
    obj$d01I_01D = model.vec[337]
    obj$d01I_01E = model.vec[338]
    obj$d01I_01F = model.vec[339]
    obj$d01I_01G = model.vec[340]
    obj$d01I_01H = model.vec[341]
    obj$d01I_01J = model.vec[342]
    
    ##Hidden State J
    #obj$s00J = model.vec[343]
    #obj$s11J = model.vec[344]
    #obj$s01J = model.vec[345]
    #obj$x00J = model.vec[346]
    #obj$x11J = model.vec[347]
    obj$s00J = model.vec[343] / (1 + model.vec[346])
    obj$s11J = model.vec[344] / (1 + model.vec[347])
    if(model.vec[345] == 0 & assume.cladogenetic==TRUE){
        obj$s01J = obj$s00J
    }else{
        obj$s01J = model.vec[345] - obj$s00J - obj$s11J
    }
    obj$x00J = (model.vec[346] * model.vec[343]) / (1 + model.vec[346])
    obj$x11J = (model.vec[347] * model.vec[344]) / (1 + model.vec[347])

    obj$d00J_11J = model.vec[348]
    obj$d00J_01J = model.vec[349]
    obj$d11J_00J = model.vec[350]
    obj$d11J_01J = model.vec[351]
    
    #This sets the extinct.frac if necessary
    if(model.vec[352]==0){
        #obj$d01J_00J = model.vec[347]
        obj$d01J_00J = obj$x11J
    }else{
        obj$d01J_00J = model.vec[352]
    }
    if(model.vec[353]==0){
        #obj$d01J_11J = model.vec[346]
        obj$d01J_11J = obj$x00J
    }else{
        obj$d01J_11J = model.vec[353]
    }
    
    obj$d00J_00A = model.vec[354]
    obj$d00J_00B = model.vec[355]
    obj$d00J_00C = model.vec[356]
    obj$d00J_00D = model.vec[357]
    obj$d00J_00E = model.vec[358]
    obj$d00J_00F = model.vec[359]
    obj$d00J_00G = model.vec[360]
    obj$d00J_00H = model.vec[361]
    obj$d00J_00I = model.vec[362]
    obj$d11J_11A = model.vec[363]
    obj$d11J_11B = model.vec[364]
    obj$d11J_11C = model.vec[365]
    obj$d11J_11D = model.vec[366]
    obj$d11J_11E = model.vec[367]
    obj$d11J_11F = model.vec[368]
    obj$d11J_11G = model.vec[369]
    obj$d11J_11H = model.vec[370]
    obj$d11J_11I = model.vec[371]
    obj$d01J_01A = model.vec[372]
    obj$d01J_01B = model.vec[373]
    obj$d01J_01C = model.vec[374]
    obj$d01J_01D = model.vec[375]
    obj$d01J_01E = model.vec[376]
    obj$d01J_01F = model.vec[377]
    obj$d01J_01G = model.vec[378]
    obj$d01J_01H = model.vec[379]
    obj$d01J_01I = model.vec[380]
    
    return(obj)
}


