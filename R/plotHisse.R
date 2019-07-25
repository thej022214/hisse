
######################################################################################################################################
######################################################################################################################################
### Plotting function for HiSSE
######################################################################################################################################
######################################################################################################################################

## Setting the functions to be full S3 and congruent with RCran S3 requirements.
plot.hisse.states <- function(x, rate.param = "net.div", type = "fan", show.tip.label = TRUE, fsize = 1.0, legend = "tips", ...) {
    hisse.results <- x
    if( inherits(hisse.results, what=c("hisse.states","list")) ){
        if( inherits(hisse.results, what="hisse.states") ){
            ## we have to make a list so we can run this generally
            if(is.null(hisse.results$aic)){
                ## If a user forgot to include the aic, then we add a random value in for them
                hisse.results$aic = 42
            }
            tmp.list <- list()
            tmp.list[[1]] <- hisse.results
            hisse.results <- tmp.list
        } else { ## Then it is a list.
            ## If x is a list we need to check if all elements have the $aic to make the model average.
            any.other.class <- any( sapply(hisse.results, function(x) !inherits(x, what="hisse.states") ) )
            if( any.other.class ) stop("All elements of the list 'x' need to be of class 'hisse.states'.")

            any.missing <- any( sapply(hisse.results, function(x) is.null(x$aic) ) )
            if( any.missing ) stop( "If x is a list, then each reconstruction need to have $aic information in order to make the model average." )

        }
    } else {
        stop( "x needs to be an object of class 'hisse.states'." )
    }

    ## Check the parameters for the functions:
    rate.param <- match.arg(rate.param, c("turnover", "net.div", "speciation", "extinction"
                                        , "extinction.fraction") )
    type <- match.arg(type, c("fan", "phylogram") )
    if( !is.logical( show.tip.label ) ) stop("show.tip.label needs TRUE or FALSE")
    legend <- match.arg(legend, c("none", "traditional", "tips", "internal", "all") )

    ## ############################################
    ## Block to get the parameters from the ... list.
    ## Will also check the parameters to see if they are valid.
    if(hasArg(do.observed.only)){
        do.observed.only <- list(...)$do.observed.only
        if( !is.logical( do.observed.only ) ) stop("do.observed.only needs TRUE or FALSE")
    } else{
        do.observed.only <- TRUE
    }
    if(hasArg(rate.colors)){
        rate.colors <- list(...)$rate.colors
    } else{
        rate.colors <- NULL
    }
    if(hasArg(state.colors)){
        state.colors <- list(...)$state.colors
    } else{
        state.colors <- NULL
    }
    if(hasArg(edge.width)){
        edge.width <- list(...)$edge.width
        if( !is.numeric( edge.width ) ) stop("edge.width needs a single numeric value")
    } else{
        edge.width <- 5
    }
    if(hasArg(width.factor)){
        width.factor <- list(...)$width.factor
        if( !is.numeric( width.factor ) ) stop("width.factor needs a single numeric value")
    } else{
        width.factor <- 0.5
    }
    if(hasArg(mar)){
        mar <- list(...)$mar
        if( !is.numeric( mar ) ) stop("mar needs a numeric vector")
    } else{
        mar <- c(0.1,0.1,0.1,0.1)
    }
    if(hasArg(outline)){
        outline <- list(...)$outline
        if( !is.logical( outline ) ) stop("outline needs TRUE or FALSE")
    } else{
        outline <- FALSE
    }
    if(hasArg(outline.color)){
        outline.color <- list(...)$outline.color
    } else{
        outline.color <- "gray"
    }
    if(hasArg(rate.range)){
        rate.range <- list(...)$rate.range
    } else{
        rate.range <- NULL
    }
    if(hasArg(swap.underscore)){
        swap.underscore <- list(...)$swap.underscore
        if( !is.logical( swap.underscore ) ) stop("swap.underscore needs TRUE or FALSE")
    } else{
        swap.underscore <- TRUE
    }
    if(hasArg(lims.percentage.correction)){
        lims.percentage.correction <- list(...)$lims.percentage.correction
        if( !is.numeric( lims.percentage.correction ) ) stop("lims.percentage.correction needs a numeric value")
    } else{
        lims.percentage.correction <- 0.001
    }
    if(hasArg(legend.position)){
        legend.position <- list(...)$legend.position
        if( !is.numeric( legend.position ) ) stop("legend.position needs a numeric vector")
    } else{
        legend.position <- c(0, 0.2, 0, 0.2)
    }
    if(hasArg(legend.cex)){
        legend.cex <- list(...)$legend.cex
        if( !is.numeric( legend.cex ) ) stop("legend.cex needs a numeric value")
    } else{
        legend.cex <- 0.4
    }
    if(hasArg(legend.kernel.rates)){
        legend.kernel.rates <- list(...)$legend.kernel.rates
    } else{
        legend.kernel.rates <- "auto"
    }
    if(hasArg(legend.kernel.states)){
        legend.kernel.states <- list(...)$legend.kernel.states
    } else{
        legend.kernel.states <- "auto"
    }
    if(hasArg(legend.bg)){
        legend.bg <- list(...)$legend.bg
    } else{
        legend.bg <- "cornsilk3"
    }
    ## ############################################


    ## Going to change par here. So need to save current par and return to previous at the end. Good practice!
    old.par <- par(no.readonly=T)
    if(is.null(rate.colors)) {
        rate.colors <- c("blue", "red")
    }
    if(is.null(state.colors)) {
        state.colors <- c("white", "black")
    }
    rates.tips <- ConvertManyToRate(hisse.results, rate.param, "tip.mat")
    rates.internal <- ConvertManyToRate(hisse.results, rate.param, "node.mat")
    states.tips <- NA
    states.internal <- NA
    if (do.observed.only) {
        states.tips <- ConvertManyToBinaryState(hisse.results, "tip.mat")
        states.internal <- ConvertManyToBinaryState(hisse.results, "node.mat")
    } else {
        stop("So far we can easily plot just the binary observed state; if you want to plot the hidden states, use a different function")
    }
    tree.to.plot <- hisse.results[[1]]$phy
    #	rate.tree <- contMapGivenAnc(tree=hisse.object$phy, x=ConvertToRate(hisse.object$tip.mat, rate.vector= rate.vector), plot=FALSE, anc.states=ConvertToRate(hisse.object$node.mat, rate.vector= rate.vector), ...)
    rate.lims <- range(c(rates.tips, rates.internal))
    if(!is.null(rate.range)) {
        if(min(rate.range) > min(rate.lims) | max(rate.range) < max(rate.lims)) {
            warning(paste("Did not override rate.lims: the specified rate.range (", rate.range[1], ", ", rate.range[2], ") did not contain all the values for the observed rates (", rate.lims[1], ", ", rate.lims[2], ")"))
        } else {
            rate.lims <- rate.range
        }
    }
    rate.lims[1] <- rate.lims[1] - lims.percentage.correction*abs(rate.lims[1])
    rate.lims[2] <- rate.lims[2] + lims.percentage.correction*abs(rate.lims[2])

    rate.tree <- contMapGivenAnc(tree= tree.to.plot, x=rates.tips, plot=FALSE, anc.states=rates.internal, lims=rate.lims)
    #change colors
    rate.colors <- colorRampPalette(rate.colors, space="Lab")(length(rate.tree$cols))
    rate.tree$cols[] <- rate.colors
    #rate.tree$cols[] <- adjustcolor(rate.tree$cols[], alpha.f=0.3)
    #state.tree <- contMapGivenAnc(tree=hisse.object$phy, x=ConvertToBinaryState(hisse.object$tip.mat, state.0.indices=state.0.indices), plot=FALSE, anc.states=ConvertToBinaryState(hisse.object$node.mat, state.0.indices=state.0.indices))

    state.lims <- range(c(states.tips, states.internal))
    state.lims[1] <- state.lims[1] - lims.percentage.correction*abs(state.lims[1])
    state.lims[2] <- state.lims[2] + lims.percentage.correction*abs(state.lims[2])

    state.tree <- contMapGivenAnc(tree=tree.to.plot, x=states.tips, plot=FALSE, anc.states=states.internal, lims=state.lims)
    #state.colors <- grey(seq(1,0,length.out=length(state.tree$cols)))
    state.colors <- colorRampPalette(state.colors, space="Lab")(length(rate.tree$cols))
    state.tree$cols[]<- state.colors

    par(fig=c(0,1, 0, 1), new=FALSE)
    plot.contMapHisse(A=rate.tree, B=state.tree, lwd.factor=width.factor, fsize=fsize,
    , add=FALSE, lwd=edge.width, type=type, mar=mar, direction="rightwards"
    , offset=NULL, xlim=NULL, ylim=NULL, hold=TRUE, swap.underscore=swap.underscore
    , outline=outline, outline.color=outline.color, show.tiplabels=show.tip.label)

    if(legend!="none") {
        par(fig=legend.position, new=TRUE)
        plot(x=c(-0.1,1.1), y=c(-1.5,1.5), xlab="", ylab="", bty="n", type="n", xaxt="n", yaxt="n")
        rect(-0.1,-1.1,1.1,1.1, border=NA, col=legend.bg)
        rates.to.plot <- c()
        states.to.plot <- c()
        if(legend=="all" | legend=="tips") {
            rates.to.plot <- append(rates.to.plot, rates.tips)
            states.to.plot <- append(states.to.plot, states.tips)
        }
        if(legend=="all" | legend=="internal") {
            rates.to.plot <- append(rates.to.plot, rates.internal)
            states.to.plot <- append(states.to.plot, states.internal)
        }

        if(legend.kernel.rates=="auto") {
            if(length(unique(rates.to.plot))<=4) {
                legend.kernel.rates <- "hist"
            } else {
                legend.kernel.rates <- "rectangular"
            }
        }
        if(legend.kernel.states=="auto") {
            if(length(unique(states.to.plot))<=4) {
                legend.kernel.states <- "hist"
            } else {
                legend.kernel.states <- "rectangular"
            }
        }

        rates.density <- GetNormalizedDensityPlot(rates.to.plot, rate.lims, legend.kernel.rates)
        states.density <- GetNormalizedDensityPlot(states.to.plot, state.lims, legend.kernel.states)
        states.density$y <- (-1) * states.density$y #so it gets drawn below the other one
        # rates.density <- c()
        # states.density <- c()
        # if (legend.kernel=="hist") {
        # rates.density<-hist(rates.to.plot, breaks=seq(from=rate.lims[1], to=rate.lims[2], length.out = max(100,nclass.Sturges(rates.to.plot)+2)), plot=FALSE)
        # states.density<-hist(states.to.plot, breaks=seq(from=state.lims[1], to=state.lims[2], length.out = max(100,nclass.Sturges(states.to.plot)+2)), plot=FALSE)
        # rates.density$x <- rates.density$mid
        # rates.density$y <- rates.density$density
        # states.density$x <- states.density$mid
        # states.density$y <- states.density$density
        # } else {
        # rates.density <- density(rates.to.plot, from=rate.lims[1], to=rate.lims[2], kernel=legend.kernel)
        # states.density <- density(states.to.plot, from=state.lims[1], to=state.lims[2], kernel=legend.kernel)
        # }
        # rates.density$x <- (rates.density$x - rate.lims[1]) / (rate.lims[2]-rate.lims[1]) #so it goes from zero to one
        # rates.density$y <- rates.density$y/max(rates.density$y)
        # states.density$y <- (-1) * states.density$y/max(states.density$y)
        # states.density$x <- (states.density$x - state.lims[1]) / (state.lims[2]-state.lims[1]) #so it goes from zero to one

        # if(legend=="rect") {
        # rates.density$y <- rep(1, length(rates.density$y))
        # states.density$y <- rep(-1, length(states.density$y))
        # }
        par(lend=1)
        segments(x0=rates.density$x, y0=rep(0, length(rates.density$y)), y1=rates.density$y, col=rate.colors[1+as.integer(round((length(rate.colors)-1)* rates.density$x))], lwd=ifelse(legend.kernel.rates=="hist",4,1))
        text(x=0, y=1.2, labels=format(rate.lims[1], digits=2), cex=legend.cex)
        text(x=1, y=1.2, labels=format(rate.lims[2], digits=2), cex=legend.cex)
        text(x=0.5, y=1.2, labels=rate.param, cex=legend.cex)
        #		lines(rates.density$x, rates.density$y, lwd=0.5, col="gray")

        segments(x0=states.density$x, y0=rep(0, length(states.density$y)), y1=states.density$y, col=state.colors[1+as.integer(round((length(state.colors)-1)* states.density$x))], lwd=ifelse(legend.kernel.states=="hist",4,1))
        text(x=0, y=-1.2, labels="0", cex=legend.cex)
        text(x=1, y=-1.2, labels="1", cex=legend.cex)
        text(x=0.5, y=-1.2, labels="State", cex=legend.cex)
        #		lines(states.density$x, states.density$y, lwd=0.5, col="gray")
        #		lines(rates.density$x, 0*rates.density$y, lwd=0.5, col="gray")

    }

    ## Return par to the previous state.
    par( old.par )
    return(list(rate.tree=rate.tree, state.tree=state.tree))
}


plot.geohisse.states <- function(x, rate.param = "net.div", type = "fan", show.tip.label = TRUE, fsize = 1.0, legend = TRUE, ...){
    hisse.results <- x
    if( inherits(hisse.results, what=c("hisse.geosse.states","list")) ){
        if( inherits(hisse.results, what="hisse.geosse.states") ){
            ## we have to make a list so we can run this generally
            if(is.null(hisse.results$aic)){
                ## If a user forgot to include the aic, then we add a random value in for them
                hisse.results$aic = 42
            }
            tmp.list <- list()
            tmp.list[[1]] <- hisse.results
            hisse.results <- tmp.list
        } else { ## Then it is a list.
            ## If x is a list we need to check if all elements have the $aic to make the model average.
            any.other.class <- any( sapply(hisse.results, function(x) !inherits(x, what="hisse.geosse.states") ) )
            if( any.other.class ) stop("All elements of the list 'x' need to be of class 'hisse.geosse.states'.")
            any.missing <- any( sapply(hisse.results, function(x) is.null(x$aic) ) )
            if( any.missing ) stop( "If x is a list, then each reconstruction need to have $aic information in order to make the model average." )

        }
    } else {
        stop( "x needs to be an object of class 'hisse.geosse.states'." )
    }

    ## Check the values for the parameters:
    rate.param <- match.arg(rate.param, c("turnover", "net.div", "speciation", "extinction"
                                        , "extinction.fraction") )
    type <- match.arg(type, c("fan", "phylogram") )
    if( !is.logical( show.tip.label ) ) stop("show.tip.label needs TRUE or FALSE")
    if( !is.numeric( fsize ) ) stop("fsize needs a numeric value")
    if( !is.logical( legend ) ) stop("legend needs TRUE or FALSE")

    ## ############################################
    ## Block to get the parameters from the ... list.
    ## Will also check the parameters to see if they are valid.
    if(hasArg(do.observed.only)){
        do.observed.only <- list(...)$do.observed.only
        if( !is.logical( do.observed.only ) ) stop("do.observed.only needs TRUE or FALSE")
    } else{
        do.observed.only <- TRUE
    }
    if(hasArg(rate.colors)){
        rate.colors <- list(...)$rate.colors
    } else{
        rate.colors <- NULL
    }
    if(hasArg(state.colors)){
        state.colors <- list(...)$state.colors
    } else{
        state.colors <- NULL
    }
    if(hasArg(edge.width)){
        edge.width <- list(...)$edge.width
        if( !is.numeric( edge.width ) ) stop("edge.width needs a single numeric value")
    } else{
        edge.width <- 5
    }
    if(hasArg(width.factor)){
        width.factor <- list(...)$width.factor
        if( !is.numeric( width.factor ) ) stop("width.factor needs a single numeric value")
    } else{
        width.factor <- 0.5
    }
    if(hasArg(mar)){
        mar <- list(...)$mar
        if( !is.numeric( mar ) ) stop("mar needs a numeric vector")
    } else{
        mar <- c(0.1,0.1,0.1,0.1)
    }
    if(hasArg(outline)){
        outline <- list(...)$outline
        if( !is.logical( outline ) ) stop("outline needs TRUE or FALSE")
    } else{
        outline <- FALSE
    }
    if(hasArg(outline.color)){
        outline.color <- list(...)$outline.color
    } else{
        outline.color <- "gray"
    }
    if(hasArg(rate.range)){
        rate.range <- list(...)$rate.range
    } else{
        rate.range <- NULL
    }
    if(hasArg(swap.underscore)){
        swap.underscore <- list(...)$swap.underscore
        if( !is.logical( swap.underscore ) ) stop("swap.underscore needs TRUE or FALSE")
    } else{
        swap.underscore <- TRUE
    }
    if(hasArg(lims.percentage.correction)){
        lims.percentage.correction <- list(...)$lims.percentage.correction
        if( !is.numeric( lims.percentage.correction ) ) stop("lims.percentage.correction needs a numeric value")
    } else{
        lims.percentage.correction <- 0.001
    }
    if(hasArg(legend.position)){
        legend.position <- list(...)$legend.position
        if( !is.numeric( legend.position ) ) stop("legend.position needs a numeric vector")
    } else{
        legend.position <- c(0, 0.2, 0, 0.2)
    }
    if(hasArg(legend.cex)){
        legend.cex <- list(...)$legend.cex
        if( !is.numeric( legend.cex ) ) stop("legend.cex needs a numeric value")
    } else{
        legend.cex <- 0.4
    }
    if(hasArg(legend.kernel)){
        legend.kernel <- list(...)$legend.kernel
    } else{
        legend.kernel <- "auto"
    }
    if(hasArg(legend.bg)){
        legend.bg <- list(...)$legend.bg
    } else{
        legend.bg <- "cornsilk3"
    }
    ## ############################################

    ## Going to change par here. So need to save current par and return to previous at the end. Good practice!
    old.par <- par(no.readonly=T)

    ## This sets the coordinates of the figure and facilitate localization.
    par(fig=c(0, 1, 0, 1), new=FALSE)
    if(is.null(rate.colors)) {
        rate.colors <- c("blue", "red")
    }
    if(is.null(state.colors)) {
        state.colors <- c("white", "black", "yellow")
        print("Using default colors: white (state 1), black (state 2), and yellow (state 0).")
    }
    rates.tips <- ConvertManyToRate(hisse.results, rate.param, "tip.mat")
    rates.internal <- ConvertManyToRate(hisse.results, rate.param, "node.mat")
    states.tips <- NA
    states.internal <- NA
    if (do.observed.only) {
        ## states.tips is a table with 3 columns.
        states.tips <- ConvertManyToMultiState(hisse.results, "tip.mat")
        states.internal <- ConvertManyToMultiState(hisse.results, "node.mat")
    } else {
        stop("So far we can easily plot just the binary observed state; if you want to plot the hidden states, use a different function")
    }
    tree.to.plot <- hisse.results[[1]]$phy
    rate.lims <- range(c(rates.tips, rates.internal))
    if(!is.null(rate.range)) {
        if(min(rate.range) > min(rate.lims) | max(rate.range) < max(rate.lims)) {
            warning(paste("Did not override rate.lims: the specified rate.range (", rate.range[1], ", ", rate.range[2], ") did not contain all the values for the observed rates (", rate.lims[1], ", ", rate.lims[2], ")"))
        } else {
            rate.lims <- rate.range
        }
    }
    rate.lims[1] <- rate.lims[1] - lims.percentage.correction*abs(rate.lims[1])
    rate.lims[2] <- rate.lims[2] + lims.percentage.correction*abs(rate.lims[2])

    rate.tree <- contMapGivenAnc(tree= tree.to.plot, x=rates.tips, plot=FALSE, anc.states=rates.internal, lims=rate.lims)
    ## change colors
    rate.colors <- colorRampPalette(rate.colors, space="Lab")(length(rate.tree$cols))
    rate.tree$cols[] <- rate.colors

    states.tips.tmp <- rowSums(states.tips %*% c(0,1,2))
    names(states.tips.tmp) <- 1:length(states.tips.tmp)
    states.internal.tmp <- rowSums(states.internal %*% c(0,1,2))
    names(states.internal.tmp) <- (length(states.tips.tmp)+1):(length(states.internal.tmp)+length(states.tips.tmp))
    state.lims <- range(c(states.tips.tmp, states.internal.tmp))
    state.lims[1] <- state.lims[1] - lims.percentage.correction*abs(state.lims[1])
    state.lims[2] <- state.lims[2] + lims.percentage.correction*abs(state.lims[2])
    state.tree <- contMapGivenAnc(tree=tree.to.plot, x=states.tips.tmp, plot=FALSE, anc.states=states.internal.tmp, lims=state.lims)
    state.colors <- colorRampPalette(state.colors, space="Lab")(length(rate.tree$cols))
    state.tree$cols[]<- state.colors

    ## Make the plot.
    ## This is NOT using the phytools version.
    ## This is a special plotting function for HiSSE.
    plot.contMapHisse(A=rate.tree, B=state.tree, lwd.factor=width.factor, fsize=fsize,
    , add=FALSE, lwd=edge.width, type=type, mar=mar, direction="rightwards"
    , offset=NULL, xlim=NULL, ylim=NULL, hold=TRUE, swap.underscore=swap.underscore
    , outline=outline, outline.color=outline.color, show.tiplabels=show.tip.label)

    ## Need to make the legend now:
    ## Supporting only the plot of the rates in this version:
    if( legend ){
        rates.to.plot <- rates.tips
        if(legend.kernel=="auto") {
            if(length(unique(rates.to.plot))<=4) {
                legend.kernel <- "hist"
            } else {
                legend.kernel <- "rectangular"
            }
        }

        ## Make the plot:
        par(fig=legend.position, new=TRUE)
        plot(x=c(-0.1, 1.1), y=c(0, 1.5), xlab="", ylab="", bty="n", type="n", xaxt="n", yaxt="n")
        rect(-0.1, 0, 1.1, 1.1, border=NA, col=legend.bg)
        par(lend=1)
        rates.density <- GetNormalizedDensityPlot(rates.to.plot, rate.lims, legend.kernel)
        segments(x0=rates.density$x, y0=rep(0, length(rates.density$y)), y1=rates.density$y
        , col=rate.colors[1+as.integer(round((length(rate.colors)-1)* rates.density$x))]
        , lwd=ifelse(legend.kernel=="hist",4,1))
        text(x=0, y=1.2, labels=format(rate.lims[1], digits=2), cex=legend.cex)
        text(x=1, y=1.2, labels=format(rate.lims[2], digits=2), cex=legend.cex)
        text(x=0.5, y=1.2, labels=rate.param, cex=legend.cex)
    }

    ## Return par to the previous state.
    par( old.par )
    ## The tiplables for the individual trees will be strange. But this seems fine.
    return(list(rate.tree=rate.tree, state.tree=state.tree))
}

plot.muhisse.states <- function(x, rate.param = "net.div", type = "fan", show.tip.label = TRUE, fsize = 1.0, legend = TRUE, ...){
    hisse.results <- x
    if( inherits(hisse.results, what=c("muhisse.states","list")) ){
        if( inherits(hisse.results, what="muhisse.states") ){
            ## we have to make a list so we can run this generally
            if(is.null(hisse.results$aic)){
                ## If a user forgot to include the aic, then we add a random value in for them
                hisse.results$aic = 42
            }
            tmp.list <- list()
            tmp.list[[1]] <- hisse.results
            hisse.results <- tmp.list
        } else { ## Then it is a list.
            ## If x is a list we need to check if all elements have the $aic to make the model average.
            any.other.class <- any( sapply(hisse.results, function(x) !inherits(x, what="muhisse.states") ) )
            if( any.other.class ) stop("All elements of the list 'x' need to be of class 'muhisse.states'.")
            any.missing <- any( sapply(hisse.results, function(x) is.null(x$aic) ) )
            if( any.missing ) stop( "If x is a list, then each reconstruction need to have $aic information in order to make the model average." )

        }
    } else {
        stop( "x needs to be an object of class 'muhisse.states'." )
    }

    ## Check the values for the parameters:
    rate.param <- match.arg(rate.param, c("turnover", "net.div", "speciation", "extinction"
                                        , "extinction.fraction") )
    type <- match.arg(type, c("fan", "phylogram") )
    if( !is.logical( show.tip.label ) ) stop("show.tip.label needs TRUE or FALSE")
    if( !is.numeric( fsize ) ) stop("fsize needs a numeric value")
    if( !is.logical( legend ) ) stop("legend needs TRUE or FALSE")

    ## ############################################
    ## Block to get the parameters from the ... list.
    ## Will also check the parameters to see if they are valid.
    if(hasArg(do.observed.only)){
        do.observed.only <- list(...)$do.observed.only
        if( !is.logical( do.observed.only ) ) stop("do.observed.only needs TRUE or FALSE")
    } else{
        do.observed.only <- TRUE
    }
    if(hasArg(rate.colors)){
        rate.colors <- list(...)$rate.colors
    } else{
        rate.colors <- NULL
    }
    if(hasArg(state.colors)){
        state.colors <- list(...)$state.colors
    } else{
        state.colors <- NULL
    }
    if(hasArg(edge.width)){
        edge.width <- list(...)$edge.width
        if( !is.numeric( edge.width ) ) stop("edge.width needs a single numeric value")
    } else{
        edge.width <- 5
    }
    if(hasArg(width.factor)){
        width.factor <- list(...)$width.factor
        if( !is.numeric( width.factor ) ) stop("width.factor needs a single numeric value")
    } else{
        width.factor <- 0.5
    }
    if(hasArg(mar)){
        mar <- list(...)$mar
        if( !is.numeric( mar ) ) stop("mar needs a numeric vector")
    } else{
        mar <- c(0.1,0.1,0.1,0.1)
    }
    if(hasArg(outline)){
        outline <- list(...)$outline
        if( !is.logical( outline ) ) stop("outline needs TRUE or FALSE")
    } else{
        outline <- FALSE
    }
    if(hasArg(outline.color)){
        outline.color <- list(...)$outline.color
    } else{
        outline.color <- "gray"
    }
    if(hasArg(rate.range)){
        rate.range <- list(...)$rate.range
    } else{
        rate.range <- NULL
    }
    if(hasArg(swap.underscore)){
        swap.underscore <- list(...)$swap.underscore
        if( !is.logical( swap.underscore ) ) stop("swap.underscore needs TRUE or FALSE")
    } else{
        swap.underscore <- TRUE
    }
    if(hasArg(lims.percentage.correction)){
        lims.percentage.correction <- list(...)$lims.percentage.correction
        if( !is.numeric( lims.percentage.correction ) ) stop("lims.percentage.correction needs a numeric value")
    } else{
        lims.percentage.correction <- 0.001
    }
    if(hasArg(legend.position)){
        legend.position <- list(...)$legend.position
        if( !is.numeric( legend.position ) ) stop("legend.position needs a numeric vector")
    } else{
        legend.position <- c(0, 0.2, 0, 0.2)
    }
    if(hasArg(legend.cex)){
        legend.cex <- list(...)$legend.cex
        if( !is.numeric( legend.cex ) ) stop("legend.cex needs a numeric value")
    } else{
        legend.cex <- 0.4
    }
    if(hasArg(legend.kernel)){
        legend.kernel <- list(...)$legend.kernel
    } else{
        legend.kernel <- "auto"
    }
    if(hasArg(legend.bg)){
        legend.bg <- list(...)$legend.bg
    } else{
        legend.bg <- "cornsilk3"
    }
    ## ############################################

    ## Going to change par here. So need to save current par and return to previous at the end. Good practice!
    old.par <- par(no.readonly=T)

    ## This sets the coordinates of the figure and facilitate localization.
    par(fig=c(0, 1, 0, 1), new=FALSE)
    if(is.null(rate.colors)) {
        rate.colors <- c("blue", "red")
    }
    if(is.null(state.colors)) {
        state.colors <- c("white", "black", "yellow", "green")
        print("Using default colors: white (state 00), black (state 01), yellow (state 10), and green (state 11).")
    }
    rates.tips <- ConvertManyToRate(hisse.results, rate.param, "tip.mat")
    rates.internal <- ConvertManyToRate(hisse.results, rate.param, "node.mat")
    states.tips <- NA
    states.internal <- NA
    if (do.observed.only) {
        ## states.tips is a table with 3 columns.
        states.tips <- ConvertManyToMultiState(hisse.results, "tip.mat")
        states.internal <- ConvertManyToMultiState(hisse.results, "node.mat")
    } else {
        stop("So far we can easily plot just the binary observed state; if you want to plot the hidden states, use a different function")
    }
    tree.to.plot <- hisse.results[[1]]$phy
    rate.lims <- range(c(rates.tips, rates.internal))
    if(!is.null(rate.range)) {
        if(min(rate.range) > min(rate.lims) | max(rate.range) < max(rate.lims)) {
            warning(paste("Did not override rate.lims: the specified rate.range (", rate.range[1], ", ", rate.range[2], ") did not contain all the values for the observed rates (", rate.lims[1], ", ", rate.lims[2], ")"))
        } else {
            rate.lims <- rate.range
        }
    }
    rate.lims[1] <- rate.lims[1] - lims.percentage.correction*abs(rate.lims[1])
    rate.lims[2] <- rate.lims[2] + lims.percentage.correction*abs(rate.lims[2])

    rate.tree <- contMapGivenAnc(tree= tree.to.plot, x=rates.tips, plot=FALSE, anc.states=rates.internal, lims=rate.lims)
    ## change colors
    rate.colors <- colorRampPalette(rate.colors, space="Lab")(length(rate.tree$cols))
    rate.tree$cols[] <- rate.colors

    states.tips.tmp <- rowSums(states.tips %*% c(1,2,3,4))
    names(states.tips.tmp) <- 1:length(states.tips.tmp)
    states.internal.tmp <- rowSums(states.internal %*% c(1,2,3,4))
    names(states.internal.tmp) <- (length(states.tips.tmp)+1):(length(states.internal.tmp)+length(states.tips.tmp))
    state.lims <- range(c(states.tips.tmp, states.internal.tmp))
    state.lims[1] <- state.lims[1] - lims.percentage.correction*abs(state.lims[1])
    state.lims[2] <- state.lims[2] + lims.percentage.correction*abs(state.lims[2])
    state.tree <- contMapGivenAnc(tree=tree.to.plot, x=states.tips.tmp, plot=FALSE, anc.states=states.internal.tmp, lims=state.lims)
    state.colors <- colorRampPalette(state.colors, space="Lab")(length(rate.tree$cols))
    state.tree$cols[]<- state.colors

    ## Make the plot.
    ## This is NOT using the phytools version.
    ## This is a special plotting function for HiSSE.
    plot.contMapHisse(A=rate.tree, B=state.tree, lwd.factor=width.factor, fsize=fsize,
    , add=FALSE, lwd=edge.width, type=type, mar=mar, direction="rightwards"
    , offset=NULL, xlim=NULL, ylim=NULL, hold=TRUE, swap.underscore=swap.underscore
    , outline=outline, outline.color=outline.color, show.tiplabels=show.tip.label)

    ## Need to make the legend now:
    ## Supporting only the plot of the rates in this version:
    if( legend ){
        rates.to.plot <- rates.tips
        if(legend.kernel=="auto") {
            if(length(unique(rates.to.plot))<=4) {
                legend.kernel <- "hist"
            } else {
                legend.kernel <- "rectangular"
            }
        }

        ## Make the plot:
        par(fig=legend.position, new=TRUE)
        plot(x=c(-0.1, 1.1), y=c(0, 1.5), xlab="", ylab="", bty="n", type="n", xaxt="n", yaxt="n")
        rect(-0.1, 0, 1.1, 1.1, border=NA, col=legend.bg)
        par(lend=1)
        rates.density <- GetNormalizedDensityPlot(rates.to.plot, rate.lims, legend.kernel)
        segments(x0=rates.density$x, y0=rep(0, length(rates.density$y)), y1=rates.density$y
        , col=rate.colors[1+as.integer(round((length(rate.colors)-1)* rates.density$x))]
        , lwd=ifelse(legend.kernel=="hist",4,1))
        text(x=0, y=1.2, labels=format(rate.lims[1], digits=2), cex=legend.cex)
        text(x=1, y=1.2, labels=format(rate.lims[2], digits=2), cex=legend.cex)
        text(x=0.5, y=1.2, labels=rate.param, cex=legend.cex)
    }

    ## Return par to the previous state.
    par( old.par )
    ## The tiplables for the individual trees will be strange. But this seems fine.
    return(list(rate.tree=rate.tree, state.tree=state.tree))
}



plot.misse.states <- function(x, rate.param = "net.div", type = "fan", show.tip.label = TRUE, fsize = 1.0, legend = "tips", ...) {
    hisse.results <- x
    if( inherits(hisse.results, what=c("misse.states","list")) ){
        if( inherits(hisse.results, what="misse.states") ){
            ## we have to make a list so we can run this generally
            if(is.null(hisse.results$aic)){
                ## If a user forgot to include the aic, then we add a random value in for them
                hisse.results$aic = 42
            }
            tmp.list <- list()
            tmp.list[[1]] <- hisse.results
            hisse.results <- tmp.list
        } else { ## Then it is a list.
            ## If x is a list we need to check if all elements have the $aic to make the model average.
            any.other.class <- any( sapply(hisse.results, function(x) !inherits(x, what="misse.states") ) )
            if( any.other.class ) stop("All elements of the list 'x' need to be of class 'misse.states'.")

            any.missing <- any( sapply(hisse.results, function(x) is.null(x$aic) ) )
            if( any.missing ) stop( "If x is a list, then each reconstruction need to have $aic information in order to make the model average." )

        }
    } else {
        stop( "x needs to be an object of class 'misse.states'." )
    }

    ## Check the parameters for the functions:
    rate.param <- match.arg(rate.param, c("turnover", "net.div", "speciation", "extinction"
    , "extinction.fraction") )
    type <- match.arg(type, c("fan", "phylogram") )
    if( !is.logical( show.tip.label ) ) stop("show.tip.label needs TRUE or FALSE")
    legend <- match.arg(legend, c("none", "traditional", "tips", "internal", "all") )

    ## ############################################
    ## Block to get the parameters from the ... list.
    ## Will also check the parameters to see if they are valid.
    if(hasArg(do.observed.only)){
        do.observed.only <- list(...)$do.observed.only
        if( !is.logical( do.observed.only ) ) stop("do.observed.only needs TRUE or FALSE")
    } else{
        do.observed.only <- TRUE
    }
    if(hasArg(rate.colors)){
        rate.colors <- list(...)$rate.colors
    } else{
        rate.colors <- NULL
    }
    if(hasArg(state.colors)){
        state.colors <- list(...)$state.colors
    } else{
        state.colors <- NULL
    }
    if(hasArg(edge.width)){
        edge.width <- list(...)$edge.width
        if( !is.numeric( edge.width ) ) stop("edge.width needs a single numeric value")
    } else{
        edge.width <- 1
    }
    if(hasArg(width.factor)){
        width.factor <- list(...)$width.factor
        if( !is.numeric( width.factor ) ) stop("width.factor needs a single numeric value")
    } else{
        width.factor <- 0
    }
    if(hasArg(mar)){
        mar <- list(...)$mar
        if( !is.numeric( mar ) ) stop("mar needs a numeric vector")
    } else{
        mar <- c(0.1,0.1,0.1,0.1)
    }
    if(hasArg(outline)){
        outline <- list(...)$outline
        if( !is.logical( outline ) ) stop("outline needs TRUE or FALSE")
    } else{
        outline <- FALSE
    }
    if(hasArg(outline.color)){
        outline.color <- list(...)$outline.color
    } else{
        outline.color <- "gray"
    }
    if(hasArg(rate.range)){
        rate.range <- list(...)$rate.range
    } else{
        rate.range <- NULL
    }
    if(hasArg(swap.underscore)){
        swap.underscore <- list(...)$swap.underscore
        if( !is.logical( swap.underscore ) ) stop("swap.underscore needs TRUE or FALSE")
    } else{
        swap.underscore <- TRUE
    }
    if(hasArg(lims.percentage.correction)){
        lims.percentage.correction <- list(...)$lims.percentage.correction
        if( !is.numeric( lims.percentage.correction ) ) stop("lims.percentage.correction needs a numeric value")
    } else{
        lims.percentage.correction <- 0.001
    }
    if(hasArg(legend.position)){
        legend.position <- list(...)$legend.position
        if( !is.numeric( legend.position ) ) stop("legend.position needs a numeric vector")
    } else{
        legend.position <- c(0, 0.2, 0, 0.2)
    }
    if(hasArg(legend.cex)){
        legend.cex <- list(...)$legend.cex
        if( !is.numeric( legend.cex ) ) stop("legend.cex needs a numeric value")
    } else{
        legend.cex <- 0.4
    }
    if(hasArg(legend.kernel.rates)){
        legend.kernel.rates <- list(...)$legend.kernel.rates
    } else{
        legend.kernel.rates <- "auto"
    }
    if(hasArg(legend.kernel.states)){
        legend.kernel.states <- list(...)$legend.kernel.states
    } else{
        legend.kernel.states <- "auto"
    }
    if(hasArg(legend.bg)){
        legend.bg <- list(...)$legend.bg
    } else{
        legend.bg <- "cornsilk3"
    }
    ## ############################################


    ## Going to change par here. So need to save current par and return to previous at the end. Good practice!
    old.par <- par(no.readonly=T)
    if(is.null(rate.colors)) {
        rate.colors <- c("blue", "red")
    }
    if(is.null(state.colors)) {
        state.colors <- c("white", "black")
    }
    rates.tips <- ConvertManyToRate(hisse.results, rate.param, "tip.mat")
    rates.internal <- ConvertManyToRate(hisse.results, rate.param, "node.mat")
    states.tips <- NA
    states.internal <- NA
    if (do.observed.only) {
        states.tips <- rep(0, dim(hisse.results[[1]][["tip.mat"]])[1])
        states.internal <- rep(0, dim(hisse.results[[1]][["node.mat"]])[1])
    }
    tree.to.plot <- hisse.results[[1]]$phy
    #	rate.tree <- contMapGivenAnc(tree=hisse.object$phy, x=ConvertToRate(hisse.object$tip.mat, rate.vector= rate.vector), plot=FALSE, anc.states=ConvertToRate(hisse.object$node.mat, rate.vector= rate.vector), ...)
    rate.lims <- range(c(rates.tips, rates.internal))
    if(!is.null(rate.range)) {
        if(min(rate.range) > min(rate.lims) | max(rate.range) < max(rate.lims)) {
            warning(paste("Did not override rate.lims: the specified rate.range (", rate.range[1], ", ", rate.range[2], ") did not contain all the values for the observed rates (", rate.lims[1], ", ", rate.lims[2], ")"))
        } else {
            rate.lims <- rate.range
        }
    }
    rate.lims[1] <- rate.lims[1] - lims.percentage.correction*abs(rate.lims[1])
    rate.lims[2] <- rate.lims[2] + lims.percentage.correction*abs(rate.lims[2])

    rate.tree <- contMapGivenAnc(tree= tree.to.plot, x=rates.tips, plot=FALSE, anc.states=rates.internal, lims=rate.lims)
    #change colors
    rate.colors <- colorRampPalette(rate.colors, space="Lab")(length(rate.tree$cols))
    rate.tree$cols[] <- rate.colors
    #rate.tree$cols[] <- adjustcolor(rate.tree$cols[], alpha.f=0.3)
    #state.tree <- contMapGivenAnc(tree=hisse.object$phy, x=ConvertToBinaryState(hisse.object$tip.mat, state.0.indices=state.0.indices), plot=FALSE, anc.states=ConvertToBinaryState(hisse.object$node.mat, state.0.indices=state.0.indices))

    #state.lims <- range(c(states.tips, states.internal))
    #state.lims[1] <- state.lims[1] - lims.percentage.correction*abs(state.lims[1])
    #state.lims[2] <- state.lims[2] + lims.percentage.correction*abs(state.lims[2])

    #state.tree <- contMapGivenAnc(tree=tree.to.plot, x=states.tips, plot=FALSE, anc.states=states.internal, lims=state.lims)
    #state.colors <- grey(seq(1,0,length.out=length(state.tree$cols)))
    #state.colors <- colorRampPalette(state.colors, space="Lab")(length(rate.tree$cols))
    #state.tree$cols[]<- state.colors

    par(fig=c(0,1, 0, 1), new=FALSE)
    plot.contMapHisse(A=rate.tree, B=rate.tree, lwd.factor=width.factor, fsize=fsize,
    , add=FALSE, lwd=edge.width, type=type, mar=mar, direction="rightwards"
    , offset=NULL, xlim=NULL, ylim=NULL, hold=TRUE, swap.underscore=swap.underscore
    , outline=outline, outline.color=outline.color, show.tiplabels=show.tip.label)

    if(legend!="none") {
        par(fig=legend.position, new=TRUE)
        plot(x=c(-0.1,1.1), y=c(-1.5,1.5), xlab="", ylab="", bty="n", type="n", xaxt="n", yaxt="n")
        rect(-0.1,-1.1,1.1,1.1, border=NA, col=legend.bg)
        rates.to.plot <- c()
        states.to.plot <- c()
        if(legend=="all" | legend=="tips") {
            rates.to.plot <- append(rates.to.plot, rates.tips)
            states.to.plot <- append(states.to.plot, states.tips)
        }
        if(legend=="all" | legend=="internal") {
            rates.to.plot <- append(rates.to.plot, rates.internal)
            states.to.plot <- append(states.to.plot, states.internal)
        }

        if(legend.kernel.rates=="auto") {
            if(length(unique(rates.to.plot))<=4) {
                legend.kernel.rates <- "hist"
            } else {
                legend.kernel.rates <- "rectangular"
            }
        }
        if(legend.kernel.states=="auto") {
            if(length(unique(states.to.plot))<=4) {
                legend.kernel.states <- "hist"
            } else {
                legend.kernel.states <- "rectangular"
            }
        }

        rates.density <- GetNormalizedDensityPlot(rates.to.plot, rate.lims, legend.kernel.rates)
        #states.density <- GetNormalizedDensityPlot(states.to.plot, state.lims, legend.kernel.states)
        #states.density$y <- (-1) * states.density$y #so it gets drawn below the other one
        # rates.density <- c()
        # states.density <- c()
        # if (legend.kernel=="hist") {
        # rates.density<-hist(rates.to.plot, breaks=seq(from=rate.lims[1], to=rate.lims[2], length.out = max(100,nclass.Sturges(rates.to.plot)+2)), plot=FALSE)
        # states.density<-hist(states.to.plot, breaks=seq(from=state.lims[1], to=state.lims[2], length.out = max(100,nclass.Sturges(states.to.plot)+2)), plot=FALSE)
        # rates.density$x <- rates.density$mid
        # rates.density$y <- rates.density$density
        # states.density$x <- states.density$mid
        # states.density$y <- states.density$density
        # } else {
        # rates.density <- density(rates.to.plot, from=rate.lims[1], to=rate.lims[2], kernel=legend.kernel)
        # states.density <- density(states.to.plot, from=state.lims[1], to=state.lims[2], kernel=legend.kernel)
        # }
        # rates.density$x <- (rates.density$x - rate.lims[1]) / (rate.lims[2]-rate.lims[1]) #so it goes from zero to one
        # rates.density$y <- rates.density$y/max(rates.density$y)
        # states.density$y <- (-1) * states.density$y/max(states.density$y)
        # states.density$x <- (states.density$x - state.lims[1]) / (state.lims[2]-state.lims[1]) #so it goes from zero to one

        # if(legend=="rect") {
        # rates.density$y <- rep(1, length(rates.density$y))
        # states.density$y <- rep(-1, length(states.density$y))
        # }
        par(lend=1)
        segments(x0=rates.density$x, y0=rep(0, length(rates.density$y)), y1=rates.density$y, col=rate.colors[1+as.integer(round((length(rate.colors)-1)* rates.density$x))], lwd=ifelse(legend.kernel.rates=="hist",4,1))
        text(x=0, y=1.2, labels=format(rate.lims[1], digits=2), cex=legend.cex)
        text(x=1, y=1.2, labels=format(rate.lims[2], digits=2), cex=legend.cex)
        text(x=0.5, y=1.2, labels=rate.param, cex=legend.cex)
        #		lines(rates.density$x, rates.density$y, lwd=0.5, col="gray")

        #segments(x0=states.density$x, y0=rep(0, length(states.density$y)), y1=states.density$y, col=state.colors[1+as.integer(round((length(state.colors)-1)* states.density$x))], lwd=ifelse(legend.kernel.states=="hist",4,1))
        #text(x=0, y=-1.2, labels="0", cex=legend.cex)
        #text(x=1, y=-1.2, labels="1", cex=legend.cex)
        #text(x=0.5, y=-1.2, labels="State", cex=legend.cex)
        #		lines(states.density$x, states.density$y, lwd=0.5, col="gray")
        #		lines(rates.density$x, 0*rates.density$y, lwd=0.5, col="gray")

    }

    ## Return par to the previous state.
    par( old.par )
    return(list(rate.tree=rate.tree))
}


GetNormalizedDensityPlot <- function(x, limits, kernel, min.breaks=100) {
    x.density <- c()

    if (kernel=="hist") {
        x.density<-hist(x, breaks=seq(from= limits[1], to= limits[2], length.out = max(min.breaks,nclass.Sturges(x)+2)), plot=FALSE)
        x.density$x <- x.density$mid
        x.density$y <- x.density$density
        x.density$x <- x.density$x[which(x.density$y>0)] #since the line is thick, do not plot it if zero
        x.density$y <- x.density$y[which(x.density$y>0)]
    } else {
        if(kernel=="traditional") {
            x.density <- density(x, from=limits[1], to=limits[2])
            x.density$y <- rep(1, length(x.density$y))
        } else {
            x.density <- density(x, from=limits[1], to=limits[2], kernel=kernel)
        }
    }
    x.density$x <- (x.density$x - limits[1]) / (limits[2]-limits[1]) #so it goes from zero to one
    if(kernel!="traditional") {
        x.density$y <- x.density$y/max(x.density$y)
    }
    return(x.density)
}


# function plots reconstructed values for ancestral characters along the edges of the tree
# Modified by Brian O'Meara, June 9, 2015
# Modified by Daniel Caetano, April 4, 2018
contMapGivenAnc <-function(tree,x,res=100,fsize=NULL,ftype=NULL,lwd=4,legend=NULL,
lims=NULL,outline=TRUE,sig=3,type="phylogram",direction="rightwards",
plot=TRUE,anc.states=NULL,...){
    if(hasArg(mar)) mar<-list(...)$mar
    else mar<-rep(0.3,4)
    if(hasArg(offset)) offset<-list(...)$offset
    else offset<-NULL
    if(hasArg(method)) method<-list(...)$method
    else method<-"fastAnc"
    if(hasArg(hold)) hold<-list(...)$hold
    else hold<-TRUE
    h<-max(nodeHeights(tree))
    steps<-0:res/res*max(h)
    H<-nodeHeights(tree)
    a <- anc.states #BCO modified
    if(is.null(a)) { #BCO put this in if loop
        if(method=="fastAnc") a<-fastAnc(tree,x)
        else {
            fit<-anc.ML(tree,x)
            a<-fit$ace
            if(!is.null(fit$missing.x)) x<-c(x,fit$missing.x)
        }
    } #end BCO if loop
    names(x) <- tree$tip.label[as.numeric(names(x))]
    y <- c(a, x[tree$tip.label])
    ## Fixed a problem here. Previous version was calling 'tree$tip' which is not an element of tree.
    names(y)[1:length(tree$tip.label)+tree$Nnode] <- 1:length(tree$tip.label)
    A<-matrix(y[as.character(tree$edge)],nrow(tree$edge),ncol(tree$edge))
    cols<-rainbow(1001,start=0,end=0.7); names(cols)<-0:1000
    if(is.null(lims)) lims<-c(min(c(a,x)),max(c(a,x))) #modified by BCO to include anc state in range for lims
    trans<-0:1000/1000*(lims[2]-lims[1])+lims[1]; names(trans)<-0:1000
    tree$maps <- list(rep(rep(NA, 2), nrow(tree$edge)))
    for(i in 1:nrow(tree$edge)){
        XX<-cbind(c(H[i,1],steps[intersect(which(steps>H[i,1]),which(steps<H[i,2]))]),
        c(steps[intersect(which(steps>H[i,1]),which(steps<H[i,2]))],H[i,2]))-H[i,1]
        YY<-rowMeans(XX)
        if(!all(YY==0)){
            b<-vector()
            for(j in 1:length(YY))
            b[j]<-(A[i,1]/YY[j]+A[i,2]/(max(XX)-YY[j]))/(1/YY[j]+1/(max(XX)-YY[j]))
        } else b<-A[i,1]
        d<-sapply(b,getState,trans=trans)
        tree$maps[[i]]<-XX[,2]-XX[,1]
        names(tree$maps[[i]])<-d
    }
    tree$mapped.edge<-makeMappedEdge(tree$edge,tree$maps)
    tree$mapped.edge<-tree$mapped.edge[,order(as.numeric(colnames(tree$mapped.edge)))]
    xx<-list(tree=tree,cols=cols,lims=lims)
    class(xx) <- "contMap"
    if(plot){
        plot.contMapHisse(xx, fsize=fsize, ftype=ftype, lwd=lwd, outline=outline
        , type=type, mar=mar, direction=direction, offset=offset
        , hold=hold, swap.underscore=TRUE, show.tiplabels=TRUE)
    }
    invisible(xx)
}


# The following function is from phytools, which is released under GPL2+
# It is not exported from phytools, so the options are phytools:::getState or copy it here
# Given that phytools could change non-exported functions in a way that breaks the above code
# I have elected to copy it unchanged.
# Brian O'Meara, June 10, 2015
getState<-function(x,trans){
    i<-1
    state <- names(trans)[1] #BCO: added to prevent error when while loop not entered
    while(x>trans[i] & i <= length(trans)){
        state<-names(trans)[i]
        i<-i+1
    }
    return(state)
}

# make a mapped edge matrix
# The following function is from phytools, which is released under GPL2+
# It is not exported from phytools, so the options are phytools:::makeMappedEdge or copy it here
# Given that phytools could change non-exported functions in a way that breaks the above code
# I have elected to copy it unchanged. Phytools remains a
# Brian O'Meara, June 10, 2015
makeMappedEdge<-function(edge,maps){
    st<-sort(unique(unlist(sapply(maps,function(x) names(x)))))
    mapped.edge<-matrix(0,nrow(edge),length(st))
    rownames(mapped.edge)<-apply(edge,1,function(x) paste(x,collapse=","))
    colnames(mapped.edge)<-st
    for(i in 1:length(maps))
    for(j in 1:length(maps[[i]]))
    mapped.edge[i,names(maps[[i]])[j]]<-mapped.edge[i,names(maps[[i]])[j]]+maps[[i]][j]
    return(mapped.edge)
}

MergeAcrossHidden <- function(x) {
    x.trimmed <- x[,-1]
    hidden.deleted <- gsub("[^0-9]", "", colnames(x.trimmed))
    all.observed.states <- unique(hidden.deleted)
    x.new <- data.frame(matrix(0, nrow=length(x[,1]), ncol=length(all.observed.states)))
    for (state.index in sequence(length(all.observed.states))) {
        state.indices <- which(hidden.deleted==all.observed.states[state.index])
        x.new[,state.index] <- rowSums(data.frame(x.trimmed[,state.indices])) #wrapping in data.frame since it drops to vector if length(state.indices)==1
    }
    x.new <- x.new/rowSums(x.new) #normalize
    rownames(x.new) <- x[,1]
    colnames(x.new) <- all.observed.states
    return(x.new)
}


#Get prob of it being 1
ConvertToBinaryState <- function(x) {
    x.trimmed <- x[,-1]
    state.0.indices <- which(grepl("0", colnames(x.trimmed)))
    state.1.indices <- which(grepl("1", colnames(x.trimmed)))
    x0.trimmed <- x.trimmed[, state.0.indices]
    if (is.null(dim(x0.trimmed)[1])) { #is a vector
        x0.trimmed <- matrix(data=x0.trimmed, nrow=length(x0.trimmed), ncol=1)
    }
    x1.trimmed <- x.trimmed[, state.1.indices]
    if (is.null(dim(x1.trimmed)[1])) { #is a vector
        x1.trimmed <- matrix(data=x1.trimmed, nrow=length(x1.trimmed), ncol=1)
    }
    result.vector.0 <- apply(x0.trimmed, 1, sum)
    result.vector.1 <- apply(x1.trimmed, 1, sum)
    result.vector <- result.vector.1 / (result.vector.0 + result.vector.1)
    names(result.vector) <- x[,1]
    return(result.vector)
}


ConvertToRate <- function(x, rate.vector) {
    x.trimmed <- x[,-1]
    ## Here is the problem!! This will cause the division by infinite!
    x.trimmed <- x.trimmed / rowSums(x.trimmed) #normalize to 1 (which it should be already)
    result.vector <- rep(0, dim(x.trimmed)[1])
    for (i in sequence(dim(x.trimmed)[2])) {
        result.vector <- result.vector + rate.vector[i] * x.trimmed[,i]	#get weighted mean
    }
    names(result.vector) <- x[,1]
    return(result.vector)
}

ConvertManyToRate <- function(hisse.results, rate.param, which.element, AIC.weights=NULL) {
    if( is.null(AIC.weights) ){
        AIC.weights <- GetAICWeights(hisse.results)
    }
    storage.matrix <- matrix(nrow=dim(hisse.results[[1]][[which.element]])[1], ncol=0)
    for (i in sequence(length(hisse.results))) {
        rate.vector <- hisse.results[[i]]$rates.mat[rate.param,]
        storage.matrix <- cbind(storage.matrix, ConvertToRate(x=hisse.results[[i]][[which.element]], rate.vector=rate.vector))
    }
    final.results <- apply(storage.matrix, 1, weighted.mean, w=AIC.weights)
    return(final.results)
}

CheckReconBounds <- function(x, n.models, AIC.weights, bound.par.matrix){
    ## Check if every column of the matrices in the list x is within the bounds set to the each of the parameters.
    ## Drop all the models that are not. Models are the columns in each of the matrices in the list x.

    check.keep.mod <- matrix(nrow = 5, ncol = n.models)
    for( i in 1:5 ){
        lim.range <- bound.par.matrix[i, ]
        check.keep.mod[i,] <- apply(x[[i]], 2, function(y) all( y >= lim.range[1] & y <= lim.range[2] ) )
    }
    ## Check if all parameters for the particular model to be averaged are within the bounds.
    keep.mod <- apply(check.keep.mod, 2, all)
    if( !all(keep.mod) ){
        warning( paste0(" Models in position ", paste(which(!keep.mod), collapse=", ")," have parameters outside the bounds defined by 'bound.matrix' argument. These will NOT be included in the reconstruction.") )
    }
    if( sum( keep.mod ) > 1 ){
        final.results <- lapply(x, function(y) apply(y[,keep.mod], 1, weighted.mean, w=AIC.weights[keep.mod]) )
    }
    if( sum( keep.mod ) == 1 ){
        final.results <- lapply(x, function(y) y[,keep.mod]) ## No need to make the weighted.mean
    }
    if( sum( keep.mod ) < 1 ) stop( "No models left to reconstruct! Check if parameter estimates for models are outside the bounds defined by 'bound.matrix'." )
    return( final.results )
}

ConvertManyToRate_ModelAve <- function(hisse.results, rate.param, which.element) {
    storage.matrix <- matrix(nrow=dim(hisse.results[[1]][[which.element]])[1], ncol=0)
    for (i in sequence(length(hisse.results))) {
        rate.vector <- hisse.results[[i]]$rates.mat[rate.param,]
        ## Some cases can have 'Inf' values for the rate.vector. This will produce NAs.
        ## Will happen in the case of problem with 0 speciation values. In the vector for extinction fraction.
        ## Need to check if the rate vector has some 'NA', 'NaN', or 'Inf'. Do something about it and throw a warning message.
        if( any( is.infinite( rate.vector ) ) ) stop( paste0("One or more values for ", rate.param, " are 'Inf'. Check model parameters." ) )
        storage.matrix <- cbind(storage.matrix, ConvertToRate(x=hisse.results[[i]][[which.element]], rate.vector=rate.vector))
    }
    return(storage.matrix)
}


ConvertManyToBinaryState <- function(hisse.results, which.element, AIC.weights=NULL) {
    if( is.null(AIC.weights) ){
        AIC.weights <- GetAICWeights(hisse.results)
    }
    storage.matrix <- matrix(nrow=dim(hisse.results[[1]][[which.element]])[1], ncol=0)
    for (i in sequence(length(hisse.results))) {
        storage.matrix <- cbind(storage.matrix, ConvertToBinaryState(x=hisse.results[[i]][[which.element]]))
    }
    final.results <- apply(storage.matrix, 1, weighted.mean, w=AIC.weights)
    return(final.results)
}

ConvertManyToMultiState <- function(hisse.results, which.element, AIC.weights=NULL) {
    if( is.null(AIC.weights) ){
        AIC.weights <- GetAICWeights(hisse.results)
    }
    storage.array <- array(dim=c(nrow(hisse.results[[1]][[which.element]]), length(unique(gsub("[^0-9]", "", colnames(hisse.results[[1]][[which.element]][,-1])))), length(hisse.results)))
    for (i in sequence(length(hisse.results))) {
        storage.array[,,i] <- as.matrix(MergeAcrossHidden(x=hisse.results[[i]][[which.element]]))
    }
    final.results <- apply(storage.array, c(1,2), weighted.mean, w=AIC.weights)
    rownames(final.results) <- rownames(hisse.results[[1]][[which.element]])
    colnames(final.results) <- unique(gsub("[^0-9]", "", colnames(hisse.results[[1]][[which.element]][,-1])))
    return(final.results)
}

GetAICWeights <- function(hisse.results, criterion="aic") {
    if(class(hisse.results)=="misse.states") {
        hisse.results <- list(hisse.results)
    }
    AIC.vector <- sapply(hisse.results, "[[", criterion)
    delta.AIC.vector <- AIC.vector - min(AIC.vector)
    rel.likelihood <- exp(-0.5 * delta.AIC.vector)
    AIC.weight.vector <- rel.likelihood / sum(rel.likelihood)
    return(AIC.weight.vector)
}


GetRateRange <- function(x, rate.param) {
    hisse.results <- x
    if(class(hisse.results)=="hisse.states") { #we have to make a list so we can run this generally
        tmp.list <- list()
        tmp.list[[1]] <- hisse.results
        hisse.results <- tmp.list
    }
    all.rates <- sapply(hisse.results, "[[", "rates.mat", simplify=FALSE)
    return(range(unname(unlist(sapply(all.rates, GetRelevantRowEntries, rate.param="turnover")))))
}


GetRelevantRowEntries <- function(x, rate.param) {
    return(x[which(rownames(x)==rate.param),])
}
