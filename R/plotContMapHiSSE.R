## IMPORTANT NOTE:
## These functions are modified from the package "phytools".
## These have been adapted to this package.
## For general use, besides internal usage in package "hisse", please use and cite "phytools".

plot.contMapHisse <- function(x, fsize=1, outline=TRUE, lwd=4, type="phylogram", mar=rep(0.3,4), direction="rightwards", offset=NULL, xlim=NULL, ylim=NULL, hold=TRUE, swap.underscore=TRUE, show.tiplabels=TRUE){
    
    ## This will plot the cont map for a HiSSE model.
    lims <- x$lims
    tree <- x$tree
    cols <- x$cols
    
    H <- nodeHeights(tree)

    if(is.null(fsize)){
        fsize <- c(1,1)
    }
    if(length(fsize)==1){
        fsize <- rep(fsize,2)
    }
    if(length(lwd)==1){
        lwd <- rep(lwd,2)
    } else if(length(lwd)>2){
        lwd <- lwd[1:2]
    }
    ## done optional arguments

    ## From here on the code comes from the "plot.densityMap" function.
    if(hold){
        set.null <- dev.hold()
    }
    
    if(type=="phylogram"){
        N<-length(tree$tip.label)
        if(is.null(ylim)){
            ylim <- NULL
        }

        ## Make the plot.
        plotSimmapHiSSE(tree, cols, pts=FALSE, lwd=lwd[1], fsize=fsize[1], mar=mar, add=outline,
                        xlim=xlim, ylim=ylim, direction=direction, offset=offset, hold=FALSE
                      , swap.underscore=swap.underscore, show.tiplabels=show.tiplabels)
        
    } else if(type=="fan"){
        invisible(
            capture.output(
                plotSimmapHiSSE(tree,cols,lwd=lwd[1],
                                mar=mar,fsize=fsize[1],add=outline,
                                type="fan",xlim=xlim,ylim=ylim,hold=FALSE,
                                swap.underscore=swap.underscore, show.tiplabels=show.tiplabels)
            )
        )
    }
    
    if(hold){
        set.null <- dev.flush()
    }
}

## Define the version of the plotSimmap function:
plotSimmapHiSSE <- function(tree,colors, fsize=1.0, lwd=2, pts=FALSE, node.numbers=FALSE
                          , mar=rep(0.1,4), add=FALSE, offset=NULL, direction="rightwards"
                          , type="phylogram", setEnv=TRUE, part=1.0, xlim=NULL, ylim=NULL
                          , nodes="intermediate", tips=NULL, maxY=NULL, hold=TRUE
                          , lend=2, asp=NA, plot=TRUE
                          , swap.underscore=TRUE, show.tiplabels){

    if(!inherits(tree,"phylo")) stop("tree should be object of class \"phylo\"")
    if(is.null(tree$maps)) stop("tree should contain mapped states on edges.")
    ftype <- "reg"
    
    ## swap out "_" character for spaces (assumes _ is a place holder)
    if( swap.underscore ){
        tree$tip.label<-gsub("_"," ",tree$tip.label)
    }
    
    if(type=="phylogram"){
        plotPhylogramHiSSE(tree = tree, colors = colors, fsize = fsize, ftype = ftype, lwd = lwd
                         , pts = pts,node.numbers = node.numbers, mar = mar,add = add, offset = offset
                         , direction = direction, xlim = xlim, ylim = ylim, placement = nodes
                         , tips = tips,lend = lend, asp = asp, plot = plot,show.tiplabels = show.tiplabels)
    } else if(type=="fan"){
        plotFanHiSSE(tree = tree, colors =  colors, fsize = fsize, ftype = ftype, lwd = lwd, mar = mar
                   , add = add, part = part, xlim = xlim, ylim = ylim, tips = tips
                   , maxY = maxY, lend = lend, plot = plot, show.tiplabels = show.tiplabels)
    }

    if(hold){
        set.null <- dev.flush()
    }
}

## #################################
## Some helping functions:

reorderSimmapHiSSE <- function(tree, order="cladewise", index.only=FALSE, ...){
    ii<-reorder.phylo(tree,order,index.only=TRUE,...)
    if(!index.only){
        if(inherits(ii,"phylo")) ii<-whichorder(ii$edge[,2],tree$edge[,2]) ## bug workaround
        tree$edge<-tree$edge[ii,]
        tree$edge.length<-tree$edge.length[ii]
        if(!is.null(tree$maps)){
            tree$maps<-tree$maps[ii]
            tree$mapped.edge<-tree$mapped.edge[ii,]
        }
        attr(tree,"order")<-order
        return(tree)
    } else return(ii)
}

## #################################

## #################################
## The rest of the plotting functions:

plotPhylogramHiSSE <- function(tree, colors, fsize, ftype, lwd, pts, node.numbers, mar, 
                               add, offset, direction, xlim, ylim, placement,
                               tips, lend, asp, plot, show.tiplabels){
    ## set offset fudge (empirically determined)
    offsetFudge <- 1.37
    ## reorder
    cw <- reorderSimmapHiSSE(tree)
    pw <- reorderSimmapHiSSE(tree, "postorder")
    ## count nodes and tips
    n<-Ntip(cw)
    m<-cw$Nnode
    ## Y coordinates for nodes
    Y<-matrix(NA, m+n, 1)
    ## first, assign y coordinates to all the tip nodes
    if(is.null(tips)){
        Y[cw$edge[cw$edge[,2]<=n,2]] <- 1:n
    } else{
        Y[cw$edge[cw$edge[,2]<=n,2]] <- ifelse(is.null(names(tips)), tips[sapply(1:Ntip(cw),function(x,y) which(y==x),y=cw$edge[cw$edge[,2]<=n,2])], tips[gsub(" ","_",cw$tip.label)] )
    }
    ## get Y coordinates of the nodes
    nodes <- unique(pw$edge[,1])
    for(i in 1:m){
        if(placement=="intermediate"){ 
            desc<-pw$edge[which(pw$edge[,1]==nodes[i]),2]
            Y[nodes[i]]<-(min(Y[desc])+max(Y[desc]))/2
        } else if(placement=="centered"){
            desc<-getDescendants(pw,nodes[i])
            desc<-desc[desc<=Ntip(pw)]
            Y[nodes[i]]<-(min(Y[desc])+max(Y[desc]))/2
        } else if(placement=="weighted"){
            desc<-pw$edge[which(pw$edge[,1]==nodes[i]),2]
            n1<-desc[which(Y[desc]==min(Y[desc]))]
            n2<-desc[which(Y[desc]==max(Y[desc]))]
            v1<-pw$edge.length[which(pw$edge[,2]==n1)]
            v2<-pw$edge.length[which(pw$edge[,2]==n2)]
            Y[nodes[i]]<-((1/v1)*Y[n1]+(1/v2)*Y[n2])/(1/v1+1/v2)
        } else if(placement=="inner"){
            desc<-getDescendants(pw,nodes[i])
            desc<-desc[desc<=Ntip(pw)]
            mm<-which(abs(Y[desc]-median(Y[1:Ntip(pw)]))==min(abs(Y[desc]-
                                                                  median(Y[1:Ntip(pw)]))))
            if(length(mm>1)) mm<-mm[which(Y[desc][mm]==min(Y[desc][mm]))]
            Y[nodes[i]]<-Y[desc][mm]
        }
    }
    ## compute node heights
    H<-nodeHeights(cw)
    ## open plot
    par(mar=mar)
    if(is.null(offset)) offset<-0.2*lwd/3+0.2/3
    if(!add) plot.new()
    ## ##
    if(is.null(xlim)){
        pp<-par("pin")[1]
        sw<-fsize*(max(strwidth(cw$tip.label,units="inches")))+
            offsetFudge*fsize*strwidth("W",units="inches")
        alp<-optimize(function(a,H,sw,pp) (a*1.04*max(H)+sw-pp)^2,H=H,sw=sw,pp=pp,
                      interval=c(0,1e6))$minimum
        xlim<-if(direction=="leftwards") c(min(H)-sw/alp,max(H)) else c(min(H),max(H)+sw/alp)
    }
    if(is.null(ylim)) ylim=range(Y)
    if(direction=="leftwards") H<-max(H)-H
    plot.window(xlim=xlim,ylim=ylim,asp=asp)
    if(plot){
        for(i in 1:m) lines(H[which(cw$edge[,1]==nodes[i]),1], Y[cw$edge[which(cw$edge[,1]==nodes[i]),2]], col=colors[names(cw$maps[[match(nodes[i],cw$edge[,1])]])[1]],lwd=lwd)
        for(i in 1:nrow(cw$edge)){
            x<-H[i,1]
            for(j in 1:length(cw$maps[[i]])){
                if(direction=="leftwards")
                    lines(c(x,x-cw$maps[[i]][j]),c(Y[cw$edge[i,2]],Y[cw$edge[i,2]]),
                          col=colors[names(cw$maps[[i]])[j]],lwd=lwd,lend=lend)
                else lines(c(x,x+cw$maps[[i]][j]),c(Y[cw$edge[i,2]],Y[cw$edge[i,2]]),
                           col=colors[names(cw$maps[[i]])[j]],lwd=lwd,lend=lend)
                if(pts) points(c(x,x+cw$maps[[i]][j]),c(Y[cw$edge[i,2]],Y[cw$edge[i,2]]),
                               pch=20,lwd=(lwd-1))
                x<-x+if(direction=="leftwards") -cw$maps[[i]][j] else cw$maps[[i]][j]
                j<-j+1
            }
        }
        if(node.numbers){
            symbols(if(direction=="leftwards") max(H) else 0,
                    mean(Y[cw$edge[cw$edge[,1]==(Ntip(cw)+1),2]]),
                    rectangles=matrix(c(1.2*fsize*strwidth(as.character(Ntip(cw)+1)),
                                        1.4*fsize*strheight(as.character(Ntip(cw)+1))),1,2),inches=FALSE,
                    bg="white",add=TRUE)
            text(if(direction=="leftwards") max(H) else 0,
                 mean(Y[cw$edge[cw$edge[,1]==(Ntip(cw)+1),2]]),Ntip(cw)+1,
                 cex=fsize)
            for(i in 1:nrow(cw$edge)){
                x<-H[i,2]
                if(cw$edge[i,2]>Ntip(cw)){
                    symbols(x,Y[cw$edge[i,2]],
                            rectangles=matrix(c(1.2*fsize*strwidth(as.character(cw$edge[i,2])),
						1.4*fsize*strheight(as.character(cw$edge[i,2]))),1,2),inches=FALSE,
                            bg="white",add=TRUE)
                    text(x,Y[cw$edge[i,2]],cw$edge[i,2],cex=fsize)
                }
            }
        }
        if(direction=="leftwards") pos<-if(par()$usr[1]>par()$usr[2]) 4 else 2
        if(direction=="rightwards") pos<-if(par()$usr[1]>par()$usr[2]) 2 else 4

        ## Option to plot the tiplabels
        if( show.tiplabels ){
            for(i in 1:n) if(ftype) text(H[which(cw$edge[,2]==i),2],Y[i],cw$tip.label[i],pos=pos,
                                         offset=offset,cex=fsize,font=ftype)
        }
    }
}

plotFanHiSSE <- function(tree,colors,fsize,ftype,lwd,mar,add,part
                        ,xlim,ylim,tips,maxY,lend,plot,show.tiplabels){
    if(!plot) cat("plot=FALSE option is not permitted for type=\"fan\". Tree will be plotted.\n")
                                        # reorder
    cw<-reorder(tree)
    pw<-reorder(tree,"pruningwise")
                                        # count nodes and tips
    n<-Ntip(cw)
    m<-cw$Nnode 
                                        # get Y coordinates on uncurved space
    Y<-vector(length=m+n)
    if(is.null(tips)) tips<-1:n
    if(part<1.0) Y[cw$edge[cw$edge[,2]<=n,2]]<-0:(n-1)
    else Y[cw$edge[cw$edge[,2]<=n,2]]<-tips
    nodes<-unique(pw$edge[,1])
    for(i in 1:m){
        desc<-pw$edge[which(pw$edge[,1]==nodes[i]),2]
        Y[nodes[i]]<-(min(Y[desc])+max(Y[desc]))/2
    }
    if(is.null(maxY)) maxY<-max(Y)
    Y<-setNames(Y/maxY*2*pi,1:(n+m))
    Y<-part*cbind(Y[as.character(cw$edge[,2])],Y[as.character(cw$edge[,2])])
    R<-nodeHeights(cw)
                                        # now put into a circular coordinate system
    x<-R*cos(Y)
    y<-R*sin(Y)
                                        # optimize x & y limits
    par(mar=mar)
    offsetFudge<-1.37 # empirically determined
    offset<-0
    pp<-par("pin")[1]
    sw<-fsize*(max(strwidth(cw$tip.label,units="inches")))+
        offsetFudge*offset*fsize*strwidth("W",units="inches") 
    alp<-optimize(function(a,H,sw,pp) (2*a*1.04*max(H)+2*sw-pp)^2,H=R,sw=sw,pp=pp,
                  interval=c(0,1e6))$minimum
    if(part<=0.25) x.lim<-y.lim<-c(0,max(R)+sw/alp)
    else if(part>0.25&&part<=0.5){ 
        x.lim<-c(-max(R)-sw/alp,max(R)+sw/alp)
        y.lim<-c(0,max(R)+sw/alp)
    } else x.lim<-y.lim<-c(-max(R)-sw/alp,max(R)+sw/alp)
    if(is.null(xlim)) xlim<-x.lim
    if(is.null(ylim)) ylim<-y.lim
                                        # plot tree
    if(!add) plot.new()
    plot.window(xlim=xlim,ylim=ylim,asp=1)
                                        # plot radial lines (edges)
    ## first, the lines emerging from the root (if there are only two):
    jj<-which(cw$edge[,1]==(Ntip(cw)+1))
    if(length(jj)==2){
        m.left<-cumsum(cw$maps[[jj[1]]])/sum(cw$maps[[jj[1]]])
        xx.left<-c(x[jj[1],1],x[jj[1],1]+(x[jj[1],2]-x[jj[1],1])*m.left)
        yy.left<-c(y[jj[1],1],y[jj[1],1]+(y[jj[1],2]-y[jj[1],1])*m.left)
        m.right<-cumsum(cw$maps[[jj[2]]])/sum(cw$maps[[jj[2]]])
        xx.right<-c(x[jj[2],1],x[jj[2],1]+(x[jj[2],2]-x[jj[2],1])*m.right)
        yy.right<-c(y[jj[2],1],y[jj[2],1]+(y[jj[2],2]-y[jj[2],1])*m.right)
        xx<-c(xx.left[length(xx.left):1],xx.right[2:length(xx.right)])
        yy<-c(yy.left[length(yy.left):1],yy.right[2:length(yy.right)])
        col<-colors[c(names(m.left)[length(m.left):1],names(m.right))]
        segments(xx[2:length(xx)-1],yy[2:length(yy)-1],xx[2:length(xx)],yy[2:length(yy)],
                 col=col,lwd=lwd,lend=lend)
    } else jj<-NULL
    for(i in 1:nrow(cw$edge)){
        if(i%in%jj==FALSE){
            maps<-cumsum(cw$maps[[i]])/sum(cw$maps[[i]])
            xx<-c(x[i,1],x[i,1]+(x[i,2]-x[i,1])*maps)
            yy<-c(y[i,1],y[i,1]+(y[i,2]-y[i,1])*maps)
            for(i in 1:(length(xx)-1)) lines(xx[i+0:1],yy[i+0:1],col=colors[names(maps)[i]],
                                             lwd=lwd,lend=lend)
        }
    }
                                        # plot circular lines
    for(i in 1:m+n){
        r<-R[match(i,cw$edge)]
        a1<-min(Y[which(cw$edge==i)])
        a2<-max(Y[which(cw$edge==i)])
        plotrix::draw.arc(0,0,r,a1,a2,lwd=lwd,col=colors[names(cw$maps[[match(i,cw$edge[,1])]])[1]])
    }
                                        # plot labels
    ## Option to plot the tiplabels
    if( show.tiplabels ){
        for(i in 1:n){
            ii<-which(cw$edge[,2]==i)
            aa<-Y[ii,2]/(2*pi)*360
            adj<-if(aa>90&&aa<270) c(1,0.25) else c(0,0.25)
            tt<-if(aa>90&&aa<270) paste(cw$tip.label[i]," ",sep="") else paste(" ",
                                                                               cw$tip.label[i],sep="")
            aa<-if(aa>90&&aa<270) 180+aa else aa
            if(ftype) text(x[ii,2],y[ii,2],tt,srt=aa,adj=adj,cex=fsize,font=ftype)
        }
    }
}
