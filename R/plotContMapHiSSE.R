## IMPORTANT NOTE:
## These functions are modified from the package "phytools".
## These have been adapted to this package.
## For general use, besides internal usage in package "hisse", please use and cite "phytools".

plot.contMapHisse <- function(A, B, lwd.factor = 0.5, fsize=1, ftype="reg", add=FALSE, lwd=4, type="phylogram", mar=rep(0.3,4), direction="rightwards", offset=NULL, xlim=NULL, ylim=NULL, hold=TRUE, swap.underscore=TRUE, outline = TRUE, outline.color = "black", show.tiplabels=TRUE){
    
    ## This will plot the cont map for a HiSSE model.
    ## A$tree
    ## A$cols
    
    H <- nodeHeights(A$tree)

    ## Translate the ftype to numeric:
    ## This follows the same nomenclature as in the phytools package.
    ftype <- which(c("reg","b","i","bi") == ftype)    

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
        N<-length(A$tree$tip.label) ## Same tree for both plots.
        if(is.null(ylim)){
            ylim <- NULL
        }

        ## Make the plot.
        plotSimmapHiSSE(treeA = A$tree, colorsA = A$cols, treeB = B$tree, colorsB = B$cols
                      , pts=FALSE, lwd=lwd[1], lwd.factor=lwd.factor, fsize=fsize[1], ftype=ftype
                      , mar=mar, add=add,xlim=xlim, ylim=ylim, direction=direction
                      , offset=offset, hold=FALSE, swap.underscore=swap.underscore
                      , show.tiplabels=show.tiplabels, outline = outline, outline.color=outline.color)
        
    } else if(type=="fan"){
        invisible(
            capture.output(
                plotSimmapHiSSE(treeA = A$tree, colorsA = A$cols, treeB = B$tree, colorsB = B$cols
                      , pts=FALSE, lwd=lwd[1], lwd.factor=lwd.factor, fsize=fsize[1], ftype=ftype
                      , mar=mar, add=add,xlim=xlim, ylim=ylim, direction=direction
                      , offset=offset, hold=FALSE, swap.underscore=swap.underscore
                      , show.tiplabels=show.tiplabels, outline = outline, outline.color=outline.color
                      ,  type="fan")
            )
        )
    }
    
    if(hold){
        set.null <- dev.flush()
    }
}

## Define the version of the plotSimmap function:
plotSimmapHiSSE <- function(treeA, colorsA, treeB, colorsB, fsize=1.0, ftype="reg", lwd=2
                          , pts=FALSE, lwd.factor=0.5
                          , mar=rep(0.1,4), add=FALSE, offset=NULL, direction="rightwards"
                          , type="phylogram", setEnv=TRUE, part=1.0, xlim=NULL, ylim=NULL
                          , nodes="intermediate", tips=NULL, maxY=NULL, hold=TRUE
                          , lend=2, asp=NA, swap.underscore=TRUE, show.tiplabels, outline, outline.color){
    ## swap out "_" character for spaces (assumes _ is a place holder)
    if( swap.underscore ){
        treeA$tip.label <- gsub("_"," ",treeA$tip.label)
        treeB$tip.label <- gsub("_"," ",treeB$tip.label)
    }
    
    if(type=="phylogram"){
        plotPhylogramHiSSE(treeA = treeA, colorsA = colorsA, treeB = treeB, colorsB = colorsB
                         , fsize = fsize, ftype = ftype, lwd = lwd, lwd.factor = lwd.factor
                         , pts = pts, mar = mar,add = add, offset = offset
                         , outline = outline, outline.color = outline.color
                         , direction = direction, xlim = xlim, ylim = ylim, placement = nodes
                         , tips = tips,lend = lend, asp = asp, show.tiplabels = show.tiplabels)
    } else if(type=="fan"){
        plotFanHiSSE(treeA = treeA, colorsA = colorsA, treeB = treeB, colorsB = colorsB
                   , fsize = fsize, ftype = ftype, lwd = lwd, lwd.factor = lwd.factor, mar = mar
                   , outline = outline, outline.color = outline.color
                   , add = add, part = part, xlim = xlim, ylim = ylim, tips = tips
                   , maxY = maxY, lend = lend, show.tiplabels = show.tiplabels)
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

plotPhylogramHiSSE <- function(treeA, colorsA, treeB, colorsB, fsize, ftype, lwd, lwd.factor,
                               pts, mar, add, offset, direction, xlim, ylim, placement,
                               tips, lend, asp, show.tiplabels, outline, outline.color){
    
    ## set offset fudge (empirically determined)
    offsetFudge <- 1.37
    ## reorder
    treeA <- reorderSimmapHiSSE(treeA)
    treeB <- reorderSimmapHiSSE(treeB)
    ptreeA <- reorderSimmapHiSSE(treeA, "postorder")
    ptreeB <- reorderSimmapHiSSE(treeB, "postorder")

    ## Prepare some objects dependent on the presence of an outline:
    if( outline ){
        ## The outline lwd need to be the base 'lwd' and the larger.
        lwdA <- lwd
        lwd <- lwd * 1.5
        lwdB <- lwdA * lwd.factor
    } else{
        lwdA <- lwd
        lwdB <- lwdA * lwd.factor
    }
    ## count nodes and tips
    n<-Ntip(treeA)
    m<-treeA$Nnode
    ## Y coordinates for nodes
    Y<-matrix(NA, m+n, 1)
    ## first, assign y coordinates to all the tip nodes
    if(is.null(tips)){
        Y[treeA$edge[treeA$edge[,2]<=n,2]] <- 1:n
    } else{
        Y[treeA$edge[treeA$edge[,2]<=n,2]] <- ifelse(is.null(names(tips)), tips[sapply(1:Ntip(treeA),function(x,y) which(y==x),y=treeA$edge[treeA$edge[,2]<=n,2])], tips[gsub(" ","_",treeA$tip.label)] )
    }
    ## get Y coordinates of the nodes
    nodes <- unique(ptreeA$edge[,1])
    for(i in 1:m){
        if(placement=="intermediate"){ 
            desc<-ptreeA$edge[which(ptreeA$edge[,1]==nodes[i]),2]
            Y[nodes[i]]<-(min(Y[desc])+max(Y[desc]))/2
        } else if(placement=="centered"){
            desc<-getDescendants(ptreeA,nodes[i])
            desc<-desc[desc<=Ntip(ptreeA)]
            Y[nodes[i]]<-(min(Y[desc])+max(Y[desc]))/2
        } else if(placement=="weighted"){
            desc<-ptreeA$edge[which(ptreeA$edge[,1]==nodes[i]),2]
            n1<-desc[which(Y[desc]==min(Y[desc]))]
            n2<-desc[which(Y[desc]==max(Y[desc]))]
            v1<-ptreeA$edge.length[which(ptreeA$edge[,2]==n1)]
            v2<-ptreeA$edge.length[which(ptreeA$edge[,2]==n2)]
            Y[nodes[i]]<-((1/v1)*Y[n1]+(1/v2)*Y[n2])/(1/v1+1/v2)
        } else if(placement=="inner"){
            desc<-getDescendants(ptreeA,nodes[i])
            desc<-desc[desc<=Ntip(ptreeA)]
            mm<-which(abs(Y[desc]-median(Y[1:Ntip(ptreeA)]))==min(abs(Y[desc]-
                                                                  median(Y[1:Ntip(ptreeA)]))))
            if(length(mm>1)) mm<-mm[which(Y[desc][mm]==min(Y[desc][mm]))]
            Y[nodes[i]]<-Y[desc][mm]
        }
    }
    ## compute node heights
    H<-nodeHeights(treeA)
    ## open plot
    par(mar=mar)
    if(is.null(offset)) offset<-0.2*lwd/3+0.2/3
    if(!add) plot.new()
    
    ## If not showing the tip labels, then need to set the size of the font to zero.
    ## This will set the space for the labels to 0 also.
    if( !show.tiplabels ){
        fsize <- 0.0
    }
    
    if(is.null(xlim)){
        pp<-par("pin")[1]
        sw<-fsize*(max(strwidth(treeA$tip.label,units="inches")))+
            offsetFudge*fsize*strwidth("W",units="inches")
        alp<-optimize(function(a,H,sw,pp) (a*1.04*max(H)+sw-pp)^2,H=H,sw=sw,pp=pp,
                      interval=c(0,1e6))$minimum
        xlim<-if(direction=="leftwards") c(min(H)-sw/alp,max(H)) else c(min(H),max(H)+sw/alp)
    }
    if(is.null(ylim)) ylim=range(Y)
    if(direction=="leftwards") H<-max(H)-H
    plot.window(xlim=xlim,ylim=ylim,asp=asp)

    ## This is the block that make the plot.
    ## We need to make the first plot and then the second one. For this we just need to duplicate this block.

    ## Plot the outline of the tree, if necessary:
    if( outline ){
        ## Make the color black
        treeOut <- treeA ## Cannot modify the original tree.
        for(i in 1:length(treeOut$maps)){
            names(treeOut$maps[[i]]) <- c("1")
        }
        
        ## Make the plot:
        for(i in 1:m){
            lines(H[which(treeOut$edge[,1]==nodes[i]),1], Y[treeOut$edge[which(treeOut$edge[,1]==nodes[i]),2]]
                , col=outline.color, lwd=lwd)
        }
        for(i in 1:nrow(treeOut$edge)){
            x<-H[i,1]
            for(j in 1:length(treeOut$maps[[i]])){
                if(direction=="leftwards")
                    lines(c(x,x-treeOut$maps[[i]][j]),c(Y[treeOut$edge[i,2]],Y[treeOut$edge[i,2]]),
                          col=outline.color, lwd=lwd, lend=lend)
                else lines(c(x,x+treeOut$maps[[i]][j]),c(Y[treeOut$edge[i,2]],Y[treeOut$edge[i,2]]),
                           col=outline.color, lwd=lwd, lend=lend)
                if(pts) points(c(x,x+treeOut$maps[[i]][j]),c(Y[treeOut$edge[i,2]],Y[treeOut$edge[i,2]]),
                               pch=20,lwd=(lwd-1))
                x<-x+if(direction=="leftwards") -treeOut$maps[[i]][j] else treeOut$maps[[i]][j]
                j<-j+1
            }
        }
    }
    
    ## Plot the first tree layer (treeA)
    for(i in 1:m){
        lines(H[which(treeA$edge[,1]==nodes[i]),1], Y[treeA$edge[which(treeA$edge[,1]==nodes[i]),2]], col=colorsA[names(treeA$maps[[match(nodes[i],treeA$edge[,1])]])[1]],lwd=lwdA)
    }
    for(i in 1:nrow(treeA$edge)){
        x<-H[i,1]
        for(j in 1:length(treeA$maps[[i]])){
            if(direction=="leftwards")
                lines(c(x,x-treeA$maps[[i]][j]),c(Y[treeA$edge[i,2]],Y[treeA$edge[i,2]]),
                      col=colorsA[names(treeA$maps[[i]])[j]],lwd=lwdA,lend=lend)
            else lines(c(x,x+treeA$maps[[i]][j]),c(Y[treeA$edge[i,2]],Y[treeA$edge[i,2]]),
                       col=colorsA[names(treeA$maps[[i]])[j]],lwd=lwdA,lend=lend)
            if(pts) points(c(x,x+treeA$maps[[i]][j]),c(Y[treeA$edge[i,2]],Y[treeA$edge[i,2]]),
                           pch=20,lwd=(lwdA-1))
            x<-x+if(direction=="leftwards") -treeA$maps[[i]][j] else treeA$maps[[i]][j]
            j<-j+1
        }
    }

    ## Plot the second tree layer (treeB)
    for(i in 1:m){
        lines(H[which(treeB$edge[,1]==nodes[i]),1], Y[treeB$edge[which(treeB$edge[,1]==nodes[i]),2]], col=colorsB[names(treeB$maps[[match(nodes[i],treeB$edge[,1])]])[1]],lwd=lwdB)
    }
    for(i in 1:nrow(treeB$edge)){
        x<-H[i,1]
        for(j in 1:length(treeB$maps[[i]])){
            if(direction=="leftwards")
                lines(c(x,x-treeB$maps[[i]][j]),c(Y[treeB$edge[i,2]],Y[treeB$edge[i,2]]),
                      col=colorsB[names(treeB$maps[[i]])[j]],lwd=lwdB,lend=lend)
            else lines(c(x,x+treeB$maps[[i]][j]),c(Y[treeB$edge[i,2]],Y[treeB$edge[i,2]]),
                       col=colorsB[names(treeB$maps[[i]])[j]],lwd=lwdB,lend=lend)
            if(pts) points(c(x,x+treeB$maps[[i]][j]),c(Y[treeB$edge[i,2]],Y[treeB$edge[i,2]]),
                           pch=20,lwd=(lwdB-1))
            x<-x+if(direction=="leftwards") -treeB$maps[[i]][j] else treeB$maps[[i]][j]
            j<-j+1
        }
    }

    if(direction=="leftwards") pos<-if(par()$usr[1]>par()$usr[2]) 4 else 2
    if(direction=="rightwards") pos<-if(par()$usr[1]>par()$usr[2]) 2 else 4

    ## Option to plot the tiplabels
    if( show.tiplabels ){
        for(i in 1:n) text(H[which(treeA$edge[,2]==i),2],Y[i],treeA$tip.label[i],pos=pos,
                           offset=offset,cex=fsize,font=ftype)
    }
}

plotFanHiSSE <- function(treeA, colorsA, treeB, colorsB, fsize, ftype, lwd, lwd.factor
                       , mar, add, part, xlim, ylim, tips, maxY, lend, show.tiplabels
                       , outline, outline.color){

    ## Reorder
    treeA <- reorderSimmapHiSSE(treeA)
    treeB <- reorderSimmapHiSSE(treeB)
    ptreeA <- reorderSimmapHiSSE(treeA, "postorder")
    ptreeB <- reorderSimmapHiSSE(treeB, "postorder")

    ## Prepare some objects dependent on the presence of an outline:
    if( outline ){
        ## The outline lwd need to be the base 'lwd' and the larger.
        lwdA <- lwd
        lwd <- lwd * 1.5
        lwdB <- lwdA * lwd.factor
    } else{
        lwdA <- lwd
        lwdB <- lwdA * lwd.factor
    }

    ## Adjust the size of the font if not showing the labels:
    if( !show.tiplabels ){
        fsize <- 0.0
    }
    
    ## count nodes and tips
    n<-Ntip(treeA)
    m<-treeA$Nnode 
    ## get Y coordinates on uncurved space
    Y<-vector(length=m+n)
    if(is.null(tips)) tips<-1:n
    if(part<1.0) Y[treeA$edge[treeA$edge[,2]<=n,2]]<-0:(n-1)
    else Y[treeA$edge[treeA$edge[,2]<=n,2]]<-tips
    nodes<-unique(ptreeA$edge[,1])
    for(i in 1:m){
        desc<-ptreeA$edge[which(ptreeA$edge[,1]==nodes[i]),2]
        Y[nodes[i]]<-(min(Y[desc])+max(Y[desc]))/2
    }
    if(is.null(maxY)) maxY<-max(Y)
    Y<-setNames(Y/maxY*2*pi,1:(n+m))
    Y<-part*cbind(Y[as.character(treeA$edge[,2])],Y[as.character(treeA$edge[,2])])
    R<-nodeHeights(treeA)
    ## now put into a circular coordinate system
    x<-R*cos(Y)
    y<-R*sin(Y)
    ## optimize x & y limits
    par(mar=mar)
    offsetFudge<-1.37 ## empirically determined
    offset <- 1
    pp<-par("pin")[1]
    sw<-fsize*(max(strwidth(treeA$tip.label,units="inches")))+
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
    ## plot tree
    if(!add) plot.new()
    plot.window(xlim=xlim,ylim=ylim,asp=1)

    ## Make the plots:
    ## Plot the outline of the tree, if necessary:
    if( outline ){
        ## Make the color black
        treeOut <- treeA ## Cannot modify the original tree.
        for(i in 1:length(treeOut$maps)){
            names(treeOut$maps[[i]]) <- c("1")
        }

        ## plot radial lines (edges)
        ## first, the lines emerging from the root (if there are only two):
        jj<-which(treeOut$edge[,1]==(Ntip(treeOut)+1))
        if(length(jj)==2){
            m.left<-cumsum(treeOut$maps[[jj[1]]])/sum(treeOut$maps[[jj[1]]])
            xx.left<-c(x[jj[1],1],x[jj[1],1]+(x[jj[1],2]-x[jj[1],1])*m.left)
            yy.left<-c(y[jj[1],1],y[jj[1],1]+(y[jj[1],2]-y[jj[1],1])*m.left)
            m.right<-cumsum(treeOut$maps[[jj[2]]])/sum(treeOut$maps[[jj[2]]])
            xx.right<-c(x[jj[2],1],x[jj[2],1]+(x[jj[2],2]-x[jj[2],1])*m.right)
            yy.right<-c(y[jj[2],1],y[jj[2],1]+(y[jj[2],2]-y[jj[2],1])*m.right)
            xx<-c(xx.left[length(xx.left):1],xx.right[2:length(xx.right)])
            yy<-c(yy.left[length(yy.left):1],yy.right[2:length(yy.right)])
            segments(xx[2:length(xx)-1],yy[2:length(yy)-1],xx[2:length(xx)],yy[2:length(yy)],
                     col=outline.color,lwd=lwd,lend=lend)
        } else jj<-NULL
        for(i in 1:nrow(treeOut$edge)){
            if(i%in%jj==FALSE){
                maps<-cumsum(treeOut$maps[[i]])/sum(treeOut$maps[[i]])
                xx<-c(x[i,1],x[i,1]+(x[i,2]-x[i,1])*maps)
                yy<-c(y[i,1],y[i,1]+(y[i,2]-y[i,1])*maps)
                for(i in 1:(length(xx)-1)) lines(xx[i+0:1],yy[i+0:1],col=outline.color,
                                                 lwd=lwd,lend=lend)
            }
        }
        ## plot circular lines
        for(i in 1:m+n){
            r<-R[match(i,treeOut$edge)]
            a1<-min(Y[which(treeOut$edge==i)])
            a2<-max(Y[which(treeOut$edge==i)])
            draw.arc(0,0,r,a1,a2,lwd=lwd,col=outline.color)
        }
        
    }

    ## Make the plot for the first layer:
    ## plot radial lines (edges)
    ## first, the lines emerging from the root (if there are only two):
    jj<-which(treeA$edge[,1]==(Ntip(treeA)+1))
    if(length(jj)==2){
        m.left<-cumsum(treeA$maps[[jj[1]]])/sum(treeA$maps[[jj[1]]])
        xx.left<-c(x[jj[1],1],x[jj[1],1]+(x[jj[1],2]-x[jj[1],1])*m.left)
        yy.left<-c(y[jj[1],1],y[jj[1],1]+(y[jj[1],2]-y[jj[1],1])*m.left)
        m.right<-cumsum(treeA$maps[[jj[2]]])/sum(treeA$maps[[jj[2]]])
        xx.right<-c(x[jj[2],1],x[jj[2],1]+(x[jj[2],2]-x[jj[2],1])*m.right)
        yy.right<-c(y[jj[2],1],y[jj[2],1]+(y[jj[2],2]-y[jj[2],1])*m.right)
        xx<-c(xx.left[length(xx.left):1],xx.right[2:length(xx.right)])
        yy<-c(yy.left[length(yy.left):1],yy.right[2:length(yy.right)])
        col<-colorsA[c(names(m.left)[length(m.left):1],names(m.right))]
        segments(xx[2:length(xx)-1],yy[2:length(yy)-1],xx[2:length(xx)],yy[2:length(yy)],
                 col=col,lwd=lwdA,lend=lend)
    } else jj<-NULL
    for(i in 1:nrow(treeA$edge)){
        if(i%in%jj==FALSE){
            maps<-cumsum(treeA$maps[[i]])/sum(treeA$maps[[i]])
            xx<-c(x[i,1],x[i,1]+(x[i,2]-x[i,1])*maps)
            yy<-c(y[i,1],y[i,1]+(y[i,2]-y[i,1])*maps)
            for(i in 1:(length(xx)-1)) lines(xx[i+0:1],yy[i+0:1],col=colorsA[names(maps)[i]],
                                             lwd=lwdA,lend=lend)
        }
    }
    ## plot circular lines
    for(i in 1:m+n){
        r<-R[match(i,treeA$edge)]
        a1<-min(Y[which(treeA$edge==i)])
        a2<-max(Y[which(treeA$edge==i)])
        draw.arc(0,0,r,a1,a2,lwd=lwdA,col=colorsA[names(treeA$maps[[match(i,treeA$edge[,1])]])[1]])
    }

    ## Make the plot for the second layer:
    ## plot radial lines (edges)
    ## first, the lines emerging from the root (if there are only two):
    jj<-which(treeB$edge[,1]==(Ntip(treeB)+1))
    if(length(jj)==2){
        m.left<-cumsum(treeB$maps[[jj[1]]])/sum(treeB$maps[[jj[1]]])
        xx.left<-c(x[jj[1],1],x[jj[1],1]+(x[jj[1],2]-x[jj[1],1])*m.left)
        yy.left<-c(y[jj[1],1],y[jj[1],1]+(y[jj[1],2]-y[jj[1],1])*m.left)
        m.right<-cumsum(treeB$maps[[jj[2]]])/sum(treeB$maps[[jj[2]]])
        xx.right<-c(x[jj[2],1],x[jj[2],1]+(x[jj[2],2]-x[jj[2],1])*m.right)
        yy.right<-c(y[jj[2],1],y[jj[2],1]+(y[jj[2],2]-y[jj[2],1])*m.right)
        xx<-c(xx.left[length(xx.left):1],xx.right[2:length(xx.right)])
        yy<-c(yy.left[length(yy.left):1],yy.right[2:length(yy.right)])
        col<-colorsB[c(names(m.left)[length(m.left):1],names(m.right))]
        segments(xx[2:length(xx)-1],yy[2:length(yy)-1],xx[2:length(xx)],yy[2:length(yy)],
                 col=col,lwd=lwdB,lend=lend)
    } else jj<-NULL
    for(i in 1:nrow(treeB$edge)){
        if(i%in%jj==FALSE){
            maps<-cumsum(treeB$maps[[i]])/sum(treeB$maps[[i]])
            xx<-c(x[i,1],x[i,1]+(x[i,2]-x[i,1])*maps)
            yy<-c(y[i,1],y[i,1]+(y[i,2]-y[i,1])*maps)
            for(i in 1:(length(xx)-1)) lines(xx[i+0:1],yy[i+0:1],col=colorsB[names(maps)[i]],
                                             lwd=lwdB,lend=lend)
        }
    }
    ## plot circular lines
    for(i in 1:m+n){
        r<-R[match(i,treeB$edge)]
        a1<-min(Y[which(treeB$edge==i)])
        a2<-max(Y[which(treeB$edge==i)])
        draw.arc(0,0,r,a1,a2,lwd=lwdB,col=colorsB[names(treeB$maps[[match(i,treeB$edge[,1])]])[1]])
    }

    ## Plot the tip labels if necessary:    
    if( show.tiplabels ){
        for(i in 1:n){
            ii<-which(treeA$edge[,2]==i)
            aa<-Y[ii,2]/(2*pi)*360
            adj<-if(aa>90&&aa<270) c(1,0.25) else c(0,0.25)
            tt<-if(aa>90&&aa<270) paste(treeA$tip.label[i]," ",sep="") else paste(" ",
                                                                               treeA$tip.label[i],sep="")
            aa<-if(aa>90&&aa<270) 180+aa else aa
            text(x[ii,2],y[ii,2],tt,srt=aa,adj=adj,cex=fsize,font=ftype)
        }
    }
}

## function whichorder
## written by Liam Revell 2011, 2013, 2015
## Replicated here to avoid dependencies problems due to RCran policy.
## Very simple function.
whichorder <- function(x,y) sapply(x, function(x,y) which(x==y), y=y)
