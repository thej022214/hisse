SimulateHisse <- function(turnover.rates, eps.rates, transition.rates, max.taxa=Inf, max.t=Inf, max.wall.time=Inf, x0, nstart=1) {
	start <- Sys.time()
	if(length(turnover.rates) != length(eps.rates)) {
		stop("need to have same number of turnover and eps rates")	
	}
	if(length(turnover.rates) != dim(transition.rates)[1]) {
		stop("need to have same number of turnover and eps rates as rows in the transition matrix")	
	}
	if(length(turnover.rates) != dim(transition.rates)[1]) {
		stop("need to have same number of turnover and eps rates as columns in the transition matrix")	
	}
	state.levels <- c(0:(length(turnover.rates)-1))
	states <- factor(as.character(state.levels), levels=state.levels)
	results <- data.table(anc=NA, id=1, state=factor(x0, state.levels), length=0, height=0, living=TRUE, descendants=FALSE)
	for(additional.row in sequence(nstart-1)) {
		results <- rbind(results, data.table(anc=NA, id=additional.row+1, state=factor(x0, state.levels), length=0, height=0, living=TRUE, descendants=FALSE))
	}
	setkey(results, id, living)
	birth.rates <- GetBirthRate(turnover.rates, eps.rates)
	death.rates <- GetDeathRate(turnover.rates, eps.rates)
	diag(transition.rates) <- NA
	all.rates <- c(birth.rates, death.rates, c(t(transition.rates))) #c(t()) so that we have first row (from 0 to 1,2,3), then second row, etc, rather than by columns
	all.rates <- all.rates[!is.na(all.rates)] #drop the diagonals
	keep.running <- TRUE
	while(keep.running) {
		tip.state.counts <- table(subset(results, living)$state)
		birth.rates.actual <- birth.rates * tip.state.counts
		death.rates.actual <- death.rates * tip.state.counts
		transition.rates.actual <- apply(transition.rates, 2, Multiply, y=tip.state.counts)
		birth.wait.times <- suppressWarnings(rexp(n=length(birth.rates.actual), birth.rates.actual))
		death.wait.times <- suppressWarnings(rexp(n=length(death.rates.actual), death.rates.actual))
		transition.wait.times <- suppressWarnings(matrix(rexp(n=length(transition.rates.actual), transition.rates.actual), nrow=dim(transition.rates.actual)[1])) #the NAs will be an issue
		min.times <- suppressWarnings(c(min(birth.wait.times, na.rm=TRUE), min(death.wait.times, na.rm=TRUE), min(transition.wait.times, na.rm=TRUE)))
		if((min(min.times, na.rm=TRUE)+max(subset(results, living)$height)) > max.t) { #gone too long
			keep.running <- FALSE	
			time.dif <- max.t - max(subset(results, living)$height)
			results[which(results$living),]$height <- results[which(results$living),]$height + time.dif
			results[which(results$living),]$length <- results[which(results$living),]$length + time.dif
		}
		if(keep.running) {
			results[which(results$living),]$height <- results[which(results$living),]$height + min(min.times)
			results[which(results$living),]$length <- results[which(results$living),]$length + min(min.times)

			if(which.min(min.times)==1) { #birth
				potential.lucky.taxa <- subset(results, living & state==states[which.min(birth.wait.times)])$id
				lucky.taxon <- potential.lucky.taxa[sample.int(length(potential.lucky.taxa), 1)]
				results[which(id==lucky.taxon),]$living <- FALSE
				results[which(id==lucky.taxon),]$descendants <- TRUE
				results <- rbind(results, data.table(anc=lucky.taxon, id=max(results$id)+1, state=subset(results, id==lucky.taxon)$state, length=0, height=subset(results, id==lucky.taxon)$height, living=TRUE, descendants=FALSE))
				results <- rbind(results, data.table(anc=lucky.taxon, id=max(results$id)+1, state=subset(results, id==lucky.taxon)$state, length=0, height=subset(results, id==lucky.taxon)$height, living=TRUE, descendants=FALSE))
			}
			if(which.min(min.times)==2) { #death
				potential.unlucky.taxa <- subset(results, living & state==states[which.min(death.wait.times)])$id
				unlucky.taxon <- potential.unlucky.taxa[sample.int(length(potential.unlucky.taxa), 1)]
				results[which(id==lucky.taxon),]$living <- FALSE
			}
			if(which.min(min.times)==3) { #transition
				from.to <- which(transition.wait.times == min(transition.wait.times, na.rm=TRUE), arr.ind=TRUE)
				if (dim(from.to)[1] > 1) {
					from.to <- from.to[sample.int(dim(from.to)[1], 1),]	
				} else {
					from.to <- from.to[1,]	
				}
				potential.changing.taxa <- subset(results, living & state==states[from.to[1]])$id
				changed.taxon <- potential.changing.taxa[sample.int(length(potential.changing.taxa), 1)]
				results[which(id==changed.taxon),]$state <- states[from.to[2]]
			}
			keep.running <- CheckKeepRunning(results, max.taxa, max.t, max.wall.time, start)
		}
	}
	return(results)
}

Multiply <- function(x, y) { #I know, this is silly. It's like the joke about Wickham's addr package
	return(x*y)	
}

GetBirthRate <- function(turnover, eps) {
	return(turnover / (1+eps))	
}

GetDeathRate <- function(turnover, eps) {
	return((turnover * eps) / (1 + eps))	
}

CheckKeepRunning <- function(results, max.taxa=Inf, max.t=Inf, max.wall.time=Inf, start=NULL) {
	keep.running <- TRUE
	if(dim(subset(results, living))[1]<1) {
		keep.running <- FALSE	
	}
	if(dim(subset(results, living))[1]>=max.taxa) {
		keep.running <- FALSE
	}
	if( max(subset(results, living)$height)>=max.t) {
		keep.running <- FALSE
	}
	if(keep.running & is.finite(max.wall.time)) { #could slow us down, so only check if needed
		if((as.numeric(Sys.time() - start)) > max.wall.time) {
			keep.running <- FALSE
		}
	}	
	return(keep.running)
}