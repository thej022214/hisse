SimulateHisse <- function(turnover.rates, eps.rates, transition.rates, max.taxa=Inf, max.t=Inf, max.walltime=Inf, x0, nstart=1) {
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
	results <- data.table(anc=NA, id=1, state=factor(x0, state.levels), length=0, height=0, living=TRUE, descendants=FALSE)
	for(additional.row in sequence(nstart-1)) {
		results <- rbind(results, data.table(anc=NA, id=additional.row+1, state=factor(x0, state.levels), length=0, height=0, living=TRUE, descendants=FALSE))
	}
	setkey(results, id, living)
	birth.rates <- GetBirthRate(turover.rates, eps.rates)
	death.rates <- GetDeathRate(turnover.rates, eps.rates)
	diag(transition.rates) <- NA
	all.rates <- c(birth.rates, death.rates, c(t(transition.rates))) #c(t()) so that we have first row (from 0 to 1,2,3), then second row, etc, rather than by columns
	all.rates <- all.rates[!is.na(all.rates)] #drop the diagonals
	keep.running <- TRUE
	while(keep.running) {
		tip.state.counts <- table(subset(results, living)$state)
		
		keep.running <- CheckKeepRunning(results, max.taxa, max.t, max.wall.time, start)
	}
	return(results)
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