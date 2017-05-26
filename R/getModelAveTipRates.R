

GetModelAveTipRates <- function(x){
	hisse.results <- x
	if(class(hisse.results)=="hisse.states") { #we have to make a list so we can run this generally
		if(is.null(hisse.results$aic)){
			#If a user forgot to include the aic, then we add a random value in for them
			hisse.results$aic = 42
		}
		tmp.list <- list()
		tmp.list[[1]] <- hisse.results
		hisse.results <- tmp.list	
	}

	rates.tips.turnover <- ConvertManyToRate(hisse.results, rate.param="turnover", "tip.mat")
	rates.tips.net.div <- ConvertManyToRate(hisse.results, rate.param="net.div", "tip.mat")
	rates.tips.speciation <- ConvertManyToRate(hisse.results, rate.param="speciation", "tip.mat")
	rates.tips.extinct.fraction <- ConvertManyToRate(hisse.results, rate.param="extinction.fraction", "tip.mat")
	rates.tips.extinction <- ConvertManyToRate(hisse.results, rate.param="extinction", "tip.mat")

	states.tips <- ConvertManyToBinaryState(hisse.results, "tip.mat")
	
	final.df <- data.frame(taxon=x[[1]]$phy$tip.label, state=states.tips, turnover=rates.tips.turnover, net.div=rates.tips.net.div, speciation=rates.tips.speciation, extinct.frac=rates.tips.extinct.fraction, extinction=rates.tips.extinction)
	return(final.df)
}



GetModelAveNodeRates <- function(x){
	hisse.results <- x
	if(class(hisse.results)=="hisse.states") { #we have to make a list so we can run this generally
		if(is.null(hisse.results$aic)){
			#If a user forgot to include the aic, then we add a random value in for them
			hisse.results$aic = 42
		}
		tmp.list <- list()
		tmp.list[[1]] <- hisse.results
		hisse.results <- tmp.list	
	}

	rates.tips.turnover <- ConvertManyToRate(hisse.results, rate.param="turnover", "node.mat")
	rates.tips.net.div <- ConvertManyToRate(hisse.results, rate.param="net.div", "node.mat")
	rates.tips.speciation <- ConvertManyToRate(hisse.results, rate.param="speciation", "node.mat")
	rates.tips.extinct.fraction <- ConvertManyToRate(hisse.results, rate.param="extinction.fraction", "node.mat")
	rates.tips.extinction <- ConvertManyToRate(hisse.results, rate.param="extinction", "node.mat")

	states.internal <- ConvertManyToBinaryState(hisse.results, "node.mat")
	
	final.df <- data.frame(id=x[[1]]$node.mat[,1], state=states.internal, turnover=rates.tips.turnover, net.div=rates.tips.net.div, speciation=rates.tips.speciation, extinct.frac=rates.tips.extinct.fraction, extinction=rates.tips.extinction)
	return(final.df)
}
