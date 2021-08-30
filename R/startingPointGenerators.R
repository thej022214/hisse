
######################################################################################################################################
######################################################################################################################################
### Generates starting parameter values for HiSSE, MuHiSSE, MiSSE without fossils
######################################################################################################################################
######################################################################################################################################

#Taken from the BiSSE code -- credit goes to Rich FitzJohn:
starting.point.tree <- function(phy, yule=FALSE) {
  if(!ape::is.binary(phy)) {
    phy <- ape::multi2di(phy)
  }
	p.yule <- c(yule(phy)$lambda, 0)
	if(yule){
		p.yule
	}else{
		suppressWarnings(c(birthdeath(phy)$para[2] / (1-birthdeath(phy)$para[1]), ((birthdeath(phy)$para[1] * birthdeath(phy)$para[2]) / (1-birthdeath(phy)$para[1]))))
	}
}


starting.point.generator <- function(phy, k, samp.freq.tree, q.div=5, yule=FALSE) {
    if(!ape::is.binary(phy)) {
        cat("Tree has polytomies:\n  Resolving randomly for initial parameter guesses only (hisse will use the tree with polytomies as given when optmizing).\n  Note that correctness of solutions with polytomies has not been established.\n")
    }
    pars.bd <- suppressWarnings(starting.point.tree(phy, yule))
    #Rescale parameters to account for sampling, if necessary, using Stadler 2013:
    pars.bd[1] = pars.bd[1] / samp.freq.tree
    pars.bd[2] = pars.bd[2] - ((pars.bd[1]*samp.freq.tree) * (1 - 1/samp.freq.tree))
    r <- if( pars.bd[1] > pars.bd[2] )
    (pars.bd[1] - pars.bd[2]) else pars.bd[1]
    p <- rep(c(pars.bd, r / q.div), c(k, k, k * (k-1)))
    names(p) <- NULL
    p
}


######################################################################################################################################
######################################################################################################################################
### Generates starting parameter values for GeoHiSSE
######################################################################################################################################
######################################################################################################################################

#Taken from the GeoSSE code modified to account for sampling -- credit goes to Emma Goldberg:
starting.point.geosse <- function(tree, eps=0.5, samp.freq.tree) {
    if (eps == 0) {
        s <- (log(Ntip(tree)) - log(2)) / max(branching.times(tree))
        s <- s/samp.freq.tree
        x <- 0
        d <- s/10
    } else {
        n <- Ntip(tree)
        r <- ( log( (n/2) * (1 - eps*eps) + 2*eps + (1 - eps)/2 *
        sqrt( n * (n*eps*eps - 8*eps + 2*n*eps + n))) - log(2)
        ) / max(branching.times(tree))
        s <- r / (1 - eps)
        x <- s * eps
        s <- s/samp.freq.tree
        x <- x - ((s * samp.freq.tree) * (1 - 1/samp.freq.tree))
        d <- x
    }
    p <- c(s, s, s, x, x, d, d)
    names(p) <- c("sA",  "sB",  "sAB", "xA" , "xB"  ,"dA"  ,"dB")
    p
}


######################################################################################################################################
######################################################################################################################################
### Generates starting parameter values for HiSSE, MuHiSSE, and MiSSE with fossils
######################################################################################################################################
######################################################################################################################################

#c1 and c2 are taken from Stadler 2010, pg. 398
c1 <- function(lambda, mu, psi){
    return(abs(sqrt(((lambda - mu - psi)^2)+(4*lambda*psi))))
}

c2 <- function(lambda, mu, psi, rho){
    return(-(lambda-mu-(2*lambda*rho)-psi)/c1(lambda, mu, psi))
}

#This is Master equation from Stalder 2010, pg. 398.
p_0 <- function(x, lambda, mu, psi, rho){
    c1_val <- c1(lambda=lambda, mu=mu, psi=psi)
    c2_val <- c2(lambda=lambda, mu=mu, psi=psi, rho=rho)
    return(log((lambda + mu + psi + c1_val * ((exp(-c1_val*x) * (1-c2_val)-(1+c2_val))/(exp(-c1_val*x) * (1-c2_val)+(1+c2_val))))/(2*lambda)))
}


#This is Master equation from Stadler 2010, pg. 398.
p_one <- function(x,lambda, mu, psi, rho){
    c1_val <- c1(lambda=lambda, mu=mu, psi=psi)
    c2_val <- c2(lambda=lambda, mu=mu, psi=psi, rho=rho)
    return(log((4*rho)/(2*(1-c2_val^2) + exp(-c1_val*x)*(1-c2_val)^2 + exp(c1_val*x) * (1+c2_val)^2)))
}


starting.point.tree.fossils <- function(x, rho, n, m, k, x_times, y_times, interval.sum){
    x <- exp(x)
    lambda <- x[1] / (1 + x[2])
    mu <- (x[2] * x[1]) / (1 + x[2])
    psi <- x[3]
    
    #Equation 5 from Stadler 2010, pg. 400 -- need to readjust the likelihood equation for logspace.
    #loglik <- log( (((lambda^(n+m-2)) * (psi^m))/(1-p_0(max(x_times),lambda,mu,psi,rho))^2) * p_one(max(x_times), lambda, mu, psi, rho) * prod(p_one(x_times, lambda,mu,psi,rho)) * prod(p_0(y_times,lambda,mu,psi,rho)/p_one(y_times,lambda,mu,psi,rho)) )
    loglik <- (((n+m-2) * log(lambda)) + ((k+m) * log(psi))) - log(1-exp(p_0(max(x_times),lambda,mu,psi=0,rho)))*2 + p_one(max(x_times), lambda, mu, psi, rho) + sum(p_one(x_times, lambda,mu,psi,rho)) + (sum(p_0(y_times,lambda,mu,psi,rho)) - sum(p_one(y_times,lambda,mu,psi,rho)))
    
    if(!is.null(interval.sum)){
        loglik <- loglik + (psi * interval.sum)
    }
    
    return(-loglik)
}


starting.point.generator.fossils <- function(n.tax, k, samp.freq.tree, q.div=5, fossil.taxa, fossil.ages, no.k.samples, split.times, interval.sum, get.likelihood=FALSE) {
    opts <- list("algorithm" = "NLOPT_LN_NELDERMEAD", "maxeval" = 100000, "ftol_rel" = .Machine$double.eps^.5)
    out <- nloptr(x0=log(c(0.1, 0.2, 0.01)), eval_f=starting.point.tree.fossils, ub=log(c(10, 0.99, 1)), lb=c(-21,-21, -21), opts=opts, rho=samp.freq.tree, n=n.tax, m=length(fossil.taxa), k=no.k.samples, x_times=split.times, y_times=fossil.ages, interval.sum=interval.sum)
    if(get.likelihood == TRUE){
        starting.rates <- exp(out$solution)
        lambda <- starting.rates[1] / (1 + starting.rates[2])
        mu <- (starting.rates[1] * starting.rates[2]) / (1 + starting.rates[2])
        psi <- starting.rates[3]
        return(c(lambda, mu, psi, out$objective))
    }else{
        starting.rates <- exp(out$solution)
        lambda <- starting.rates[1] / (1 + starting.rates[2])
        mu <- (starting.rates[1] * starting.rates[2]) / (1 + starting.rates[2])
        psi <- starting.rates[3]
        r <- if( lambda > mu )
        (lambda - mu) else lambda
        p <- c(rep(c(c(lambda, mu), r / q.div), c(k, k, k * (k-1))), psi)
        names(p) <- NULL
        return(p)
    }
}





