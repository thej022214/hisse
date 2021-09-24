
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
p_0 <- function(x, lambda, mu, psi, rho, log=TRUE){
    c1_val <- c1(lambda=lambda, mu=mu, psi=psi)
    c2_val <- c2(lambda=lambda, mu=mu, psi=psi, rho=rho)
    if(log == TRUE){
        return(log((lambda + mu + psi + c1_val * ((exp(-c1_val*x) * (1-c2_val)-(1+c2_val))/(exp(-c1_val*x) * (1-c2_val)+(1+c2_val))))/(2*lambda)))
    }else{
        return((lambda + mu + psi + c1_val * ((exp(-c1_val*x) * (1-c2_val)-(1+c2_val))/(exp(-c1_val*x) * (1-c2_val)+(1+c2_val))))/(2*lambda))
    }
}


#This is Master equation from Stadler 2010, pg. 398.
p_one <- function(x,lambda, mu, psi, rho){
    c1_val <- c1(lambda=lambda, mu=mu, psi=psi)
    c2_val <- c2(lambda=lambda, mu=mu, psi=psi, rho=rho)
    return(log((4*rho)/(2*(1-c2_val^2) + exp(-c1_val*x)*(1-c2_val)^2 + exp(c1_val*x) * (1+c2_val)^2)))
}


starting.point.tree.fossils <- function(x, rho, n, m, k, x_times, y_times){
    x <- exp(x)
    lambda <- x[1] / (1 + x[2])
    mu <- (x[2] * x[1]) / (1 + x[2])
    psi <- x[3]
    
    #Equation 5 from Stadler 2010, pg. 400 -- need to readjust the likelihood equation for logspace.
    #loglik <- log( (((lambda^(n+m-2)) * (psi^m))/(1-p_0(max(x_times),lambda,mu,psi,rho))^2) * p_one(max(x_times), lambda, mu, psi, rho) * prod(p_one(x_times, lambda,mu,psi,rho)) * prod(p_0(y_times,lambda,mu,psi,rho)/p_one(y_times,lambda,mu,psi,rho)) )
    loglik <- (((n+m-2) * log(lambda)) + ((k+m) * log(psi))) - log(1-p_0(max(x_times),lambda,mu,psi=0,rho, log=FALSE))*2 + p_one(max(x_times), lambda, mu, psi, rho) + sum(p_one(x_times, lambda,mu,psi,rho)) + (sum(p_0(y_times,lambda,mu,psi,rho)) - sum(p_one(y_times,lambda,mu,psi,rho)))
    
    return(-loglik)
}


starting.point.generator.fossils <- function(n.tax, k, samp.freq.tree, q.div=5, fossil.taxa, fossil.ages, no.k.samples, split.times, get.likelihood=FALSE) {
    opts <- list("algorithm" = "NLOPT_LN_NELDERMEAD", "maxeval" = 100000, "ftol_rel" = .Machine$double.eps^.5)
    out <- nloptr(x0=log(c(0.1, 0.2, 0.01)), eval_f=starting.point.tree.fossils, ub=log(c(10, 0.99, 1)), lb=c(-21,-21, -21), opts=opts, rho=samp.freq.tree, n=n.tax, m=length(fossil.taxa), k=no.k.samples, x_times=split.times, y_times=fossil.ages)
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


######################################################################################################################################
######################################################################################################################################
### Generates starting parameter values for HiSSE, MuHiSSE, and MiSSE with stratigraphic intervals
######################################################################################################################################
######################################################################################################################################

#This is Eq. 4 from Stadler et al 2018, pg. 46
q_t <- function(t, lambda, mu, psi, rho){
    c1_val <- c1(lambda=lambda, mu=mu, psi=psi)
    c2_val <- c2(lambda=lambda, mu=mu, psi=psi, rho=rho)
    qt_val <- (4*exp(-c1_val*t)) / ((exp(-c1_val*t) * (1 - c2_val) + (1 + c2_val))^2)
    return(qt_val)
}

#This is Eq. from Stadler et al 2018, pg. 51
q_sym_tilde <- function(t, lambda, mu, psi, rho){
    q_sym_val <- exp(-(lambda + mu + psi)*t)
    return(q_sym_val)
}

#This is Eq. from Stadler et al 2018, pg. 51
q_sym_ratio <- function(si, ei, lambda, mu, psi, rho){
    q_sym_val1 <- exp(-(lambda + mu + psi)*si)
    q_sym_val2 <- exp(-(lambda + mu + psi)*ei)
    q_sym_hat_val <- q_sym_val1 / q_sym_val2
    return(q_sym_hat_val)
}

#This is Eq. from Stadler et al 2018, pg. 51
q_ratio <- function(si, ei, lambda, mu, psi, rho){
    q_val1 <- q_t(t=si, lambda=lambda, mu=mu, psi=psi, rho=rho)
    q_val2 <- q_t(t=ei, lambda=lambda, mu=mu, psi=psi, rho=rho)
    q_sym_hat_val <- q_val1 / q_val2
    return(q_sym_hat_val)
}

#This is Eq. from Stadler et al 2018, pg. 51
GetB_i <- function(x, lambda, mu, psi, rho){
    res <- numeric(dim(x)[1])
    for(index in 1:dim(x)[1]){
        if(x$branch.type[index] == 3){
            res[index] <- q_sym_ratio(si=x$RootwardAge[index], ei=x$TipwardAge[index], lambda=lambda, mu=mu, psi=psi, rho=rho)
        }else{
            res[index] <- q_ratio(si=x$RootwardAge[index], ei=x$TipwardAge[index], lambda=lambda, mu=mu, psi=psi, rho=rho)
        }
    }
    return(res)
}

#This is right hand most product from corollary 13 from Stadler et al 2018, pg. 51
GetUnobsSpecEventProbs <- function(x, lambda, mu, psi, rho){
    res <- numeric(dim(x)[1])
    for(index  in 1:dim(x)[1]){
        first.ratio <- q_t(t=x[index,4], lambda=lambda, mu=mu, psi=psi, rho=rho) / q_sym_tilde(t=x[index,4],lambda=lambda, mu=mu, psi=psi, rho=rho)
        second.ratio <- q_sym_tilde(t=x[index,3], lambda=lambda, mu=mu, psi=psi, rho=rho)/ q_t(t=x[index,3],lambda=lambda, mu=mu, psi=psi, rho=rho)
        res[index] <- 1 - (first.ratio * second.ratio)
    }
    return(res)
}


starting.point.tree.intervals <- function(x, n.tax, rho, seg_map, x_times, y_times, strat.cache){
    x <- exp(x)
    lambda <- x[1] / (1 + x[2])
    mu <- (x[2] * x[1]) / (1 + x[2])
    psi <- x[3]
    k <- strat.cache$k
    intervening.intervals <- strat.cache$intervening.intervals
    l_s <- strat.cache$l_s
    
    #STEP 1: Get branch-wise probs:
    branch.segs <- GetB_i(x=seg_map, lambda=lambda, mu=mu, psi=psi, rho=rho)
    #STEP 2: Get extinct starting probs:
    extinct.starts <- p_0(x=y_times, lambda=lambda, mu=mu, psi=psi, rho=rho, log=FALSE)
    #STEP 3: Get double interval segs probs:
    if(!is.null(intervening.intervals)){
        unobs.spec <- GetUnobsSpecEventProbs(x=intervening.intervals, lambda=lambda, mu=mu, psi=psi, rho=rho)
    }else{
        unobs.spec <- 1
    }
    
    #STEP 4: Get likelihood:
    loglik <- (l_s * psi) + sum(log(unobs.spec)) + sum(log(branch.segs)) + sum(log(extinct.starts*psi)) + (n.tax * log(rho)) + (length(x_times) * log(lambda)) + k*log(psi) - (log(1-p_0(x=max(x_times),lambda=lambda,mu=mu,psi=0,rho=rho,log=FALSE))*2 + log(lambda))
    
    return(-loglik)
}


starting.point.generator.intervals <- function(k, samp.freq.tree, q.div=5, n.tax, seg_map, split.times, fossil.ages, strat.cache, get.likelihood=FALSE) {
    opts <- list("algorithm" = "NLOPT_LN_NELDERMEAD", "maxeval" = 100000, "ftol_rel" = .Machine$double.eps^.5)
    out <- nloptr(x0=log(c(0.1, 0.2, 0.01)), eval_f=starting.point.tree.intervals, ub=log(c(10, 0.99, 2)), lb=c(-21,-21, -21), opts=opts, n.tax=n.tax, rho=samp.freq.tree, seg_map=seg_map, x_times=split.times, y_times=fossil.ages, strat.cache=strat.cache)
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





