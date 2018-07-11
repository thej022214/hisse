test_that("BiSSE_HiSSE_test",{
	skip_on_cran()
    
    library(diversitree)
	pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
	set.seed(4)
	phy <- tree.bisse(pars, max.t=30, x0=0)
	lik <- make.bisse(phy, phy$tip.state, sampling.f=c(.4,.6))
	diversitree.full <- lik(pars)
	
	hidden.states=FALSE
	states <- data.frame(phy$tip.state, phy$tip.state, row.names=names(phy$tip.state))
	states <- states[phy$tip.label,]
	pars.bisse <- c(0.1+0.03, 0.2+0.03, 0, 0, 0.03/0.1, 0.03/0.2, 0, 0, 0.01, 0, 0, 0.01, 0, 0, 0, 0, 0, 0, 0, 0)
	model.vec = c(pars.bisse, rep(1,36))
	phy$node.label = NULL
	cache = ParametersToPass(phy, states[,1], model.vec, f=c(.4,.6), timeslice=NULL, hidden.states=hidden.states)
	cache$turnover.beta.factor0 = 1 / dbeta(0.1, 1, 1)
	cache$eps.beta.factor0 = 1 / dbeta(0.1, 1, 1)
	cache$turnover.beta.factor1 = 1 / dbeta(0.1, 1, 1)
	cache$eps.beta.factor1 = 1 / dbeta(0.1, 1, 1)
	cache$turnover.beta.factorA = 1 / dbeta(0.1, 1, 1)
	cache$eps.beta.factorA = 1 / dbeta(0.1, 1, 1)
	cache$turnover.beta.factorB = 1 / dbeta(0.1, 1, 1)
	cache$eps.beta.factorB = 1 / dbeta(0.1, 1, 1)
    hisse.full <- DownPass(phy, cache, hidden.states=TRUE, root.type="madfitz", condition.on.survival=TRUE)
	comparison <- identical(round(hisse.full,4), round(diversitree.full,4))
	expect_true(comparison)
})


test_that("MuSSE_HiSSE_test1", {
    skip_on_cran()

	library(diversitree)
	pars <- c(.1,  .15,  .2,  # lambda 1, 2, 3
			  .03, .045, .06, # mu 1, 2, 3
			  .05, 0,         # q12, q13
			  .05, .05,       # q21, q23
			  0,   .05)       # q31, q32
	set.seed(2)
	phy <- tree.musse(pars, 30, x0=1)
	states <- phy$tip.state
	lik <- make.musse(phy, states, 3)
	lik.base <- constrain(lik, lambda2 ~ lambda1, lambda3 ~ lambda1,
						  mu2 ~ mu1, mu3 ~ mu1,
						  q13 ~ 0, q21 ~ q12, q23 ~ q12, q31 ~ 0, q32 ~ q12)
	diversitree.constrained = lik.base(c(.1, .03, .05))
	diversitree.full = lik(pars)
	
	hidden.states="TEST"
	states <- data.frame(phy$tip.state, phy$tip.state, row.names=names(phy$tip.state))
	states <- states[phy$tip.label,]
	states[states[,1]==3,] = 4
	pars.hisse <- c(0.1+0.03, 0.1+0.03, 0, 0.1+0.03, 0.03/0.1, 0.03/0.1, 0, 0.03/0.1, 0.05, 0, 0, 0.05, 0, 0.05,0, 0, 0, 0, 0.05, 0)
	model.vec = c(pars.hisse, rep(1,36))
	phy$node.label = NULL
	cache <- ParametersToPass(phy, states[,1], model.vec, f=c(1,1), timeslice=NULL, hidden.states=hidden.states)
	cache$turnover.beta.factor0 = 1 / dbeta(0.1, 1, 1)
	cache$eps.beta.factor0 = 1 / dbeta(0.1, 1, 1)
	cache$turnover.beta.factor1 = 1 / dbeta(0.1, 1, 1)
	cache$eps.beta.factor1 = 1 / dbeta(0.1, 1, 1)
	cache$turnover.beta.factorA = 1 / dbeta(0.1, 1, 1)
	cache$eps.beta.factorA = 1 / dbeta(0.1, 1, 1)
	cache$turnover.beta.factorB = 1 / dbeta(0.1, 1, 1)
	cache$eps.beta.factorB = 1 / dbeta(0.1, 1, 1)
	hisse.constrained <- DownPass(phy, cache, hidden.states=TRUE, root.type="madfitz", condition.on.survival=TRUE)
	comparison <- identical(round(hisse.constrained,4), round(diversitree.constrained,4))
	expect_true(comparison)
})


test_that("MuSSE_HiSSE_test2", {
    skip_on_cran()

	library(diversitree)
	pars <- c(.1,  .15,  .2,  # lambda 1, 2, 3
			  .03, .045, .06, # mu 1, 2, 3
			  .05, 0,         # q12, q13
			  .05, .05,       # q21, q23
			  0,   .05)       # q31, q32
	set.seed(2)
	phy <- tree.musse(pars, 30, x0=1)
	states <- phy$tip.state
	lik <- make.musse(phy, states, 3)
	diversitree.full = lik(pars)
	
	hidden.states="TEST"
	states <- data.frame(phy$tip.state, phy$tip.state, row.names=names(phy$tip.state))
	states <- states[phy$tip.label,]
	states[states[,1]==3,] = 4
	pars.hisse <- c(0.1+0.03, 0.15+0.045, 0, 0.2+0.06, 0.03/0.1, 0.045/0.15, 0, 0.06/0.2, 0.05, 0, 0, 0.05, 0, 0.05,0, 0, 0, 0, 0.05, 0)
	model.vec = c(pars.hisse, rep(1,36))
	phy$node.label = NULL
	cache = ParametersToPass(phy, states[,1], model.vec, f=c(1,1), timeslice=NULL, hidden.states=hidden.states)
	cache$turnover.beta.factor0 = 1 / dbeta(0.1, 1, 1)
	cache$eps.beta.factor0 = 1 / dbeta(0.1, 1, 1)
	cache$turnover.beta.factor1 = 1 / dbeta(0.1, 1, 1)
	cache$eps.beta.factor1 = 1 / dbeta(0.1, 1, 1)
	cache$turnover.beta.factorA = 1 / dbeta(0.1, 1, 1)
	cache$eps.beta.factorA = 1 / dbeta(0.1, 1, 1)
	cache$turnover.beta.factorB = 1 / dbeta(0.1, 1, 1)
	cache$eps.beta.factorB = 1 / dbeta(0.1, 1, 1)
	hisse.full <- DownPass(phy, cache, hidden.states=TRUE, root.type="madfitz", condition.on.survival=TRUE)
	comparison <- identical(round(hisse.full,4), round(diversitree.full,4))
	expect_true(comparison)
})


test_that("HiSSE_Null_Four_test", {
    skip_on_cran()

	library(diversitree)
	pars <- c(0.1, 0.1, 0.03, 0.03, 0.01, 0.01)
	set.seed(4)
	phy <- tree.bisse(pars, max.t=30, x0=0)
	lik <- make.bisse(phy, phy$tip.state)
	diversitree.full <- lik(pars)
	
	states <- data.frame(phy$tip.state, phy$tip.state, row.names=names(phy$tip.state))
	states <- states[phy$tip.label,]
	pars.hisse.null <- c(rep(0.1+0.03,8), rep(0.03/.1, 8), rep(0.01, 32))
	model.vec = pars.hisse.null
	phy$node.label = NULL
	cache <- ParametersToPassNull(phy, states[,1], model.vec, f=c(1,1))
	hisse.full <- DownPassNull(phy, cache, root.type="madfitz", condition.on.survival=TRUE)
	comparison <- identical(round(hisse.full,4), round(diversitree.full,4))
	expect_true(comparison)
})


test_that("GeoHiSSE_test1", {
    skip_on_cran()
    
    library(diversitree)
    pars <- c(1.5, 0.5, 1.0, 0.7, 0.7, 2.5, 0.5)
    names(pars) <- diversitree:::default.argnames.geosse()
    set.seed(5)
    phy <- tree.geosse(pars, max.t=4, x0=0)
    lik <- make.geosse(phy, phy$tip.state)
    diversitree.full <- lik(pars)
    
    states <- data.frame(phy$tip.state, phy$tip.state, row.names=names(phy$tip.state))
    states <- states[phy$tip.label,]
    names(pars) <- NULL
    model.vec <- numeric(115)
    model.vec[1:11] <- c(pars[1:3], pars[4:5], 0, pars[6], 0, pars[7], pars[5:4])
    phy$node.label <- NULL
    cache <- ParametersToPassGeoHiSSE(phy, states[,1], f=c(1,1,1), model.vec, hidden.states=FALSE)
    geohisse.full <- DownPassGeoHisse(phy=phy, cache=cache, hidden.states=FALSE, bad.likelihood=-1000000, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL)
    comparison <- identical(round(geohisse.full,4), round(diversitree.full,4))
    expect_true(comparison)
})


test_that("GeoHiSSE_test2", {
    skip_on_cran()
    
    library(diversitree)
    pars <- c(1.5, 0.5, 1.0, 0.7, 0.7, 2.5, 0.5)
    names(pars) <- diversitree:::default.argnames.geosse()
    set.seed(5)
    phy <- tree.geosse(pars, max.t=4, x0=0)
    lik <- make.geosse(phy, phy$tip.state)
    diversitree.full <- lik(pars)
    
    states <- data.frame(phy$tip.state, phy$tip.state, row.names=names(phy$tip.state))
    states <- states[phy$tip.label,]
    names(pars) <- NULL
    model.vec <- numeric(115)
    model.vec[1:11] <- c(pars[1:3], pars[4:5], 0, pars[6], 0, pars[7], pars[5:4])
    phy$node.label <- NULL
    cache <- ParametersToPassGeoHiSSE(phy, states[,1], f=c(1,1,1), model.vec, hidden.states=TRUE)
    geohisse.full <- DownPassGeoHisse(phy=phy, cache=cache, hidden.states=TRUE, bad.likelihood=-1000000, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL)
    comparison <- identical(round(geohisse.full,4), round(diversitree.full,4))
    expect_true(comparison)
})


test_that("GeoHiSSE_test3", {
    skip_on_cran()
    
    library(diversitree)
    pars <- c(1.5, 0.5, 1.0, 0.7, 0.4, 2.5, 0.5)
    names(pars) <- diversitree:::default.argnames.geosse()
    set.seed(9)
    phy <- tree.geosse(pars, max.t=4, x0=0)
    lik <- make.geosse(phy, phy$tip.state)
    diversitree.full <- lik(pars)
    
    states <- data.frame(phy$tip.state, phy$tip.state, row.names=names(phy$tip.state))
    states <- states[phy$tip.label,]
    names(pars) <- NULL
    model.vec <- numeric(115)
    model.vec[1:11] <- c(pars[1:3], pars[4:5], 0, pars[6], 0, pars[7], pars[5:4])
    phy$node.label <- NULL
    cache <- ParametersToPassGeoHiSSE(phy, states[,1], f=c(1,1,1), model.vec, hidden.states=FALSE)
    geohisse.full <- DownPassGeoHisse(phy=phy, cache=cache, hidden.states=TRUE, bad.likelihood=-1000000, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL)
    comparison <- identical(round(geohisse.full,4), round(diversitree.full,4))
    expect_true(comparison)
})


test_that("GeoHiSSE_test4", {
    skip_on_cran()
    
    library(diversitree)
    pars <- c(1.5, 0.5, 1.0, 0.7, 0.4, 2.5, 0.5)
    names(pars) <- diversitree:::default.argnames.geosse()
    set.seed(9)
    phy <- tree.geosse(pars, max.t=4, x0=0)
    lik <- make.geosse(phy, phy$tip.state)
    diversitree.full <- lik(pars)
    
    states <- data.frame(phy$tip.state, phy$tip.state, row.names=names(phy$tip.state))
    states <- states[phy$tip.label,]
    names(pars) <- NULL
    model.vec <- numeric(115)
    model.vec[1:11] <- c(pars[1:3], pars[4:5], 0, pars[6], 0, pars[7], pars[5:4])
    phy$node.label <- NULL
    cache <- ParametersToPassGeoHiSSE(phy, states[,1], f=c(1,1,1), model.vec, hidden.states=TRUE)
    geohisse.full <- DownPassGeoHisse(phy=phy, cache=cache, hidden.states=TRUE, bad.likelihood=-1000000, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL)
    comparison <- identical(round(geohisse.full,4), round(diversitree.full,4))
    expect_true(comparison)
})

test_that("GeoHiSSE_test5", {
    skip_on_cran()
    
    library(diversitree)
    
    ## Define some parameters:
    ##         s01,  s0,  s1,  x0,  x1,  d0,  d1
    pars1 <- c(1.5, 0.5, 1.0, 0.7, 0.7, 2.5, 0.5)
    pars2 <- c(1.1, 1.6, 0.4, 0.6, 0.5, 0.6, 2.7)
    rates.q <- 0.1 ## Transition between hidden rates.
    
    ## Check the likelihood for the hidden null model.
    sim.pars <- SimulateGeoHiSSE(hidden.areas = 1, return.GeoHiSSE_pars = TRUE)
    sim.pars$model.pars[,1] <- pars1
    sim.pars$model.pars[,2] <- pars2
    sim.pars$q.01[2,1] <- rates.q
    sim.pars$q.01[1,2] <- rates.q
    sim.pars$q.0[2,1] <- rates.q
    sim.pars$q.0[1,2] <- rates.q
    sim.pars$q.1[2,1] <- rates.q
    sim.pars$q.1[1,2] <- rates.q

    sim.get.par <- SimulateGeoHiSSE(pars=sim.pars, hidden.areas = 1, max.taxa = 500)
    classe.pars <- sim.get.par$classe.pars
    sim.data <- NULL
    while( is.null(sim.data) ) sim.data <- tree.classe(pars=classe.pars, max.taxa = 500)
   
    ## Get the likelihood for the ClaSSE model:
    ## Here we assume root value is EQUAL
    sim.data$node.label <- NULL
    lik.classe <- make.classe(tree=sim.data, states=sim.data$tip.state, k=6)
    classe.full <- lik.classe(pars=classe.pars, root=ROOT.FLAT, condition.surv = TRUE, root.p = NULL)
    
    ## Now the lik for the GeoHiSSE model.
    ## Our function starts from 0.
    states <- sim.data$tip.state - 1
    states.mat <- data.frame(states, states, row.names=names(states))
    order.mat <- match(rownames(states.mat), sim.data$tip.label)
    states.mat <- states.mat[order.mat,]
    model.vec <- numeric(115)
    order.pars1 <- c(pars1[c(2,3,1,4:5)],0,pars1[6],0,pars1[c(7,5:4)], rates.q,0,0,0, rates.q,0,0,0
                   , rates.q,0,0,0)
    order.pars2 <- c(pars2[c(2,3,1,4:5)],0,pars2[6],0,pars2[c(7,5:4)], rates.q,0,0,0, rates.q,0,0,0
                   , rates.q,0,0,0)
    order.pars <- c(order.pars1, order.pars2)
    model.vec[1:46] <- order.pars
    
    cache <- ParametersToPassGeoHiSSE(sim.data, states.mat[,1], f=c(1,1,1), model.vec, hidden.states="TEST")
    geohisse.full <- DownPassGeoHisse(phy=sim.data, cache=cache, hidden.states=TRUE, bad.likelihood=-1000000
                                    , condition.on.survival=TRUE, root.type="equal", root.p=NULL)
    comparison <- identical(round(geohisse.full,4), round(classe.full,4))
    expect_true(comparison)
})

test_that("MuSSE_test1", {
    skip_on_cran()
    
    library(diversitree)
    pars <- c(.1,  .15,  .2,  # lambda 1, 2, 3
    .03, .045, .06, # mu 1, 2, 3
    .05, 0,         # q12, q13
    .05, .05,       # q21, q23
    0,   .05)       # q31, q32
    set.seed(2)
    phy <- tree.musse(pars, 30, x0=1)
    states <- phy$tip.state
    lik <- make.musse(phy, states, 3)
    lik.base <- constrain(lik, lambda2 ~ lambda1, lambda3 ~ lambda1,
    mu2 ~ mu1, mu3 ~ mu1,
    q13 ~ 0, q21 ~ q12, q23 ~ q12, q31 ~ 0.03, q32 ~ q12)
    diversitree.constrained = lik.base(c(.1, .03, .05))
    diversitree.full = lik(pars)
    
    states <- data.frame(phy$tip.state, phy$tip.state, row.names=names(phy$tip.state))
    states <- states[phy$tip.label,]
    states[states[,1]==3,] = 4
    pars.hisse <- c(0.1, 0.1, 0.1, 0.03, 0.03, 0.03, 0.05, 0, 0.05, 0.05, 0, 0.05)
    model.vec = rep(0,120)
    model.vec[1:12] = pars.hisse
    phy$node.label = NULL
    cache <- hisse:::ParametersToPassMuSSE(phy, states[,1], model.vec, f=c(1,1,1), hidden.states="TEST1")
    hisse.constrained <- hisse:::DownPassMusse(phy, cache, hidden.states=FALSE, root.type="madfitz", condition.on.survival=TRUE)
    comparison <- identical(round(hisse.constrained,4), round(diversitree.constrained,4))
    expect_true(comparison)
})


test_that("MuSSE_test2", {
    skip_on_cran()

    library(diversitree)
    pars <- c(.1,  .15,  .2,  # lambda 1, 2, 3
    .03, .045, .06, # mu 1, 2, 3
    .05, 0,         # q12, q13
    .05, .05,       # q21, q23
    0,   .05)       # q31, q32
    set.seed(2)
    phy <- tree.musse(pars, 30, x0=1)
    states <- phy$tip.state
    lik <- make.musse(phy, states, 3)
    lik.base <- constrain(lik,
    q13 ~ 0, q21 ~ q12, q23 ~ q12, q31 ~ 0.04, q32 ~ q12)
    diversitree.constrained = lik.base(c(.1, .2, .3, .03,.04,.05, .05))
    #diversitree.full = lik(pars)
    
    states <- data.frame(phy$tip.state, phy$tip.state, row.names=names(phy$tip.state))
    states <- states[phy$tip.label,]
    states[states[,1]==3,] = 4
    pars.hisse <- c(.1, .2, .3, .03, .04, .05, 0.05, 0, 0.05, 0.05, 0, 0.05)
    model.vec = rep(0,120)
    model.vec[25:36] = pars.hisse
    phy$node.label = NULL
    cache <- ParametersToPassMuSSE(phy, states[,1], model.vec, f=c(1,1,1), hidden.states="TEST2")
    hisse.constrained <- DownPassMusse(phy, cache, hidden.states=TRUE, root.type="madfitz", condition.on.survival=TRUE)
    comparison <- identical(round(hisse.constrained,4), round(diversitree.constrained,4))
    expect_true(comparison)
})

