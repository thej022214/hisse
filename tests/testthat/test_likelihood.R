test_that("BiSSE_HiSSE_test",{
	skip_on_cran()

    library(diversitree)
	pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
	set.seed(4)
        phy <- NULL
        while( is.null( phy ) ){
            phy <- tree.bisse(pars, max.t=30, x0=0)
        }
	lik <- make.bisse(phy, phy$tip.state, sampling.f=c(.4,.6))
	diversitree.full <- lik(pars)

	hidden.states=FALSE
	states <- data.frame(phy$tip.state, phy$tip.state, row.names=names(phy$tip.state))
	states <- states[phy$tip.label,]
	pars.bisse <- c(0.1+0.03, 0.2+0.03, 0, 0, 0.03/0.1, 0.03/0.2, 0, 0, 0.01, 0, 0, 0.01, 0, 0, 0, 0, 0, 0, 0, 0)
	model.vec = c(pars.bisse, rep(1,36))
	phy$node.label = NULL
	cache = hisse:::ParametersToPass(phy, states[,1], model.vec, f=c(.4,.6), timeslice=NULL, hidden.states=hidden.states)
	cache$turnover.beta.factor0 = 1 / dbeta(0.1, 1, 1)
	cache$eps.beta.factor0 = 1 / dbeta(0.1, 1, 1)
	cache$turnover.beta.factor1 = 1 / dbeta(0.1, 1, 1)
	cache$eps.beta.factor1 = 1 / dbeta(0.1, 1, 1)
	cache$turnover.beta.factorA = 1 / dbeta(0.1, 1, 1)
	cache$eps.beta.factorA = 1 / dbeta(0.1, 1, 1)
	cache$turnover.beta.factorB = 1 / dbeta(0.1, 1, 1)
	cache$eps.beta.factorB = 1 / dbeta(0.1, 1, 1)
    hisse.full <- hisse:::DownPass(phy, cache, hidden.states=TRUE, root.type="madfitz", condition.on.survival=TRUE, root.p=NULL)
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
    phy <- NULL
    while( is.null( phy ) ){
        phy <- tree.musse(pars, 30, x0=1)
    }
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
    cache <- hisse:::ParametersToPass(phy, states[,1], model.vec, f=c(1,1), timeslice=NULL, hidden.states=hidden.states)
	cache$turnover.beta.factor0 = 1 / dbeta(0.1, 1, 1)
	cache$eps.beta.factor0 = 1 / dbeta(0.1, 1, 1)
	cache$turnover.beta.factor1 = 1 / dbeta(0.1, 1, 1)
	cache$eps.beta.factor1 = 1 / dbeta(0.1, 1, 1)
	cache$turnover.beta.factorA = 1 / dbeta(0.1, 1, 1)
	cache$eps.beta.factorA = 1 / dbeta(0.1, 1, 1)
	cache$turnover.beta.factorB = 1 / dbeta(0.1, 1, 1)
	cache$eps.beta.factorB = 1 / dbeta(0.1, 1, 1)
    hisse.constrained <- hisse:::DownPass(phy, cache, hidden.states=TRUE, root.type="madfitz", condition.on.survival=TRUE, root.p=NULL)
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
    phy <- NULL
    while( is.null( phy ) ){
        phy <- tree.musse(pars, 30, x0=1)
    }
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
    cache = hisse:::ParametersToPass(phy, states[,1], model.vec, f=c(1,1), timeslice=NULL, hidden.states=hidden.states)
	cache$turnover.beta.factor0 = 1 / dbeta(0.1, 1, 1)
	cache$eps.beta.factor0 = 1 / dbeta(0.1, 1, 1)
	cache$turnover.beta.factor1 = 1 / dbeta(0.1, 1, 1)
	cache$eps.beta.factor1 = 1 / dbeta(0.1, 1, 1)
	cache$turnover.beta.factorA = 1 / dbeta(0.1, 1, 1)
	cache$eps.beta.factorA = 1 / dbeta(0.1, 1, 1)
	cache$turnover.beta.factorB = 1 / dbeta(0.1, 1, 1)
	cache$eps.beta.factorB = 1 / dbeta(0.1, 1, 1)
    hisse.full <- hisse:::DownPass(phy, cache, hidden.states=TRUE, root.type="madfitz", condition.on.survival=TRUE, root.p=NULL)
	comparison <- identical(round(hisse.full,4), round(diversitree.full,4))
	expect_true(comparison)
})


test_that("HiSSE_Null_Four_test", {
    skip_on_cran()

	library(diversitree)
	pars <- c(0.1, 0.1, 0.03, 0.03, 0.01, 0.01)
    set.seed(4)
    phy <- NULL
    while( is.null( phy ) ){
        phy <- tree.bisse(pars, max.t=30, x0=0)
    }
	lik <- make.bisse(phy, phy$tip.state)
	diversitree.full <- lik(pars)

	states <- data.frame(phy$tip.state, phy$tip.state, row.names=names(phy$tip.state))
	states <- states[phy$tip.label,]
	pars.hisse.null <- c(rep(0.1+0.03,8), rep(0.03/.1, 8), rep(0.01, 32))
	model.vec = pars.hisse.null
	phy$node.label = NULL
    cache <- hisse:::ParametersToPassNull(phy, states[,1], model.vec, f=c(1,1))
    hisse.full <- hisse:::DownPassNull(phy, cache, root.type="madfitz", condition.on.survival=TRUE, root.p=NULL)
	comparison <- identical(round(hisse.full,4), round(diversitree.full,4))
	expect_true(comparison)
})


test_that("GeoHiSSE_test1", {
    skip_on_cran()

    library(diversitree)
    pars <- c(1.5, 0.5, 1.0, 0.7, 0.7, 2.5, 0.5)
    names(pars) <- diversitree:::default.argnames.geosse()
    set.seed(5)
    phy <- NULL
    while( is.null( phy ) ){
        phy <- tree.geosse(pars, max.t=4, x0=0)
    }
    lik <- make.geosse(phy, phy$tip.state)
    diversitree.full <- lik(pars)

    states <- data.frame(phy$tip.state, phy$tip.state, row.names=names(phy$tip.state))
    states <- states[phy$tip.label,]
    names(pars) <- NULL
    model.vec <- numeric(115)
    model.vec[1:11] <- c(pars[1:3], pars[4:5], 0, pars[6], 0, pars[7], pars[5:4])
    phy$node.label <- NULL
    cache <- hisse:::ParametersToPassGeoHiSSE(phy, states[,1], f=c(1,1,1), model.vec, hidden.states=FALSE)
    geohisse.full <- hisse:::DownPassGeoHisse(phy=phy, cache=cache, hidden.states=FALSE, bad.likelihood=-1000000, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL)
    comparison <- identical(round(geohisse.full,4), round(diversitree.full,4))
    expect_true(comparison)
})


test_that("GeoHiSSE_test2", {
    skip_on_cran()

    library(diversitree)
    pars <- c(1.5, 0.5, 1.0, 0.7, 0.7, 2.5, 0.5)
    names(pars) <- diversitree:::default.argnames.geosse()
    set.seed(5)
    phy <- NULL
    while( is.null( phy ) ){
        phy <- tree.geosse(pars, max.t=4, x0=0)
    }
    lik <- make.geosse(phy, phy$tip.state)
    diversitree.full <- lik(pars)

    states <- data.frame(phy$tip.state, phy$tip.state, row.names=names(phy$tip.state))
    states <- states[phy$tip.label,]
    names(pars) <- NULL
    model.vec <- numeric(115)
    model.vec[1:11] <- c(pars[1:3], pars[4:5], 0, pars[6], 0, pars[7], pars[5:4])
    phy$node.label <- NULL
    cache <- hisse:::ParametersToPassGeoHiSSE(phy, states[,1], f=c(1,1,1), model.vec, hidden.states=TRUE)
    geohisse.full <- hisse:::DownPassGeoHisse(phy=phy, cache=cache, hidden.states=TRUE, bad.likelihood=-1000000, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL)
    comparison <- identical(round(geohisse.full,4), round(diversitree.full,4))
    expect_true(comparison)
})


test_that("GeoHiSSE_test3", {
    skip_on_cran()

    library(diversitree)
    pars <- c(1.5, 0.5, 1.0, 0.7, 0.4, 2.5, 0.5)
    names(pars) <- diversitree:::default.argnames.geosse()
    set.seed(9)
    phy <- NULL
    while( is.null( phy ) ){
        phy <- tree.geosse(pars, max.t=4, x0=0)
    }
    lik <- make.geosse(phy, phy$tip.state)
    diversitree.full <- lik(pars)

    states <- data.frame(phy$tip.state, phy$tip.state, row.names=names(phy$tip.state))
    states <- states[phy$tip.label,]
    names(pars) <- NULL
    model.vec <- numeric(115)
    model.vec[1:11] <- c(pars[1:3], pars[4:5], 0, pars[6], 0, pars[7], pars[5:4])
    phy$node.label <- NULL
    cache <- hisse:::ParametersToPassGeoHiSSE(phy, states[,1], f=c(1,1,1), model.vec, hidden.states=FALSE)
    geohisse.full <- hisse:::DownPassGeoHisse(phy=phy, cache=cache, hidden.states=TRUE, bad.likelihood=-1000000, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL)
    comparison <- identical(round(geohisse.full,4), round(diversitree.full,4))
    expect_true(comparison)
})


test_that("GeoHiSSE_test4", {
    skip_on_cran()

    library(diversitree)
    pars <- c(1.5, 0.5, 1.0, 0.7, 0.4, 2.5, 0.5)
    names(pars) <- diversitree:::default.argnames.geosse()
    set.seed(9)
    phy <- NULL
    while( is.null( phy ) ){
        phy <- tree.geosse(pars, max.t=4, x0=0)
    }
    lik <- make.geosse(phy, phy$tip.state)
    diversitree.full <- lik(pars)

    states <- data.frame(phy$tip.state, phy$tip.state, row.names=names(phy$tip.state))
    states <- states[phy$tip.label,]
    names(pars) <- NULL
    model.vec <- numeric(115)
    model.vec[1:11] <- c(pars[1:3], pars[4:5], 0, pars[6], 0, pars[7], pars[5:4])
    phy$node.label <- NULL
    cache <- hisse:::ParametersToPassGeoHiSSE(phy, states[,1], f=c(1,1,1), model.vec, hidden.states=TRUE)
    geohisse.full <- hisse:::DownPassGeoHisse(phy=phy, cache=cache, hidden.states=TRUE, bad.likelihood=-1000000, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL)
    comparison <- identical(round(geohisse.full,4), round(diversitree.full,4))
    expect_true(comparison)
})

#test_that("GeoHiSSE_test5", {
#    skip_on_cran()

#    library(diversitree)

    ## Define some parameters:
    ##         s01,  s0,  s1,  x0,  x1,  d0,  d1
#    pars1 <- c(1.5, 0.5, 1.0, 0.7, 0.7, 2.5, 0.5)
#    pars2 <- c(1.1, 1.6, 0.4, 0.6, 0.5, 0.6, 2.7)
#    rates.q <- 0.1 ## Transition between hidden rates.

    ## Check the likelihood for the hidden null model.
#    sim.pars <- SimulateGeoHiSSE(hidden.areas = 1, return.GeoHiSSE_pars = TRUE)
#    sim.pars$model.pars[,1] <- pars1
#    sim.pars$model.pars[,2] <- pars2
#    sim.pars$q.01[2,1] <- rates.q
#    sim.pars$q.01[1,2] <- rates.q
#    sim.pars$q.0[2,1] <- rates.q
#    sim.pars$q.0[1,2] <- rates.q
#    sim.pars$q.1[2,1] <- rates.q
#    sim.pars$q.1[1,2] <- rates.q

#    sim.get.par <- SimulateGeoHiSSE(pars=sim.pars, hidden.areas = 1, max.taxa = 500)
#    classe.pars <- sim.get.par$classe.pars
#    sim.data <- NULL
#    while( is.null(sim.data) ) sim.data <- tree.classe(pars=classe.pars, max.taxa = 500)

    ## Get the likelihood for the ClaSSE model:
    ## Here we assume root value is EQUAL
#    sim.data$node.label <- NULL
#    lik.classe <- make.classe(tree=sim.data, states=sim.data$tip.state, k=6)
#    classe.full <- lik.classe(pars=classe.pars, root=ROOT.FLAT, condition.surv = TRUE, root.p = NULL)

    ## Now the lik for the GeoHiSSE model.
    ## Our function starts from 0.
#    states <- sim.data$tip.state - 1
#    states.mat <- data.frame(states, states, row.names=names(states))
#    order.mat <- match(rownames(states.mat), sim.data$tip.label)
#    states.mat <- states.mat[order.mat,]
#    model.vec <- numeric(115)
#    order.pars1 <- c(pars1[c(2,3,1,4:5)],0,pars1[6],0,pars1[c(7,5:4)], rates.q,0,0,0, rates.q,0,0,0
#                   , rates.q,0,0,0)
#    order.pars2 <- c(pars2[c(2,3,1,4:5)],0,pars2[6],0,pars2[c(7,5:4)], rates.q,0,0,0, rates.q,0,0,0
#                   , rates.q,0,0,0)
#    order.pars <- c(order.pars1, order.pars2)
#    model.vec[1:46] <- order.pars

#    cache <- hisse:::ParametersToPassGeoHiSSE(sim.data, states.mat[,1], f=c(1,1,1), model.vec, hidden.states="TEST")
#    geohisse.full <- hisse:::DownPassGeoHisse(phy=sim.data, cache=cache, hidden.states=TRUE, bad.likelihood=-1000000
#                                   , condition.on.survival=TRUE, root.type="madfitz", root.p=c(1/3,1/3,1/3))
#    comparison <- identical(round(geohisse.full,4), round(classe.full,4))
#    expect_true(comparison)
#})

test_that("MuSSE_test1", {
    skip_on_cran()

    library(diversitree)
    pars <- c(.1,  .15,  .2,  # lambda 1, 2, 3
    .03, .045, .06, # mu 1, 2, 3
    .05, 0,         # q12, q13
    .05, .05,       # q21, q23
    0,   .05)       # q31, q32
    set.seed(2)
    phy <- NULL
    while( is.null( phy ) ){
        phy <- tree.musse(pars, 30, x0=1)
    }
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
    hisse.constrained <- hisse:::DownPassMusse(phy, cache, hidden.states=FALSE, root.type="madfitz", condition.on.survival=TRUE, root.p=NULL)
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
    phy <- NULL
    while( is.null( phy ) ){
        phy <- tree.musse(pars, 30, x0=1)
    }
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
    cache <- hisse:::ParametersToPassMuSSE(phy, states[,1], model.vec, f=c(1,1,1), hidden.states="TEST2")
    hisse.constrained <- hisse:::DownPassMusse(phy, cache, hidden.states=TRUE, root.type="madfitz", condition.on.survival=TRUE, root.p=NULL)
    comparison <- identical(round(hisse.constrained,4), round(diversitree.constrained,4))
    expect_true(comparison)
})


test_that("MuHiSSE_test1", {
    skip_on_cran()

    library(diversitree)
    pars <- c(.1,  .15,  .2, .1, # lambda 1, 2, 3, 4
    .03, .045, .06, 0.03, # mu 1, 2, 3, 4
    .05, .05, .00,        # q12, q13, q14
    .05, .00, .05,     # q21, q23, q24
    .05, .00, .05,     # q31, q32, q34
    .00, .05, .05)
    set.seed(2)
    phy <- NULL
    while( is.null( phy ) ){
        phy <- tree.musse(pars, 30, x0=1)
    }
    states <- phy$tip.state
    lik <- make.musse(phy, states, 4)
    #lik <- make.musse(phy, states, 3)
    diversitree.free = lik(pars)

    states <- data.frame(phy$tip.state, phy$tip.state,
    row.names=names(phy$tip.state))
    states <- states[phy$tip.label,]
    states.trans <- states
    for(i in 1:Ntip(phy)){
        if(states[i,1] == 1){
            states.trans[i,1] = 0
            states.trans[i,2] = 0
        }
        if(states[i,1] == 2){
            states.trans[i,1] = 0
            states.trans[i,2] = 1
        }
        if(states[i,1] == 3){
            states.trans[i,1] = 1
            states.trans[i,2] = 0
        }
        if(states[i,1] == 4){
            states.trans[i,1] = 1
            states.trans[i,2] = 1
        }
    }
    pars.hisse <- c(pars[1]+pars[5],pars[2]+pars[6],pars[3]+pars[7],pars[4]+pars[8],pars[5]/pars[1],pars[6]/pars[2],pars[7]/pars[3],pars[8]/pars[4], 0.05,0.05,0, 0.05,0,0.05, 0.05,0,.05, 0,0.05,.05)
    model.vec = rep(0,384)
    model.vec[1:20] = pars.hisse
    phy$node.label = NULL
    cache <- hisse:::ParametersToPassMuHiSSE(model.vec=model.vec, hidden.states=TRUE, nb.tip=Ntip(phy), nb.node=Nnode(phy), bad.likelihood=exp(-300), f=c(1,1,1,1), ode.eps=0)
    cache$psi <- 0
    gen <- hisse:::FindGenerations(phy)
    dat.tab <- hisse:::OrganizeData(states.trans, phy, f=c(1,1,1,1), hidden.states=TRUE)
    hisse.constrained <- hisse:::DownPassMuHisse(dat.tab, gen=gen, cache=cache, root.type="madfitz", condition.on.survival=TRUE, root.p=NULL)
    comparison <- identical(round(hisse.constrained,4), round(diversitree.free,4))
    expect_true(comparison)
})


test_that("MuHiSSE_test2", {
    skip_on_cran()

    library(diversitree)
    pars <- c(.1,  .15,  .2,  # lambda 1, 2, 3
    .03, .045, .06, # mu 1, 2, 3
    .05, 0,         # q12, q13
    .05, .05,       # q21, q23
    0,   .05)       # q31, q32
    set.seed(2)
    phy <- NULL
    while( is.null( phy ) ){
        phy <- tree.musse(pars, 30, x0=1)
    }
    states <- phy$tip.state
    lik <- make.musse(phy, states, 3)
    lik.base <- constrain(lik, lambda2 ~ lambda1, lambda3 ~ lambda1,
    mu2 ~ mu1, mu3 ~ mu1,
    q13 ~ 0, q21 ~ q12, q23 ~ q12, q31 ~ 0.03, q32 ~ q12)
    diversitree.constrained = lik.base(c(.1, .03, .05))

    states <- data.frame(phy$tip.state, phy$tip.state, row.names=names(phy$tip.state))
    states <- states[phy$tip.label,]
    states[states[,1]==3,] = 4
    states.trans <- states
    for(i in 1:Ntip(phy)){
        if(states[i,1] == 1){
            states.trans[i,1] = 0
            states.trans[i,2] = 0
        }
        if(states[i,1] == 2){
            states.trans[i,1] = 0
            states.trans[i,2] = 1
        }
        if(states[i,1] == 4){
            states.trans[i,1] = 1
            states.trans[i,2] = 0
        }
    }
    pars.hisse <- c(0.1+0.03,0.1+0.03,0.1+0.03,0,0.03/0.1,0.03/0.1,0.03/0.1,0,0.05,0,0, 0.05,0.05,0, 0.03,0.05,0, 0,0,0)
    model.vec = rep(0,384)
    model.vec[1:20] = pars.hisse
    phy$node.label = NULL
    cache <- hisse:::ParametersToPassMuHiSSE(model.vec=model.vec, hidden.states=TRUE, nb.tip=Ntip(phy), nb.node=Nnode(phy), bad.likelihood=exp(-300), f=c(1,1,1,0), ode.eps=0)
    cache$psi <- 0
    gen <- hisse:::FindGenerations(phy)
    dat.tab <- hisse:::OrganizeData(states.trans, phy, f=c(1,1,1,0), hidden.states=TRUE)
    muhisse.constrained <- hisse:::DownPassMuHisse(dat.tab, gen=gen, cache=cache, root.type="madfitz", condition.on.survival=TRUE, root.p=NULL)
    comparison <- identical(round(muhisse.constrained,4), round(diversitree.constrained,4))
    expect_true(comparison)
})


test_that("MiSSE_test1", {
    skip_on_cran()

    phy <- read.tree("whales_Steemanetal2009.tre")
    gen <- hisse:::FindGenerations(phy)
    dat.tab <- hisse:::OrganizeDataMiSSE(phy=phy, f=1, hidden.states=1)
    nb.tip <- Ntip(phy)
    nb.node <- phy$Nnode
    model.vec <- c(0.103624, 5.207178e-09, rep(0,52), 0)
    cache = hisse:::ParametersToPassMiSSE(model.vec=model.vec, hidden.states=1, fixed.eps=NULL, nb.tip=nb.tip, nb.node=nb.node, bad.likelihood=exp(-300), ode.eps=0)
    logl.one.rate <- hisse:::DownPassMisse(dat.tab=dat.tab, cache=cache, gen=gen, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL)
    right.logl <- -277.6942
    comparison <- round(logl.one.rate,4) == round(right.logl,4)
    expect_true(comparison)
})


test_that("MiSSE_test2", {
    skip_on_cran()

    phy <- read.tree("whales_Steemanetal2009.tre")
    gen <- hisse:::FindGenerations(phy)
    dat.tab <- hisse:::OrganizeDataMiSSE(phy=phy, f=1, hidden.states=3)
    nb.tip <- Ntip(phy)
    nb.node <- phy$Nnode
    model.vec <- c(0.103624, 5.207178e-09, 0.103624, 5.207178e-09, 0.103624, 5.207178e-09, rep(0,46), 1, 0, 0)
    cache = hisse:::ParametersToPassMiSSE(model.vec=model.vec, hidden.states=3, fixed.eps=NULL, nb.tip=nb.tip, nb.node=nb.node, bad.likelihood=exp(-300), ode.eps=0)
    logl.three.rate <- hisse:::DownPassMisse(dat.tab=dat.tab, cache=cache, gen=gen, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL)
    right.logl <- -277.6942
    comparison <- round(logl.three.rate,4) == round(right.logl,4)
    expect_true(comparison)
})


test_that("MiSSE_test3", {
    skip_on_cran()

    load("CID2_MiSSE_check.Rsave")
    
    fix.type <- NULL
    nb.tip <- Ntip(phy)
    nb.node <- phy$Nnode
    gen <- hisse:::FindGenerations(phy)
    
    dat.tab <- hisse:::OrganizeDataMiSSE(phy=phy, f=1, hidden.states=2)
    model.vec <- c(pp.hisse$solution[1], pp.hisse$solution[3], pp.hisse$solution[13], pp.hisse$solution[15], rep(0,48),  pp.hisse$solution[17],0)
    cache = hisse:::ParametersToPassMiSSE(model.vec=model.vec, hidden.states=2, fixed.eps=NULL, nb.tip=nb.tip, nb.node=nb.node, bad.likelihood=exp(-300), ode.eps=0)#
    gen <- hisse:::FindGenerations(phy)
    MiSSE.logL <- hisse:::DownPassMisse(dat.tab=dat.tab, cache=cache, gen=gen, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, node=NULL, fix.type=NULL)
    
    #Part 2: corHMM loglik:
    library(corHMM)
    data <- data.frame(taxon=names(phy$tip.state), phy$tip.state, stringsAsFactors=FALSE)
    char.logL <- corHMM(phy, data, rate.cat=1, model = "ER", node.states = "none", fixed.nodes=FALSE, p=pp.hisse$solution[17], root.p="maddfitz")
    tot.logL <- char.logL$loglik + MiSSE.logL
    
    comparison <- identical(round(pp.hisse$loglik,3), round(tot.logL,3))
    expect_true(comparison)
})


test_that("BiSSE_fHiSSE_test1", {
    skip_on_cran()

    library(diversitree)
    pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
    set.seed(4)
    phy <- NULL
    while( is.null( phy ) ){
        phy <- tree.bisse(pars, max.t=30, x0=0)
    }
    lik <- make.bisse(phy, phy$tip.state, sampling.f=c(.4,.6))
    diversitree.full <- lik(pars)

    hidden.states=FALSE
    states <- data.frame(phy$tip.state, phy$tip.state, row.names=names(phy$tip.state))
    states <- states[phy$tip.label,]
    gen <- hisse:::FindGenerations(phy)
    dat.tab <- hisse:::OrganizeDataHiSSE(states, phy=phy, f=c(.4,.6), hidden.states=FALSE)
    pars.bisse <- c(0.1+0.03, 0.2+0.03, 0.03/0.1, 0.03/0.2, 0.01, 0.01)
    model.vec <- numeric(48)
    model.vec[1:6] = pars.bisse
    phy$node.label = NULL
    cache <- hisse:::ParametersToPassfHiSSE(model.vec, hidden.states=hidden.states, nb.tip=Ntip(phy), nb.node=Nnode(phy),  bad.likelihood=-300, f=c(.4,.6), ode.eps=0)
    cache$psi <- 0
    hisse.full <- hisse:::DownPassHiSSE(dat.tab, gen, cache, root.type="madfitz", condition.on.survival=TRUE, root.p=NULL, fossil.taxa=NULL)
    comparison <- identical(round(hisse.full,4), round(diversitree.full,4))
    expect_true(comparison)

})


test_that("BiSSE_fHiSSE_test2", {
    skip_on_cran()

    library(diversitree)
    pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
    set.seed(4)
    phy <- NULL
    while( is.null( phy ) ){
        phy <- tree.bisse(pars, max.t=30, x0=0)
    }
    lik <- make.bisse(phy, phy$tip.state, sampling.f=c(.4,.6))
    diversitree.full <- lik(pars)

    hidden.states=TRUE
    states <- data.frame(phy$tip.state, phy$tip.state, row.names=names(phy$tip.state))
    states <- states[phy$tip.label,]
    gen <- hisse:::FindGenerations(phy)
    dat.tab <- hisse:::OrganizeDataHiSSE(states, phy=phy, f=c(.4,.6), hidden.states=hidden.states)
    pars.hisse <- c(0.1+0.03, 0.2+0.03, 0.03/0.1, 0.03/0.2, 0.01, 0.01, 0.1, rep(0,5), 0.1+0.03, 0.2+0.03, 0.03/0.1, 0.03/0.2, 0.01, 0.01, 0.1, rep(0,5))
    model.vec <- numeric(48)
    model.vec[1:24] = pars.hisse
    phy$node.label = NULL
    cache = hisse:::ParametersToPassfHiSSE(model.vec, hidden.states=hidden.states, nb.tip=Ntip(phy), nb.node=Nnode(phy), bad.likelihood=-300, f=c(.4,.6), ode.eps=0)
    cache$psi <- 0
    hisse.full <- hisse:::DownPassHiSSE(dat.tab, gen, cache, root.type="madfitz", condition.on.survival=TRUE, root.p=NULL)
    comparison <- identical(round(hisse.full,4), round(diversitree.full,4))
    expect_true(comparison)
})


test_that("HiSSE_Null_Four_fHiSSE_test", {
    skip_on_cran()
    
    library(diversitree)
    pars <- c(0.1, 0.1, 0.03, 0.03, 0.01, 0.01)
    set.seed(4)
    phy <- NULL
    while( is.null( phy ) ){
        phy <- tree.bisse(pars, max.t=30, x0=0)
    }
    lik <- make.bisse(phy, phy$tip.state)
    diversitree.full <- lik(pars)
    
    states <- data.frame(phy$tip.state, phy$tip.state, row.names=names(phy$tip.state))
    states <- states[phy$tip.label,]
    pars.hisse.null <- c(0.1+0.03, 0.2+0.03, 0.1+0.05, 0.2+0.05, 0.1+0.03, 0.2+0.03, 0.1+0.05, 0.2+0.05, 0.03/.1, 0.03/.2, 0.05/.1, 0.05/.2, 0.03/.1, 0.03/.2, 0.05/.1, 0.05/.2, rep(0.01, 32))
    model.vec = pars.hisse.null
    phy$node.label = NULL
    cache <- hisse:::ParametersToPassNull(phy, states[,1], model.vec, f=c(1,1))
    hisse.nullOG.full <- hisse:::DownPassNull(phy, cache, root.type="madfitz", condition.on.survival=TRUE, root.p=NULL)
    
    hidden.states=TRUE
    states <- data.frame(phy$tip.state, phy$tip.state, row.names=names(phy$tip.state))
    states <- states[phy$tip.label,]
    gen <- hisse:::FindGenerations(phy)
    dat.tab <- hisse:::OrganizeDataHiSSE(states, phy=phy, f=c(1,1), hidden.states=hidden.states)
    pars.hisse <- c(0.1+0.03, 0.1+0.03, 0.03/0.1, 0.03/0.1, 0.01, 0.01, 0.01, rep(0.01,5), 0.2+0.03, 0.2+0.03, 0.03/.2, 0.03/.2, 0.01, 0.01, 0.01, rep(0.01,5), 0.1+0.05, 0.1+0.05, 0.05/.1, 0.05/.1, 0.01, 0.01, 0.01, rep(0.01,5), 0.2+0.05, 0.2+0.05, 0.05/.2, 0.05/.2, 0.01, 0.01, 0.01, rep(0.01,5))
    model.vec <- numeric(48)
    model.vec[1:48] = pars.hisse
    phy$node.label = NULL
    cache = hisse:::ParametersToPassfHiSSE(model.vec, hidden.states=hidden.states, nb.tip=Ntip(phy), nb.node=Nnode(phy), bad.likelihood=-300, f=c(1,1), ode.eps=0)
    cache$psi <- 0
    hisse.null.full <- hisse:::DownPassHiSSE(dat.tab, gen, cache, root.type="madfitz", condition.on.survival=TRUE, root.p=NULL)
    comparison <- identical(round(hisse.nullOG.full,4), round(hisse.null.full,4))
    
    expect_true(comparison)
})


test_that("GeoSSE_fGeoSSSE_test1", {
    skip_on_cran()

    library(diversitree)
    pars <- c(1.5, 0.5, 1.0, 0.7, 0.4, 2.5, 0.5)
    names(pars) <- diversitree:::default.argnames.geosse()
    set.seed(9)
    phy <- tree.geosse(pars, max.t=4, x0=0)
    lik <- make.geosse(phy, phy$tip.state, sampling.f=c(.2,.2,.2))
    diversitree.full <- lik(pars)
    
    states <- data.frame(phy$tip.state, phy$tip.state, row.names=names(phy$tip.state))
    states <- states[phy$tip.label,]
    names(pars) <- NULL
    model.vec <- numeric(380)
    model.vec[1:11] <- c(1.5+0.7, 0.5+0.4, sum(1.5, 0.5, 1.0), pars[4]/pars[1], pars[5]/pars[2], 0, pars[6], 0, pars[7], 0,0)
    phy$node.label <- NULL
    cache.slim <- hisse:::ParametersToPassGeoHiSSEfast(model.vec=model.vec, hidden.states=FALSE, nb.tip=Ntip(phy), nb.node=Nnode(phy), bad.likelihood=exp(-300), ode.eps=0)
    gen <- hisse:::FindGenerations(phy)
    dat.tab <- hisse:::OrganizeDataGeo(states[,1], phy, c(.2,.2,.2), hidden.states=FALSE)
    geohisse.new <- hisse:::DownPassGeoHissefast(dat.tab, gen=gen, cache=cache.slim, root.type="madfitz", condition.on.survival=TRUE, root.p=NULL)
    
    model.vec[1:11] <- c(1.5, 0.5, 1.0, pars[4], pars[5], 0, pars[6], 0, pars[7], pars[5], pars[4])
    cache <- hisse:::ParametersToPassGeoHiSSE(phy, states[,1], f=c(.2,.2,.2), model.vec, hidden.states=FALSE)
    geohisse.full <- hisse:::DownPassGeoHisse(phy=phy, cache=cache, hidden.states=TRUE, bad.likelihood=-1000000, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL)
    #comparison <- identical(round(geohisse.full,4), round(diversitree.full,4), round(geohisse.new,4))
    comparison <- identical(round(geohisse.new,4), round(diversitree.full,4))
    
    expect_true(comparison)
})


test_that("GeoSSE_fGeoSSSE_test2", {
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
    model.vec <- numeric(380)
    model.vec[1:76] <- c(0.07612170+0.02869795, 0.07612170+0.02869795, 0, 0.02869795/0.07612170, 0.02869795/0.07612170, 0, 0.02727061, 0, 0.03946784, 0.02869795, 0.02869795, 0.01055818, rep(0, 8), 0.01055818, rep(0, 8), 0.01055818, rep(0, 8), 0.18861747+0.02869795, 0.18861747+0.02869795, 0, 0.02869795/0.18861747, 0.02869795/0.18861747, 0, 0.02727061, 0, 0.03946784, 0.02869795, 0.02869795, 0.01055818, rep(0, 8), 0.01055818, rep(0, 8), 0.01055818, rep(0, 8))
    phy$node.label <- NULL
    cache.slim <- hisse:::ParametersToPassGeoHiSSEfast(model.vec=model.vec, hidden.states=TRUE, nb.tip=Ntip(phy), nb.node=Nnode(phy), bad.likelihood=exp(-300), ode.eps=0)
    gen <- hisse:::FindGenerations(phy)
    dat.tab <- hisse:::OrganizeDataGeo(states[,1], phy, c(1,1,1), hidden.states=TRUE)
    geohisse.new <- hisse:::DownPassGeoHissefast(dat.tab, gen=gen, cache=cache.slim, root.type="madfitz", root.p=NULL, condition.on.survival=TRUE)
    
    model.vec <- c(0.07612170, 0.07612170, 0.07612170, 0.02869795, 0.02869795, 0.00000000, 0.02727061,
    0.00000000, 0.03946784, 0.00000000, 0.00000000, 0.01055818, 0.00000000, 0.00000000,
    0.00000000, 0.01055818, 0.00000000, 0.00000000, 0.00000000, 0.01055818, 0.00000000,
    0.00000000, 0.00000000, 0.18861747, 0.18861747, 0.18861747, 0.02869795, 0.02869795,
    0.00000000, 0.02727061, 0.00000000, 0.03946784, 0.00000000, 0.00000000, 0.01055818,
    0.00000000, 0.00000000, 0.00000000, 0.01055818, 0.00000000, 0.00000000, 0.00000000,
    0.01055818, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000,
    0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000,
    0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000,
    0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000,
    0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000,
    0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000,
    0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000,
    0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000,
    0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000,
    0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000, 0.00000000,
    0.00000000, 0.00000000, 0.00000000)
    cache <- hisse:::ParametersToPassGeoHiSSE(phy, states[,1], f=c(1,1,1), model.vec, hidden.states=TRUE)
    geohisse.full <- hisse:::DownPassGeoHisse(phy=phy, cache=cache, hidden.states=TRUE, bad.likelihood=-1000000, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL)
    comparison <- identical(round(geohisse.full,4), round(geohisse.new,4))
    
    expect_true(comparison)
})


test_that("GeoSSE_fGeoSSSE_test3", {
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
    lik.base <- constrain(lik, lambda2 ~ lambda1, lambda3 ~ 0,
    mu2 ~ mu1, mu3 ~ 0,
    q13 ~ 0, q21 ~ q12, q23 ~ q12, q31 ~ 0.03, q32 ~ q12)
    diversitree.constrained = lik.base(c(.1, .03, .05))
    diversitree.full = lik(pars)

    states <- data.frame(phy$tip.state, phy$tip.state, row.names=names(phy$tip.state))
    states <- states[phy$tip.label,]
    states[states[,1]==3,] = 4
    pars.hisse <- c(0.1, 0.1, 0, 0.03, 0.03, 0.0, 0.05, 0, 0.05, 0.05, 0, 0.05)
    model.vec = rep(0,120)
    model.vec[1:12] = pars.hisse
    phy$node.label = NULL
    cache <- hisse:::ParametersToPassMuSSE(phy, states[,1], model.vec, f=c(1,1,1), hidden.states="TEST1")
    hisse.constrained <- hisse:::DownPassMusse(phy, cache, hidden.states=FALSE, root.type="madfitz", condition.on.survival=TRUE, root.p=NULL)
    #comparison <- identical(round(hisse.constrained,4), round(diversitree.constrained,4))

    states.new <- numeric(Ntip(phy))
    for(i in 1:Ntip(phy)){
        if(states[i,1]==1){states.new[i]=1}
        if(states[i,1]==2){states.new[i]=2}
        if(states[i,1]==4){states.new[i]=0}
    }

    states <- cbind(states, states.new)
    pars.hisse <- c(0.1+0.03, 0.1+0.03, 0, 0.03/0.1, 0.03/0.1, 0.05, 0, 0.05, 0.05, 0, 0.05)
    model.vec = rep(0,380)
    model.vec[1:11] = pars.hisse
    phy$node.label = NULL
    cache.slim <- hisse:::ParametersToPassGeoHiSSEfast(model.vec=model.vec, hidden.states=FALSE, assume.cladogenetic=FALSE, nb.tip=Ntip(phy), nb.node=Nnode(phy), bad.likelihood=exp(-300), ode.eps=0)
    gen <- hisse:::FindGenerations(phy)
    dat.tab <- hisse:::OrganizeDataGeo(states[,3], phy, c(1,1,1), hidden.states=FALSE)
    geohisse.new <- hisse:::DownPassGeoHissefast(dat.tab, gen=gen, cache=cache.slim, root.type="madfitz", condition.on.survival=TRUE, root.p=NULL)
    #comparison <- identical(round(hisse.constrained,4), round(diversitree.constrained,4), round(geohisse.new,4))
    comparison <- identical(round(geohisse.new,4), round(diversitree.constrained,4))
    
    expect_true(comparison)
})


test_that("MiSSE_fossil_test1", {
    skip_on_cran()

    #Tests the loglikehood for the starting values code against Stadler (2010)
    set.seed(42)
    phy <- TreeSim::sim.bd.taxa(n = 100, numbsim = 1, lambda = 0.3, mu = 0.2)[[1]]
    f <- hisse:::GetFossils(phy, psi=0.05)
    pp <- hisse:::ProcessSimSample(phy, f)
    
    dat.tab <- hisse:::OrganizeDataMiSSE(phy=pp$phy, f=1, hidden.states=1)
    edge_details <- hisse:::GetEdgeDetails(phy, includes.intervals=FALSE, intervening.intervals=NULL)
    fossil.taxa <- edge_details$tipward_node[which(edge_details$type == "extinct_tip")]
    #This is for starting values:
    fossil.ages <- dat.tab$TipwardAge[which(dat.tab$DesNode %in% fossil.taxa)]
    #Drop all k.samples to get split times:
    k.sample.tip.no <- grep("Ksamp*", x=phy$tip.label)
    phy.no.k <- drop.tip(pp$phy, k.sample.tip.no)
    split.times <- paleotree:::dateNodes(phy.no.k, rootAge=max(node.depth.edgelength(phy.no.k)))[-c(1:Ntip(phy.no.k))]
    n <- Ntip(phy.no.k)-length(fossil.taxa)
    m <- length(fossil.taxa)
    x_times <- split.times
    y_times <- fossil.ages
    k <- dim(pp$k.samples)[1]
    
    starting.point.code <- hisse:::starting.point.generator.fossils(n.tax=n, k=1, samp.freq.tree=1, q.div=5, fossil.taxa=fossil.taxa, fossil.ages=fossil.ages, no.k.samples=k, split.times=split.times, get.likelihood=TRUE)
    
    rho=1
    lambda <- starting.point.code[1]
    mu <-  starting.point.code[2]
    psi <- starting.point.code[3]
    
    logLikLogSpace <- (((n+m-2) * log(lambda)) + ((k+m) * log(psi))) - log(1-exp(hisse:::p_0(max(x_times),lambda,mu,psi=0,rho)))*2 + hisse:::p_one(max(x_times), lambda, mu, psi, rho) + sum(hisse:::p_one(x_times, lambda,mu,psi,rho)) + (sum(hisse:::p_0(y_times,lambda,mu,psi,rho)) - sum(hisse:::p_one(y_times,lambda,mu,psi,rho)))

    comparison <- identical(round(-starting.point.code[4], 4), round(logLikLogSpace, 4))

    expect_true(comparison)
})


test_that("MiSSE_fossil_test2", {
    skip_on_cran()

    #Tests the loglikehood for MiSSE when there are only m+k samples
    set.seed(42)
    phy <- TreeSim::sim.bd.taxa(n = 100, numbsim = 1, lambda = 0.3, mu = 0.2)[[1]]
    f <- hisse:::GetFossils(phy, psi=0.05)
    pp <- hisse:::ProcessSimSample(phy, f)

    dat.tab <- hisse:::OrganizeDataMiSSE(phy=pp$phy, f=1, hidden.states=1)
    edge_details <- hisse:::GetEdgeDetails(phy=pp$phy, includes.intervals=FALSE, intervening.intervals=NULL)
    fossil.taxa <- edge_details$tipward_node[which(edge_details$type == "extinct_tip")]
    #This is for starting values:
    fossil.ages <- dat.tab$TipwardAge[which(dat.tab$DesNode %in% fossil.taxa)]
    #Drop all k.samples to get split times:
    k.sample.tip.no <- grep("Ksamp*", x=phy$tip.label)
    phy.no.k <- drop.tip(pp$phy, k.sample.tip.no)
    split.times <- paleotree:::dateNodes(phy.no.k, rootAge=max(node.depth.edgelength(phy.no.k)))[-c(1:Ntip(phy.no.k))]
    n <- Ntip(phy.no.k)-length(fossil.taxa)
    m <- length(fossil.taxa)
    x_times <- split.times
    y_times <- fossil.ages
    k <- dim(pp$k.samples)[1]

    starting.point.code <- hisse:::starting.point.generator.fossils(n.tax=n, k=1, samp.freq.tree=1, q.div=5, fossil.taxa=fossil.taxa, fossil.ages=fossil.ages, no.k.samples=k, split.times=split.times, get.likelihood=TRUE)

    rho=1
    lambda <- starting.point.code[1]
    mu <-  starting.point.code[2]
    psi <- starting.point.code[3]

    logLikLogSpace <- (((n+m-2) * log(lambda)) + ((k+m) * log(psi))) - log(1-exp(hisse:::p_0(max(x_times),lambda,mu,psi=0,rho)))*2 + hisse:::p_one(max(x_times), lambda, mu, psi, rho) + sum(hisse:::p_one(x_times, lambda,mu,psi,rho)) + (sum(hisse:::p_0(y_times,lambda,mu,psi,rho)) - sum(hisse:::p_one(y_times,lambda,mu,psi,rho)))

    phy <- hisse:::AddKNodes(pp$phy, pp$k.samples)
    nb.tip <- Ntip(phy)
    nb.node <- phy$Nnode
    dat.tab <- hisse:::OrganizeDataMiSSE(phy=phy, f=1, hidden.states=1)
    model.vec <- c(lambda+mu, mu/lambda, rep(0,51))
    cache = hisse:::ParametersToPassMiSSE(model.vec=model.vec, hidden.states=1, fixed.eps=NULL, nb.tip=nb.tip, nb.node=nb.node, bad.likelihood=exp(-300), ode.eps=0)#
    cache$psi = psi
    gen <- hisse:::FindGenerations(phy)
    k.samples <- hisse:::GetKSampleMRCA(phy, k.samples)
    MiSSE.logL <- hisse:::DownPassMisse(dat.tab=dat.tab, cache=cache, gen=gen, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, fossil.taxa=fossil.taxa, node=k.samples$node, fix.type=k.samples$type)
    comparison <- identical(round(logLikLogSpace,3), round(MiSSE.logL,3))

    expect_true(comparison)
})


test_that("MiSSE_fossil_test3", {
    skip_on_cran()

    #Tests the loglikehood when there are only m samples
    set.seed(42)
    phy <- TreeSim::sim.bd.taxa(n = 100, numbsim = 1, lambda = 0.3, mu = 0.2)[[1]]
    f <- hisse:::GetFossils(phy, psi=0.05)
    pp <- hisse:::ProcessSimSample(phy, f)

    dat.tab <- hisse:::OrganizeDataMiSSE(phy=pp$phy, f=1, hidden.states=2)
    edge_details <- hisse:::GetEdgeDetails(phy=pp$phy, intervening.intervals=strat.cache$intervening.intervals)
    fossil.taxa <- edge_details$tipward_node[which(edge_details$type == "extinct_tip")]
    #This is for starting values:
    fossil.ages <- dat.tab$TipwardAge[which(dat.tab$DesNode %in% fossil.taxa)]
    #Drop all k.samples to get split times:
    k.sample.tip.no <- grep("Ksamp*", x=phy$tip.label)
    phy.no.k <- drop.tip(pp$phy, k.sample.tip.no)
    split.times <- paleotree:::dateNodes(phy.no.k, rootAge=max(node.depth.edgelength(phy.no.k)))[-c(1:Ntip(phy.no.k))]
    n <- Ntip(phy.no.k)-length(fossil.taxa)
    m <- length(fossil.taxa)
    x_times <- split.times
    y_times <- fossil.ages
    k <- 0

    starting.point.code <- hisse:::starting.point.generator.fossils(n.tax=n, k=1, samp.freq.tree=1, q.div=5, fossil.taxa=fossil.taxa, fossil.ages=fossil.ages, no.k.samples=k, split.times=split.times, get.likelihood=TRUE)

    rho=1
    lambda <- starting.point.code[1]
    mu <-  starting.point.code[2]
    psi <- starting.point.code[3]

    logLikLogSpace <- (((n+m-2) * log(lambda)) + ((k+m) * log(psi))) - log(1-exp(hisse:::p_0(max(x_times),lambda,mu,psi=0,rho)))*2 + hisse:::p_one(max(x_times), lambda, mu, psi, rho) + sum(hisse:::p_one(x_times, lambda,mu,psi,rho)) + (sum(hisse:::p_0(y_times,lambda,mu,psi,rho)) - sum(hisse:::p_one(y_times,lambda,mu,psi,rho)))

    phy <- pp$phy
    nb.tip <- Ntip(phy)
    nb.node <- phy$Nnode
    dat.tab <- hisse:::OrganizeDataMiSSE(phy=phy, f=1, hidden.states=2)
    model.vec <- c(lambda+mu, mu/lambda, lambda+mu, mu/lambda, rep(0,48), 0.01, psi)
    cache = hisse:::ParametersToPassMiSSE(model.vec=model.vec, hidden.states=1, fixed.eps=NULL, nb.tip=nb.tip, nb.node=nb.node, bad.likelihood=exp(-300), ode.eps=0)#
    gen <- hisse:::FindGenerations(phy)
    k.samples <- NULL
    MiSSE.logL <- hisse:::DownPassMisse(dat.tab=dat.tab, cache=cache, gen=gen, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, fossil.taxa=fossil.taxa, node=NULL, fix.type=NULL)
    comparison <- identical(round(logLikLogSpace,3), round(MiSSE.logL,3))

    expect_true(comparison)
})


test_that("Placing_m_k_fossils", {
    skip_on_cran()
    
    set.seed(42)
    phy <- TreeSim::sim.bd.taxa(n = 200, numbsim = 1, lambda = .3, mu = .2)[[1]]
    f <- hisse:::GetFossils(phy, psi=0.1)
    pp <- hisse:::ProcessSimSample(phy, f)

    extinct.samples <- f[which(f$fossiltype_long=="extinct_terminal" | f$fossiltype_long=="extinct_internal"),]
    extinct.samples <- extinct.samples[which(extinct.samples$has_sampled_descendant == FALSE),]
    k.samples = pp$k.samples

    phy.k <- hisse:::AddKNodes(pp$phy, pp$k.samples)
    dat.tab <- hisse:::OrganizeDataMiSSE(phy=phy.k, f=1, hidden.states=1)
    edge_details <- hisse:::GetEdgeDetails(phy=phy.k, intervening.intervals=strat.cache$intervening.intervals)
    fossil.taxa <- edge_details$tipward_node[which(edge_details$type == "extinct_tip")]
    fossil.ages <- dat.tab$TipwardAge[which(dat.tab$branch.type == 1)]
    k.ages <- dat.tab$TipwardAge[which(dat.tab$branch.type == 2)]
    split.times <- paleotree::dateNodes(phy.k, rootAge=max(node.depth.edgelength(phy.k)))

    tot.diff.extinct <- sum(as.numeric(extinct.samples[,5]) - fossil.ages)
    tot.diff.ksample <- sum(as.numeric(k.samples[,3]) - k.ages)
    
    comparison1 <- identical(round(tot.diff.extinct,10), 0)
    comparison2 <- identical(round(tot.diff.ksample,10), 0)
    
    expect_true(comparison1)
    expect_true(comparison2)

})


test_that("HiSSE_fossil_test1", {
    skip_on_cran()
    
    library(diversitree)
    pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
    set.seed(4)
    phy <- NULL
    while( is.null( phy ) ){
        phy <- tree.bisse(pars, max.t=30, x0=0, include.extinct=TRUE)
    }
    #h <- history.from.sim.discrete(phy, 0:1)
    #plot(h, phy)
    k.samples <- data.frame(taxon1="sp12", taxon2="sp12", timefrompresent=3.164384, state=1, stringsAsFactors=FALSE)
    
    hidden.states=FALSE
    
    phy.k <- hisse:::AddKNodes(phy, k.samples)
    fix.type <- hisse:::GetKSampleMRCA(phy.k, k.samples)
    nb.tip <- Ntip(phy.k)
    nb.node <- phy.k$Nnode
    gen <- hisse:::FindGenerations(phy.k)
    
    data <- data.frame(taxon=names(phy$tip.state), phy$tip.state, stringsAsFactors=FALSE)
    data <- hisse:::AddKData(data, k.samples)
    data.new <- data.frame(data[,2], data[,2], row.names=data[,1])
    data.new <- data.new[phy.k$tip.label,]
    
    dat.tab <- hisse:::OrganizeDataHiSSE(data.new, phy=phy.k, f=c(1,1), hidden.states=FALSE)
    edge_details <- hisse:::GetEdgeDetails(phy.k, includes.intervals=FALSE, intervening.intervals=NULL)
    fossil.taxa <- edge_details$tipward_node[which(edge_details$type == "extinct_tip")]
    pars.bisse <- c(0.1+0.03, 0.1+0.03, 0.03/0.1, 0.03/0.1, 0.01, 0.01)
 
    model.vec <- numeric(48)
    model.vec[1:6] = pars.bisse
    phy$node.label = NULL
    cache <- hisse:::ParametersToPassfHiSSE(model.vec, hidden.states=hidden.states, nb.tip=Ntip(phy.k), nb.node=Nnode(phy.k), bad.likelihood=-300, f=c(1,1), ode.eps=0)
    cache$psi <- 0.01
    hisse.full <- hisse:::DownPassHiSSE(dat.tab, gen, cache, root.type="madfitz", condition.on.survival=TRUE, root.p=NULL, node=fix.type$node, state=fix.type$state, fossil.taxa=fossil.taxa, fix.type=fix.type$type)
    
    ## Trait independent model should be loglik_tree + loglik_character ##
    
    #Part 1: MiSSE loglik:
    dat.tab <- hisse:::OrganizeDataMiSSE(phy=phy.k, f=1, hidden.states=1)
    model.vec <- c(0.1+0.03, 0.03/0.1, rep(0,51))
    cache = hisse:::ParametersToPassMiSSE(model.vec=model.vec, hidden.states=1, fixed.eps=NULL, nb.tip=nb.tip, nb.node=nb.node, bad.likelihood=exp(-300), ode.eps=0)#
    cache$psi <- 0.01
    edge_details <- hisse:::GetEdgeDetails(phy=phy.k, includes.intervals=FALSE, intervening.intervals=NULL)
    fossil.taxa <- edge_details$tipward_node[which(edge_details$type == "extinct_tip")]
    #This is for starting values:
    gen <- hisse:::FindGenerations(phy.k)
    MiSSE.logL <- hisse:::DownPassMisse(dat.tab=dat.tab, cache=cache, gen=gen, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, fossil.taxa=fossil.taxa, node=fix.type$node, fix.type=fix.type$type)
    
    #Part 2: corHMM loglik:
    library(corHMM)
    char.logL <- corHMM(phy.k, data, rate.cat=1, model = "ER", node.states = "none", fixed.nodes=FALSE, p=0.01, root.p="maddfitz")
    tot.logL <- char.logL$loglik + MiSSE.logL
    
    comparison1 <- identical(round(hisse.full,3), round(tot.logL,3))
    expect_true(comparison1)

    #Same test, but fossils removed:
    extinct.tip.no <- grep("ex*", x=phy$tip.label)
    phy.extant <- drop.tip(phy, phy$tip.label[extinct.tip.no])
    hidden.states=FALSE
    gen <- hisse:::FindGenerations(phy.extant)
    
    dat.tab <- hisse:::OrganizeDataHiSSE(data.new[c(9:15,17:18),], phy=phy.extant, f=c(1,1), hidden.states=FALSE)
    pars.bisse <- c(0.1+0.03, 0.1+0.03, 0.03/0.1, 0.03/0.1, 0.01, 0.01)
    model.vec <- numeric(48)
    model.vec[1:6] = pars.bisse
    phy$node.label = NULL
    cache <- hisse:::ParametersToPassfHiSSE(model.vec, hidden.states=hidden.states, nb.tip=Ntip(phy.extant), nb.node=Nnode(phy.extant), bad.likelihood=-300, f=c(1,1), ode.eps=0)
    cache$psi <- 0
    hisse.full <- hisse:::DownPassHiSSE(dat.tab, gen, cache, root.type="madfitz", condition.on.survival=TRUE, root.p=NULL, fossil.taxa=NULL)

    dat.tab <- hisse:::OrganizeDataMiSSE(phy=phy.extant, f=1, hidden.states=1)
    model.vec <- c(0.1+0.03, 0.03/0.1, rep(0,51))
    cache = hisse:::ParametersToPassMiSSE(model.vec=model.vec, hidden.states=1, fixed.eps=NULL, nb.tip=Ntip(phy.extant), nb.node=Nnode(phy.extant), bad.likelihood=exp(-300), ode.eps=0)#
    cache$psi <- 0.0
    gen <- hisse:::FindGenerations(phy.extant)
    MiSSE.logL <- hisse:::DownPassMisse(dat.tab=dat.tab, cache=cache, gen=gen, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, fossil.taxa=NULL, node=NULL, fix.type=NULL)

    library(corHMM)
    char.logL <- corHMM(phy.extant, data[c(9:15,17:18),], rate.cat=1, model = "ER", node.states = "none", fixed.nodes=FALSE, p=0.01, root.p="maddfitz")
    tot.logL <- char.logL$loglik + MiSSE.logL

    comparison2 <- identical(round(hisse.full,3), round(tot.logL,3))
    expect_true(comparison2)
})


test_that("HiSSE_fossil_test2", {
    skip_on_cran()
    
    #Tests when there are no ksamples
    library(diversitree)
    pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
    set.seed(4)
    phy <- NULL
    while( is.null( phy ) ){
        phy <- tree.bisse(pars, max.t=30, x0=0, include.extinct=TRUE)
    }
    #h <- history.from.sim.discrete(phy, 0:1)
    #plot(h, phy)
    
    hidden.states=FALSE
    
    fix.type <- NULL
    nb.tip <- Ntip(phy)
    nb.node <- phy$Nnode
    gen <- hisse:::FindGenerations(phy)
    
    data <- data.frame(taxon=names(phy$tip.state), phy$tip.state, stringsAsFactors=FALSE)
    data.new <- data.frame(data[,2], data[,2], row.names=data[,1])
    data.new <- data.new[phy$tip.label,]
    
    dat.tab <- hisse:::OrganizeDataHiSSE(data.new, phy=phy, f=c(1,1), hidden.states=FALSE)
    edge_details <- hisse:::GetEdgeDetails(phy, includes.intervals=FALSE, intervening.intervals=NULL)
    fossil.taxa <- edge_details$tipward_node[which(edge_details$type == "extinct_tip")]
    pars.bisse <- c(0.1+0.03, 0.1+0.03, 0.03/0.1, 0.03/0.1, 0.01, 0.01)
    
    model.vec <- numeric(48)
    model.vec[1:6] = pars.bisse
    phy$node.label = NULL
    cache <- hisse:::ParametersToPassfHiSSE(model.vec, hidden.states=hidden.states, nb.tip=Ntip(phy), nb.node=Nnode(phy), bad.likelihood=-300, f=c(1,1), ode.eps=0)
    cache$psi <- 0.01
    hisse.full <- hisse:::DownPassHiSSE(dat.tab, gen, cache, root.type="madfitz", condition.on.survival=TRUE, root.p=NULL, node=NULL, state=NULL, fossil.taxa=fossil.taxa, fix.type=NULL)
    
    ## Trait independent model should be loglik_tree + loglik_character ##
    
    #Part 1: MiSSE loglik:
    dat.tab <- hisse:::OrganizeDataMiSSE(phy=phy, f=1, hidden.states=1)
    model.vec <- c(0.1+0.03, 0.03/0.1, rep(0,51))
    cache = hisse:::ParametersToPassMiSSE(model.vec=model.vec, hidden.states=1, fixed.eps=NULL, nb.tip=nb.tip, nb.node=nb.node, bad.likelihood=exp(-300), ode.eps=0)#
    cache$psi <- 0.01
    gen <- hisse:::FindGenerations(phy)
    MiSSE.logL <- hisse:::DownPassMisse(dat.tab=dat.tab, cache=cache, gen=gen, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, fossil.taxa=fossil.taxa, node=NULL, fix.type=NULL)
    
    #Part 2: corHMM loglik:
    library(corHMM)
    char.logL <- corHMM(phy, data, rate.cat=1, model = "ER", node.states = "none", fixed.nodes=FALSE, p=0.01, root.p="maddfitz")
    tot.logL <- char.logL$loglik + MiSSE.logL
    
    comparison <- identical(round(hisse.full,3), round(tot.logL,3))
    expect_true(comparison)
})


test_that("MuHiSSE_fossil_test1", {
    skip_on_cran()

    library(diversitree)
    pars <- c(.1,  .15,  .2, .1, # lambda 1, 2, 3, 4
            .03, .045, .06, 0.03, # mu 1, 2, 3, 4
            .05, .05, .00,        # q12, q13, q14
            .05, .00, .05,     # q21, q23, q24
            .05, .00, .05,     # q31, q32, q34
            .00, .05, .05)
    set.seed(2)
    phy <- NULL
    while( is.null( phy ) ){
        phy <- tree.musse(pars, 30, x0=1, include.extinct=TRUE)
    }
    #f <- hisse:::GetFossils(phy, psi=0.01)
    #h <- history.from.sim.discrete(phy, 1:4)
    #plot(h, phy)
    k.samples <- data.frame(taxon1="sp20", taxon2="sp37", timefrompresent=8.54554, state1=0, state2=1, stringsAsFactors=FALSE)
    
    phy.k <- hisse:::AddKNodes(phy, k.samples)
    fix.type <- hisse:::GetKSampleMRCA(phy.k, k.samples)
    nb.tip <- Ntip(phy.k)
    nb.node <- phy.k$Nnode
    gen <- hisse:::FindGenerations(phy.k)

    states <- phy$tip.state
    states <- data.frame(phy$tip.state, phy$tip.state,
    row.names=names(phy$tip.state))
    states <- states[phy$tip.label,]
    states.trans <- states
    for(i in 1:Ntip(phy)){
        if(states[i,1] == 1){
            states.trans[i,1] = 0
            states.trans[i,2] = 0
        }
        if(states[i,1] == 2){
            states.trans[i,1] = 0
            states.trans[i,2] = 1
        }
        if(states[i,1] == 3){
            states.trans[i,1] = 1
            states.trans[i,2] = 0
        }
        if(states[i,1] == 4){
            states.trans[i,1] = 1
            states.trans[i,2] = 1
        }
    }
    
    data <- data.frame(taxon=names(phy$tip.state), states.trans[,1], states.trans[,2], stringsAsFactors=FALSE)
    data <- hisse:::AddKData(data, k.samples, muhisse=TRUE)
    data.new <- data.frame(data[,2], data[,3], row.names=data[,1])
    data.new <- data.new[phy.k$tip.label,]
    
    pars.muhisse <- c(rep(0.1+0.03,4), rep(0.03/.1, 4), 0.05,0.05,0, 0.05,0,0.05, 0.05,0,.05, 0,0.05,.05)
    model.vec = rep(0,384)
    model.vec[1:20] = pars.muhisse
    cache <- hisse:::ParametersToPassMuHiSSE(model.vec=model.vec, hidden.states=FALSE, nb.tip=Ntip(phy.k), nb.node=Nnode(phy.k), bad.likelihood=exp(-300), f=c(1,1,1,1), ode.eps=0)
    cache$psi <- 0.01
    gen <- hisse:::FindGenerations(phy.k)
    dat.tab <- hisse:::OrganizeData(data.new, phy.k, f=c(1,1,1,1), hidden.states=FALSE, includes.fossils=TRUE)
    fossil.taxa <- which(dat.tab$branch.type == 1)
    
    muhisse.full <- hisse:::DownPassMuHisse(dat.tab, gen=gen, cache=cache, root.type="madfitz", condition.on.survival=TRUE, root.p=NULL, node=fix.type$node, state=fix.type$state, fossil.taxa=fossil.taxa, fix.type=fix.type$type)

    ## Trait independent model should be loglik_tree + loglik_character ##

    #Part 1: MiSSE loglik:
    dat.tab <- hisse:::OrganizeDataMiSSE(phy=phy.k, f=1, hidden.states=1)
    model.vec <- c(0.1+0.03, 0.03/0.1, rep(0,51))
    cache = hisse:::ParametersToPassMiSSE(model.vec=model.vec, hidden.states=1, fixed.eps=NULL, nb.tip=nb.tip, nb.node=nb.node, bad.likelihood=exp(-300), ode.eps=0)#
    cache$psi <- 0.01
    edge_details <- hisse:::GetEdgeDetails(phy, intervening.intervals=strat.cache$intervening.intervals)
    fossil.taxa <- edge_details$tipward_node[which(edge_details$type == "extinct_tip")]
    gen <- hisse:::FindGenerations(phy.k)
    MiSSE.logL <- hisse:::DownPassMisse(dat.tab=dat.tab, cache=cache, gen=gen, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, fossil.taxa=fossil.taxa, node=fix.type$node, fix.type=fix.type$type)
    
    #Part 2: corHMM loglik:
    library(corHMM)
    char.logL <- corHMM(phy.k, data, rate.cat=1, model = "ER", node.states = "none", fixed.nodes=FALSE, p=0.05, root.p="maddfitz")
    tot.logL <- char.logL$loglik + MiSSE.logL
    
    comparison1 <- identical(round(muhisse.full,3), round(tot.logL,3))
    expect_true(comparison1)
    
    #Same test, but fossils removed:
    extinct.tip.no <- grep("ex*", x=phy$tip.label)
    k.samp.tip.no <- grep("Ksamp*", x=data[,1])
    to.delete <- c(extinct.tip.no, k.samp.tip.no)
    phy.extant <- drop.tip(phy, phy$tip.label[extinct.tip.no])
    data.tmp <- data
    data.tmp <- data.tmp[-c(to.delete),]
    data.new <- data.frame(data.tmp[,2], data.tmp[,3], row.names=data.tmp[,1])
    data.new <- data.new[phy.extant$tip.label,]
    
    gen <- hisse:::FindGenerations(phy.extant)
    
    dat.tab <- hisse:::OrganizeData(data.new, phy=phy.extant, f=c(1,1,1,1), hidden.states=FALSE)
    pars.muhisse <- c(rep(0.1+0.03,4), rep(0.03/.1, 4), 0.05,0.05,0, 0.05,0,0.05, 0.05,0,.05, 0,0.05,.05)
    model.vec = rep(0,384)
    model.vec[1:20] = pars.muhisse
    cache <- hisse:::ParametersToPassMuHiSSE(model.vec, hidden.states=FALSE, nb.tip=Ntip(phy.extant), nb.node=Nnode(phy.extant), bad.likelihood=-300, f=c(1,1), ode.eps=0)
    cache$psi <- 0
    muhisse.full <- hisse:::DownPassMuHisse(dat.tab, gen, cache, root.type="madfitz", condition.on.survival=TRUE, root.p=NULL, fossil.taxa=NULL)
    
    dat.tab <- hisse:::OrganizeDataMiSSE(phy=phy.extant, f=1, hidden.states=1)
    model.vec <- c(0.1+0.03, 0.03/0.1, rep(0,51))
    cache = hisse:::ParametersToPassMiSSE(model.vec=model.vec, hidden.states=1, fixed.eps=NULL, nb.tip=Ntip(phy.extant), nb.node=Nnode(phy.extant), bad.likelihood=exp(-300), ode.eps=0)#
    cache$psi <- 0.0
    edge_details <- hisse:::GetEdgeDetails(phy, intervening.intervals=strat.cache$intervening.intervals)
    fossil.taxa <- edge_details$tipward_node[which(edge_details$type == "extinct_tip")]
    gen <- hisse:::FindGenerations(phy.extant)
    MiSSE.logL <- hisse:::DownPassMisse(dat.tab=dat.tab, cache=cache, gen=gen, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, fossil.taxa=NULL, node=NULL, fix.type=NULL)
    
    library(corHMM)
    char.logL <- corHMM(phy.extant, data, rate.cat=1, model = "ER", node.states = "none", fixed.nodes=FALSE, p=0.05, root.p="maddfitz")
    tot.logL <- char.logL$loglik + MiSSE.logL
    
    comparison2 <- identical(round(muhisse.full,3), round(tot.logL,3))
    expect_true(comparison2)
})


test_that("MuHiSSE_fossil_test2", {
    skip_on_cran()
    
    library(diversitree)
    pars <- c(.1,  .15,  .2, .1, # lambda 1, 2, 3, 4
    .03, .045, .06, 0.03, # mu 1, 2, 3, 4
    .05, .05, .00,        # q12, q13, q14
    .05, .00, .05,     # q21, q23, q24
    .05, .00, .05,     # q31, q32, q34
    .00, .05, .05)
    set.seed(2)
    phy <- NULL
    while( is.null( phy ) ){
        phy <- tree.musse(pars, 30, x0=1, include.extinct=TRUE)
    }
    
    fix.type <- NULL
    nb.tip <- Ntip(phy)
    nb.node <- phy$Nnode
    gen <- hisse:::FindGenerations(phy)
    
    states <- phy$tip.state
    states <- data.frame(phy$tip.state, phy$tip.state,
    row.names=names(phy$tip.state))
    states <- states[phy$tip.label,]
    states.trans <- states
    for(i in 1:Ntip(phy)){
        if(states[i,1] == 1){
            states.trans[i,1] = 0
            states.trans[i,2] = 0
        }
        if(states[i,1] == 2){
            states.trans[i,1] = 0
            states.trans[i,2] = 1
        }
        if(states[i,1] == 3){
            states.trans[i,1] = 1
            states.trans[i,2] = 0
        }
        if(states[i,1] == 4){
            states.trans[i,1] = 1
            states.trans[i,2] = 1
        }
    }
    
    data <- data.frame(taxon=names(phy$tip.state), states.trans[,1], states.trans[,2], stringsAsFactors=FALSE)
    data.new <- data.frame(data[,2], data[,3], row.names=data[,1])
    data.new <- data.new[phy$tip.label,]
    
    pars.muhisse <- c(rep(0.1+0.03,4), rep(0.03/.1, 4), 0.05,0.05,0, 0.05,0,0.05, 0.05,0,.05, 0,0.05,.05)
    model.vec = rep(0,384)
    model.vec[1:20] = pars.muhisse
    cache <- hisse:::ParametersToPassMuHiSSE(model.vec=model.vec, hidden.states=FALSE, nb.tip=Ntip(phy), nb.node=Nnode(phy), bad.likelihood=exp(-300), f=c(1,1,1,1), ode.eps=0)
    cache$psi <- 0.01
    gen <- hisse:::FindGenerations(phy)
    dat.tab <- hisse:::OrganizeData(data.new, phy, f=c(1,1,1,1), hidden.states=FALSE, includes.fossils=TRUE)
    fossil.taxa <- which(dat.tab$branch.type == 1)
    
    muhisse.full <- hisse:::DownPassMuHisse(dat.tab, gen=gen, cache=cache, root.type="madfitz", condition.on.survival=TRUE, root.p=NULL, node=fix.type$node, state=fix.type$state, fossil.taxa=fossil.taxa, fix.type=fix.type$type)
    
    ## Trait independent model should be loglik_tree + loglik_character ##
    
    #Part 1: MiSSE loglik:
    dat.tab <- hisse:::OrganizeDataMiSSE(phy=phy, f=1, hidden.states=1)
    model.vec <- c(0.1+0.03, 0.03/0.1, rep(0,51))
    cache = hisse:::ParametersToPassMiSSE(model.vec=model.vec, hidden.states=1, fixed.eps=NULL, nb.tip=nb.tip, nb.node=nb.node, bad.likelihood=exp(-300), ode.eps=0)#
    cache$psi <- 0.01
    edge_details <- GetEdgeDetails(phy, intervening.intervals=strat.cache$intervening.intervals)
    fossil.taxa <- edge_details$tipward_node[which(edge_details$type == "extinct_tip")]
    gen <- hisse:::FindGenerations(phy)
    MiSSE.logL <- hisse:::DownPassMisse(dat.tab=dat.tab, cache=cache, gen=gen, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, fossil.taxa=fossil.taxa, node=fix.type$node, fix.type=fix.type$type)
    
    #Part 2: corHMM loglik:
    library(corHMM)
    char.logL <- corHMM(phy, data, rate.cat=1, model = "ER", node.states = "none", fixed.nodes=FALSE, p=0.05, root.p="maddfitz")
    tot.logL <- char.logL$loglik + MiSSE.logL
    
    comparison <- identical(round(muhisse.full,3), round(tot.logL,3))
    expect_true(comparison)

})


test_that("MiSSE_fossil_test4", {
    skip_on_cran()
    
    library(diversitree)
    pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
    set.seed(4)
    phy <- NULL
    while( is.null( phy ) ){
        phy <- tree.bisse(pars, max.t=30, x0=0, include.extinct=TRUE)
    }
    #h <- history.from.sim.discrete(phy, 0:1)
    #plot(h, phy)
    k.samples <- data.frame(taxon1="sp12", taxon2="sp12", timefrompresent=3.164384, state=1, stringsAsFactors=FALSE)
    
    hidden.states=TRUE
    
    phy.k <- hisse:::AddKNodes(phy, k.samples)
    fix.type <- hisse:::GetKSampleMRCA(phy.k, k.samples)
    nb.tip <- Ntip(phy.k)
    nb.node <- phy.k$Nnode
    gen <- hisse:::FindGenerations(phy.k)
    
    data <- data.frame(taxon=names(phy$tip.state), phy$tip.state, stringsAsFactors=FALSE)
    data <- hisse:::AddKData(data, k.samples)
    data.new <- data.frame(data[,2], data[,2], row.names=data[,1])
    data.new <- data.new[phy.k$tip.label,]
    
    dat.tab <- hisse:::OrganizeDataHiSSE(data.new, phy=phy.k, f=c(1,1), hidden.states=TRUE)
    edge_details <- hisse:::GetEdgeDetails(phy=phy.k, includes.intervals=FALSE, intervening.intervals=NULL)
    fossil.taxa <- edge_details$tipward_node[which(edge_details$type == "extinct_tip")]
    pars.hisse <- c(0.1+0.03, 0.1+0.03, 0.03/0.1, 0.03/0.1, 0.01, 0.01, 0.01, rep(0,5), 0.2+0.03, 0.2+0.03, 0.03/0.2, 0.03/0.2, 0.01, 0.01, 0.01, rep(0,5))
    pars.hisse[10] <- 0.01
    pars.hisse[22] <- 0.01
    model.vec <- numeric(48)
    model.vec[1:24] = pars.hisse
    
    phy$node.label = NULL
    cache <- hisse:::ParametersToPassfHiSSE(model.vec, hidden.states=hidden.states, nb.tip=Ntip(phy.k), nb.node=Nnode(phy.k), bad.likelihood=-300, f=c(1,1), ode.eps=0)
    cache$psi <- 0.02
    hisse.full <- hisse:::DownPassHiSSE(dat.tab, gen, cache, root.type="madfitz", condition.on.survival=TRUE, root.p=NULL, node=fix.type$node, state=fix.type$state, fossil.taxa=fossil.taxa, fix.type=fix.type$type)
    
    ## Trait independent model should be loglik_tree + loglik_character ##
    
    #Part 1: MiSSE loglik:
    dat.tab <- hisse:::OrganizeDataMiSSE(phy=phy.k, f=1, hidden.states=2)
    model.vec <- c(0.1+0.03, 0.03/0.1, 0.2+0.03, 0.03/0.2, rep(0,48), 0.01, 0.02)
    cache = hisse:::ParametersToPassMiSSE(model.vec=model.vec, hidden.states=2, fixed.eps=NULL, nb.tip=nb.tip, nb.node=nb.node, bad.likelihood=exp(-300), ode.eps=0)#
    gen <- hisse:::FindGenerations(phy.k)
    MiSSE.logL <- hisse:::DownPassMisse(dat.tab=dat.tab, cache=cache, gen=gen, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, fossil.taxa=fossil.taxa, node=fix.type$node, fix.type=fix.type$type)
    
    #Part 2: corHMM loglik:
    library(corHMM)
    char.logL <- corHMM(phy.k, data, rate.cat=1, model = "ER", node.states = "none", fixed.nodes=FALSE, p=0.01, root.p="maddfitz")
    tot.logL <- char.logL$loglik + MiSSE.logL
    
    comparison1 <- identical(round(hisse.full,3), round(tot.logL,3))
    expect_true(comparison1)
    
})


test_that("interval_branch_test1", {
    skip_on_cran()

    compE <- numeric(26)
    compD <- numeric(26)
    compD[1] <- 1
    
    yini <- c(E0A = compE[1], E0B = compE[2], E0C = compE[3], E0D = compE[4], E0E = compE[5], E0F = compE[6], E0G = compE[7], E0H = compE[8], E0I = compE[9], E0J = compE[10], E0K = compE[11], E0L = compE[12], E0M = compE[13], E0N = compE[14], E0O = compE[15], E0P = compE[16], E0Q = compE[17], E0R = compE[18], E0S = compE[19], E0T = compE[20], E0U = compE[21], E0V = compE[22], E0W = compE[23], E0X = compE[24], E0Y = compE[25], E0Z = compE[26], D0A = compD[1], D0B = compD[2], D0C = compD[3], D0D = compD[4], D0E = compD[5], D0F = compD[6], D0G = compD[7], D0H = compD[8], D0I = compD[9], D0J = compD[10], D0K = compD[11], D0L = compD[12], D0M = compD[13], D0N = compD[14], D0O = compD[15], D0P = compD[16], D0Q = compD[17], D0R = compD[18], D0S = compD[19], D0T = compD[20], D0U = compD[21], D0V = compD[22], D0W = compD[23], D0X = compD[24], D0Y = compD[25], D0Z = compD[26])
    
    times=c(0, 10)
    
    pars <- numeric(55)
    pars[55] <- 0.05
    pars[1] <- .2
    pars[2] <- .1
    pars[54] <- 1
        
    stad.1 <- hisse:::q_sym_ratio(si=10, ei=0, lambda=.2, mu=.1, psi=0.05, rho=1)
    
    ## Now do same calculation using SingleChildCode
    rho = 1
    turnover = .2 + .1
    eps = .1 / .2
    psi = 0.05
    model.vec <- c(turnover, eps, rep(0,51))
    cache <- hisse:::ParametersToPassMiSSE(model.vec=model.vec, hidden.states=1, fixed.eps=NULL, nb.tip=1, nb.node=1, bad.likelihood=exp(-300), ode.eps=0)#
    cache$psi <- psi
    strat.range.calc <- hisse:::SingleChildProbMiSSE(cache, pars, compD, compE,  0, 10, 3)
    comparison <- identical(round(unname(strat.range.calc[27]),4), round(stad.1,4))
    
    expect_true(comparison)
})


test_that("interval_branch_test2", {
    skip_on_cran()

    compE <- numeric(26)
    compD <- numeric(26)
    compD[1] <- 1
    
    yini <- c(E0A = compE[1], E0B = compE[2], E0C = compE[3], E0D = compE[4], E0E = compE[5], E0F = compE[6], E0G = compE[7], E0H = compE[8], E0I = compE[9], E0J = compE[10], E0K = compE[11], E0L = compE[12], E0M = compE[13], E0N = compE[14], E0O = compE[15], E0P = compE[16], E0Q = compE[17], E0R = compE[18], E0S = compE[19], E0T = compE[20], E0U = compE[21], E0V = compE[22], E0W = compE[23], E0X = compE[24], E0Y = compE[25], E0Z = compE[26], D0A = compD[1], D0B = compD[2], D0C = compD[3], D0D = compD[4], D0E = compD[5], D0F = compD[6], D0G = compD[7], D0H = compD[8], D0I = compD[9], D0J = compD[10], D0K = compD[11], D0L = compD[12], D0M = compD[13], D0N = compD[14], D0O = compD[15], D0P = compD[16], D0Q = compD[17], D0R = compD[18], D0S = compD[19], D0T = compD[20], D0U = compD[21], D0V = compD[22], D0W = compD[23], D0X = compD[24], D0Y = compD[25], D0Z = compD[26])
    
    times=c(0, 10)
    
    pars <- numeric(55)
    pars[55] <- 0.05
    pars[1] <- .2
    pars[2] <- .1
    pars[54] <- 1
    
    stad.1 <- hisse:::q_sym_ratio(si=10, ei=0, lambda=.2, mu=.1, psi=0.05, rho=1)
    
    ## Now do same calculation using SingleChildCode
    rho = 1
    turnover = .2 + .1
    eps = .1 / .2
    psi = 0.05
    model.vec <- c(turnover, eps, rep(0,51))
    cache <- hisse:::ParametersToPassMiSSE(model.vec=model.vec, hidden.states=1, fixed.eps=NULL, nb.tip=1, nb.node=1, bad.likelihood=exp(-300), ode.eps=0)#
    cache$psi <- psi
    strat.range.calc <- hisse:::SingleChildProbMiSSE(cache, pars, compD, compE,  0, 10, 3)
    comparison <- identical(round(unname(strat.range.calc[27]),4), round(stad.1,4))
    
    ## Now test an intervening branch calculations:
    first.ratio <- hisse:::q_t(t=10, lambda=.2, mu=.1, psi=0.05, rho=1) / hisse:::q_sym_tilde(t=10,lambda=.2, mu=.1, psi=0.05, rho=1)
    second.ratio <- hisse:::q_sym_tilde(t=12, lambda=.2, mu=.1, psi=0.05, rho=1)/ hisse:::q_t(t=12,lambda=.2, mu=.1, psi=0.05, rho=1)
    stad.unobs <- 1 - (first.ratio * second.ratio)
    stad.2 <- stad.unobs * hisse:::q_ratio(12,10,lambda=.2, mu=.1, psi=0.05, rho=1) * stad.1
    
    strat.intervene.calc <- hisse:::SingleChildProbMiSSE(cache, pars, strat.range.calc[27:52], strat.range.calc[1:26], 10, 12, 4)
    comparison <- identical(round(unname(strat.intervene.calc[27]),4), round(stad.2,4))
    
    expect_true(comparison)
})


test_that("MiSSE_interval_test1", {
    skip_on_cran()
    
    s <- "(A:20,(ex_1:4,C:12):8):0;"
    cat(s, file = "ex.tre", sep = "\n")
    phy <- read.tree(s, file="ex.tre", sep="\n")
    strat.intervals <- data.frame(taxon1=c("A", "ex_1"), taxon2=c("A", "C"), timefrompresentroot=c(10, 17), timefrompresenttip=c(0, 15), type=c("R", "R"), stringsAsFactors=FALSE)
    
    rho = 1
    turnover = .2 + .1
    eps = .1 / .2
    psi = 0.05
    
    phy.og <- phy
    
    split.times.all <- paleotree::dateNodes(phy, rootAge=max(node.depth.edgelength(phy)))
    split.times <- split.times.all[-c(1:Ntip(phy))]
    strat.cache <- hisse:::GetStratInfo(strat.intervals)
    k.samples <- hisse:::GetIntervalToK(strat.intervals, intervening.intervals=strat.cache$intervening.intervals)
    phy <- hisse:::AddKNodes(phy, k.samples)
    fix.type <- hisse:::GetKSampleMRCA(phy, k.samples, strat.intervals=TRUE)
    edge_details <- hisse:::GetEdgeDetails(phy, includes.intervals=TRUE, intervening.intervals=strat.cache$intervening.intervals)
    
    gen <- hisse:::FindGenerations(phy)
    dat.tab <- hisse:::OrganizeDataMiSSE(phy=phy, f=1, hidden.states=1, includes.intervals=TRUE, intervening.intervals=strat.cache$intervening.intervals)
    
    #These are all inputs for generating starting values:
    fossil.taxa <- edge_details$tipward_node[which(edge_details$type == "extinct_tip" | edge_details$type == "k_extinct_interval")]
    fossil.ages <- dat.tab$TipwardAge[which(dat.tab$DesNode %in% fossil.taxa)]
    
    cols <- c("FocalNode","DesNode", "RootwardAge", "TipwardAge", "branch.type")
    seg.map <- dat.tab[, cols, with=FALSE]
    #remove k tips -- we do not do anything with them.
    data.table::setkey(seg.map, branch.type)
    seg.map <- seg.map[branch.type != 2]
    
    logLikLogSpace <- -hisse:::starting.point.tree.intervals(x=log(c(turnover, eps, psi)), n.tax=3, rho=1, seg_map=seg.map, x_times=split.times, y_times=fossil.ages, strat.cache=strat.cache)
    
    nb.tip <- Ntip(phy)
    nb.node <- phy$Nnode
    model.vec <- c(turnover, eps, rep(0,51))
    cache <- hisse:::ParametersToPassMiSSE(model.vec=model.vec, hidden.states=1, fixed.eps=NULL, nb.tip=nb.tip, nb.node=nb.node, bad.likelihood=exp(-300), ode.eps=0)#
    cache$psi <- psi
    
    MiSSE.logL <- hisse:::DownPassMisse(dat.tab=dat.tab, cache=cache, gen=gen, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, fossil.taxa=fossil.taxa, node=fix.type$node, fix.type=fix.type$type) + (strat.cache$k*log(psi)) + (psi*strat.cache$l_s)
    
    comparison <- identical(round(logLikLogSpace,4), round(MiSSE.logL,4))
    
    expect_true(comparison)
})


test_that("MiSSE_interval_test2", {
    skip_on_cran()
    
    s <- "(A:20,(ex_1:4,C:12):8):0;"
    cat(s, file = "ex.tre", sep = "\n")
    phy <- read.tree(s, file="ex.tre", sep="\n")
    strat.intervals <- data.frame(taxon1=c("A", "A", "ex_1"), taxon2=c("A", "A", "C"), timefrompresentroot=c(10, 14, 17), timefrompresenttip=c(0, 14, 15), type=c("R", "S","R"), stringsAsFactors=FALSE)
    
    rho = 1
    turnover = .2 + .1
    eps = .1 / .2
    psi = 0.05
        
    split.times.all <- paleotree::dateNodes(phy, rootAge=max(node.depth.edgelength(phy)))
    split.times <- split.times.all[-c(1:Ntip(phy))]
    strat.cache <- hisse:::GetStratInfo(strat.intervals)
    k.samples <- hisse:::GetIntervalToK(strat.intervals, intervening.intervals=strat.cache$intervening.intervals)
    phy <- hisse:::AddKNodes(phy, k.samples)
    fix.type <- hisse:::GetKSampleMRCA(phy, k.samples, strat.intervals=TRUE)
    edge_details <- hisse:::GetEdgeDetails(phy, includes.intervals=TRUE, intervening.intervals=strat.cache$intervening.intervals)
    
    gen <- hisse:::FindGenerations(phy)
    dat.tab <- hisse:::OrganizeDataMiSSE(phy=phy, f=1, hidden.states=1, includes.intervals=TRUE, intervening.intervals=strat.cache$intervening.intervals)
    
    #These are all inputs for generating starting values:
    fossil.taxa <- edge_details$tipward_node[which(edge_details$type == "extinct_tip" | edge_details$type == "k_extinct_interval")]
    fossil.ages <- dat.tab$TipwardAge[which(dat.tab$DesNode %in% fossil.taxa)]
    
    cols <- c("FocalNode","DesNode", "RootwardAge", "TipwardAge", "branch.type")
    seg.map <- dat.tab[, cols, with=FALSE]
    #remove k tips -- we do not do anything with them.
    data.table::setkey(seg.map, branch.type)
    seg.map <- seg.map[branch.type != 2]
    
    logLikLogSpace <- -hisse:::starting.point.tree.intervals(x=log(c(turnover, eps, psi)), n.tax=3, rho=1, seg_map=seg.map, x_times=split.times, y_times=fossil.ages, strat.cache=strat.cache)
    
    nb.tip <- Ntip(phy)
    nb.node <- phy$Nnode
    model.vec <- c(turnover, eps, rep(0,51))
    cache <- hisse:::ParametersToPassMiSSE(model.vec=model.vec, hidden.states=1, fixed.eps=NULL, nb.tip=nb.tip, nb.node=nb.node, bad.likelihood=exp(-300), ode.eps=0)#
    cache$psi <- psi
    
    MiSSE.logL <- hisse:::DownPassMisse(dat.tab=dat.tab, cache=cache, gen=gen, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, fossil.taxa=fossil.taxa, node=fix.type$node, fix.type=fix.type$type) + (strat.cache$k*log(psi)) + (psi*strat.cache$l_s)
    
    comparison <- identical(round(logLikLogSpace,4), round(MiSSE.logL,4))

    expect_true(comparison)
})


test_that("MiSSE_interval_test3", {
    skip_on_cran()
    
    s <- "(A:20,(ex_1:4,C:12):8):0;"
    cat(s, file = "ex.tre", sep = "\n")
    phy <- read.tree(s, file="ex.tre", sep="\n")
    strat.intervals <- data.frame(taxon1=c("A", "A", "ex_1"), taxon2=c("A", "A", "C"), timefrompresentroot=c(10, 14, 17), timefrompresenttip=c(0, 12, 15), type=c("R", "R", "R"), stringsAsFactors=FALSE)
    
    rho = 1
    turnover = .2 + .1
    eps = .1 / .2
    psi = 0.05
    
    phy.og <- phy
    
    split.times.all <- paleotree::dateNodes(phy, rootAge=max(node.depth.edgelength(phy)))
    split.times <- split.times.all[-c(1:Ntip(phy))]
    strat.cache <- hisse:::GetStratInfo(strat.intervals)
    k.samples <- hisse:::GetIntervalToK(strat.intervals, intervening.intervals=strat.cache$intervening.intervals)
    phy <- hisse:::AddKNodes(phy, k.samples)
    fix.type <- hisse:::GetKSampleMRCA(phy, k.samples, strat.intervals=TRUE)
    edge_details <- hisse:::GetEdgeDetails(phy, includes.intervals=TRUE, intervening.intervals=strat.cache$intervening.intervals)
    
    gen <- hisse:::FindGenerations(phy)
    dat.tab <- hisse:::OrganizeDataMiSSE(phy=phy, f=1, hidden.states=1, includes.intervals=TRUE, intervening.intervals=strat.cache$intervening.intervals)
    
    #These are all inputs for generating starting values:
    fossil.taxa <- edge_details$tipward_node[which(edge_details$type == "extinct_tip" | edge_details$type == "k_extinct_interval")]
    fossil.ages <- dat.tab$TipwardAge[which(dat.tab$DesNode %in% fossil.taxa)]
    
    cols <- c("FocalNode","DesNode", "RootwardAge", "TipwardAge", "branch.type")
    seg.map <- dat.tab[, cols, with=FALSE]
    #remove k tips -- we do not do anything with them.
    data.table::setkey(seg.map, branch.type)
    seg.map <- seg.map[branch.type != 2]
    
    logLikLogSpace <- -hisse:::starting.point.tree.intervals(x=log(c(turnover, eps, psi)), n.tax=3, rho=1, seg_map=seg.map, x_times=split.times, y_times=fossil.ages, strat.cache=strat.cache)
    
    nb.tip <- Ntip(phy)
    nb.node <- phy$Nnode
    model.vec <- c(turnover, eps, rep(0,51))
    cache <- hisse:::ParametersToPassMiSSE(model.vec=model.vec, hidden.states=1, fixed.eps=NULL, nb.tip=nb.tip, nb.node=nb.node, bad.likelihood=exp(-300), ode.eps=0)#
    cache$psi <- psi
    
    MiSSE.logL <- hisse:::DownPassMisse(dat.tab=dat.tab, cache=cache, gen=gen, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, fossil.taxa=fossil.taxa, node=fix.type$node, fix.type=fix.type$type) + (strat.cache$k*log(psi)) + (psi*strat.cache$l_s)
    
    comparison <- identical(round(logLikLogSpace,4), round(MiSSE.logL,4))
    
    expect_true(comparison)
})


test_that("MiSSE_interval_test4", {
    skip_on_cran()
    
    s <- "(A:20,(ex_1:4,C:12):8):0;"
    cat(s, file = "ex.tre", sep = "\n")
    phy <- read.tree(s, file="ex.tre", sep="\n")
    strat.intervals <- data.frame(taxon1=c("A", "ex_1", "A", "ex_1"), taxon2=c("A", "ex_1", "A", "C"), timefrompresentroot=c(10, 11, 14, 17), timefrompresenttip=c(0, 8, 12, 15), type=c("R", "R", "R", "R"), stringsAsFactors=FALSE)
    
    rho = 1
    turnover = .2 + .1
    eps = .1 / .2
    psi = 0.05
    
    phy.og <- phy
    
    split.times.all <- paleotree::dateNodes(phy, rootAge=max(node.depth.edgelength(phy)))
    split.times <- split.times.all[-c(1:Ntip(phy))]
    strat.cache <- hisse:::GetStratInfo(strat.intervals)
    k.samples <- hisse:::GetIntervalToK(strat.intervals, intervening.intervals=strat.cache$intervening.intervals)
    k.samples <- k.samples[-1,]
    phy <- hisse:::AddKNodes(phy, k.samples)
    fix.type <- hisse:::GetKSampleMRCA(phy, k.samples, strat.intervals=TRUE)
    edge_details <- hisse:::GetEdgeDetails(phy, includes.intervals=TRUE, intervening.intervals=strat.cache$intervening.intervals)
    
    gen <- hisse:::FindGenerations(phy)
    dat.tab <- hisse:::OrganizeDataMiSSE(phy=phy, f=1, hidden.states=1, includes.intervals=TRUE, intervening.intervals=strat.cache$intervening.intervals)
    
    #These are all inputs for generating starting values:
    fossil.taxa <- edge_details$tipward_node[which(edge_details$type == "extinct_tip" | edge_details$type == "k_extinct_interval")]
    fossil.ages <- dat.tab$TipwardAge[which(dat.tab$DesNode %in% fossil.taxa)]
    
    cols <- c("FocalNode","DesNode", "RootwardAge", "TipwardAge", "branch.type")
    seg.map <- dat.tab[, cols, with=FALSE]
    #remove k tips -- we do not do anything with them.
    data.table::setkey(seg.map, branch.type)
    seg.map <- seg.map[branch.type != 2]
    
    logLikLogSpace <- -hisse:::starting.point.tree.intervals(x=log(c(turnover, eps, psi)), n.tax=3, rho=1, seg_map=seg.map, x_times=split.times, y_times=fossil.ages, strat.cache=strat.cache)
    
    nb.tip <- Ntip(phy)
    nb.node <- phy$Nnode
    model.vec <- c(turnover, eps, rep(0,51))
    cache <- hisse:::ParametersToPassMiSSE(model.vec=model.vec, hidden.states=1, fixed.eps=NULL, nb.tip=nb.tip, nb.node=nb.node, bad.likelihood=exp(-300), ode.eps=0)#
    cache$psi <- psi
    
    MiSSE.logL <- hisse:::DownPassMisse(dat.tab=dat.tab, cache=cache, gen=gen, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, fossil.taxa=fossil.taxa, node=fix.type$node, fix.type=fix.type$type) + (strat.cache$k*log(psi)) + (psi*strat.cache$l_s)
    
    comparison <- identical(round(logLikLogSpace,4), round(MiSSE.logL,4))
    
    expect_true(comparison)
})


test_that("MiSSE_interval_test5", {
    skip_on_cran()
    
    ntax=200
    true.psi=0.1
    set.seed(4)
    try(sim.tab <- hisse::SimulateHisse(turnover=c(0.25,0.25), eps=rep(0.75,2), max.taxa=ntax, x0=0, transition.rates=matrix(c(NA, 0.005, 0.005, NA), nrow=2), nstart=2))
    phy <- SimToPhylo(sim.tab, include.extinct=TRUE)
    f <- hisse:::GetFossils(phy, psi=true.psi)
    
    ss <- hisse:::ProcessSimStrat(phy, f)
    split.times.all <- paleotree::dateNodes(ss$phy, rootAge=max(node.depth.edgelength(ss$phy)))
    split.times <- split.times.all[-c(1:Ntip(ss$phy))]
    strat.cache <- hisse:::GetStratInfo(ss$strat.intervals)
    k.samples <- hisse:::GetIntervalToK(ss$strat.intervals, intervening.intervals=strat.cache$intervening.intervals)
    k.samples <- k.samples[-which(round(k.samples$timefrompresent,8) %in% round(split.times.all[c(1:Ntip(ss$phy))],8)),]
    phy <- hisse:::AddKNodes(ss$phy, k.samples)
    fix.type <- hisse:::GetKSampleMRCA(phy, k.samples, strat.intervals=TRUE)
    edge_details <- hisse:::GetEdgeDetails(phy, includes.intervals=TRUE, intervening.intervals=strat.cache$intervening.intervals)
    
    gen <- hisse:::FindGenerations(phy)
    dat.tab <- hisse:::OrganizeDataMiSSE(phy=phy, f=1, hidden.states=1, includes.intervals=TRUE, intervening.intervals=strat.cache$intervening.intervals)
    fossil.taxa <- edge_details$tipward_node[which(edge_details$type == "extinct_tip" | edge_details$type == "k_extinct_interval")]
    fossil.ages <- dat.tab$TipwardAge[which(dat.tab$DesNode %in% fossil.taxa)]
    
    cols <- c("FocalNode","DesNode", "RootwardAge", "TipwardAge", "branch.type")
    seg.map <- dat.tab[, cols, with=FALSE]
    #remove k tips -- we do not do anything with them.
    data.table::setkey(seg.map, branch.type)
    #drop the k.tips because we do not do calculation on these zero length edges:
    seg.map <- seg.map[branch.type != 2]
    #######################
    
    start.points <- hisse:::starting.point.generator.intervals(k=1, n=ntax, samp.freq=1, seg_map=seg.map, split.times=split.times, fossil.ages=fossil.ages, strat.cache=strat.cache, get.likelihood=TRUE)
    
    turnover <- start.points[1] + start.points[2]
    eps <- start.points[2]/start.points[1]
    psi <- start.points[3]
    
    nb.tip <- Ntip(phy)
    nb.node <- phy$Nnode
    model.vec <- c(turnover, eps, rep(0,51))
    cache <- hisse:::ParametersToPassMiSSE(model.vec=model.vec, hidden.states=1, fixed.eps=NULL, nb.tip=nb.tip, nb.node=nb.node, bad.likelihood=exp(-300), ode.eps=0)#
    cache$psi <- psi
    
    MiSSE.logL <- hisse:::DownPassMisse(dat.tab=dat.tab, cache=cache, gen=gen, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, fossil.taxa=fossil.taxa, node=fix.type$node, fix.type=fix.type$type) + (strat.cache$k*log(psi)) + (psi*strat.cache$l_s)
    
    comparison <- identical(round(-start.points[4],3), round(MiSSE.logL,3))

    expect_true(comparison)
})

