FormatOutHiGeoSSE <- function(fit.geosse){
    ## Function to repackage the output from HiGeoSSE.
    ## This seems to be more readable.

    ## Create a list to store everything.
    out <- list()
    ## List the model fit information.
    out$fit <- list()
    out$fit$loglik <- fit.geosse$loglik
    out$fit$AIC <- fit.geosse$AIC
    out$fit$AICc <- fit.geosse$AICc
    ## List with the imput data for the analysis.
    imput <- list()
    out$imput$phy <- fit.geosse$phy
    out$imput$data <- fit.geosse$data
    ## List with the control for the analysis.
    control <- list()
    out$control$index.par <- fit.geosse$index.par
    out$control$f <- fit.geosse$f
    out$control$hidden.areas <- fit.geosse$hidden.areas
    out$control$condition.on.survival <- fit.geosse$condition.on.survival
    out$control$root.type <- fit.geosse$root.type
    out$control$root.p <- fit.geosse$root.p
    out$control$trans.matrix <- fit.geosse$trans.matrix
    out$control$max.tol <- fit.geosse$max.tol
    out$control$starting.vals <- fit.geosse$starting.vals
    out$control$upper.bounds <- fit.geosse$upper.bounds
    out$control$lower.bounds <- fit.geosse$lower.bounds
    out$control$ode.eps <- fit.geosse$ode.eps

    ## Create edited format to report parameters back:
    
    ## The solution vector is hard coded because of the C code.
    ## The position of the names will always be the same.
    parnames <- c("s0A", "s1A", "s01A", "x0A", "x1A", "d0A", "d1A"
                , "q0A_0B", "q0A_0C", "q0A_0D", "q0A_0E"
                , "q1A_1B", "q1A_1C", "q1A_1D", "q1A_1E"
                , "q01A_01B", "q01A_01C", "q01A_01D", "q01A_01E"
                , "s0B", "s1B", "s01B", "x0B", "x1B", "d0B", "d1B"
                , "q0B_0A", "q0B_0C", "q0B_0D", "q0B_0E"
                , "q1B_1A", "q1B_1C", "q1B_1D", "q1B_1E"
                , "q01B_01A", "q01B_01C", "q01B_01D", "q01B_01E"
                , "s0C", "s1C", "s01C", "x0C", "x1C", "d0C", "d1C"
                , "q0C_0A", "q0C_0B", "q0C_0D", "q0C_0E"
                , "q1C_1A", "q1C_1B", "q1C_1D", "q1C_1E"
                , "q01C_01A", "q01C_01B", "q01C_01D", "q01C_01E"
                , "s0D", "s1D", "s01D", "x0D", "x1D", "d0D", "d1D"
                , "q0D_0A", "q0D_0B", "q0D_0C", "q0D_0E"
                , "q1D_1A", "q1D_1B", "q1D_1C", "q1D_1E"
                , "q01D_01A", "q01D_01B", "q01D_01C", "q01D_01E"
                , "s0E", "s1E", "s01E", "x0E", "x1E", "d0E", "d1E"
                , "q0E_0A", "q0E_0B", "q0E_0C", "q0E_0D"
                , "q1E_1A", "q1E_1B", "q1E_1C", "q1E_1ED"
                , "q01E_01A", "q01E_01B", "q01E_01C", "q01E_01D")
    ## Add names to the 'index par' control. This will help A LOT!
    names( out$control$index.par ) <- parnames
    solution <- fit.geosse$solution
    names( solution ) <- parnames
    ## solution[out$control$index.par == max(out$control$index.par)] <- NA
    ## First group of parameters:
    geosse.pars <- matrix(nrow=7, ncol=5)
    colnames(geosse.pars) <- LETTERS[1:5]
    rownames(geosse.pars) <- c("s0","s1","s01","x0","x1","d0","d1")
    for(i in 1:5) geosse.pars[,(1*i)] <- solution[(1+(19*(i-1))):(7+(19*(i-1)))]
    ## Second group of parameters:
    q0.pars <- matrix(nrow=5, ncol=5)
    colnames(q0.pars) <- rownames(q0.pars) <- paste0("0", LETTERS[1:5])
    for(i in 1:5) q0.pars[(1*i),(1:5)[-i]] <- solution[(8+(19*(i-1))):(11+(19*(i-1)))]
    q1.pars <- matrix(nrow=5, ncol=5)
    colnames(q1.pars) <- rownames(q1.pars) <- paste0("1", LETTERS[1:5])
    for(i in 1:5) q1.pars[(1*i),(1:5)[-i]] <- solution[(12+(19*(i-1))):(15+(19*(i-1)))]
    q01.pars <- matrix(nrow=5, ncol=5)
    colnames(q01.pars) <- rownames(q01.pars) <- paste0("01", LETTERS[1:5])
    for(i in 1:5) q01.pars[(1*i),(1:5)[-i]] <- solution[(16+(19*(i-1))):(19+(19*(i-1)))]
    ## Drop the columns and lines of the result that do not belong to the model.
    n.areas <- ncol( out$control$trans.matrix ) / 3
    geosse.pars <- geosse.pars[,1:n.areas]
    q0.pars <- q0.pars[1:n.areas,1:n.areas]
    q1.pars <- q1.pars[1:n.areas,1:n.areas]
    q01.pars <- q01.pars[1:n.areas,1:n.areas]

    ## List with the parameter estimates for the model.
    out$solution
    out$solution$model.pars <- geosse.pars
    out$solution$q.0 <- q0.pars
    out$solution$q.1 <- q1.pars
    out$solution$q.01 <- q01.pars

    return(out)
}
