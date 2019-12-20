## Exporting the likelihood for the GeoHiSSE function.
## This is based on the fGeoHiSSE version of the function.

makeGeoHiSSELikelihood <- function(phy, data, hidden.areas=0, f = c(1,1,1), assume.cladogenetic=TRUE, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, dt.threads=1, ode.eps = 0, bad.likelihood = exp(-500)){
    
    ## ########################################################
    ## Block for basic checks.
        
    ## Temporary fix for the current BUG:
    if( !is.null(phy$node.label) ) phy$node.label <- NULL

    ## The maximum number of hidden areas is 9
    if( hidden.areas > 9 ){
        stop("Maximum of hidden areas is 9.")
    }    
    
    if(!is.null(root.p)) {
        ## The vector of the root.p need to be as long as the speciation vector.
        ## Check the length and do the correction to sum to 1.
        ## All good:
        root.p <- root.p / sum(root.p)
    }

    ## Needed for the data.table package.
    setDTthreads(threads = dt.threads)

    if(!root.type == "madfitz" & !root.type == "herr_als"){
        stop("Check that you specified a proper root.type option. Options are 'madfitz' or 'herr_als'. See help for more details.", call.=FALSE)
    }
        
    ## Return error message if the data is not in the correct format.
    if( !inherits(data, what = c("matrix","data.frame")) ){
        stop("'data' needs to be a matrix or data.frame with 2 columns. See help.")
    }
    if( !ncol( data ) == 2 ){
        stop("'data' needs to be a matrix or data.frame with 2 columns. See help.")
    }
    ## Check if the states are %in% 0:2:
    states.check <- all( as.numeric(data[,2]) %in% 0:2 )
    if( !states.check ){
        stop("states need to be one of 0, 1, or 2. See help.")
    }

    ## END: Block for basic checks.
    ## ########################################################

    ## This is used to scale starting values to account for sampling:
    if( length(f) == 3 ){
        ## All good.
    } else{    
        if(length(f) == Ntip(phy)){
            stop("This is functionality has been temporarily removed.")
            ## samp.freq.tree <- Ntip(phy) / sum(table(data.new[,1]) / mean(f[as.numeric(names(freqs))+1]))
        }else{
            stop("The vector of sampling frequencies does not match the number of tips in the tree.")
        }
    }
    
    ## Generate the parameter key:
    par_key <- c("tau00A","tau11A","tau01A","ef00A","ef11A","d00A_11A","d00A_01A","d11A_00A","d11A_01A","d01A_00A","d01A_11A","d00A_00B","d00A_00C","d00A_00D","d00A_00E","d00A_00F","d00A_00G","d00A_00H","d00A_00I","d00A_00J","d11A_11B","d11A_11C","d11A_11D","d11A_11E","d11A_11F","d11A_11G","d11A_11H","d11A_11I","d11A_11J","d01A_01B","d01A_01C","d01A_01D","d01A_01E","d01A_01F","d01A_01G","d01A_01H","d01A_01I","d01A_01J","tau00B","tau11B","tau01B","ef00B","ef11B","d00B_11B","d00B_01B","d11B_00B","d11B_01B","d01B_00B","d01B_11B","d00B_00A","d00B_00C","d00B_00D","d00B_00E","d00B_00F","d00B_00G","d00B_00H","d00B_00I","d00B_00J","d11B_11A","d11B_11C","d11B_11D","d11B_11E","d11B_11F","d11B_11G","d11B_11H","d11B_11I","d11B_11J","d01B_01A","d01B_01C","d01B_01D","d01B_01E","d01B_01F","d01B_01G","d01B_01H","d01B_01I","d01B_01J","tau00C","tau11C","tau01C","ef00C","ef11C","d00C_11C","d00C_01C","d11C_00C","d11C_01C","d01C_00C","d01C_11C","d00C_00A","d00C_00B","d00C_00D","d00C_00E","d00C_00F","d00C_00G","d00C_00H","d00C_00I","d00C_00J","d11C_11A","d11C_11B","d11C_11D","d11C_11E","d11C_11F","d11C_11G","d11C_11H","d11C_11I","d11C_11J","d01C_01A","d01C_01B","d01C_01D","d01C_01E","d01C_01F","d01C_01G","d01C_01H","d01C_01I","d01C_01J","tau00D","tau11D","tau01D","ef00D","ef11D","d00D_11D","d00D_01D","d11D_00D","d11D_01D","d01D_00D","d01D_11D","d00D_00A","d00D_00B","d00D_00C","d00D_00E","d00D_00F","d00D_00G","d00D_00H","d00D_00I","d00D_00J","d11D_11A","d11D_11B","d11D_11C","d11D_11E","d11D_11F","d11D_11G","d11D_11H","d11D_11I","d11D_11J","d01D_01A","d01D_01B","d01D_01C","d01D_01E","d01D_01F","d01D_01G","d01D_01H","d01D_01I","d01D_01J","tau00E","tau11E","tau01E","ef00E","ef11E","d00E_11E","d00E_01E","d11E_00E","d11E_01E","d01E_00E","d01E_11E","d00E_00A","d00E_00B","d00E_00C","d00E_00D","d00E_00F","d00E_00G","d00E_00H","d00E_00I","d00E_00J","d11E_11A","d11E_11B","d11E_11C","d11E_11D","d11E_11F","d11E_11G","d11E_11H","d11E_11I","d11E_11J","d01E_01A","d01E_01B","d01E_01C","d01E_01D","d01E_01F","d01E_01G","d01E_01H","d01E_01I","d01E_01J","tau00F","tau11F","tau01F","ef00F","ef11F","d00F_11F","d00F_01F","d11F_00F","d11F_01F","d01F_00F","d01F_11F","d00F_00A","d00F_00B","d00F_00C","d00F_00D","d00F_00E","d00F_00G","d00F_00H","d00F_00I","d00F_00J","d11F_11A","d11F_11B","d11F_11C","d11F_11D","d11F_11E","d11F_11G","d11F_11H","d11F_11I","d11F_11J","d01F_01A","d01F_01B","d01F_01C","d01F_01D","d01F_01E","d01F_01G","d01F_01H","d01F_01I","d01F_01J","tau00G","tau11G","tau01G","ef00G","ef11G","d00G_11G","d00G_01G","d11G_00G","d11G_01G","d01G_00G","d01G_11G","d00G_00A","d00G_00B","d00G_00C","d00G_00D","d00G_00E","d00G_00F","d00G_00H","d00G_00I","d00G_00J","d11G_11A","d11G_11B","d11G_11C","d11G_11D","d11G_11E","d11G_11F","d11G_11H","d11G_11I","d11G_11J","d01G_01A","d01G_01B","d01G_01C","d01G_01D","d01G_01E","d01G_01F","d01G_01H","d01G_01I","d01G_01J","tau00H","tau11H","tau01H","ef00H","ef11H","d00H_11H","d00H_01H","d11H_00H","d11H_01H","d01H_00H","d01H_11H","d00H_00A","d00H_00B","d00H_00C","d00H_00D","d00H_00E","d00H_00F","d00H_00G","d00H_00I","d00H_00J","d11H_11A","d11H_11B","d11H_11C","d11H_11D","d11H_11E","d11H_11F","d11H_11G","d11H_11I","d11H_11J","d01H_01A","d01H_01B","d01H_01C","d01H_01D","d01H_01E","d01H_01F","d01H_01G","d01H_01I","d01H_01J","tau00I","tau11I","tau01I","ef00I","ef11I","d00I_11I","d00I_01I","d11I_00I","d11I_01I","d01I_00I","d01I_11I","d00I_00A","d00I_00B","d00I_00C","d00I_00D","d00I_00E","d00I_00F","d00I_00G","d00I_00H","d00I_00J","d11I_11A","d11I_11B","d11I_11C","d11I_11D","d11I_11E","d11I_11F","d11I_11G","d11I_11H","d11I_11J","d01I_01A","d01I_01B","d01I_01C","d01I_01D","d01I_01E","d01I_01F","d01I_01G","d01I_01H","d01I_01J","tau00J","tau11J","tau01J","ef00J","ef11J","d00J_11J","d00J_01J","d11J_00J","d11J_01J","d01J_00J","d01J_11J","d00J_00A","d00J_00B","d00J_00C","d00J_00D","d00J_00E","d00J_00F","d00J_00G","d00J_00H","d00J_00I","d11J_11A","d11J_11B","d11J_11C","d11J_11D","d11J_11E","d11J_11F","d11J_11G","d11J_11H","d11J_11I","d01J_01A","d01J_01B","d01J_01C","d01J_01D","d01J_01E","d01J_01F","d01J_01G","d01J_01H","d01J_01I")    

    ## Need a boolean for the presence of hidden states:
    if( hidden.areas > 0 ){ 
        hidden.states <- TRUE
    } else{
        hidden.states <- FALSE
    }
    
    data.new <- data.frame(data[,2], data[,2], row.names=data[,1])
    data.new <- data.new[phy$tip.label,]
    
    ## Some new prerequisites ##
    gen <- FindGenerations(phy)
    dat.tab <- OrganizeDataGeo(data=data.new[,1], phy=phy, f=f, hidden.states=hidden.states)
    nb.tip <- Ntip(phy)
    nb.node <- phy$Nnode
    
    if( hidden.areas == 0 ){
        names.pars <- paste("tau", c("00A","11A","01A"), sep="")
        names.pars <- c(names.pars, paste("ef", c("00A","11A"), sep="") )
        names.pars <- c(names.pars, "d00A_11A","d00A_01A","d11A_00A","d11A_01A","d01A_00A","d01A_11A")
        ## Make an empty numeric vector to return to the user.
        initial.pars <- setNames(object = rep(0.0, times = length(names.pars))
                               , nm = names.pars)
        ## Generate the translation key for the parameters of this model.
        trans_key <- unname( sapply(names.pars, function(x) which( par_key == x ) ) )
        category.rates.unique <- 0 ## Unclear why this is necessary.
        
        loglik <- function(p){
            ## Make the log.lik function for this case here.
            p.full <- rep(0.0, times = length(par_key) )
            p.full[ trans_key ] <- p
            cache <- ParametersToPassGeoHiSSEfast(model.vec = p.full, hidden.states = hidden.states, assume.cladogenetic = assume.cladogenetic
                                                , nb.tip = nb.tip, nb.node = nb.node, bad.likelihood = bad.likelihood, ode.eps = ode.eps)
            if(any(cache$s01A<0, cache$s01B<0, cache$s01C<0, cache$s01D<0, cache$s01E<0, cache$s01F<0, cache$s01G<0, cache$s01H<0, cache$s01I<0, cache$s01J<0)){
                ## If this condition is meet than the likelihood will blow, the function will return the bad value for the likelihood.
                return(log(cache$bad.likelihood)^13)
            }else{
                logl <- DownPassGeoHissefast(dat.tab = dat.tab, cache = cache, gen = gen, condition.on.survival = condition.on.survival, root.type=root.type, root.p=root.p)
                return(logl)
            }
        }

        return( list( log.lik = loglik, pars = initial.pars ) )
    }

    ## If here, then we have n hidden states.
    rate.cats <- hidden.areas + 1
        ## names.pars <- paste("tau", c("00A","11A","01A", "00B","11B","01B"), sep="")
        names.pars <- paste("tau", c("00","11","01"), rep(LETTERS[1:rate.cats], each = 3), sep="")
        names.pars <- c(names.pars, paste("ef", c("00","11"), rep(LETTERS[1:rate.cats], each = 2), sep="") )
        base_string <- c("d00x_11x","d00x_01x","d11x_00x","d11x_01x","d01x_00x","d01x_11x")
        names.pars <- c(names.pars, c( sapply(LETTERS[1:rate.cats], function(str) gsub(pattern = "x", replacement = str, x = base_string) ) ) )
        get_string <- function( base ){
            base_string <- c("d00x_00y","d11x_11y","d01x_01y")
            sub_from <- gsub(pattern = "x", replacement = base[1], x = base_string)
            sub_to <- gsub(pattern = "y", replacement = base[2], x = sub_from)
            return( sub_to )
        }
        comb_mat <- expand.grid(LETTERS[1:rate.cats], LETTERS[1:rate.cats], stringsAsFactors = FALSE)
        comb_mat <- comb_mat[!comb_mat[,1] == comb_mat[,2],]
        names.pars <- c(names.pars, c( apply(comb_mat, 1, get_string) ) )
        ## Make an empty numeric vector to return to the user.
        initial.pars <- setNames(object = rep(0.0, times = length(names.pars))
                               , nm = names.pars)
        ## Generate the translation key for the parameters of this model.
        trans_key <- unname( sapply(names.pars, function(x) which( par_key == x ) ) )
        
        loglik <- function(p){
            ## Make the log.lik function for this case here.
            p.full <- rep(0.0, times = length(par_key) )
            p.full[ trans_key ] <- p
            cache <- ParametersToPassGeoHiSSEfast(model.vec = p.full, hidden.states = hidden.states, assume.cladogenetic = assume.cladogenetic
                                                , nb.tip = nb.tip, nb.node = nb.node, bad.likelihood = exp(-500), ode.eps = ode.eps)
            if(any(cache$s01A<0, cache$s01B<0, cache$s01C<0, cache$s01D<0, cache$s01E<0, cache$s01F<0, cache$s01G<0, cache$s01H<0, cache$s01I<0, cache$s01J<0)){
                ## If this condition is meet than the likelihood will blow, the function will return the bad value for the likelihood.
                return(log(cache$bad.likelihood)^13)
            }else{
                logl <- DownPassGeoHissefast(dat.tab = dat.tab, cache = cache, gen = gen, condition.on.survival = condition.on.survival, root.type=root.type, root.p=root.p)
                return(logl)
            }
        }
    
    return( list( log.lik = loglik, pars = initial.pars ) )
}

