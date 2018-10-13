## Common function used by multiple methods.

GenerateValues <- function(par, lower, upper, scale.int, max.tries=100, expand.prob=0, examined.max, examined.min) {
    pass=FALSE
    tries=0
    while(!pass && tries<=max.tries) {
        tries <- tries+1
        pass=TRUE
        new.vals <- rep(NA, length(par))
        for(i in sequence(length(par))) {
            examined.max[i] <- max(0.001, examined.max[i])
            min.val <- min(max(lower[i], (1-scale.int)*examined.min[i]), examined.max[i]) #just in case min is greater than max
            max.val <- max(min(upper[i], (1+scale.int)*examined.max[i]), examined.min[i])
            if(isTRUE(all.equal(min.val, max.val))) {
                min.val <- min.val * 0.9999
                max.val <- max.val * 1.0001
            }
            new.vals[i] <- runif(1, min.val, max.val)
            if(new.vals[i]<lower[i]) {
                pass=FALSE
            }
            if(new.vals[i]>upper[i]) {
                pass=FALSE
            }
        }
    }
    if(tries>max.tries) {
        return(NA)
    }
    return(new.vals)
}
