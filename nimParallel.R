

library(parallel)
this_cluster <- makeCluster(2) # Use only 2 for simplicity

set.seed(10120)
# Simulate some data
myData <- rgamma(1000, shape = 0.4, rate = 0.8)

# Option 1: make global lists and use the X input from parLapply
# to look up elements.
run_MCMC_allcode <- function(X, data) {
    library(nimble)

    myCode <- nimbleCode({
        a ~ dunif(0, 100)
        b ~ dnorm(0, 100)

        for (i in 1:length_y) {
            y[i] ~ dgamma(shape = a, rate = b)
        }
    })

    inputs <- inputsList[[X]]
    inits <- inputs$inits
    seed <- inputs$seed

    myModel <- nimbleModel(code = myCode,
                           data = list(y = data),
                           constants = list(length_y = 1000),
                           inits = inits)

    CmyModel <- compileNimble(myModel)

    myMCMC <- buildMCMC(CmyModel)
    CmyMCMC <- compileNimble(myMCMC)

    results <- runMCMC(CmyMCMC, niter = 10000, setSeed = seed)

    return(list(results = results, seed = seed))
}

inputsList <- list(
    list(seed = 123, # The seeds are arbitrary numbers
         inits = list(a = 0.5, b = 0.5)),
    list(seed = 234,
         inits = list(a = 0.3, b = 0.3)))
clusterExport(this_cluster, c("inputsList"))
chain_output <- parLapply(cl = this_cluster, X = 1:2,
                          fun = run_MCMC_allcode,
                          data = myData)
plot(chain_output[[1]]$results[, 'a'])

# Option 2: Make the X input to parLapply be a list
# including elements such as seed and inits
#
run_MCMC_allcode <- function(X, data) {
    library(nimble)

    myCode <- nimbleCode({
        a ~ dunif(0, 100)
        b ~ dnorm(0, 100)

        for (i in 1:length_y) {
            y[i] ~ dgamma(shape = a, rate = b)
        }
    })

    inits <- X$inits
    myModel <- nimbleModel(code = myCode,
                           data = list(y = data),
                           constants = list(length_y = 1000),
                           inits = inits)

    CmyModel <- compileNimble(myModel)

    myMCMC <- buildMCMC(CmyModel)
    CmyMCMC <- compileNimble(myMCMC)

    seed <- X$seed
    results <- runMCMC(CmyMCMC, niter = 10000, setSeed = seed)

    return(list(results = results, seed = seed))
}

# Xlist is the same as inputsList above.
# I'm  just naming it Xlist to make clear it will be
# the value passed as the X argument.
Xlist <- list(
    list(seed = 123,
         inits = list(a = 0.5, b = 0.5)),
    list(seed = 234,
         inits = list(a = 0.3, b = 0.3)))
chain_output <- parLapply(cl = this_cluster, X = Xlist,
                          fun = run_MCMC_allcode,
                          data = myData)

plot(chain_output[[1]]$results[, 'a'])

stopCluster(this_cluster)



