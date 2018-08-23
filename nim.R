library(nimble)

code <- nimbleCode({
    a ~ dnorm(0, 1)
})
constants <- list()
data <- list()
inits <- list(a = 0)

Rmodel <- nimbleModel(code, constants, data, inits)
Rmodel$calculate()

conf <- configureMCMC(Rmodel)
conf$printSamplers()
Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)#, showCompilerOutput = TRUE)
##compiledList <- compileNimble(list(model=Rmodel, mcmc=Rmcmc))
##Cmodel <- compiledList$model; Cmcmc <- compiledList$mcmc

set.seed(0)
samples <- runMCMC(Cmcmc, 10000)

colnames(samples)
apply(samples, 2, mean)

samplesSummary(samples)
samplesPlot(samples)

library(coda)
apply(samples, 2, effectiveSize)


nfDef <- nimbleFunction(
    setup = function() {},
    run = function() {
        returnType()
    }
)

Rnf <- nfDef()
Cnf <- compileNimble(Rnf)#, showCompilerOutput = TRUE)

Rnf$run()
Cnf$run()


Rnf <- nimbleFunction(
    run = function() {
        returnType()
    }
)

Cnf <- compileNimble(Rnf)#, showCompilerOutput = TRUE)

Rnf()
Cnf()

